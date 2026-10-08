"""Shared one-parameter conditional composite fit for genomic exon spans."""
from __future__ import annotations

import math
import logging
from collections import OrderedDict
from threading import Lock
from time import perf_counter
import numpy as np
from scipy.optimize import minimize

from ..structure.edits import EDIT_KINDS
from .configuration_model import RateModel, evaluate_model
from .configuration_fit_cache import CommonRateKernels
from .kernel_cache import KernelCache
from .locus_rates import _local_curvature


logger = logging.getLogger("intraphy")


def _model(mu, foreground=frozenset(), foreground_multiplier=1.):
    return RateModel({kind: float(mu) for kind in EDIT_KINDS}, scale=1.,
                     foreground=frozenset(foreground),
                     foreground_multiplier=float(foreground_multiplier))


class FamilyRateLikelihood:
    """Memoized likelihood workspace for one fixed family roster."""

    def __init__(self, units, branch_length_mode="supplied", max_origins=None, *,
                 unit_map=None, family_id=None):
        self.units = tuple(units)
        self.branch_length_mode = branch_length_mode
        self.max_origins = max_origins
        self.family_id = family_id
        if branch_length_mode not in {"supplied", "unit"}:
            raise ValueError("branch_length_mode must be 'supplied' or 'unit'")
        tree = self.units[0]["tree"] if self.units else None
        self.exposure = (sum(1. for _ in tree.edges()) if tree is not None and branch_length_mode == "unit" else
                         sum(tree.branch_length(child) for _, child in tree.edges())
                         if tree is not None else 0.)
        template_cache = KernelCache()
        cache_lock = Lock()
        self.templates = tuple(CommonRateKernels(unit["space"], template_cache, cache_lock)
                               for unit in self.units)
        self.ordered_map = map if unit_map is None else unit_map
        self._values = OrderedDict()

    def log_likelihood(self, x, foreground=frozenset(), foreground_multiplier=1.):
        key = (float(x), tuple(sorted(foreground)), float(foreground_multiplier))
        if key in self._values:
            self._values.move_to_end(key)
            return self._values[key]
        if self.exposure <= 0:
            return -math.inf
        started = perf_counter()
        model = _model(float(x) / self.exposure, foreground, foreground_multiplier)

        def evaluate_unit(item):
            index, (unit, template) = item
            unit_id = unit.get("unit_id", f"unit-{index}")
            unit_started = perf_counter()
            logger.info("Family unit likelihood started: family=%s unit=%s x=%.17g foreground_multiplier=%.17g",
                        self.family_id, unit_id, x, foreground_multiplier)
            try:
                value = float(evaluate_model(unit["space"], unit["tree"], unit["tips"], model,
                    max_origins=self.max_origins, posterior=False, counts=False,
                    branch_length_mode=self.branch_length_mode, backend="sparse",
                    sparse_templates=template)["log_likelihood"])
            except BaseException:
                logger.exception("Family unit likelihood failed: family=%s unit=%s x=%.17g",
                                 self.family_id, unit_id, x)
                raise
            logger.info("Family unit likelihood completed: family=%s unit=%s elapsed_seconds=%.6f",
                         self.family_id, unit_id, perf_counter() - unit_started)
            return value

        values = list(self.ordered_map(evaluate_unit, enumerate(zip(self.units, self.templates))))
        if any(math.isnan(value) or value == math.inf for value in values):
            raise ArithmeticError("Family likelihood contains NaN or positive infinity")
        if not values or any(value == -math.inf for value in values):
            result = -math.inf
        else:
            try:
                result = math.fsum(values)
            except OverflowError as exc:
                raise ArithmeticError("Family likelihood sum overflowed") from exc
            if not math.isfinite(result):
                raise ArithmeticError("Family likelihood sum is nonfinite")
        self._values[key] = result
        if len(self._values) > 128:
            self._values.popitem(last=False)
        logger.info("Family likelihood completed: family=%s x=%.17g foreground_multiplier=%.17g elapsed_seconds=%.6f",
                    self.family_id, x, foreground_multiplier, perf_counter() - started)
        return result


def fit_family_rate(units, branch_length_mode="supplied", max_origins=None, *,
                    unit_map=None, family_id=None, foreground=frozenset(),
                    foreground_multiplier=1., likelihood_workspace=None):
    """Fit a nonnegative rate per opportunity and branch length unit."""
    if not units:
        return {"status": "no_estimable_units", "converged": False, "mu": None}
    foreground = frozenset(foreground)
    workspace = likelihood_workspace or FamilyRateLikelihood(
        units, branch_length_mode, max_origins, unit_map=unit_map, family_id=family_id)
    if (len(units) != len(workspace.units)
            or any(unit is not cached for unit, cached in zip(units, workspace.units))
            or branch_length_mode != workspace.branch_length_mode
            or max_origins != workspace.max_origins):
        raise ValueError("likelihood_workspace does not match units or fit settings")
    base_exposure = workspace.exposure
    if not math.isfinite(base_exposure):
        raise ArithmeticError("Tree rate exposure is nonfinite")
    if base_exposure <= 0:
        return {"status": "flat_or_nonidentified", "converged": False, "mu": None,
                "rate_exposure": base_exposure,
                "reason": "tree_has_no_positive_branch_length"}
    multiplier = float(foreground_multiplier)
    if not math.isfinite(multiplier) or multiplier < 0:
        raise ValueError("foreground_multiplier must be finite and nonnegative")
    exposure = base_exposure
    if foreground and multiplier != 1.:
        tree = workspace.units[0]["tree"]
        measure = (lambda child: 1. if branch_length_mode == "unit"
                   else tree.branch_length(child))
        try:
            foreground_exposure = math.fsum(measure(child) for child in foreground)
            background_exposure = math.fsum(measure(child) for _, child in tree.edges()
                                            if child not in foreground)
            exposure = math.fsum((background_exposure,
                                  multiplier * foreground_exposure))
        except OverflowError as exc:
            raise ArithmeticError("Weighted rate exposure overflowed") from exc
        if not math.isfinite(exposure):
            raise ArithmeticError("Weighted rate exposure is nonfinite")
    if exposure <= 0:
        return {"status": "flat_or_nonidentified", "converged": False, "mu": None,
                "rate_exposure": exposure,
                "reason": "weighted_tree_has_no_positive_rate_exposure"}
    fit_scope = ("best_converged_optimizer_result_among_tree_exposure_scaled_starts"
        if exposure == base_exposure else
        "best_converged_optimizer_result_among_weighted_rate_exposure_scaled_starts")

    def log_likelihood(x):
        mu = float(x) / exposure
        workspace_x = (float(x) if exposure == base_exposure
                       else mu * base_exposure)
        if not math.isfinite(workspace_x):
            raise ArithmeticError("Workspace rate coordinate is nonfinite")
        return workspace.log_likelihood(workspace_x, foreground, multiplier)

    def objective(point):
        value = log_likelihood(float(point[0]))
        return 1e300 if value == -math.inf else -value

    # x = mu * rate exposure gives comparable starts across profile values.
    starts = (0., .01, .1, 1., 10., 100.)
    sampled = [(x, log_likelihood(x)) for x in starts]
    finite = [value for _, value in sampled if np.isfinite(value)]
    if not finite:
        return {"status": "impossible_observation_under_model", "converged": False,
                "mu": None, "rate_exposure": exposure,
                "starting_log_likelihoods": sampled}
    if len(finite) == len(starts) and max(finite)-min(finite) <= 1e-10:
        return {"status": "no_rate_information", "converged": False, "mu": None,
                "log_likelihood": max(finite), "rate_exposure": exposure,
                "starting_log_likelihoods": sampled}

    candidates = []
    for x0 in starts[1:]:
        result = minimize(objective, np.array([x0]), method="L-BFGS-B",
                          bounds=[(0., None)])
        value = log_likelihood(float(result.x[0]))
        if np.isfinite(value):
            candidates.append({"x": float(result.x[0]), "log_likelihood": value,
                "success": bool(result.success), "message": str(result.message),
                "status_code": int(result.status), "iterations": int(result.nit),
                "start": x0})
    zero_ll = log_likelihood(0.)
    if np.isfinite(zero_ll):
        candidates.append({"x": 0., "log_likelihood": zero_ll, "success": True,
            "message": "explicit_zero_boundary", "status_code": 0, "iterations": 0,
            "start": 0.})
    successful = [candidate for candidate in candidates if candidate["success"]]
    if not successful:
        return {"status": "optimizer_failed", "converged": False, "mu": None,
                "rate_exposure": exposure, "starting_log_likelihoods": sampled,
                "optimizer_candidates": candidates}
    best = max(successful, key=lambda candidate: candidate["log_likelihood"])
    tolerance = 1e-10 * max(1., abs(best["log_likelihood"]))
    higher_failed = [candidate for candidate in candidates
                     if not candidate["success"] and candidate["log_likelihood"] > best["log_likelihood"]+tolerance]
    if higher_failed:
        return {"status": "optimizer_unresolved_higher_failed_candidate", "converged": False,
                "mu": None, "rate_exposure": exposure,
                "log_likelihood": best["log_likelihood"],
                "selected_candidate": best, "higher_failed_candidates": higher_failed,
                "optimizer_candidates": candidates, "starting_log_likelihoods": sampled}

    xhat = best["x"]
    x_probe = 2. * max(xhat, max(starts))
    probe_ll = log_likelihood(x_probe)
    if not np.isfinite(probe_ll):
        probe_reason = "nonfinite_probe"
    elif probe_ll > best["log_likelihood"] + tolerance:
        probe_reason = "higher_probe"
    elif probe_ll >= best["log_likelihood"] - tolerance:
        probe_reason = "numerically_indistinguishable_probe"
    else:
        probe_reason = "lower_probe"
    upper_tail_diagnostic = {
        "method": "single_finite_probe_at_twice_max_selected_x_or_start",
        "reason": probe_reason,
        "x_probe": x_probe,
        "log_likelihood": probe_ll if np.isfinite(probe_ll) else None,
        "best_log_likelihood": best["log_likelihood"],
        "comparison_tolerance": tolerance,
        "scope": "finite_probe_only_not_an_asymptotic_or_global_optimum_test",
    }
    if probe_reason != "lower_probe":
        return {"status": "upper_tail_unresolved", "converged": False, "mu": None,
            "dimensionless_rate_x": xhat, "rate_exposure": exposure,
            "log_likelihood": best["log_likelihood"], "selected_candidate": best,
            "upper_tail_diagnostic": upper_tail_diagnostic,
            "optimizer_candidates": candidates, "starting_log_likelihoods": sampled,
            "rate_unit": "per_eligible_elementary_edit_opportunity_per_branch_length_unit",
            "objective": "conditional_composite_log_likelihood_across_local_units_within_family",
            "fit_scope": fit_scope}

    at_zero = xhat <= 1e-8
    curvature = _local_curvature(objective, np.array([xhat]), np.array([at_zero]))
    eigenvalues, rank, dimension, positive, negative, flat = curvature
    mu = xhat / exposure
    if at_zero:
        status = "zero_boundary"
    elif flat or rank < dimension or not positive or negative:
        status, mu = "flat_or_nonidentified", None
    else:
        status = "estimated_conditional_composite_rate"
    return {"status": status, "converged": True, "mu": mu,
        "dimensionless_rate_x": xhat, "rate_exposure": exposure,
        "log_likelihood": best["log_likelihood"], "selected_candidate": best,
        "optimizer_candidates": candidates, "starting_log_likelihoods": sampled,
        "curvature_eigenvalues": list(eigenvalues), "curvature_rank": rank,
        "curvature_dimension": dimension, "curvature_positive_definite": positive,
        "curvature_has_negative_eigenvalue": negative, "locally_flat": flat,
        "rate_unit": "per_eligible_elementary_edit_opportunity_per_branch_length_unit",
        "objective": "conditional_composite_log_likelihood_across_local_units_within_family",
        "fit_scope": fit_scope,
        "upper_tail_diagnostic": upper_tail_diagnostic,
        "zero_boundary_log_likelihood": zero_ll,
        "optimizer": "scipy.optimize.minimize(method='L-BFGS-B', bounds=[(0, None)])"}


def fixed_family_fit(units, rates, branch_length_mode="supplied", max_origins=None, *,
                     unit_map=None, family_id=None):
    ordered_map = map if unit_map is None else unit_map
    def evaluate_unit(item):
        index, unit = item
        unit_id = unit.get("unit_id", f"unit-{index}")
        started = perf_counter()
        logger.info("Family unit likelihood started: family=%s unit=%s x=fixed",
                    family_id, unit_id)
        try:
            value = float(evaluate_model(unit["space"], unit["tree"], unit["tips"], rates,
                max_origins=max_origins, posterior=False, counts=False,
                branch_length_mode=branch_length_mode, backend="sparse")["log_likelihood"])
        except BaseException:
            logger.exception("Family unit likelihood failed: family=%s unit=%s x=fixed elapsed_seconds=%.6f",
                             family_id, unit_id, perf_counter() - started)
            raise
        logger.info("Family unit likelihood completed: family=%s unit=%s x=fixed elapsed_seconds=%.6f",
                     family_id, unit_id, perf_counter() - started)
        return value
    values = list(ordered_map(evaluate_unit, enumerate(units)))
    ll = float(sum(values)) if values and all(np.isfinite(values)) else -math.inf
    return {"status": "fixed_parameters", "converged": True,
        "log_likelihood": ll if np.isfinite(ll) else None,
        "unit_log_likelihoods": values,
        "rate_parameters": {"rates": rates.rates, "scale": rates.scale,
            "foreground_multiplier": rates.foreground_multiplier,
            "foreground": sorted(rates.foreground), "origin_root_weight": rates.origin_root_weight},
        "rate_unit": "per_eligible_elementary_edit_opportunity_per_branch_length_unit",
        "objective": "conditional_composite_log_likelihood_across_local_units_within_family"}
