"""Shared one-parameter conditional composite fit for genomic exon spans."""
from __future__ import annotations

import math
import numpy as np
from scipy.optimize import minimize

from ..structure.edits import EDIT_KINDS
from .configuration_model import RateModel, evaluate_model
from .locus_rates import _local_curvature


def _model(mu):
    return RateModel({kind: float(mu) for kind in EDIT_KINDS}, scale=1.)


def fit_family_rate(units, branch_length_mode="supplied", max_origins=None):
    """Fit a nonnegative rate per opportunity and branch length unit."""
    if not units:
        return {"status": "no_estimable_units", "converged": False, "mu": None}
    tree = units[0]["tree"]
    exposure = (sum(1. for _ in tree.edges()) if branch_length_mode == "unit" else
                sum(tree.branch_length(child) for _, child in tree.edges()))
    if exposure <= 0:
        return {"status": "flat_or_nonidentified", "converged": False, "mu": None,
                "reason": "tree_has_no_positive_branch_length"}

    def log_likelihood(x):
        model = _model(float(x)/exposure)
        values = [float(evaluate_model(unit["space"], unit["tree"], unit["tips"], model,
            max_origins=max_origins, posterior=False, counts=False,
            branch_length_mode=branch_length_mode)["log_likelihood"]) for unit in units]
        return float(sum(values)) if all(np.isfinite(values)) else -math.inf

    def objective(point):
        value = log_likelihood(float(point[0]))
        return 1e300 if value == -math.inf else -value

    # x = mu * total tree exposure gives comparable starts across tree units.
    starts = (0., .01, .1, 1., 10., 100.)
    sampled = [(x, log_likelihood(x)) for x in starts]
    finite = [value for _, value in sampled if np.isfinite(value)]
    if not finite:
        return {"status": "impossible_observation_under_model", "converged": False,
                "mu": None, "starting_log_likelihoods": sampled}
    if len(finite) == len(starts) and max(finite)-min(finite) <= 1e-10:
        return {"status": "no_rate_information", "converged": False, "mu": None,
                "log_likelihood": max(finite), "starting_log_likelihoods": sampled}

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
                "starting_log_likelihoods": sampled, "optimizer_candidates": candidates}
    best = max(successful, key=lambda candidate: candidate["log_likelihood"])
    tolerance = 1e-10 * max(1., abs(best["log_likelihood"]))
    higher_failed = [candidate for candidate in candidates
                     if not candidate["success"] and candidate["log_likelihood"] > best["log_likelihood"]+tolerance]
    if higher_failed:
        return {"status": "optimizer_unresolved_higher_failed_candidate", "converged": False,
                "mu": None, "log_likelihood": best["log_likelihood"],
                "selected_candidate": best, "higher_failed_candidates": higher_failed,
                "optimizer_candidates": candidates, "starting_log_likelihoods": sampled}

    xhat = best["x"]
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
        "fit_scope": "best_converged_optimizer_result_among_tree_exposure_scaled_starts",
        "zero_boundary_log_likelihood": zero_ll,
        "optimizer": "scipy.optimize.minimize(method='L-BFGS-B', bounds=[(0, None)])"}


def fixed_family_fit(units, rates, branch_length_mode="supplied", max_origins=None):
    values = [float(evaluate_model(unit["space"], unit["tree"], unit["tips"], rates,
        max_origins=max_origins, posterior=False, counts=False,
        branch_length_mode=branch_length_mode)["log_likelihood"]) for unit in units]
    ll = float(sum(values)) if values and all(np.isfinite(values)) else -math.inf
    return {"status": "fixed_parameters", "converged": True,
        "log_likelihood": ll if np.isfinite(ll) else None,
        "unit_log_likelihoods": values,
        "rate_parameters": {"rates": rates.rates, "scale": rates.scale,
            "foreground_multiplier": rates.foreground_multiplier,
            "foreground": sorted(rates.foreground), "origin_root_weight": rates.origin_root_weight},
        "rate_unit": "per_eligible_elementary_edit_opportunity_per_branch_length_unit",
        "objective": "conditional_composite_log_likelihood_across_local_units_within_family"}
