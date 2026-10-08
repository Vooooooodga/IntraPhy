"""Profile a shared foreground rate multiplier with family-specific rates."""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
import math

import numpy as np
from scipy.optimize import minimize

from ..structure.tree_context import canonical_tree
from .genomic_exon_rates import FamilyRateLikelihood, fit_family_rate
from .locus_rates import _local_curvature


def _tree_signature(tree):
    return tuple(sorted((node, tree.parent[node], tree.label[node], tree.length[node])
                        for node in tree.parent))


def _validate_roster(family_units, foreground, branch_length_mode):
    if branch_length_mode not in {"supplied", "unit"}:
        raise ValueError("branch_length_mode must be 'supplied' or 'unit'")
    if not isinstance(family_units, dict):
        raise ValueError("family_units must map family IDs to unit lists")
    if isinstance(foreground, (str, bytes)):
        raise ValueError("foreground must be an iterable of canonical child IDs")
    foreground = frozenset(str(node) for node in foreground)
    families = {}
    for family, units in family_units.items():
        label = str(family)
        if not label or label in families:
            raise ValueError("family IDs must be nonempty and unique after string conversion")
        families[label] = tuple(units)
    reference = next((unit for units in families.values() for unit in units), None)
    if reference is None:
        return families, None, foreground, 0., 0., 0, 0

    tree = reference["tree"]
    base_context = canonical_tree(tree)
    base_tree = base_context.tree
    base_signature = _tree_signature(base_tree)
    children = set(base_tree.parent) - {base_tree.root}
    if not foreground or not foreground <= children:
        raise ValueError("foreground must select canonical non-root branch children")
    selected_context = canonical_tree(tree, foreground)
    if _tree_signature(selected_context.tree) != base_signature:
        raise ValueError("foreground selection changes the canonical tree or origin opportunities")
    if selected_context.foreground != foreground:
        raise ValueError("foreground must use canonical biological branch identifiers")
    tips_expected = set(base_tree.leaf_by_label)

    for family, units in families.items():
        for unit in units:
            other = unit.get("tree")
            if other is None or _tree_signature(canonical_tree(other).tree) != base_signature:
                raise ValueError("all units must retain the same canonical tree and branch lengths")
            tips = unit.get("tips")
            if not isinstance(tips, dict) or set(tips) != tips_expected:
                raise ValueError("all units must retain the same complete tip roster")

    def exposure(children_subset):
        if branch_length_mode == "unit":
            return float(len(children_subset))
        return math.fsum(base_tree.branch_length(child) for child in children_subset)

    foreground_exposure = exposure(foreground)
    background_exposure = exposure(children - foreground)
    if foreground_exposure <= 0 or background_exposure <= 0:
        raise ValueError("foreground and background must both have positive exposure")
    return (families, base_tree, foreground, foreground_exposure,
            background_exposure, len(foreground), len(children - foreground))


def _family_record(fit, multiplier):
    return {"mu": fit.get("mu"), "log_likelihood": fit.get("log_likelihood"),
            "status": fit.get("status", "unknown"),
            "diagnostic": {key: fit[key] for key in (
                "reason", "upper_tail_diagnostic", "selected_candidate",
                "optimizer_candidates", "starting_log_likelihoods") if key in fit},
            "foreground_multiplier": float(multiplier)}


def compare_foreground(family_units, foreground, *, branch_length_mode="supplied",
                       max_origins=None, profile_multipliers=(), threads=1):
    """Compare null multiplier 1 with a shared, nonnegative foreground multiplier.

    Every family keeps its input units at every profile point and has its own
    nonnegative background rate. The likelihood remains the native conditional
    composite likelihood, marginalized over the existing origin assignments.
    """
    if type(threads) is not int or threads < 1:
        raise ValueError("threads must be a positive integer")
    requested = []
    for value in profile_multipliers:
        number = float(value)
        if not math.isfinite(number) or number < 0:
            raise ValueError("profile multipliers must be finite and nonnegative")
        requested.append(number)
    (families, tree, foreground, foreground_exposure, background_exposure,
     foreground_branches, background_branches) = _validate_roster(
        family_units, foreground, branch_length_mode)

    counts = {"requested_families": len(families),
        "eligible_families": sum(bool(units) for units in families.values()),
        "excluded_families": [family for family, units in sorted(families.items()) if not units],
        "units": sum(len(units) for units in families.values()),
        "informative_units": sum(len(units) for units in families.values()),
        "foreground_branches": foreground_branches,
        "background_branches": background_branches,
        "foreground_exposure": foreground_exposure,
        "background_exposure": background_exposure}

    if tree is None:
        def empty_family_records(multiplier):
            return {family: {"mu": None, "log_likelihood": None,
                "status": "no_estimable_units", "diagnostic": {},
                "foreground_multiplier": multiplier} for family in sorted(families)}
        def empty_point(multiplier):
            return {"foreground_multiplier": (None if multiplier is None else float(multiplier)),
                "log_likelihood": None, "status": "no_estimable_units",
                "families": empty_family_records(multiplier),
                "unresolved_families": list(sorted(families))}
        points = sorted({0., 1., *requested, .01, .1, 10., 100.})
        null = empty_point(1.)
        return {"schema": "intraphy.exon-foreground-comparison/1",
            "model": "exon-structure-ctmc", "scope": "no_estimable_families",
            "foreground": {"canonical_children": sorted(foreground),
                "branch_length_mode": branch_length_mode}, "counts": counts,
            "null": null,
            "alternative": {**empty_point(None), "families": empty_family_records(None),
                "best_evaluated_multiplier": None},
            "profile": [empty_point(point) for point in points],
            "comparison_statistic": None, "p_value": None,
            "calibration_status": "not_calibrated_composite_likelihood"}

    executor = ThreadPoolExecutor(max_workers=threads) if threads > 1 else None
    unit_map = map if executor is None else executor.map
    workspaces = {family: FamilyRateLikelihood(units, branch_length_mode, max_origins,
        unit_map=unit_map, family_id=family) for family, units in families.items() if units}
    evaluations = {}

    def evaluate_multiplier(multiplier):
        key = float(multiplier)
        if key in evaluations:
            return evaluations[key]
        family_records = {}
        failed = []
        impossible = []
        terms = []
        for family, units in sorted(families.items()):
            if not units:
                family_records[family] = {"mu": None, "log_likelihood": None,
                    "status": "no_estimable_units", "diagnostic": {},
                    "foreground_multiplier": key}
                continue
            try:
                fit = fit_family_rate(units, branch_length_mode, max_origins,
                    family_id=family, foreground=foreground,
                    foreground_multiplier=key, likelihood_workspace=workspaces[family])
            except (ValueError, ArithmeticError) as exc:
                reason = str(exc)
                status = ("origin_scenarios_incomplete" if "origin_scenarios_incomplete" in reason
                          else "likelihood_evaluation_failed")
                fit = {"status": status, "converged": False,
                       "mu": None, "reason": f"{type(exc).__name__}: {exc}"}
            record = _family_record(fit, key)
            family_records[family] = record
            if fit.get("status") == "no_rate_information" and fit.get("log_likelihood") is not None:
                terms.append(float(fit["log_likelihood"]))
            elif fit.get("status") == "impossible_observation_under_model":
                impossible.append(family)
            elif fit.get("converged") and fit.get("log_likelihood") is not None and fit.get("mu") is not None:
                terms.append(float(fit["log_likelihood"]))
            else:
                failed.append(family)
        eligible = any(families.values())
        if failed:
            total = None
            status = "profile_incomplete"
        elif impossible:
            total = None
            status = "impossible_observation_under_model"
        elif not eligible:
            total = None
            status = "no_estimable_units"
        else:
            total = float(sum(terms))
            status = "conditional_composite_log_likelihood"
        result = {"foreground_multiplier": key, "log_likelihood": total,
            "status": status, "families": family_records,
            "unresolved_families": failed, "impossible_families": impossible}
        evaluations[key] = result
        return result

    try:
        profile_points = {0., 1., *requested, .01, .1, 10., 100.}
        for point in sorted(profile_points):
            evaluate_multiplier(point)

        starts = (.01, .1, 1., 10., 100.)
        optimization_evaluations = set(starts)
        candidates = []
        for start in starts:
            def objective(point):
                multiplier = float(point[0])
                optimization_evaluations.add(multiplier)
                log_likelihood = evaluate_multiplier(multiplier)["log_likelihood"]
                return 1e300 if log_likelihood is None else -log_likelihood
            result = minimize(objective, np.array([start]), method="L-BFGS-B",
                              bounds=[(0., None)])
            multiplier = float(result.x[0])
            point = evaluate_multiplier(multiplier)
            candidates.append({"foreground_multiplier": multiplier,
                "success": bool(result.success), "message": str(result.message),
                "status_code": int(result.status), "iterations": int(result.nit),
                "log_likelihood": point["log_likelihood"], "profile_status": point["status"]})
        zero = evaluate_multiplier(0.)
        if zero["log_likelihood"] is not None:
            candidates.append({"foreground_multiplier": 0., "success": True,
                "message": "explicit_zero_boundary", "status_code": 0,
                "iterations": 0, "log_likelihood": zero["log_likelihood"]})
        null_point = evaluate_multiplier(1.)
        if null_point["log_likelihood"] is not None:
            candidates.append({"foreground_multiplier": 1., "success": True,
                "message": "explicit_nested_null", "status_code": 0,
                "iterations": 0, "log_likelihood": null_point["log_likelihood"]})
        successful = [candidate for candidate in candidates
                      if candidate["success"] and candidate["log_likelihood"] is not None]
        best = max(successful, key=lambda item: item["log_likelihood"]) if successful else None
        alternative = {"foreground_multiplier": None, "best_evaluated_multiplier": None,
            "log_likelihood": None, "status": "profile_incomplete",
            "families": evaluate_multiplier(1.)["families"],
            "optimizer_candidates": candidates, "resolution_reasons": ["no_successful_candidate"],
            "incomplete_profile_points": []}
        best_observed_multiplier = None
        if best is not None:
            rho_hat = best["foreground_multiplier"]
            best_point = evaluate_multiplier(rho_hat)
            positive_explored = [value for value in evaluations if value > 0]
            probe = max(2. * max(rho_hat, max(positive_explored)), 1.)
            optimization_evaluations.add(probe)
            probe_point = evaluate_multiplier(probe)
            tolerance = 1e-10 * max(1., abs(best_point["log_likelihood"]))
            if probe_point["status"] == "profile_incomplete":
                tail_reason = "probe_fit_incomplete"
            elif probe_point["status"] == "impossible_observation_under_model":
                tail_reason = "lower_probe"
            elif probe_point["log_likelihood"] > best_point["log_likelihood"] + tolerance:
                tail_reason = "higher_probe"
            elif probe_point["log_likelihood"] >= best_point["log_likelihood"] - tolerance:
                tail_reason = "numerically_indistinguishable_probe"
            else:
                tail_reason = "lower_probe"
            diagnostic_points = {value for value in (.01, .1, 1., 10., 100., rho_hat, probe)
                                 if value > 0}
            positive_ll = [evaluate_multiplier(value)["log_likelihood"]
                for value in diagnostic_points
                if evaluate_multiplier(value)["log_likelihood"] is not None]
            flat_profile = (len(positive_ll) >= 2 and
                max(positive_ll) - min(positive_ll) <= tolerance)
            if rho_hat > 1e-8:
                def curvature_objective(point):
                    multiplier = float(point[0])
                    optimization_evaluations.add(multiplier)
                    log_likelihood = evaluate_multiplier(multiplier)["log_likelihood"]
                    return 1e300 if log_likelihood is None else -log_likelihood
                try:
                    curvature = _local_curvature(curvature_objective,
                        np.array([rho_hat]), np.array([False]))
                    _, rank, dimension, positive, negative, flat = curvature
                except (ValueError, ArithmeticError):
                    rank, dimension, positive, negative, flat = 0, 1, False, False, True
            else:
                rank, dimension, positive, negative, flat = None, None, None, None, False
            curvature_unresolved = (rho_hat > 1e-8 and
                (flat or rank < dimension or not positive or negative))
            higher_evaluated = [point for point in evaluations.values()
                if point["log_likelihood"] is not None
                and point["log_likelihood"] > best_point["log_likelihood"] + tolerance]
            failed_higher = [candidate for candidate in candidates
                if not candidate["success"] and candidate["log_likelihood"] is not None
                and candidate["log_likelihood"] > best_point["log_likelihood"] + tolerance]
            null_exceeded = (null_point["log_likelihood"] is not None and
                best_point["log_likelihood"] < null_point["log_likelihood"] - tolerance)
            incomplete_search = sorted(value for value in optimization_evaluations
                if evaluate_multiplier(value)["status"] == "profile_incomplete")
            unidentified = (tail_reason == "numerically_indistinguishable_probe"
                            or flat_profile or curvature_unresolved)
            optimization_unresolved = bool(higher_evaluated or failed_higher or null_exceeded)
            tail_unresolved = tail_reason == "higher_probe"
            null_incomplete = null_point["status"] == "profile_incomplete"
            reasons = []
            if incomplete_search or null_incomplete or tail_reason == "probe_fit_incomplete":
                reasons.append("profile_incomplete")
            if tail_unresolved:
                reasons.append("upper_tail_unresolved")
            if optimization_unresolved:
                reasons.append("optimization_unresolved")
            if unidentified:
                reasons.append("foreground_multiplier_unidentified")
            if not reasons:
                reasons.append("resolved")
            best_observed = max((point for point in evaluations.values()
                if point["log_likelihood"] is not None),
                key=lambda point: point["log_likelihood"])
            best_observed_multiplier = best_observed["foreground_multiplier"]
            if reasons[0] != "resolved":
                status = reasons[0]
                reported_rho = None
            elif rho_hat <= 1e-8:
                status = "zero_boundary"
                reported_rho = 0.
            else:
                status = "estimated_conditional_composite_multiplier"
                reported_rho = rho_hat
            alternative = {"foreground_multiplier": reported_rho,
                "best_evaluated_multiplier": (best_observed["foreground_multiplier"]
                    if reported_rho is None else rho_hat),
                "log_likelihood": (best_observed["log_likelihood"] if reported_rho is None
                    else best_point["log_likelihood"]), "status": status,
                "families": best_observed["families"] if reported_rho is None else best_point["families"],
                "optimizer_candidates": candidates,
                "resolution_reasons": reasons,
                "incomplete_profile_points": incomplete_search,
                "upper_tail_diagnostic": {"method": "finite_profile_probe",
                    "probe_multiplier": probe, "probe_log_likelihood": probe_point["log_likelihood"],
                    "best_log_likelihood": (best_observed["log_likelihood"] if reported_rho is None
                        else best_point["log_likelihood"]),
                    "reason": tail_reason,
                    "scope": "finite_probe_only_not_an_asymptotic_or_global_optimum_test"},
                "profile_curvature": {"rank": rank, "dimension": dimension,
                    "positive_definite": positive, "has_negative_eigenvalue": negative,
                    "locally_flat": flat}}
        else:
            fallback = max((point for point in evaluations.values()
                if point["log_likelihood"] is not None),
                key=lambda point: point["log_likelihood"], default=evaluate_multiplier(1.))
            incomplete_search = sorted(value for value in optimization_evaluations
                if evaluate_multiplier(value)["status"] == "profile_incomplete")
            status = "profile_incomplete"
            alternative = {"foreground_multiplier": None,
                "best_evaluated_multiplier": (fallback["foreground_multiplier"]
                    if fallback["log_likelihood"] is not None else None),
                "log_likelihood": fallback["log_likelihood"], "status": status,
                "families": fallback["families"], "optimizer_candidates": candidates,
                "resolution_reasons": [status],
                "incomplete_profile_points": incomplete_search}

        null = evaluate_multiplier(1.)
        profile_keys = {0., 1., *requested, .01, .1, 10., 100.}
        if best is not None:
            profile_keys.add(best["foreground_multiplier"])
            profile_keys.add(probe)
            profile_keys.add(best_observed_multiplier)
        profile = [evaluations[key] for key in sorted(profile_keys) if key in evaluations]
        statistic = None
        if (null["log_likelihood"] is not None and alternative["log_likelihood"] is not None
                and alternative["status"] in {"estimated_conditional_composite_multiplier", "zero_boundary"}):
            statistic = 2. * (alternative["log_likelihood"] - null["log_likelihood"])
        scope = "single_gene" if len([units for units in families.values() if units]) == 1 else "multi_gene"
        return {"schema": "intraphy.exon-foreground-comparison/1",
            "model": "exon-structure-ctmc", "scope": scope,
            "foreground": {"canonical_children": sorted(foreground),
                "branch_length_mode": branch_length_mode,
                "foreground_exposure": foreground_exposure,
                "background_exposure": background_exposure,
                "origin_prior_unchanged": True},
            "counts": counts, "null": null, "alternative": alternative,
            "profile": profile, "comparison_statistic": statistic, "p_value": None,
            "calibration_status": "not_calibrated_composite_likelihood"}
    finally:
        if executor is not None:
            executor.shutdown(wait=True)
