"""Family-shared two-state DNA presence CTMC fitting."""
from __future__ import annotations

import math
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from scipy.optimize import minimize

from intraphy.inference.ctmc import _ascertainment_log_probability
from intraphy.inference.ctmc import _compress_patterns
from intraphy.inference.ctmc import _pattern_log_likelihood
from intraphy.inference.locus_rates import _local_curvature


def _site_log_likelihood(tree, observations, gain, loss, root_presence):
    value = _pattern_log_likelihood(tree, observations, gain, loss, 1.0, frozenset(), root_presence)
    denominator = _ascertainment_log_probability(
        tree, observations, gain, loss, 1.0, frozenset(), root_presence,
        "observed-at-least-one",
    )
    if not math.isfinite(value) or not math.isfinite(denominator):
        return -math.inf
    return value - denominator


def _weighted_site_log_likelihood(args):
    tree, observations, weight, gain, loss, root_presence = args
    return int(weight) * _site_log_likelihood(tree, observations, gain, loss, root_presence)


def _fit_family(tree, patterns, root_frequency, root_presence, fixed_rates, workers):
    fixed = dict(fixed_rates or {})
    free = tuple(rate for rate in ("gain", "loss") if rate not in fixed)

    def rates_at(point):
        rates = dict(fixed)
        rates.update(zip(free, (float(value) for value in point)))
        return rates

    executor = ThreadPoolExecutor(max_workers=workers) if workers > 1 and len(patterns) > 1 else None

    def log_likelihood(point):
        rates = rates_at(point)
        rho = None if root_frequency == "stationary" else root_presence
        tasks = ((tree, observations, weight, rates["gain"], rates["loss"], rho)
                 for observations, weight in patterns)
        terms = list(executor.map(_weighted_site_log_likelihood, tasks)) if executor else [
            _weighted_site_log_likelihood(task) for task in tasks
        ]
        return math.fsum(terms) if all(math.isfinite(term) for term in terms) else -math.inf

    try:
        return _optimize_rates(tree, patterns, fixed, free, rates_at, log_likelihood)
    finally:
        if executor:
            executor.shutdown(wait=True, cancel_futures=True)


def _optimize_rates(tree, patterns, fixed, free, rates_at, log_likelihood):
    if not free:
        ll = log_likelihood(())
        return {
            "rates": rates_at(()), "log_likelihood": ll,
            "optimizer_success": math.isfinite(ll), "optimizer_status": 0,
            "optimizer_message": "all rates fixed" if math.isfinite(ll) else "fixed rates give zero likelihood",
            "iterations": 0, "boundary_rates": [name for name, value in fixed.items() if value == 0],
            "locally_flat": None, "curvature_eigenvalues": [], "selected_start": None,
        }

    tree_scale = math.fsum(tree.branch_length(child) for _parent, child in tree.edges())
    base = 1.0 / tree_scale if tree_scale > 0.0 else 1.0
    starts = [base * scale for scale in (0.1, 1.0, 10.0)]
    candidates = []
    records = []

    def objective(point):
        value = log_likelihood(point)
        return -value if math.isfinite(value) else 1e300

    for index, scale in enumerate(starts):
        start = np.full(len(free), scale, dtype=float)
        initial_ll = log_likelihood(start)
        result = minimize(
            objective, start, method="L-BFGS-B", bounds=[(0.0, None)] * len(free),
            options={"maxiter": 1000, "ftol": 1e-9, "gtol": 1e-6},
        )
        fitted_ll = log_likelihood(result.x)
        tolerance = 1e-10 * max(1.0, abs(initial_ll)) if math.isfinite(initial_ll) else 0.0
        success = bool(result.success and math.isfinite(fitted_ll) and fitted_ll >= initial_ll - tolerance)
        record = {
            "start_index": index, "start_rate": scale,
            "start_log_likelihood": initial_ll if math.isfinite(initial_ll) else None,
            "optimizer_success": success, "optimizer_status": int(result.status),
            "optimizer_message": str(result.message), "iterations": int(result.nit),
            "log_likelihood": fitted_ll if math.isfinite(fitted_ll) else None,
            "rates": rates_at(result.x),
        }
        records.append(record)
        if math.isfinite(fitted_ll):
            candidates.append((index, result, fitted_ll, success, record))

    if not candidates:
        return {
            "rates": {}, "log_likelihood": None, "optimizer_success": False,
            "optimizer_status": -1, "optimizer_message": "all starts had zero or non-finite likelihood",
            "iterations": 0, "boundary_rates": [], "locally_flat": None,
            "curvature_eigenvalues": [], "selected_start": None, "starts": records,
        }
    successful = [item for item in candidates if item[3]]
    selected = max(successful if successful else candidates, key=lambda item: item[2])
    index, result, fitted_ll, success, _record = selected
    point = np.asarray(result.x, dtype=float)
    boundary_mask = point <= 1e-8
    boundary_rates = [name for name, flag in zip(free, boundary_mask) if flag]
    eigenvalues, _rank, _dimension, _positive, _negative, flat = _local_curvature(
        objective, point, boundary_mask
    )
    higher_failed = any(not item[3] and item[2] > fitted_ll + 1e-10 * max(1.0, abs(fitted_ll)) for item in candidates)
    message = str(result.message)
    if higher_failed:
        message += "; a higher-likelihood start did not converge"
    return {
        "rates": rates_at(point), "log_likelihood": fitted_ll,
        "optimizer_success": success and not higher_failed,
        "optimizer_status": int(result.status), "optimizer_message": message,
        "iterations": int(result.nit), "boundary_rates": boundary_rates,
        "locally_flat": bool(flat), "curvature_eigenvalues": list(eigenvalues),
        "selected_start": index, "starts": records,
    }


def fit_dna_family(tree, observations_by_site, *, root_frequency, root_presence, fixed_rates, workers=1):
    """Fit rates on the ascertainment-conditioned independent-site likelihood."""
    labels = tuple(sorted(tree.leaf_by_label))
    patterns = _compress_patterns(list(observations_by_site.values()), labels)
    return _fit_family(tree, patterns, root_frequency, root_presence, fixed_rates, workers)
