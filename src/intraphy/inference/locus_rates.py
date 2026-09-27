"""Nonnegative maximum-likelihood fitting for independent locus units."""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
import math
from typing import Mapping, Sequence

import numpy as np
from scipy.optimize import minimize

from intraphy.inference.locus_likelihood import evaluate_locus


@dataclass(frozen=True)
class LocusFitUnit:
    """One connected locus and its observations on a fixed tree."""
    process: object
    tree: object
    tips: Mapping[str, object]
    root_prior: object
    name: str = ""
    feature_sensitivity: Mapping[str, float] | None = None
    feature_specificity: Mapping[str, float] | None = None


@dataclass(frozen=True)
class LocusFitResult:
    initial_log_likelihood: float
    fitted_log_likelihood: float
    rates: Mapping[str, float]
    optimizer_success: bool
    optimizer_status: int
    optimizer_message: str
    iterations: int
    boundary_rates: tuple[str, ...]
    curvature_eigenvalues: tuple[float, ...]
    curvature_rank: int
    curvature_dimension: int
    curvature_positive_definite: bool
    curvature_has_negative_eigenvalue: bool
    locally_flat: bool
    workers: int
    start_results: tuple[Mapping, ...] = ()
    selected_start: int = 0
    curvature_evaluated: bool = True
    unresolved_higher_likelihood: bool = False


_WORKER_UNITS = None


def _install_units(units):
    global _WORKER_UNITS
    _WORKER_UNITS = units


def _unit_log_likelihood(task):
    index, rates = task
    unit = _WORKER_UNITS[index]
    return evaluate_locus(unit.process, unit.tree, unit.tips, unit.root_prior,
                          rates, posterior=False,
                          feature_sensitivity=unit.feature_sensitivity,
                          feature_specificity=unit.feature_specificity).log_likelihood


def _evaluate_units(units, rates, pool=None):
    if pool is None:
        values = [evaluate_locus(u.process, u.tree, u.tips, u.root_prior,
                                 rates, posterior=False,
                                 feature_sensitivity=u.feature_sensitivity,
                                 feature_specificity=u.feature_specificity).log_likelihood for u in units]
    else:
        values = list(pool.map(_unit_log_likelihood,
                               ((index, dict(rates)) for index in range(len(units)))))
    if any(value == -math.inf for value in values):
        return -math.inf
    total = math.fsum(values)
    if not math.isfinite(total):
        raise FloatingPointError("Joint log likelihood is non-finite")
    return total


def _local_curvature(objective, point, boundary):
    """Finite-difference observed Hessian; a local diagnostic only."""
    n = len(point)
    if n == 0:
        return (), 0, 0, True, False, False
    steps = np.maximum(1e-5, np.maximum(np.abs(point), 1.0) * 2e-4)
    one_sided = np.asarray(boundary, dtype=bool) | (point <= steps)
    base = objective(point)
    hessian = np.zeros((n, n), dtype=float)
    for i in range(n):
        h = steps[i]
        ei = np.zeros(n); ei[i] = h
        if one_sided[i]:
            f1 = objective(point + ei)
            f2 = objective(point + 2.0 * ei)
            hessian[i, i] = (f2 - 2.0 * f1 + base) / (h * h)
        else:
            hessian[i, i] = (objective(point + ei) - 2.0 * base
                             + objective(point - ei)) / (h * h)
        for j in range(i):
            hj = steps[j]
            ej = np.zeros(n); ej[j] = hj
            if one_sided[i] and one_sided[j]:
                value = (objective(point + ei + ej) - objective(point + ei)
                         - objective(point + ej) + base) / (h * hj)
            elif one_sided[i]:
                value = (objective(point + ei + ej) - objective(point + ei - ej)
                         - objective(point + ej) + objective(point - ej)) / (2 * h * hj)
            elif one_sided[j]:
                value = (objective(point + ei + ej) - objective(point - ei + ej)
                         - objective(point + ei) + objective(point - ei)) / (2 * h * hj)
            else:
                value = (objective(point + ei + ej) - objective(point + ei - ej)
                         - objective(point - ei + ej) + objective(point - ei - ej)) / (4 * h * hj)
            hessian[i, j] = hessian[j, i] = value
    if not np.all(np.isfinite(hessian)):
        return (), 0, n, False, False, True
    eigenvalues = np.linalg.eigvalsh((hessian + hessian.T) / 2.0)
    scale = max(float(np.max(np.abs(eigenvalues))), 1.0)
    rank = int(np.count_nonzero(np.abs(eigenvalues) > scale * 1e-7))
    positive = bool(np.all(eigenvalues > scale * 1e-7))
    negative = bool(np.any(eigenvalues < -scale * 1e-7))
    return tuple(float(x) for x in eigenvalues), rank, n, positive, negative, rank < n


def fit_locus_rates(units, initial_rates: Mapping[str, float], fixed_rates=None,
                    *, maxiter=1000, workers=1, ftol=1e-9, gtol=1e-6,
                    start_scales=(1.0,), compute_curvature=True,
                    _explicit_start_rates: Sequence[Mapping[str, float]] = ()):
    """Fit shared rate groups by ML over conditionally independent units.

    Equal edge rate_group labels tie rates across all units. Independent
    likelihood contributions are summed; units must represent independent
    observations under the declared model. No global scale parameter or
    arbitrary upper rate bound is introduced.
    """
    units = tuple(units)
    start_scales = tuple(start_scales)
    _explicit_start_rates = tuple(_explicit_start_rates)
    if not units:
        raise ValueError("At least one locus unit is required")
    if type(workers) is not int or workers < 1:
        raise ValueError("workers must be a positive integer")
    effective_workers = min(workers, len(units))
    if type(maxiter) is not int or maxiter <= 0:
        raise ValueError("maxiter must be a positive integer")
    if type(compute_curvature) is not bool:
        raise ValueError("compute_curvature must be Boolean")
    if not isinstance(start_scales, (tuple, list)) or not start_scales:
        raise ValueError("start_scales must be a nonempty sequence of finite positive numbers")
    checked_scales = []
    for value in start_scales:
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, float, np.number)):
            raise ValueError("start_scales must contain finite positive numbers, not booleans")
        value = float(value)
        if not math.isfinite(value) or value <= 0:
            raise ValueError("start_scales must contain finite positive numbers")
        checked_scales.append(value)
    fixed = dict(fixed_rates or {})
    initial = dict(initial_rates)
    if set(fixed) & set(initial):
        raise ValueError("A rate group cannot be both fixed and fitted")
    required = {event.rate_group for unit in units
                for event in unit.process.catalogue.opportunities}
    required.update(edge.rate_group for unit in units for edge in unit.process.edges)
    if set(fixed) | set(initial) != required:
        missing = sorted(required - (set(fixed) | set(initial)))
        extra = sorted((set(fixed) | set(initial)) - required)
        raise ValueError(f"Rate groups mismatch; missing={missing}, extra={extra}")
    for group, value in {**fixed, **initial}.items():
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, float, np.number)) or not math.isfinite(float(value)) or float(value) < 0:
            raise ValueError(f"Rate {group!r} must be finite and nonnegative")
    free_groups = tuple(sorted(initial))
    if any(float(initial[group]) <= 0 for group in free_groups):
        raise ValueError("Fitted rates require positive explicit starting values")
    scales = np.array([float(initial[group]) for group in free_groups], dtype=float)
    pool = None
    try:
        if effective_workers > 1:
            pool = ProcessPoolExecutor(max_workers=effective_workers,
                                       initializer=_install_units, initargs=(units,))
        def rates_at(point):
            rates = dict(fixed)
            rates.update(zip(free_groups, (float(scale * x) for scale, x in zip(scales, point))))
            return rates

        def log_likelihood(point):
            return _evaluate_units(units, rates_at(point), pool)

        x0 = np.ones(len(free_groups), dtype=float)
        initial_ll = log_likelihood(x0)
        if initial_ll == -math.inf:
            return LocusFitResult(initial_ll, initial_ll, rates_at(x0), False, -1,
                                  "observations have zero likelihood under the declared reachable process",
                                  0, (), (), 0, len(free_groups), False, False, False,
                                  effective_workers, (), 0, False)
        if not free_groups:
            return LocusFitResult(initial_ll, initial_ll, dict(fixed), True, 0,
                                  "all rates fixed", 0, (), (), 0, 0, compute_curvature, False, False,
                                  effective_workers, (), 0, compute_curvature)
        def objective(point):
            value = log_likelihood(point)
            return 1e300 if value == -math.inf else -value
        starts = [(scale, np.full(len(free_groups), scale), "scale") for scale in checked_scales]
        for explicit in _explicit_start_rates:
            values = dict(explicit)
            if set(values) != set(free_groups):
                raise ValueError("Each explicit start must specify exactly the free rate groups")
            point = []
            for group in free_groups:
                value = values[group]
                if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, float, np.number)):
                    raise ValueError("Explicit starting rates must be finite positive numbers")
                value = float(value)
                if not math.isfinite(value) or value <= 0:
                    raise ValueError("Explicit starting rates must be finite positive numbers")
                point.append(value / float(initial[group]))
            starts.append((None, np.asarray(point, dtype=float), "explicit_rates"))

        start_records = []
        candidates = []
        for start_index, (start_scale, start_point, start_source) in enumerate(starts):
            start_ll = log_likelihood(start_point)
            result = minimize(objective, start_point, method="L-BFGS-B",
                              bounds=[(0.0, None)] * len(free_groups),
                              options={"maxiter": maxiter, "ftol": ftol, "gtol": gtol})
            fitted_ll = log_likelihood(result.x)
            tolerance = 1e-10 * max(1.0, abs(start_ll)) if math.isfinite(start_ll) else 0.0
            valid = bool(result.success and math.isfinite(fitted_ll) and fitted_ll >= start_ll - tolerance)
            fitted_rates = rates_at(result.x)
            record = {"start_index": start_index, "start_source": start_source,
                      "start_scale": start_scale, "start_rates": rates_at(start_point),
                      "start_log_likelihood": start_ll, "optimizer_success": valid,
                      "optimizer_status": int(result.status), "optimizer_message": str(result.message),
                      "iterations": int(result.nit), "fitted_log_likelihood": fitted_ll,
                      "fitted_rates": fitted_rates}
            start_records.append(record)
            if math.isfinite(fitted_ll):
                candidates.append((start_index, result, fitted_rates, fitted_ll, valid, record))

        successful = [item for item in candidates if item[4]]
        selected = max(successful or candidates, key=lambda item: item[3]) if candidates else None
        if selected is None:
            return LocusFitResult(initial_ll, -math.inf, rates_at(x0), False, -1,
                                  "all optimizer starts returned non-finite likelihoods", 0, (), (),
                                  0, len(free_groups), False, False, False, effective_workers,
                                  tuple(start_records), 0, False)
        selected_index, result, fitted_rates, fitted_ll, selected_success, selected_record = selected
        boundary_mask = np.asarray(result.x <= 1e-8 * np.maximum(1.0, np.abs(np.ones(len(free_groups)))))
        boundary = tuple(group for group, flag in zip(free_groups, boundary_mask) if flag)
        if compute_curvature:
            eigenvalues, rank, dimension, positive, negative, flat = _local_curvature(
                objective, np.asarray(result.x, dtype=float), boundary_mask)
        else:
            eigenvalues, rank, dimension, positive, negative, flat = (), 0, len(free_groups), False, False, False
        message = str(result.message)
        higher_failed = [item for item in candidates if not item[4] and item[3] > fitted_ll + 1e-10 * max(1.0, abs(fitted_ll))]
        if higher_failed:
            message += "; a higher-likelihood start did not converge; selected result is not claimed as global optimum"
        return LocusFitResult(initial_ll, fitted_ll, fitted_rates, selected_success,
                              int(result.status), message, int(result.nit), boundary,
                              eigenvalues, rank, dimension, positive, negative, flat,
                              effective_workers, tuple(start_records), selected_index, compute_curvature,
                              bool(higher_failed))
    finally:
        if pool is not None:
            pool.shutdown(wait=True, cancel_futures=True)
