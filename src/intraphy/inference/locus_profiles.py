"""Conditional fixed-rate profiles using the locus likelihood MLE API."""
from __future__ import annotations

from dataclasses import dataclass
import math
from numbers import Real
from typing import Mapping

from .locus_rates import LocusFitResult, fit_locus_rates


@dataclass(frozen=True)
class LocusProfilePoint:
    value: float
    fit: LocusFitResult
    log_likelihood_difference: float | None
    status: str


@dataclass(frozen=True)
class LocusProfileResult:
    rate_group: str
    reference_fit: LocusFitResult
    points: tuple[LocusProfilePoint, ...]
    reference_inadequate: bool
    status: str


def _nonnegative_finite(value, name):
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be finite and nonnegative, not Boolean")
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be finite and nonnegative") from error
    if not math.isfinite(number) or number < 0:
        raise ValueError(f"{name} must be finite and nonnegative")
    return number


def profile_locus_rate(units, initial_rates: Mapping[str, float], fixed_rates: Mapping[str, float],
                       rate_group: str, values, *, start_scales=(0.2, 1.0, 5.0),
                       workers=1, maxiter=1000) -> LocusProfileResult:
    """Profile one free rate while refitting all other free rate groups.

    The supplied values are evaluated as given; no upper bound, interval,
    confidence limit, or reference-distribution calibration is imposed. A
    fixed, zero-likelihood point with no nuisance rates is retained as
    ``zero_likelihood`` rather than treated as an optimization failure.
    """
    units = tuple(units)
    start_scales = tuple(start_scales)
    profile_values = tuple(_nonnegative_finite(value, "profile value") for value in values)
    initial = dict(initial_rates)
    fixed = dict(fixed_rates)
    if rate_group not in initial or rate_group in fixed:
        raise ValueError("rate_group must be declared as a free rate in initial_rates")
    if not profile_values:
        raise ValueError("values must contain at least one rate value")
    reference = fit_locus_rates(units, initial, fixed, start_scales=start_scales,
                                workers=workers, maxiter=maxiter, compute_curvature=False)
    nuisance_initial = {key: value for key, value in initial.items() if key != rate_group}
    points = []
    reference_inadequate = False
    for value in profile_values:
        point_fixed = dict(fixed)
        point_fixed[rate_group] = value
        fit = fit_locus_rates(units, nuisance_initial, point_fixed, start_scales=start_scales,
                              workers=workers, maxiter=maxiter, compute_curvature=False)
        if (not nuisance_initial and fit.fitted_log_likelihood == -math.inf
                and reference.optimizer_success and not reference.unresolved_higher_likelihood):
            difference, status = -math.inf, "zero_likelihood"
        elif not fit.optimizer_success or fit.unresolved_higher_likelihood:
            difference, status = None, "profile_fit_unsuccessful"
        elif (not reference.optimizer_success or reference.unresolved_higher_likelihood
              or not math.isfinite(reference.fitted_log_likelihood)):
            difference, status = None, "reference_fit_unsuccessful"
        else:
            difference = fit.fitted_log_likelihood - reference.fitted_log_likelihood
            tolerance = 1e-10 * max(1.0, abs(reference.fitted_log_likelihood))
            if difference > tolerance:
                reference_inadequate = True
                status = "profile_exceeds_reference"
            else:
                status = "evaluated"
        points.append(LocusProfilePoint(value, fit, difference, status))
    if reference_inadequate:
        status = "reference_inadequate"
    elif not reference.optimizer_success or reference.unresolved_higher_likelihood:
        status = "reference_fit_unsuccessful"
    elif any(point.status not in {"evaluated", "zero_likelihood"} for point in points):
        status = "partial_profile"
    else:
        status = "evaluated"
    return LocusProfileResult(rate_group, reference, tuple(points), reference_inadequate, status)
