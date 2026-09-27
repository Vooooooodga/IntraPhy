"""Fixed-zero nested comparisons for identical evidence-conditioned loci."""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Mapping

from .locus_rates import LocusFitResult, fit_locus_rates


@dataclass(frozen=True)
class LocusRateComparison:
    null_fit: LocusFitResult
    full_fit: LocusFitResult
    statistic: float | None
    status: str
    raw_log_likelihood_difference: float | None
    likelihood_tolerance: float
    interpretation_limit: str


def compare_locus_rates(units, initial_rates: Mapping[str, float], fixed_rates: Mapping[str, float],
                        null_rate_groups, *, start_scales=(0.2, 1.0, 5.0),
                        workers=1, maxiter=1000, compute_curvature=True) -> LocusRateComparison:
    """Compare full MLE with a same-data nested model fixing selected rates at 0.

    No P value or chi-square calibration is returned. Likelihood differences
    within the declared numerical tolerance are represented by statistic 0.0
    in either direction while the raw difference is retained. The comparison is
    conditional on the caller supplying identical units, tree, tips, root prior,
    state space, and finite opportunity catalogue to both fits.
    """
    units = tuple(units)
    start_scales = tuple(start_scales)
    initial = dict(initial_rates)
    fixed = dict(fixed_rates)
    null_groups = tuple(null_rate_groups)
    if not null_groups or len(set(null_groups)) != len(null_groups):
        raise ValueError("null_rate_groups must contain unique free-rate group names")
    if any(group not in initial or group in fixed for group in null_groups):
        raise ValueError("Every null rate group must be free in initial_rates")

    null_initial = {key: value for key, value in initial.items() if key not in null_groups}
    null_fixed = dict(fixed)
    null_fixed.update({group: 0.0 for group in null_groups})
    null_fit = fit_locus_rates(units, null_initial, null_fixed, start_scales=start_scales,
                               workers=workers, maxiter=maxiter, compute_curvature=compute_curvature)

    # Use the null nuisance estimates as an additional full-model start. Keep
    # all tested rates positive so they remain valid free-parameter starts.
    smallest_scale = min(float(value) for value in start_scales)
    null_seed = {}
    for group, initial_value in initial.items():
        if group in null_groups:
            null_seed[group] = float(initial_value) * smallest_scale
        else:
            value = float(null_fit.rates.get(group, initial_value))
            null_seed[group] = value if math.isfinite(value) and value > 0 else float(initial_value) * smallest_scale
    full_fit = fit_locus_rates(units, initial, fixed, start_scales=start_scales,
                               workers=workers, maxiter=maxiter, compute_curvature=compute_curvature,
                               _explicit_start_rates=(null_seed,))
    limit = ("Conditional likelihood comparison for the same supplied locus, tree, root prior, observations, "
             "state space, and finite opportunity catalogue. No default P value or chi-square calibration is supplied.")
    finite_fits = (math.isfinite(null_fit.fitted_log_likelihood)
                   and math.isfinite(full_fit.fitted_log_likelihood))
    raw = (full_fit.fitted_log_likelihood - null_fit.fitted_log_likelihood) if finite_fits else None
    if (not null_fit.optimizer_success or not full_fit.optimizer_success
            or null_fit.unresolved_higher_likelihood or full_fit.unresolved_higher_likelihood):
        return LocusRateComparison(null_fit, full_fit, None, "one_or_both_fits_unsuccessful",
                                   raw, 0.0, limit)
    tolerance = 1e-10 * max(1.0, abs(null_fit.fitted_log_likelihood), abs(full_fit.fitted_log_likelihood))
    if raw < -tolerance:
        return LocusRateComparison(null_fit, full_fit, None, "full_fit_below_null_beyond_tolerance",
                                   raw, tolerance, limit)
    statistic = 0.0 if abs(raw) <= tolerance else 2.0 * raw
    return LocusRateComparison(null_fit, full_fit, statistic,
                               "evaluated_within_numerical_tolerance", raw, tolerance, limit)
