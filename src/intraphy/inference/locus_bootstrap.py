"""Fixed-catalogue parametric bootstrap for DNA-only nested rate comparisons."""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
import math
from typing import Mapping

import numpy as np

from ..structure.locus_observations import observation_emission
from ..structure.locus_sampling import sample_history, sample_tip_observation
from .locus_comparison import compare_locus_rates
from .locus_rates import LocusFitUnit


@dataclass(frozen=True)
class LocusBootstrapResult:
    observed_comparison: object
    null_rates: Mapping[str, float]
    replicates: tuple[Mapping, ...]
    p_value: float | None
    tail_count: int
    resolution: float
    requested_replicates: int
    completed_replicates: int
    seed: int
    workers: int
    sampling_design: str
    status: str
    limitation: str = "Plug-in parametric-bootstrap calibration conditional on the fixed catalogue, tree, root, and detection/missingness design; not an exact test for composite nulls."


def validate_dna_bootstrap_design(bundle):
    """Validate fixed-catalogue DNA-only inputs without fitting or sampling.

    Latent catalogue splice features are permitted. Tip feature/path
    observations and feature-level detection metadata are not.
    """
    raw_units = bundle.model_record.get("units", ())
    if len(raw_units) != len(bundle.units):
        raise ValueError("model record and compiled locus units differ in length")
    empty_values = {
        "features": {}, "surveyed_features": [], "sensitivity": {},
        "specificity": {}, "observed_paths": [],
    }
    for unit, raw in zip(bundle.units, raw_units):
        if (raw.get("family"), raw.get("unit")) != (unit.family, unit.unit):
            raise ValueError("raw and compiled family/unit ordering does not match")
        raw_material_ids = tuple(entry.get("id") for entry in
                                 raw.get("catalogue", {}).get("material", ()))
        compiled_material_ids = tuple(tract.id for tract in unit.process.catalogue.material)
        if raw_material_ids != compiled_material_ids:
            raise ValueError(f"raw and compiled material ID order differs for {unit.family}/{unit.unit}")
        observations = raw["observations"]["tips"]
        if set(observations) != set(unit.tree.leaf_by_label):
            raise ValueError("raw observation tips do not match the fixed tree")
        for label, tip in observations.items():
            for field, empty in empty_values.items():
                if tip.get(field, empty) != empty:
                    raise ValueError(f"DNA-only bootstrap requires empty tip field {field!r} ({label})")


def _fit_rates(bundle):
    initial = {group: float(record["value"])
               for group, record in bundle.rates.items() if record["mode"] == "fit"}
    fixed = {group: float(record["value"])
             for group, record in bundle.rates.items() if record["mode"] == "fixed"}
    return initial, fixed


def _simulated_units(bundle, rates, rng):
    sampled_units = []
    for unit, raw in zip(bundle.units, bundle.model_record["units"]):
        latent = sample_history(unit.process, unit.tree, unit.root_prior, rates, rng)
        tips = {}
        for label, node in unit.tree.leaf_by_label.items():
            record = raw["observations"]["tips"][label]
            observation = sample_tip_observation(unit.process.catalogue, record, latent[node], rng)
            tips[label] = observation_emission(
                unit.process.catalogue, unit.process.states, observation,
                material_sensitivity=record.get("material_sensitivity") or None,
                material_specificity=record.get("material_specificity") or None,
            )
        sampled_units.append(LocusFitUnit(unit.process, unit.tree, tips, unit.root_prior,
                                          name=f"{unit.family}/{unit.unit}"))
    return tuple(sampled_units)


def _usable(comparison):
    statistic = comparison.statistic
    return (comparison.null_fit.optimizer_success and comparison.full_fit.optimizer_success
            and statistic is not None and math.isfinite(float(statistic)))


_WORKER_CONTEXT = None


def _init_replicate_worker(context):
    global _WORKER_CONTEXT
    _WORKER_CONTEXT = context


def _run_replicate(args):
    replicate, context = args
    bundle, null_rates, initial, fixed, null_rate_groups, seed, start_scales, maxiter = context
    row = {"replicate": replicate, "status": "failed", "comparison_status": "",
           "statistic": None, "null_log_likelihood": None, "full_log_likelihood": None,
           "error": ""}
    try:
        rng = np.random.default_rng(np.random.SeedSequence([seed, replicate]))
        units = _simulated_units(bundle, null_rates, rng)
        result = compare_locus_rates(units, initial, fixed, null_rate_groups,
                                     start_scales=start_scales, workers=1,
                                     maxiter=maxiter, compute_curvature=False)
        row["comparison_status"] = result.status
        row["null_log_likelihood"] = float(result.null_fit.fitted_log_likelihood)
        row["full_log_likelihood"] = float(result.full_fit.fitted_log_likelihood)
        if not _usable(result):
            row["error"] = "null/full fit unsuccessful or statistic non-finite"
        else:
            row["status"] = "success"
            row["statistic"] = float(result.statistic)
    except Exception as exc:
        row["error"] = f"{type(exc).__name__}: {exc}"
    return row


def _run_replicate_from_worker(args):
    replicate, _unused = args
    return _run_replicate((replicate, _WORKER_CONTEXT))


def bootstrap_locus_comparison(bundle, null_rate_groups, replicates, seed, *,
                               sampling_design, start_scales=(0.2, 1.0, 5.0),
                               workers=1, maxiter=1000):
    """Calibrate a nested DNA-rate statistic under a declared fixed catalogue.

    ``sampling_design`` must explicitly be ``"fixed_catalogue"``. Missing masks,
    detection parameters, tree, root prior, and opportunity catalogue stay fixed;
    no ascertainment correction or catalogue resampling is attempted. The result's
    ``workers`` records the effective process-worker count; it is zero when the
    observed comparison fails and no replicate is attempted.
    """
    if sampling_design != "fixed_catalogue":
        raise ValueError("sampling_design must explicitly equal 'fixed_catalogue'")
    for value, name in ((replicates, "replicates"), (workers, "workers"), (maxiter, "maxiter")):
        if type(value) is not int or value < 1:
            raise ValueError(f"{name} must be a positive integer")
    if type(seed) is not int or seed < 0:
        raise ValueError("seed must be a nonnegative integer")
    null_rate_groups = tuple(null_rate_groups)
    start_scales = tuple(start_scales)
    if not start_scales or any(not math.isfinite(float(scale)) or float(scale) <= 0 for scale in start_scales):
        raise ValueError("start_scales must contain finite positive values")
    validate_dna_bootstrap_design(bundle)
    initial, fixed = _fit_rates(bundle)
    observed = compare_locus_rates(bundle.fit_units(), initial, fixed, null_rate_groups,
                                   start_scales=start_scales, workers=1,
                                   maxiter=maxiter, compute_curvature=False)
    if not _usable(observed):
        return LocusBootstrapResult(observed, {}, (), None, 0, 1.0 / (replicates + 1),
            replicates, 0, seed, 0, sampling_design, "observed_comparison_failed")
    null_rates = {key: float(value) for key, value in observed.null_fit.rates.items()}
    if set(null_rates) != set(bundle.rates):
        raise ValueError("null fitted rates do not exactly cover process rate groups")
    observed_statistic = float(observed.statistic)

    count_workers = min(workers, replicates)
    context = (bundle, null_rates, initial, fixed, null_rate_groups, seed, start_scales, maxiter)
    if count_workers == 1:
        rows = [_run_replicate((i, context)) for i in range(replicates)]
    else:
        with ProcessPoolExecutor(max_workers=count_workers,
                                 initializer=_init_replicate_worker,
                                 initargs=(context,)) as pool:
            rows = list(pool.map(_run_replicate_from_worker,
                                 ((i, None) for i in range(replicates))))
    successful = [row for row in rows if row["status"] == "success"]
    tail = sum(float(row["statistic"]) >= observed_statistic for row in successful)
    complete = len(successful) == replicates
    p_value = (1 + tail) / (replicates + 1) if complete else None
    return LocusBootstrapResult(observed, null_rates, tuple(rows), p_value, tail,
        1.0 / (replicates + 1), replicates, len(successful), seed, count_workers,
        sampling_design, "completed" if complete else "incomplete")
