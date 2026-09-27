"""Run conditional nested-rate statistics on a preflight-loaded locus bundle."""
from __future__ import annotations

from dataclasses import asdict
import math
from pathlib import Path

from ..run_result import RunResult
from ..storage.tabular import write_tsv
from ..structure.serialization import json_safe, write_json
from .locus_comparison import compare_locus_rates
from .locus_profiles import profile_locus_rate


def _rates(bundle):
    initial = {key: float(value["value"]) for key, value in bundle.rates.items()
               if value["mode"] == "fit"}
    fixed = {key: float(value["value"]) for key, value in bundle.rates.items()
             if value["mode"] == "fixed"}
    return initial, fixed


def _usable_comparison(result):
    return (result.status == "evaluated_within_numerical_tolerance"
            and result.statistic is not None and math.isfinite(float(result.statistic))
            and result.null_fit.optimizer_success and result.full_fit.optimizer_success
            and not result.null_fit.unresolved_higher_likelihood
            and not result.full_fit.unresolved_higher_likelihood)


def run_locus_statistics(bundle, output_dir, *, null_rate_groups,
                         start_scales=(0.2, 1.0, 5.0), profile_rate_group=None,
                         profile_values=None, bootstrap_replicates=0, seed=None,
                         sampling_design=None, workers=1, maxiter=1000):
    """Compare rate groups, optionally profiling and bootstrapping, then serialize."""
    null_rate_groups = tuple(null_rate_groups)
    start_scales = tuple(start_scales)
    profile_values = None if profile_values is None else tuple(profile_values)
    target = Path(output_dir)
    target.mkdir(parents=True, exist_ok=True)
    paths = (target / "locus_statistics.json", target / "species_tree.tsv",
             target / "run_result.json")
    if any(path.exists() or path.is_symlink() for path in paths):
        raise ValueError("Refusing to overwrite existing locus-statistics output")

    initial, fixed = _rates(bundle)
    bootstrap = None
    if bootstrap_replicates:
        from .locus_bootstrap import bootstrap_locus_comparison

        bootstrap = bootstrap_locus_comparison(
            bundle, null_rate_groups, bootstrap_replicates, seed,
            sampling_design=sampling_design, start_scales=start_scales,
            workers=workers, maxiter=maxiter)
        comparison = bootstrap.observed_comparison
    else:
        comparison = compare_locus_rates(
            bundle.fit_units(), initial, fixed, null_rate_groups,
            start_scales=start_scales, workers=workers, maxiter=maxiter)

    profile = None
    if profile_rate_group is not None:
        profile = profile_locus_rate(
            bundle.fit_units(), initial, fixed, profile_rate_group, profile_values,
            start_scales=start_scales, workers=workers, maxiter=maxiter)

    profile_ok = (profile is None or
                  (profile.status == "evaluated"
                   and profile.reference_fit.optimizer_success
                   and not profile.reference_fit.unresolved_higher_likelihood
                   and all(point.status == "zero_likelihood"
                           or (point.status == "evaluated" and point.fit.optimizer_success
                               and not point.fit.unresolved_higher_likelihood)
                           for point in profile.points)))
    bootstrap_ok = bootstrap is None or bootstrap.status == "completed"
    status_ok = _usable_comparison(comparison) and profile_ok and bootstrap_ok
    failure_reasons = []
    if not _usable_comparison(comparison):
        failure_reasons.append(f"comparison:{comparison.status}")
    if not profile_ok:
        failure_reasons.append(f"profile:{profile.status}")
    if not bootstrap_ok:
        failure_reasons.append(f"bootstrap:{bootstrap.status}")

    payload = {
        "schema": "intraphy.locus-statistics/2",
        "status": "completed" if status_ok else "failed",
        "model": "exon-locus-ctmc",
        "analysis_scope": "conditional rate comparison for supplied genomic DNA copy units",
        "state_models": {f"{unit.family}/{unit.unit}": unit.process.catalogue.state_model
                         for unit in bundle.units},
        "material_state_semantics": {
            f"{unit.family}/{unit.unit}": ({"0": "absent", "1": "present"}
                                         if unit.process.catalogue.state_model == "binary" else
                                         {"0": "unintroduced", "1": "present", "2": "deleted"})
            for unit in bundle.units},
        "model_record": bundle.model_record,
        "tree": list(bundle.tree_rows),
        "tree_provenance": bundle.tree_provenance,
        "branch_length_unit": bundle.branch_length_unit,
        "model_provenance": bundle.provenance,
        "observation_provenance": [unit.observation_provenance for unit in bundle.units],
        "settings": {"null_rate_groups": list(null_rate_groups),
                     "start_scales": list(start_scales),
                     "profile_rate_group": profile_rate_group,
                     "profile_values": profile_values,
                     "bootstrap_replicates": bootstrap_replicates,
                     "seed": seed, "sampling_design": sampling_design,
                     "workers": workers, "maxiter": maxiter},
        "comparison": asdict(comparison),
        "profile": None if profile is None else asdict(profile),
        "bootstrap": None if bootstrap is None else asdict(bootstrap),
        "failure_reasons": failure_reasons,
        "interpretation_limit": comparison.interpretation_limit,
    }
    write_json(paths[0], json_safe(payload))
    write_tsv(paths[1], bundle.tree_rows,
              ["node_id", "parent_id", "label", "branch_length"])
    RunResult("exon-locus-ctmc", "conditional-locus-statistics", paths[1].name,
              (paths[0].name, paths[1].name),
              status="completed" if status_ok else "failed").write(target)
    if not status_ok:
        raise RuntimeError("Locus statistics contain unsuccessful fit/profile/bootstrap diagnostics; see locus_statistics.json")
