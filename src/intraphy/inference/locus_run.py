"""Run and serialize finite evidence-conditioned genomic copy histories."""
from __future__ import annotations

from pathlib import Path
import numpy as np

from ..run_result import RunResult
from ..structure.serialization import json_safe, write_json
from ..storage.tabular import write_tsv
from .locus_likelihood import evaluate_locus
from .locus_rates import fit_locus_rates


def _state_record(state, catalogue):
    return {
        "material": {tract.id: state.material[i] for i, tract in enumerate(catalogue.material)},
    }


def _node_marginals(process, posterior):
    catalogue = process.catalogue
    rows = []
    for node_id, probabilities in posterior.items():
        values = np.asarray(probabilities, dtype=float)
        material = {}
        for i, tract in enumerate(catalogue.material):
            material[tract.id] = {str(status): float(sum(values[j] for j, state in enumerate(process.states)
                                                          if state.material[i] == status))
                                  for status in catalogue.material_states}
        copies = {copy.id: float(sum(values[j] for j, state in enumerate(process.states)
                                    if all(state.material[next(k for k, tract in enumerate(catalogue.material)
                                                               if tract.id == material_id)] == 1
                                           for material_id in copy.material_ids)))
                  for copy in catalogue.copies}
        rows.append({"node_id": node_id, "material_state": material,
                     "copy_intact": copies})
    return rows


def _fit_record(result, mode, bundle):
    estimated_groups = tuple(sorted(key for key, value in bundle.rates.items()
                                    if mode == "fit" and value["mode"] == "fit"))
    if not estimated_groups:
        parameter_status = "fixed_parameter_evaluation"
    elif not result.optimizer_success:
        parameter_status = "optimizer_unsuccessful"
    elif result.boundary_rates:
        parameter_status = "boundary_rate_estimate"
    elif result.locally_flat:
        parameter_status = "locally_flat_objective"
    elif result.curvature_has_negative_eigenvalue:
        parameter_status = "negative_local_curvature"
    elif not result.curvature_positive_definite:
        parameter_status = "curvature_not_positive_definite"
    else:
        parameter_status = "converged_local_solution"
    return {
        "schema": "intraphy.exon-locus-fit/2", "model": "exon-locus-ctmc",
        "parameter_mode": mode, "log_likelihood_initial": result.initial_log_likelihood,
        "log_likelihood_fitted": result.fitted_log_likelihood,
        "optimizer_success": result.optimizer_success, "optimizer_status": result.optimizer_status,
        "optimizer_message": result.optimizer_message, "iterations": result.iterations,
        "rates": dict(result.rates), "boundary_rates": list(result.boundary_rates),
        "curvature_eigenvalues": list(result.curvature_eigenvalues),
        "curvature_rank": result.curvature_rank, "curvature_dimension": result.curvature_dimension,
        "curvature_positive_definite": result.curvature_positive_definite,
        "curvature_has_negative_eigenvalue": result.curvature_has_negative_eigenvalue,
        "locally_flat": result.locally_flat, "workers": result.workers,
        "estimated_rate_groups": list(estimated_groups),
        "initial_rates": {key: bundle.rates[key]["value"] for key in estimated_groups},
        "curvature_coordinates": "rate divided by supplied initial rate, in sorted estimated_rate_groups order",
        "rate_provenance": bundle.rate_provenance,
        "branch_length_unit": bundle.branch_length_unit,
        "conditional_independence_provenance": bundle.independence_provenance,
        "state_models": {f"{unit.family}/{unit.unit}": unit.process.catalogue.state_model
                         for unit in bundle.units},
        "parameter_status": parameter_status,
        "scope": "genomic DNA material and copy histories conditional on labelled positions, supplied finite opportunities, root distributions, fixed tree, and declared independent units",
        "model_specification": bundle.model_record,
        "tree": list(bundle.tree_rows),
        "tree_provenance": bundle.tree_provenance,
    }


def analyze_locus(bundle, output_dir, *, parameter_mode="fit", expected_counts=False, workers=1):
    """Fit or evaluate declared global rate groups and write posterior summaries."""
    if parameter_mode not in {"fit", "fixed"}:
        raise ValueError("parameter_mode must be 'fit' or 'fixed'")
    if type(expected_counts) is not bool:
        raise TypeError("expected_counts must be Boolean")
    target = Path(output_dir)
    target.mkdir(parents=True, exist_ok=True)
    history_path = target / "locus_history.json"
    fit_path = target / "locus_fit.json"
    tree_path = target / "species_tree.tsv"
    result_manifest = target / "run_result.json"
    for path in (history_path, fit_path, tree_path, result_manifest):
        if path.exists() or path.is_symlink() or path.with_name(path.name + ".tmp").exists():
            raise ValueError(f"Refusing existing locus result: {path}")
    write_tsv(tree_path, bundle.tree_rows, ["node_id", "parent_id", "label", "branch_length"])

    initial = {key: record["value"] for key, record in bundle.rates.items()
               if parameter_mode == "fit" and record["mode"] == "fit"}
    fixed = {key: record["value"] for key, record in bundle.rates.items()
             if parameter_mode == "fixed" or record["mode"] == "fixed"}
    fit = fit_locus_rates(bundle.fit_units(), initial, fixed_rates=fixed,
                          workers=min(workers, len(bundle.units)))
    write_json(fit_path, json_safe(_fit_record(fit, parameter_mode, bundle)))
    if not fit.optimizer_success or not np.isfinite(fit.fitted_log_likelihood):
        RunResult("exon-locus-ctmc", "single-locus-copy-slots", tree_path.name,
                  (fit_path.name, tree_path.name), status="failed").write(target)
        raise RuntimeError("The declared exon-locus likelihood could not be fitted/evaluated; see locus_fit.json.")
    units = []
    for unit in bundle.units:
        result = evaluate_locus(unit.process, unit.tree, unit.tips, unit.root_prior,
                                fit.rates, posterior=True, counts=expected_counts)
        catalogue = unit.process.catalogue
        branches = []
        if result.branch_event_counts is not None:
            for child, counts in result.branch_event_counts.items():
                branches.append({"child_node_id": child, "parent_node_id": unit.tree.parent[child],
                                 "branch_length": unit.tree.branch_length(child),
                                 "expected_event_counts": dict(counts)})
        units.append({
            "family": unit.family, "unit": unit.unit,
            "state_model": catalogue.state_model,
            "material_state_semantics": ({"0": "absent", "1": "present"}
                                         if catalogue.state_model == "binary" else
                                         {"0": "unintroduced", "1": "present", "2": "deleted"}),
            "state_count": len(unit.process.states),
            "states": [_state_record(state, catalogue) for state in unit.process.states],
            "root_distribution": {str(i): float(value) for i, value in enumerate(unit.root_prior) if value > 0},
            "root_provenance": unit.root_provenance,
            "observation_provenance": unit.observation_provenance,
            "log_likelihood": result.log_likelihood,
            "node_marginals": _node_marginals(unit.process, result.node_posteriors),
            "branch_event_counts": branches,
            "event_opportunities": [{"id": item.id, "outcome_id": item.outcome_id,
                                     "kind": item.kind, "rate_group": item.rate_group,
                                     "weight": item.weight}
                                    for item in catalogue.opportunities],
            "catalogue_provenance": catalogue.provenance,
        })
    history = {
        "schema": "intraphy.exon-locus-history/2", "model": "exon-locus-ctmc",
        "status": "completed", "parameter_mode": parameter_mode,
        "parameter_status": _fit_record(fit, parameter_mode, bundle)["parameter_status"],
        "log_likelihood": fit.fitted_log_likelihood,
        "rates": dict(fit.rates), "branch_length_unit": bundle.branch_length_unit,
        "tree_provenance": bundle.tree_provenance,
        "model_provenance": bundle.provenance,
        "rate_provenance": bundle.rate_provenance,
        "conditional_independence_provenance": bundle.independence_provenance,
        "model_specification_ref": fit_path.name,
        "finite_supplied_opportunity_catalogue_only": True,
        "units": units,
    }
    write_json(history_path, json_safe(history))
    RunResult("exon-locus-ctmc", "single-locus-copy-slots", tree_path.name,
              (history_path.name, fit_path.name, tree_path.name)).write(target)
