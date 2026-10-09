"""Conditional simulation for the native genomic-exon foreground comparison."""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import lru_cache
import json
import logging
import math
from multiprocessing import get_context
from pathlib import Path
from time import perf_counter

import numpy as np

from ..structure.serialization import json_safe, write_json
from ..structure.space import enumerate_space
from ..topology import SpeciesTree
from .configuration_model import RateModel
from .exon_rates import InferenceUnit
from .exon_resampling import sample_unit_history
from .genomic_exon_calibration import _catalogues
from .genomic_exon_comparison import compare_foreground
from ..verification.exon_comparison_scoring import (
    independent_truth_pairs,
    native_endpoint_scores,
    paired_panel_summary,
    panel_summary,
    summarize_calibration,
)
from ..structure.edits import EDIT_KINDS


def balanced_tree(taxa):
    """Return a balanced power-of-two tree with unit root-to-tip height."""
    if type(taxa) is not int or taxa < 4 or taxa & (taxa - 1):
        raise ValueError("taxa must be a power of two and at least four")
    labels = tuple(f"taxon_{index:03d}" for index in range(1, taxa + 1))
    depth = int(math.log2(taxa))
    edge_length = 1.0 / depth
    rows = [{"node_id": "root", "parent_id": "", "label": "root",
             "branch_length": 0.0}]

    def descend(group, parent, level, serial):
        if len(group) == 1:
            return
        midpoint = len(group) // 2
        for side, child_group in enumerate((group[:midpoint], group[midpoint:])):
            child = (child_group[0] if len(child_group) == 1 else
                     f"node_{level + 1}_{serial * 2 + side}")
            rows.append({"node_id": child, "parent_id": parent, "label": child,
                         "branch_length": edge_length})
            descend(child_group, child, level + 1, serial * 2 + side)

    descend(labels, "root", 0, 0)
    return tuple(rows), labels


def foreground_clade(tree_rows, labels):
    """Select the stem and descendant edges subtending the first taxon quarter."""
    target = set(labels[:len(labels) // 4])
    parent = {row["node_id"]: row["parent_id"] for row in tree_rows}
    children = {}
    for row in tree_rows:
        if row["parent_id"]:
            children.setdefault(row["parent_id"], []).append(row["node_id"])
    descendants = {}

    def collect(node):
        if node not in children:
            descendants[node] = {node}
        else:
            descendants[node] = set().union(*(collect(child) for child in children[node]))
        return descendants[node]

    collect("root")
    roots = [node for node, values in descendants.items()
             if node != "root" and values == target]
    if len(roots) != 1:
        raise ValueError("Could not resolve the prespecified foreground clade")
    clade = roots[0]
    selected = set()
    for node in parent:
        if node == "root":
            continue
        current = node
        while current and current != "root":
            if current == clade:
                selected.add(node)
                break
            current = parent[current]
    return frozenset(selected), clade


@lru_cache(maxsize=2)
def _catalogue_spaces(scenario):
    catalogues = _catalogues(scenario)
    spaces = tuple(enumerate_space(catalogue) for catalogue in catalogues)
    if any(not space.complete for space in spaces):
        raise ValueError("Calibration catalogue must enumerate completely")
    return catalogues, spaces


def _unit_for_gene(catalogue, space, tree, labels, gene_id):
    reference = catalogue.observations[0].configurations[0]
    visible = (tuple((span.start, span.end) for span in reference.exons),
               tuple(int(value == 1) for value in reference.material))
    tips = {label: np.asarray([
        float(_visible_state(state) == visible)
        for state in space.states], dtype=float) for label in labels}
    return InferenceUnit(gene_id, catalogue.unit, space, tree, tips)


def _unit_record(unit):
    return {"family_id": unit.family, "unit_id": unit.unit,
            "space": unit.space, "tree": unit.tree, "tips": unit.tips}


def _mask_unit(unit, missing_tip):
    tips = {label: np.asarray(values, dtype=float).copy()
            for label, values in unit.tips.items()}
    if missing_tip is not None:
        tips[missing_tip] = np.ones(len(unit.space.states), dtype=float)
    from dataclasses import replace
    return replace(unit, tips=tips)


def _validate_tip_support(unit):
    """Accept exact calls and whole-tip missingness only."""
    for label, values in unit.tips.items():
        vector = np.asarray(values, dtype=float)
        all_unknown = np.all(vector == 1.)
        exact = False
        if np.all((vector == 0.) | (vector == 1.)) and np.any(vector > 0):
            signatures = {_visible_state(unit.space.states[index])
                          for index in np.flatnonzero(vector > 0)}
            exact = len(signatures) == 1
        if not (all_unknown or exact):
            raise ValueError(f"Partial candidate observation is unsupported for tip {label}")


def _visible_state(state):
    return (tuple((span.start, span.end) for span in state.exons),
            tuple(int(value == 1) for value in state.material))


def _mean(values):
    return float(np.mean(values)) if values else None


def _profile_record(comparison, multiplier):
    return next((point for point in comparison.get("profile", ())
                 if point.get("foreground_multiplier") == multiplier), None)


def _scoreable_native_fit(gene_fit):
    mu = gene_fit.get("mu")
    likelihood = gene_fit.get("log_likelihood")
    return (gene_fit.get("status") in {
                "estimated_conditional_composite_rate", "zero_boundary"}
            and mu is not None and math.isfinite(float(mu)) and float(mu) >= 0
            and likelihood is not None and math.isfinite(float(likelihood)))


def _replicate(task):
    replicate_started = perf_counter()
    replicate, scenario, rho, mask, taxa, gene_rates, seed = task
    row = {"replicate": replicate, "status": "recorded", "simulation_status": "incomplete",
           "generated_rho": rho, "observation_mask": mask, "genes": [],
           "foreground_comparison": {"status": "not_run", "p_value": None}}
    tree_rows, labels = balanced_tree(taxa)
    tree = SpeciesTree(tree_rows)
    foreground, clade = foreground_clade(tree_rows, labels)
    background_labels = [label for label in labels
                         if label not in set(labels[:len(labels) // 4])]
    missing_tip = background_labels[-1] if mask == "missing-tip" else None
    catalogue_data, spaces = _catalogue_spaces(scenario)
    rngs = [np.random.default_rng(child)
            for child in np.random.SeedSequence([seed, replicate]).spawn(4)]
    family_units, observed_units, latent_by_gene = {}, {}, {}
    gene_scores = []
    row["genes"] = gene_scores
    for index, (gene_rate, rng) in enumerate(zip(gene_rates, rngs)):
        gene_id = f"gene_{index + 1:02d}"
        gene_result = {"gene_id": gene_id, "generating_mu": float(gene_rate),
            "joint_history": {"status": "not_scored"},
            "independent_history": {"status": "not_scored"}}
        gene_scores.append(gene_result)
        try:
            template_index = index % len(catalogue_data)
            unit = _unit_for_gene(catalogue_data[template_index], spaces[template_index],
                                  tree, labels, gene_id)
            truth_model = RateModel({kind: float(gene_rate) for kind in EDIT_KINDS},
                foreground=foreground, foreground_multiplier=float(rho))
            sampled, assigned, origins = sample_unit_history(unit, truth_model, rng)
            observed = _mask_unit(sampled, missing_tip)
            _validate_tip_support(observed)
        except (ValueError, ArithmeticError, FloatingPointError, OverflowError,
                np.linalg.LinAlgError) as exc:
            gene_result["simulation_status"] = "failed"
            gene_result["simulation_error"] = {"type": type(exc).__name__, "message": str(exc)}
            continue
        gene_result["simulation_status"] = "complete"
        family_units[gene_id] = [_unit_record(observed)]
        observed_units[gene_id] = observed
        latent_by_gene[gene_id] = (assigned, origins, sampled)
    row["simulation_status"] = "complete" if len(observed_units) == 4 else "incomplete"

    comparison = None
    if len(observed_units) == 4:
        foreground_started = perf_counter()
        try:
            comparison = compare_foreground(family_units, foreground,
                branch_length_mode="supplied", max_origins=None,
                profile_multipliers=(.5, 1., 2.), threads=1)
            alt = comparison.get("alternative", {})
            row["foreground_comparison"] = {
                "status": alt.get("status"),
                "foreground_multiplier": alt.get("foreground_multiplier"),
                "best_evaluated_multiplier": alt.get("best_evaluated_multiplier"),
                "comparison_statistic": comparison.get("comparison_statistic"),
                "p_value": None, "null_status": comparison.get("null", {}).get("status"),
                "null_log_likelihood": comparison.get("null", {}).get("log_likelihood"),
                "alternative_log_likelihood": alt.get("log_likelihood"),
                "resolution_reasons": alt.get("resolution_reasons"),
                "profile_curvature": alt.get("profile_curvature"),
                "upper_tail_diagnostic": alt.get("upper_tail_diagnostic"),
                "profile": comparison.get("profile", []),
                "calibration_status": comparison.get("calibration_status")}
        except (ValueError, ArithmeticError, FloatingPointError, OverflowError,
                np.linalg.LinAlgError) as exc:
            row["foreground_comparison"] = {"status": "failed", "p_value": None,
                "error_type": type(exc).__name__, "error": str(exc)}
        row["foreground_comparison"]["elapsed_seconds"] = perf_counter() - foreground_started
    profile = _profile_record(comparison, rho) if comparison else None
    from ..verification.exon_independent import (
        fit_independent_exon_characters,
        score_independent_endpoint_posteriors,
    )
    for gene_id, gene_result in zip(
            (f"gene_{index + 1:02d}" for index in range(4)), gene_scores):
        observed = observed_units.get(gene_id)
        latent = latent_by_gene.get(gene_id)
        if observed is None or latent is None:
            gene_result["joint_history"] = {"status": "simulation_failed", "branch_count": 0}
            gene_result["independent_history"] = {"status": "simulation_failed", "branch_count": 0}
            continue
        assigned, _origins, sampled = latent
        gene_fit = (profile.get("families", {}).get(gene_id, {})
                    if profile is not None else {})
        mu = gene_fit.get("mu")
        gene_result["joint_fit_status"] = gene_fit.get("status", "comparison_unavailable")
        gene_result["joint_fitted_mu"] = mu
        gene_result["joint_fit_diagnostic"] = gene_fit.get("diagnostic")
        try:
            if profile is None:
                raise ValueError("Fixed-rho native profile is unavailable")
            if not _scoreable_native_fit(gene_fit):
                raise ValueError("No estimable fixed-rho gene rate for native history score")
            joint = native_endpoint_scores(observed, assigned, float(mu), foreground, rho)
            logs = [entry["true_visible_pair_log_score"] for entry in joint]
            zeros = sum(entry["zero_probability_true_pair"] for entry in joint)
            gene_result["joint_history"] = {"status": "scored",
                "branch_count": len(joint), "mean_endpoint_brier": _mean(
                    [entry["brier"] for entry in joint]),
                "mean_true_visible_pair_log_score": None if zeros else _mean(logs),
                "log_score_status": "negative_infinity" if zeros else "finite",
                "zero_probability_true_pairs": zeros}
        except (ValueError, ArithmeticError, FloatingPointError, OverflowError,
                np.linalg.LinAlgError) as exc:
            gene_result["joint_history"] = {"status": "failed", "branch_count": 0,
                "mean_endpoint_brier": None, "mean_true_visible_pair_log_score": None,
                "error_type": type(exc).__name__, "error": str(exc)}
        try:
            fit_input = _unit_record(observed)
            independent_fit = fit_independent_exon_characters(
                fit_input, foreground, foreground_multiplier=rho,
                branch_length_mode="supplied")
            gene_result["independent_fit"] = {
                "status": independent_fit.get("status"),
                "rate": independent_fit.get("rate"),
                "log_likelihood": independent_fit.get("log_likelihood"),
                "diagnostic": independent_fit.get("diagnostic"),
                "root_prior": independent_fit.get("root_prior"),
                "rate_units": independent_fit.get("rate_units")}
            truth_pairs = independent_truth_pairs(independent_fit, sampled, assigned, tree)
            score_result = score_independent_endpoint_posteriors(independent_fit, truth_pairs)
            independent_scores = score_result.get("edges", [])
            independent_briers = []
            for entry in independent_scores:
                probability = float(entry["p_any_exon_geometry_endpoint_change"])
                parent, child = entry["parent_node"], entry["child_node"]
                before, after = sampled.space.states[assigned[parent]], sampled.space.states[assigned[child]]
                independent_briers.append((probability - float(before.exons != after.exons)) ** 2)
            independent_logs = [entry.get("true_full_endpoint_log_probability")
                                for entry in independent_scores]
            zero_logs = sum(value is None or not math.isfinite(float(value))
                            for value in independent_logs)
            gene_result["independent_history"] = {"status": score_result.get(
                "status", independent_fit.get("status")),
                "branch_count": len(independent_scores),
                "mean_endpoint_brier": _mean(independent_briers) if independent_scores else None,
                "mean_true_visible_pair_log_score": None if zero_logs else _mean(independent_logs),
                "log_score_status": "negative_infinity" if zero_logs else "finite",
                "zero_probability_true_pairs": zero_logs,
                "node_illegal_mass": independent_fit.get("node_illegal_mass")}
        except (ValueError, ArithmeticError, FloatingPointError, OverflowError,
                np.linalg.LinAlgError) as exc:
            gene_result["independent_history"] = {"status": "failed", "branch_count": 0,
                "mean_endpoint_brier": None, "mean_true_visible_pair_log_score": None,
                "error_type": type(exc).__name__, "error": str(exc)}
            gene_result.setdefault("independent_fit", {
                "status": "failed", "rate": None,
                "diagnostic": {"type": type(exc).__name__, "message": str(exc)},
                "root_prior": None})

    expected_branches = sum(1 for _ in tree.edges())
    row["joint_history_panel"] = panel_summary(gene_scores, "joint_history", expected_branches)
    row["independent_history_panel"] = panel_summary(
        gene_scores, "independent_history", expected_branches)
    row["paired_panel"] = paired_panel_summary(
        row["joint_history_panel"], row["independent_history_panel"])
    for key, panel in (("native_history", row["joint_history_panel"]),
                       ("independent", row["independent_history_panel"])):
        row[f"{key}_replicate_brier"] = panel.get("mean_endpoint_brier")
        row[f"{key}_replicate_log_score"] = panel.get("mean_true_visible_pair_log_score")
    row["design"] = {"foreground_clade_node": clade,
        "foreground_children": sorted(foreground), "missing_tip": missing_tip}
    row["elapsed_seconds"] = perf_counter() - replicate_started
    return row


def calibrate_genomic_exon_comparison(output_dir, scenario, foreground_multiplier,
        observation_mask, gene_rates, taxa=8, replicates=20, seed=101, threads=1):
    """Run one prespecified, fixed-catalogue foreground calibration condition."""
    validate_calibration_parameters(scenario, foreground_multiplier, observation_mask,
        gene_rates, taxa, replicates, seed, threads)
    tree_rows, labels = balanced_tree(taxa)
    tree = SpeciesTree(tree_rows)
    foreground, clade = foreground_clade(tree_rows, labels)
    background = [label for label in labels if label not in set(labels[:len(labels) // 4])]
    missing_tip = background[-1] if observation_mask == "missing-tip" else None
    metadata = {"schema": "intraphy.exon-comparison-calibration/1",
        "scenario": scenario, "observation_unit": "genomic_exon_spans",
        "generating_foreground_multiplier": float(foreground_multiplier),
        "observation_mask": observation_mask,
        "missing_tip": missing_tip, "taxa": taxa, "taxon_labels": labels,
        "tree": tree_rows, "root_to_tip_height": 1.0,
        "foreground_clade_node": clade, "foreground_children": sorted(foreground),
        "background_children": sorted(set(tree.parent) - {tree.root} - set(foreground)),
        "gene_rates": {f"gene_{index + 1:02d}": float(rate)
                       for index, rate in enumerate(gene_rates)},
        "genes_per_replicate": 4, "local_units_per_gene": 1,
        "replicates_requested": replicates, "seed": seed,
        "seed_scheme": "numpy SeedSequence([seed, replicate_index]).spawn(4)",
        "catalogue_source": "existing fixed geometry/shared-deletion calibration catalogues",
        "catalogue_discovery": "independent_catalogue",
        "generation": "sample_unit_history with equal elementary edit rates per gene and shared foreground multiplier",
        "fit_inputs": "observed state spaces, tree, and tip likelihood vectors only; latent assignments excluded",
        "history_score_scope": "conditional on generating foreground multiplier; not a general joint-model advantage test",
        "posterior_backend": "dense full finite-state and full origin mixture; no state/origin cap",
        "independent_model": {
            "process": "symmetric binary character CTMC",
            "rate_sharing": "one per-gene rate shared across that gene's features",
            "root_marginals": "matched to native feature marginals; joint root prior differs",
            "reversibility": "DNA gains are allowed; this differs from native material process",
            "evaluation_scope": "conditional on native-model generated histories only"},
        "partial_candidate_observations": "rejected; only exact calls and whole-tip unknown masking are supported",
        "failed_replicates_retained": True, "p_values_reported": False,
        "confidence_intervals_reported": False}
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    worker_count = min(threads, replicates)
    metadata["worker_count"] = worker_count
    logger = logging.getLogger("intraphy")
    logger.info("Starting exon comparison calibration: %s replicates, %s workers",
                replicates, worker_count)
    write_json(out / "calibration_metadata.json", metadata)
    tasks = [(index, scenario, float(foreground_multiplier), observation_mask,
              taxa, tuple(float(value) for value in gene_rates), seed)
             for index in range(1, replicates + 1)]
    records = []
    path = out / "replicates.jsonl"
    if worker_count == 1:
        generated = map(_replicate, tasks)
        with path.open("w", encoding="utf-8") as handle:
            for record in generated:
                safe = json_safe(record)
                handle.write(json.dumps(safe, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n")
                handle.flush()
                records.append(record)
                logger.info("Completed exon comparison calibration replicate %d/%d",
                            record["replicate"], replicates)
    else:
        with ProcessPoolExecutor(max_workers=worker_count,
                mp_context=get_context("spawn")) as executor, \
             path.open("w", encoding="utf-8") as handle:
            futures = [executor.submit(_replicate, task) for task in tasks]
            for future in as_completed(futures):
                record = future.result()
                safe = json_safe(record)
                handle.write(json.dumps(safe, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n")
                handle.flush()
                records.append(record)
                logger.info("Completed exon comparison calibration replicate %d/%d",
                            record["replicate"], replicates)
    records.sort(key=lambda record: record["replicate"])
    summary = summarize_calibration(records, replicates)
    write_json(out / "summary.json", summary)
    return summary


def validate_calibration_parameters(scenario, foreground_multiplier, observation_mask,
                                    gene_rates, taxa, replicates, seed, threads):
    if scenario not in {"geometry", "shared-deletion"}:
        raise ValueError("scenario must be 'geometry' or 'shared-deletion'")
    rho = float(foreground_multiplier)
    if not math.isfinite(rho) or rho < 0:
        raise ValueError("foreground multiplier must be finite and nonnegative")
    if rho not in {.5, 1., 2.}:
        raise ValueError("foreground multiplier must be one of the prespecified values 0.5, 1, or 2")
    if observation_mask not in {"complete", "missing-tip"}:
        raise ValueError("observation mask must be 'complete' or 'missing-tip'")
    if len(gene_rates) != 4:
        raise ValueError("gene-rates must contain exactly four per-gene rates")
    if any(not math.isfinite(float(value)) or float(value) < 0 for value in gene_rates):
        raise ValueError("gene rates must be finite and nonnegative")
    balanced_tree(taxa)
    if type(replicates) is not int or replicates < 1:
        raise ValueError("replicates must be a positive integer")
    if type(seed) is not int or seed < 0:
        raise ValueError("seed must be a nonnegative integer")
    if type(threads) is not int or threads < 1:
        raise ValueError("threads must be a positive integer")
