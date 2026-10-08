"""Per-family conditional genomic exon-span CTMC fitting and exports."""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import numpy as np

from ..run_result import RunResult
from ..storage.tabular import write_tsv
from ..structure.prepare import prepare_configurations
from ..structure.serialization import read_catalogues, write_catalogues
from ..structure.tree_context import canonical_tree, tree_rows as normalized_rows
from ..structure.validation import validate_collection
from ..inputs.species_tree import read_species_tree_rows, select_species_tree_rows
from ..topology import SpeciesTree
from .configuration_model import RateModel
from .genomic_exon_family import run_family_tasks
from .genomic_exon_output import assemble_rows, write_outputs
from .genomic_exon_branches import branch_change_rows
from .genomic_exon_inputs import _diagnose_catalogue, _write_state_space_progress


MODEL = "exon-structure-ctmc"


def infer_genomic_exons(input_dir, output_dir, *, configurations=None, rates=None,
                        parameter_mode="fit", species_tree=None, max_states=None,
                        max_origins=None, branch_length_mode="supplied", expected_edits=False,
                        alignment_timeout=600, max_locus_bases=100000, exon_identity=.7,
                        anchor_bases=12, anchor_identity=.8, threads=1,
                        alignment_evidence_dir=None, state_space_only=False):
    """Infer physical exon-span histories with one shared rate per family."""
    if configurations and alignment_evidence_dir:
        raise ValueError("alignment_evidence_dir and exon configurations are mutually exclusive")
    if parameter_mode not in {"fit", "fixed"}:
        raise ValueError("parameter_mode must be 'fit' or 'fixed'")
    if parameter_mode == "fixed" and not rates:
        raise ValueError("Fixed parameter_mode requires an explicit exon rate file")
    if parameter_mode == "fit" and rates:
        raise ValueError("Rate files are accepted only when parameter_mode='fixed'")
    if state_space_only and max_states is not None:
        raise ValueError("--state-space-only enumerates complete spaces; omit --max-states")
    if branch_length_mode not in {"supplied", "unit"}:
        raise ValueError("branch_length_mode must be 'supplied' or 'unit'")
    if type(threads) is not int or threads < 1:
        raise ValueError("threads must be a positive integer")

    root, out = Path(input_dir), Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    RunResult(MODEL, "genomic-exon-spans", "species_tree.tsv", (), "running").write(out)
    if configurations:
        catalogues = read_catalogues(configurations)
        if any(c.observation_unit != "genomic_exon_spans" for c in catalogues):
            raise ValueError("Explicit catalogues must declare observation_unit='genomic_exon_spans'; regenerate transcript catalogues")
        preparation_reports = []
        write_catalogues(out/"exon_configurations.jsonl", catalogues)
    else:
        catalogues, preparation_reports = prepare_configurations(root, out,
            timeout=alignment_timeout, max_locus_bases=max_locus_bases,
            minimum_identity=exon_identity, anchor_bases=anchor_bases,
            anchor_identity=anchor_identity, observation_unit="genomic_exon_spans",
            threads=threads, alignment_evidence_dir=alignment_evidence_dir)
    catalogues = validate_collection(tuple(catalogues), allow_empty=True)
    if any(c.observation_unit != "genomic_exon_spans" for c in catalogues):
        raise ValueError("Genomic exon inference requires genomic_exon_spans catalogues")

    prepared_rows = read_species_tree_rows(root/"species_tree.tsv")
    prepared_tree = SpeciesTree(prepared_rows)
    tree_source = species_tree if species_tree is not None else root/"species_tree.tsv"
    input_tree_rows = read_species_tree_rows(tree_source)
    input_tree = SpeciesTree(input_tree_rows)
    if set(input_tree.leaf_by_label) != set(prepared_tree.leaf_by_label):
        raise ValueError("species_tree override must have exactly the prepared tip set")
    # The existing selector validates and retains the exact input tip panel and
    # every original node path; canonical_tree then collapses only redundant unary nodes.
    exact_rows = select_species_tree_rows(tree_source, sorted(prepared_tree.leaf_by_label), prune=False)
    rate_model = None
    if parameter_mode == "fixed":
        from .configuration_run import load_rates
        rate_model = load_rates(rates)
    foreground = rate_model.foreground if rate_model else frozenset()
    context = canonical_tree(SpeciesTree(exact_rows), foreground)
    tree = context.tree
    if rate_model is not None:
        from dataclasses import replace
        rate_model = replace(rate_model, foreground=context.foreground)
    taxa = tuple(sorted(tree.leaf_by_label))
    unexpected = {o.species for c in catalogues for o in c.observations} - set(taxa)
    if unexpected:
        raise ValueError("Catalogue contains species absent from the tree: " + ", ".join(sorted(unexpected)))
    tree_table = normalized_rows(tree)
    write_tsv(out/"species_tree.tsv", tree_table,
              ["node_id", "parent_id", "label", "branch_length"])

    unresolved, diagnostic_units, family_units = [], [], defaultdict(list)
    unit_details = {}
    with (out/"state_space_diagnostics.jsonl").open("w", encoding="utf-8") as progress:
        for catalogue in catalogues:
            key = (catalogue.family, catalogue.unit)
            diag, unit_detail, unit_input = _diagnose_catalogue(
                catalogue, taxa, tree, None if state_space_only else max_states,
                state_space_only, progress)
            if unit_input is not None:
                family_units[catalogue.family].append(unit_input)
                unit_details[key] = unit_detail
            if diag["status"] == "eligible_with_unknown_tips":
                unresolved.append({"family_id": catalogue.family,
                    "unit_id": catalogue.unit, "reason": "unknown_tip_observations"})
            elif diag["status"] not in {"eligible_conditional_unit", "eligible_with_unknown_tips"}:
                unresolved.append({"family_id": catalogue.family, "unit_id": catalogue.unit,
                                   "reason": diag["status"]})
            diagnostic_units.append(diag)

    if state_space_only:
        artifacts = ["exon_configurations.jsonl", "species_tree.tsv",
                     "state_space_diagnostics.jsonl"]
        if not configurations:
            artifacts.extend(["exon_correspondence.tsv", "exon_coordinates.tsv",
                "exon_preparation_summary.tsv", "annotation_structure_candidates.tsv",
                "exon_evidence_policy.json", "native_cds_consequences.json"])
        RunResult(MODEL, "genomic-exon-state-space-diagnostic", "species_tree.tsv",
                  tuple(artifacts), "state_space_only").write(out)
        return diagnostic_units

    worker_parameters = {"parameter_mode": parameter_mode, "rates": rate_model,
        "branch_length_mode": branch_length_mode, "max_origins": max_origins,
        "expected_edits": expected_edits}
    tasks = [(family, units, worker_parameters)
             for family, units in sorted(family_units.items())]
    family_results, effective_threads, effective_unit_workers = run_family_tasks(tasks, threads)
    fit_by_family = {item["family_id"]: item["fit"] for item in family_results}
    result_by_unit = {(result["family_id"], result["unit_id"]): result
                      for family_result in family_results for result in family_result["units"]}
    for diagnostic in diagnostic_units:
        unit_key = (diagnostic["family_id"], diagnostic["unit_id"])
        if diagnostic["status"] in {"eligible_conditional_unit", "eligible_with_unknown_tips"}:
            diagnostic["family_rate_status"] = fit_by_family.get(unit_key[0], {}).get(
                "status", "no_estimable_units")
            diagnostic["status"] = result_by_unit.get(unit_key, {}).get("status",
                diagnostic["family_rate_status"])
            diagnostic["probability_status"] = ("conditional_on_declared_model"
                if result_by_unit.get(unit_key, {}).get("ctmc") is not None
                and result_by_unit[unit_key].get("log_likelihood") is not None
                and np.isfinite(result_by_unit[unit_key]["log_likelihood"])
                and result_by_unit[unit_key]["ctmc"].get("nodes") else "not_available")
    all_families = {c.family for c in catalogues} | {r.get("family_id") for r in preparation_reports
                   if r.get("family_id")}
    if not tasks:
        unresolved.append({"family_id": "NA", "unit_id": "NA", "reason": "no_estimable_units"})
    fit_by_family, details, summaries, ancestral, branches = assemble_rows(
        family_results, unit_details, diagnostic_units, unresolved, all_families)
    branch_exon_changes = branch_change_rows(family_results, unit_details, tree)
    fit_record = {"model": MODEL, "parameter_mode": parameter_mode,
        "objective_scope": "conditional_composite_likelihood_within_family_across_local_units",
        "families": fit_by_family,
        "requested_threads": threads, "effective_family_workers": effective_threads,
        "effective_unit_workers": effective_unit_workers}
    completed_status = "completed_with_unresolved" if unresolved else "completed"
    diagnostics = {"model": MODEL,
        "method_scope": "constrained_elementary_edit_graph_conditional_composite_likelihood",
        "ascertainment_correction": "not_applied",
        "observation_unit": "distinct_nonoverlapping_physical_exon_span_union_per_species",
        "transcript_usage_modelled": False, "branch_length_mode": branch_length_mode,
        "rate_unit": "per_eligible_elementary_edit_opportunity_per_branch_length_unit",
        "rate_mode": parameter_mode,
        "rates": ({"rates": rate_model.rates, "scale": rate_model.scale,
                   "foreground_multiplier": rate_model.foreground_multiplier,
                   "foreground": sorted(rate_model.foreground),
                   "origin_root_weight": rate_model.origin_root_weight} if rate_model else
                  {family: fit_by_family[family].get("mu") for family in fit_by_family}),
        "root_structure_prior": "uniform_valid_exon_geometries_given_material_origin_scenario",
        "origin_root_weight": rate_model.origin_root_weight if rate_model else 1.0,
        "branch_origin_opportunity_weight": 1.0,
        "origin_prior": "declared_root_weight_and_unit_weight_per_canonical_branch_opportunity_before_observation",
        "tree_normalization": context.diagnostics(), "input_species_tree": input_tree_rows,
        "deletion_reversible": False, "threads": {"requested": threads,
            "effective_family_workers": effective_threads,
            "effective_unit_workers": effective_unit_workers},
        "preparation": preparation_reports, "units": diagnostic_units,
        "unresolved": unresolved}
    artifacts = write_outputs(out, tree=tree_table, details=details, summaries=summaries,
        ancestral=ancestral, branches=branches, branch_exon_changes=branch_exon_changes,
        expected_edits=expected_edits,
        fit_record=fit_record, diagnostics=diagnostics,
        preparation_artifacts=not bool(configurations))
    artifacts = tuple(artifacts) + ("state_space_diagnostics.jsonl",)
    RunResult(MODEL, "genomic-exon-spans", "species_tree.tsv", tuple(artifacts),
              completed_status).write(out)
    return summaries
