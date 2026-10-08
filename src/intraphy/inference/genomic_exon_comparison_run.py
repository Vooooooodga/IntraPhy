"""Prepare native genomic exon units and export a shared-foreground comparison."""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

from ..inputs.foreground import read_canonical_foreground
from ..inputs.species_tree import read_species_tree_rows
from ..run_result import RunResult
from ..storage.tabular import write_tsv
from ..structure.serialization import json_safe, read_catalogues, write_catalogues, write_json
from ..structure.tree_context import canonical_tree, tree_rows as normalized_rows
from ..structure.validation import validate_collection
from ..topology import SpeciesTree
from .genomic_exon_run import MODEL
from .genomic_exon_inputs import _diagnose_catalogue
from .genomic_exon_comparison import compare_foreground


def _json_cell(value):
    safe = json_safe(value)
    return json.dumps(safe, sort_keys=True, separators=(",", ":")) if isinstance(
        safe, (dict, list, tuple)) else safe


def compare_genomic_exon_foreground(input_dir, configurations, foreground_branches,
                                    output_dir, *, branch_length_mode="supplied",
                                    max_states=None, max_origins=None,
                                    profile_multipliers=(), threads=1):
    root, out = Path(input_dir), Path(output_dir)
    catalogues = tuple(c for path in configurations for c in read_catalogues(path))
    if not catalogues:
        raise ValueError("At least one exon catalogue record is required")
    catalogues = validate_collection(catalogues, allow_empty=False)
    if any(c.observation_unit != "genomic_exon_spans" for c in catalogues):
        raise ValueError("Catalogues must declare observation_unit='genomic_exon_spans'")

    input_tree_rows = read_species_tree_rows(root / "species_tree.tsv")
    input_tree = SpeciesTree(input_tree_rows)
    context = canonical_tree(input_tree)
    tree = context.tree
    taxa = tuple(sorted(tree.leaf_by_label))
    unexpected = {o.species for c in catalogues for o in c.observations} - set(taxa)
    if unexpected:
        raise ValueError("Catalogue contains species absent from the prepared tree: "
                         + ", ".join(sorted(unexpected)))
    foreground, foreground_provenance, foreground_exposure, background_exposure = \
        read_canonical_foreground(foreground_branches, context, branch_length_mode)

    out.mkdir(parents=True, exist_ok=True)
    RunResult(MODEL, "genomic-exon-foreground-comparison", "species_tree.tsv", (),
              "running").write(out)
    write_catalogues(out / "exon_configurations.jsonl", catalogues)
    write_tsv(out / "species_tree.tsv", normalized_rows(tree),
              ["node_id", "parent_id", "label", "branch_length"])

    diagnostics = []
    family_units = defaultdict(list)
    with (out / "state_space_diagnostics.jsonl").open("w", encoding="utf-8") as progress:
        for catalogue in catalogues:
            diag, _detail, unit_input = _diagnose_catalogue(
                catalogue, taxa, tree, max_states, False, progress,
                include_detail=False)
            diagnostics.append(diag)
            if unit_input is not None:
                family_units[catalogue.family].append(unit_input)
            else:
                family_units[catalogue.family]

    comparison = compare_foreground(
        dict(family_units), foreground, branch_length_mode=branch_length_mode,
        max_origins=max_origins, profile_multipliers=profile_multipliers,
        threads=threads,
    )
    input_scope = {
        "prepared_input_dir": str(root.resolve()),
        "species_tree_file": "species_tree.tsv",
        "catalogue_files": [str(Path(path).resolve()) for path in configurations],
        "foreground_file": str(Path(foreground_branches).resolve()),
        "catalogue_count": len(catalogues),
        "requested_roster": [{"family_id": c.family, "unit_id": c.unit,
                              "catalogue_status": c.status} for c in catalogues],
        "eligible_roster": [{"family_id": d["family_id"], "unit_id": d["unit_id"],
                             "status": d["status"]} for d in diagnostics
                            if d["status"] in {"eligible_conditional_unit",
                                                "eligible_with_unknown_tips"}],
        "excluded_units": [{"family_id": d["family_id"], "unit_id": d["unit_id"],
                            "status": d["status"]} for d in diagnostics
                           if d["status"] not in {"eligible_conditional_unit",
                                                   "eligible_with_unknown_tips"}],
        "requested_families": sorted({c.family for c in catalogues}),
        "excluded_families": sorted(family for family, units in family_units.items() if not units),
        "discovery_scopes": [{"family_id": c.family, "unit_id": c.unit,
                              "discovery": c.discovery} for c in catalogues],
        "root_structure_prior": "uniform_valid_exon_geometries_given_material_origin_scenario",
        "material_origin_prior": "declared_root_weight_and_unit_weight_per_canonical_branch_opportunity_before_observation",
        "root_origin_weight": 1.0,
        "branch_origin_opportunity_weight": 1.0,
        "tree_normalization": context.diagnostics(),
        "input_species_tree": input_tree_rows,
        "foreground_branch_provenance": foreground_provenance,
        "foreground_exposure": foreground_exposure,
        "background_exposure": background_exposure,
    }
    result = {**comparison, "input_scope": input_scope,
              "observation_unit": "genomic_exon_spans",
              "method_scope": "conditional_composite_likelihood_within_family_across_local_units",
              "transcript_usage_modelled": False,
              "biological_event_identification": "not_estimated"}
    write_json(out / "exon_foreground_comparison.json", result)

    profile_rows = []
    family_fit_rows = []
    for point in comparison.get("profile", []):
        multiplier = point.get("foreground_multiplier")
        profile_rows.append({"foreground_multiplier": multiplier,
            "log_likelihood": point.get("log_likelihood"), "status": point.get("status"),
            "unresolved_families": _json_cell(point.get("unresolved_families", [])),
            "impossible_families": _json_cell(point.get("impossible_families", []))})
    for fit_role in ("null", "alternative"):
        fit = comparison.get(fit_role, {})
        for family, detail in sorted(fit.get("families", {}).items()):
            used_multiplier = fit.get("foreground_multiplier")
            if used_multiplier is None:
                used_multiplier = detail.get("foreground_multiplier")
            best_multiplier = fit.get("best_evaluated_multiplier")
            if best_multiplier is None:
                best_multiplier = used_multiplier
            family_fit_rows.append({"fit_role": fit_role,
                "foreground_multiplier": used_multiplier,
                "best_evaluated_multiplier": best_multiplier,
                "family_id": family, "mu": detail.get("mu"),
                "log_likelihood": detail.get("log_likelihood"), "status": detail.get("status"),
                "diagnostic": _json_cell(detail.get("diagnostic", {}))})
    for point in comparison.get("profile", []):
        for family, detail in sorted(point.get("families", {}).items()):
            family_fit_rows.append({"fit_role": "profile",
                "foreground_multiplier": point.get("foreground_multiplier"),
                "best_evaluated_multiplier": point.get("foreground_multiplier"),
                "family_id": family, "mu": detail.get("mu"),
                "log_likelihood": detail.get("log_likelihood"), "status": detail.get("status"),
                "diagnostic": _json_cell(detail.get("diagnostic", {}))})
    write_tsv(out / "foreground_profile.tsv", profile_rows,
              ["foreground_multiplier", "log_likelihood", "status", "unresolved_families",
               "impossible_families"])
    write_tsv(out / "foreground_family_fits.tsv", family_fit_rows,
              ["fit_role", "foreground_multiplier", "best_evaluated_multiplier", "family_id", "mu",
               "log_likelihood", "status", "diagnostic"])
    write_tsv(out / "foreground_unit_diagnostics.tsv",
              [{key: _json_cell(value) for key, value in row.items()} for row in diagnostics],
              ["family_id", "unit_id", "catalogue_status", "status", "probability_status",
               "state_count", "state_space_complete", "informative_tips", "unknown_tips",
               "state_space_reason", "reason"])
    write_tsv(out / "foreground_branches.tsv",
              [{**row, "input_children_on_canonical_edge": _json_cell(
                  row["input_children_on_canonical_edge"]), "input_row": _json_cell(row["input_row"])}
               for row in foreground_provenance],
              ["row", "branch_scope", "parent_node", "child_node",
               "input_children_on_canonical_edge", "input_row"])

    artifacts = ("exon_foreground_comparison.json", "foreground_profile.tsv",
        "foreground_family_fits.tsv", "foreground_unit_diagnostics.tsv",
        "foreground_branches.tsv", "state_space_diagnostics.jsonl",
        "species_tree.tsv", "exon_configurations.jsonl")
    successful_alternative = comparison.get("alternative", {}).get("status") in {
        "estimated_conditional_composite_multiplier", "zero_boundary"}
    successful_null = comparison.get("null", {}).get("status") == \
        "conditional_composite_log_likelihood"
    status = ("completed" if successful_alternative and successful_null
              and not input_scope["excluded_units"] else "completed_with_unresolved")
    RunResult(MODEL, "genomic-exon-foreground-comparison", "species_tree.tsv",
              artifacts, status).write(out)
    return result
