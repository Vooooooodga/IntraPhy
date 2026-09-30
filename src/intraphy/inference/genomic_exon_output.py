"""TSV and JSON organization for genomic exon-span inference results."""
from __future__ import annotations

import json
import numpy as np

from ..storage.tabular import write_tsv
from ..structure.serialization import json_safe, write_json


def state_rows(space):
    return [{"state_id": f"C{i:04d}", "exons": [[e.start, e.end] for e in state.exons],
             "material": list(state.material), "key": state.key}
            for i, state in enumerate(space.states)]


def _table(directory, name, rows, columns):
    normalized = []
    for row in rows:
        safe = {key: json_safe(value) for key, value in row.items()}
        normalized.append({key: json.dumps(value, separators=(",", ":"))
                           if isinstance(value, (dict, list, tuple)) else value
                           for key, value in safe.items()})
    fields = list(dict.fromkeys(columns + [key for row in normalized for key in row]))
    write_tsv(directory/name, normalized, fields)


def assemble_rows(family_results, unit_details, diagnostic_units, unresolved, all_families):
    details, summaries, ancestral, branches = [], [], [], []
    fit_by_family = {item["family_id"]: item["fit"] for item in family_results}
    for family_result in family_results:
        family = family_result["family_id"]
        for result in family_result["units"]:
            key = (family, result["unit_id"])
            details.append({**unit_details[key], **result})
            summaries.append({"family_id": family, "unit_id": result["unit_id"],
                "status": result["status"], "log_likelihood": result["log_likelihood"],
                "unknown_tips": unit_details[key]["unknown_tips"]})
            ctmc = result["ctmc"]
            if ctmc and result["log_likelihood"] is not None and np.isfinite(result["log_likelihood"]):
                for node, probabilities in ctmc["nodes"].items():
                    for state, probability in enumerate(probabilities):
                        ancestral.append({"family_id": family, "unit_id": result["unit_id"],
                            "node": node, "state_id": f"C{state:04d}", "state_index": state,
                            "probability": probability})
                for branch in ctmc["branches"]:
                    branches.append({"family_id": family, "unit_id": result["unit_id"],
                        **{key: value for key, value in branch.items()
                           if key != "endpoint_probabilities"}})
            if result["status"] in {"unresolved_fit", "impossible_observation_under_model",
                    "origin_scenarios_incomplete", "zero_probability_under_supplied_or_fitted_parameters"}:
                unresolved.append({"family_id": family, "unit_id": result["unit_id"],
                                   "reason": result["status"]})
        fit = family_result["fit"]
        if fit.get("mu") is None and fit.get("status") != "fixed_parameters":
            unresolved.append({"family_id": family, "unit_id": "family_fit", "reason": fit["status"]})
    for family in sorted(all_families):
        if family not in fit_by_family:
            fit_by_family[family] = {"status": "no_estimable_units", "converged": False,
                                     "mu": None, "unit_count": 0}
            summaries.append({"family_id": family, "unit_id": "family_fit",
                "status": "no_estimable_units", "log_likelihood": None, "unknown_tips": []})
            unresolved.append({"family_id": family, "unit_id": "family_fit", "reason": "no_estimable_units"})
    summary_lookup = {(row["family_id"], row["unit_id"]): row for row in summaries}
    for diagnostic in diagnostic_units:
        key = (diagnostic["family_id"], diagnostic["unit_id"])
        if key not in summary_lookup:
            summaries.append({"family_id": key[0], "unit_id": key[1],
                "status": diagnostic["status"], "log_likelihood": None,
                "unknown_tips": diagnostic.get("unknown_tips", [])})
            summary_lookup[key] = summaries[-1]
        summary_lookup[key]["status"] = diagnostic["status"]
    for row in unresolved:
        key = (row["family_id"], row["unit_id"])
        if key in summary_lookup and summary_lookup[key]["status"] in {
                "eligible_conditional_unit", "eligible_with_unknown_tips"}:
            summary_lookup[key]["status"] = row["reason"]
    return fit_by_family, details, summaries, ancestral, branches


def write_outputs(directory, *, tree, details, summaries, ancestral, branches,
                  branch_exon_changes, expected_edits, fit_record, diagnostics,
                  preparation_artifacts=False):
    write_json(directory/"exon_structure_fit.json", fit_record)
    write_json(directory/"exon_history.json", {"model": "exon-structure-ctmc",
        "tree": tree, "observation_unit": "genomic_exon_spans", "units": details})
    _table(directory, "exon_structure_summary.tsv", summaries,
           ["family_id", "unit_id", "status", "log_likelihood", "unknown_tips"])
    _table(directory, "ancestral_exon_states.tsv", ancestral,
           ["family_id", "unit_id", "node", "state_id", "state_index", "probability"])
    _table(directory, "branch_exon_changes.tsv", branch_exon_changes,
           ["family_id", "unit_id", "parent_node_id", "child_node_id",
           "descendant_species", "probability_exon_structure_change",
            "probability_dna_presence_change", "exon_count_pair_probabilities",
            "probability_exon_count_increase", "probability_exon_count_decrease",
            "probability_exon_count_unchanged", "dna_presence_scope",
            "joint_configuration_probability",
            "parent_exons", "child_exons", "parent_dna_presence", "child_dna_presence",
            "material_ids", "changed_material_tracts",
            "change_classification", "coordinate_system", "alignment_offset",
            "interpretation_scope"])
    artifacts = ["exon_structure_fit.json", "exon_history.json", "exon_structure_summary.tsv",
        "ancestral_exon_states.tsv", "branch_exon_changes.tsv", "model_diagnostics.json",
        "species_tree.tsv", "exon_configurations.jsonl"]
    if preparation_artifacts:
        artifacts.extend(["exon_correspondence.tsv", "exon_coordinates.tsv",
                          "exon_preparation_summary.tsv", "annotation_structure_candidates.tsv",
                          "exon_evidence_policy.json", "native_cds_consequences.json"])
    if expected_edits:
        _table(directory, "branch_exon_events.tsv", branches,
            ["family_id", "unit_id", "parent", "child", "probability_different_endpoints",
             "probability_at_least_one_edit", "expected_edits"])
        artifacts.append("branch_exon_events.tsv")
    write_json(directory/"model_diagnostics.json", diagnostics)
    return tuple(artifacts)
