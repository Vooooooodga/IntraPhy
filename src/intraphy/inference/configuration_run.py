"""V19 configuration inference orchestration. Does not rewrite correspondence."""
from __future__ import annotations
from collections import defaultdict
from dataclasses import asdict, replace
import json
from pathlib import Path
import shutil
import tempfile
import numpy as np

from ..run_result import RunResult
from ..storage.tabular import read_tsv, write_tsv
from ..topology import SpeciesTree
from ..structure import MODEL_VERSION
from ..structure.space import enumerate_space
from ..structure.validation import validate_collection
from ..structure.tree_context import canonical_tree, tree_rows as normalized_rows
from ..structure.observations import observation_scenarios, observation_scenario_count
from ..structure.serialization import read_catalogues, write_catalogues, write_json, json_safe
from ..structure.prepare import prepare_configurations
from ..structure.edits import EDIT_KINDS
from .configuration_history import ReconstructionWorkspace
from .configuration_model import RateModel, evaluate_model


def load_rates(path: str | Path) -> RateModel:
    data = json.loads(Path(path).read_text())
    if data.get("schema") != "intraphy.exon-rates/1":
        raise ValueError("Expected intraphy.exon-rates/1; no implicit or estimated rate defaults")
    if not data.get("provenance"):
        raise ValueError("A rate file must describe its source in provenance")
    return RateModel(data["rates"], float(data.get("scale", 1.)),
                     float(data.get("foreground_multiplier", 1.)), frozenset(data.get("foreground_children", ())),
                     float(data.get("origin_root_weight", 1.)))


def _state_rows(space):
    return [{"state_id": f"C{i:04d}", "exons": [[e.start, e.end] for e in s.exons],
             "material": list(s.material), "key": s.key} for i, s in enumerate(space.states)]


def _robust_events(histories):
    if not histories:
        return []
    maps = [{(e["parent"], e["child"], e["edit_key"]): e for e in h.get("events", ())} for h in histories]
    keys = set().union(*(m.keys() for m in maps))
    output = []
    for key in sorted(keys):
        examples = [m[key] for m in maps if key in m]
        record = dict(examples[0])
        lower = min(m.get(key, {}).get("minimum_count", 0) for m in maps)
        upper = max(m.get(key, {}).get("maximum_count", 0) for m in maps)
        record.update(minimum_count=lower, maximum_count=upper,
                      support="required" if lower else "possible",
                      observation_scenarios=len(histories))
        record["consequences"] = sorted({v for e in examples for v in e["consequences"]})
        record["affected_spans"] = sorted({tuple(v) for e in examples for v in e["affected_spans"]})
        output.append(record)
    return output


def _table(directory, name, rows, first):
    # Nested objects have an unambiguous JSON encoding in TSV cells.
    normalized = [{k: json.dumps(json_safe(v), separators=(",", ":")) if isinstance(v, (dict, list, tuple)) else v
                   for k, v in r.items()} for r in rows]
    fields = list(dict.fromkeys(first+[key for row in normalized for key in row]))
    write_tsv(directory/name, normalized, fields)


def infer_configurations(input_dir, output_dir, *, model="exon-parsimony", configurations=None,
                         rates=None, observation_view="evidence", max_states=None, max_origins=None,
                         max_observation_scenarios=None, branch_length_mode="supplied", expected_edits=False,
                         alignment_timeout=600, max_locus_bases=100000,
                         exon_identity=.7, anchor_bases=12, anchor_identity=.8, origin_root_sensitivity=()):
    if (max_observation_scenarios is not None and
            (not isinstance(max_observation_scenarios, int) or isinstance(max_observation_scenarios, bool)
             or max_observation_scenarios <= 0)):
        raise ValueError("max_observation_scenarios must be a positive integer or None")
    if any(not np.isfinite(w) or w <= 0 for w in origin_root_sensitivity):
        raise ValueError("Origin sensitivity weights must be finite and positive")
    if origin_root_sensitivity and model != "exon-ctmc":
        raise ValueError("Origin sensitivity requires exon-ctmc and explicit rates")
    if model not in {"exon-parsimony", "exon-ctmc"}:
        raise ValueError("Use the explicit V19 exon-configuration model identifier")
    if model == "exon-ctmc" and not rates:
        raise ValueError("exon-ctmc requires --exon-rates with explicit externally specified parameters")
    if model != "exon-ctmc" and rates:
        raise ValueError("--exon-rates is only used by exon-ctmc")
    if observation_view not in {"annotation", "evidence"}:
        raise ValueError("Unknown observation view")
    root, out = Path(input_dir), Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    RunResult(model, "exon-configurations", "species_tree.tsv", (), "running").write(out)
    input_tree_rows = read_tsv(root/"species_tree.tsv")
    rate_model = load_rates(rates) if rates else None
    tree_context = canonical_tree(SpeciesTree(input_tree_rows), rate_model.foreground if rate_model else frozenset())
    tree, tree_rows = tree_context.tree, normalized_rows(tree_context.tree)
    if rate_model:
        rate_model = replace(rate_model, foreground=tree_context.foreground)
    write_tsv(out/"species_tree.tsv", tree_rows, ["node_id", "parent_id", "label", "branch_length"])
    write_json(out/"input_tree.json", input_tree_rows)
    write_json(out/"tree_normalization.json", tree_context.diagnostics())
    preparation_reports = []
    if configurations:
        catalogues = read_catalogues(configurations)
        write_catalogues(out/"exon_configurations.jsonl", catalogues)
    else:
        catalogues, preparation_reports = prepare_configurations(root, out,
            timeout=alignment_timeout, max_locus_bases=max_locus_bases,
            minimum_identity=exon_identity, anchor_bases=anchor_bases, anchor_identity=anchor_identity)
    validate_collection(catalogues, allow_empty=True)
    for c in catalogues:
        unexpected = {o.species for o in c.observations} - set(tree.leaf_by_label)
        if unexpected:
            raise ValueError("Configuration contains species absent from the tree: " + ", ".join(sorted(unexpected)))
    taxa = tuple(sorted(tree.leaf_by_label))
    summary, events, witnesses, ancestors, posterior_branches, diagnostic_units, detailed = [], [], [], [], [], [], []
    summary.extend({"family_id": r["family_id"], "unit_id": "not_resolved",
        "status": r["status"], "reason": r["reason"], "minimum_structural_edits": None}
        for r in preparation_reports if not r["units"])
    for c in catalogues:
        prefix = {"family_id": c.family, "unit_id": c.unit}
        info = {**prefix, "status": c.status, "reasons": list(c.reasons), "discovery": c.discovery,
                "input_taxa": len(taxa), "known_structure_taxa": sum(o.kind in {"observed", "partial"} for o in c.observations),
                "coexisting_taxa": sum(o.kind == "coexisting" for o in c.observations),
                "unknown_taxa": sum(o.kind in {"unknown", "excluded"} for o in c.observations)
                    + len(set(taxa)-{o.species for o in c.observations})}
        space = enumerate_space(c, max_states)
        info.update(state_count=len(space.states), state_space_complete=space.complete,
                    state_space_scope="closed_declared_candidate_catalogue_not_all_possible_historical_exons",
                    state_space_reason=space.reason, state_limit=max_states, state_space_estimate=space.diagnostics)
        detail = {**prefix, "catalogue": asdict(c), "states": _state_rows(space), "views": {}}
        theoretical_scenarios = observation_scenario_count(space, taxa)
        info["observation_scenarios"] = {
            view: {"theoretical_count": theoretical_scenarios, "evaluated_count": 0,
                   "resolved_count": 0, "unresolved_count": 0, "limit": max_observation_scenarios,
                   "status": "not_attempted", "scope": "exact_declared_cartesian_coexisting_scenarios"}
            for view in ("evidence", "annotation")}
        if c.status != "qualified" or not space.complete:
            info["inference_status"] = "unresolved" if c.status != "qualified" else "state_space_incomplete"
            info["probability_status"] = "not_available"
            selected_counts = info["observation_scenarios"][observation_view]
            summary.append({**prefix, "status": info["inference_status"], "minimum_structural_edits": None,
                "observation_view": observation_view,
                "observation_scenario_theoretical_count": selected_counts["theoretical_count"],
                "observation_scenario_evaluated_count": selected_counts["evaluated_count"],
                "observation_scenario_resolved_count": selected_counts["resolved_count"],
                "observation_scenario_unresolved_count": selected_counts["unresolved_count"],
                "observation_scenario_limit": selected_counts["limit"],
                "observation_scenario_status": selected_counts["status"],
                "observation_scenario_scope": selected_counts["scope"]})
            diagnostic_units.append(info)
            detailed.append(detail)
            continue
        views = {}
        reconstruction_workspace = ReconstructionWorkspace(space, tree)
        try:
            for view in ("evidence", "annotation"):
                scenarios = observation_scenarios(space, taxa, view, max_observation_scenarios)
                info["observation_scenarios"][view]["status"] = "processing"
                histories = []
                for label, tips in scenarios:
                    informative = sum(not np.all(v == v[0]) for v in tips.values())
                    if informative < 2:
                        histories.append({"status": "insufficient_observed_taxa", "minimum_cost": None,
                                          "minimum_structural_edits": None, "events": [], "witness": [],
                                          "observation_scenario": label})
                    else:
                        result = reconstruction_workspace.reconstruct(tips, max_origins=max_origins)
                        result["observation_scenario"] = label
                        histories.append(result)
                    if histories[-1].get("minimum_structural_edits") is None:
                        info["observation_scenarios"][view]["unresolved_count"] += 1
                    else:
                        info["observation_scenarios"][view]["resolved_count"] += 1
                    info["observation_scenarios"][view]["evaluated_count"] += 1
                info["observation_scenarios"][view]["status"] = "complete"
                views[view] = {"histories": histories, "events": _robust_events(histories)}
                info["observation_scenarios"][view]["workspace"] = reconstruction_workspace.diagnostics()
        except ValueError as exc:
            if "incomplete" not in str(exc):
                raise
            info["observation_scenarios"][view]["status"] = "blocked"
            info["observation_scenarios"][view]["workspace"] = reconstruction_workspace.diagnostics()
            info["inference_status"] = str(exc)
            info["probability_status"] = "not_available"
            selected_counts = info["observation_scenarios"][observation_view]
            summary.append({**prefix, "status": str(exc), "minimum_structural_edits": None,
                "observation_view": observation_view,
                "observation_scenario_theoretical_count": selected_counts["theoretical_count"],
                "observation_scenario_evaluated_count": selected_counts["evaluated_count"],
                "observation_scenario_resolved_count": selected_counts["resolved_count"],
                "observation_scenario_unresolved_count": selected_counts["unresolved_count"],
                "observation_scenario_limit": selected_counts["limit"],
                "observation_scenario_status": selected_counts["status"],
                "observation_scenario_scope": selected_counts["scope"]})
            diagnostic_units.append(info)
            detailed.append(detail)
            continue
        detail["views"] = views
        chosen = views[observation_view]
        minimums = [h.get("minimum_structural_edits") for h in chosen["histories"]]
        valid_values = [v for v in minimums if v is not None]
        info["inference_status"] = "conditional_structure_inference" if valid_values else "insufficient_observed_taxa"
        annotation_values = [h.get("minimum_structural_edits") for h in views["annotation"]["histories"]]
        evidence_values = [h.get("minimum_structural_edits") for h in views["evidence"]["histories"]]
        def event_signature(view):
            return {(e["parent"], e["child"], e["edit_key"], e["minimum_count"], e["maximum_count"])
                    for e in views[view]["events"]}
        info["annotation_sensitivity"] = (annotation_values != evidence_values or
            event_signature("annotation") != event_signature("evidence"))
        summary.append({**prefix, "status": info["inference_status"], "observation_view": observation_view,
            "observation_scenario_theoretical_count": info["observation_scenarios"][observation_view]["theoretical_count"],
            "observation_scenario_evaluated_count": info["observation_scenarios"][observation_view]["evaluated_count"],
            "observation_scenario_resolved_count": info["observation_scenarios"][observation_view]["resolved_count"],
            "observation_scenario_unresolved_count": info["observation_scenarios"][observation_view]["unresolved_count"],
            "observation_scenario_limit": info["observation_scenarios"][observation_view]["limit"],
            "observation_scenario_status": "complete",
            "observation_scenario_scope": info["observation_scenarios"][observation_view]["scope"],
            "minimum_structural_edits": min(valid_values) if valid_values else None,
            "maximum_conditional_minimum": max(valid_values) if valid_values else None,
            "required_edit_placements": sum(e["support"] == "required" for e in chosen["events"]),
            "possible_edit_placements": sum(e["support"] == "possible" for e in chosen["events"]),
            "possible_placements_are_not_additive": True,
            "annotation_conditional_minima": annotation_values, "evidence_compatible_minima": evidence_values,
            "coexisting_structures": info["coexisting_taxa"], "molecular_mutation_count": "not_identified"})
        for i, event in enumerate(chosen["events"], 1):
            events.append({**prefix, "event_id": f"{c.family}/{c.unit}/E{i:05d}", "observation_view": observation_view,
                "interpretation": "required_under_declared_observation_conditions_not_confirmed_mutation",
                "robust_to_allowed_annotation_alternatives": (event["support"] == "required" and
                    any(e["edit_key"] == event["edit_key"] and e["parent"] == event["parent"] and
                        e["child"] == event["child"] and e["support"] == "required"
                        for e in views["evidence"]["events"])), **event})
        for number, history in enumerate(chosen["histories"], 1):
            for witness in history.get("witness", ()):
                witnesses.append({**prefix, "conditional_scenario": number, "representative_only": True, **witness})
            for node, states in history.get("nodes", {}).items():
                ancestors.append({**prefix, "conditional_scenario": number, "node": node, "state_indices": states,
                                  "method": "all_global_parsimony_optima", "probability": None})
        info["probability_status"] = "not_requested"
        if rate_model:
            if any(o.kind == "coexisting" for o in c.observations):
                info["probability_status"] = "coexisting_structures_are_not_a_probability_mixture"
            elif not valid_values:
                info["probability_status"] = "insufficient_observed_taxa"
            else:
                _, tips = observation_scenarios(space, taxa, observation_view, 1)[0]
                result = evaluate_model(space, tree, tips, rate_model, max_origins=max_origins,
                                        counts=expected_edits, branch_length_mode=branch_length_mode)
                detail["ctmc"] = json_safe(result)
                if not np.isfinite(result["log_likelihood"]):
                    info["probability_status"] = "zero_probability_observation_under_supplied_model"
                else:
                    info["probability_status"] = "conditional_on_fixed_parameters_and_declared_catalogue"
                    info["log_likelihood"] = result["log_likelihood"]
                    detail["probability_conditions"] = {"origin_scenarios": result["origins"],
                        "root_prior": result["root_prior"], "origin_prior": result["origin_prior"],
                        "probability_status": info["probability_status"]}
                    detail["origin_prior_sensitivity"] = []
                    for weight in sorted(set(origin_root_sensitivity)):
                        alternate = evaluate_model(space, tree, tips,
                            replace(rate_model, origin_root_weight=weight), max_origins=max_origins,
                            counts=False, branch_length_mode=branch_length_mode)
                        detail["origin_prior_sensitivity"].append({
                            "origin_root_weight": weight, "log_likelihood": alternate["log_likelihood"],
                            "root_probabilities": json_safe(alternate["nodes"].get(tree.root)),
                            "conditional_sensitivity_not_model_selection": True})
                    for node, values in result["nodes"].items():
                        for state, value in enumerate(values):
                            ancestors.append({**prefix, "conditional_scenario": 1, "node": node, "state_indices": [state],
                                              "method": "fixed_parameter_ctmc", "probability": float(value)})
                    for branch in result["branches"]:
                        posterior_branches.append({**prefix, **{k: v for k, v in branch.items() if k != "endpoint_probabilities"}})
        diagnostic_units.append(info)
        detailed.append(detail)
    grouped_summaries = defaultdict(list)
    for row in summary:
        grouped_summaries[row["family_id"]].append(row)
    gene_summary = []
    for family, records in sorted(grouped_summaries.items()):
        resolved = [r for r in records if r.get("minimum_structural_edits") is not None]
        gene_summary.append({"family_id": family, "analyses_reported": len(records),
            "resolved_local_configurations": len(resolved), "unresolved_analyses": len(records)-len(resolved),
            "minimum_edits_within_resolved_units": sum(r["minimum_structural_edits"] for r in resolved) if resolved else None,
            "whole_gene_mutation_count": "not_identified", "scope": "local_conditional_histories_not_a_joint_ancestral_transcript"})
    # Stage final inferential products; publish their owner manifest last.
    with tempfile.TemporaryDirectory(prefix=".v19-results-", dir=out) as tmp:
        staging = Path(tmp)
        write_json(staging/"exon_history.json", {"model": MODEL_VERSION, "engine": model,
            "tree": tree_rows, "observation_view": observation_view, "units": detailed})
        _table(staging, "structural_history.tsv", events, ["family_id", "unit_id", "event_id", "operation", "support"])
        _table(staging, "representative_structural_history.tsv", witnesses, ["family_id", "unit_id", "conditional_scenario", "parent", "child", "order", "operation"])
        _table(staging, "gene_structure_summary.tsv", gene_summary, ["family_id", "resolved_local_configurations", "unresolved_analyses", "minimum_edits_within_resolved_units"])
        _table(staging, "exon_structure_summary.tsv", summary, ["family_id", "unit_id", "status", "minimum_structural_edits"])
        _table(staging, "exon_ancestral_states.tsv", ancestors, ["family_id", "unit_id", "node", "method", "probability"])
        _table(staging, "exon_branch_posteriors.tsv", posterior_branches, ["family_id", "unit_id", "parent", "child", "probability_different_endpoints"])
        write_json(staging/"model_diagnostics.json", {"model_version": MODEL_VERSION, "engine": model,
            "purpose": "exon_structure_evolution_only", "biological_accuracy_calibrated": False,
            "finite_sample_significance_calibrated": False, "formal_pvalues": "not_produced_by_single_family_analysis",
            "observation_view": observation_view, "parameters": asdict(rate_model) if rate_model else None,
            "rate_unit": "per_eligible_edit_opportunity_per_branch_length_unit", "branch_length_mode": branch_length_mode,
            "tree_normalization": tree_context.diagnostics(),
            "root_structure_prior": "uniform_valid_configurations_given_material_origin_scenario",
            "material_origin_prior": "declared_root_weight_and_unit_branch_opportunities_on_canonical_tree",
            "origin_root_sensitivity_weights": list(origin_root_sensitivity),
            "deletion_reversible": False, "transcript_usage_modelled": False,
            "preparation": preparation_reports, "units": diagnostic_units})
        for path in staging.iterdir():
            path.replace(out/path.name)
    artifacts = ("exon_configurations.jsonl", "exon_history.json", "structural_history.tsv",
                 "representative_structural_history.tsv", "exon_structure_summary.tsv", "gene_structure_summary.tsv", "exon_ancestral_states.tsv",
                 "exon_branch_posteriors.tsv", "model_diagnostics.json")
    RunResult(model, "exon-configurations", "species_tree.tsv", artifacts).write(out)
    return summary
