"""Serialization for DNA presence CTMC results."""
from __future__ import annotations

import json
from pathlib import Path

from intraphy.storage.tabular import write_tsv


MODEL = "dna-presence-ctmc"


def write_dna_outputs(output_dir, tree, fit_rows, history_rows):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    tree_path = output_dir / "species_tree.tsv"
    fit_path = output_dir / "dna_fit.json"
    history_path = output_dir / "dna_history.json"
    tree_rows = [{
        "node_id": node,
        "parent_id": tree.parent[node],
        "label": tree.label[node],
        "branch_length": "NA" if tree.length[node] is None else tree.length[node],
    } for node in tree.parent]
    write_tsv(tree_path, tree_rows, ["node_id", "parent_id", "label", "branch_length"])
    fit_data = {
        "schema_version": 1, "model": MODEL,
        "analysis_status": "no_candidate_sites" if not fit_rows else "completed_with_diagnostics",
        "assumptions": {
            "states": {"0": "DNA absent at the homologous position", "1": "DNA present at the homologous position"},
            "gain_interpretation": "structural material gain; source and mechanism are unresolved",
            "rate_sharing": "gain and loss rates shared among eligible homologous sites within each family",
            "ascertainment": "conditioned per site on at least one observed present state",
            "missingness": "observation masks are conditioned on; nonrandom annotation, mapping, and discovery biases are not fully corrected",
            "independence": "site likelihoods are treated as independent; for linked sites their sum is a composite likelihood",
            "event_interpretation": "expected CTMC site-state changes are not counts of molecular lesions",
            "inference_limit": "rates are conditional estimates; no likelihood-ratio tests, confidence intervals, or guaranteed global-optimum claim are reported",
        },
        "families": fit_rows,
    }
    history_data = {
        "schema_version": 1, "model": MODEL,
        "interpretation": "node and edge quantities describe DNA presence-state histories conditional on fitted rates and observed tip states",
        "families": history_rows,
    }
    _write_json(fit_path, fit_data)
    _write_json(history_path, history_data)
    artifacts = [tree_path.name, fit_path.name, history_path.name]
    run_data = {
        "schema_version": 1,
        "status": "completed_with_unresolved" if not fit_rows or any(
            row.get("status") not in {"estimated", "fixed_rates"} for row in fit_rows
        ) else "completed",
        "model": MODEL, "analysis_scope": "single-copy-dna-presence",
        "tree_file": tree_path.name, "artifacts": artifacts,
        "diagnostics": ["no_candidate_sites"] if not fit_rows else [],
    }
    _write_json(output_dir / "run_result.json", run_data)
    return fit_path, history_path


def _write_json(path, data):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)
