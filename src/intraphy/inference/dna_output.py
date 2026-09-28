"""Serialization for separate binary-presence CTMC result domains."""
from __future__ import annotations

import json
from pathlib import Path

from intraphy.storage.tabular import write_tsv


MODEL = "dna-presence-ctmc"


def write_binary_presence_outputs(output_dir, tree, fit_rows, history_rows, domain):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    tree_path = output_dir / "species_tree.tsv"
    fit_path = output_dir / f"{domain['file_prefix']}_fit.json"
    history_path = output_dir / f"{domain['file_prefix']}_history.json"
    tree_rows = [{
        "node_id": node, "parent_id": tree.parent[node], "label": tree.label[node],
        "branch_length": "NA" if tree.length[node] is None else tree.length[node],
    } for node in tree.parent]
    write_tsv(tree_path, tree_rows, ["node_id", "parent_id", "label", "branch_length"])
    fit_data = {
        "schema_version": 1, "model": domain["model"],
        "observation_type": domain["observation_type"], "analysis_scope": domain["analysis_scope"],
        "analysis_status": "no_candidate_sites" if not fit_rows else "completed_with_diagnostics",
        "assumptions": {
            "states": {"0": domain["state_0"], "1": domain["state_1"]},
            "transition_0_to_1": domain["transition_01"],
            "transition_1_to_0": domain["transition_10"],
            "rate_sharing": "gain and loss rates shared among eligible homologous sites within each family",
            "ascertainment": "conditioned per site on at least one observed state 1 and on that site's fixed observation mask",
            "root_distribution": "stationary distribution from fitted rates or fixed root_presence, as declared per family",
            "tree_conditioning": "the supplied rooted tree and finite nonnegative branch lengths are fixed conditions",
            "unobserved_all_zero_sites": "not modeled; no estimate of total possible positions or ancestral site density",
            "missingness": "observation masks are conditioned on; nonrandom annotation, mapping, and discovery biases are not fully corrected",
            "rate_heterogeneity": "no branch-specific or site-specific rate heterogeneity is modeled",
            "independence": "site likelihoods are treated as independent; for linked sites their sum is a composite likelihood",
            "event_interpretation": domain["event_interpretation"],
            "inference_limit": "rates are conditional estimates; no likelihood-ratio tests, confidence intervals, or guaranteed global-optimum claim are reported",
        },
        "families": fit_rows,
    }
    if domain["model"] == MODEL:
        fit_data["assumptions"]["gain_interpretation"] = domain["transition_01"]
    history_data = {
        "schema_version": 1, "model": domain["model"],
        "observation_type": domain["observation_type"], "analysis_scope": domain["analysis_scope"],
        "interpretation": domain["history_interpretation"], "families": history_rows,
    }
    _write_json(fit_path, fit_data)
    _write_json(history_path, history_data)
    artifacts = [tree_path.name, fit_path.name, history_path.name]
    run_data = {
        "schema_version": 1,
        "status": "completed_with_unresolved" if not fit_rows or any(
            row.get("status") not in {"estimated", "fixed_rates"} for row in fit_rows
        ) else "completed",
        "model": domain["model"], "observation_type": domain["observation_type"],
        "analysis_scope": domain["analysis_scope"], "tree_file": tree_path.name,
        "artifacts": artifacts,
        "diagnostics": ["no_candidate_sites"] if not fit_rows else [],
    }
    _write_json(output_dir / "run_result.json", run_data)
    return fit_path, history_path


def write_dna_outputs(output_dir, tree, fit_rows, history_rows):
    """Backward-compatible DNA output function."""
    from intraphy.inference.binary_presence import DNA_DOMAIN
    return write_binary_presence_outputs(output_dir, tree, fit_rows, history_rows, DNA_DOMAIN)


def _write_json(path, data):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)
