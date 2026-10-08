#!/usr/bin/env python3
"""Exercise the current default exon-structure CLI on fixed synthetic cases.

This is a software/output-contract regression, not a biological accuracy test.
Unmarked omission is interpreted only relative to supplied annotation paths; it
cannot be distinguished from genuine non-exonic annotation or auto-corrected.
Each fixture's truth record is moved outside every analysis input before inference.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time


def _read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _number(value, field, context):
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise AssertionError(f"{context}: {field} is not numeric: {value!r}") from exc
    if not math.isfinite(parsed) or parsed < -1e-12 or parsed > 1 + 1e-12:
        raise AssertionError(f"{context}: {field} is not a finite probability in [0, 1]: {value!r}")
    return parsed


def _validate_native_result(directory):
    directory = Path(directory)
    required_artifacts = {"exon_structure_fit.json", "exon_history.json",
        "exon_structure_summary.tsv", "ancestral_exon_states.tsv",
        "branch_exon_changes.tsv", "model_diagnostics.json"}
    run_record = json.loads((directory / "run_result.json").read_text(encoding="utf-8"))
    history = json.loads((directory / "exon_history.json").read_text(encoding="utf-8"))
    diagnostics = json.loads((directory / "model_diagnostics.json").read_text(encoding="utf-8"))
    if run_record.get("status") not in {"completed", "completed_with_unresolved"}:
        raise AssertionError(f"Unexpected native completion status: {run_record}")
    if run_record.get("model") != "exon-structure-ctmc" or not required_artifacts.issubset(run_record.get("artifacts", [])):
        raise AssertionError(f"Unexpected native run manifest: {run_record}")
    for artifact in required_artifacts:
        path = directory / artifact
        if not path.is_file() or path.stat().st_size == 0:
            raise AssertionError(f"Native artifact is missing or empty: {path}")
    if history.get("model") != "exon-structure-ctmc" or diagnostics.get("model") != "exon-structure-ctmc":
        raise AssertionError("Native result labels do not identify exon-structure-ctmc")
    if history.get("observation_unit") != "genomic_exon_spans":
        raise AssertionError("Unexpected native observation unit")
    fit = json.loads((directory / "exon_structure_fit.json").read_text(encoding="utf-8"))
    if fit.get("model") != "exon-structure-ctmc" or fit.get("parameter_mode") != "fit":
        raise AssertionError("Native fit record does not identify the default fitted exon-structure model")

    summary = _read_tsv(directory / "exon_structure_summary.tsv")
    branches = _read_tsv(directory / "branch_exon_changes.tsv")
    ancestral = _read_tsv(directory / "ancestral_exon_states.tsv")
    if not summary or not branches or not ancestral:
        raise AssertionError("Native unit, branch, and node-state tables must be nonempty")
    if not {"family_id", "unit_id", "status", "log_likelihood", "unknown_tips"}.issubset(summary[0]):
        raise AssertionError("Native unit summary has an unexpected schema")
    if not {"family_id", "unit_id", "parent_node_id", "child_node_id",
            "probability_exon_structure_change", "probability_dna_presence_change",
            "joint_configuration_probability", "exon_count_pair_probabilities"}.issubset(branches[0]):
        raise AssertionError("Native branch table has an unexpected schema")
    if not {"family_id", "unit_id", "node", "state_id", "probability"}.issubset(ancestral[0]):
        raise AssertionError("Native node table has an unexpected schema")

    for row in summary:
        if row["log_likelihood"] not in {"", "NA", "None"}:
            if not math.isfinite(float(row["log_likelihood"])):
                raise AssertionError("Native summary contains a nonfinite log-likelihood")
        unknown = json.loads(row["unknown_tips"] or "[]")
        if not isinstance(unknown, list):
            raise AssertionError("Native whole-tip unknown_tips field must be a JSON list")

    node_sums = {}
    for row in ancestral:
        key = (row["family_id"], row["unit_id"], row["node"])
        node_sums[key] = node_sums.get(key, 0.0) + _number(row["probability"], "probability", key)
    if any(not math.isclose(total, 1.0, rel_tol=1e-9, abs_tol=1e-12)
           for total in node_sums.values()):
        raise AssertionError(f"Native node-state posterior does not normalize: {node_sums}")

    probability_columns = ("probability_exon_structure_change", "probability_dna_presence_change",
        "joint_configuration_probability", "probability_exon_count_increase",
        "probability_exon_count_decrease", "probability_exon_count_unchanged")
    required_probabilities = ("probability_exon_structure_change",
        "probability_dna_presence_change", "joint_configuration_probability",
        "exon_count_pair_probabilities")
    for row in branches:
        context = (row.get("family_id"), row.get("unit_id"), row.get("parent_node_id"), row.get("child_node_id"))
        if any(row.get(field) in {None, "", "NA", "None"} for field in required_probabilities):
            raise AssertionError(f"{context}: required native branch probabilities are missing")
        for field in probability_columns:
            if row.get(field) not in {None, "", "NA", "None"}:
                _number(row[field], field, context)
        if row.get("exon_count_pair_probabilities") not in {None, "", "NA", "None"}:
            pairs = json.loads(row["exon_count_pair_probabilities"])
            if not isinstance(pairs, list):
                raise AssertionError(f"{context}: exon-count probabilities must be a JSON list")
            pair_total = 0.0
            for pair in pairs:
                pair_total += _number(pair["probability"], "exon_count_pair_probabilities.probability", context)
            if not math.isclose(pair_total, 1.0, rel_tol=1e-9, abs_tol=1e-12):
                raise AssertionError(f"{context}: exon-count-pair probabilities do not normalize")
    return {"run_record": run_record, "summary": summary,
            "branches": branches, "ancestral": ancestral}


def _unique_rows(rows, fields, label):
    indexed = {}
    for row in rows:
        key = tuple(row.get(field, "") for field in fields)
        if key in indexed:
            raise AssertionError(f"{label}: duplicate native row key {key}")
        indexed[key] = row
    return indexed


def _compare_native_results(first, second):
    summary_fields = ("family_id", "unit_id")
    first_summary = _unique_rows(first["summary"], summary_fields, "default summary")
    second_summary = _unique_rows(second["summary"], summary_fields, "reload summary")
    if first_summary.keys() != second_summary.keys():
        raise AssertionError("Reload summary rows differ from the default result")
    for key in first_summary:
        left, right = first_summary[key], second_summary[key]
        if (left["status"], left["unknown_tips"]) != (right["status"], right["unknown_tips"]):
            raise AssertionError(f"Reload summary status differs for {key}")
        left_has_loglik = left["log_likelihood"] not in {"", "NA", "None"}
        right_has_loglik = right["log_likelihood"] not in {"", "NA", "None"}
        if left_has_loglik != right_has_loglik:
            raise AssertionError(f"Reload log-likelihood availability differs for {key}")
        if left_has_loglik:
            if not math.isclose(float(left["log_likelihood"]), float(right["log_likelihood"]),
                                rel_tol=1e-9, abs_tol=1e-12):
                raise AssertionError(f"Reload log-likelihood differs for {key}")

    ancestral_fields = ("family_id", "unit_id", "node", "state_id")
    first_nodes = _unique_rows(first["ancestral"], ancestral_fields, "default ancestral states")
    second_nodes = _unique_rows(second["ancestral"], ancestral_fields, "reload ancestral states")
    if first_nodes.keys() != second_nodes.keys():
        raise AssertionError("Reload ancestral-state rows differ from the default result")
    for key in first_nodes:
        if not math.isclose(float(first_nodes[key]["probability"]), float(second_nodes[key]["probability"]),
                            rel_tol=1e-9, abs_tol=1e-12):
            raise AssertionError(f"Reload ancestral probability differs for {key}")

    branch_fields = ("family_id", "unit_id", "parent_node_id", "child_node_id",
        "parent_exons", "child_exons", "parent_dna_presence", "child_dna_presence",
        "material_ids", "changed_material_tracts")
    probability_fields = ("probability_exon_structure_change", "probability_dna_presence_change",
        "joint_configuration_probability", "probability_exon_count_increase",
        "probability_exon_count_decrease", "probability_exon_count_unchanged")
    first_branches = _unique_rows(first["branches"], branch_fields, "default branch rows")
    second_branches = _unique_rows(second["branches"], branch_fields, "reload branch rows")
    if first_branches.keys() != second_branches.keys():
        raise AssertionError("Reload branch/configuration rows differ from the default result")
    for key in first_branches:
        for field in probability_fields:
            left, right = first_branches[key].get(field), second_branches[key].get(field)
            if left in {None, "", "NA", "None"} or right in {None, "", "NA", "None"}:
                if left != right:
                    raise AssertionError(f"Reload {field} availability differs for {key}")
            elif not math.isclose(float(left), float(right), rel_tol=1e-9, abs_tol=1e-12):
                raise AssertionError(f"Reload {field} differs for {key}")
        left_pairs = { (item["parent_count"], item["child_count"]): float(item["probability"])
                      for item in json.loads(first_branches[key]["exon_count_pair_probabilities"]) }
        right_pairs = { (item["parent_count"], item["child_count"]): float(item["probability"])
                       for item in json.loads(second_branches[key]["exon_count_pair_probabilities"]) }
        if left_pairs.keys() != right_pairs.keys():
            raise AssertionError(f"Reload exon-count-pair rows differ for {key}")
        for pair in left_pairs:
            if not math.isclose(left_pairs[pair], right_pairs[pair], rel_tol=1e-9, abs_tol=1e-12):
                raise AssertionError(f"Reload exon-count-pair probability differs for {key}, {pair}")


def _localized_unknown_observations(catalogues, species):
    rows = [observation for catalogue in catalogues
            for observation in catalogue.get("observations", ())
            if observation.get("species") == species
            and observation.get("kind") in {"partial", "unknown"}
            and observation.get("unknown_intervals")]
    if not rows:
        raise AssertionError(f"No localized unknown observation intervals found for {species}")
    return rows


def _run_raw_fixture(root, scenario, invoke):
    raw = root / f"raw_{scenario}"
    result = root / f"analyze_{scenario}"
    invoke("example-exons", "--scenario", scenario, "--seed", "19",
           "--output-dir", raw)
    truth = raw / "truth.json"
    if not truth.is_file():
        raise AssertionError(f"{scenario} fixture lacks its separate evaluation truth record")
    truth.rename(root / f"{scenario}_evaluation_truth.json")
    if truth.exists():
        raise AssertionError(f"{scenario} truth remained inside analysis inputs")

    # Omit --model to exercise the production default analyze route.
    invoke("analyze", "--fasta", raw, "--gff", raw,
           "--species-tree", raw / "species_tree.nwk", "--output-dir", result)
    catalogues = [json.loads(line) for line in
                  (result / "exon_configurations.jsonl").read_text(encoding="utf-8").splitlines()
                  if line.strip()]
    if not catalogues or any(c.get("observation_unit") != "genomic_exon_spans"
                             or c.get("schema") != "intraphy.exon-configurations/2"
                             for c in catalogues):
        raise AssertionError(f"{scenario} catalogue has an unexpected schema or observation unit")
    return result, catalogues, _validate_native_result(result)


def _validate_annotation_dropout(catalogues, species):
    supported = [
        (catalogue, observation)
        for catalogue in catalogues
        for observation in catalogue.get("observations", ())
        if observation.get("species") == species
        and "annotation_supported_nonexonic_interval" in observation.get("reasons", ())
    ]
    if not supported:
        raise AssertionError(f"No annotation-conditioned omission observation found for {species}")
    for catalogue, observation in supported:
        if observation.get("kind") != "observed":
            raise AssertionError(f"{species} omission was not observed under the supplied annotation")
        candidate_spans = catalogue.get("spans", ())
        configurations = observation.get("configurations", ())
        if not candidate_spans or not configurations:
            raise AssertionError(f"{species} omission lacks candidate spans or observed configurations")
        for configuration in configurations:
            for candidate in candidate_spans:
                if any(candidate["start"] < exon["end"] and exon["start"] < candidate["end"]
                       for exon in configuration.get("exons", ())):
                    raise AssertionError(f"{species} omitted candidate span appears in its exon configuration")
    return supported


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--isolated", action="store_true",
                        help="Run child CLI commands with Python isolated mode (installed-wheel smoke).")
    args = parser.parse_args()
    root = args.output_dir.resolve()
    root.mkdir(parents=True, exist_ok=True)
    if any(root.iterdir()):
        raise SystemExit("Validation output must be empty")
    reload_result = root / "explicit_model_reload"
    commands = []
    direct_env = {key: value for key, value in os.environ.items()
                  if "proxy" not in key.casefold()}

    def invoke(*argv):
        command = [sys.executable, *( ["-I"] if args.isolated else [] ), "-m", "intraphy",
                    *map(str, argv)]
        started = time.monotonic()
        completed = subprocess.run(command, cwd=root, env=direct_env,
                                   capture_output=True, text=True, timeout=900)
        item = {"argv": command, "returncode": completed.returncode,
                "seconds": round(time.monotonic() - started, 3)}
        commands.append(item)
        index = len(commands)
        (root / f"command-{index:02d}.stdout").write_text(completed.stdout, encoding="utf-8")
        (root / f"command-{index:02d}.stderr").write_text(completed.stderr, encoding="utf-8")
        (root / "commands.json").write_text(json.dumps(commands, indent=2) + "\n", encoding="utf-8")
        if completed.returncode:
            raise RuntimeError(f"CLI command failed ({completed.returncode}): {command}\n{completed.stderr}")

    default_result, _, default_native = _run_raw_fixture(root, "split_insertion", invoke)
    _, dropout_catalogues, _ = _run_raw_fixture(root, "annotation_dropout", invoke)
    _validate_annotation_dropout(dropout_catalogues, "Species_D")
    _, gap_catalogues, _ = _run_raw_fixture(root, "assembly_gap", invoke)
    _localized_unknown_observations(gap_catalogues, "Species_D")

    # Re-enter through analyze with the explicit current model and native catalogue.
    invoke("analyze", "--model", "exon-structure-ctmc", "--input-dir",
           default_result / "prepared_inputs", "--exon-configurations",
           default_result / "exon_configurations.jsonl", "--output-dir", reload_result)
    reload_catalogues = [json.loads(line) for line in
        (reload_result / "exon_configurations.jsonl").read_text(encoding="utf-8").splitlines() if line.strip()]
    if not reload_catalogues or any(c.get("observation_unit") != "genomic_exon_spans"
                                    or c.get("schema") != "intraphy.exon-configurations/2"
                                    for c in reload_catalogues):
        raise AssertionError("Reload catalogue has an unexpected schema or observation unit")
    reload_native = _validate_native_result(reload_result)
    _compare_native_results(default_native, reload_native)

    record = {"status": "passed", "model": "exon-structure-ctmc",
              "route": "default analyze plus explicit-model analyze reload",
              "fixtures": ["split_insertion", "annotation_dropout", "assembly_gap"],
              "truth_handling": "each fixture truth record moved outside raw inputs before analyze",
              "annotation_dropout_scope": "unmarked omission cannot be distinguished from genuine non-exonic annotation under supplied paths; no automatic dropout correction or accuracy claim",
              "unknown_scope_checked": "assembly-gap localized unknown intervals are distinct from whole-tip unknown_tips",
              "artifact_contract": "nonempty native summary/branch/node tables; finite bounded probabilities; normalized node posteriors; reload agreement",
              "biological_accuracy_benchmark": False,
              "statistical_calibration": False,
              "isolated_cli": bool(args.isolated), "commands": len(commands)}
    (root / "validation.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(record, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
