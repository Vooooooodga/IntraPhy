#!/usr/bin/env python3
"""Smoke the installed wheel outside the checkout without mixing model outputs."""
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
import xml.etree.ElementTree as ET


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    root = Path(args.output_dir).resolve()
    root.mkdir(parents=True, exist_ok=True)
    if any(root.iterdir()):
        raise SystemExit("Smoke output must be empty")

    commands = []
    direct_env = {key: value for key, value in os.environ.items()
                  if "proxy" not in key.casefold()}

    def invoke(*argv):
        started = time.monotonic()
        command = [sys.executable, "-I", "-m", "intraphy", *map(str, argv)]
        result = subprocess.run(command, cwd=root, env=direct_env,
                                capture_output=True, text=True, timeout=900)
        index = len(commands) + 1
        (root / f"command-{index:02d}.stdout").write_text(result.stdout, encoding="utf-8")
        (root / f"command-{index:02d}.stderr").write_text(result.stderr, encoding="utf-8")
        commands.append({"argv": command, "returncode": result.returncode,
                         "seconds": round(time.monotonic() - started, 3)})
        (root / "commands.json").write_text(json.dumps(commands, indent=2) + "\n", encoding="utf-8")
        if result.returncode:
            raise RuntimeError(f"Installed command failed: {command}\n{result.stderr}")

    invoke("--version")
    invoke("inspect-aligners")

    # Reuse the same default-analyze CLI regression as source/CI validation.
    regression = Path(__file__).resolve().with_name("validate_exon_structure.py")
    started = time.monotonic()
    command = [sys.executable, "-I", str(regression), "--isolated",
               "--output-dir", str(root / "current_model")]
    result = subprocess.run(command, cwd=root, env=direct_env,
                            capture_output=True, text=True, timeout=1800)
    index = len(commands) + 1
    (root / f"command-{index:02d}.stdout").write_text(result.stdout, encoding="utf-8")
    (root / f"command-{index:02d}.stderr").write_text(result.stderr, encoding="utf-8")
    commands.append({"argv": command, "returncode": result.returncode,
                     "seconds": round(time.monotonic() - started, 3)})
    (root / "commands.json").write_text(json.dumps(commands, indent=2) + "\n", encoding="utf-8")
    if result.returncode:
        raise RuntimeError(f"Installed current-model regression failed: {command}\n{result.stderr}")

    # Renderer compatibility is checked only with a separate explicit legacy result.
    invoke("example-exons", "--scenario", "split_insertion", "--output-dir", "legacy_raw")
    (root / "legacy_raw/truth.json").unlink(missing_ok=True)
    raw = ("--fasta", "legacy_raw", "--gff", "legacy_raw",
           "--species-tree", "legacy_raw/species_tree.nwk")
    invoke("check", *raw)
    invoke("extract-loci", *raw, "--flank", "50", "--output-dir", "portable")
    invoke("analyze", "--model", "exon-parsimony", *raw, "--output-dir", "legacy_parsimony")
    invoke("infer-phylogeny", "--model", "exon-parsimony", "--input-dir",
           "legacy_parsimony/prepared_inputs", "--exon-configurations",
           "legacy_parsimony/exon_configurations.jsonl", "--output-dir", "legacy_reload")
    invoke("exon-rate-template", "--output", "legacy_rates.json", "--rate", ".1")
    invoke("infer-phylogeny", "--model", "exon-ctmc", "--exon-rates", "legacy_rates.json",
           "--expected-edits", "--input-dir",
           "legacy_parsimony/prepared_inputs", "--exon-configurations",
           "legacy_parsimony/exon_configurations.jsonl", "--output-dir", "legacy_ctmc")
    invoke("visualize", "--input-dir", "legacy_parsimony/prepared_inputs",
           "--result-dir", "legacy_parsimony", "--output-dir", "legacy_figures")
    invoke("explain", "--output-dir", "guide")

    parsimony_record = json.loads((root / "legacy_parsimony/run_result.json").read_text(encoding="utf-8"))
    assert parsimony_record["model"] == "exon-parsimony"
    parsimony_history = json.loads((root / "legacy_parsimony/exon_history.json").read_text(encoding="utf-8"))
    assert parsimony_history["model"] == "exon_configuration_v2"
    assert parsimony_history["engine"] == "exon-parsimony"
    required = [event for unit in parsimony_history["units"]
                for event in unit.get("views", {}).get("evidence", {}).get("events", ())
                if event.get("support") == "required"]
    assert len(required) == 1 and required[0]["operation"] == "dna_insertion", required
    assert "exon_split" in required[0]["consequences"], required

    ctmc_record = json.loads((root / "legacy_ctmc/run_result.json").read_text(encoding="utf-8"))
    assert ctmc_record["model"] == "exon-ctmc"
    ctmc_diagnostics = json.loads((root / "legacy_ctmc/model_diagnostics.json").read_text(encoding="utf-8"))
    assert ctmc_diagnostics["units"] and all(row.get("probability_status") ==
               "conditional_on_fixed_parameters_and_declared_catalogue"
               for row in ctmc_diagnostics["units"]), ctmc_diagnostics
    with (root / "legacy_ctmc/exon_branch_posteriors.tsv").open(encoding="utf-8", newline="") as handle:
        ctmc_branches = list(csv.DictReader(handle, delimiter="\t"))
    if not ctmc_branches or "expected_edits" not in ctmc_branches[0]:
        raise AssertionError("Legacy expected-edits smoke did not emit branch expected counts")
    for row in ctmc_branches:
        value = float(row["expected_edits"])
        if not math.isfinite(value) or value < 0:
            raise AssertionError(f"Invalid legacy expected-edits value: {row}")

    figure_groups = {"legacy_figures": list((root / "legacy_figures").rglob("*.svg")),
                     "guide": list((root / "guide").rglob("*.svg"))}
    for label, paths in figure_groups.items():
        if not paths:
            raise AssertionError(f"{label} did not produce its own SVG output")
        for path in paths:
            ET.parse(path)
    svgs = [path for paths in figure_groups.values() for path in paths]
    summary = {"status": "passed", "commands": len(commands), "svg_files": len(svgs),
               "current_model_cli_regression": "shared validate_exon_structure.py under isolated installed interpreter",
               "visualization_model": "exon-parsimony (separate explicit legacy result)",
               "ctmc_model": "exon-ctmc (fixed rate file)",
               "legacy_contracts": ["check", "extract-loci", "required DNA insertion with exon_split consequence",
                   "fixed-parameter CTMC conditional diagnostics and expected edits", "separate renderer/guide SVGs"],
               "validation_scope": "installed_wheel_cli_contracts_not_biological_accuracy_or_calibration"}
    (root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
