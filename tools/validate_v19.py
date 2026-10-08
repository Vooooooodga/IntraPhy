#!/usr/bin/env python3
"""Retained raw FASTA/GFF validation, including observation damage (not calibration)."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
import time

# Executing this developer utility intentionally validates the checkout.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"src"))
from intraphy.commands.parser import build_parser
from intraphy.commands.exons import dispatch_analyze
from intraphy.commands.environment import environment_report
from intraphy.inputs.selection import resolve_inputs
from intraphy.verification.exon_cases import SCENARIOS, write_exon_example
from intraphy.structure.serialization import write_json


def expected_check(name, rows, details):
    total = sum(r.get("minimum_structural_edits") or 0 for r in rows)
    events = [e for u in details["units"] for e in u.get("views", {}).get("evidence", {}).get("events", [])]
    required = [e for e in events if e["support"] == "required"]
    positive = {"split_insertion", "fusion_phase0", "fusion_phase1", "fusion_phase2",
                "exon_deletion", "multi_exon_deletion", "negative_strand"}
    if name in positive:
        assert total == 1, (name, total)
        assert len(required) == 1, (name, required)
    elif name in {"duplication", "inversion"}:
        assert any(r["status"] == "unresolved" for r in rows), (name, rows)
    else:
        assert total == 0 and not required, (name, total, required)
    if name == "intronization":
        annotation = [e for u in details["units"] for e in u.get("views", {}).get("annotation", {}).get("events", [])]
        assert any(e["operation"] == "split" and e["support"] == "required" for e in annotation)
    if name == "multi_exon_deletion":
        assert "deleted_exons:2" in required[0]["consequences"]
    if name.startswith("fusion_phase"):
        assert required[0]["operation"] == "dna_deletion" and "exon_fusion" in required[0]["consequences"]
    return {"minimum_edits_in_resolved_units": total, "required_edits": len(required),
            "unresolved_units": sum(r["status"] == "unresolved" for r in rows)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--scenario", action="append", choices=SCENARIOS)
    args = parser.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if any(out.iterdir()):
        raise SystemExit("Validation output must be empty")
    write_json(out/"environment.json", environment_report())
    results = []
    for scenario in args.scenario or SCENARIOS:
        begin = time.monotonic()
        try:
            raw = write_exon_example(out/scenario/"raw", scenario)
            # Retain truth OUTSIDE the input, never as an analysis parameter.
            (raw/"truth.json").replace(out/scenario/"evaluation_truth.json")
            command = build_parser().parse_args(["analyze", "--model", "exon-parsimony",
                "--fasta", str(raw), "--gff", str(raw),
                "--species-tree", str(raw/"species_tree.nwk"), "--output-dir", str(out/scenario/"result")])
            command._input_selection = resolve_inputs(command)
            rows = dispatch_analyze(command)
            details = json.loads((out/scenario/"result/exon_history.json").read_text())
            counts = expected_check(scenario, rows, details)
            item = {"scenario": scenario, "status": "passed", **counts}
        except (Exception, SystemExit) as exc:
            item = {"scenario": scenario, "status": "failed", "error_type": type(exc).__name__, "error": str(exc)}
        item["seconds"] = round(time.monotonic()-begin, 3)
        results.append(item)
        print(json.dumps(item), flush=True)
        write_json(out/"validation.json", {"model_contract": "exon-parsimony (explicit legacy regression)",
            "cases": results, "scope": "synthetic_raw_input_contracts",
            "biological_accuracy_benchmark": False, "statistical_calibration": False})
    return int(any(r["status"] != "passed" for r in results))


if __name__ == "__main__":
    raise SystemExit(main())
