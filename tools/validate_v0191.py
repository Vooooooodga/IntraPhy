#!/usr/bin/env python3
"""Keep all raw audit cases, outputs and evaluation metadata. No empirical calibration claim."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
import time
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from intraphy.commands.parser import build_parser
from intraphy.commands.exons import dispatch_analyze
from intraphy.inputs.selection import resolve_inputs
from intraphy.structure.serialization import write_json
from intraphy.commands.environment import environment_report
from intraphy.verification.exon_audit_cases import AUDIT_SCENARIOS, write_audit_example, check_audit_result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--output-dir",required=True,type=Path);args=p.parse_args()
    out=args.output_dir.resolve();out.mkdir(parents=True,exist_ok=True)
    if any(out.iterdir()): raise SystemExit("Validation output must be empty")
    write_json(out/"environment.json",environment_report());results=[]
    for name in AUDIT_SCENARIOS:
        start=time.monotonic()
        try:
            raw=write_audit_example(out/name/"raw",name)
            (raw/"truth.json").replace(out/name/"evaluation_truth.json")
            result=out/name/"result"
            cmd=build_parser().parse_args(["analyze","--model","exon-parsimony",
                "--fasta",str(raw),"--gff",str(raw),
                "--species-tree",str(raw/"species_tree.nwk"),"--output-dir",str(result)])
            cmd._input_selection=resolve_inputs(cmd);summary=dispatch_analyze(cmd)
            details=json.loads((result/"exon_history.json").read_text())
            item={"scenario":name,"status":"passed",**check_audit_result(name,summary,details,result)}
        except (Exception,SystemExit) as exc:
            item={"scenario":name,"status":"failed","error_type":type(exc).__name__,"error":str(exc)}
        item["seconds"]=round(time.monotonic()-start,3);results.append(item);print(json.dumps(item),flush=True)
        write_json(out/"validation.json",{"version":"0.19.1",
            "model_contract":"exon-parsimony (explicit legacy regression)","cases":results,
            "scope":"synthetic_raw_audit_regressions_not_biological_accuracy_or_statistical_calibration"})
    return int(any(r["status"]!="passed" for r in results))


if __name__=="__main__": raise SystemExit(main())
