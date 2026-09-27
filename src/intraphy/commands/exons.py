"""CLI surface for genomic DNA-copy and legacy exon models."""
from __future__ import annotations
from pathlib import Path
from .file_inputs import add_file_inputs
from ..verification.exon_cases import SCENARIOS

MODELS = ("exon-locus-ctmc", "exon-parsimony", "exon-ctmc")


def add_configuration_options(command):
    command.add_argument("--exon-configurations", help="Explicit exon-configuration JSONL catalogue.")
    command.add_argument("--exon-rates", help="Explicit fixed-rate JSON for exon-ctmc; no invented default estimates.")
    command.add_argument("--observation-view", choices=("evidence", "annotation"), default="evidence")
    command.add_argument("--max-states", type=int, default=None, help="Optional finite candidate-state limit; exceeding it stops inference.")
    command.add_argument("--max-origin-scenarios", type=int, default=None,
                         help="Optional origin-scenario limit; exceeding it stops inference.")
    command.add_argument("--max-observation-scenarios", type=int, default=None,
                         help="Optional positive scenario limit; default enumerates all scenarios lazily.")
    command.add_argument("--max-locus-bases", type=int, default=100000, help="Genomic MSA budget; sequence is not silently cropped.")
    command.add_argument("--alignment-timeout", type=int, default=600)
    command.add_argument("--exon-identity", type=float, default=.7, help="Explicit, uncalibrated nucleotide evidence threshold.")
    command.add_argument("--anchor-bases", type=int, default=12)
    command.add_argument("--anchor-identity", type=float, default=.8)
    command.add_argument("--origin-root-sensitivity", nargs="+", type=float, default=[],
                         help="Explicit alternative root-opportunity weights; conditional sensitivity, not model selection.")
    command.add_argument("--expected-edits", action="store_true",
                         help="Compute marked CTMC branch counts; can be substantially slower.")


def add_exon_commands(sub):
    from .exon_statistics import add_statistics_commands
    add_statistics_commands(sub)
    analyze = sub.add_parser("analyze", help="Fit supplied genomic DNA-copy histories or prepare raw genomic inputs.")
    add_file_inputs(analyze, required=False)
    from .locus import add_locus_options
    add_locus_options(analyze)
    analyze.add_argument("--output-dir", required=True)
    analyze.add_argument("--threads", type=int, default=1)
    analyze.add_argument("--flank", type=int, default=1000)
    analyze.add_argument("--max-extension", type=int, default=10000)
    analyze.add_argument("--model", choices=MODELS, default="exon-locus-ctmc")
    analyze.add_argument("--branch-length-mode", choices=("supplied", "unit"), default="supplied")
    add_configuration_options(analyze)
    cesar = sub.add_parser("realign-exons", help="Optional CESAR2 coding gene-mode prediction, separate from observations.")
    cesar.add_argument("--input-dir", required=True)
    cesar.add_argument("--output-dir", required=True)
    cesar.add_argument("--family-id", required=True)
    cesar.add_argument("--reference-species", required=True)
    cesar.add_argument("--query-species", required=True)
    cesar.add_argument("--reference-transcript")
    cesar.add_argument("--cesar", default="cesar")
    cesar.add_argument("--profile-dir", required=True, help="Explicit first/last/acceptor/donor profiles, with no clade defaults.")
    cesar.add_argument("--codon-matrix", required=True)
    cesar.add_argument("--timeout", type=int, default=600)
    example = sub.add_parser("example-exons", help="Synthetic structural and observation-error test cases.")
    example.add_argument("--output-dir", required=True)
    example.add_argument("--scenario", choices=SCENARIOS, default="split_insertion")
    example.add_argument("--seed", type=int, default=19)
    rate = sub.add_parser("exon-rate-template", help="Write explicit illustrative rate parameters, not estimates.")
    rate.add_argument("--output", required=True)
    rate.add_argument("--rate", type=float, default=.1)


def configuration_arguments(args):
    return {"model": args.model, "configurations": args.exon_configurations, "rates": args.exon_rates,
            "observation_view": args.observation_view, "max_states": args.max_states,
            "max_origins": args.max_origin_scenarios, "max_observation_scenarios": args.max_observation_scenarios,
            "branch_length_mode": args.branch_length_mode, "expected_edits": args.expected_edits,
            "origin_root_sensitivity": args.origin_root_sensitivity,
            "alignment_timeout": args.alignment_timeout, "max_locus_bases": args.max_locus_bases,
            "exon_identity": args.exon_identity, "anchor_bases": args.anchor_bases, "anchor_identity": args.anchor_identity}


def dispatch_analyze(args):
    from ..case import build_case
    from ..inference.configuration_run import infer_configurations
    out = Path(args.output_dir)
    prepared = out/"prepared_inputs"
    args._input_selection.write(prepared)
    build_case(None, prepared, species_tree=args.species_tree, target_rows=args._input_selection.rows,
               flank=args.flank, max_extension=args.max_extension, threads=args.threads, allow_unannotated=True)
    return infer_configurations(prepared, out, **configuration_arguments(args))
