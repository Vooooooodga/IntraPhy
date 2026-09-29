"""CLI surface for genomic DNA-copy and legacy exon models."""
from __future__ import annotations
import argparse
from pathlib import Path
from .file_inputs import add_file_inputs
from ..verification.exon_cases import SCENARIOS

MODELS = ("exon-locus-ctmc", "exon-parsimony", "exon-ctmc")
EXON_STRUCTURE_MODEL = "exon-structure-ctmc"
DNA_PRESENCE_MODEL = "dna-presence-ctmc"
INTRON_POSITION_MODEL = "intron-position-ctmc"


class _LocusModelAction(argparse.Action):
    """Keep the historical explicit-model parse behavior for /2 inputs."""

    def __call__(self, parser, namespace, value, option_string=None):
        setattr(namespace, self.dest, value)
        if not getattr(namespace, "_model_explicit", False):
            namespace.model = "exon-locus-ctmc"
        if not getattr(namespace, "_root_frequency_explicit", False):
            namespace.root_frequency = "estimated"


class _AnalyzeModelAction(argparse.Action):
    def __call__(self, parser, namespace, value, option_string=None):
        setattr(namespace, self.dest, value)
        namespace._model_explicit = True
        if not getattr(namespace, "_root_frequency_explicit", False):
            namespace.root_frequency = "stationary" if value in {DNA_PRESENCE_MODEL, INTRON_POSITION_MODEL} else (None if value == EXON_STRUCTURE_MODEL else "estimated")


class _RootFrequencyAction(argparse.Action):
    def __call__(self, parser, namespace, value, option_string=None):
        setattr(namespace, self.dest, value)
        namespace._root_frequency_explicit = True


class _RootPresenceAction(argparse.Action):
    def __call__(self, parser, namespace, value, option_string=None):
        setattr(namespace, self.dest, value)
        namespace._root_presence_explicit = True


class _TrackedValueAction(argparse.Action):
    def __call__(self, parser, namespace, value, option_string=None):
        setattr(namespace, self.dest, value)
        setattr(namespace, f"_{self.dest}_explicit", True)


def add_dna_observation_options(command):
    command.add_argument("--survey-min-identity", type=float, default=.7, action=_TrackedValueAction,
                         help="Minimum identity for a qualified DNA position-homology observation.")
    command.add_argument("--survey-min-coverage", type=float, default=.8, action=_TrackedValueAction,
                         help="Minimum aligned coverage for a qualified DNA position-homology observation.")
    command.add_argument("--survey-max-dp-cells", type=int, default=250000, action=_TrackedValueAction,
                         help="Pairwise alignment cell budget; over-budget candidates remain unknown.")


def add_intron_observation_options(command):
    command.add_argument("--intron-anchor-window", type=int, default=15, action=_TrackedValueAction,
                         help="Maximum aligned amino-acid columns on each flank used to identify corresponding intron positions.")
    command.add_argument("--intron-min-anchor-pairs", type=int, default=8, action=_TrackedValueAction,
                         help="Minimum paired alignment columns required to qualify an intron position.")


def add_configuration_options(command):
    command.add_argument("--exon-configurations", help="Explicit exon-configuration JSONL catalogue.")
    command.add_argument("--exon-rates", help="Explicit rate JSON; exon-structure-ctmc requires --parameter-mode fixed.")
    command.add_argument("--observation-view", choices=("evidence", "annotation"), default="evidence",
                         action=_TrackedValueAction)
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
    analyze = sub.add_parser("analyze", help="Infer exon-structure histories from genomic exon spans, or select another supported model.")
    add_file_inputs(analyze, required=False)
    analyze.add_argument("--input-dir", help="Prepared genomic case directory containing species_tree.tsv.")
    analyze.add_argument("--genomic-evidence-dir", help="Reuse staged observations from prepare-genomic-evidence for the selected character type.")
    analyze.add_argument("--locus-model", action=_LocusModelAction,
                         help="Advanced DNA-only labelled-copy JSON model (intraphy.exon-locus-model/2).")
    analyze.add_argument("--species-tree", help="Rooted Newick or TSV tree for raw inputs or --locus-model; with --input-dir, replace its tree only with the same tip panel.")
    analyze.add_argument("--parameter-mode", choices=("fit", "fixed"), default="fit",
                         help="Fit the default exon-structure-ctmc rate, or evaluate an explicit --exon-rates file in fixed mode.")
    analyze.add_argument("--output-dir", required=True)
    analyze.add_argument("--threads", type=int, default=1)
    analyze.add_argument("--flank", type=int, default=1000, action=_TrackedValueAction)
    analyze.add_argument("--max-extension", type=int, default=10000, action=_TrackedValueAction)
    analyze.add_argument("--model", choices=(DNA_PRESENCE_MODEL, INTRON_POSITION_MODEL, EXON_STRUCTURE_MODEL, *MODELS), default=EXON_STRUCTURE_MODEL,
                         action=_AnalyzeModelAction,
                         help="Genomic exon-structure CTMC is the default; DNA presence, intron position, and supplied locus models are separate explicit routes.")
    analyze.add_argument("--root-frequency", choices=("stationary", "fixed", "estimated"), default=None,
                         action=_RootFrequencyAction,
                         help="Root presence distribution for the selected binary genomic character model.")
    analyze.add_argument("--root-presence", type=float, default=.5, action=_RootPresenceAction,
                         help="Root presence probability when --root-frequency fixed is selected.")
    analyze.add_argument("--gain-rate", "--dna-gain-rate", dest="dna_gain_rate", type=float,
                         help="Fixed gain rate for the selected binary character; supply with --loss-rate. --dna-gain-rate remains an alias.")
    analyze.add_argument("--loss-rate", "--dna-loss-rate", dest="dna_loss_rate", type=float,
                         help="Fixed loss rate for the selected binary character; supply with --gain-rate. --dna-loss-rate remains an alias.")
    analyze.add_argument("--branch-length-mode", choices=("supplied", "unit"), default="supplied")
    add_configuration_options(analyze)
    add_dna_observation_options(analyze)
    add_intron_observation_options(analyze)
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
