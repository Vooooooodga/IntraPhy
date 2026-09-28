"""Stable argparse interface, separate from command execution."""
import argparse
from intraphy import __version__
from intraphy.aligners.types import MAX_INTERNAL_DP_CELLS
from intraphy.verification.calibration import DEFAULT_SCENARIOS


def _positive_int_or_unlimited(value):
    if value == "unlimited":
        return None
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a positive integer or 'unlimited'") from exc
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be a positive integer or 'unlimited'")
    return parsed


def build_parser():
    parser = argparse.ArgumentParser(
        prog="intraphy", description="Infer genomic DNA-presence or intron-position histories on a supplied rooted species tree.",
        epilog="Primary workflow: analyze --fasta GENOMES --gff ANNOTATIONS --species-tree TREE. Advanced labelled-copy inference uses --locus-model. See docs/quickstart.md.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    from .exons import add_exon_commands, add_configuration_options
    add_exon_commands(sub)
    from .exons import add_dna_observation_options
    prepare_dna = sub.add_parser(
        "prepare-genomic-evidence",
        help="Derive qualified DNA-presence or intron-position observations from prepared genomic inputs.",
    )
    prepare_dna.add_argument("--character-type", choices=("dna-presence", "intron-position"), default="dna-presence",
                             help="Observation domain to derive; staged tables remain separate by character type.")
    prepare_dna.add_argument("--input-dir", required=True,
                             help="Prepared build-case directory, optionally already containing derive-tables outputs.")
    prepare_dna.add_argument("--output-dir", required=True)
    prepare_dna.add_argument("--threads", type=int, default=1)
    add_dna_observation_options(prepare_dna)
    from .exons import add_intron_observation_options
    add_intron_observation_options(prepare_dna)
    from .repertoires import add_repertoire_commands
    add_repertoire_commands(sub)
    from .locus_statistics import add_locus_statistics_commands
    add_locus_statistics_commands(sub)

    example = sub.add_parser("example", help="Write a small synthetic raw FASTA/GFF3 dataset.")
    example.add_argument("--output-dir", required=True)
    example.add_argument("--scenario", choices=["conserved", "splice_difference", "annotation_dropout"],
                         default="splice_difference")
    example.add_argument("--seed", type=int, default=18)

    explain = sub.add_parser("explain", help="Render the synthetic structure/tree methods guide (no input data).")
    explain.add_argument("--output-dir", required=True)
    explain.add_argument("--legacy-v18", action="store_true", help="Render the previous three-layer teaching guide, explicitly as legacy.")
    explain.add_argument("--png", action="store_true", help="Also render PNG previews using optional CairoSVG.")

    extract = sub.add_parser("extract-gene")
    extract.add_argument("--genome", required=True)
    extract.add_argument("--annotation", required=True)
    extract.add_argument("--gene-id", required=True)
    extract.add_argument("--family-id", required=True)
    extract.add_argument("--species", required=True)
    extract.add_argument("--gene-copy-id", required=True)
    extract.add_argument("--output-dir", required=True)
    extract.add_argument("--append", action="store_true")
    extract.add_argument(
        "--transcript-policy",
        choices=["canonical", "all"],
        default="all",
        help="Transcript structures to extract; default keeps all annotated transcript paths.",
    )
    extract.add_argument("--canonical-rule", choices=["longest_cds", "longest_span"], default="longest_cds")
    extract.add_argument("--source-label", default="unknown_source")
    extract.add_argument("--copy-role", choices=["source", "background", "derived", "candidate"], default="candidate")
    extract.add_argument("--flank", type=int, default=1000)
    extract.add_argument("--max-extension", type=int, default=10000)

    derive = sub.add_parser("derive-tables")
    derive.add_argument("--input-dir", required=True)
    derive.add_argument("--output-dir")
    derive.add_argument("--identity-threshold", type=float, default=0.7)
    derive.add_argument("--distance-table")
    derive.add_argument("--aligner", choices=["auto", "internal", "mafft", "minimap2", "lastz"], default="mafft", help="Exon-pair backend: mafft/auto uses overlap projection; others use local alignment.")
    derive.add_argument("--context-aligner", choices=["internal", "minimap2", "lastz"], default="minimap2", help="Local backend for pairs involving non-exon sequence.")
    derive.add_argument(
        "--coding-msa-mode",
        choices=["linsi", "einsi"],
        default="linsi",
        help="MAFFT strategy for the family-level protein alignment.",
    )
    derive.add_argument(
        "--short-context-max-length",
        type=int,
        default=300,
        help="Maximum anchor-bounded nucleotide interval length for enumerating short-alignment candidates.",
    )
    derive.add_argument(
        "--short-alignment-max-dp-cells",
        type=_positive_int_or_unlimited,
        default=MAX_INTERNAL_DP_CELLS,
        metavar="N|unlimited",
        help=(
            f"Maximum DP cells for one complete short alignment (default: {MAX_INTERNAL_DP_CELLS}); "
            "exceeding the budget fails the alignment without truncation. Use 'unlimited' "
            "to remove this DP guard; actual system resources still limit the alignment."
        ),
    )
    derive.add_argument("--threads", type=int, default=1)
    derive.add_argument("--min-size-ratio", type=float, default=0.25)

    sim = sub.add_parser(
        "simulate",
        description="Legacy experimental-multicopy simulator; it does not validate formal single-copy statistics.",
    )
    sim.add_argument("--output-dir", required=True)
    sim.add_argument("--seed", type=int, default=7)
    sim.add_argument(
        "--scenario",
        choices=[
            "compound",
            "exonization",
            "source_join",
            "tandem_duplication",
            "processed_copy_or_intron_loss",
            "segment_split",
            "segment_fusion",
            "segment_split_fusion",
            "splice_boundary_shift",
            "te_exonization",
            "gene_conversion",
            "negative_control",
            "annotation_dropout",
        ],
        default="compound",
    )

    bench = sub.add_parser("benchmark")
    bench.add_argument("--input-dir", required=True)
    bench.add_argument("--output-dir", required=True)

    cal = sub.add_parser(
        "calibrate",
        description="Legacy experimental-multicopy calibration wrapper; it cannot validate formal single-copy CTMC results.",
    )
    cal.add_argument("--output-dir", required=True)
    cal.add_argument("--scenario", action="append", choices=DEFAULT_SCENARIOS)
    cal.add_argument("--replicates", type=int, default=10)
    cal.add_argument("--bootstrap-replicates", type=int, default=0)
    cal.add_argument("--stochastic-maps", type=int, default=0)
    cal.add_argument("--seed", type=int, default=101)

    viz = sub.add_parser("visualize")
    viz.add_argument("--layout", choices=["target-groups", "legacy-overview"], default="target-groups",
                     help="One target per group on unchanged native annotation (default).")
    viz.add_argument("--target", action="append", default=[],
                     help="Exact family/layer/site ID (repeatable); unambiguous site IDs also accepted.")
    viz.add_argument("--target-manifest", help="TSV of targets and optional exact native intervals; see docs/target_views.md.")
    viz.add_argument("--input-dir", required=True)
    viz.add_argument("--result-dir", required=True)
    viz.add_argument("--output-dir", required=True)
    viz.add_argument("--correspondence-encoding", choices=["pattern", "color"], default=None)

    sub.add_parser("inspect-aligners")

    inspect = sub.add_parser("inspect-annotation")
    inspect.add_argument("--annotation", required=True)
    inspect.add_argument("--output-dir", required=True)
    inspect.add_argument("--query", action="append")
    inspect.add_argument("--alias-file")
    inspect.add_argument("--species", default="NA")
    inspect.add_argument("--case-id", default="case")

    case = sub.add_parser("build-case")
    from .file_inputs import add_file_inputs
    add_file_inputs(case)
    case.add_argument("--output-dir", required=True)
    case.add_argument("--identity-threshold", type=float, default=0.7)
    case.add_argument("--species-tree")
    case.add_argument("--copy-tree")
    case.add_argument("--gene-tree")
    case.add_argument(
        "--transcript-policy",
        choices=["canonical", "all"],
        default="all",
        help="Transcript structures to extract; default keeps all annotated transcript paths.",
    )
    case.add_argument("--canonical-rule", choices=["longest_cds", "longest_span"], default="longest_cds")
    case.add_argument("--aligner", choices=["auto", "internal", "mafft", "minimap2", "lastz"], default="mafft", help="Exon-pair backend: mafft/auto uses overlap projection; others use local alignment.")
    case.add_argument("--context-aligner", choices=["internal", "minimap2", "lastz"], default="minimap2", help="Local backend for pairs involving non-exon sequence.")
    case.add_argument(
        "--coding-msa-mode",
        choices=["linsi", "einsi"],
        default="linsi",
        help="MAFFT strategy for the family-level protein alignment.",
    )
    case.add_argument(
        "--short-context-max-length",
        type=int,
        default=300,
        help="Maximum anchor-bounded nucleotide interval length for enumerating short-alignment candidates.",
    )
    case.add_argument("--threads", type=int, default=1)
    case.add_argument("--min-size-ratio", type=float, default=0.25)
    case.add_argument("--flank", type=int, default=1000)
    case.add_argument("--max-extension", type=int, default=10000)

    orthofinder = sub.add_parser(
        "import-orthofinder",
        description=(
            "Import one upstream orthogroup: every selected species must map to exactly "
            "one annotated gene locus. Multiple isoform members of that locus are accepted "
            "and their source IDs retained. Ambiguous or unresolved mappings error by "
            "default; --on-unresolved exclude records and omits those species."
        ),
    )
    orthofinder.add_argument(
        "--orthofinder-dir", required=True,
        help="Completed results directory containing Orthogroups.tsv or Orthogroups.txt; uses WorkingDirectory/SequenceIDs.txt when available.",
    )
    orthofinder.add_argument("--orthogroup", required=True)
    orthofinder.add_argument(
        "--genome-manifest", required=True,
        help="TSV selecting species and genome_fasta/annotation_file resources; may declare an exact member_id_prefix per species.",
    )
    orthofinder.add_argument(
        "--species-tree",
        help="Supplied TSV or Newick tree; the full tree is retained unless --prune-species-tree is set.",
    )
    orthofinder.add_argument(
        "--on-unresolved",
        choices=["error", "exclude"],
        default="error",
        help="Action for species whose members do not resolve to one annotated locus (default: error).",
    )
    orthofinder.add_argument(
        "--prune-species-tree",
        action="store_true",
        help="With --species-tree, retain only tips with uniquely mapped loci (default keeps the full tree).",
    )
    orthofinder.add_argument("--output-dir", required=True)

    hidden = sub.add_parser("scan-hidden-segments")
    hidden.add_argument("--source-fasta", required=True)
    hidden.add_argument("--target-fasta", required=True)
    hidden.add_argument("--output-dir", required=True)
    hidden.add_argument("--family-id", default="NA")
    hidden.add_argument("--species", default="NA")
    hidden.add_argument("--gene-copy-id", default="NA")
    hidden.add_argument("--min-identity", type=float, default=0.75)
    hidden.add_argument("--min-coverage", type=float, default=0.5)
    hidden.add_argument("--aligner", choices=["internal", "minimap2", "miniprot"], default="internal")
    hidden.add_argument("--threads", type=int, default=1)

    sensitivity = sub.add_parser(
        "candidate-sensitivity",
        description=(
            "Re-evaluate saved nucleotide candidates while retaining each species "
            "pair's original acceptance threshold."
        ),
    )
    sensitivity.add_argument("--segment-matches", required=True)
    sensitivity.add_argument("--output", required=True)
    sensitivity.add_argument(
        "--identity",
        action="append",
        type=float,
        required=True,
        help="Short-fragment minimum nucleotide identity; repeat for a threshold grid.",
    )
    sensitivity.add_argument(
        "--coverage",
        action="append",
        type=float,
        required=True,
        help="Short-fragment minimum query coverage; repeat for a threshold grid.",
    )

    for name in ["complete-annotation", "segment-correspondence", "compare-baselines"]:
        cmd = sub.add_parser(name)
        cmd.add_argument("--input-dir", required=True)
        cmd.add_argument("--output-dir", required=True)
    for name in ["infer-phylogeny", "run"]:
        cmd = sub.add_parser(name)
        cmd.add_argument("--input-dir")
        cmd.add_argument("--output-dir", required=True)
        cmd.add_argument(
            "--analysis-scope",
            choices=["single-copy", "experimental-multicopy"],
            default="single-copy",
        )
        cmd.add_argument("--model", choices=["exon-locus-ctmc", "exon-parsimony", "exon-ctmc", "parsimony", "er-ard", "foreground"],
                         default="exon-locus-ctmc", help="Explicit genomic DNA-copy CTMC by default; other model choices select their named compatibility path.")
        from .locus import add_locus_options
        add_locus_options(cmd)
        add_configuration_options(cmd)
        cmd.add_argument("--branch-length-mode", choices=["supplied", "unit"], default="supplied")
        cmd.add_argument(
            "--ascertainment",
            choices=["observed-at-least-one", "complete-universe", "variable-only"],
            default="observed-at-least-one",
        )
        cmd.add_argument(
            "--root-frequency",
            choices=["estimated", "stationary", "fixed"],
            default="estimated",
        )
        cmd.add_argument("--root-presence", type=float, default=0.5)
        cmd.add_argument(
            "--annotation-view",
            choices=["repertoire", "canonical"],
            default="repertoire",
            help="Summarize all supplied transcript paths or the explicitly marked canonical path.",
        )
        cmd.add_argument("--threads", type=int, default=1)
        cmd.add_argument("--bootstrap-replicates", type=int, default=0)
        cmd.add_argument("--stochastic-maps", type=int, default=0)
        cmd.add_argument("--seed", type=int, default=7)
        cmd.add_argument("--foreground-branches")
        cmd.add_argument(
            "--structural-site-matrix",
            help="Frozen schema-v3 structural observation matrix shared across model fits.",
        )
        if name == "run":
            cmd.add_argument(
                "--evidence-aligner",
                choices=["internal", "mafft", "minimap2", "lastz", "miniprot"],
                default="minimap2",
                help="Aligner for sequence/protein evidence projection.",
            )
            cmd.add_argument(
                "--short-context-max-length",
                type=int,
                default=300,
                help="Maximum anchor-bounded interval length for the short local DNA route.",
            )

    for command in (sub.choices["run"], sub.choices["infer-phylogeny"]):
        command.add_argument("--analysis-range", choices=["all", "high-coverage"], default="all",
                             help="Fit all callable observations, or an explicit character subset. Never trims DNA.")
        command.add_argument("--min-callable-fraction", type=float, default=0.70,
                             help="Known 0/1 divided by full tree panel; also sets the high-coverage report view.")

    check = sub.add_parser("check", help="Validate genomic FASTA, GFF and tree inputs without running alignments.")
    add_file_inputs(check)
    check.add_argument("--species-tree")
    loci = sub.add_parser("extract-loci", help="Export paired genomic FASTA/GFF loci with automatic flanks.")
    add_file_inputs(loci)
    loci.add_argument("--species-tree", required=True)
    loci.add_argument("--output-dir", required=True)
    loci.add_argument("--flank", type=int, default=1000)
    normalize = sub.add_parser("normalize-annotation", help="Explicitly normalize GFF/GTF with optional AGAT.")
    normalize.add_argument("--gff", nargs="+", required=True)
    normalize.add_argument("--output-dir", required=True)
    normalize.add_argument("--config", help="Optional AGAT YAML configuration.")
    normalize.add_argument("--timeout", type=int, default=600)
    for command in sub.choices.values():
        command.formatter_class = argparse.ArgumentDefaultsHelpFormatter
        command.add_argument("--quiet", action="store_true", help="Suppress console progress; retain logs.")
        command.add_argument("--debug", action="store_true", help="Show a Python traceback on errors.")
    for name in ("build-case", "run", "infer-phylogeny", "visualize", "import-orthofinder", "extract-loci", "normalize-annotation", "explain", "analyze", "prepare-genomic-evidence", "fit-exon-rates"):
        sub.choices[name].add_argument("--force", action="store_true",
            help="Preserve an existing IntraPhy output directory as a backup, then rerun. Never overwrites inputs.")
    return parser
