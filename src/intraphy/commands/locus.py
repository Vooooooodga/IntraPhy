"""Shared command-line options and dispatch for the exon-locus CTMC."""


def add_locus_options(command):
    """Register the explicit model/tree interface on supported commands."""
    command.add_argument(
        "--locus-model",
        help="Complete exon-locus JSON model (intraphy.exon-locus-model/1).",
    )
    command.add_argument(
        "--species-tree",
        help="Fixed rooted species tree for --locus-model; Newick or TSV.",
    )
    command.add_argument(
        "--parameter-mode",
        choices=("fit", "fixed"),
        default="fit",
        help="Fit declared global rate groups or evaluate their supplied values.",
    )


def dispatch_locus(args):
    """Run the preflight-loaded locus bundle through the common result writer."""
    from ..inference.locus_run import analyze_locus

    analyze_locus(
        args._locus_bundle,
        args.output_dir,
        parameter_mode=args.parameter_mode,
        expected_counts=args.expected_edits,
        workers=args.threads,
    )
