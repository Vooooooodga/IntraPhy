"""CLI path for the exon-centered genomic structure model."""
from pathlib import Path


def dispatch_analyze(args):
    from ..case import build_case
    from ..inference.genomic_exon_run import infer_genomic_exons

    output = Path(args.output_dir)
    if getattr(args, "input_dir", None):
        prepared = Path(args.input_dir)
    else:
        prepared = output / "prepared_inputs"
        args._input_selection.write(prepared)
        build_case(
            None, prepared, species_tree=args.species_tree,
            target_rows=args._input_selection.rows, flank=args.flank,
            max_extension=args.max_extension, threads=args.threads,
            allow_unannotated=True, derive_correspondence=False,
        )

    return infer_genomic_exons(
        prepared, output, configurations=args.exon_configurations,
        rates=args.exon_rates, parameter_mode=args.parameter_mode,
        max_states=args.max_states, max_origins=args.max_origin_scenarios,
        branch_length_mode=args.branch_length_mode,
        expected_edits=args.expected_edits,
        alignment_timeout=args.alignment_timeout,
        max_locus_bases=args.max_locus_bases,
        exon_identity=args.exon_identity, anchor_bases=args.anchor_bases,
        anchor_identity=args.anchor_identity, threads=args.threads,
        species_tree=getattr(args, "species_tree", None),
        alignment_evidence_dir=getattr(args, "_alignment_evidence_dir",
                                       getattr(args, "alignment_evidence_dir", None)),
    )
