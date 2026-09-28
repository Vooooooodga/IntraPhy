"""CLI dispatch for qualified genomic DNA-presence observations."""
from __future__ import annotations

from pathlib import Path

from ..storage.tabular import read_tsv


def prepare_evidence(input_dir, output_dir, *, min_identity=.7, min_coverage=.8,
                     max_dp_cells=250000, threads=1):
    from ..genomic.observations import prepare_dna_observations

    return prepare_dna_observations(
        input_dir, output_dir, min_identity=min_identity, min_coverage=min_coverage,
        max_dp_cells=max_dp_cells, threads=threads,
    )


def dispatch_prepare_evidence(args):
    return prepare_evidence(
        args.input_dir, args.output_dir,
        min_identity=args.survey_min_identity,
        min_coverage=args.survey_min_coverage,
        max_dp_cells=args.survey_max_dp_cells,
        threads=args.threads,
    )


def dispatch_analyze(args):
    """Prepare or reuse genomic observations, then fit their binary histories."""
    from ..inference.dna_presence import analyze_dna_presence

    output_dir = Path(args.output_dir)
    if getattr(args, "input_dir", None):
        prepared = Path(args.input_dir)
    else:
        from ..case import build_case

        prepared = output_dir / "prepared_inputs"
        args._input_selection.write(prepared)
        build_case(
            None, prepared, species_tree=args.species_tree,
            target_rows=args._input_selection.rows, flank=args.flank,
            max_extension=args.max_extension, threads=args.threads,
            allow_unannotated=True,
        )

    evidence_dir = getattr(args, "genomic_evidence_dir", None)
    if evidence_dir:
        rows = read_tsv(Path(evidence_dir) / "dna_observations.tsv")
    else:
        rows = prepare_evidence(
            prepared, output_dir,
            min_identity=args.survey_min_identity,
            min_coverage=args.survey_min_coverage,
            max_dp_cells=args.survey_max_dp_cells,
            threads=args.threads,
        )
    fixed_rates = None
    if args.dna_gain_rate is not None:
        fixed_rates = {"gain": args.dna_gain_rate, "loss": args.dna_loss_rate}
    return analyze_dna_presence(
        rows, args._species_tree, output_dir,
        root_frequency=args.root_frequency,
        root_presence=args.root_presence,
        fixed_rates=fixed_rates,
        branch_length_unit="supplied_tree_units",
        expected_counts=args.expected_edits,
        workers=args.threads,
    )
