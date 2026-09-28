"""CLI dispatch for qualified genomic DNA-presence observations."""
from __future__ import annotations

import json
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
    if getattr(args, "character_type", "dna-presence") == "intron-position":
        from ..genomic.introns import prepare_intron_observations
        return prepare_intron_observations(
            args.input_dir, args.output_dir, threads=args.threads,
            anchor_window=args.intron_anchor_window,
            min_anchor_pairs=args.intron_min_anchor_pairs,
        )
    return prepare_evidence(
        args.input_dir, args.output_dir,
        min_identity=args.survey_min_identity,
        min_coverage=args.survey_min_coverage,
        max_dp_cells=args.survey_max_dp_cells,
        threads=args.threads,
    )


def dispatch_analyze_introns(args):
    """Prepare/reuse intron-position observations, then fit the separate binary model."""
    from ..inference.intron_presence import analyze_intron_positions

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
            allow_unannotated=True, derive_correspondence=False,
        )

    evidence_dir = getattr(args, "genomic_evidence_dir", None)
    if evidence_dir:
        rows = read_tsv(Path(evidence_dir) / "intron_observations.tsv")
    else:
        from ..genomic.introns import prepare_intron_observations
        rows = prepare_intron_observations(
            prepared, output_dir, threads=args.threads,
            anchor_window=args.intron_anchor_window,
            min_anchor_pairs=args.intron_min_anchor_pairs,
        )
    fixed_rates = None
    if args.dna_gain_rate is not None:
        fixed_rates = {"gain": args.dna_gain_rate, "loss": args.dna_loss_rate}
    result = analyze_intron_positions(
        rows, args._species_tree, output_dir,
        root_frequency=args.root_frequency,
        root_presence=args.root_presence,
        fixed_rates=fixed_rates,
        branch_length_unit="supplied_tree_units",
        expected_counts=args.expected_edits,
        workers=args.threads,
    )
    if not evidence_dir:
        _include_intron_preparation_artifacts(output_dir)
    return result


def _include_intron_preparation_artifacts(output_dir):
    """Record preparation files generated in this run, without copying staged inputs."""
    result_path = Path(output_dir) / "run_result.json"
    if not result_path.is_file():
        return
    record = json.loads(result_path.read_text(encoding="utf-8"))
    for name in ("intron_families.tsv", "intron_positions.tsv", "intron_observations.tsv", "intron_alignment.tsv"):
        if (Path(output_dir) / name).is_file() and name not in record.get("artifacts", []):
            record.setdefault("artifacts", []).append(name)
    temporary = result_path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(result_path)


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
