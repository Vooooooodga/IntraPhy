"""The IntraPhy command-line interface."""

from .commands.parser import build_parser
import json
from pathlib import Path

from . import __version__
from .alignment import available_alignment_backends
from .annotation import complete_annotation, generate_sequence_evidence
from .baseline import evaluate_baselines
from .benchmark import benchmark_events
from .calibration import DEFAULT_SCENARIOS, calibrate_simulations
from .case import build_case, inspect_annotation, scan_hidden_segments
from .correspondence import infer_correspondence
from .orthofinder import import_orthofinder
from .phylogeny import infer_phylogeny
from .preprocess import assess_short_candidate_thresholds, derive_tables, extract_gene
from .simulate import simulate_dataset
from .visualize import visualize_results


from .workflow import run_all, _record_evidence_aligner



def _dispatch(args):
    if args.command == "prepare-genomic-evidence":
        from .commands.genomic import dispatch_prepare_evidence
        dispatch_prepare_evidence(args)
        return
    if args.command in {"locus-statistics", "prepare-locus-evidence"}:
        from .commands.locus_statistics import (
            dispatch_locus_statistics, dispatch_prepare_locus_evidence,
        )
        if args.command == "locus-statistics":
            dispatch_locus_statistics(args)
        else:
            dispatch_prepare_locus_evidence(args)
        return
    if args.command in {"infer-exon-repertoires", "fit-exon-repertoire-rates"}:
        from .inference.repertoire_run import infer_repertoires, fit_repertoire_rates
        if args.command == "infer-exon-repertoires":
            infer_repertoires(args._repertoire_bundle, args.output_dir, expected_edits=args.expected_edits)
        else:
            fit_repertoire_rates(args._repertoire_bundle, args.output_dir, log_scale_bounds=args.log_scale_bounds, collection_provenance=args.collection_provenance)
        return
    if args.command == "realign-exons":
        from .aligners.cesar_adapter import realign
        realign(args)
    elif args.command == "fit-exon-rates":
        from .commands.exon_statistics import fit_command
        fit_command(args)
    elif args.command == "analyze":
        if args.model == "exon-locus-ctmc":
            from .commands.locus import dispatch_locus
            dispatch_locus(args)
        elif args.model == "dna-presence-ctmc":
            from .commands.genomic import dispatch_analyze
            dispatch_analyze(args)
        elif args.model == "intron-position-ctmc":
            from .commands.genomic import dispatch_analyze_introns
            dispatch_analyze_introns(args)
        else:
            from .commands.exons import dispatch_analyze
            dispatch_analyze(args)
    elif args.command == "example-exons":
        from .verification.exon_cases import write_exon_example
        write_exon_example(args.output_dir, args.scenario, args.seed)
    elif args.command == "exon-rate-template":
        from .structure.edits import EDIT_KINDS
        from .structure.serialization import write_json
        import math
        if not math.isfinite(args.rate) or args.rate < 0:
            raise ValueError("Rate must be finite and nonnegative")
        if Path(args.output).exists():
            raise ValueError("Rate output already exists")
        write_json(args.output, {"schema": "intraphy.exon-rates/1", "rates": {k: args.rate for k in EDIT_KINDS},
                               "provenance": "user-requested illustrative fixed rates; not fitted biological estimates"})
    elif args.command in {"run", "infer-phylogeny"} and args.model == "exon-locus-ctmc":
        from .commands.locus import dispatch_locus
        dispatch_locus(args)
    elif args.command in {"run", "infer-phylogeny"} and args.model.startswith("exon-"):
        from .commands.exons import configuration_arguments
        from .inference.configuration_run import infer_configurations
        infer_configurations(args.input_dir, args.output_dir, **configuration_arguments(args))
    elif args.command == "example":
        from .verification.native_cases import build_native_example
        manifest = build_native_example(args.output_dir, args.scenario, args.seed)
        print(f"Synthetic FASTA/GFF/tree inputs written: {Path(manifest).parent}")
    elif args.command == "explain":
        if args.legacy_v18:
            from .reporting.methods import render_guide
        else:
            from .reporting.exon_guide import render_guide
        render_guide(args.output_dir, png=args.png)
    elif args.command == "check":
        from .commands.environment import environment_report
        selection = args._input_selection
        print(json.dumps({"status": "valid", "input_mode": selection.source_mode,
                          "species": list(selection.tree_species),
                          "families": sorted({r["family_id"] for r in selection.rows}),
                          "selected_gene_loci": len(selection.rows),
                          "gene_only_unknown_loci": sum(r.get("annotation_structure_status") == "gene_only_unknown" for r in selection.rows),
                          "environment": environment_report()}, indent=2))
    elif args.command == "extract-loci":
        from .inputs.loci import export_loci
        export_loci(args._input_selection, args.species_tree, args.output_dir, args.flank)
    elif args.command == "normalize-annotation":
        from .inputs.agat import normalize_annotations
        normalize_annotations(args.gff, args.output_dir, args.config, args.timeout)
    elif args.command == "extract-gene":
        extract_gene(args.genome, args.annotation, args.gene_id, args.family_id, args.species, args.gene_copy_id, args.output_dir, args.append, args.transcript_policy, args.canonical_rule, args.source_label, args.copy_role, flank=args.flank, max_extension=args.max_extension)
    elif args.command == "derive-tables":
        derive_tables(
            args.input_dir,
            args.output_dir,
            args.identity_threshold,
            args.distance_table,
            args.aligner,
            args.threads,
            args.min_size_ratio,
            context_aligner=args.context_aligner,
            coding_msa_mode=args.coding_msa_mode,
            short_context_max_length=args.short_context_max_length,
            short_alignment_max_dp_cells=args.short_alignment_max_dp_cells,
        )
    elif args.command == "simulate":
        simulate_dataset(args.output_dir, args.seed, args.scenario)
    elif args.command == "benchmark":
        benchmark_events(args.input_dir, args.output_dir)
    elif args.command == "calibrate":
        calibrate_simulations(args.output_dir, args.scenario, args.replicates, args.bootstrap_replicates, args.stochastic_maps, args.seed)
    elif args.command == "visualize":
        from .run_result import result_model
        model_name = str(result_model(args.result_dir))
        if model_name == "exon-locus-ctmc":
            raise ValueError("The exon-locus-ctmc result currently provides auditable JSON posterior and event-count outputs; no locus-specific graphic renderer is available.")
        if model_name == "dna-presence-ctmc":
            raise ValueError("The dna-presence-ctmc result contains DNA presence histories; no structural graphic renderer is available.")
        if model_name == "intron-position-ctmc":
            raise ValueError("The intron-position-ctmc result contains intron-position histories; no structural graphic renderer is available.")
        if model_name.startswith("exon-"):
            from .reporting.exon_results import render_exon_results
            render_exon_results(args.result_dir, args.output_dir)
            return
        encoding = args.correspondence_encoding or "color"
        visualize_results(args.input_dir, args.result_dir, args.output_dir, encoding,
                          layout=args.layout, targets=args.target,
                          target_manifest=args.target_manifest)
    elif args.command == "inspect-aligners":
        from .commands.environment import inspect_tools
        print("tool\tavailable\tpath\tversion")
        for row in inspect_tools():
            print(f"{row['tool']}\t{row['available']}\t{row['path'] or 'NA'}\t{row['version'] or 'NA'}")
    elif args.command == "inspect-annotation":
        inspect_annotation(args.annotation, args.output_dir, args.query, args.alias_file, args.species, args.case_id)
    elif args.command == "build-case":
        args._input_selection.write(Path(args.output_dir))
        build_case(
            args.manifest,
            args.output_dir,
            args.identity_threshold,
            args.species_tree,
            args.transcript_policy,
            args.canonical_rule,
            args.aligner,
            args.threads,
            args.min_size_ratio,
            args.copy_tree,
            args.gene_tree,
            flank=args.flank,
            max_extension=args.max_extension,
            context_aligner=args.context_aligner,
            coding_msa_mode=args.coding_msa_mode,
            short_context_max_length=args.short_context_max_length,
            target_rows=args._input_selection.rows,
            allow_unannotated=args.allow_unannotated_loci,
        )
    elif args.command == "import-orthofinder":
        import_orthofinder(args.orthofinder_dir, args.orthogroup, args.genome_manifest, args.output_dir, args.species_tree)
    elif args.command == "scan-hidden-segments":
        scan_hidden_segments(args.source_fasta, args.target_fasta, args.output_dir, args.family_id, args.species, args.gene_copy_id, args.min_identity, args.min_coverage, args.aligner, args.threads)
    elif args.command == "candidate-sensitivity":
        assess_short_candidate_thresholds(
            args.segment_matches,
            args.output,
            args.identity,
            args.coverage,
        )
    elif args.command == "complete-annotation":
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
        complete_annotation(args.input_dir, args.output_dir)
    elif args.command == "segment-correspondence":
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
        infer_correspondence(args.input_dir, args.output_dir)
    elif args.command == "infer-phylogeny":
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
        infer_phylogeny(
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            bootstrap_replicates=args.bootstrap_replicates,
            stochastic_maps=args.stochastic_maps,
            seed=args.seed,
            foreground_branches=args.foreground_branches,
            analysis_scope=args.analysis_scope,
            model=args.model,
            branch_length_mode=args.branch_length_mode,
            ascertainment=args.ascertainment,
            threads=args.threads,
            root_frequency=args.root_frequency,
            root_presence=args.root_presence,
            structural_site_matrix_path=args.structural_site_matrix,
            analysis_range=args.analysis_range,
            min_callable_fraction=args.min_callable_fraction,
            annotation_view=args.annotation_view,
        )
    elif args.command == "compare-baselines":
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
        evaluate_baselines(args.input_dir, args.output_dir)
    elif args.command == "run":
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
        run_all(
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            bootstrap_replicates=args.bootstrap_replicates,
            stochastic_maps=args.stochastic_maps,
            seed=args.seed,
            foreground_branches=args.foreground_branches,
            analysis_scope=args.analysis_scope,
            model=args.model,
            branch_length_mode=args.branch_length_mode,
            ascertainment=args.ascertainment,
            threads=args.threads,
            root_frequency=args.root_frequency,
            root_presence=args.root_presence,
            evidence_aligner=args.evidence_aligner,
            short_context_max_length=args.short_context_max_length,
            annotation_view=args.annotation_view,
            structural_site_matrix_path=args.structural_site_matrix,
            analysis_range=args.analysis_range,
            min_callable_fraction=args.min_callable_fraction,
        )


from .commands.session import command_session


def main(argv=None):
    import sys
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        with command_session(args):
            _dispatch(args)
        return 0
    except KeyboardInterrupt:
        print("IntraPhy interrupted; see execution.json for the failed run.", file=sys.stderr)
        return 130
    except (OSError, ValueError, RuntimeError, SystemExit) as exc:
        if getattr(args, "debug", False):
            raise
        if isinstance(exc, SystemExit) and exc.code in {None, 0}:
            return 0
        print(f"IntraPhy error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
