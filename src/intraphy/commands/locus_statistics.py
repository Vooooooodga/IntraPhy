"""CLI registration, validation, and dispatch for conditional locus statistics."""
from __future__ import annotations

import math


def add_locus_statistics_commands(sub):
    stats = sub.add_parser("locus-statistics", help="Compare declared locus-rate groups conditionally.")
    stats.add_argument("--locus-model", required=True)
    stats.add_argument("--species-tree", required=True)
    stats.add_argument("--output-dir", required=True)
    stats.add_argument("--null-rate-group", action="append", required=True)
    stats.add_argument("--start-scales", nargs="+", type=float, default=[0.2, 1.0, 5.0])
    stats.add_argument("--profile-rate-group")
    stats.add_argument("--profile-values", nargs="+", type=float)
    stats.add_argument("--bootstrap-replicates", type=int, default=0)
    stats.add_argument("--seed", type=int)
    stats.add_argument("--sampling-design", choices=("fixed_catalogue",))
    stats.add_argument("--threads", type=int, default=1)
    stats.add_argument("--maxiter", type=int, default=1000)
    stats.set_defaults(model="exon-locus-ctmc")

    prepare = sub.add_parser("prepare-locus-evidence", help="Apply a complete tip-by-material DNA evidence table.")
    prepare.add_argument("--locus-model", required=True, help="Template exon-locus JSON model.")
    prepare.add_argument("--species-tree", required=True)
    prepare.add_argument("--evidence-json", required=True, help="JSON array of complete per-tip material evidence rows.")
    prepare.add_argument("--output-dir", required=True)
    prepare.set_defaults(model="exon-locus-ctmc")


def validate_locus_statistics_arguments(args):
    if not args.null_rate_group or len(set(args.null_rate_group)) != len(args.null_rate_group):
        raise ValueError("--null-rate-group must be supplied at least once with unique values")
    if args.threads < 1 or args.maxiter < 1:
        raise ValueError("--threads and --maxiter must be positive integers")
    if args.bootstrap_replicates < 0:
        raise ValueError("--bootstrap-replicates must be nonnegative")
    if any(not math.isfinite(value) or value <= 0 for value in args.start_scales):
        raise ValueError("--start-scales must contain finite positive values")
    if (args.profile_rate_group is None) != (args.profile_values is None):
        raise ValueError("--profile-rate-group and --profile-values must be supplied together")
    if args.profile_values is not None and any(not math.isfinite(value) or value < 0
                                               for value in args.profile_values):
        raise ValueError("--profile-values must contain finite nonnegative values")
    if args.seed is not None and args.seed < 0:
        raise ValueError("--seed must be nonnegative")
    if args.bootstrap_replicates > 0:
        if args.seed is None:
            raise ValueError("--seed is required when --bootstrap-replicates is positive")
        if args.sampling_design != "fixed_catalogue":
            raise ValueError("--sampling-design fixed_catalogue is required for bootstrap replicates")


def validate_loaded_locus_statistics(args, bundle):
    free_groups = {name for name, record in bundle.rates.items() if record["mode"] == "fit"}
    invalid = set(args.null_rate_group) - free_groups
    if invalid:
        raise ValueError(f"Null rate groups must be declared free in the model: {sorted(invalid)}")
    if args.profile_rate_group is not None and args.profile_rate_group not in free_groups:
        raise ValueError("--profile-rate-group must name a free rate group in the model")
    if args.bootstrap_replicates:
        from ..inference.locus_bootstrap import validate_dna_bootstrap_design
        validate_dna_bootstrap_design(bundle)


def dispatch_locus_statistics(args):
    from ..inference.locus_statistics_run import run_locus_statistics

    run_locus_statistics(args._locus_bundle, args.output_dir,
                         null_rate_groups=args.null_rate_group,
                         start_scales=args.start_scales,
                         profile_rate_group=args.profile_rate_group,
                         profile_values=args.profile_values,
                         bootstrap_replicates=args.bootstrap_replicates,
                         seed=args.seed, sampling_design=args.sampling_design,
                         workers=args.threads, maxiter=args.maxiter)


def dispatch_prepare_locus_evidence(args):
    from pathlib import Path

    from ..run_result import RunResult
    from ..storage.tabular import write_tsv
    from ..structure.serialization import write_json

    output = args.output_dir
    targets = [Path(output) / name for name in ("locus_model.json", "species_tree.tsv", "run_result.json")]
    if any(path.exists() or path.is_symlink() for path in targets):
        raise ValueError("Refusing to overwrite existing prepared-locus-evidence output")
    write_json(f"{output}/locus_model.json", args._prepared_model_record)
    write_tsv(f"{output}/species_tree.tsv", args._prepared_locus_bundle.tree_rows,
              ["node_id", "parent_id", "label", "branch_length"])
    RunResult("exon-locus-ctmc", "evidence-preparation", "species_tree.tsv",
              ("locus_model.json", "species_tree.tsv")).write(output)
