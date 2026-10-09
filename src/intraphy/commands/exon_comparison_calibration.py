"""Argument validation and dispatch for native foreground calibration."""

from ..inference.genomic_exon_comparison_calibration import (
    calibrate_genomic_exon_comparison,
    validate_calibration_parameters,
)


def validate_calibration_arguments(args):
    validate_calibration_parameters(args.scenario, args.foreground_multiplier,
        args.observation_mask, args.gene_rates, args.taxa, args.replicates,
        args.seed, args.threads)


def dispatch_calibration(args):
    return calibrate_genomic_exon_comparison(args.output_dir, args.scenario,
        args.foreground_multiplier, args.observation_mask, args.gene_rates,
        taxa=args.taxa, replicates=args.replicates, seed=args.seed,
        threads=args.threads)
