"""Explicit pooled statistical commands, separate from single-family histories."""
from __future__ import annotations
from pathlib import Path
from ..storage.tabular import read_tsv, write_tsv
from ..topology import SpeciesTree
from ..structure.serialization import write_json
from ..inference.configuration_run import load_rates
from ..inference.exon_rates import load_units, fit_scale, gene_bootstrap
from ..inference.exon_resampling import foreground_bootstrap


def add_statistics_commands(sub):
    cmd = sub.add_parser("fit-exon-rates", help="Fit a shared rate scale across declared independent genes.")
    cmd.add_argument("--exon-configurations", nargs="+", required=True)
    cmd.add_argument("--species-tree", required=True, help="Prepared species_tree.tsv")
    cmd.add_argument("--exon-rates", required=True, help="Fixed relative operation rates and their provenance.")
    cmd.add_argument("--output-dir", required=True)
    cmd.add_argument("--max-states", type=int, default=None,
                     help="Optional explicit state-space limit; default enumerates the full declared closure.")
    cmd.add_argument("--observation-view", choices=("annotation", "evidence"), default="evidence")
    cmd.add_argument("--gene-bootstrap", type=int, default=0)
    cmd.add_argument("--foreground-child", action="append", default=[])
    cmd.add_argument("--parametric-bootstrap", type=int, default=0)
    cmd.add_argument("--seed", type=int, default=19)


def fit_command(args):
    tree = SpeciesTree(read_tsv(args.species_tree))
    units, excluded = load_units(args.exon_configurations, tree, max_states=args.max_states,
                                 observation_view=args.observation_view)
    base = load_rates(args.exon_rates)
    if base.foreground or base.foreground_multiplier != 1:
        raise ValueError("A pooled fitting template must not preapply a foreground multiplier")
    fit = fit_scale(units, base)
    out = Path(args.output_dir)
    write_json(out/"rate_fit.json", {**fit, "excluded_units": excluded,
        "discovery_scopes": sorted({u.space.catalogue.discovery for u in units}),
        "scope": "conditional_scale_estimation_not_genomic_error_calibration"})
    if fit.get("valid_for_resampling"):
        write_json(out/"fitted_exon_rates.json", {"schema": "intraphy.exon-rates/1", "rates": base.rates,
            "scale": fit["scale"], "origin_root_weight": base.origin_root_weight, "provenance": "conditional pooled scale fit; see rate_fit.json",
            "gene_collection": sorted({u.family for u in units})})
    if args.gene_bootstrap:
        write_json(out/"gene_bootstrap.json", gene_bootstrap(units, base, args.gene_bootstrap, args.seed))
    if args.foreground_child:
        test = foreground_bootstrap(units, base, args.foreground_child, args.parametric_bootstrap, args.seed)
        write_json(out/"foreground_comparison.json", test)
    elif args.parametric_bootstrap:
        raise ValueError("--parametric-bootstrap requires --foreground-child")
    return fit
