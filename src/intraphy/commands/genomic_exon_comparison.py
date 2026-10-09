"""Preflight and input contracts for the native exon foreground comparison."""
from __future__ import annotations

import math
from pathlib import Path

from ..inputs.foreground import read_canonical_foreground
from ..inputs.species_tree import read_species_tree_rows
from ..structure.serialization import read_catalogues
from ..structure.prepare import read_preparation_summary, validate_preparation_catalogues
from ..structure.tree_context import canonical_tree
from ..structure.validation import validate_collection
from ..topology import SpeciesTree


def validate_comparison_arguments(args):
    if args.threads < 1:
        raise ValueError("--threads must be positive")
    for name in ("max_states", "max_origin_scenarios"):
        value = getattr(args, name)
        if value is not None and value < 1:
            raise ValueError(f"--{name.replace('_', '-')} must be positive when supplied")
    multipliers = args.profile_multipliers or []
    if any(not math.isfinite(value) or value < 0 for value in multipliers):
        raise ValueError("--profile-multipliers values must be finite and nonnegative")
    return multipliers


def validate_comparison_inputs(args):
    input_dir = Path(args.input_dir)
    if not input_dir.is_dir():
        raise FileNotFoundError(f"--input-dir does not exist: {input_dir}")
    tree_path = input_dir / "species_tree.tsv"
    if not tree_path.is_file():
        raise FileNotFoundError(f"Prepared input lacks {tree_path.name}")
    catalog_paths = [Path(value) for value in args.exon_configurations]
    for path in catalog_paths:
        if not path.is_file():
            raise FileNotFoundError(f"--exon-configurations file does not exist: {path}")
    foreground_path = Path(args.foreground_branches)
    if not foreground_path.is_file():
        raise FileNotFoundError(f"--foreground-branches file does not exist: {foreground_path}")
    rows = read_species_tree_rows(tree_path)
    input_tree = SpeciesTree(rows)
    context = canonical_tree(input_tree)
    read_canonical_foreground(foreground_path, context, args.branch_length_mode)
    summary_path = getattr(args, "preparation_summary", None)
    preparation_roster = None
    if summary_path is not None:
        summary_path = Path(summary_path)
        if not summary_path.is_file():
            raise FileNotFoundError(f"--preparation-summary file does not exist: {summary_path}")
        preparation_roster = read_preparation_summary(summary_path)
    catalogues = tuple(c for path in catalog_paths for c in read_catalogues(
        path, allow_empty=preparation_roster is not None))
    if not catalogues and preparation_roster is None:
        raise ValueError("At least one exon catalogue record is required")
    catalogues = validate_collection(catalogues, allow_empty=preparation_roster is not None)
    if any(c.observation_unit != "genomic_exon_spans" for c in catalogues):
        raise ValueError("Catalogues must declare observation_unit='genomic_exon_spans'")
    validate_preparation_catalogues(preparation_roster, catalogues)
    unexpected = {o.species for c in catalogues for o in c.observations} - set(context.tree.leaf_by_label)
    if unexpected:
        raise ValueError("Catalogue contains species absent from the prepared tree: "
                         + ", ".join(sorted(unexpected)))
