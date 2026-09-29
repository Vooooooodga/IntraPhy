"""Preflight contracts for exon-centered genomic analyses."""
import json
import math
from pathlib import Path
from ..inputs.species_tree import read_species_tree_rows
from ..topology import SpeciesTree


def validate_analysis_options(args):
    if getattr(args, "command", None) != "analyze" or getattr(args, "model", None) != "exon-structure-ctmc":
        return
    if getattr(args, "locus_model", None):
        raise ValueError("--locus-model requires --model exon-locus-ctmc.")
    if getattr(args, "genomic_evidence_dir", None):
        raise ValueError("--genomic-evidence-dir contains DNA-presence or intron-position evidence; exon-structure-ctmc derives exon-span observations.")
    if getattr(args, "input_dir", None) and any(getattr(args, name, None) for name in
                                                  ("fasta", "manifest", "gff", "orthologs")):
        raise ValueError("Use either --input-dir prepared inputs or raw FASTA/GFF inputs.")
    if getattr(args, "exon_configurations", None) and (
            not getattr(args, "input_dir", None)
            or any(getattr(args, name, None) for name in ("fasta", "manifest", "gff", "orthologs"))):
        raise ValueError("--exon-configurations reuse requires prepared --input-dir inputs and their species tree.")
    if (getattr(args, "_root_frequency_explicit", False)
            or getattr(args, "_root_presence_explicit", False)):
        raise ValueError("Binary root-frequency controls do not apply to exon-structure-ctmc.")
    if getattr(args, "dna_gain_rate", None) is not None or getattr(args, "dna_loss_rate", None) is not None:
        raise ValueError("Binary gain/loss rates do not apply to exon-structure-ctmc.")
    if getattr(args, "exon_rates", None) and getattr(args, "parameter_mode", "fit") != "fixed":
        raise ValueError("--exon-rates requires --parameter-mode fixed.")
    if not getattr(args, "exon_rates", None) and getattr(args, "parameter_mode", "fit") == "fixed":
        raise ValueError("--parameter-mode fixed requires --exon-rates.")
    for name, default in (("max_observation_scenarios", None),
                          ("root_presence", .5),
                          ("survey_min_identity", .7), ("survey_min_coverage", .8),
                          ("survey_max_dp_cells", 250000), ("intron_anchor_window", 15),
                          ("intron_min_anchor_pairs", 8)):
        if getattr(args, name, default) != default:
            raise ValueError(f"--{name.replace('_', '-')} is unsupported for exon-structure-ctmc.")
    if getattr(args, "_observation_view_explicit", False):
        raise ValueError("--observation-view is unsupported for exon-structure-ctmc.")
    for name in ("survey_min_identity", "survey_min_coverage", "survey_max_dp_cells",
                 "intron_anchor_window", "intron_min_anchor_pairs"):
        if getattr(args, f"_{name}_explicit", False):
            raise ValueError(f"--{name.replace('_', '-')} is unsupported for exon-structure-ctmc.")
    if getattr(args, "origin_root_sensitivity", []):
        raise ValueError("--origin-root-sensitivity is unsupported for exon-structure-ctmc.")
    if getattr(args, "branch_length_mode", "supplied") not in {"supplied", "unit"}:
        raise ValueError("Unsupported branch-length mode for exon-structure-ctmc.")
    for name in ("alignment_timeout", "max_locus_bases", "anchor_bases"):
        if getattr(args, name, 1) < 1:
            raise ValueError(f"--{name.replace('_', '-')} must be positive.")
    for name in ("max_states", "max_origin_scenarios"):
        value = getattr(args, name, None)
        if value is not None and value < 1:
            raise ValueError(f"--{name.replace('_', '-')} must be positive when supplied.")
    for name in ("exon_identity", "anchor_identity"):
        value = getattr(args, name, 0.0)
        if not math.isfinite(value) or value < 0 or value > 1:
            raise ValueError(f"--{name.replace('_', '-')} must lie in [0, 1].")
    if getattr(args, "species_tree", None) is None and not getattr(args, "input_dir", None):
        raise ValueError("Raw exon-structure analysis requires --species-tree.")
    if getattr(args, "exon_configurations", None):
        markers = []
        with Path(args.exon_configurations).open(encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, 1):
                if line.strip():
                    try:
                        markers.append(json.loads(line).get("observation_unit"))
                    except (ValueError, AttributeError) as exc:
                        raise ValueError(f"Invalid exon catalogue record at line {line_number}.") from exc
        if not markers or any(marker != "genomic_exon_spans" for marker in markers):
            raise ValueError("--exon-configurations must contain only catalogues marked observation_unit='genomic_exon_spans'.")
    if getattr(args, "exon_rates", None):
        from ..inference.configuration_run import load_rates
        load_rates(args.exon_rates)


def validate_input_paths(args):
    from ..inputs.selection import resolve_inputs

    if getattr(args, "input_dir", None):
        directory = Path(args.input_dir)
        if not directory.is_dir():
            raise FileNotFoundError(f"--input-dir: directory does not exist: {args.input_dir}")
        required = ["species_tree.tsv"]
        if not getattr(args, "exon_configurations", None):
            required.extend(("gene_loci.tsv", "gene_loci.fasta", "transcript_paths.tsv"))
        for name in required:
            if not (directory / name).is_file():
                raise FileNotFoundError(f"Prepared exon-structure input lacks {name}")
        prepared_rows = read_species_tree_rows(directory / "species_tree.tsv")
        prepared_tree = SpeciesTree(prepared_rows)
        selected_rows = prepared_rows
        selected_tree = prepared_tree
        if getattr(args, "species_tree", None):
            selected_rows = read_species_tree_rows(args.species_tree)
            selected_tree = SpeciesTree(selected_rows)
            if set(selected_tree.leaf_by_label) != set(prepared_tree.leaf_by_label):
                raise ValueError("--species-tree override must have exactly the prepared tree tip panel.")
        if getattr(args, "branch_length_mode", "supplied") == "supplied":
            for _, child in selected_tree.edges():
                selected_tree.branch_length(child)
        args._prepared_input_dir = directory
        args._species_tree = selected_tree
        args._species_tree_rows = selected_rows
        if getattr(args, "exon_configurations", None) and any(
                getattr(args, name, None) for name in ("fasta", "manifest", "gff", "orthologs")):
            raise ValueError("--exon-configurations reuse requires --input-dir prepared inputs.")
        if getattr(args, "_flank_explicit", False) or getattr(args, "_max_extension_explicit", False):
            raise ValueError("--flank and --max-extension apply only while preparing raw genomic inputs.")
        return

    if not getattr(args, "species_tree", None):
        raise ValueError("Raw exon-structure analysis requires --species-tree.")
    tree_rows = read_species_tree_rows(args.species_tree)
    args._species_tree = SpeciesTree(tree_rows)
    if getattr(args, "branch_length_mode", "supplied") == "supplied":
        for _, child in args._species_tree.edges():
            args._species_tree.branch_length(child)
    args._species_tree_rows = tree_rows
    if not (getattr(args, "fasta", None) or getattr(args, "manifest", None)):
        raise ValueError("Raw exon-structure analysis requires --fasta or --manifest.")
    args._input_selection = resolve_inputs(args)
    if set(args._input_selection.tree_species) != set(args._species_tree.leaf_by_label):
        raise ValueError("Resolved raw input species do not match the supplied tree tips.")
