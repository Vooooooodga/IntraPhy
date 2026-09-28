"""Validation specific to prepared genomic intron-position observations."""
from pathlib import Path
import math
import json

from ..storage.tabular import read_tsv


REQUIRED_PREPARED = (
    "gene_loci.tsv", "segment_occurrences.tsv", "segment_sequences.fasta",
    "transcript_paths.tsv", "intron_sites.tsv",
)


def validate_analysis_options(args):
    """Reject controls from other domains and validate the binary intron model."""
    if getattr(args, "model", None) != "intron-position-ctmc":
        if (getattr(args, "intron_anchor_window", 15) != 15
                or getattr(args, "intron_min_anchor_pairs", 8) != 8):
            raise ValueError("Intron anchor controls apply only to --model intron-position-ctmc.")
        return
    if getattr(args, "root_frequency", "stationary") not in {"stationary", "fixed"}:
        raise ValueError("intron-position-ctmc supports --root-frequency stationary or fixed.")
    if getattr(args, "locus_model", None) or getattr(args, "exon_configurations", None) or getattr(args, "exon_rates", None):
        raise ValueError("The intron-position model uses genomic observations, not exon configurations or a locus model.")
    if getattr(args, "input_dir", None) and any(getattr(args, name, None) for name in
                                                  ("fasta", "manifest", "gff", "orthologs")):
        raise ValueError("Use either --input-dir prepared inputs or raw FASTA/GFF inputs.")
    if getattr(args, "genomic_evidence_dir", None) and not getattr(args, "input_dir", None):
        raise ValueError("--genomic-evidence-dir reuse requires --input-dir prepared inputs.")
    if getattr(args, "branch_length_mode", "supplied") != "supplied":
        raise ValueError("intron-position-ctmc uses supplied tree branch lengths.")
    if getattr(args, "parameter_mode", "fit") != "fit":
        raise ValueError("--parameter-mode applies only to exon-locus-ctmc.")
    if (getattr(args, "dna_gain_rate", None) is None) != (getattr(args, "dna_loss_rate", None) is None):
        raise ValueError("--gain-rate and --loss-rate must be supplied together.")
    for field in ("dna_gain_rate", "dna_loss_rate"):
        value = getattr(args, field, None)
        if value is not None and (not math.isfinite(value) or value < 0):
            raise ValueError("Fixed gain and loss rates must be finite and nonnegative.")
    if getattr(args, "root_frequency", "stationary") == "stationary" and getattr(args, "root_presence", .5) != .5:
        raise ValueError("--root-presence is used only with --root-frequency fixed.")
    if getattr(args, "root_frequency", "stationary") == "fixed" and not 0 <= float(getattr(args, "root_presence", .5)) <= 1:
        raise ValueError("--root-presence must lie in [0, 1] for a fixed root distribution.")
    if getattr(args, "intron_anchor_window", 15) < 1 or getattr(args, "intron_min_anchor_pairs", 8) < 1:
        raise ValueError("Intron anchor window and minimum paired anchors must be positive.")
    if getattr(args, "intron_min_anchor_pairs", 8) > getattr(args, "intron_anchor_window", 15):
        raise ValueError("--intron-min-anchor-pairs cannot exceed --intron-anchor-window.")
    legacy_controls = (
        ("observation_view", "evidence"), ("max_locus_bases", 100000),
        ("alignment_timeout", 600), ("exon_identity", .7),
        ("anchor_bases", 12), ("anchor_identity", .8),
        ("survey_min_identity", .7), ("survey_min_coverage", .8),
        ("survey_max_dp_cells", 250000),
    )
    if any(getattr(args, name, default) != default for name, default in legacy_controls) or any(
            getattr(args, name, None) not in (None, [], ()) for name in
            ("max_states", "max_origin_scenarios", "max_observation_scenarios", "origin_root_sensitivity")):
        raise ValueError("DNA-survey and exon-configuration controls do not apply to intron-position-ctmc.")


def validate_prepare_options(args):
    if getattr(args, "threads", 1) < 1:
        raise ValueError("--threads must be at least 1")
    if getattr(args, "character_type", "dna-presence") == "intron-position":
        if getattr(args, "intron_anchor_window", 15) < 1 or getattr(args, "intron_min_anchor_pairs", 8) < 1:
            raise ValueError("Intron anchor window and minimum paired anchors must be positive.")
        if getattr(args, "intron_min_anchor_pairs", 8) > getattr(args, "intron_anchor_window", 15):
            raise ValueError("--intron-min-anchor-pairs cannot exceed --intron-anchor-window.")
        if (getattr(args, "survey_min_identity", .7) != .7
                or getattr(args, "survey_min_coverage", .8) != .8
                or getattr(args, "survey_max_dp_cells", 250000) != 250000):
            raise ValueError("DNA survey thresholds do not apply to --character-type intron-position.")
        return
    if (getattr(args, "intron_anchor_window", 15) != 15
            or getattr(args, "intron_min_anchor_pairs", 8) != 8):
        raise ValueError("Intron anchor controls apply only to --character-type intron-position.")
    if getattr(args, "survey_max_dp_cells", 250000) < 1:
        raise ValueError("--survey-max-dp-cells must be positive.")
    for field in ("survey_min_identity", "survey_min_coverage"):
        value = getattr(args, field, None)
        if value is not None and (not math.isfinite(value) or not 0 <= value <= 1):
            raise ValueError(f"--{field.replace('_', '-')} must be finite and in [0, 1]")


def validate_prepared_input(directory, tree):
    """Require extracted genomic structure and an exact family/tree roster."""
    directory = Path(directory)
    for name in REQUIRED_PREPARED:
        if not (directory / name).is_file():
            raise FileNotFoundError(f"Prepared intron-position input lacks {name}")
    occurrences = read_tsv(directory / "segment_occurrences.tsv", ["family_id", "species", "gene_copy_id"])
    from ..genomic import family_species_roster
    family_species, _ = family_species_roster(directory, {str(i): row for i, row in enumerate(occurrences)})
    tips = set(tree.leaf_by_label)
    for family, species in family_species.items():
        if species != tips:
            raise ValueError(f"Prepared family {family!r} does not match the supplied tree tips.")


def validate_reusable_evidence(evidence_dir, input_dir, tree, *, anchor_window, min_anchor_pairs):
    """Match staged intron sites to the prepared native intervals and tree roster."""
    from ..genomic.intron_coordinates import physical_key

    evidence_dir, input_dir = Path(evidence_dir), Path(input_dir)
    occurrences = read_tsv(input_dir / "segment_occurrences.tsv", ["family_id", "species", "gene_copy_id"])
    from ..genomic import family_species_roster
    family_species, _ = family_species_roster(input_dir, {str(i): row for i, row in enumerate(occurrences)})
    tips = set(tree.leaf_by_label)
    if any(species != tips for species in family_species.values()):
        raise ValueError("Prepared families do not match the supplied tree tips.")
    prepared_introns = read_tsv(input_dir / "intron_sites.tsv", optional=True)
    expected_physical = {physical_key(row) for row in prepared_introns}
    positions = read_tsv(evidence_dir / "intron_positions.tsv", [
        "family_id", "site_id", "site_kind", "eligible", "reason", "physical_intervals",
        "left_msa_nt0", "right_msa_nt0",
        "anchor_window", "min_anchor_pairs", "alignment_mode", "alignment_scope", "observation_type",
    ])
    observations = read_tsv(evidence_dir / "intron_observations.tsv", [
        "family_id", "site_id", "species", "state", "eligible", "physical_intervals",
        "left_msa_nt0", "right_msa_nt0", "evidence",
        "anchor_window", "min_anchor_pairs", "alignment_mode", "alignment_scope", "observation_type",
    ])
    alignments = read_tsv(evidence_dir / "intron_alignment.tsv", [
        "family_id", "record_id", "aligned_protein", "transcript_keys", "residue_columns0",
        "msa_mode", "alignment_scope",
    ])
    scope = "family_protein_msa_conditional_position_not_orthology_certification"
    summaries = read_tsv(evidence_dir / "intron_families.tsv", [
        "family_id", "species", "n_annotated_introns", "n_projected_sites",
        "n_eligible_sites", "status", "reasons", "anchor_window", "min_anchor_pairs",
        "alignment_mode", "alignment_scope", "observation_type",
    ])
    if (len(summaries) != len(family_species)
            or len({row["family_id"] for row in summaries}) != len(summaries)
            or {row["family_id"] for row in summaries} != set(family_species)):
        raise ValueError("Staged intron family summary does not match the prepared target roster.")
    for row in summaries:
        if (set(filter(None, row["species"].split(";"))) != family_species[row["family_id"]]
                or int(row["anchor_window"]) != anchor_window
                or int(row["min_anchor_pairs"]) != min_anchor_pairs
                or row["alignment_mode"] != "linsi"
                or row["alignment_scope"] != scope
                or row["observation_type"] != "intron_position_presence"):
            raise ValueError("Staged intron family summary has a different species roster.")
    expected_ids = set()
    observed_physical = set()
    for row in positions:
        family = row["family_id"]
        if family not in family_species:
            raise ValueError("Staged intron positions contain a family absent from prepared targets.")
        site_id = row["site_id"]
        if (row["observation_type"] != "intron_position_presence"
                or row["alignment_scope"] != scope or row["alignment_mode"] != "linsi"
                or int(row["anchor_window"]) != anchor_window
                or int(row["min_anchor_pairs"]) != min_anchor_pairs):
            raise ValueError("Staged intron positions have incompatible character type, alignment scope, or thresholds.")
        key = (family, site_id)
        if key in expected_ids:
            raise ValueError("Staged intron position catalogue has duplicate family/site IDs.")
        expected_ids.add(key)
        try:
            for physical in json.loads(row["physical_intervals"]):
                observed_physical.add(tuple(physical))
        except (TypeError, ValueError, json.JSONDecodeError) as exc:
            raise ValueError("Staged intron physical intervals are invalid JSON.") from exc
    if not expected_physical.issubset(observed_physical):
        raise ValueError("Staged intron evidence omits one or more prepared annotated intron intervals.")
    summary_by_family = {row["family_id"]: row for row in summaries}
    for family, summary in summary_by_family.items():
        family_physical = {physical for physical in expected_physical if physical[0] == family}
        family_positions = [row for row in positions if row["family_id"] == family]
        projected = sum(row["site_kind"] == "annotated_intron_boundary" for row in family_positions)
        eligible = sum(row["site_kind"] == "annotated_intron_boundary"
                       and row["eligible"] in {"1", "true", "True"} for row in family_positions)
        if (int(summary["n_annotated_introns"]) != len(family_physical)
                or int(summary["n_projected_sites"]) != projected
                or int(summary["n_eligible_sites"]) != eligible):
            raise ValueError("Staged intron summary counts disagree with the prepared physical sites.")
    expected_rows = {(family, site, species) for family, site in expected_ids
                     for species in family_species[family]}
    actual_rows = set()
    for row in observations:
        key = (row["family_id"], row["site_id"], row["species"])
        if key in actual_rows:
            raise ValueError("Staged intron observations contain duplicate family/site/species rows.")
        actual_rows.add(key)
        if (key not in expected_rows or row["observation_type"] != "intron_position_presence"
                or row["alignment_scope"] != scope or row["alignment_mode"] != "linsi"
                or int(row["anchor_window"]) != anchor_window
                or int(row["min_anchor_pairs"]) != min_anchor_pairs
                or row["state"] not in {"0", "1", "unknown"}
                or row["eligible"] not in {"0", "1", "true", "false", "True", "False"}):
            raise ValueError("Staged intron observations have incompatible type, roster, state, or thresholds.")
    if actual_rows != expected_rows:
        raise ValueError("Staged intron observations do not match the prepared site and tree-tip roster.")
    for row in alignments:
        if row["alignment_scope"] != scope or row["msa_mode"] != "linsi":
            raise ValueError("Staged intron alignment has incompatible MSA provenance.")
    from .intron_reuse import validate_observation_details, validate_saved_protein_alignments
    indexes = validate_saved_protein_alignments(input_dir, alignments, family_species)
    validate_observation_details(input_dir, positions, observations, family_species, indexes)
