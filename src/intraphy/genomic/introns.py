"""Prepare physical intron-boundary observations from DNA and annotation."""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from intraphy.coding_correspondence import CodingProjectionIndex
from intraphy.genomic.groups import family_species_roster
from intraphy.genomic.intron_coordinates import physical_key
from intraphy.genomic.intron_sites import ALIGNMENT_FIELDS, OBSERVATION_FIELDS, POSITION_FIELDS
from intraphy.genomic.intron_sites import family_intron_rows, unresolved_intron_rows
from intraphy.preparation.annotation_index import parse_attributes
from intraphy.storage.fasta import parse_fasta
from intraphy.storage.tabular import read_tsv, write_tsv


FAMILY_FIELDS = ("family_id", "species", "n_annotated_introns", "n_projected_sites",
                 "n_eligible_sites", "status", "reasons", "anchor_window",
                 "min_anchor_pairs", "alignment_mode", "alignment_scope", "observation_type")


def prepare_intron_observations(input_dir, output_dir, *, threads=1,
                                anchor_window=15, min_anchor_pairs=8) -> list[dict]:
    """Write intron-position DNA calls and their conditional family-MSA mapping."""
    if type(threads) is not int or threads < 1:
        raise ValueError("threads must be a positive integer")
    if type(anchor_window) is not int or anchor_window < 1:
        raise ValueError("anchor_window must be a positive integer")
    if type(min_anchor_pairs) is not int or min_anchor_pairs < 1:
        raise ValueError("min_anchor_pairs must be a positive integer")
    if min_anchor_pairs > anchor_window:
        raise ValueError("min_anchor_pairs cannot exceed anchor_window")
    input_dir, output_dir = Path(input_dir), Path(output_dir)
    required_inputs = ("transcript_paths.tsv", "intron_sites.tsv", "gene_loci.tsv",
                       "segment_sequences.fasta")
    missing_inputs = [name for name in required_inputs if not (input_dir / name).is_file()]
    if missing_inputs:
        raise FileNotFoundError("missing prepared intron inputs: " + ", ".join(missing_inputs))
    occurrences = read_tsv(input_dir / "segment_occurrences.tsv", required=(
        "occurrence_id", "family_id", "species", "gene_copy_id", "contig", "start", "end", "strand"))
    paths = read_tsv(input_dir / "transcript_paths.tsv", optional=True)
    introns = read_tsv(input_dir / "intron_sites.tsv", optional=True)
    raw = [dict(row, attrs=parse_attributes(row.get("attrs", "")))
           for row in read_tsv(input_dir / "raw_gene_features.tsv", optional=True)]
    loci = defaultdict(list)
    for row in read_tsv(input_dir / "gene_loci.tsv", optional=True):
        loci[(row["species"], row["gene_copy_id"])].append(row)
    sequences = parse_fasta(input_dir / "segment_sequences.fasta")
    family_species, family_copies = family_species_roster(input_dir, {
        str(index): row for index, row in enumerate(occurrences)
    })
    aliases, introns_by_family = defaultdict(list), defaultdict(list)
    for row in introns:
        try:
            physical = physical_key(row)
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("invalid intron_sites.tsv physical coordinates") from exc
        if int(row["start"]) < 1 or int(row["end"]) < int(row["start"]) or row["strand"] not in {"+", "-"}:
            raise ValueError(f"invalid intron interval: {physical}")
        if not aliases[physical]:
            introns_by_family[row["family_id"]].append(row)
        aliases[physical].append(row)

    positions, observations, alignments, family_summaries = [], [], [], []
    for family in sorted(family_species):
        family_occurrences = [row for row in occurrences if row["family_id"] == family]
        family_paths = [row for row in paths if row.get("family_id") == family]
        family_features = [row for row in raw if row.get("family_id") == family]
        index = (CodingProjectionIndex(family_occurrences, sequences, family_paths, family_features,
                 threads=threads, msa_mode="linsi") if family_paths else None)
        in_family = [row for row in introns if row["family_id"] == family]
        aligned, family_positions, family_observations, errors, candidates = family_intron_rows(
            family, index, in_family, aliases, sequences, loci, family_species[family], family_copies,
            anchor_window=anchor_window, min_anchor_pairs=min_anchor_pairs,
        )
        projected = {item["physical"] for items in candidates.values() for item in items
                     if item["kind"] == "intron"}
        unresolved_positions, unresolved_observations = unresolved_intron_rows(
            family, introns_by_family.get(family, ()), projected, aliases, family_species[family], errors,
            anchor_window=anchor_window, min_anchor_pairs=min_anchor_pairs,
            site_ordinal_start=len(family_positions),
        )
        alignments.extend(aligned)
        positions.extend(family_positions)
        positions.extend(unresolved_positions)
        observations.extend(family_observations)
        observations.extend(unresolved_observations)
        annotated_count = len(introns_by_family.get(family, ()))
        projected_count = len(family_positions)
        eligible_count = sum(int(row["eligible"]) for row in family_positions)
        if not annotated_count:
            status, reasons = "no_annotated_introns", "no_annotated_intron_candidates"
        elif not projected_count:
            status = "no_usable_coding_projection"
            reasons = ";".join(sorted({row["reason"] for row in unresolved_positions})) or "no_projected_intron_boundary"
        elif not eligible_count:
            status = "sites_excluded"
            reasons = ";".join(sorted({row["reason"] for row in family_positions})) or "no_eligible_intron_positions"
        else:
            status, reasons = "observations_prepared", "NA"
        family_summaries.append({"family_id": family,
            "species": ";".join(sorted(family_species[family])),
            "n_annotated_introns": annotated_count, "n_projected_sites": projected_count,
            "n_eligible_sites": eligible_count, "status": status, "reasons": reasons,
            "anchor_window": anchor_window, "min_anchor_pairs": min_anchor_pairs,
            "alignment_mode": "linsi",
            "alignment_scope": "family_protein_msa_conditional_position_not_orthology_certification",
            "observation_type": "intron_position_presence"})

    write_tsv(output_dir / "intron_positions.tsv", positions, POSITION_FIELDS)
    write_tsv(output_dir / "intron_observations.tsv", observations, OBSERVATION_FIELDS)
    write_tsv(output_dir / "intron_alignment.tsv", alignments, ALIGNMENT_FIELDS)
    write_tsv(output_dir / "intron_families.tsv", family_summaries, FAMILY_FIELDS)
    return observations
