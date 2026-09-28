"""Coordinate and local-anchor helpers for genomic intron observations."""
from __future__ import annotations

from pathlib import Path

from intraphy.aligners.columns import revcomp
from intraphy.coding.projection import _blosum62_score
from intraphy.storage.fasta import read_fasta_interval


def physical_key(row):
    return (row["family_id"], row["species"], row["contig"], row["strand"],
            int(row["start"]), int(row["end"]))


def fasta_locus(row, loci):
    resources = {}
    for locus in loci.get((row["species"], row["gene_copy_id"]), ()):
        if locus.get("contig") != row["contig"] or locus.get("strand") != row["strand"]:
            continue
        try:
            lo = int(locus.get("search_start", locus.get("linked_start")))
            hi = int(locus.get("search_end", locus.get("linked_end")))
        except (TypeError, ValueError):
            continue
        fasta = locus.get("genome_fasta")
        if not fasta or fasta == "NA":
            continue
        key = (str(Path(fasta).resolve()), row["contig"], row["strand"], lo, hi)
        resources[key] = locus
    return next(iter(resources.values())) if len(resources) == 1 else None


def multiple_physical_copy_species(family_copies, family, loci):
    physical = {}
    for species, copy_id in family_copies.get(family, set()):
        rows = loci.get((species, copy_id), ())
        signatures = set()
        for row in rows:
            path = row.get("genome_fasta")
            try:
                start = int(row.get("annotation_start", row.get("linked_start")))
                end = int(row.get("annotation_end", row.get("linked_end")))
            except (TypeError, ValueError):
                start, end = -1, -1
            signatures.add((str(Path(path).resolve()) if path else "NA", row.get("contig"),
                            row.get("strand"), start, end))
        if not signatures:
            signatures.add(("unresolved_copy", copy_id))
        physical.setdefault(species, set()).update(signatures)
    return {species for species, signatures in physical.items() if len(signatures) > 1}


def genomic_evidence(row, left_base, right_base, locus, left_nt, right_nt):
    if locus is None:
        return None, "genome_resource_ambiguous_or_unavailable"
    left_pos = left_base.genome_interval.start0 + 1
    right_pos = right_base.genome_interval.start0 + 1
    lo, hi = min(left_pos, right_pos), max(left_pos, right_pos)
    try:
        locus_start = int(locus.get("search_start", locus.get("linked_start")))
        locus_end = int(locus.get("search_end", locus.get("linked_end")))
    except (TypeError, ValueError):
        return None, "genomic_locus_bounds_unavailable"
    if lo < locus_start or hi > locus_end:
        return None, "genomic_boundary_outside_locus_bounds"
    intron = None
    if hi - lo > 1:
        intron = (lo + 1, hi - 1)
        if intron != (int(row["start"]), int(row["end"])):
            return None, "coding_gap_does_not_match_annotated_intron_interval"
    elif hi - lo != 1:
        return None, "coding_bases_not_consecutive_or_intronic"
    try:
        left_genome = read_fasta_interval(locus["genome_fasta"], row["contig"], left_pos, left_pos)
        right_genome = read_fasta_interval(locus["genome_fasta"], row["contig"], right_pos, right_pos)
    except SystemExit as exc:
        # An absent declared resource is operational failure. A missing contig or
        # out-of-range endpoint is unresolved biological coordinate evidence.
        if not Path(locus["genome_fasta"]).is_file():
            raise FileNotFoundError(f"declared genome FASTA is unavailable: {locus['genome_fasta']}") from exc
        message = str(exc)
        if (message.startswith("FASTA record not found:")
                or message.startswith("FASTA interval exceeds record length")
                or message.startswith("invalid FASTA interval for")):
            return None, "genomic_boundary_sequence_unavailable"
        raise RuntimeError(f"failed reading declared genome FASTA {locus['genome_fasta']}: {message}") from exc
    except OSError as exc:
        raise OSError(f"failed reading declared genome FASTA {locus['genome_fasta']}") from exc
    if left_genome not in {"A", "C", "G", "T"} or right_genome not in {"A", "C", "G", "T"}:
        return None, "genomic_boundary_endpoint_non_acgt"
    if left_base.strand == "-":
        left_genome, right_genome = revcomp(left_genome), revcomp(right_genome)
    if left_genome != left_nt or right_genome != right_nt:
        return None, "genomic_boundary_bases_disagree_with_coding_projection"
    interval = (int(row["start"]), int(row["end"])) if hi - lo > 1 else None
    return {
        "left_base_coordinate": left_pos, "right_base_coordinate": right_pos,
        "left_base": left_genome, "right_base": right_genome,
        "intron_length": interval[1] - interval[0] + 1 if interval else 0,
        "genome_fasta": str(Path(locus["genome_fasta"]).resolve()),
        "evidence_scope": "coding_endpoints_and_annotated_interval;intron_interior_unassessed",
    }, "annotated_intron_interval_and_coding_endpoints_verified" if intron else "consecutive_genomic_coding_bases_verified"


def coding_boundaries(transcript, columns, introns, sequences):
    bases = []
    for residue0, codon in enumerate(transcript.codon_sources):
        if codon is None or transcript.codons[residue0] is None:
            continue
        for offset, base in enumerate(codon):
            source = sequences.get(base.occurrence_id, "").upper()
            occurrence_offset = base.occurrence_interval.start0
            nucleotide = source[occurrence_offset] if occurrence_offset < len(source) else ""
            bases.append((base, 3 * columns[residue0] + offset, offset,
                          nucleotide))
    found = []
    for (left, left_msa, _left_offset, left_nt), (right, right_msa, right_offset, right_nt) in zip(bases, bases[1:]):
        if right.cds_interval.start0 != left.cds_interval.end0 or right_msa != left_msa + 1:
            continue
        if left.contig != right.contig or left.strand != right.strand:
            continue
        left_pos, right_pos = left.genome_interval.start0 + 1, right.genome_interval.start0 + 1
        lo, hi = min(left_pos, right_pos), max(left_pos, right_pos)
        matching = [row for row in introns
                    if row["contig"] == left.contig and row["strand"] == left.strand
                    and int(row["start"]) == lo + 1 and int(row["end"]) == hi - 1
                    and hi - lo > 1]
        if hi - lo == 1:
            found.append((left_msa, right_msa, "continuous", None, left, right,
                          right_offset, left_nt, right_nt))
        else:
            found.extend((left_msa, right_msa, "intron", intron, left, right,
                          right_offset, left_nt, right_nt)
                         for intron in matching)
    return found


def prepared_base_matches(base, nucleotide, sequences):
    sequence = sequences.get(base.occurrence_id, "").upper()
    position = base.occurrence_interval.start0
    return 0 <= position < len(sequence) and sequence[position] == nucleotide


def anchor_metrics(projection, transcript_key, left_nt0, right_nt0, window, minimum):
    aligned = projection.aligned_for(transcript_key)
    focal_left, focal_right = left_nt0 // 3, right_nt0 // 3
    supported = []
    for key in sorted(projection.record_by_transcript):
        if key == transcript_key or key[1] == transcript_key[1]:
            continue
        other = projection.aligned_for(key)
        sides = []
        for columns in (range(max(0, focal_left - window), focal_left),
                        range(focal_right + 1, min(len(aligned), focal_right + window + 1))):
            scores = [_blosum62_score(aligned[col], other[col])
                      for col in columns if aligned[col] not in "-X" and other[col] not in "-X"]
            sides.append((len(scores), sum(scores)))
        if all(count >= minimum and score > 0 for count, score in sides):
            supported.append({"transcript_key": list(key), "left_pairs": sides[0][0],
                              "right_pairs": sides[1][0], "left_score": sides[0][1],
                              "right_score": sides[1][1]})
    return supported
