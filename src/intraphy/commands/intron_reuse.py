"""Read-only checks that staged intron observations belong to prepared coding inputs."""
import json
from pathlib import Path


def validate_saved_protein_alignments(input_dir, alignment_rows, family_species):
    """Check stored aligned proteins and aliases against current prepared CDS translations."""
    from ..coding.types import FamilyCodingProjection
    from ..coding_correspondence import CodingProjectionIndex
    from ..preparation.annotation_index import parse_attributes
    from ..storage.fasta import parse_fasta
    from ..storage.tabular import read_tsv

    input_dir = Path(input_dir)
    occurrences = read_tsv(input_dir / "segment_occurrences.tsv", [
        "occurrence_id", "family_id", "species", "gene_copy_id", "contig", "start", "end", "strand",
    ])
    paths = read_tsv(input_dir / "transcript_paths.tsv", optional=True)
    raw = [dict(row, attrs=parse_attributes(row.get("attrs", "")))
           for row in read_tsv(input_dir / "raw_gene_features.tsv", optional=True)]
    sequences = parse_fasta(input_dir / "segment_sequences.fasta")
    rows_by_family = {}
    for row in alignment_rows:
        if row["family_id"] not in family_species:
            raise ValueError("Staged intron alignment contains a family absent from prepared targets.")
        rows_by_family.setdefault(row["family_id"], []).append(row)
    indexes = {}
    for family in family_species:
        family_occ = [row for row in occurrences if row["family_id"] == family]
        family_paths = [row for row in paths if row.get("family_id") == family]
        family_raw = [row for row in raw if row.get("family_id") == family]
        index = CodingProjectionIndex(family_occ, sequences, family_paths, family_raw,
                                      threads=1, msa_mode="linsi")
        indexes[family] = index
        usable = {key for key, transcript in index.transcripts.items()
                  if key[0] == family and not transcript.unavailable_reason and transcript.protein}
        saved = rows_by_family.get(family, [])
        if not usable:
            if saved:
                raise ValueError("Staged intron alignment has records but prepared family has no usable CDS translations.")
            continue
        if not saved:
            raise ValueError("Staged intron alignment omits usable prepared CDS translations.")
        aligned_records, aliases_by_record, residue_columns, record_by_transcript = {}, {}, {}, {}
        lengths = set()
        for row in saved:
            record_id = row["record_id"]
            if not record_id or record_id in aligned_records:
                raise ValueError("Staged intron alignment has duplicate or empty record IDs.")
            sequence = row["aligned_protein"]
            try:
                aliases = tuple(tuple(key) for key in json.loads(row["transcript_keys"])
                                if isinstance(key, list) and len(key) == 4)
                columns = tuple(json.loads(row["residue_columns0"]))
            except (TypeError, ValueError, json.JSONDecodeError) as exc:
                raise ValueError("Staged intron alignment alias or residue-column data is invalid.") from exc
            if (not sequence or not aliases or len(aliases) != len(json.loads(row["transcript_keys"]))
                    or any(type(column) is not int for column in columns)
                    or tuple(sorted(set(columns))) != columns):
                raise ValueError("Staged intron alignment has incomplete transcript aliases.")
            if columns != tuple(i for i, aa in enumerate(sequence) if aa != "-"):
                raise ValueError("Staged intron alignment residue columns disagree with its protein sequence.")
            if any(key in record_by_transcript for key in aliases):
                raise ValueError("Staged intron alignment repeats a transcript alias.")
            aligned_records[record_id] = sequence
            aliases_by_record[record_id] = aliases
            residue_columns[record_id] = columns
            lengths.add(len(sequence))
            for key in aliases:
                record_by_transcript[key] = record_id
        if len(lengths) != 1 or set(record_by_transcript) != usable:
            raise ValueError("Staged intron alignment length or transcript alias roster differs from prepared CDS.")
        projection = FamilyCodingProjection(family, "linsi", aligned_records,
                                            aliases_by_record, record_by_transcript,
                                            residue_columns)
        # This checks residue identity and projections without invoking an MSA backend.
        index.use_family_projection(projection)
    return indexes


def validate_observation_details(input_dir, positions, observations, family_species, indexes):
    """Check per-tip catalogue masks, physical membership, and genomic provenance."""
    from ..storage.fasta import parse_fasta
    from ..storage.tabular import read_tsv
    from ..genomic.intron_coordinates import physical_key

    input_dir = Path(input_dir)
    prepared_introns = read_tsv(input_dir / "intron_sites.tsv", optional=True)
    expected_introns = {physical_key(row) for row in prepared_introns}
    occurrence_rows = read_tsv(input_dir / "segment_occurrences.tsv", [
        "occurrence_id", "family_id", "species", "gene_copy_id", "contig", "start", "end", "strand",
    ])
    occurrences = {row["occurrence_id"]: row for row in occurrence_rows}
    sequences = parse_fasta(input_dir / "segment_sequences.fasta")
    loci = read_tsv(input_dir / "gene_loci.tsv", ["species", "gene_copy_id", "genome_fasta"])
    locus_by_copy = {}
    for row in loci:
        locus_by_copy.setdefault((row["species"], row["gene_copy_id"]), []).append(row)
    sites = {(row["family_id"], row["site_id"]): row for row in positions}
    for row in observations:
        site = sites.get((row["family_id"], row["site_id"]))
        if site is None:
            raise ValueError("Staged intron observation has no matching catalogue site.")
        site_eligible = site["eligible"] in {"1", "true", "True"}
        observation_eligible = row["eligible"] in {"1", "true", "True"}
        if site_eligible != observation_eligible:
            raise ValueError("Staged intron observation eligibility disagrees with its site catalogue.")
        try:
            physical_intervals = {tuple(value) for value in json.loads(row["physical_intervals"])}
            evidence = json.loads(row["evidence"])
            site_intervals = {tuple(value) for value in json.loads(site["physical_intervals"])}
        except (TypeError, ValueError, json.JSONDecodeError) as exc:
            raise ValueError("Staged intron physical interval or evidence data is invalid JSON.") from exc
        if any(len(key) != 6 or key[0] != row["family_id"] or key[1] != row["species"]
               for key in physical_intervals):
            raise ValueError("Staged observation physical intervals do not belong to its family/species.")
        if not isinstance(evidence, list) or any(not isinstance(item, dict) for item in evidence):
            raise ValueError("Staged intron endpoint evidence must be a JSON array of records.")
        if not physical_intervals.issubset(site_intervals):
            raise ValueError("Staged observation physical intervals are absent from its catalogue site.")
        if row["state"] != "unknown" and not evidence:
            raise ValueError("Known intron-position state lacks genomic endpoint evidence.")
        for item in evidence:
            physical = tuple(item.get("physical", ()))
            dna = item.get("evidence") or {}
            if (physical not in physical_intervals or item.get("species") != row["species"]
                    or not dna.get("genome_fasta")):
                raise ValueError("Known intron-position evidence lacks matching physical and genome provenance.")
            if row["state"] == "1" and (item.get("kind") != "intron" or physical not in expected_introns):
                raise ValueError("Present intron state does not refer to a prepared annotated physical intron.")
            if row["state"] == "0" and item.get("kind") != "continuous":
                raise ValueError("Absent intron state lacks evidence of genomic coding continuity.")
            if item.get("kind") == "continuous" and physical in expected_introns:
                raise ValueError("Continuous coding evidence conflicts with an annotated intron interval.")
            if item.get("kind") == "intron" and physical not in expected_introns:
                raise ValueError("Intron evidence interval is absent from prepared annotation.")
            left_coordinate = dna.get("left_base_coordinate")
            right_coordinate = dna.get("right_base_coordinate")
            lo, hi = sorted((left_coordinate, right_coordinate)) if type(left_coordinate) is int and type(right_coordinate) is int else (None, None)
            if item.get("kind") == "intron":
                if (lo is None or physical[4:] != (lo + 1, hi - 1)):
                    raise ValueError("Annotated intron interval does not exactly span the mapped coding cut.")
            elif item.get("kind") == "continuous":
                if lo is None or hi != lo + 1 or physical[4:] != (lo, hi):
                    raise ValueError("State-zero evidence does not map adjacent coding bases on one genomic strand.")
            left_id = dna.get("left_occurrence_id")
            right_id = dna.get("right_occurrence_id")
            if any(type(dna.get(name)) is not int for name in (
                    "left_occurrence_offset0", "right_occurrence_offset0",
                    "left_base_coordinate", "right_base_coordinate")):
                raise ValueError("Staged intron endpoint coordinates and offsets must be integers.")
            if any(not isinstance(dna.get(name), str) for name in ("left_base", "right_base")):
                raise ValueError("Staged intron endpoint bases must be strings.")
            transcript_key = tuple(item.get("transcript_key", ()))
            index = indexes.get(row["family_id"])
            if (index is None or transcript_key not in index.transcripts
                    or transcript_key not in index.family_alignments[row["family_id"]].record_by_transcript):
                raise ValueError("Staged intron evidence refers to a transcript absent from its saved MSA.")
            try:
                left_msa_nt = int(site["left_msa_nt0"])
                right_msa_nt = int(site["right_msa_nt0"])
                if (int(row["left_msa_nt0"]) != left_msa_nt
                        or int(row["right_msa_nt0"]) != right_msa_nt):
                    raise ValueError("Observation MSA cut differs from catalogue site.")
            except (KeyError, TypeError, ValueError) as exc:
                raise ValueError("Staged intron evidence has no projected MSA cut.") from exc
            if right_msa_nt != left_msa_nt + 1:
                raise ValueError("Staged intron cut does not join adjacent projected coding nucleotides.")
            transcript = index.transcripts[transcript_key]
            columns = index.family_alignments[row["family_id"]].columns_for(transcript_key)
            bases = []
            for position in (left_msa_nt, right_msa_nt):
                residue_column = position // 3
                residue = columns.index(residue_column) if residue_column in columns else None
                if (residue is None or transcript.codons[residue] is None
                        or transcript.codon_sources[residue] is None):
                    bases.append(None)
                else:
                    bases.append(transcript.codon_sources[residue][position % 3])
            left_base, right_base = bases
            if left_base is None or right_base is None:
                raise ValueError("Staged intron cut is not present in the prepared MSA-to-CDS projection.")
            if left_base.cds_interval.end0 != right_base.cds_interval.start0:
                raise ValueError("Staged intron cut does not join adjacent prepared coding nucleotides.")
            for base, occurrence_id, offset, coordinate_name, base_name in (
                    (left_base, left_id, dna.get("left_occurrence_offset0"), "left_base_coordinate", "left_base"),
                    (right_base, right_id, dna.get("right_occurrence_offset0"), "right_base_coordinate", "right_base")):
                if (base.occurrence_id != occurrence_id or base.occurrence_interval.start0 != offset
                        or base.genome_interval.start0 + 1 != dna.get(coordinate_name)
                        or sequences.get(occurrence_id, "").upper()[offset:offset + 1] != dna.get(base_name)):
                    raise ValueError("Staged intron cut endpoints differ from the current coding projection.")
            for occurrence_id, offset_name in ((left_id, "left_occurrence_offset0"),
                                               (right_id, "right_occurrence_offset0")):
                occurrence = occurrences.get(occurrence_id)
                offset = dna.get(offset_name)
                if (occurrence is None or occurrence["family_id"] != row["family_id"]
                        or occurrence["species"] != row["species"]
                        or occurrence["contig"] != physical[2] or occurrence["strand"] != physical[3]
                        or type(offset) is not int or occurrence_id not in sequences
                        or not 0 <= offset < len(sequences[occurrence_id])):
                    raise ValueError("Staged intron endpoint evidence does not match prepared coding occurrences.")
                expected_coordinate = (int(occurrence["start"]) + offset if occurrence["strand"] == "+"
                                       else int(occurrence["end"]) - offset)
                if (expected_coordinate != dna.get("left_base_coordinate" if offset_name.startswith("left")
                                                   else "right_base_coordinate")
                        or sequences[occurrence_id][offset].upper() != dna.get(
                            "left_base" if offset_name.startswith("left") else "right_base")):
                    raise ValueError("Staged intron endpoint bases or coordinates disagree with prepared occurrences.")
                resources = [entry.get("genome_fasta") for entry in locus_by_copy.get(
                    (occurrence["species"], occurrence["gene_copy_id"]), [])
                    if entry.get("contig") == physical[2] and entry.get("strand") == physical[3]
                    and entry.get("genome_fasta") not in {None, "", "NA"}]
                if not resources or not Path(dna["genome_fasta"]).is_file() or str(Path(dna["genome_fasta"]).resolve()) not in {
                        str(Path(path).resolve()) for path in resources}:
                    raise ValueError("Staged intron evidence genome FASTA differs from prepared loci.")
