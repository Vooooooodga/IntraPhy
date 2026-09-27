"""mapping / short alignment: explicit implementation ownership."""
from __future__ import annotations

from intraphy.aligners.short import _validate_max_dp_cells
from intraphy.aligners.short import anchored_short_alignment
from intraphy.aligners.types import AlignmentBackendError
from intraphy.aligners.types import MAX_INTERNAL_DP_CELLS
from intraphy.coordinates import ClosedInterval1
from intraphy.coordinates import CoordinateBlock
from intraphy.coordinates import Interval0
from intraphy.coordinates import genome_interval_to_local
from intraphy.coordinates import local_interval_to_genome
from intraphy.mapping.candidate_coordinates import _coordinate_block0
from intraphy.mapping.candidate_coordinates import _format_alignment_blocks
from intraphy.mapping.candidate_coordinates import _format_genomic_blocks
from intraphy.mapping.candidate_coordinates import _genomic_blocks0
from intraphy.mapping.candidate_records import _candidate_record
from intraphy.mapping.candidate_serialization import _public_candidate_record
from intraphy.mapping.candidate_serialization import _public_gap_blocks
from intraphy.mapping.candidate_serialization import _public_interval
from intraphy.mapping.candidate_serialization import _transpose_candidate_record
from intraphy.mapping.candidate_serialization import _transpose_gap_blocks
from intraphy.mapping.fields import DEFAULT_CORRESPONDENCE_CRITERIA
from intraphy.mapping.flank_context import _anchor_bounded_gap_blocks
from intraphy.mapping.flank_context import _candidate_blocks_for_copy
from intraphy.mapping.flank_context import _context_flank_pairs
from intraphy.mapping.flank_context import _interval_between_transcript_flanks
from intraphy.mapping.flank_context import _owner_occurrence_for_copy
from intraphy.mapping.flank_context import _value_tokens
from intraphy.mapping.match_records import _format_contract_value
from intraphy.mapping.match_records import _format_optional_number
from intraphy.mapping.policies import _candidate_sequence_accepted
from intraphy.mapping.policies import occurrence_copy_key
from intraphy.storage.values import to_float
import json


def _rerun_anchor_bounded_short_candidates(
    rows,
    occurrence_by_id,
    seqs,
    gene_loci,
    transcript_paths=None,
    *,
    short_alignment_max_dp_cells=MAX_INTERNAL_DP_CELLS,
):
    short_alignment_max_dp_cells = _validate_max_dp_cells(short_alignment_max_dp_cells)
    candidate_owner = {}
    candidate_record = {}
    for owner in rows:
        for record in owner.get("_candidate_records", ()):
            candidate_id = record.get("candidate_id")
            if candidate_id not in {None, "", "NA"}:
                candidate_owner[candidate_id] = owner
                candidate_record[candidate_id] = record

    changed = False
    for row in rows:
        if row.get("short_context_route") != "feature_bounded_candidate":
            continue
        row["membership_edge_eligible"] = 0
        row["position_edge_eligible"] = 0
        row["_membership_edge_eligible"] = False
        row["_position_edge_eligible"] = False
        input_transposed = row.get("alignment_input_transposed") in {1, "1", True}
        bounded_side = "query" if input_transposed else "target"
        source_side = "target" if input_transposed else "query"
        bounded_id = (
            row["query_occurrence_id"] if bounded_side == "query"
            else row["subject_occurrence_id"]
        )
        source_id = (
            row["query_occurrence_id"] if source_side == "query"
            else row["subject_occurrence_id"]
        )
        bounded = occurrence_by_id.get(bounded_id, {})
        source = occurrence_by_id.get(source_id, {})
        retained = [
            record for record in row.get("_candidate_records", ())
            if record.get("candidate_id") in _value_tokens(row.get("retained_candidate_ids"))
        ]
        flank_pairs = {
            (
                tuple(sorted(_value_tokens(record.get("left_anchor_id")))),
                tuple(sorted(_value_tokens(record.get("right_anchor_id")))),
            )
            for record in retained
        }
        flank_pairs.discard((tuple(), tuple()))
        inferred_pairs = _context_flank_pairs(
            row,
            source=source,
            bounded=bounded,
            candidate_owner=candidate_owner,
            candidate_record=candidate_record,
            occurrence_by_id=occurrence_by_id,
            transcript_paths=transcript_paths,
        )
        if transcript_paths:
            flank_pairs = inferred_pairs
        if len(flank_pairs) != 1:
            row["match_status"] = "candidate_unanchored"
            row["candidate_resolution"] = "candidate"
            row["incomplete_reason"] = "unique_same_context_flank_pair_unavailable"
            continue
        left_ids, right_ids = next(iter(flank_pairs))
        if len(left_ids) != 1 or len(right_ids) != 1:
            row["match_status"] = "candidate_ambiguous"
            row["candidate_resolution"] = "candidate"
            row["incomplete_reason"] = "competing_same_context_flank_pairs"
            continue

        bounded_copy = occurrence_copy_key(bounded)
        locus = gene_loci.get((bounded.get("species"), bounded.get("gene_copy_id")))
        if (
            locus is None
            or locus.get("contig") != bounded.get("contig")
            or locus.get("strand") != bounded.get("strand")
        ):
            row["match_status"] = "candidate_unanchored"
            row["candidate_resolution"] = "candidate"
            row["incomplete_reason"] = "target_gene_locus_sequence_unavailable"
            continue

        left_id, right_id = left_ids[0], right_ids[0]
        left_record = candidate_record.get(left_id)
        right_record = candidate_record.get(right_id)
        left_owner = candidate_owner.get(left_id)
        right_owner = candidate_owner.get(right_id)
        if not left_record or not right_record or not left_owner or not right_owner:
            row["match_status"] = "candidate_unanchored"
            row["candidate_resolution"] = "candidate"
            row["incomplete_reason"] = "flank_candidate_coordinates_unavailable"
            continue
        flank_occurrences = (
            _owner_occurrence_for_copy(
                left_owner, occurrence_copy_key(source), occurrence_by_id,
            ),
            _owner_occurrence_for_copy(
                right_owner, occurrence_copy_key(source), occurrence_by_id,
            ),
            _owner_occurrence_for_copy(
                left_owner, bounded_copy, occurrence_by_id,
            ),
            _owner_occurrence_for_copy(
                right_owner, bounded_copy, occurrence_by_id,
            ),
        )
        expected_geometry = (
            (source.get("contig"), source.get("strand")),
            (source.get("contig"), source.get("strand")),
            (bounded.get("contig"), bounded.get("strand")),
            (bounded.get("contig"), bounded.get("strand")),
        )
        if any(
            occurrence is None
            or (occurrence.get("contig"), occurrence.get("strand")) != expected
            for occurrence, expected in zip(flank_occurrences, expected_geometry)
        ):
            row["match_status"] = "candidate_ambiguous"
            row["candidate_resolution"] = "candidate"
            row["incomplete_reason"] = "flanks_are_not_on_matching_contigs_and_gene_strands"
            continue
        left_blocks = _candidate_blocks_for_copy(
            left_record, left_owner, bounded_copy, occurrence_by_id,
        )
        right_blocks = _candidate_blocks_for_copy(
            right_record, right_owner, bounded_copy, occurrence_by_id,
        )
        search_interval = _interval_between_transcript_flanks(
            left_blocks, right_blocks, bounded.get("strand"),
        )
        if (
            search_interval is None
            or search_interval.length == 0
            or search_interval.start0 < locus["interval"].start0
            or search_interval.end0 > locus["interval"].end0
        ):
            row["match_status"] = "candidate_ambiguous"
            row["candidate_resolution"] = "candidate"
            row["incomplete_reason"] = "flanks_do_not_define_one_contained_genome_interval"
            continue
        local_search = genome_interval_to_local(
            search_interval, locus["interval"], locus["strand"],
        )
        target_sequence = locus["sequence"][local_search.start0:local_search.end0]
        query_sequence = seqs.get(source_id, "")
        search_metadata = {
            "coordinate_system": "0-based-half-open",
            "contig": locus["contig"],
            "start0": search_interval.start0,
            "end0": search_interval.end0,
            "strand": locus["strand"],
        }
        try:
            candidate_set = anchored_short_alignment(
                query_sequence,
                target_sequence,
                mode="local",
                query_occurrence_id=source_id,
                target_occurrence_id=bounded_id,
                query_transcript_id=source.get("transcript_id"),
                target_transcript_id=bounded.get("transcript_id"),
                left_anchor_id=left_id,
                right_anchor_id=right_id,
                search_interval=search_metadata,
                max_dp_cells=short_alignment_max_dp_cells,
            )
        except AlignmentBackendError as error:
            row["match_status"] = "candidate_unanchored"
            row["candidate_resolution"] = "candidate"
            row["short_context_route"] = "anchor_bounded_unavailable"
            row["incomplete_reason"] = str(error)
            continue

        records = []
        try:
            bounded_locus = ClosedInterval1(
                int(bounded["start"]), int(bounded["end"]),
            ).to_interval0()
        except (KeyError, TypeError, ValueError):
            bounded_locus = None
        for rank, candidate in enumerate(candidate_set.candidates, start=1):
            record = _candidate_record(
                candidate, rank, candidate.backend, candidate.score_scheme,
            )
            source_genomic = _genomic_blocks0(
                source, record["aligned_blocks"], "query",
            )
            bounded_genomic = tuple(
                local_interval_to_genome(
                    block.target, search_interval, locus["strand"],
                )
                for block in record["aligned_blocks"]
            )
            try:
                if bounded_locus is None:
                    raise ValueError("bounded parent feature is unavailable")
                bounded_local = tuple(
                    genome_interval_to_local(block, bounded_locus, bounded.get("strand"))
                    for block in bounded_genomic
                )
                parent_gap_blocks = _anchor_bounded_gap_blocks(
                    record.get("gap_blocks", ()),
                    search_interval,
                    bounded_locus,
                    bounded.get("strand"),
                )
            except (KeyError, TypeError, ValueError):
                bounded_local = tuple()
                parent_gap_blocks = []
            if (
                len(source_genomic) == len(record["aligned_blocks"])
                and len(bounded_local) == len(record["aligned_blocks"])
            ):
                row_blocks = tuple(
                    CoordinateBlock(source_block.query, target_local)
                    for source_block, target_local in zip(
                        record["aligned_blocks"], bounded_local,
                    )
                )
                if input_transposed:
                    row_blocks = tuple(
                        CoordinateBlock(block.target, block.query) for block in row_blocks
                    )
                    parent_gap_blocks = _transpose_gap_blocks(parent_gap_blocks)
            else:
                row_blocks = tuple()
                parent_gap_blocks = []
            if input_transposed:
                record = _transpose_candidate_record(record)
                record["query_genomic_blocks0"] = bounded_genomic
                record["target_genomic_blocks0"] = source_genomic
            else:
                record["query_genomic_blocks0"] = source_genomic
                record["target_genomic_blocks0"] = bounded_genomic
            record["aligned_blocks"] = row_blocks
            record["gap_blocks"] = parent_gap_blocks
            if row_blocks:
                record["query_start0"] = min(block.query.start0 for block in row_blocks)
                record["query_end0"] = max(block.query.end0 for block in row_blocks)
                record["target_start0"] = min(block.target.start0 for block in row_blocks)
                record["target_end0"] = max(block.target.end0 for block in row_blocks)
            record["candidate_id"] = (
                f"{row['match_id']}.anchor_bounded_candidate_{rank:03d}"
            )
            record["left_anchor_id"] = left_id
            record["right_anchor_id"] = right_id
            record["search_interval"] = search_metadata
            record["search_interval_side"] = bounded_side
            record["source"] = "nucleotide_alignment"
            record["short_sequence_coverage"] = candidate.query_coverage
            record["accepted"] = int(bool(
                row_blocks
                and _candidate_sequence_accepted(
                    record, to_float(row.get("threshold"), 0.0), short_context=True,
                    anchored_microexon=(
                        3 <= len(query_sequence) < DEFAULT_CORRESPONDENCE_CRITERIA.short_min_aligned_pairs
                        and candidate_set.enumeration_complete
                        and len(candidate_set.candidates) == 1
                        and bool(left_id) and bool(right_id) and left_id != right_id
                        and all(base in "ACGT" for base in query_sequence.upper())
                    ),
                )
            ))
            record["acceptance_threshold"] = row.get("threshold", "NA")
            records.append(record)

        final_candidate_ids = [record["candidate_id"] for record in records]
        for record in records:
            record["hit_count"] = len(records)
            record["alternative_candidate_ids"] = tuple(
                candidate_id for candidate_id in final_candidate_ids
                if candidate_id != record["candidate_id"]
            )

        row["_candidate_records"] = records
        row["short_context_route"] = "anchor_bounded_local"
        row["flanking_anchor_status"] = "ordered_double_sided_homologous_flanks_same_path"
        row["left_anchor_id"] = left_id
        row["right_anchor_id"] = right_id
        row["search_interval"] = _format_contract_value(_public_interval(search_metadata))
        row["search_interval_side"] = bounded_side
        row["local_boundary_range"] = row["search_interval"]
        row["enumeration_complete"] = int(candidate_set.enumeration_complete)
        row["candidate_enumeration_status"] = (
            "complete" if candidate_set.enumeration_complete else "incomplete"
        )
        row["incomplete_reason"] = candidate_set.incomplete_reason or "NA"
        row["hit_count"] = len(records)
        row["ambiguous_hit_count"] = max(0, len(records) - 1)
        accepted = [record for record in records if record.get("accepted") == 1]
        row["match_status"] = "mapped" if accepted else "candidate_low_similarity"
        row["candidate_resolution"] = "unassessed"
        primary = records[0] if records else None
        if primary is not None:
            row["alignment_score"] = f"{to_float(primary.get('identity'), 0.0):.6g}"
            row["coverage_score"] = f"{to_float(primary.get('coverage'), 0.0):.6g}"
            row["sequence_score"] = f"{(0.70 * to_float(primary.get('identity'), 0.0) + 0.30 * to_float(primary.get('coverage'), 0.0)):.6g}"
            row["correspondence_score"] = row["sequence_score"]
            row["candidate_id"] = primary["candidate_id"]
            row["candidate_ids"] = ";".join(record["candidate_id"] for record in records)
            row["alternative_candidate_ids"] = ";".join(
                record["candidate_id"] for record in records[1:]
            ) or "NA"
            row["matched_blocks"] = _format_alignment_blocks(primary["aligned_blocks"])
            row["aligned_blocks"] = row["matched_blocks"]
            row["projected_reference_blocks"] = row["matched_blocks"]
            query_blocks = tuple(
                _coordinate_block0(block).query for block in primary["aligned_blocks"]
            )
            target_blocks = tuple(
                _coordinate_block0(block).target for block in primary["aligned_blocks"]
            )
            if query_blocks:
                public_query = ClosedInterval1.from_interval0(Interval0(
                    min(block.start0 for block in query_blocks),
                    max(block.end0 for block in query_blocks),
                ))
                public_target = ClosedInterval1.from_interval0(Interval0(
                    min(block.start0 for block in target_blocks),
                    max(block.end0 for block in target_blocks),
                ))
                row["query_alignment_start"] = public_query.start
                row["query_alignment_end"] = public_query.end
                row["target_alignment_start"] = public_target.start
                row["target_alignment_end"] = public_target.end
            row["query_genomic_matched_blocks"] = _format_genomic_blocks(
                occurrence_by_id[row["query_occurrence_id"]].get("contig", "NA"),
                occurrence_by_id[row["query_occurrence_id"]].get("strand", "NA"),
                primary.get("query_genomic_blocks0", ()),
            )
            row["subject_genomic_matched_blocks"] = _format_genomic_blocks(
                occurrence_by_id[row["subject_occurrence_id"]].get("contig", "NA"),
                occurrence_by_id[row["subject_occurrence_id"]].get("strand", "NA"),
                primary.get("target_genomic_blocks0", ()),
            )
            row["gap_blocks"] = json.dumps(
                _public_gap_blocks(primary.get("gap_blocks", ())),
                sort_keys=True, separators=(",", ":"),
            )
            row["alignment_cigar"] = primary.get("cigar", "NA")
            row["raw_alignment_score"] = _format_optional_number(primary.get("score"))
            row["raw_score"] = row["raw_alignment_score"]
            row["short_sequence_coverage"] = primary.get("short_sequence_coverage", "NA")
        public_records = [_public_candidate_record(record) for record in records]
        row["candidate_assessments"] = json.dumps(
            public_records, sort_keys=True, separators=(",", ":"),
        )
        row["dna_candidate_assessments"] = row["candidate_assessments"]
        row["alternative_hits"] = json.dumps(
            public_records[1:], sort_keys=True, separators=(",", ":"),
        )
        changed = True
    return changed
