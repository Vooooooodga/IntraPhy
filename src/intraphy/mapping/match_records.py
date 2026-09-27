"""mapping / match_records: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from collections import defaultdict
from intraphy.candidate_chain import DEFAULT_CHAIN_CONFIGURATION
from intraphy.coordinates import ClosedInterval1
from intraphy.coordinates import Interval0
from intraphy.coordinates import parse_legacy_blocks
from intraphy.mapping.candidate_codec import _coordinate_block0
from intraphy.mapping.candidate_codec import _covered_bases
from intraphy.mapping.candidate_codec import _format_alignment_blocks
from intraphy.mapping.candidate_codec import _format_genomic_blocks
from intraphy.mapping.candidate_codec import _genomic_blocks0
from intraphy.mapping.candidate_codec import _public_candidate_record
from intraphy.mapping.candidate_codec import _public_interval
from intraphy.mapping.candidate_records import _scope_candidate_record_ids
from intraphy.mapping.policies import occurrence_copy_key
import json


def _projection_interval(evidence, side):
    if side not in {"query", "target"}:
        raise ValueError(f"invalid projection side: {side}")
    blocks = ()
    if (
        "annotated_CDS_protein" in str(evidence.get("correspondence_basis", ""))
        and evidence.get("protein_hard_observation_eligible") in {1, "1", True}
    ):
        blocks = tuple(
            _coordinate_block0(block)
            for block in evidence.get("protein_projected_blocks", ())
        )
    if not blocks:
        accepted = [
            record
            for record in evidence.get("candidate_records", ())
            if record.get("accepted", 1) in {1, "1", True}
        ]
        if accepted:
            blocks = tuple(
                _coordinate_block0(block)
                for block in accepted[0].get("aligned_blocks", ())
            )
    if blocks:
        intervals = [getattr(block, side) for block in blocks]
        return Interval0(
            min(interval.start0 for interval in intervals),
            max(interval.end0 for interval in intervals),
        )
    if side == "query":
        start = evidence.get("query_alignment_start", "NA")
        end = evidence.get("query_alignment_end", "NA")
    else:
        start = evidence.get("target_alignment_start", "NA")
        end = evidence.get("target_alignment_end", "NA")
    if start in {"NA", None, ""} or end in {"NA", None, ""}:
        return None
    start, end = int(start), int(end)
    if end < start:
        return None
    return ClosedInterval1(start, end).to_interval0()


def _projection_record(evidence, side):
    interval = _projection_interval(evidence, side)
    if interval is None:
        return None
    return {
        "interval": interval,
        "strand": (
            "+"
            if "annotated_CDS_protein" in str(evidence.get("correspondence_basis", ""))
            and evidence.get("protein_position_eligible") in {1, "1", True}
            else evidence.get("alignment_strand", "NA")
        ),
    }


def _ordered_projection_compatible(left, right, left_ref_interval, right_ref_interval):
    if not left_ref_interval or not right_ref_interval:
        return False
    if not _disjoint_reference_intervals(left_ref_interval, right_ref_interval):
        return False
    if left.get("contig") != right.get("contig") or left.get("strand") != right.get("strand"):
        return False
    left_order = int(left["start"])
    right_order = int(right["start"])
    if left.get("strand") == "-":
        left_order, right_order = -left_order, -right_order
    copy_order = -1 if left_order < right_order else 1
    ref_order = -1 if left_ref_interval.start0 < right_ref_interval.start0 else 1
    return copy_order == ref_order


def _disjoint_reference_intervals(left_ref_interval, right_ref_interval):
    if not left_ref_interval or not right_ref_interval:
        return False
    return not left_ref_interval.overlaps(right_ref_interval)


def _projection_compatibility(projection_by_occ_ref, occurrence_by_id):
    same_copy_compatible = set()
    cross_copy_compatible = set()
    by_reference = defaultdict(list)
    for (occ_id, ref_id), projection in projection_by_occ_ref.items():
        if projection:
            by_reference[ref_id].append((occ_id, projection))
    for projected in by_reference.values():
        for left_index, (left_id, left_projection) in enumerate(projected):
            left = occurrence_by_id.get(left_id, {})
            for right_id, right_projection in projected[left_index + 1 :]:
                right = occurrence_by_id.get(right_id, {})
                pair = frozenset((left_id, right_id))
                if occurrence_copy_key(left) == occurrence_copy_key(right):
                    if (
                        left_projection.get("strand") == "+"
                        and right_projection.get("strand") == "+"
                        and _ordered_projection_compatible(
                            left,
                            right,
                            left_projection["interval"],
                            right_projection["interval"],
                        )
                    ):
                        same_copy_compatible.add(pair)
                elif _disjoint_reference_intervals(left_projection["interval"], right_projection["interval"]):
                    cross_copy_compatible.add(pair)
    return same_copy_compatible, cross_copy_compatible


def _candidate_records_for_match(evidence, match_id):
    records = [dict(record) for record in evidence.get("candidate_records", ())]
    if not records and evidence.get("projected_reference_blocks") not in {None, "", "NA"}:
        sequence_kind = (
            "amino_acid"
            if str(evidence.get("score_scheme", "")).startswith("blosum")
            else evidence.get("sequence_kind", "nucleotide")
        )
        blocks = list(parse_legacy_blocks(evidence["projected_reference_blocks"]))
        records.append(
            {
                "rank": 1,
                "identity": evidence.get("alignment_score", "NA"),
                "coverage": evidence.get("coverage_score", "NA"),
                "query_start0": min((block.query.start0 for block in blocks), default="NA"),
                "query_end0": max((block.query.end0 for block in blocks), default="NA"),
                "target_start0": min((block.target.start0 for block in blocks), default="NA"),
                "target_end0": max((block.target.end0 for block in blocks), default="NA"),
                "strand": evidence.get("alignment_strand", "+"),
                "mapping_quality": evidence.get("mapping_quality", "NA"),
                "is_secondary": 0,
                "score": evidence.get("raw_alignment_score", evidence.get("alignment_score", "NA")),
                "cigar": evidence.get("alignment_cigar", "NA"),
                "aligned_blocks": blocks,
                "query_transcript_id": evidence.get(
                    "protein_best_query_transcript",
                    evidence.get("query_transcript_ids", "NA"),
                ),
                "target_transcript_id": evidence.get(
                    "protein_best_target_transcript",
                    evidence.get("subject_transcript_ids", "NA"),
                ),
                "backend": evidence.get("alignment_backend", "NA"),
                "score_scheme": evidence.get("score_scheme", "unspecified"),
                "source": (
                    "protein_msa_projection"
                    if str(evidence.get("score_scheme", "")).startswith("blosum")
                    else "nucleotide_alignment"
                ),
                "query_coverage": evidence.get("query_coverage", "NA"),
                "target_coverage": evidence.get("target_coverage", "NA"),
                "aligned_pairs": evidence.get("aligned_pairs", 0),
                "known_aligned_pairs": evidence.get(
                    "protein_known_aa_pairs" if sequence_kind == "amino_acid" else "known_aligned_pairs",
                    evidence.get("aligned_pairs", 0),
                ),
                "unknown_aligned_pairs": evidence.get("unknown_aligned_pairs", 0),
                "query_covered_bases": _covered_bases(blocks, "query"),
                "target_covered_bases": _covered_bases(blocks, "target"),
                "query_length": evidence.get("query_length", "NA"),
                "target_length": evidence.get("target_length", "NA"),
                "gap_blocks": evidence.get("gap_blocks", []),
                "sequence_kind": sequence_kind,
                "backend_version": evidence.get("backend_version", "NA"),
                "raw_score": evidence.get("raw_score", evidence.get("raw_alignment_score", "NA")),
                "nt_identity": (
                    evidence.get("nt_identity", evidence.get("alignment_score", "NA"))
                    if sequence_kind == "nucleotide" else "NA"
                ),
                "aa_identity": (
                    evidence.get("aa_identity", evidence.get("protein_aa_identity", "NA"))
                    if sequence_kind == "amino_acid" else "NA"
                ),
                "relative_strand": evidence.get("relative_strand", evidence.get("alignment_strand", "+")),
                "mapq": evidence.get("mapq", evidence.get("mapping_quality", "NA")),
                "search_interval": evidence.get("search_interval", evidence.get("local_boundary_range", "NA")),
                "search_interval_side": evidence.get("search_interval_side", "target"),
                "accepted": int(bool(evidence.get("candidate_accepted", True))),
            }
        )
    _scope_candidate_record_ids(records, match_id, "candidate")
    return records


def _format_optional_number(value):
    if value in {None, "", "NA"}:
        return "NA"
    return f"{float(value):.6g}"


def _format_contract_value(value):
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    return value


def _match_row(left, right, evidence, score, threshold, distance_class, status, match_index):
    match_id = f"match_{match_index:05d}"
    candidate_records = _candidate_records_for_match(evidence, match_id)
    dna_candidate_records = [
        dict(record)
        for record in evidence.get("dna_candidate_assessments", evidence.get("candidate_records", ()))
    ]
    _scope_candidate_record_ids(dna_candidate_records, match_id, "dna_candidate")
    for record in candidate_records + dna_candidate_records:
        blocks = record.get("aligned_blocks", ())
        record.setdefault("query_genomic_blocks0", _genomic_blocks0(left, blocks, "query"))
        record.setdefault("target_genomic_blocks0", _genomic_blocks0(right, blocks, "target"))
        record["query_genomic_matched_blocks"] = _format_genomic_blocks(
            left.get("contig", "NA"), left.get("strand", "NA"),
            record["query_genomic_blocks0"],
        )
        record["subject_genomic_matched_blocks"] = _format_genomic_blocks(
            right.get("contig", "NA"), right.get("strand", "NA"),
            record["target_genomic_blocks0"],
        )
        record["query_parent_feature_ids"] = evidence.get(
            "query_parent_feature_ids", left.get("source_feature_id", "NA"),
        )
        record["subject_parent_feature_ids"] = evidence.get(
            "subject_parent_feature_ids", right.get("source_feature_id", "NA"),
        )
        record["query_transcript_ids"] = evidence.get(
            "query_transcript_ids", left.get("transcript_id", "NA"),
        )
        record["subject_transcript_ids"] = evidence.get(
            "subject_transcript_ids", right.get("transcript_id", "NA"),
        )
    primary = candidate_records[0] if candidate_records else {}
    public_candidates = [_public_candidate_record(record) for record in candidate_records]
    public_dna_candidates = [
        _public_candidate_record(record) for record in dna_candidate_records
    ]
    public_primary = public_candidates[0] if public_candidates else {}
    row = {
        "match_id": match_id,
        "query_occurrence_id": left["occurrence_id"],
        "subject_occurrence_id": right["occurrence_id"],
        "alignment_score": f"{evidence['alignment_score']:.6g}",
        "coverage_score": f"{evidence['coverage_score']:.6g}",
        "sequence_score": f"{evidence.get('sequence_score', 0.0):.6g}",
        "structural_context_score": f"{evidence.get('structural_context_score', 0.5):.6g}",
        "query_coverage": f"{evidence.get('query_coverage', 0.0):.6g}",
        "target_coverage": f"{evidence.get('target_coverage', 0.0):.6g}",
        "aligned_pairs": evidence.get("aligned_pairs", 0),
        "query_alignment_start": evidence.get("query_alignment_start", "NA"),
        "query_alignment_end": evidence.get("query_alignment_end", "NA"),
        "target_alignment_start": evidence.get("target_alignment_start", "NA"),
        "target_alignment_end": evidence.get("target_alignment_end", "NA"),
        "query_mapped_contig": evidence.get("query_mapped_contig", "NA"),
        "query_mapped_start": evidence.get("query_mapped_start", "NA"),
        "query_mapped_end": evidence.get("query_mapped_end", "NA"),
        "query_mapped_strand": evidence.get("query_mapped_strand", "NA"),
        "query_mapped_length": evidence.get("query_mapped_length", "NA"),
        "subject_mapped_contig": evidence.get("subject_mapped_contig", "NA"),
        "subject_mapped_start": evidence.get("subject_mapped_start", "NA"),
        "subject_mapped_end": evidence.get("subject_mapped_end", "NA"),
        "subject_mapped_strand": evidence.get("subject_mapped_strand", "NA"),
        "subject_mapped_length": evidence.get("subject_mapped_length", "NA"),
        "alignment_strand": evidence.get("alignment_strand", "NA"),
        "projected_reference_occurrence_id": evidence.get("projected_reference_occurrence_id", "NA"),
        "projected_reference_start": evidence.get("projected_reference_start", "NA"),
        "projected_reference_end": evidence.get("projected_reference_end", "NA"),
        "projected_reference_blocks": evidence.get("projected_reference_blocks", "NA"),
        "matched_blocks": evidence.get("matched_blocks", evidence.get("projected_reference_blocks", "NA")),
        "query_genomic_matched_blocks": evidence.get("query_genomic_matched_blocks", "NA"),
        "subject_genomic_matched_blocks": evidence.get("subject_genomic_matched_blocks", "NA"),
        "query_parent_feature_ids": evidence.get("query_parent_feature_ids", left.get("source_feature_id", "NA")),
        "subject_parent_feature_ids": evidence.get("subject_parent_feature_ids", right.get("source_feature_id", "NA")),
        "query_transcript_ids": evidence.get("query_transcript_ids", left.get("transcript_id", "NA")),
        "subject_transcript_ids": evidence.get("subject_transcript_ids", right.get("transcript_id", "NA")),
        "candidate_id": primary.get("candidate_id", "NA"),
        "alternative_candidate_ids": ";".join(
            record["candidate_id"] for record in public_candidates[1:]
        ) or "NA",
        "aligned_blocks": _format_alignment_blocks(primary.get("aligned_blocks", ())),
        "gap_blocks": json.dumps(public_primary.get("gap_blocks", []), sort_keys=True, separators=(",", ":")),
        "sequence_kind": primary.get("sequence_kind", evidence.get("sequence_kind", "nucleotide")),
        "backend": primary.get("backend", evidence.get("backend", evidence.get("alignment_backend", "NA"))),
        "backend_version": primary.get("backend_version", evidence.get("backend_version", "NA")),
        "raw_score": _format_optional_number(primary.get("raw_score", primary.get("score", evidence.get("raw_score")))),
        "nt_identity": _format_optional_number(primary.get("nt_identity", evidence.get("nt_identity"))),
        "aa_identity": _format_optional_number(primary.get("aa_identity", evidence.get("aa_identity"))),
        "known_aligned_pairs": primary.get("known_aligned_pairs", evidence.get("known_aligned_pairs", "NA")),
        "unknown_aligned_pairs": primary.get("unknown_aligned_pairs", evidence.get("unknown_aligned_pairs", "NA")),
        "query_covered_bases": primary.get("query_covered_bases", evidence.get("query_covered_bases", "NA")),
        "target_covered_bases": primary.get("target_covered_bases", evidence.get("target_covered_bases", "NA")),
        "query_length": primary.get("query_length", evidence.get("query_length", "NA")),
        "target_length": primary.get("target_length", evidence.get("target_length", "NA")),
        "relative_strand": primary.get("relative_strand", evidence.get("relative_strand", evidence.get("alignment_strand", "NA"))),
        "is_secondary": primary.get("is_secondary", 0),
        "left_anchor_id": primary.get("left_anchor_id", "NA"),
        "right_anchor_id": primary.get("right_anchor_id", "NA"),
        "search_interval": _format_contract_value(
            _public_interval(primary.get(
                "search_interval", evidence.get("search_interval", evidence.get("local_boundary_range", "NA")),
            ))
        ),
        "search_interval_side": primary.get(
            "search_interval_side", evidence.get("search_interval_side", "NA"),
        ),
        "short_sequence_coverage": primary.get("short_sequence_coverage", evidence.get("short_sequence_coverage", "NA")),
        "alignment_input_transposed": primary.get("alignment_input_transposed", evidence.get("alignment_input_transposed", 0)),
        "mapq": primary.get("mapq", evidence.get("mapping_quality", "NA")),
        "mapping_quality": evidence.get("mapping_quality", "NA"),
        "hit_count": evidence.get("hit_count", len(candidate_records)),
        "ambiguous_hit_count": evidence.get("ambiguous_hit_count", max(0, len(candidate_records) - 1)),
        "alternative_hits": json.dumps(public_candidates[1:], sort_keys=True, separators=(",", ":")),
        "candidate_assessments": json.dumps(public_candidates, sort_keys=True, separators=(",", ":")),
        "dna_candidate_assessments": json.dumps(
            public_dna_candidates,
            sort_keys=True,
            separators=(",", ":"),
        ),
        "candidate_ids": ";".join(record["candidate_id"] for record in candidate_records) or "NA",
        "score_scheme": evidence.get("score_scheme", "unspecified"),
        "raw_alignment_score": _format_optional_number(evidence.get("raw_alignment_score")),
        "enumeration_complete": int(bool(evidence.get("enumeration_complete", False))),
        "candidate_enumeration_status": evidence.get("candidate_enumeration_status", "unassessed"),
        "incomplete_reason": evidence.get("incomplete_reason", "NA") or "NA",
        "short_context_route": evidence.get("short_context_route", "not_used"),
        "local_boundary_range": _format_contract_value(_public_interval(
            evidence.get("local_boundary_range", "not_evaluated")
        )),
        "candidate_resolution": "unassessed",
        "retained_candidate_ids": "NA",
        "best_path_candidate_ids": "NA",
        "chain_best_score": "NA",
        "chain_score_delta": "NA",
        "chain_configuration": DEFAULT_CHAIN_CONFIGURATION.name,
        "chain_delta_rule": DEFAULT_CHAIN_CONFIGURATION.delta_rule,
        "chain_local_mode": "NA",
        "chain_status": "unassessed",
        "chain_ambiguity": "unassessed",
        "chain_start_anchor_ids": "NA",
        "chain_end_anchor_ids": "NA",
        "chain_retained_edges": "NA",
        "chain_best_path_count_capped": 0,
        "chain_near_optimal_path_count_capped": 0,
        "flanking_anchor_status": evidence.get("flanking_anchor_status", "not_evaluated_pre_chain"),
        "membership_edge_eligible": 0,
        "membership_edge_reason": "candidate_chain_not_evaluated",
        "position_edge_eligible": 0,
        "position_edge_reason": "candidate_chain_not_evaluated",
        "true_absence_eligible": 0,
        "true_absence_evidence_status": evidence.get("true_absence_evidence_status", "insufficient_evidence"),
        "true_absence_reason": evidence.get("true_absence_reason", "ordered_double_flanks_not_evaluated;anchor_interval_sequence_not_extracted;assembly_continuity_unassessed;ambiguous_base_status_unassessed;query_only_deletion_gap_unassessed;alternative_alignment_concordance_unassessed"),
        "left_context_score": f"{evidence['left_context_score']:.6g}",
        "right_context_score": f"{evidence['right_context_score']:.6g}",
        "boundary_score": f"{evidence['boundary_score']:.6g}",
        "phase_score": f"{evidence['phase_score']:.6g}",
        "order_score": f"{evidence['order_score']:.6g}",
        "strand_score": f"{evidence['strand_score']:.6g}",
        "splice_score": f"{evidence['splice_score']:.6g}",
        "size_ratio": f"{evidence['size_ratio']:.6g}",
        "total_score": f"{evidence['total_score']:.6g}",
        "distance_class": distance_class,
        "threshold": f"{threshold:.6g}",
        "alignment_cigar": evidence["alignment_cigar"],
        "alignment_backend": evidence["alignment_backend"],
        "match_status": status,
        "alignment_mode": evidence["alignment_mode"],
        "alignment_meaning": evidence["alignment_meaning"],
        "alignment_requested_backend": evidence["alignment_requested_backend"],
        "dna_match_status": evidence.get("dna_match_status", status),
        "correspondence_basis": evidence.get("correspondence_basis", "DNA"),
        "correspondence_score": f"{score:.6g}",
        "protein_status": evidence.get("protein_status", "not_requested"),
        "protein_unavailable_reason": evidence.get("protein_unavailable_reason", "NA"),
        "protein_metrics_scope": evidence.get("protein_metrics_scope", "NA"),
        "protein_best_query_transcript": evidence.get("protein_best_query_transcript", "NA"),
        "protein_best_target_transcript": evidence.get("protein_best_target_transcript", "NA"),
        "protein_supporting_transcripts": evidence.get("protein_supporting_transcripts", "NA"),
        "protein_projected_blocks": _format_alignment_blocks(evidence.get("protein_projected_blocks", ())),
        "protein_mapping_status": evidence.get("protein_mapping_status", "uncovered"),
        "protein_candidate_evidence_available": int(bool(
            evidence.get("protein_candidate_evidence_available", False)
        )),
        "protein_membership_eligible": int(bool(evidence.get("protein_membership_eligible", False))),
        "protein_position_eligible": int(bool(evidence.get("protein_position_eligible", False))),
        "protein_hard_observation_eligible": int(bool(evidence.get("protein_hard_observation_eligible", False))),
        "protein_terminal_side": evidence.get("protein_terminal_side", "NA"),
        "protein_msa_column_start0": evidence.get("protein_msa_column_start0", "NA"),
        "protein_msa_column_end0": evidence.get("protein_msa_column_end0", "NA"),
        "protein_msa_column_interval": evidence.get("protein_msa_column_interval", "NA"),
        "protein_msa_mode": evidence.get("protein_msa_mode", "NA"),
        "protein_competing_occurrences": evidence.get("protein_competing_occurrences", "NA"),
        "protein_candidate_coordinate_consensus": int(bool(evidence.get("protein_candidate_coordinate_consensus", False))),
        "protein_candidate_details": evidence.get("protein_candidate_details", "NA"),
        "protein_query_source_features": evidence.get("protein_query_source_features", "NA"),
        "protein_target_source_features": evidence.get("protein_target_source_features", "NA"),
        **{
            field: f"{evidence[field]:.6g}" if field in evidence else "NA"
            for field in (
                "protein_aa_identity",
                "protein_query_cds_coverage",
                "protein_target_cds_coverage",
                "protein_blosum62_score",
                "protein_gap_fraction",
                "protein_left_anchor_score",
                "protein_right_anchor_score",
            )
        },
        **{
            field: evidence.get(field, "NA")
            for field in (
                "protein_known_aa_pairs",
                "protein_left_anchor_pairs",
                "protein_right_anchor_pairs",
                "protein_left_anchor_supported",
                "protein_right_anchor_supported",
                "protein_candidate_mapping_count",
            )
        },
    }
    row["_candidate_records"] = candidate_records
    return row
