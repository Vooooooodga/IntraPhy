"""Occurrence-pair evidence construction over a coding projection index."""
from __future__ import annotations

from intraphy.coding.evidence import finalize_evidence
from intraphy.coding.projection import _anchor_metrics, _coordinate_blocks0
from intraphy.coding.types import CodingProjectionCandidate


def evidence(index, query_occurrence, target_occurrence):
    query_keys = index.by_occurrence.get(query_occurrence, [])
    target_keys = index.by_occurrence.get(target_occurrence, [])
    result = {
        "protein_status": "unavailable",
        "protein_mapping_status": "uncovered",
        "protein_unavailable_reason": "no_CDS_transcript_path",
        "protein_membership_eligible": False,
        "protein_position_eligible": False,
        "protein_hard_observation_eligible": False,
        "protein_candidate_evidence_available": False,
    }
    unavailable, candidates = set(), []
    for query_key in query_keys:
        query = index.transcripts[query_key]
        for target_key in target_keys:
            target = index.transcripts[target_key]
            for side, transcript in (("query", query), ("target", target)):
                if transcript.unavailable_reason:
                    unavailable.add(f"{side}:{transcript.key[-1]}:{transcript.unavailable_reason}")
            if query.unavailable_reason or target.unavailable_reason:
                continue
            if not query.coding_lengths.get(query_occurrence) or not target.coding_lengths.get(target_occurrence):
                unavailable.add("no_CDS_in_requested_occurrence")
                continue
            result["protein_status"] = "no_aligned_CDS"
            pairs, known_columns, aligned_query, aligned_target, _inverted = index._pair(
                query_key, target_key,
            )
            record = pairs.get((query_occurrence, target_occurrence))
            if not record:
                continue
            positions = set(record["positions0"])
            anchors = _anchor_metrics(
                record, known_columns, aligned_query, aligned_target,
                query_occurrence, target_occurrence,
            )
            ordered = sorted(positions)
            query_positions0 = {query0 for query0, _target0 in positions}
            target_positions0 = {target0 for _query0, target0 in positions}
            candidates.append(CodingProjectionCandidate(
                query_transcript_key=query_key,
                target_transcript_key=target_key,
                coordinate_blocks=_coordinate_blocks0(positions),
                aa_identity=record["aa_matches"] / record["aa_pairs"],
                query_cds_coverage=len(positions) / query.coding_lengths[query_occurrence],
                target_cds_coverage=len(positions) / target.coding_lengths[target_occurrence],
                known_aa_pairs=record["aa_pairs"],
                blosum62_score=record["blosum62_score"],
                gap_fraction=anchors["gap_fraction"],
                left_anchor_pairs=anchors["left_pairs"],
                right_anchor_pairs=anchors["right_pairs"],
                left_anchor_score=anchors["left_score"],
                right_anchor_score=anchors["right_score"],
                left_anchor_supported=anchors["left_supported"],
                right_anchor_supported=anchors["right_supported"],
                terminal_side=anchors["terminal_side"],
                msa_column_interval=anchors["column_interval"],
                anchor_resolved=anchors["resolved"],
                position_monotonic=all(
                    right[0] > left[0] and right[1] > left[1]
                    for left, right in zip(ordered, ordered[1:])
                ),
                competing_occurrences=index._competing_occurrences(
                    query_key, target_key, query_occurrence, target_occurrence,
                    query_positions0, target_positions0,
                ),
            ))
    result["_unavailable"] = unavailable
    return finalize_evidence(result, candidates, index.transcripts, index.msa_mode)
