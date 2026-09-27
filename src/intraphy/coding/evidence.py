"""Pure finalization of typed transcript-pair candidates into public evidence."""
from __future__ import annotations

import json

from intraphy.coding.projection import _candidate_json
from intraphy.coding.types import CodingProjectionCandidateSet


def finalize_evidence(result, candidates, transcripts, msa_mode):
    """Apply the legacy candidate ordering, winner rule, and evidence schema."""
    result["protein_unavailable_reason"] = ";".join(sorted(result.pop("_unavailable", ()))) or (
        "NA" if candidates or result["protein_status"] != "unavailable" else "no_CDS_transcript_path")
    if not candidates:
        return result
    candidates = tuple(sorted(candidates, key=lambda candidate: (
        candidate.query_transcript_key, candidate.target_transcript_key,
        candidate.coordinate_blocks)))
    candidate_set = CodingProjectionCandidateSet(candidates)
    best = max(candidates, key=lambda candidate: (
        candidate.position_eligible, candidate.anchor_resolved,
        candidate.known_aa_pairs, candidate.blosum62_score,
        candidate.aa_identity,
        max(candidate.query_cds_coverage, candidate.target_cds_coverage)))
    ambiguous = (not candidate_set.coordinate_consensus
        or bool(candidate_set.competing_occurrences)
        or any(not candidate.position_monotonic for candidate in candidates))
    mapping_status = ("ambiguous_mapping" if ambiguous else
        "resolved_local" if candidate_set.position_eligible else "supported_unanchored")
    interval = best.msa_column_interval
    result.update(
        protein_status="ambiguous_transcript_projection" if ambiguous else "supported",
        protein_mapping_status=mapping_status,
        protein_membership_eligible=candidate_set.membership_eligible,
        protein_position_eligible=candidate_set.position_eligible,
        protein_hard_observation_eligible=candidate_set.position_eligible,
        protein_candidate_evidence_available=candidate_set.candidate_evidence_available,
        protein_aa_identity=best.aa_identity,
        protein_query_cds_coverage=best.query_cds_coverage,
        protein_target_cds_coverage=best.target_cds_coverage,
        protein_known_aa_pairs=best.known_aa_pairs,
        protein_blosum62_score=best.blosum62_score,
        protein_gap_fraction=best.gap_fraction,
        protein_left_anchor_pairs=best.left_anchor_pairs,
        protein_right_anchor_pairs=best.right_anchor_pairs,
        protein_left_anchor_score=best.left_anchor_score,
        protein_right_anchor_score=best.right_anchor_score,
        protein_left_anchor_supported=best.left_anchor_supported,
        protein_right_anchor_supported=best.right_anchor_supported,
        protein_terminal_side=best.terminal_side,
        protein_msa_column_start0=interval.start0,
        protein_msa_column_end0=interval.end0,
        protein_msa_column_interval=f"{interval.start0}:{interval.end0}",
        protein_msa_mode=msa_mode,
        protein_candidate_mapping_count=len(candidates),
        protein_candidate_coordinate_consensus=candidate_set.coordinate_consensus,
        protein_competing_occurrences=";".join(candidate_set.competing_occurrences) or "NA",
        protein_candidate_details=json.dumps([_candidate_json(c) for c in candidates],
            sort_keys=True, separators=(",", ":")),
        protein_candidate_set=candidate_set,
        protein_query_source_features=";".join(transcripts[best.query_transcript_key].source_ids) or "NA",
        protein_target_source_features=";".join(transcripts[best.target_transcript_key].source_ids) or "NA",
        protein_metrics_scope="best_transcript_pair",
        protein_best_query_transcript=best.query_transcript_id,
        protein_best_target_transcript=best.target_transcript_id,
        protein_supporting_transcripts=";".join(
            f"{candidate.query_transcript_id}>{candidate.target_transcript_id}" for candidate in candidates),
    )
    if not ambiguous:
        result["protein_projected_coordinate_blocks"] = candidates[0].coordinate_blocks
        result["protein_projected_blocks"] = candidates[0].coordinate_blocks
    return result
