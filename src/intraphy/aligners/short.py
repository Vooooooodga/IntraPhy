"""aligners / short: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from collections.abc import Mapping
from intraphy.aligners.columns import _stats_from_alignment_columns
from intraphy.aligners.legacy import _candidate_from_legacy_stats
from intraphy.aligners.runner import _backend_version
from intraphy.aligners.types import AlignmentBackendError
from intraphy.aligners.types import AlignmentCandidateSet
from intraphy.aligners.types import AlignmentGap
from intraphy.aligners.types import KNOWN_NT
from intraphy.aligners.types import MAX_INTERNAL_DP_CELLS
from intraphy.aligners.types import MAX_OPTIMAL_ALIGNMENTS
from intraphy.aligners.types import NT_BLASTN_V1
from intraphy.aligners.types import NT_BLASTN_V1_GAP_EXTEND
from intraphy.aligners.types import NT_BLASTN_V1_GAP_OPEN
from intraphy.aligners.types import NT_BLASTN_V1_MATCH
from intraphy.aligners.types import NT_BLASTN_V1_MISMATCH
from intraphy.coordinates import ClosedInterval1
from intraphy.coordinates import Interval0
from typing import Optional


def _validate_max_dp_cells(max_dp_cells):
    if max_dp_cells is not None and (type(max_dp_cells) is not int or max_dp_cells < 1):
        raise ValueError("max_dp_cells must be None or a positive integer")
    return max_dp_cells


def _normalized_nucleotide_sequence(sequence: str) -> str:
    normalized = []
    for base in (sequence or "").upper():
        base = "T" if base == "U" else base
        normalized.append(base if base in KNOWN_NT else "N")
    return "".join(normalized)


def _pairwise_alignment_columns(seq_a: str, seq_b: str, coordinates) -> tuple[str, str]:
    query_cols = []
    target_cols = []
    for index in range(coordinates.shape[1] - 1):
        a_start, a_end = int(coordinates[0, index]), int(coordinates[0, index + 1])
        b_start, b_end = int(coordinates[1, index]), int(coordinates[1, index + 1])
        a_span = a_end - a_start
        b_span = b_end - b_start
        if a_span and b_span:
            if a_span != b_span:
                raise AlignmentBackendError("PairwiseAligner returned an unequal paired coordinate step")
            query_cols.append(seq_a[a_start:a_end])
            target_cols.append(seq_b[b_start:b_end])
        elif a_span:
            query_cols.append(seq_a[a_start:a_end])
            target_cols.append("-" * a_span)
        elif b_span:
            query_cols.append("-" * b_span)
            target_cols.append(seq_b[b_start:b_end])
    return "".join(query_cols), "".join(target_cols)


def _gap_blocks_from_coordinates(coordinates) -> tuple[AlignmentGap, ...]:
    gaps = []
    for index in range(coordinates.shape[1] - 1):
        query = Interval0(int(coordinates[0, index]), int(coordinates[0, index + 1]))
        target = Interval0(int(coordinates[1, index]), int(coordinates[1, index + 1]))
        if query.length == 0 or target.length == 0:
            gaps.append(AlignmentGap(query, target))
    return tuple(gaps)


def _normalized_search_interval(search_interval, target_length: int) -> dict:
    """Adapt a bounded target interval to one internal coordinate convention."""

    if search_interval is None:
        return {
            "coordinate_system": "0-based-half-open",
            "start0": 0,
            "end0": target_length,
        }
    if isinstance(search_interval, Interval0):
        interval = search_interval
        metadata = {}
    elif isinstance(search_interval, Mapping):
        metadata = {
            key: value
            for key, value in search_interval.items()
            if key not in {"coordinate_system", "start", "end", "start0", "end0"}
        }
        if search_interval.get("start0") not in {None, "", "NA"}:
            interval = Interval0(
                int(search_interval["start0"]), int(search_interval["end0"]),
            )
        elif search_interval.get("start") not in {None, "", "NA"}:
            interval = ClosedInterval1(
                int(search_interval["start"]), int(search_interval["end"]),
            ).to_interval0()
        else:
            raise ValueError("search interval requires start0/end0 or start/end")
    else:
        raise TypeError("search interval must be Interval0 or a coordinate mapping")
    if interval.length != target_length:
        raise ValueError(
            "search interval length does not match the bounded target sequence: "
            f"{interval.length} != {target_length}"
        )
    return {
        "coordinate_system": "0-based-half-open",
        **metadata,
        "start0": interval.start0,
        "end0": interval.end0,
    }


def anchored_short_alignment(
    query: str,
    target: str,
    mode: str = "global",
    max_alignments: int = MAX_OPTIMAL_ALIGNMENTS,
    *,
    query_occurrence_id: Optional[str] = None,
    target_occurrence_id: Optional[str] = None,
    query_transcript_id: Optional[str] = None,
    target_transcript_id: Optional[str] = None,
    left_anchor_id: Optional[str] = None,
    right_anchor_id: Optional[str] = None,
    search_interval: Optional[dict] = None,
    max_dp_cells: Optional[int] = MAX_INTERNAL_DP_CELLS,
) -> AlignmentCandidateSet:
    """Return distinct optimal alignments for one anchor-bounded DNA interval.

    Candidate coordinates are 0-based and half-open. Only optimal alignments
    reported by Biopython are enumerated; near-optimal candidates are outside
    this adapter's search contract.
    """
    max_dp_cells = _validate_max_dp_cells(max_dp_cells)
    mode = (mode or "global").lower()
    if mode not in {"global", "local"}:
        raise AlignmentBackendError(f"unsupported bounded short alignment mode: {mode}")
    if max_alignments < 1:
        raise ValueError("max_alignments must be at least 1")

    query = _normalized_nucleotide_sequence(query)
    target = _normalized_nucleotide_sequence(target)
    backend_version = _backend_version("internal")
    bounded_interval = _normalized_search_interval(search_interval, len(target))
    if not query or not target:
        return AlignmentCandidateSet(
            alignment_mode=mode,
            backend_version=backend_version,
            query_length=len(query),
            target_length=len(target),
            query_occurrence_id=query_occurrence_id,
            target_occurrence_id=target_occurrence_id,
            query_transcript_id=query_transcript_id,
            target_transcript_id=target_transcript_id,
            left_anchor_id=left_anchor_id,
            right_anchor_id=right_anchor_id,
            search_interval=bounded_interval,
        )
    cell_count = len(query) * len(target)
    if max_dp_cells is not None and cell_count > max_dp_cells:
        raise AlignmentBackendError(
            f"bounded internal {mode} alignment rejected by the DP resource guard: "
            f"query_occurrence_id={query_occurrence_id!r} "
            f"(length={len(query)}), target_occurrence_id={target_occurrence_id!r} "
            f"(length={len(target)}), cells={cell_count}, budget={max_dp_cells}; "
            "the complete alignment was rejected without truncating DP states"
        )

    from Bio.Align import PairwiseAligner

    aligner = PairwiseAligner()
    aligner.mode = mode
    aligner.match_score = NT_BLASTN_V1_MATCH
    aligner.mismatch_score = NT_BLASTN_V1_MISMATCH
    aligner.open_gap_score = NT_BLASTN_V1_GAP_OPEN
    aligner.extend_gap_score = NT_BLASTN_V1_GAP_EXTEND
    aligner.wildcard = "N"

    candidates: list[AlignmentCandidate] = []
    signatures = set()
    enumeration_complete = True
    incomplete_reason = ""
    for alignment in aligner.align(query, target):
        coordinates = alignment.coordinates
        signature = tuple(tuple(int(value) for value in row) for row in coordinates)
        if signature in signatures:
            continue
        if len(candidates) >= max_alignments:
            enumeration_complete = False
            incomplete_reason = "optimal_alignment_limit_reached"
            break
        signatures.add(signature)
        query_cols, target_cols = _pairwise_alignment_columns(query, target, coordinates)
        stats = _stats_from_alignment_columns(
            query_cols,
            target_cols,
            len(query),
            len(target),
            query_start0=int(coordinates[0, 0]),
            target_start0=int(coordinates[1, 0]),
            score=float(alignment.score),
            backend="internal",
            mode=mode,
            alignment_mode=mode,
            alignment_meaning=f"bounded short nucleotide {mode} alignment",
        )
        stats.score_scheme = NT_BLASTN_V1
        candidates.append(_candidate_from_legacy_stats(stats, _gap_blocks_from_coordinates(coordinates)))
    candidates.sort(
        key=lambda candidate: (
            candidate.target_interval.start0,
            candidate.target_interval.end0,
            candidate.query_interval.start0,
            candidate.query_interval.end0,
            tuple(
                (block.query.start0, block.query.end0, block.target.start0, block.target.end0)
                for block in candidate.aligned_blocks
            ),
            tuple(
                (gap.query.start0, gap.query.end0, gap.target.start0, gap.target.end0)
                for gap in candidate.gap_blocks
            ),
            candidate.cigar,
        )
    )
    candidate_ids = [
        f"{query_occurrence_id or 'query'}->{target_occurrence_id or 'target'}.candidate_{index:03d}"
        for index in range(1, len(candidates) + 1)
    ]
    for index, candidate in enumerate(candidates):
        candidate.candidate_id = candidate_ids[index]
        candidate.query_occurrence_id = query_occurrence_id
        candidate.target_occurrence_id = target_occurrence_id
        candidate.query_transcript_id = query_transcript_id
        candidate.target_transcript_id = target_transcript_id
        candidate.backend_version = backend_version
        candidate.raw_score = candidate.score
        candidate.nt_identity = candidate.identity
        candidate.known_aligned_pairs = candidate.matches + candidate.mismatches
        candidate.query_covered_bases = sum(block.query.length for block in candidate.aligned_blocks)
        candidate.target_covered_bases = sum(block.target.length for block in candidate.aligned_blocks)
        candidate.query_length = len(query)
        candidate.target_length = len(target)
        candidate.relative_strand = candidate.strand
        candidate.is_secondary = index > 0
        candidate.hit_count = len(candidates)
        candidate.alternative_candidate_ids = tuple(
            candidate_id for candidate_id in candidate_ids if candidate_id != candidate.candidate_id
        )
        candidate.left_anchor_id = left_anchor_id
        candidate.right_anchor_id = right_anchor_id
        candidate.search_interval = bounded_interval
        candidate.enumeration_complete = enumeration_complete
        candidate.incomplete_reason = incomplete_reason
    return AlignmentCandidateSet(
        candidates=candidates,
        enumeration_complete=enumeration_complete,
        incomplete_reason=incomplete_reason,
        score_scheme=NT_BLASTN_V1,
        alignment_mode=mode,
        backend_version=backend_version,
        query_length=len(query),
        target_length=len(target),
        query_occurrence_id=query_occurrence_id,
        target_occurrence_id=target_occurrence_id,
        query_transcript_id=query_transcript_id,
        target_transcript_id=target_transcript_id,
        left_anchor_id=left_anchor_id,
        right_anchor_id=right_anchor_id,
        search_interval=bounded_interval,
    )
