"""mapping / pairwise matches: explicit implementation ownership."""
from __future__ import annotations

from intraphy.aligners.pairwise import local_alignment_stats
from intraphy.aligners.pairwise import overlap_alignment_stats
from intraphy.aligners.short import _validate_max_dp_cells
from intraphy.aligners.short import anchored_short_alignment
from intraphy.aligners.types import MAX_INTERNAL_DP_CELLS
from intraphy.aligners.types import AlignmentStats
from intraphy.coordinates import ClosedInterval1
from intraphy.coordinates import CoordinateBlock
from intraphy.coordinates import Interval0
from intraphy.mapping.candidate_coordinates import _alignment_blocks0
from intraphy.mapping.candidate_coordinates import _covered_bases
from intraphy.mapping.candidate_coordinates import _explicit_bounded_target
from intraphy.mapping.candidate_coordinates import _format_alignment_blocks
from intraphy.mapping.candidate_coordinates import _genomic_matched_blocks
from intraphy.mapping.candidate_coordinates import _valid_local_boundary_range
from intraphy.mapping.candidate_serialization import _alignment_candidate_records
from intraphy.mapping.candidate_serialization import _mapped_genomic_interval0
from intraphy.mapping.candidate_serialization import _set_explicit_alignment_score
from intraphy.mapping.candidate_serialization import _transpose_candidate_record
from intraphy.mapping.candidate_serialization import _transpose_cigar
from intraphy.mapping.fields import EXON_LIKE_ROLES
from intraphy.mapping.match_context import context_score
from intraphy.mapping.match_context import phase_score
from intraphy.mapping.match_context import role_boundary_score
from intraphy.mapping.match_context import segment_length
from intraphy.storage.values import to_float


def match_evidence(left, right, seqs, context, aligner="auto", threads=1, context_aligner="minimap2", short_context_max_length=300, *, short_alignment_max_dp_cells=MAX_INTERNAL_DP_CELLS):
    short_alignment_max_dp_cells = _validate_max_dp_cells(short_alignment_max_dp_cells)
    left_seq = seqs.get(left["occurrence_id"], "")
    right_seq = seqs.get(right["occurrence_id"], "")
    exon_pair = left.get("role") in EXON_LIKE_ROLES and right.get("role") in EXON_LIKE_ROLES
    requested_backend = aligner if exon_pair else context_aligner
    candidate_set = None
    left_is_short_query = (
        len(left_seq) <= int(short_context_max_length)
        and _valid_local_boundary_range(left, left_seq)
        and _explicit_bounded_target(right, right_seq)
    )
    right_is_short_query = (
        not left_is_short_query
        and len(right_seq) <= int(short_context_max_length)
        and _valid_local_boundary_range(right, right_seq)
        and _explicit_bounded_target(left, left_seq)
    )
    bounded_short_context = left_is_short_query or right_is_short_query
    alignment_input_transposed = bool(right_is_short_query)
    bounded_row = left if alignment_input_transposed else right
    # Keep left/query and right/target order for blocks, coverage, and genomic projection.
    if bounded_short_context:
        adapter_query_row, adapter_target_row = (
            (right, left) if alignment_input_transposed else (left, right)
        )
        adapter_query, adapter_target = (
            (right_seq, left_seq) if alignment_input_transposed else (left_seq, right_seq)
        )
        bounded_interval = {
            "coordinate_system": "0-based-half-open",
            "contig": bounded_row.get("contig", "NA"),
            "start0": int(bounded_row["start"]) - 1,
            "end0": int(bounded_row["end"]),
            "strand": bounded_row.get("strand", "NA"),
        }
        candidate_set = anchored_short_alignment(
            adapter_query,
            adapter_target,
            mode="local",
            query_occurrence_id=adapter_query_row.get("occurrence_id"),
            target_occurrence_id=adapter_target_row.get("occurrence_id"),
            query_transcript_id=adapter_query_row.get("transcript_id"),
            target_transcript_id=adapter_target_row.get("transcript_id"),
            search_interval=bounded_interval,
            max_dp_cells=short_alignment_max_dp_cells,
        )
        aln = candidate_set.primary
        requested_backend = "anchored_short_alignment"
        if aln is None:
            aln = AlignmentStats(
                0.0,
                0.0,
                0.0,
                backend="internal",
                alignment_mode="local",
                alignment_meaning="bounded short nucleotide local alignment",
                score_scheme=candidate_set.score_scheme,
                enumeration_complete=candidate_set.enumeration_complete,
                incomplete_reason=candidate_set.incomplete_reason,
            )
    elif exon_pair and aligner in {"mafft", "auto"}:
        aln = overlap_alignment_stats(left_seq, right_seq, backend="mafft", threads=threads)
    else:
        if not exon_pair and context_aligner not in {"internal", "minimap2", "lastz"}:
            raise ValueError("context_aligner must be internal, minimap2, or lastz")
        aln = local_alignment_stats(left_seq, right_seq, backend=requested_backend, threads=threads)
    _set_explicit_alignment_score(aln)
    identity = aln.identity
    aligned_pairs = int(getattr(aln, "aligned_pairs", 0) or 0)
    coverage = aligned_pairs / max(1, min(len(left_seq), len(right_seq)))
    blocks = _alignment_blocks0(aln)
    if alignment_input_transposed:
        blocks = tuple(CoordinateBlock(block.target, block.query) for block in blocks)
    query_interval = (
        Interval0(
            min(block.query.start0 for block in blocks),
            max(block.query.end0 for block in blocks),
        ) if blocks else None
    )
    target_interval = (
        Interval0(
            min(block.target.start0 for block in blocks),
            max(block.target.end0 for block in blocks),
        ) if blocks else None
    )
    public_query = (
        ClosedInterval1.from_interval0(query_interval)
        if query_interval is not None else None
    )
    public_target = (
        ClosedInterval1.from_interval0(target_interval)
        if target_interval is not None else None
    )
    qstart = public_query.start if public_query is not None else "NA"
    qend = public_query.end if public_query is not None else "NA"
    tstart = public_target.start if public_target is not None else "NA"
    tend = public_target.end if public_target is not None else "NA"
    q_gen_start, q_gen_end, q_gen_len = _mapped_genomic_interval0(left, query_interval)
    t_gen_start, t_gen_end, t_gen_len = _mapped_genomic_interval0(right, target_interval)
    score_scheme = getattr(aln, "score_scheme", "unspecified")
    if score_scheme in {None, "", "unspecified"}:
        score_scheme = f"{getattr(aln, 'backend', 'unknown')}_{getattr(aln, 'alignment_mode', 'alignment')}_raw_score"
        aln.score_scheme = score_scheme
    candidate_records = _alignment_candidate_records(aln, candidate_set)
    if alignment_input_transposed:
        candidate_records = [
            _transpose_candidate_record(record) for record in candidate_records
        ]
    for record in candidate_records:
        record["coverage"] = (
            int(record.get("aligned_pairs", 0) or 0)
            / max(1, min(len(left_seq), len(right_seq)))
        )
        if record.get("short_sequence_coverage") in {None, "", "NA"}:
            record["short_sequence_coverage"] = (
                record.get("query_coverage", "NA") if bounded_short_context else "NA"
            )
        record.setdefault("alignment_input_transposed", int(alignment_input_transposed))
        if record.get("query_length") in {None, "", "NA", 0, "0"}:
            record["query_length"] = len(left_seq)
        if record.get("target_length") in {None, "", "NA", 0, "0"}:
            record["target_length"] = len(right_seq)
        if record.get("backend_version") in {None, "", "NA"}:
            record["backend_version"] = getattr(aln, "backend_version", None) or "NA"
        if record.get("raw_score") in {None, "", "NA"}:
            record["raw_score"] = getattr(aln, "raw_score", None)
            if record["raw_score"] is None:
                record["raw_score"] = getattr(aln, "score", "NA")
        if (
            record.get("sequence_kind", "nucleotide") == "nucleotide"
            and record.get("nt_identity") in {None, "", "NA"}
        ):
            record["nt_identity"] = record.get("identity", identity)
        if (
            record.get("sequence_kind") == "amino_acid"
            and record.get("aa_identity") in {None, "", "NA"}
        ):
            record["aa_identity"] = record.get("identity", "NA")
        if record.get("query_covered_bases") in {None, "", "NA", 0, "0"}:
            record["query_covered_bases"] = _covered_bases(record.get("aligned_blocks", ()), "query")
        if record.get("target_covered_bases") in {None, "", "NA", 0, "0"}:
            record["target_covered_bases"] = _covered_bases(record.get("aligned_blocks", ()), "target")
        search_interval = record.get("search_interval")
        if not bounded_short_context:
            record["search_interval"] = "NA"
        elif search_interval is None or search_interval == "" or search_interval == "NA":
            record["search_interval"] = bounded_interval
        record.setdefault(
            "search_interval_side", "query" if alignment_input_transposed else "target",
        )
        if record.get("left_anchor_id") in {None, "", "NA"}:
            record["left_anchor_id"] = "NA"
        if record.get("right_anchor_id") in {None, "", "NA"}:
            record["right_anchor_id"] = "NA"
    primary_record = candidate_records[0] if candidate_records else {}
    unknown_pairs = int(
        getattr(aln, "unknown_aligned_pairs", getattr(aln, "unknown_bases", 0)) or 0
    )
    if candidate_set is not None or getattr(aln, "backend", "") == "internal":
        enumeration_complete = (
            candidate_set.enumeration_complete
            if candidate_set is not None
            else bool(getattr(aln, "enumeration_complete", True))
        )
        enumeration_status = "complete" if enumeration_complete else "incomplete"
        incomplete_reason = (
            candidate_set.incomplete_reason
            if candidate_set is not None
            else getattr(aln, "incomplete_reason", "")
        )
    else:
        enumeration_complete = False
        enumeration_status = "unassessed"
        incomplete_reason = "external_backend_candidate_enumeration_unassessed"
    left_ctx = context.get(left["occurrence_id"], {})
    right_ctx = context.get(right["occurrence_id"], {})
    order = 1.0 - abs(to_float(left_ctx.get("scaled_index"), 0.5) - to_float(right_ctx.get("scaled_index"), 0.5))
    left_context = context_score(left_ctx, right_ctx, "left")
    right_context = context_score(left_ctx, right_ctx, "right")
    boundary = role_boundary_score(left, right)
    phase = phase_score(left, right)
    alignment_strand = getattr(aln, "strand", "+")
    strand = 1.0 if alignment_strand == "+" else 0.0
    splice = 1.0 - abs(to_float(left.get("splice_motif_score"), 0.5) - to_float(right.get("splice_motif_score"), 0.5))
    sequence_score = 0.70 * identity + 0.30 * coverage
    structural_context_score = 0.5 * left_context + 0.5 * right_context
    total = (
        0.34 * identity
        + 0.14 * coverage
        + 0.10 * left_context
        + 0.10 * right_context
        + 0.10 * boundary
        + 0.08 * phase
        + 0.06 * order
        + 0.04 * strand
        + 0.04 * splice
    )
    return {
        "alignment_score": identity,
        "coverage_score": coverage,
        "sequence_score": sequence_score,
        "structural_context_score": structural_context_score,
        "query_coverage": aln.target_coverage if alignment_input_transposed else aln.query_coverage,
        "target_coverage": aln.query_coverage if alignment_input_transposed else aln.target_coverage,
        "aligned_pairs": aligned_pairs,
        "sequence_kind": primary_record.get("sequence_kind", "nucleotide"),
        "backend": primary_record.get("backend", getattr(aln, "backend", "NA")),
        "backend_version": primary_record.get("backend_version", getattr(aln, "backend_version", "NA")),
        "raw_score": primary_record.get("raw_score", getattr(aln, "score", "NA")),
        "nt_identity": primary_record.get("nt_identity", identity),
        "aa_identity": primary_record.get("aa_identity", "NA"),
        "known_aligned_pairs": primary_record.get(
            "known_aligned_pairs", getattr(aln, "known_aligned_pairs", aligned_pairs),
        ),
        "unknown_aligned_pairs": primary_record.get("unknown_aligned_pairs", unknown_pairs),
        "query_covered_bases": primary_record.get("query_covered_bases", _covered_bases(blocks, "query")),
        "target_covered_bases": primary_record.get("target_covered_bases", _covered_bases(blocks, "target")),
        "query_length": len(left_seq),
        "target_length": len(right_seq),
        "relative_strand": primary_record.get("relative_strand", getattr(aln, "strand", "+")),
        "gap_blocks": primary_record.get("gap_blocks", []),
        "search_interval": bounded_interval if bounded_short_context else "NA",
        "search_interval_side": (
            "query" if alignment_input_transposed else "target"
        ) if bounded_short_context else "NA",
        "alignment_input_transposed": int(alignment_input_transposed),
        "short_sequence_coverage": primary_record.get(
            "short_sequence_coverage",
            aln.query_coverage if bounded_short_context else "NA",
        ),
        "query_alignment_start": qstart,
        "query_alignment_end": qend,
        "target_alignment_start": tstart,
        "target_alignment_end": tend,
        "query_mapped_contig": left.get("contig", "NA"),
        "query_mapped_start": q_gen_start,
        "query_mapped_end": q_gen_end,
        "query_mapped_strand": left.get("strand", "NA"),
        "query_mapped_length": q_gen_len,
        "subject_mapped_contig": right.get("contig", "NA"),
        "subject_mapped_start": t_gen_start,
        "subject_mapped_end": t_gen_end,
        "subject_mapped_strand": right.get("strand", "NA"),
        "subject_mapped_length": t_gen_len,
        "alignment_strand": alignment_strand,
        "projected_reference_occurrence_id": right["occurrence_id"],
        "projected_reference_start": tstart,
        "projected_reference_end": tend,
        "projected_reference_blocks": _format_alignment_blocks(blocks),
        "matched_blocks": _format_alignment_blocks(blocks),
        "query_genomic_matched_blocks": _genomic_matched_blocks(left, blocks, "query"),
        "subject_genomic_matched_blocks": _genomic_matched_blocks(right, blocks, "subject"),
        "query_parent_feature_ids": left.get("source_feature_id", "NA"),
        "subject_parent_feature_ids": right.get("source_feature_id", "NA"),
        "query_transcript_ids": left.get("transcript_id", "NA"),
        "subject_transcript_ids": right.get("transcript_id", "NA"),
        "left_context_score": left_context,
        "right_context_score": right_context,
        "boundary_score": boundary,
        "phase_score": phase,
        "order_score": order,
        "strand_score": strand,
        "splice_score": splice,
        "size_ratio": min(segment_length(left), segment_length(right)) / max(segment_length(left), segment_length(right)),
        "alignment_cigar": _transpose_cigar(aln.cigar) if alignment_input_transposed else aln.cigar,
        "alignment_backend": aln.backend,
        "alignment_mode": getattr(aln, "alignment_mode", "local"),
        "alignment_meaning": getattr(aln, "alignment_meaning", "bounded short nucleotide local alignment"),
        "alignment_requested_backend": requested_backend,
        "mapping_quality": (
            getattr(aln, "mapping_quality", "NA")
            if getattr(aln, "backend", "internal") not in {"internal", "mafft"}
            else "NA"
        ),
        "hit_count": (
            len(candidate_set.candidates)
            if candidate_set is not None
            else getattr(aln, "hit_count", 0)
        ),
        "ambiguous_hit_count": (
            max(0, len(candidate_set.candidates) - 1)
            if candidate_set is not None
            else getattr(aln, "ambiguous_hit_count", 0)
        ),
        "alternative_hits": candidate_records[1:],
        "score_scheme": getattr(aln, "score_scheme", "unspecified"),
        "raw_alignment_score": getattr(aln, "score", "NA"),
        "enumeration_complete": enumeration_complete,
        "candidate_enumeration_status": enumeration_status,
        "incomplete_reason": incomplete_reason or "NA",
        "short_context_route": "bounded_local" if bounded_short_context else "not_used",
        "local_boundary_range": bounded_interval if bounded_short_context else "not_evaluated",
        "flanking_anchor_status": "not_evaluated_pre_chain",
        "true_absence_eligible": 0,
        "true_absence_evidence_status": "insufficient_evidence",
        "true_absence_reason": "ordered_double_flanks_not_evaluated;anchor_interval_sequence_not_extracted;assembly_continuity_unassessed;ambiguous_base_status_unassessed;query_only_deletion_gap_unassessed;alternative_alignment_concordance_unassessed",
        "candidate_records": candidate_records,
        "total_score": total,
    }
