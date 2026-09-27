"""mapping / candidate records: explicit implementation ownership."""
from __future__ import annotations

from intraphy.coordinates import ClosedInterval1
from intraphy.mapping.candidate_coordinates import _alignment_blocks0
from intraphy.mapping.candidate_coordinates import _alignment_gap_blocks
from intraphy.mapping.candidate_coordinates import _candidate_value
from intraphy.mapping.candidate_coordinates import _coordinate_block0
from intraphy.mapping.candidate_coordinates import _covered_bases
from intraphy.mapping.candidate_coordinates import _unknown_pair_count


def _scope_candidate_record_ids(records, match_id, namespace):
    """Give each record a row-unique ID while retaining backend IDs."""
    source_ids = []
    seen = set()
    for record in records:
        source_id = record.get("candidate_id")
        if source_id in {None, "", "NA"}:
            source_ids.append(None)
            continue
        source_id = str(source_id)
        if source_id in seen:
            raise ValueError(
                f"ambiguous duplicate source candidate ID within {namespace} "
                f"for {match_id}: {source_id}"
            )
        seen.add(source_id)
        source_ids.append(source_id)

    scoped_by_source = {
        source_id: f"{match_id}.{namespace}.{source_id}"
        for source_id in source_ids if source_id is not None
    }
    for index, (record, source_id) in enumerate(zip(records, source_ids), start=1):
        record["source_candidate_id"] = source_id or "NA"
        record["candidate_id"] = (
            scoped_by_source[source_id]
            if source_id is not None
            else f"{match_id}.{namespace}_{index:03d}"
        )
        alternatives = record.get("alternative_candidate_ids")
        if isinstance(alternatives, str):
            if alternatives in {"", "NA"}:
                continue
            tokens = alternatives.split(";")
            record["alternative_candidate_ids"] = ";".join(
                scoped_by_source.get(token, token) for token in tokens
            )
        elif alternatives is not None:
            remapped = [
                scoped_by_source.get(str(candidate_id), candidate_id)
                for candidate_id in alternatives
            ]
            record["alternative_candidate_ids"] = (
                tuple(remapped) if isinstance(alternatives, tuple) else remapped
            )


def _candidate_record(candidate, rank, fallback_backend, fallback_scheme):
    if isinstance(candidate, dict):
        record = dict(candidate)
        blocks = tuple(_coordinate_block0(block) for block in record.get("aligned_blocks", ()))
        record["aligned_blocks"] = blocks
        if blocks:
            record["query_start0"] = min(block.query.start0 for block in blocks)
            record["query_end0"] = max(block.query.end0 for block in blocks)
            record["target_start0"] = min(block.target.start0 for block in blocks)
            record["target_end0"] = max(block.target.end0 for block in blocks)
        elif record.get("query_start") not in {None, "", "NA"}:
            query = ClosedInterval1(
                int(record["query_start"]), int(record["query_end"]),
            ).to_interval0()
            target = ClosedInterval1(
                int(record["target_start"]), int(record["target_end"]),
            ).to_interval0()
            record["query_start0"], record["query_end0"] = query.start0, query.end0
            record["target_start0"], record["target_end0"] = target.start0, target.end0
        for field in ("query_start", "query_end", "target_start", "target_end"):
            record.pop(field, None)
        record.setdefault("rank", rank)
        record.setdefault("backend", fallback_backend)
        record.setdefault("query_coverage", record.get("coverage", "NA"))
        record.setdefault("target_coverage", record.get("coverage", "NA"))
        record.setdefault("aligned_pairs", sum(block.query.length for block in blocks))
        record["gap_blocks"] = _alignment_gap_blocks(record)
        if record.get("score_scheme") in {None, "", "unspecified"}:
            record["score_scheme"] = fallback_scheme
        sequence_kind = record.setdefault("sequence_kind", "nucleotide")
        record.setdefault("backend_version", "NA")
        record.setdefault("raw_score", record.get("score", "NA"))
        record.setdefault("nt_identity", record.get("identity", "NA") if sequence_kind == "nucleotide" else "NA")
        record.setdefault("aa_identity", record.get("identity", "NA") if sequence_kind == "amino_acid" else "NA")
        unknown = _unknown_pair_count(record)
        record.setdefault("unknown_aligned_pairs", unknown)
        record.setdefault("known_aligned_pairs", int(record.get("aligned_pairs", 0) or 0))
        record.setdefault("query_covered_bases", _covered_bases(blocks, "query"))
        record.setdefault("target_covered_bases", _covered_bases(blocks, "target"))
        record.setdefault("query_length", "NA")
        record.setdefault("target_length", "NA")
        record.setdefault("relative_strand", record.get("strand", "+"))
        record.setdefault("mapq", record.get("mapping_quality", "NA"))
        if record.get("mapq") is None:
            record["mapq"] = "NA"
        record.setdefault("search_interval_side", "target")
        return record
    adapter_fields = dict(vars(candidate)) if hasattr(candidate, "__dict__") else {}
    for field in (
        "aligned_blocks", "gap_blocks", "alternative_hits",
        "query_interval", "target_interval",
    ):
        adapter_fields.pop(field, None)
    blocks = _alignment_blocks0(candidate)
    aligned_pairs = int(_candidate_value(candidate, "aligned_pairs", 0) or 0)
    unknown_pairs = _unknown_pair_count(candidate)
    sequence_kind = _candidate_value(candidate, "sequence_kind", "nucleotide")
    mapq = _candidate_value(
        candidate, "mapq", _candidate_value(candidate, "mapping_quality", "NA"),
    )
    if mapq is None:
        mapq = "NA"
    return {
        **adapter_fields,
        "candidate_id": _candidate_value(candidate, "candidate_id", "NA"),
        "rank": rank,
        "identity": _candidate_value(candidate, "identity", "NA"),
        "coverage": _candidate_value(candidate, "coverage", "NA"),
        "query_coverage": _candidate_value(candidate, "query_coverage", "NA"),
        "target_coverage": _candidate_value(candidate, "target_coverage", "NA"),
        "aligned_pairs": aligned_pairs,
        "query_start0": min((block.query.start0 for block in blocks), default="NA"),
        "query_end0": max((block.query.end0 for block in blocks), default="NA"),
        "target_start0": min((block.target.start0 for block in blocks), default="NA"),
        "target_end0": max((block.target.end0 for block in blocks), default="NA"),
        "strand": _candidate_value(candidate, "strand", "+"),
        "mapping_quality": _candidate_value(candidate, "mapping_quality", "NA"),
        "is_secondary": int(bool(_candidate_value(candidate, "is_secondary", rank > 1))),
        "score": _candidate_value(candidate, "score", "NA"),
        "cigar": _candidate_value(candidate, "cigar", "NA"),
        "aligned_blocks": blocks,
        "gap_blocks": _alignment_gap_blocks(candidate),
        "sequence_kind": sequence_kind,
        "backend": _candidate_value(candidate, "backend", fallback_backend),
        "backend_version": _candidate_value(candidate, "backend_version", "NA"),
        "score_scheme": _candidate_value(candidate, "score_scheme", fallback_scheme),
        "raw_score": _candidate_value(candidate, "raw_score", _candidate_value(candidate, "score", "NA")),
        "nt_identity": _candidate_value(
            candidate,
            "nt_identity",
            _candidate_value(candidate, "identity", "NA") if sequence_kind == "nucleotide" else "NA",
        ),
        "aa_identity": _candidate_value(
            candidate,
            "aa_identity",
            _candidate_value(candidate, "identity", "NA") if sequence_kind == "amino_acid" else "NA",
        ),
        "known_aligned_pairs": _candidate_value(
            candidate, "known_aligned_pairs", aligned_pairs,
        ),
        "unknown_aligned_pairs": unknown_pairs,
        "query_covered_bases": _candidate_value(
            candidate, "query_covered_bases", _covered_bases(blocks, "query"),
        ),
        "target_covered_bases": _candidate_value(
            candidate, "target_covered_bases", _covered_bases(blocks, "target"),
        ),
        "query_length": _candidate_value(candidate, "query_length", "NA"),
        "target_length": _candidate_value(candidate, "target_length", "NA"),
        "relative_strand": _candidate_value(
            candidate, "relative_strand", _candidate_value(candidate, "strand", "+"),
        ),
        "mapq": mapq,
        "left_anchor_id": _candidate_value(candidate, "left_anchor_id", "NA"),
        "right_anchor_id": _candidate_value(candidate, "right_anchor_id", "NA"),
        "search_interval": _candidate_value(candidate, "search_interval", "NA"),
        "search_interval_side": "target",
        "enumeration_complete": _candidate_value(
            candidate, "enumeration_complete", "NA",
        ),
        "incomplete_reason": _candidate_value(candidate, "incomplete_reason", "NA"),
    }
