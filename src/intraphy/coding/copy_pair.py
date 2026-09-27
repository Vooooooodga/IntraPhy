"""Stable public facade for copy-pair keys and immutable evidence artifacts."""
from __future__ import annotations


def expected_copy_pair_keys(occurrences, query_copy_key, target_copy_key):
    """Return every directed exon-like occurrence pair in both orientations."""
    from intraphy.coding.transcripts import _copy_key
    from intraphy.mapping.fields import EXON_LIKE_ROLES

    query_copy_key, target_copy_key = tuple(query_copy_key), tuple(target_copy_key)
    if len(query_copy_key) != 3 or len(target_copy_key) != 3 or query_copy_key == target_copy_key:
        raise ValueError("two distinct copy keys are required")
    if query_copy_key[0] != target_copy_key[0]:
        raise ValueError("copy pair must belong to one family")
    query = sorted(row["occurrence_id"] for row in occurrences
                   if _copy_key(row) == query_copy_key and row.get("role") in EXON_LIKE_ROLES)
    target = sorted(row["occurrence_id"] for row in occurrences
                    if _copy_key(row) == target_copy_key and row.get("role") in EXON_LIKE_ROLES)
    return tuple((q, t) for q in query for t in target) + tuple((t, q) for q in query for t in target)


def expected_all_copy_pair_keys(occurrences):
    """Return exact bidirectional keys for all distinct within-family copy pairs."""
    from intraphy.coding.transcripts import _copy_key

    by_family = {}
    for row in occurrences:
        key = _copy_key(row)
        by_family.setdefault(key[0], set()).add(key)
    result = []
    for family_id in sorted(by_family):
        copies = sorted(by_family[family_id])
        for index, left in enumerate(copies):
            for right in copies[index + 1:]:
                result.extend(expected_copy_pair_keys(occurrences, left, right))
    return tuple(result)


from intraphy.coding.copy_pair_codec import (
    _decode_candidate, _encode_candidate, _decode_evidence, _encode_evidence,
)
from intraphy.coding.copy_pair_artifact import (
    COPY_PAIR_SCHEMA, CopyPairEvidenceProvider, open_copy_pair_evidence,
    write_copy_pair_evidence,
)
from intraphy.coding.family_projection_io import (
    FAMILY_PROJECTION_SCHEMA, read_family_projection, write_family_projection,
)
