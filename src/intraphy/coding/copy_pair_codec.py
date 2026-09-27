"""Typed JSON codec for candidate and copy-pair evidence records."""
from __future__ import annotations

from dataclasses import fields
import json

from intraphy.coding.types import CodingProjectionCandidate, CodingProjectionCandidateSet
from intraphy.coordinates import CoordinateBlock, Interval0


def _key_json(key):
    if not isinstance(key, (tuple, list)) or len(key) != 4 or any(not isinstance(x, str) or not x for x in key):
        raise ValueError("copy/transcript key must contain four nonempty strings")
    return list(key)


def encode_candidate(candidate):
    if not isinstance(candidate, CodingProjectionCandidate):
        raise TypeError("typed CodingProjectionCandidate required")
    result = {}
    for field in fields(CodingProjectionCandidate):
        value = getattr(candidate, field.name)
        if field.name in {"query_transcript_key", "target_transcript_key"}:
            value = _key_json(value)
        elif field.name == "coordinate_blocks":
            value = [[b.query.start0, b.query.end0, b.target.start0, b.target.end0] for b in value]
        elif field.name == "msa_column_interval":
            value = [value.start0, value.end0]
        elif field.name == "competing_occurrences":
            value = list(value)
        result[field.name] = value
    return result


def decode_candidate(value):
    names = {field.name for field in fields(CodingProjectionCandidate)}
    if not isinstance(value, dict) or set(value) != names:
        raise ValueError("invalid typed candidate record")
    if any(not isinstance(block, list) or len(block) != 4 for block in value["coordinate_blocks"]):
        raise ValueError("invalid candidate coordinate block")
    blocks = tuple(CoordinateBlock(Interval0(*block[:2]), Interval0(*block[2:]))
                   for block in value["coordinate_blocks"])
    interval = Interval0(*value["msa_column_interval"])
    data = dict(value)
    data["query_transcript_key"] = tuple(_key_json(data["query_transcript_key"]))
    data["target_transcript_key"] = tuple(_key_json(data["target_transcript_key"]))
    data["coordinate_blocks"] = blocks
    data["msa_column_interval"] = interval
    data["competing_occurrences"] = tuple(data["competing_occurrences"])
    return CodingProjectionCandidate(**data)


def encode_evidence(evidence):
    if not isinstance(evidence, dict):
        raise TypeError("copy-pair evidence must be a dictionary")
    result = dict(evidence)
    candidate_set = result.pop("protein_candidate_set", None)
    if candidate_set is not None:
        if not isinstance(candidate_set, CodingProjectionCandidateSet):
            raise TypeError("protein_candidate_set must retain its declared typed class")
        result["__candidate_set__"] = [encode_candidate(candidate) for candidate in candidate_set.candidates]
    for key in ("protein_projected_blocks", "protein_projected_coordinate_blocks"):
        if key in result and isinstance(result[key], (tuple, list)):
            result[key] = [[block.query.start0, block.query.end0,
                            block.target.start0, block.target.end0] for block in result[key]]
    return json.dumps(result, sort_keys=True, separators=(",", ":"), allow_nan=False)


def decode_evidence(raw):
    result = json.loads(raw)
    if "__candidate_set__" in result:
        candidates = tuple(decode_candidate(value) for value in result.pop("__candidate_set__"))
        result["protein_candidate_set"] = CodingProjectionCandidateSet(candidates)
    for key in ("protein_projected_blocks", "protein_projected_coordinate_blocks"):
        if key in result and isinstance(result[key], list):
            raw_blocks = result[key]
            if any(not isinstance(block, list) or len(block) != 4 for block in raw_blocks):
                raise ValueError("invalid projected coordinate block")
            result[key] = tuple(CoordinateBlock(Interval0(*block[:2]), Interval0(*block[2:]))
                                for block in raw_blocks)
    return result


_encode_candidate = encode_candidate
_decode_candidate = decode_candidate
_encode_evidence = encode_evidence
_decode_evidence = decode_evidence
