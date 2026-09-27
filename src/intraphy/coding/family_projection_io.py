"""Strict JSON persistence for full-family protein MSA projections."""
from __future__ import annotations

import json
import os
from pathlib import Path

from intraphy.coding.copy_pair_codec import _key_json
from intraphy.coding.types import FamilyCodingProjection


FAMILY_PROJECTION_SCHEMA = "intraphy.family-coding-projection/1"


def _identity(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"nonempty {name} is required")
    return value


def _atomic_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(path)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.partial")
    if temporary.exists():
        raise FileExistsError(temporary)
    try:
        with temporary.open("x", encoding="utf-8") as handle:
            json.dump(payload, handle, sort_keys=True, separators=(",", ":"), allow_nan=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.link(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def write_family_projection(path, projection, *, source_identity, input_identity):
    if not isinstance(projection, FamilyCodingProjection):
        raise TypeError("FamilyCodingProjection required")
    records = {}
    for record_id, sequence in projection.aligned_records.items():
        if not isinstance(record_id, str) or not record_id or not isinstance(sequence, str) or not sequence:
            raise ValueError("family projection contains an empty record")
        aliases = projection.aliases_by_record.get(record_id)
        if not isinstance(aliases, (tuple, list)) or not aliases:
            raise ValueError(f"family projection has no aliases for {record_id}")
        records[record_id] = {"sequence": sequence, "aliases": [_key_json(key) for key in aliases],
                              "residue_columns": list(projection.residue_columns[record_id])}
    if set(projection.aliases_by_record) != set(records) or set(projection.residue_columns) != set(records):
        raise ValueError("family projection record indexes disagree")
    transcript_records = {}
    for key, record_id in projection.record_by_transcript.items():
        key = json.dumps(_key_json(key), separators=(",", ":"))
        if record_id not in records or key in transcript_records:
            raise ValueError("family projection transcript aliases are duplicate or dangling")
        transcript_records[key] = record_id
    payload = {"schema": FAMILY_PROJECTION_SCHEMA, "family_id": projection.family_id,
        "mode": projection.mode, "source_identity": _identity(source_identity, "source_identity"),
        "input_identity": _identity(input_identity, "input_identity"),
        "records": records, "record_by_transcript": transcript_records}
    _atomic_json(path, payload)


def read_family_projection(path, *, expected_family_id, expected_mode, source_identity, input_identity):
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    required = {"schema", "family_id", "mode", "source_identity", "input_identity", "records", "record_by_transcript"}
    if not isinstance(payload, dict) or set(payload) != required or payload["schema"] != FAMILY_PROJECTION_SCHEMA:
        raise ValueError("invalid family projection artifact schema")
    if (payload["family_id"] != expected_family_id or payload["mode"] != expected_mode
            or payload["source_identity"] != _identity(source_identity, "source_identity")
            or payload["input_identity"] != _identity(input_identity, "input_identity")):
        raise ValueError("family projection artifact identity mismatch")
    records = payload["records"]
    aliases_by_record, aligned_records, residue_columns, record_by_transcript = {}, {}, {}, {}
    if not isinstance(records, dict):
        raise ValueError("family projection records must be an object")
    lengths = set()
    for record_id, value in records.items():
        if not isinstance(value, dict) or set(value) != {"sequence", "aliases", "residue_columns"}:
            raise ValueError("invalid family projection record")
        sequence = value["sequence"]
        aliases = tuple(tuple(_key_json(key)) for key in value["aliases"])
        columns = tuple(value["residue_columns"])
        if (not isinstance(sequence, str) or not sequence or not aliases or len(set(aliases)) != len(aliases)
                or len(columns) != len(sequence.replace("-", ""))
                or any(type(col) is not int or col < 0 or col >= len(sequence) for col in columns)
                or tuple(sorted(columns)) != columns):
            raise ValueError("invalid family projection sequence, aliases, or residue columns")
        if columns != tuple(i for i, aa in enumerate(sequence) if aa != "-"):
            raise ValueError("family projection residue columns disagree with aligned sequence")
        lengths.add(len(sequence))
        aligned_records[record_id] = sequence
        aliases_by_record[record_id] = aliases
        residue_columns[record_id] = columns
        for key in aliases:
            if key in record_by_transcript:
                raise ValueError("duplicate transcript alias in family projection")
            record_by_transcript[key] = record_id
    mapped = {}
    if not isinstance(payload["record_by_transcript"], dict):
        raise ValueError("record_by_transcript must be an object")
    for raw_key, record_id in payload["record_by_transcript"].items():
        key = tuple(_key_json(json.loads(raw_key)))
        if record_id not in aligned_records or key in mapped:
            raise ValueError("dangling or duplicate family transcript alias")
        mapped[key] = record_id
    if lengths and len(lengths) != 1:
        raise ValueError("family MSA records have unequal aligned lengths")
    if mapped != record_by_transcript:
        raise ValueError("family transcript aliases disagree with record aliases")
    return FamilyCodingProjection(payload["family_id"], payload["mode"], aligned_records,
        aliases_by_record, mapped, residue_columns)
