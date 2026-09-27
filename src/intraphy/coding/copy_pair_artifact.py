"""Immutable SQLite writer and read-only provider for pairwise evidence."""
from __future__ import annotations

import json
import os
from pathlib import Path
import sqlite3
from collections.abc import Mapping

from intraphy.coding.copy_pair_codec import decode_evidence, encode_evidence

COPY_PAIR_SCHEMA = "intraphy.copy-pair-evidence/1"


def _identity(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"nonempty {name} is required")
    return value


def _copy_key_json(key):
    if not isinstance(key, (tuple, list)) or len(key) != 3 or any(not isinstance(x, str) or not x for x in key):
        raise ValueError("copy key must contain three nonempty strings")
    return list(key)


def write_copy_pair_evidence(path, *, family_id, query_copy_key, target_copy_key,
                             source_identity, input_identity, expected_pairs, records):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(path)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.partial")
    if temporary.exists():
        raise FileExistsError(temporary)
    left, right = tuple(_copy_key_json(query_copy_key)), tuple(_copy_key_json(target_copy_key))
    if left == right or left[0] != family_id or right[0] != family_id:
        raise ValueError("copy-pair keys must be distinct members of the declared family")
    source_identity = _identity(source_identity, "source_identity")
    input_identity = _identity(input_identity, "input_identity")
    connection = None
    try:
        connection = sqlite3.connect(temporary)
        connection.execute("PRAGMA journal_mode=DELETE")
        connection.execute("PRAGMA synchronous=FULL")
        connection.execute("CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL) WITHOUT ROWID")
        connection.execute("CREATE TABLE expected (query_occurrence_id TEXT NOT NULL, target_occurrence_id TEXT NOT NULL, PRIMARY KEY(query_occurrence_id,target_occurrence_id)) WITHOUT ROWID")
        connection.execute("CREATE TABLE evidence (query_occurrence_id TEXT NOT NULL, target_occurrence_id TEXT NOT NULL, payload TEXT NOT NULL, PRIMARY KEY(query_occurrence_id,target_occurrence_id)) WITHOUT ROWID")
        metadata = {"schema": COPY_PAIR_SCHEMA, "status": "complete", "family_id": family_id,
                    "query_copy_key": json.dumps(left, separators=(",", ":")),
                    "target_copy_key": json.dumps(right, separators=(",", ":")),
                    "source_identity": source_identity, "input_identity": input_identity}
        connection.executemany("INSERT INTO metadata(key,value) VALUES (?,?)", sorted(metadata.items()))
        expected_rows = []
        for key in expected_pairs:
            if not isinstance(key, (tuple, list)) or len(key) != 2 or any(not isinstance(x, str) or not x for x in key):
                raise ValueError("expected directed occurrence pair keys must be two nonempty strings")
            expected_rows.append((key[0], key[1]))
        if len(expected_rows) != len(set(expected_rows)):
            raise ValueError("duplicate expected directed occurrence pair")
        connection.executemany("INSERT INTO expected VALUES (?,?)", expected_rows)
        count = 0
        try:
            for key, evidence in records:
                if not isinstance(key, (tuple, list)) or len(key) != 2 or any(not isinstance(x, str) or not x for x in key):
                    raise ValueError("evidence key must be a directed occurrence pair")
                connection.execute("INSERT INTO evidence VALUES (?,?,?)",
                                   (key[0], key[1], encode_evidence(evidence)))
                count += 1
        except Exception:
            close = getattr(records, "close", None)
            if close is not None:
                close()
            raise
        expected_count = len(expected_rows)
        if count != expected_count:
            raise ValueError(f"copy-pair evidence count mismatch: expected {expected_count}, observed {count}")
        missing, extra = connection.execute(
            "SELECT (SELECT COUNT(*) FROM expected e LEFT JOIN evidence v USING(query_occurrence_id,target_occurrence_id) WHERE v.query_occurrence_id IS NULL), "
            "(SELECT COUNT(*) FROM evidence v LEFT JOIN expected e USING(query_occurrence_id,target_occurrence_id) WHERE e.query_occurrence_id IS NULL)"
        ).fetchone()
        if missing or extra:
            raise ValueError(f"copy-pair evidence key mismatch: missing={missing}, unexpected={extra}")
        connection.commit()
        connection.close()
        connection = None
        with temporary.open("rb") as handle:
            os.fsync(handle.fileno())
        os.link(temporary, path)
    finally:
        if connection is not None:
            connection.close()
        temporary.unlink(missing_ok=True)


class CopyPairEvidenceProvider:
    """Read-only indexed lookup over validated immutable copy-pair artifacts."""
    def __init__(self, connections, key_to_connection):
        self._connections = tuple(connections)
        self._key_to_connection = key_to_connection

    def get(self, query_occurrence_id, target_occurrence_id):
        key = (query_occurrence_id, target_occurrence_id)
        connection = self._key_to_connection.get(key)
        if connection is None:
            raise KeyError(key)
        row = connection.execute(
            "SELECT payload FROM evidence WHERE query_occurrence_id=? AND target_occurrence_id=?", key
        ).fetchone()
        if row is None:
            raise KeyError(key)
        return decode_evidence(row[0])

    @property
    def keys(self):
        return frozenset(self._key_to_connection)

    def close(self):
        for connection in self._connections:
            connection.close()


def open_copy_pair_evidence(paths, *, expected_source_identity, expected_input_identity, expected_pairs):
    paths = tuple(paths)
    expected = set()
    for key in expected_pairs:
        if not isinstance(key, (tuple, list)) or len(key) != 2 or any(not isinstance(x, str) or not x for x in key):
            raise ValueError("expected directed occurrence pair keys must be two nonempty strings")
        normalized = (key[0], key[1])
        if normalized in expected:
            raise ValueError(f"duplicate expected directed evidence key: {normalized}")
        expected.add(normalized)
    connections, key_to_connection = [], {}
    observed_copy_pairs = set()
    try:
        for path in paths:
            uri = Path(path).resolve().as_uri() + "?mode=ro&immutable=1"
            connection = sqlite3.connect(uri, uri=True)
            connections.append(connection)
            metadata = dict(connection.execute("SELECT key,value FROM metadata"))
            required = {"schema", "status", "family_id", "query_copy_key", "target_copy_key", "source_identity", "input_identity"}
            if set(metadata) != required or metadata["schema"] != COPY_PAIR_SCHEMA or metadata["status"] != "complete":
                raise ValueError(f"invalid or incomplete copy-pair artifact: {path}")
            query_copy = tuple(_copy_key_json(json.loads(metadata["query_copy_key"])))
            target_copy = tuple(_copy_key_json(json.loads(metadata["target_copy_key"])))
            if query_copy == target_copy or query_copy[0] != metadata["family_id"] or target_copy[0] != metadata["family_id"]:
                raise ValueError(f"invalid copy-pair artifact metadata: {path}")
            if isinstance(expected_input_identity, Mapping):
                pair = (query_copy, target_copy)
                if pair not in expected_input_identity:
                    raise ValueError(f"no expected input identity for copy pair: {pair}")
                expected_pair_identity = expected_input_identity[pair]
            else:
                if len(paths) != 1:
                    raise ValueError("per-copy-pair expected_input_identity mapping required for multiple artifacts")
                expected_pair_identity = expected_input_identity
            if (metadata["source_identity"] != _identity(expected_source_identity, "source_identity")
                    or metadata["input_identity"] != _identity(expected_pair_identity, "input_identity")):
                raise ValueError(f"copy-pair artifact identity mismatch: {path}")
            copy_pair = (query_copy, target_copy)
            if copy_pair in observed_copy_pairs:
                raise ValueError(f"duplicate copy-pair artifact: {copy_pair}")
            observed_copy_pairs.add(copy_pair)
            artifact_expected = {(q, t) for q, t in connection.execute(
                "SELECT query_occurrence_id,target_occurrence_id FROM expected")}
            stored = {(q, t) for q, t in connection.execute(
                "SELECT query_occurrence_id,target_occurrence_id FROM evidence")}
            if stored != artifact_expected or not stored <= expected:
                raise ValueError(f"copy-pair artifact key mismatch: {path}")
            for key in stored:
                if key in key_to_connection:
                    raise ValueError(f"duplicate directed evidence key across artifacts: {key}")
                key_to_connection[key] = connection
        if set(key_to_connection) != expected:
            missing = len(expected - set(key_to_connection))
            unexpected = len(set(key_to_connection) - expected)
            raise ValueError(f"copy-pair provider key mismatch: missing={missing}, unexpected={unexpected}")
        if isinstance(expected_input_identity, Mapping) and set(expected_input_identity) != observed_copy_pairs:
            raise ValueError("expected input identity copy-pair mapping does not match artifact set")
        return CopyPairEvidenceProvider(connections, key_to_connection)
    except Exception:
        for connection in connections:
            connection.close()
        raise
