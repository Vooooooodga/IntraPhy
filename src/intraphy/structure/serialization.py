"""Versioned configuration JSONL, with strict validation and no stale schema alias."""
from __future__ import annotations
from dataclasses import asdict
import json
import math
from pathlib import Path
import numpy as np
from . import MODEL_VERSION, SCHEMA_VERSION
from .types import Catalogue, ExonConfiguration, ExonInstance, ExonSpan, Material, ObservationEvidence, ConfigurationAlternative, InsertionPayload


def json_safe(value):
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, np.generic):
        return json_safe(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (set, frozenset)):
        return sorted(json_safe(v) for v in value)
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    return value


def write_json(path: str | Path, value):
    path = Path(path)
    temporary = path.with_name(path.name+".tmp")
    temporary.write_text(json.dumps(json_safe(value), indent=2, sort_keys=True, allow_nan=False)+"\n")
    temporary.replace(path)


def encode_catalogue(c: Catalogue) -> dict:
    return {"schema": SCHEMA_VERSION, "model": MODEL_VERSION, **asdict(c)}


def decode_catalogue(record: dict) -> Catalogue:
    data = dict(record)
    if data.pop("schema", None) != SCHEMA_VERSION or data.pop("model", None) != MODEL_VERSION:
        raise ValueError("Configuration schema/model mismatch; regenerate 0.19.0 catalogues from original inputs; do not relabel their schema")
    def span(item):
        return ExonSpan(**item) if isinstance(item, dict) else ExonSpan(*item)
    def config(item):
        return ExonConfiguration(tuple(span(e) for e in item["exons"]), tuple(item.get("material", ())))
    data["spans"] = tuple(span(e) for e in data["spans"])
    data["boundary_candidates"] = tuple(span(e) for e in data.get("boundary_candidates", ()))
    data["junctions"] = tuple(tuple(v) for v in data["junctions"])
    data["material"] = tuple(Material(**v) for v in data.get("material", ()))
    data["insertion_payloads"] = tuple(InsertionPayload(
        v["material_id"], tuple(span(e) for e in v["exons"]), v["source_id"], v["evidence"])
        for v in data.get("insertion_payloads", ()))
    observations = []
    for item in data.get("observations", ()):
        o = dict(item)
        o["configurations"] = tuple(config(c) for c in o["configurations"])
        for key in ("unknown_intervals", "alternative_exons"):
            o[key] = tuple(span(e) for e in o.get(key, ()))
        for key in ("material_presence", "reasons", "native_exon_ids"):
            o[key] = tuple(o.get(key, ()))
        o["alternatives"] = tuple(ConfigurationAlternative(
            config(a["configuration"]), a["replaces_key"], a["source_species"],
            a["source_transcript"], a["evidence"], float(a["identity"]))
            for a in o.get("alternatives", ()))
        observations.append(ObservationEvidence(**o))
    data["observations"] = tuple(observations)
    exons = []
    for item in data.get("exon_instances", ()):
        item = dict(item)
        item["transcripts"] = tuple(item["transcripts"])
        item["cds"] = tuple(tuple(p) for p in item.get("cds", ()))
        exons.append(ExonInstance(**item))
    data["exon_instances"] = tuple(exons)
    data["reasons"] = tuple(data.get("reasons", ()))
    return Catalogue(**data)


def write_catalogues(path: str | Path, catalogues):
    from .validation import validate_collection
    catalogues = validate_collection(tuple(catalogues), allow_empty=True)
    path = Path(path)
    temporary = path.with_name(path.name+".tmp")
    with temporary.open("w") as handle:
        for c in catalogues:
            handle.write(json.dumps(encode_catalogue(c), sort_keys=True, allow_nan=False)+"\n")
    temporary.replace(path)


def read_catalogues(path: str | Path, *, allow_empty=False) -> tuple[Catalogue, ...]:
    def unique_keys(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"Duplicate JSON key: {key}")
            result[key] = value
        return result
    records = []
    with Path(path).open() as handle:
        for line_number, line in enumerate(handle, 1):
            if line.strip():
                try:
                    records.append(decode_catalogue(json.loads(line, object_pairs_hook=unique_keys)))
                except (ValueError, TypeError, KeyError) as exc:
                    raise ValueError(f"Invalid configuration at {path}:{line_number}: {exc}") from exc
    from .validation import validate_collection
    return validate_collection(tuple(records), allow_empty=allow_empty)
