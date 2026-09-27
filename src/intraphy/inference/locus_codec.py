"""Strict JSON record decoding for evidence-qualified locus models."""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from ..structure.locus_types import (
    CopySlot, EventOpportunity, LocusCatalogue, LocusObservation,
    LocusState, MaterialTract, SpliceFeature,
)
from ..structure.locus_process import build_locus_process
from ..structure.locus_observations import observation_emission


def object_fields(value, keys, name):
    if not isinstance(value, dict) or set(value) != set(keys):
        raise ValueError(f"{name} must contain exactly these fields: {', '.join(sorted(keys))}")
    return value


def identifier(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty string")
    return value.strip()


def finite(value, name, *, minimum=None):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite number")
    value = float(value)
    if not math.isfinite(value) or (minimum is not None and value < minimum):
        raise ValueError(f"{name} must be finite and at least {minimum}")
    return value


def strict_json(path):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError(f"duplicate JSON key: {key}")
            result[key] = value
        return result

    def parse_float(value):
        result = float(value)
        if not math.isfinite(result):
            raise ValueError("JSON numbers must be finite")
        return result

    return json.loads(Path(path).read_text(encoding="utf-8"),
                      object_pairs_hook=pairs, parse_float=parse_float,
                      parse_constant=lambda x: (_ for _ in ()).throw(ValueError("JSON numbers must be finite")))


def list_value(value, name):
    if not isinstance(value, list):
        raise ValueError(f"{name} must be an array")
    return value


def pair_values(value, name):
    pairs = list_value(value, name)
    if any(not isinstance(pair, list) or len(pair) != 2 for pair in pairs):
        raise ValueError(f"{name} entries must be two-item arrays")
    return tuple(tuple(pair) for pair in pairs)


def catalogue(record):
    allowed = {"provenance", "material", "copies", "features", "opportunities", "mxe_groups"}
    required = {"provenance", "material", "copies", "features", "opportunities"}
    if not isinstance(record, dict) or set(record) - allowed or not required <= set(record):
        raise ValueError("catalogue has missing required or unknown fields")
    if not all(isinstance(record[name], list) for name in ("material", "copies", "features", "opportunities")):
        raise ValueError("catalogue collections must be arrays")
    material = []
    for raw in record["material"]:
        fields = {"id", "start", "end", "evidence"}
        if not isinstance(raw, dict) or set(raw) - fields or not {"id", "start", "end"} <= set(raw):
            raise ValueError("material tract has missing required or unknown fields")
        material.append(MaterialTract(raw["id"], raw["start"], raw["end"],
                                      tuple(list_value(raw.get("evidence", []), "material evidence"))))
    copies = []
    fields = {"id", "material_ids", "homologous_to", "collinear_with", "orientation", "evidence"}
    for raw in record["copies"]:
        if not isinstance(raw, dict) or set(raw) - fields or not {"id", "material_ids"} <= set(raw):
            raise ValueError("copy slot has missing required or unknown fields")
        copies.append(CopySlot(raw["id"], tuple(list_value(raw["material_ids"], "copy material_ids")),
                               tuple(list_value(raw.get("homologous_to", []), "copy homologous_to")),
                               tuple(list_value(raw.get("collinear_with", []), "copy collinear_with")),
                               raw.get("orientation", "+"),
                               tuple(list_value(raw.get("evidence", []), "copy evidence"))))
    features = []
    fields = {"id", "kind", "required_material", "prerequisites", "copy_id", "start", "end", "donor", "acceptor", "mxe_group", "evidence"}
    for raw in record["features"]:
        if not isinstance(raw, dict) or set(raw) - fields or not {"id", "kind"} <= set(raw):
            raise ValueError("splice feature has missing required or unknown fields")
        features.append(SpliceFeature(raw["id"], raw["kind"],
            tuple(list_value(raw.get("required_material", []), "feature required_material")),
            tuple(list_value(raw.get("prerequisites", []), "feature prerequisites")), raw.get("copy_id"),
            raw.get("start"), raw.get("end"), raw.get("donor"), raw.get("acceptor"),
            raw.get("mxe_group"), tuple(list_value(raw.get("evidence", []), "feature evidence"))))
    opportunities = []
    event_fields = {"id", "outcome_id", "kind", "rate_group", "weight", "preconditions", "required_features", "forbidden_features", "material_gains", "material_deletions", "feature_on", "feature_off", "source_copy", "target_copy", "feature_map", "context_features", "graft_features", "interval"}
    defaults = {"weight": 1.0, "preconditions": [], "required_features": [],
                "forbidden_features": [], "material_gains": [], "material_deletions": [],
                "feature_on": [], "feature_off": [], "source_copy": None,
                "target_copy": None, "feature_map": [], "context_features": [],
                "graft_features": [], "interval": None}
    array_fields = ("preconditions", "required_features", "forbidden_features", "material_gains",
                    "material_deletions", "feature_on", "feature_off", "feature_map",
                    "context_features", "graft_features")
    for raw in record["opportunities"]:
        required_event = {"id", "outcome_id", "kind", "rate_group"}
        if not isinstance(raw, dict) or set(raw) - event_fields or not required_event <= set(raw):
            raise ValueError("event opportunity has missing required or unknown fields")
        x = {key: raw.get(key, defaults.get(key)) for key in event_fields}
        for key in array_fields:
            list_value(x[key], f"event opportunity {key}")
        preconditions = pair_values(x["preconditions"], "event preconditions")
        feature_map = pair_values(x["feature_map"], "event feature_map")
        interval = None if x["interval"] is None else tuple(list_value(x["interval"], "event interval"))
        opportunities.append(EventOpportunity(
            id=x["id"], outcome_id=x["outcome_id"], kind=x["kind"], rate_group=x["rate_group"],
            weight=finite(x["weight"], "opportunity weight", minimum=0),
            preconditions=preconditions, required_features=tuple(x["required_features"]),
            forbidden_features=tuple(x["forbidden_features"]), material_gains=tuple(x["material_gains"]),
            material_deletions=tuple(x["material_deletions"]), feature_on=tuple(x["feature_on"]),
            feature_off=tuple(x["feature_off"]), source_copy=x["source_copy"], target_copy=x["target_copy"],
            feature_map=feature_map, context_features=tuple(x["context_features"]),
            graft_features=tuple(x["graft_features"]), interval=interval))
    groups = record.get("mxe_groups", {})
    if not isinstance(groups, dict):
        raise ValueError("mxe_groups must be an object")
    mxe_groups = tuple((key, tuple(list_value(value, f"MXE group {key}")))
                       for key, value in sorted(groups.items()))
    return LocusCatalogue(tuple(material), tuple(copies), tuple(features), tuple(opportunities),
                          mxe_groups, identifier(record["provenance"], "catalogue provenance"))


def root_states(record):
    record = object_fields(record, {"entries", "provenance"}, "root distribution")
    if not isinstance(record["entries"], list) or not record["entries"]:
        raise ValueError("root entries must be a nonempty array")
    states, entries = [], []
    for entry in record["entries"]:
        entry = object_fields(entry, {"state", "weight"}, "root entry")
        if finite(entry["weight"], "root weight", minimum=0) <= 0:
            raise ValueError("Root support entries require positive weights; omit zero-weight states")
        state_record = entry["state"]
        if not isinstance(state_record, dict) or set(state_record) - {"material", "active_features"} or "material" not in state_record:
            raise ValueError("root state requires material and permits only active_features as optional")
        states.append(LocusState(tuple(list_value(state_record["material"], "root material")),
                                 frozenset(list_value(state_record.get("active_features", []), "root active_features"))))
        entries.append(entry)
    return tuple(states), tuple(entries)


def root_distribution(record, process):
    states, entries = root_states(record)
    weights = np.zeros(len(process.states), dtype=float)
    used = set()
    for entry, state in zip(entries, states):
        index = process.index.get(state)
        if index is None or index in used:
            raise ValueError("root state is outside the compiled process or duplicated")
        used.add(index)
        weights[index] = finite(entry["weight"], "root weight", minimum=0)
    if not np.isclose(weights.sum(), 1.0, rtol=0.0, atol=1e-12):
        raise ValueError("Root probabilities must sum to one")
    return weights, identifier(record["provenance"], "root provenance")


def observations(record, catalogue_value, process, tree):
    record = object_fields(record, {"provenance", "tips"}, "observations")
    provenance = identifier(record["provenance"], "observation provenance")
    if not isinstance(record["tips"], dict) or set(record["tips"]) != set(tree.leaf_by_label):
        raise ValueError("Observations must provide exactly one record for every tree tip label")
    tips = {}
    tip_fields = {"material", "features", "surveyed_features", "sensitivity", "specificity",
                  "observed_paths", "surveyed_material", "material_sensitivity",
                  "material_specificity", "material_evidence"}
    for label, raw in record["tips"].items():
        if not isinstance(raw, dict) or set(raw) - tip_fields or "material" not in raw:
            raise ValueError(f"Tip {label} requires material and has unknown fields")
        material = tuple(list_value(raw["material"], f"tip {label} material"))
        if len(material) != len(catalogue_value.material) or any(v not in (None, 0, 1) for v in material):
            raise ValueError(f"Tip {label} material observations must be binary 0/1 or null and match the catalogue")
        feature_values = raw.get("features", {})
        if not isinstance(feature_values, dict) or any(v not in (None, 0, 1) for v in feature_values.values()):
            raise ValueError(f"Tip {label} feature values must be 0/1 or null")
        surveyed = list_value(raw.get("surveyed_features", []), f"tip {label} surveyed_features")
        paths = raw.get("observed_paths", [])
        if not isinstance(paths, list):
            raise ValueError(f"Tip {label} observed_paths must be an array")
        paths = tuple(tuple(list_value(path, f"tip {label} path")) for path in paths)
        sensitivity, specificity = raw.get("sensitivity", {}), raw.get("specificity", {})
        if not isinstance(sensitivity, dict) or not isinstance(specificity, dict):
            raise ValueError(f"Tip {label} sensitivity and specificity must be objects")
        sensitivity = {key: finite(value, "feature sensitivity", minimum=0) for key, value in sensitivity.items()}
        specificity = {key: finite(value, "feature specificity", minimum=0) for key, value in specificity.items()}
        if any(value > 1 for value in (*sensitivity.values(), *specificity.values())):
            raise ValueError("Feature sensitivity and specificity must be in [0, 1]")
        surveyed_material = list_value(raw.get("surveyed_material", []), f"tip {label} surveyed_material")
        if any(not isinstance(value, str) or not value.strip() for value in surveyed_material):
            raise ValueError(f"Tip {label} surveyed_material entries must be nonempty strings")
        if len(set(surveyed_material)) != len(surveyed_material):
            raise ValueError(f"Tip {label} surveyed_material contains duplicates")
        material_sensitivity = raw.get("material_sensitivity", {})
        material_specificity = raw.get("material_specificity", {})
        material_evidence = raw.get("material_evidence", {})
        if not isinstance(material_sensitivity, dict) or not isinstance(material_specificity, dict):
            raise ValueError(f"Tip {label} material sensitivity and specificity must be objects")
        if not isinstance(material_evidence, dict):
            raise ValueError(f"Tip {label} material_evidence must be an object")
        material_sensitivity = {key: finite(value, "material sensitivity", minimum=0)
                                for key, value in material_sensitivity.items()}
        material_specificity = {key: finite(value, "material specificity", minimum=0)
                                for key, value in material_specificity.items()}
        if any(value > 1 for value in (*material_sensitivity.values(), *material_specificity.values())):
            raise ValueError("Material sensitivity and specificity must be in [0, 1]")
        material_ids = {tract.id for tract in catalogue_value.material}
        if not set(material_evidence) <= material_ids:
            raise ValueError(f"Tip {label} material_evidence refers to an undeclared material ID")
        if any(not isinstance(value, str) or not value.strip() for value in material_evidence.values()):
            raise ValueError(f"Tip {label} material_evidence values must be nonempty strings")
        if set(material_sensitivity) != set(material_specificity):
            raise ValueError(f"Tip {label} material sensitivity and specificity must have matching IDs")
        if not set(material_sensitivity) <= set(surveyed_material):
            raise ValueError(f"Tip {label} material detection parameters require surveyed material IDs")
        if any(material_id not in material_evidence for material_id in material_sensitivity):
            raise ValueError(f"Tip {label} material detection parameters require source evidence")
        observation = LocusObservation(material, tuple(sorted(feature_values.items())),
                                       frozenset(surveyed), paths, frozenset(surveyed_material))
        emission = observation_emission(catalogue_value, process.states, observation,
                                        sensitivity=sensitivity or None,
                                        specificity=specificity or None,
                                        material_sensitivity=material_sensitivity or None,
                                        material_specificity=material_specificity or None)
        if not np.any(emission):
            raise ValueError(f"Tip {label!r} has no compatible state in the root-seeded reachable closure")
        tips[label] = emission
    return tips, provenance
