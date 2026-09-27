"""Apply caller-supplied per-tip DNA evidence to a locus-model record."""
from __future__ import annotations

from copy import deepcopy
import math
from collections.abc import Mapping


_REQUIRED = {"family", "unit", "tip", "material_id", "call", "evidence", "survey"}
_OPTIONAL = {"sensitivity", "specificity"}
_DNA_FIELDS = ("surveyed_material", "material_sensitivity",
               "material_specificity", "material_evidence")


def _text(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty string")
    if value != value.strip():
        raise ValueError(f"{name} must not have surrounding whitespace")
    return value


def _probability(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite number in [0, 1]")
    value = float(value)
    if not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError(f"{name} must be a finite number in [0, 1]")
    return value


def apply_material_evidence(model_record, rows):
    """Return a deep-copied model with a complete DNA evidence table applied.

    ``rows`` must cover every declared family/unit/tip/material combination
    exactly once. Calls are ``present``, ``absent`` or ``unknown``. Detection
    probabilities are optional fixed caller-supplied values; this function
    never estimates them or infers homology, roots, or event opportunities.
    Existing DNA observation/detection fields are replaced from the complete
    table; all supplied catalogue, tree, rate, and root records are retained.
    """
    if isinstance(model_record, dict) and model_record.get("schema") == "intraphy.exon-locus-model/1":
        raise ValueError("Schema /1 includes splice and transcript-path state; supply a DNA-only intraphy.exon-locus-model/2 record")
    if not isinstance(model_record, dict) or model_record.get("schema") != "intraphy.exon-locus-model/2" or not isinstance(model_record.get("units"), list):
        raise ValueError("model_record must contain a units array")
    result = deepcopy(model_record)
    expected = {}
    seen_units = set()
    for unit in result["units"]:
        if not isinstance(unit, dict):
            raise ValueError("Every locus unit must be an object")
        family = _text(unit.get("family"), "unit family")
        unit_id = _text(unit.get("unit"), "unit id")
        key = (family, unit_id)
        if key in seen_units:
            raise ValueError(f"Duplicate declared family/unit: {family}/{unit_id}")
        seen_units.add(key)
        catalogue = unit.get("catalogue")
        observations = unit.get("observations")
        if not isinstance(catalogue, dict) or not isinstance(catalogue.get("material"), list):
            raise ValueError(f"Unit {family}/{unit_id} requires a catalogue material array")
        if set(catalogue) != {"provenance", "state_model", "material", "copies", "opportunities"}:
            raise ValueError(f"Unit {family}/{unit_id} requires a DNA-only /2 catalogue")
        material_ids = [_text(item.get("id"), "material id")
                        for item in catalogue["material"] if isinstance(item, dict)]
        if len(material_ids) != len(catalogue["material"]) or len(set(material_ids)) != len(material_ids):
            raise ValueError(f"Unit {family}/{unit_id} has malformed or duplicate material IDs")
        if not isinstance(observations, dict) or not isinstance(observations.get("tips"), dict):
            raise ValueError(f"Unit {family}/{unit_id} requires an observations tips object")
        for tip, tip_record in observations["tips"].items():
            tip = _text(tip, "tip id")
            if not isinstance(tip_record, dict) or not isinstance(tip_record.get("material"), list):
                raise ValueError(f"Tip {tip!r} requires a material observation array")
            if set(tip_record) - set(_DNA_FIELDS) - {"material"}:
                raise ValueError(f"Tip {tip!r} contains unsupported non-DNA observations")
            if len(tip_record["material"]) != len(material_ids):
                raise ValueError(f"Tip {tip!r} material vector does not match its catalogue")
            expected[(family, unit_id, tip)] = material_ids

    if not isinstance(rows, (list, tuple)):
        rows = list(rows)
    by_key = {}
    for row in rows:
        if not isinstance(row, Mapping) or set(row) - (_REQUIRED | _OPTIONAL) or not _REQUIRED <= set(row):
            raise ValueError("Each evidence row requires the declared fields and permits only sensitivity/specificity")
        family = _text(row["family"], "row family")
        unit_id = _text(row["unit"], "row unit")
        tip = _text(row["tip"], "row tip")
        material_id = _text(row["material_id"], "row material_id")
        key = (family, unit_id, tip, material_id)
        if key in by_key:
            raise ValueError(f"Duplicate DNA evidence row: {key}")
        tip_key = (family, unit_id, tip)
        if tip_key not in expected or material_id not in expected[tip_key]:
            raise ValueError(f"Evidence row refers to an undeclared unit, tip, or material: {key}")
        if type(row["survey"]) is not bool:
            raise ValueError("survey must be a Boolean")
        call = row["call"]
        if not isinstance(call, str) or call not in {"present", "absent", "unknown"}:
            raise ValueError("call must be present, absent, or unknown")
        evidence = _text(row["evidence"], "evidence source")
        has_sensitivity, has_specificity = "sensitivity" in row, "specificity" in row
        if has_sensitivity != has_specificity:
            raise ValueError("sensitivity and specificity must be supplied together")
        if call == "unknown" and (has_sensitivity or has_specificity):
            raise ValueError("Unknown calls cannot have detection probabilities")
        if call != "unknown" and (not row["survey"] or not evidence):
            raise ValueError("Known calls require a surveyed material and evidence source")
        parsed = dict(call=call, survey=row["survey"], evidence=evidence)
        if has_sensitivity:
            if not row["survey"] or call == "unknown":
                raise ValueError("Detection probabilities require a known surveyed call")
            parsed["sensitivity"] = _probability(row["sensitivity"], "sensitivity")
            parsed["specificity"] = _probability(row["specificity"], "specificity")
        by_key[key] = parsed

    required_keys = {(family, unit_id, tip, material_id)
                     for (family, unit_id, tip), material_ids in expected.items()
                     for material_id in material_ids}
    if set(by_key) != required_keys:
        missing = sorted(required_keys - set(by_key))
        extra = sorted(set(by_key) - required_keys)
        raise ValueError(f"DNA evidence table must cover every declared row exactly once; missing={missing}, extra={extra}")

    units = {(_text(unit["family"], "unit family"), _text(unit["unit"], "unit id")): unit
             for unit in result["units"]}
    for (family, unit_id, tip), material_ids in expected.items():
        tip_record = units[(family, unit_id)]["observations"]["tips"][tip]
        tip_record["material"] = [None] * len(material_ids)
        surveyed, sensitivity, specificity, evidence = [], {}, {}, {}
        for index, material_id in enumerate(material_ids):
            record = by_key[(family, unit_id, tip, material_id)]
            call = record["call"]
            tip_record["material"][index] = (1 if call == "present" else 0) if call != "unknown" else None
            if record["survey"]:
                surveyed.append(material_id)
            evidence[material_id] = record["evidence"]
            if "sensitivity" in record:
                sensitivity[material_id] = record["sensitivity"]
                specificity[material_id] = record["specificity"]
        for field in _DNA_FIELDS:
            tip_record.pop(field, None)
        tip_record["surveyed_material"] = surveyed
        tip_record["material_evidence"] = evidence
        tip_record["material_sensitivity"] = sensitivity
        tip_record["material_specificity"] = specificity
    return result
