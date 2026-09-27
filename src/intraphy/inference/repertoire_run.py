"""Production orchestration for explicitly declared repertoire models."""
from __future__ import annotations
import math
from pathlib import Path
import numpy as np
from ..structure.serialization import write_json, json_safe
from ..structure.serialization import encode_catalogue
from ..storage.tabular import read_tsv
from ..topology import SpeciesTree
from .repertoire_model import evaluate_repertoire_model
from .repertoire_rates import fit_repertoire_scale
from .repertoire_inputs import read_repertoire_model, encode_repertoire

def load_repertoire_bundle(catalogues, species_tree, repertoire_model):
    tree = species_tree if isinstance(species_tree, SpeciesTree) else SpeciesTree(read_tsv(species_tree))
    return read_repertoire_model(repertoire_model, catalogues, tree)

def _unit_inputs(unit, rates, rate_provenance, branch_length_unit, observation_metadata):
    index = unit.process.index
    record = {"family": unit.process.catalogue.family, "unit": unit.process.catalogue.unit,
              "states": [encode_repertoire(s) for s in unit.process.states],
              "events": [{"source": index[e.source], "target": index[e.target], "kind": e.kind, "opportunity": e.opportunity, "weight": e.weight} for e in unit.process.events],
              "rates": dict(rates), "root_prior": np.asarray(unit.root_prior).tolist(),
              "root_prior_provenance": unit.root_prior_provenance, "observation_provenance": unit.observation_provenance,
              "rate_provenance": rate_provenance, "branch_length_unit": branch_length_unit, "process_provenance": unit.process.provenance,
              "tips": {str(k): np.asarray(v).tolist() for k, v in unit.tips.items()},
              "observation_metadata": observation_metadata,
              "catalogue": encode_catalogue(unit.process.catalogue),
              "tree": [{"node_id": n, "parent_id": unit.tree.parent[n], "label": unit.tree.label[n], "branch_length": unit.tree.length[n]} for n in unit.tree.parent],
              "finite_declared_event_table_only": True}
    return record

def _unit_record(unit, result, rates, rate_provenance, branch_length_unit, observation_metadata):
    record = _unit_inputs(unit, rates, rate_provenance, branch_length_unit, observation_metadata)
    record.update({"scale": result["scale"], "nodes": {str(k): np.asarray(v).tolist() for k, v in result["nodes"].items()}, "branches": []})
    for (parent, child), value in result["branches"].items():
        item = {"parent": parent, "child": child}
        item.update({k: (np.asarray(v).tolist() if isinstance(v, np.ndarray) else v) for k, v in value.items()})
        record["branches"].append(item)
    ll = result["log_likelihood"]
    if math.isfinite(ll): record.update(log_likelihood=ll, log_likelihood_status="finite")
    elif ll == -math.inf: record.update(log_likelihood=None, log_likelihood_status="impossible_observation")
    else: raise ArithmeticError("Invalid repertoire log likelihood")
    return record

def _target(output_dir, name):
    target = Path(output_dir) / name
    if target.exists() or target.is_symlink() or target.with_suffix(target.suffix + ".tmp").exists() or target.with_suffix(target.suffix + ".tmp").is_symlink(): raise ValueError(f"Refusing existing output: {target}")
    return target

def infer_repertoires(bundle, output_dir, *, expected_edits=False):
    target = _target(output_dir, "joint_repertoire_history.json")
    records = []
    for unit in bundle.units:
        key = (unit.process.catalogue.family, unit.process.catalogue.unit)
        result = evaluate_repertoire_model(unit.process, unit.tree, unit.tips, rates=bundle.rates[key], root_prior=unit.root_prior, root_prior_provenance=unit.root_prior_provenance, observation_provenance=unit.observation_provenance, rate_provenance=bundle.rate_provenance[key], branch_length_unit=bundle.branch_length_unit, scale=bundle.scale, posterior=True, counts=expected_edits)
        records.append(_unit_record(unit, result, bundle.rates[key], bundle.rate_provenance[key], bundle.branch_length_unit, bundle.observation_metadata.get(key, {})))
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    target = _target(output_dir, "joint_repertoire_history.json")
    write_json(target, json_safe({"schema":"intraphy.exon-repertoire-history/1", "model":"exon-repertoire-ctmc", "units":records, "provenance":bundle.provenance, "model_specification":bundle.model_record}))

def fit_repertoire_rates(bundle, output_dir, *, log_scale_bounds, collection_provenance):
    target = _target(output_dir, "joint_repertoire_rates.json")
    kinds = set().union(*(set(v) for v in bundle.rates.values())); rates = {}
    for kind in kinds:
        values = {r[kind] for r in bundle.rates.values() if kind in r}
        if len(values) != 1: raise ValueError(f"Conflicting relative rate for {kind}")
        rates[kind] = values.pop()
    if bundle.scale != 1.: raise ValueError("Rate fitting requires an unscaled rate template")
    provenance = "; ".join(f"{k[0]}/{k[1]}: {v}" for k, v in sorted(bundle.rate_provenance.items()))
    result = fit_repertoire_scale(bundle.units, rates=rates, rate_provenance=provenance, branch_length_unit=bundle.branch_length_unit, collection_provenance=collection_provenance, log_bounds=tuple(log_scale_bounds))
    target = _target(output_dir, "joint_repertoire_rates.json")
    payload = {"schema":"intraphy.exon-repertoire-rates/1", "model":"exon-repertoire-ctmc", "fit":result, "provenance":bundle.provenance, "model_specification":bundle.model_record, "rate_provenance":[{"family":k[0], "unit":k[1], "provenance":v} for k,v in bundle.rate_provenance.items()]}
    payload["units"] = [_unit_inputs(u, bundle.rates[(u.process.catalogue.family,u.process.catalogue.unit)], bundle.rate_provenance[(u.process.catalogue.family,u.process.catalogue.unit)], bundle.branch_length_unit, bundle.observation_metadata.get((u.process.catalogue.family,u.process.catalogue.unit), {})) for u in bundle.units]
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    write_json(target, json_safe(payload))
