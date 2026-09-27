"""Read the complete exon-locus model and its fixed rooted species tree."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..storage.tabular import read_tsv
from ..topology import SpeciesTree
from ..structure.locus_process import build_locus_process
from .locus_codec import (
    catalogue as decode_catalogue, finite, identifier,
    observations as decode_observations, object_fields, root_distribution, root_states,
    strict_json,
)
from .locus_rates import LocusFitUnit

SCHEMA = "intraphy.exon-locus-model/2"


@dataclass(frozen=True)
class LocusUnitInput:
    family: str
    unit: str
    process: object
    tree: SpeciesTree
    tips: dict
    root_prior: np.ndarray
    root_provenance: str
    observation_provenance: str


@dataclass(frozen=True)
class LocusModelBundle:
    units: tuple[LocusUnitInput, ...]
    tree: SpeciesTree
    tree_source: str
    tree_rows: tuple[dict, ...]
    branch_length_unit: str
    tree_provenance: str
    provenance: str
    independence_provenance: str
    rates: dict[str, dict]
    rate_provenance: str
    model_record: dict

    def fit_units(self):
        return tuple(LocusFitUnit(u.process, u.tree, u.tips, u.root_prior,
                                  name=f"{u.family}/{u.unit}") for u in self.units)


def _tree(path):
    source = Path(path)
    if source.suffix.lower() == ".tsv":
        tree = SpeciesTree(read_tsv(source))
    else:
        from Bio import Phylo
        phylogeny = Phylo.read(str(source), "newick")
        clades = list(phylogeny.find_clades(order="preorder"))
        ids, used, parents = {}, set(), {}
        for parent in clades:
            for child in parent.clades:
                parents[id(child)] = parent
        for index, clade in enumerate(clades):
            label = str(clade.name or "").strip()
            base = label if label and label not in used else f"node_{index}"
            node_id, suffix = base, 1
            while node_id in used:
                node_id = f"{base}_{suffix}"
                suffix += 1
            ids[id(clade)] = node_id
            used.add(node_id)
        rows = []
        for clade in clades:
            parent = parents.get(id(clade))
            node_id = ids[id(clade)]
            rows.append({"node_id": node_id, "parent_id": ids[id(parent)] if parent is not None else "",
                         "label": str(clade.name or node_id),
                         "branch_length": "" if clade.branch_length is None else str(clade.branch_length)})
        tree = SpeciesTree(rows)
    for _, child in tree.edges():
        tree.branch_length(child)
    rows = tuple({"node_id": node, "parent_id": tree.parent[node], "label": tree.label[node],
                  "branch_length": "" if tree.length[node] is None else tree.length[node]}
                 for node in tree.parent)
    return tree, rows


def load_locus_model(model_path, tree_path):
    """Validate all declared objects before command output is reserved."""
    data = strict_json(model_path)
    return locus_model_from_record(data, tree_path)


def locus_model_from_record(data, tree_path):
    """Validate an already-decoded model record against a fixed species tree."""
    if isinstance(data, dict) and data.get("schema") == "intraphy.exon-locus-model/1":
        raise ValueError("Schema /1 includes splice and transcript-path state; supply a DNA-only intraphy.exon-locus-model/2 record")
    fields = {"schema", "model", "branch_length_unit", "tree_provenance", "provenance",
              "independence_provenance", "rates", "units"}
    data = object_fields(data, fields, "locus model")
    if data["schema"] != SCHEMA or data["model"] != "exon-locus-ctmc":
        raise ValueError(f"Expected schema {SCHEMA!r} and model 'exon-locus-ctmc'")
    tree, tree_rows = _tree(tree_path)
    branch_unit = identifier(data["branch_length_unit"], "branch_length_unit")
    tree_provenance = identifier(data["tree_provenance"], "tree_provenance")
    provenance = identifier(data["provenance"], "model provenance")
    independence = identifier(data["independence_provenance"], "independence provenance")
    rate_record = object_fields(data["rates"], {"provenance", "groups"}, "global rates")
    rate_provenance = identifier(rate_record["provenance"], "rate provenance")
    if not isinstance(rate_record["groups"], dict):
        raise ValueError("global rate groups must be an object")
    rates = {}
    for key, raw in rate_record["groups"].items():
        if not isinstance(raw, dict) or raw.get("mode") not in {"fixed", "fit"}:
            raise ValueError(f"Rate group {key!r} must declare mode 'fixed' or 'fit'")
        group = object_fields(raw, {"mode", "value"} if raw["mode"] == "fixed" else {"mode", "initial"}, f"rate group {key}")
        mode = group["mode"]
        value = finite(group["value"] if mode == "fixed" else group["initial"], f"rate group {key}", minimum=0)
        if mode == "fit" and value <= 0:
            raise ValueError(f"Fitted rate group {key!r} requires a positive starting value")
        rates[identifier(key, "rate group")] = {"mode": mode, "value": value}

    if not isinstance(data["units"], list) or not data["units"]:
        raise ValueError("locus model requires at least one unit")
    result, keys, used_groups = [], set(), set()
    family_material_ids, family_intervals = {}, {}
    unit_fields = {"family", "unit", "catalogue", "root", "observations"}
    for raw in data["units"]:
        unit = object_fields(raw, unit_fields, "locus unit")
        family, unit_id = identifier(unit["family"], "family"), identifier(unit["unit"], "unit")
        if (family, unit_id) in keys:
            raise ValueError(f"Duplicate family/unit key: {family}/{unit_id}")
        keys.add((family, unit_id))
        catalogue = decode_catalogue(unit["catalogue"])
        material_ids = {entry.id for entry in catalogue.material}
        prior_ids = family_material_ids.setdefault(family, set())
        overlap = prior_ids & material_ids
        if overlap:
            raise ValueError(f"Independent units in family {family!r} reuse material IDs: {sorted(overlap)}")
        prior_ids.update(material_ids)
        intervals = family_intervals.setdefault(family, [])
        current = [(entry.start, entry.end) for entry in catalogue.material]
        if any(max(a, c) < min(b, d) for a, b in intervals for c, d in current):
            raise ValueError(f"Independent units in family {family!r} overlap on the shared coordinate axis")
        intervals.extend(current)

        root_record = object_fields(unit["root"], {"entries", "provenance"}, "root distribution")
        seeds, _ = root_states(root_record)
        process = build_locus_process(catalogue, seeds)
        root_prior, root_provenance = root_distribution(root_record, process)
        tips, observation_provenance = decode_observations(unit["observations"], catalogue, process, tree)
        groups = {event.rate_group for event in catalogue.opportunities}
        if not groups <= set(rates):
            raise ValueError(f"Rate groups missing values: {sorted(groups - set(rates))}")
        used_groups.update(groups)
        result.append(LocusUnitInput(family, unit_id, process, tree, tips, root_prior,
                                     root_provenance, observation_provenance))
    if used_groups != set(rates):
        raise ValueError(f"Global rates must exactly cover used rate groups; unused={sorted(set(rates) - used_groups)}")
    return LocusModelBundle(tuple(result), tree, str(Path(tree_path).resolve()), tree_rows,
                            branch_unit, tree_provenance, provenance, independence,
                            rates, rate_provenance, data)
