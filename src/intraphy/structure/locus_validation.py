"""Cross-reference and biological consistency checks for DNA catalogues."""
from __future__ import annotations

import math

from .locus_types import LocusCatalogue


def validate_locus_catalogue(catalogue: LocusCatalogue) -> None:
    if not isinstance(catalogue.provenance, str) or not catalogue.provenance.strip():
        raise ValueError("Catalogue provenance is required")
    if catalogue.state_model not in {"binary", "irreversible"}:
        raise ValueError("state_model must be binary or irreversible")
    material_ids = [tract.id for tract in catalogue.material]
    copy_ids = [copy.id for copy in catalogue.copies]
    if not material_ids or not copy_ids or len(set(material_ids)) != len(material_ids) or len(set(copy_ids)) != len(copy_ids):
        raise ValueError("A catalogue requires unique material and copy identifiers")
    ordered = sorted(catalogue.material, key=lambda tract: tract.start)
    if any(a.end > b.start for a, b in zip(ordered, ordered[1:])):
        raise ValueError("Material tracts must not overlap on the declared coordinate axis")
    mats, copies = set(material_ids), set(copy_ids)
    material_index = {tract.id: tract for tract in catalogue.material}
    copy_index = {copy.id: copy for copy in catalogue.copies}
    owned = set()
    for copy in catalogue.copies:
        if not set(copy.material_ids) <= mats or not set(copy.homologous_to + copy.collinear_with) <= copies:
            raise ValueError(f"Copy {copy.id} refers to undeclared material or copy")
        if copy.id in copy.homologous_to or copy.id in copy.collinear_with:
            raise ValueError("A copy cannot refer to itself as homologous or collinear")
        if owned & set(copy.material_ids):
            raise ValueError("A physical material tract cannot belong to two labelled copy slots")
        owned.update(copy.material_ids)
    if len({(event.id, event.outcome_id) for event in catalogue.opportunities}) != len(catalogue.opportunities):
        raise ValueError("Opportunity outcome IDs must be unique within each opportunity")
    by_opportunity = {}
    for event in catalogue.opportunities:
        if not set(event.material_gains + event.material_deletions) <= mats:
            raise ValueError(f"Opportunity {event.id} refers to undeclared material")
        if not {mid for mid, _ in event.preconditions} <= mats:
            raise ValueError(f"Opportunity {event.id} has undeclared material preconditions")
        if any(state not in catalogue.material_states for _, state in event.preconditions):
            raise ValueError(f"Opportunity {event.id} has a precondition outside the selected state model")
        if event.source_copy is not None and event.source_copy not in copies or event.target_copy is not None and event.target_copy not in copies:
            raise ValueError(f"Opportunity {event.id} refers to undeclared copy")
        by_opportunity.setdefault(event.id, []).append(event)
        _validate_event(event, material_index, copy_index)
    for opportunity_id, outcomes in by_opportunity.items():
        common = {(e.kind, e.rate_group, e.preconditions, e.source_copy, e.target_copy) for e in outcomes}
        if len(common) != 1:
            raise ValueError(f"Opportunity {opportunity_id} outcomes must share all applicability conditions")
        if not math.isclose(sum(e.weight for e in outcomes), 1.0, rel_tol=0.0, abs_tol=1e-12):
            raise ValueError(f"Opportunity {opportunity_id} alternative outcome weights must sum to one")


def _validate_event(event, material_index, copy_index):
    if event.kind == "dna_deletion":
        if event.interval is None or not event.material_deletions or event.material_gains or event.source_copy is not None or event.target_copy is not None:
            raise ValueError("DNA deletion requires one physical interval and only its intersected material IDs")
        start, end = event.interval
        intersecting = [tract for tract in material_index.values() if tract.start < end and start < tract.end]
        if any(start > tract.start or end < tract.end for tract in intersecting):
            raise ValueError("Split material tracts at deletion boundaries so an interval never deletes outside its extent")
        if not intersecting or set(event.material_deletions) != {tract.id for tract in intersecting}:
            raise ValueError("Deletion must include every supplied material tract intersected by its interval")
    elif event.kind == "copy_duplication":
        if event.interval is not None or event.source_copy is None or event.target_copy is None or event.source_copy == event.target_copy:
            raise ValueError("Copy duplication requires distinct labelled source and target copy slots")
        if event.material_deletions:
            raise ValueError("Copy duplication cannot delete material")
        source, target = copy_index[event.source_copy], copy_index[event.target_copy]
        if source.id not in target.homologous_to and target.id not in source.homologous_to:
            raise ValueError("Duplication requires a declared source-target homology relation")
        if not source.evidence or not target.evidence:
            raise ValueError("Duplication requires evidence on both labelled copy slots")
        if set(source.material_ids) & set(target.material_ids) or set(event.material_gains) != set(target.material_ids):
            raise ValueError("Duplication must introduce the disjoint target copy material exactly")
        preconditions = dict(event.preconditions)
        if any(preconditions.get(mid) != 1 for mid in source.material_ids) or any(preconditions.get(mid) != 0 for mid in target.material_ids):
            raise ValueError("Duplication preconditions must require present source and absent target")
    else:
        raise ValueError("Supported event kinds are dna_deletion and copy_duplication")
