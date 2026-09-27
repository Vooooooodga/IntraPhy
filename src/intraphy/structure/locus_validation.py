"""Cross-reference and biological consistency checks for locus catalogues."""
from __future__ import annotations

import math

from .locus_types import LocusCatalogue


def validate_locus_catalogue(catalogue: LocusCatalogue) -> None:
    if not catalogue.provenance.strip():
        raise ValueError("Catalogue provenance is required")
    material_ids = [x.id for x in catalogue.material]
    copy_ids = [x.id for x in catalogue.copies]
    feature_ids = [x.id for x in catalogue.features]
    if (len(set(material_ids)) != len(material_ids) or len(set(copy_ids)) != len(copy_ids)
            or len(set(feature_ids)) != len(feature_ids)):
        raise ValueError("Material, copy, and feature identifiers must be unique")
    if not catalogue.material or not catalogue.copies or not catalogue.features:
        raise ValueError("A locus catalogue requires material, labelled copies, and features")
    ordered = sorted(catalogue.material, key=lambda x: x.start)
    if any(a.end > b.start for a, b in zip(ordered, ordered[1:])):
        raise ValueError("Material tracts must not overlap on the declared coordinate axis")
    mats, copies, features = set(material_ids), set(copy_ids), set(feature_ids)
    material_index = {x.id: x for x in catalogue.material}
    copy_index = {x.id: x for x in catalogue.copies}
    feature_index = {x.id: x for x in catalogue.features}

    owned = {}
    for copy in catalogue.copies:
        if (not set(copy.material_ids) <= mats or not set(copy.homologous_to) <= copies
                or not set(copy.collinear_with) <= copies):
            raise ValueError(f"Copy {copy.id} refers to undeclared material or copy")
        if copy.id in copy.homologous_to or copy.id in copy.collinear_with:
            raise ValueError("A copy cannot be homologous or collinear with itself")
        for material_id in copy.material_ids:
            if material_id in owned:
                raise ValueError("A physical material tract cannot belong to two labelled copy slots")
            owned[material_id] = copy.id

    prerequisites = {f.id: set(f.prerequisites) for f in catalogue.features}
    if any(not deps <= features for deps in prerequisites.values()):
        raise ValueError("Feature prerequisite refers to an undeclared feature")
    resolved = set()
    while True:
        ready = {node for node, deps in prerequisites.items() if node not in resolved and deps <= resolved}
        if not ready:
            break
        resolved.update(ready)
    if resolved != features:
        raise ValueError("Feature prerequisite structure must be acyclic")

    for feature in catalogue.features:
        if not set(feature.required_material) <= mats:
            raise ValueError(f"Feature {feature.id} has undeclared material requirements")
        if feature.copy_id is not None and feature.copy_id not in copies:
            raise ValueError(f"Feature {feature.id} refers to undeclared copy")
        if feature.kind == "exon":
            if not feature.required_material or feature.copy_id is None:
                raise ValueError("An exon feature requires supported material and a labelled copy")
            tract_ids = set(copy_index[feature.copy_id].material_ids)
            intersecting = sorted((material_index[mid] for mid in tract_ids
                                   if material_index[mid].start < feature.end
                                   and feature.start < material_index[mid].end), key=lambda x: x.start)
            cursor = feature.start
            for tract in intersecting:
                if tract.start > cursor:
                    break
                cursor = max(cursor, tract.end)
            if not intersecting or cursor < feature.end:
                raise ValueError(f"Exon feature {feature.id} is not fully covered by its copy material")
            if not set(feature.required_material) <= tract_ids:
                raise ValueError(f"Exon feature {feature.id} requires material outside its copy")
            if any(tract.id not in feature.required_material for tract in intersecting):
                raise ValueError(f"Exon feature {feature.id} must require every material tract it intersects")
        elif feature.kind == "splice":
            if any(feature_index[p].kind != "exon" for p in feature.prerequisites):
                raise ValueError("Splice feature prerequisites must identify exon features")
            upstream = [feature_index[p] for p in feature.prerequisites
                        if feature_index[p].end == feature.donor]
            downstream = [feature_index[p] for p in feature.prerequisites
                          if feature_index[p].start == feature.acceptor]
            if not any(a.id != b.id and a.end == feature.donor < feature.acceptor == b.start
                       for a in upstream for b in downstream):
                raise ValueError("Splice endpoints must connect two ordered prerequisite exon boundaries")

    group_map = dict(catalogue.mxe_groups)
    if len(group_map) != len(catalogue.mxe_groups):
        raise ValueError("MXE group identifiers must be unique")
    for group_id, members in catalogue.mxe_groups:
        if not group_id or not members or len(set(members)) != len(members) or not set(members) <= features:
            raise ValueError("MXE groups require unique declared feature members")
        if any(feature_index[x].kind != "exon" or feature_index[x].mxe_group != group_id for x in members):
            raise ValueError("MXE groups may contain only exon features annotated with that group")
    for feature in catalogue.features:
        if feature.mxe_group is not None and (feature.mxe_group not in group_map
                                              or feature.id not in group_map[feature.mxe_group]):
            raise ValueError("MXE feature annotation must match a declared group membership")

    if len({(e.id, e.outcome_id) for e in catalogue.opportunities}) != len(catalogue.opportunities):
        raise ValueError("Opportunity outcome IDs must be unique within each opportunity")
    by_opportunity = {}
    for event in catalogue.opportunities:
        if not set(event.material_gains + event.material_deletions) <= mats:
            raise ValueError(f"Opportunity {event.id} refers to undeclared material")
        if not set(event.feature_on + event.feature_off + event.context_features +
                   event.required_features + event.forbidden_features + event.graft_features) <= features:
            raise ValueError(f"Opportunity {event.id} refers to undeclared feature")
        if any(a not in features or b not in features for a, b in event.feature_map):
            raise ValueError(f"Opportunity {event.id} has an undeclared feature mapping")
        if ((event.source_copy is not None and event.source_copy not in copies) or
                (event.target_copy is not None and event.target_copy not in copies)):
            raise ValueError(f"Opportunity {event.id} refers to undeclared copy")
        preconditions = dict(event.preconditions)
        if len(preconditions) != len(event.preconditions) or not set(preconditions) <= mats:
            raise ValueError(f"Opportunity {event.id} has invalid material preconditions")
        if set(event.required_features) & set(event.forbidden_features):
            raise ValueError(f"Opportunity {event.id} both requires and forbids one feature")
        by_opportunity.setdefault(event.id, []).append(event)

    for opportunity_id, outcomes in by_opportunity.items():
        common = {(e.kind, e.rate_group, e.preconditions, e.required_features,
                   e.forbidden_features, e.context_features, e.source_copy, e.target_copy)
                  for e in outcomes}
        if len(common) != 1:
            raise ValueError(f"Opportunity {opportunity_id} outcomes must share all applicability conditions")
        if not math.isclose(sum(e.weight for e in outcomes), 1.0, rel_tol=0.0, abs_tol=1e-12):
            raise ValueError(f"Opportunity {opportunity_id} alternative outcome weights must sum to one")
        for event in outcomes:
            _validate_event(event, material_index, copy_index, feature_index)


def _validate_event(event, material_index, copy_index, feature_index):
    if event.kind == "dna_deletion":
        if (event.interval is None or not event.material_deletions or event.material_gains
                or event.feature_on or event.feature_off or event.feature_map or event.graft_features
                or event.source_copy is not None or event.target_copy is not None):
            raise ValueError("DNA deletion requires one physical interval and only its intersected material IDs")
        start, end = event.interval
        intersecting = [m for m in material_index.values() if m.start < end and start < m.end]
        if any(start > m.start or end < m.end for m in intersecting):
            raise ValueError("Split material tracts at deletion boundaries so an interval never deletes outside its extent")
        intersected = {m.id for m in intersecting}
        if not intersected or set(event.material_deletions) != intersected:
            raise ValueError("Deletion must include every supplied material tract intersected by its interval")
    elif event.kind == "copy_duplication":
        if event.interval is not None or event.source_copy is None or event.target_copy is None or event.source_copy == event.target_copy:
            raise ValueError("Copy duplication requires distinct labelled source and target copy slots")
        if event.material_deletions or event.feature_on or event.feature_off:
            raise ValueError("Copy duplication gains target DNA and mapped features; it has no unrelated edits")
        source, target = copy_index[event.source_copy], copy_index[event.target_copy]
        if (source.id not in target.homologous_to and target.id not in source.homologous_to):
            raise ValueError("Duplication requires a declared source-target homology relation")
        if not source.evidence or not target.evidence:
            raise ValueError("Duplication requires evidence on both labelled copy slots")
        if set(source.material_ids) & set(target.material_ids) or set(event.material_gains) != set(target.material_ids):
            raise ValueError("Duplication must introduce the disjoint target copy material exactly")
        preconditions = dict(event.preconditions)
        if any(preconditions.get(mid) != 1 for mid in source.material_ids) or any(preconditions.get(mid) != 0 for mid in target.material_ids):
            raise ValueError("Duplication preconditions must require a present source and unintroduced target slot")
        if not event.feature_map:
            raise ValueError("Duplication requires explicit source-to-target feature mappings")
        mapping = dict(event.feature_map)
        if len(mapping) != len(event.feature_map) or len(set(mapping.values())) != len(mapping):
            raise ValueError("Duplication feature mapping must be one-to-one")
        for src, dst in event.feature_map:
            a, b = feature_index[src], feature_index[dst]
            if (a.copy_id != source.id or b.copy_id != target.id or a.kind != b.kind
                    or not a.evidence or not b.evidence):
                raise ValueError("Duplication mappings require evidenced homologous feature pairs on source/target copies")
            source_deps = {p for p in a.prerequisites if feature_index[p].copy_id == source.id}
            target_deps = {p for p in b.prerequisites if feature_index[p].copy_id == target.id}
            if not source_deps <= mapping.keys() or {mapping[p] for p in source_deps} != target_deps:
                raise ValueError("Duplication feature mapping must preserve copy-local prerequisite structure")
        for feature_id in event.graft_features:
            feature = feature_index[feature_id]
            if feature.copy_id == target.id:
                continue
            target_dependencies = {p for p in feature.prerequisites
                                   if feature_index[p].copy_id == target.id}
            external_dependencies = set(feature.prerequisites) - target_dependencies
            if (feature.kind != "splice" or not target_dependencies or not external_dependencies
                    or not external_dependencies <= set(event.context_features)):
                raise ValueError("A graft outside the target copy must be a splice edge connecting target and declared context features")
    elif event.kind == "splice_change":
        if (event.interval is not None or event.material_gains or event.material_deletions
                or event.source_copy is not None or event.target_copy is not None
                or event.feature_map or event.graft_features):
            raise ValueError("Splice changes alter feature availability while DNA material remains unchanged")
        if not (event.feature_on or event.feature_off) or set(event.feature_on) & set(event.feature_off):
            raise ValueError("Splice changes require distinct feature availability changes")
    else:
        raise ValueError("Supported event kinds are dna_deletion, copy_duplication, and splice_change")
