"""Apply one declared copy, material, or splice-feature opportunity outcome."""
from __future__ import annotations

from .locus_types import EventOpportunity, LocusCatalogue, LocusState


class OpportunityUnavailable(Exception):
    """The whole opportunity is inapplicable in this state."""


def state_is_valid(catalogue: LocusCatalogue, state: LocusState) -> bool:
    if len(state.material) != len(catalogue.material):
        return False
    feature_index = {feature.id: feature for feature in catalogue.features}
    material_index = {tract.id: i for i, tract in enumerate(catalogue.material)}
    if not state.active_features <= feature_index.keys():
        return False
    for feature_id in state.active_features:
        feature = feature_index[feature_id]
        if any(state.material[material_index[mid]] != 1 for mid in feature.required_material):
            return False
        if not set(feature.prerequisites) <= state.active_features:
            return False
    return True


def _clean_features(catalogue: LocusCatalogue, material: tuple[int, ...], active: set[str]) -> frozenset[str]:
    """Drop unavailable features and their dependants after a DNA deletion."""
    material_index = {tract.id: i for i, tract in enumerate(catalogue.material)}
    feature_index = {feature.id: feature for feature in catalogue.features}
    changed = True
    while changed:
        changed = False
        for feature_id in tuple(active):
            feature = feature_index[feature_id]
            required = set(feature.required_material)
            if (any(material[material_index[mid]] != 1 for mid in required) or
                    not set(feature.prerequisites) <= active):
                active.remove(feature_id)
                changed = True
    return frozenset(active)


def apply_outcome(catalogue: LocusCatalogue, source: LocusState,
                  event: EventOpportunity) -> LocusState:
    """Apply an outcome, raising if its opportunity is inapplicable.

    Alternative outcomes are handled collectively by ``build_locus_process``;
    this operation never renormalizes their weights.
    """
    material_index = {tract.id: i for i, tract in enumerate(catalogue.material)}
    features = set(source.active_features)
    for material_id, expected in event.preconditions:
        if source.material[material_index[material_id]] != expected:
            raise OpportunityUnavailable(f"material precondition failed: {material_id}")
    if not set(event.required_features) <= features or set(event.forbidden_features) & features:
        raise OpportunityUnavailable("feature precondition failed")
    if not set(event.context_features) <= features:
        raise OpportunityUnavailable("event context feature is unavailable")
    material = list(source.material)

    if event.kind == "dna_deletion":
        present = [mid for mid in event.material_deletions if material[material_index[mid]] == 1]
        if not present:
            raise OpportunityUnavailable("deletion interval contains no present material")
        for material_id in present:
            material[material_index[material_id]] = 2
        features = set(_clean_features(catalogue, tuple(material), features))

    elif event.kind == "copy_duplication":
        source_copy = next(copy for copy in catalogue.copies if copy.id == event.source_copy)
        target_copy = next(copy for copy in catalogue.copies if copy.id == event.target_copy)
        if any(material[material_index[mid]] != 1 for mid in source_copy.material_ids):
            raise OpportunityUnavailable("duplication source copy is not present")
        if any(material[material_index[mid]] != 0 for mid in target_copy.material_ids):
            raise OpportunityUnavailable("duplication target slot is not unintroduced")
        for material_id in event.material_gains:
            material[material_index[material_id]] = 1
        for source_feature, target_feature in event.feature_map:
            if source_feature in features:
                features.add(target_feature)
        feature_index = {feature.id: feature for feature in catalogue.features}
        # Grafts are a conditional projection of the new DNA context. A copy
        # event remains valid when an exon/splice feature lacks active inputs.
        pending = set(event.graft_features)
        while pending:
            ready = {feature_id for feature_id in pending
                     if set(feature_index[feature_id].prerequisites) <= features
                     and all(material[material_index[mid]] == 1
                             for mid in feature_index[feature_id].required_material)}
            if not ready:
                break
            features.update(ready)
            pending.difference_update(ready)
        features = set(_clean_features(catalogue, tuple(material), features))

    elif event.kind == "splice_change":
        feature_index = {feature.id: feature for feature in catalogue.features}
        additions = set(event.feature_on) - features
        removals = set(event.feature_off) & features
        proposed_features = features.difference(removals).union(additions)
        for feature_id in additions:
            feature = feature_index[feature_id]
            if not set(feature.prerequisites) <= proposed_features:
                raise OpportunityUnavailable("splice feature prerequisites are unavailable")
            if any(material[material_index[mid]] != 1 for mid in feature.required_material):
                raise OpportunityUnavailable("splice feature DNA material is unavailable")
        features = proposed_features
        features = set(_clean_features(catalogue, tuple(material), features))

    result = LocusState(tuple(material), frozenset(features))
    if not state_is_valid(catalogue, result):
        raise ValueError(f"Opportunity {event.id}/{event.outcome_id} produces an invalid state")
    return result
