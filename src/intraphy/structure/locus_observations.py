"""Tip evidence constraints for DNA copy lifecycle and splice-feature availability."""
from __future__ import annotations

import numpy as np

from .locus_types import LocusCatalogue, LocusObservation, LocusState


def compatible_splice_path(catalogue: LocusCatalogue, state: LocusState,
                           path: tuple[str, ...]) -> bool:
    """Check one observed directed exon path and its MXE constraint.

    ``path`` lists ordered exon-feature IDs. The state may contain multiple
    available MXE alternatives; a single path can select at most one member
    from each group and each adjacent exon pair needs an active splice edge.
    """
    feature_index = {feature.id: feature for feature in catalogue.features}
    if any(feature_id not in feature_index or feature_index[feature_id].kind != "exon"
           for feature_id in path):
        raise ValueError("A splice path must list declared exon features")
    if len(set(path)) != len(path):
        raise ValueError("A splice path cannot repeat an exon feature")
    if any(feature_id not in state.active_features for feature_id in path):
        return False
    selected_groups = [feature_index[x].mxe_group for x in path if feature_index[x].mxe_group]
    if len(set(selected_groups)) != len(selected_groups):
        return False
    exons = [feature_index[x] for x in path]
    if any(left.start >= right.start for left, right in zip(exons, exons[1:])):
        return False
    splice_edges = [feature for feature in catalogue.features if feature.kind == "splice"
                    and feature.id in state.active_features]
    for left, right in zip(exons, exons[1:]):
        if not any(left.id in edge.prerequisites and right.id in edge.prerequisites
                   and edge.donor == left.end and edge.acceptor == right.start
                   for edge in splice_edges):
            return False
    return True


def observation_emission(catalogue: LocusCatalogue, states, observation: LocusObservation,
                         *, sensitivity: dict[str, float] | None = None,
                         specificity: dict[str, float] | None = None,
                         material_sensitivity: dict[str, float] | None = None,
                         material_specificity: dict[str, float] | None = None) -> np.ndarray:
    """Return a likelihood vector over states without enumerating transcripts.

    DNA observations constrain present versus absent; an observed absence is
    compatible with both unintroduced and deleted latent histories.
    Unsurveyed feature absences and ``None`` values are uninformative. A
    positive feature observation is evidence of availability. Optional
    Bernoulli sensitivity/specificity values apply only to explicitly surveyed
    feature observations. Material detection parameters follow the same rule
    and are keyed by material ID. Both parameter pairs must be supplied
    together for each surveyed binary call. Without DNA detection parameters,
    binary material calls retain their hard-constraint behavior.
    No feature likelihood estimates tissue-specific usage or isoform frequency.
    """
    if not isinstance(catalogue, LocusCatalogue) or not isinstance(observation, LocusObservation):
        raise TypeError("A validated catalogue and LocusObservation are required")
    states = tuple(states)
    if len(observation.material) != len(catalogue.material):
        raise ValueError("Observed material vector does not match catalogue")
    material_ids = {tract.id for tract in catalogue.material}
    if not observation.surveyed_material <= material_ids:
        raise ValueError("Observation surveys an undeclared material tract")
    feature_ids = {feature.id for feature in catalogue.features}
    observed = dict(observation.features)
    if not set(observed) <= feature_ids or not observation.surveyed_features <= feature_ids:
        raise ValueError("Observation refers to an undeclared splice feature")
    for path in observation.observed_paths:
        if any(feature_id not in feature_ids or next(f for f in catalogue.features if f.id == feature_id).kind != "exon"
               for feature_id in path):
            raise ValueError("Observed paths must contain declared exon feature IDs")
    sensitivity = {} if sensitivity is None else dict(sensitivity)
    specificity = {} if specificity is None else dict(specificity)
    material_sensitivity = {} if material_sensitivity is None else dict(material_sensitivity)
    material_specificity = {} if material_specificity is None else dict(material_specificity)
    if set(sensitivity) != set(specificity):
        raise ValueError("Sensitivity and specificity must be supplied for the same features")
    if not (set(sensitivity) | set(specificity)) <= observation.surveyed_features:
        raise ValueError("Detection parameters are permitted only for explicitly surveyed features")
    if set(material_sensitivity) != set(material_specificity):
        raise ValueError("Material sensitivity and specificity must be supplied for the same tracts")
    if not set(material_sensitivity) <= material_ids:
        raise ValueError("Material detection parameters refer to an undeclared tract")
    if not set(material_sensitivity) <= observation.surveyed_material:
        raise ValueError("Material detection parameters require explicitly surveyed tracts")
    material_index = {tract.id: index for index, tract in enumerate(catalogue.material)}
    if any(observation.material[material_index[material_id]] is None
           for material_id in material_sensitivity):
        raise ValueError("Material detection parameters require an observed binary result")
    if any(observed.get(feature_id) is None for feature_id in sensitivity):
        raise ValueError("Detection parameters require an observed binary result")
    for parameter in (sensitivity, specificity, material_sensitivity, material_specificity):
        if any(isinstance(value, bool) or not isinstance(value, (int, float)) or not np.isfinite(value) or not 0 <= value <= 1
               for value in parameter.values()):
            raise ValueError("Detection sensitivity/specificity must lie in [0,1]")

    result = np.ones(len(states), dtype=float)
    for i, state in enumerate(states):
        if not isinstance(state, LocusState) or len(state.material) != len(catalogue.material):
            raise ValueError("State vector does not match locus catalogue")
        for k, value in enumerate(observation.material):
            if value is None:
                continue
            tract_id = catalogue.material[k].id
            present = state.material[k] == 1
            if tract_id in material_sensitivity:
                sens = material_sensitivity[tract_id]
                spec = material_specificity[tract_id]
                probability = (sens if present else 1.0 - spec) if value == 1 else (
                    1.0 - sens if present else spec)
                result[i] *= probability
            elif present != (value == 1):
                result[i] = 0.0
                break
        if result[i] == 0:
            continue
        for feature_id, value in observed.items():
            if value is None:
                continue
            present = feature_id in state.active_features
            if feature_id in observation.surveyed_features and feature_id in sensitivity:
                probability = (sensitivity[feature_id] if present else 1.0 - specificity[feature_id]) if value == 1 else (
                    1.0 - sensitivity[feature_id] if present else specificity[feature_id])
                result[i] *= probability
            elif value == 1:
                if not present:
                    result[i] = 0.0
                    break
            elif feature_id in observation.surveyed_features and present:
                result[i] = 0.0
                break
        if result[i] > 0 and any(not compatible_splice_path(catalogue, state, path)
                                for path in observation.observed_paths):
            result[i] = 0.0
    return result


emission_vectors = observation_emission
