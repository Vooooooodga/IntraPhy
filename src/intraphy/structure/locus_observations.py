"""DNA material observations and optional detection-error likelihoods."""
from __future__ import annotations

import numpy as np

from .locus_events import state_is_valid
from .locus_types import LocusCatalogue, LocusObservation, LocusState


def observation_emission(catalogue: LocusCatalogue, states, observation: LocusObservation,
                         *, material_sensitivity: dict[str, float] | None = None,
                         material_specificity: dict[str, float] | None = None) -> np.ndarray:
    """Return a DNA likelihood vector; an unknown call is uninformative."""
    if not isinstance(catalogue, LocusCatalogue) or not isinstance(observation, LocusObservation):
        raise TypeError("A validated catalogue and LocusObservation are required")
    states = tuple(states)
    if len(observation.material) != len(catalogue.material):
        raise ValueError("Observed material vector does not match catalogue")
    material_ids = {tract.id for tract in catalogue.material}
    if not observation.surveyed_material <= material_ids:
        raise ValueError("Observation surveys an undeclared material tract")
    sensitivity = {} if material_sensitivity is None else dict(material_sensitivity)
    specificity = {} if material_specificity is None else dict(material_specificity)
    if set(sensitivity) != set(specificity):
        raise ValueError("Material sensitivity and specificity must be supplied for the same tracts")
    if not set(sensitivity) <= material_ids:
        raise ValueError("Material detection parameters refer to an undeclared tract")
    if not set(sensitivity) <= observation.surveyed_material:
        raise ValueError("Material detection parameters require explicitly surveyed tracts")
    material_index = {tract.id: i for i, tract in enumerate(catalogue.material)}
    if any(observation.material[material_index[mid]] is None for mid in sensitivity):
        raise ValueError("Material detection parameters require an observed binary result")
    for parameter in (sensitivity, specificity):
        if any(isinstance(value, bool) or not isinstance(value, (int, float)) or not np.isfinite(value)
               or not 0 <= value <= 1 for value in parameter.values()):
            raise ValueError("Detection sensitivity/specificity must lie in [0,1]")
    result = np.ones(len(states), dtype=float)
    for i, state in enumerate(states):
        if not isinstance(state, LocusState) or not state_is_valid(catalogue, state):
            raise ValueError("State vector does not match locus catalogue")
        for k, value in enumerate(observation.material):
            if value is None:
                continue
            tract_id = catalogue.material[k].id
            present = state.material[k] == 1
            if tract_id in sensitivity:
                sens, spec = sensitivity[tract_id], specificity[tract_id]
                result[i] *= (sens if present else 1.0 - spec) if value == 1 else (
                    1.0 - sens if present else spec)
            elif present != (value == 1):
                result[i] = 0.0
                break
    return result


emission_vectors = observation_emission
