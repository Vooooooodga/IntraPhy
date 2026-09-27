"""Apply declared DNA duplication and continuous-interval deletion outcomes."""
from __future__ import annotations

from .locus_types import EventOpportunity, LocusCatalogue, LocusState


class OpportunityUnavailable(Exception):
    """The whole opportunity is inapplicable in this state."""


def state_is_valid(catalogue: LocusCatalogue, state: LocusState) -> bool:
    return (len(state.material) == len(catalogue.material)
            and all(value in catalogue.material_states for value in state.material))


def apply_outcome(catalogue: LocusCatalogue, source: LocusState,
                  event: EventOpportunity) -> LocusState:
    """Apply one event outcome without renormalizing alternative weights."""
    if not state_is_valid(catalogue, source):
        raise ValueError("Source state is invalid for the selected state model")
    material_index = {tract.id: i for i, tract in enumerate(catalogue.material)}
    for material_id, expected in event.preconditions:
        if source.material[material_index[material_id]] != expected:
            raise OpportunityUnavailable(f"material precondition failed: {material_id}")
    material = list(source.material)
    if event.kind == "dna_deletion":
        present = [mid for mid in event.material_deletions if material[material_index[mid]] == 1]
        if not present:
            raise OpportunityUnavailable("deletion interval contains no present material")
        for material_id in present:
            material[material_index[material_id]] = 0 if catalogue.state_model == "binary" else 2
    elif event.kind == "copy_duplication":
        copy_index = {copy.id: copy for copy in catalogue.copies}
        source_copy = copy_index[event.source_copy]
        target_copy = copy_index[event.target_copy]
        if any(material[material_index[mid]] != 1 for mid in source_copy.material_ids):
            raise OpportunityUnavailable("duplication source copy is not present")
        if any(material[material_index[mid]] != 0 for mid in target_copy.material_ids):
            raise OpportunityUnavailable("duplication target slot is not absent")
        for material_id in event.material_gains:
            material[material_index[material_id]] = 1
    else:
        raise ValueError(f"Unsupported DNA event kind {event.kind!r}")
    return LocusState(tuple(material))
