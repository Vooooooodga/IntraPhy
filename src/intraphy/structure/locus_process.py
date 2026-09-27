"""Compile the complete reachable closure of a finite evidence-conditioned locus."""
from __future__ import annotations

from collections import defaultdict, deque

from .locus_events import OpportunityUnavailable, apply_outcome, state_is_valid
from .locus_types import LocusCatalogue, LocusProcess, LocusState, ProcessEdge


def build_locus_process(catalogue: LocusCatalogue, root_states) -> LocusProcess:
    """Build all states reachable from explicitly supplied root support.

    The finite closure is conditional on the supplied labelled DNA catalogue
    and event opportunities. There is no default state cap or tip-seeded state
    discovery. A partial set of applicable alternatives fails rather than
    redistributing outcome weights.
    """
    if not isinstance(catalogue, LocusCatalogue):
        raise TypeError("catalogue must be a validated LocusCatalogue")
    roots = tuple(root_states)
    if not roots or any(not isinstance(state, LocusState) for state in roots):
        raise ValueError("Explicit nonempty root_states are required")
    if len(set(roots)) != len(roots):
        raise ValueError("Duplicate root states are not allowed")
    for state in roots:
        if not state_is_valid(catalogue, state):
            raise ValueError("A root state is invalid for the supplied catalogue")

    by_id = defaultdict(list)
    for opportunity in catalogue.opportunities:
        by_id[opportunity.id].append(opportunity)
    ordered_opportunities = tuple((key, tuple(sorted(value, key=lambda e: e.outcome_id)))
                                  for key, value in sorted(by_id.items()))
    states, index = [], {}
    edges = []
    queue = deque()
    for root in roots:
        index[root] = len(states)
        states.append(root)
        queue.append(root)

    while queue:
        source = queue.popleft()
        source_index = index[source]
        for opportunity_id, outcomes in ordered_opportunities:
            probe = outcomes[0]
            # Preconditions are shared by all outcomes in one opportunity.
            try:
                # Apply each outcome below; this first check only evaluates the
                # common gate without changing state.
                _check_gate(catalogue, source, probe)
            except OpportunityUnavailable:
                continue
            targets, unavailable = [], []
            for event in outcomes:
                try:
                    targets.append((event, apply_outcome(catalogue, source, event)))
                except OpportunityUnavailable as error:
                    unavailable.append((event, error))
            if unavailable and targets:
                raise ValueError(f"Opportunity {opportunity_id} has partially applicable outcomes at a source state")
            if unavailable:
                continue
            for event, target in targets:
                if target == source:
                    # A declared null outcome retains its share of the
                    # opportunity weight and contributes no jump or mark.
                    continue
                if target not in index:
                    index[target] = len(states)
                    states.append(target)
                    queue.append(target)
                edges.append(ProcessEdge(source_index, index[target], event.id,
                                         event.outcome_id, event.kind,
                                         event.rate_group, event.weight))

    return LocusProcess(catalogue, tuple(states), tuple(edges))


def _check_gate(catalogue, state, event):
    """Test only opportunity-level preconditions before outcome resolution."""
    material_index = {tract.id: i for i, tract in enumerate(catalogue.material)}
    if any(state.material[material_index[mid]] != required for mid, required in event.preconditions):
        raise OpportunityUnavailable("material precondition failed")
