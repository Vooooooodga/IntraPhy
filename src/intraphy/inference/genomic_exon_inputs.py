"""Shared catalogue diagnostics and eligible-unit preparation for genomic-exon models."""
from __future__ import annotations

from dataclasses import asdict
import json

import numpy as np

from ..structure.observations import observation_scenarios
from ..structure.space import enumerate_space
from .genomic_exon_output import state_rows
from ..structure.serialization import json_safe


def _write_state_space_progress(handle, record):
    handle.write(json.dumps(json_safe(record), sort_keys=True, separators=(",", ":")) + "\n")
    handle.flush()


def _diagnose_catalogue(catalogue, taxa, tree, max_states, state_space_only, progress,
                        *, include_detail=True):
    _write_state_space_progress(progress, {"family_id": catalogue.family,
        "unit_id": catalogue.unit,
        "stage": "enumeration_started" if catalogue.status == "qualified" else "unit_started"})
    diag = {"family_id": catalogue.family, "unit_id": catalogue.unit,
        "catalogue_status": catalogue.status, "catalogue_reasons": list(catalogue.reasons),
        "observation_unit": catalogue.observation_unit}
    unit_detail = None
    unit_input = None
    if catalogue.status != "qualified":
        diag.update(status="unqualified_catalogue", probability_status="not_available",
                    state_count=None, state_space_complete=None,
                    state_space_reason="not_enumerated_unqualified_catalogue")
    else:
        space = enumerate_space(catalogue, max_states)
        diag.update(state_count=len(space.states), state_space_complete=space.complete,
                    state_space_reason=space.reason, state_space_diagnostics=space.diagnostics)
        if not space.complete:
            diag.update(status="state_space_incomplete", probability_status="not_available")
        else:
            try:
                scenarios = observation_scenarios(space, taxa, "annotation")
                if len(scenarios) != 1:
                    diag.update(status="coexisting_structures_unresolved", probability_status="not_available")
                else:
                    label, tips = scenarios[0]
                    informative = sum(not np.all(values == values[0]) for values in tips.values())
                    unknown_tips = [species for species, values in tips.items()
                                    if np.all(values == 1)]
                    diag.update(observation_scenario=label, informative_tips=informative,
                                unknown_tips=unknown_tips)
                    if informative < 2:
                        diag.update(status="fewer_than_two_informative_tips", probability_status="not_available")
                    else:
                        diag.update(status="eligible_conditional_unit",
                            probability_status="not_computed" if state_space_only else "pending")
                        if not state_space_only and include_detail:
                            unit_detail = {"catalogue": asdict(catalogue),
                                "states": state_rows(space), "observation_scenario": label,
                                "unknown_tips": unknown_tips}
                        if not state_space_only:
                            unit_input = {"unit_id": catalogue.unit, "space": space,
                                "tree": tree, "tips": tips}
                        if unknown_tips:
                            diag["status"] = "eligible_with_unknown_tips"
            except ValueError as exc:
                diag.update(status="observation_unresolved", probability_status="not_available",
                            reason=str(exc))
    _write_state_space_progress(progress, {**diag,
        "stage": ("enumeration_and_eligibility_complete" if catalogue.status == "qualified"
                  else "unit_assessment_complete")})
    return diag, unit_detail, unit_input
