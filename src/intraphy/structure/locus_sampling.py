"""Parametric sampling for fixed-catalogue DNA-only locus observations."""
from __future__ import annotations

from typing import Mapping

import numpy as np

from .locus_types import LocusObservation


def sample_history(process, tree, root_prior, rates: Mapping[str, float], rng):
    """Sample one latent root state and evolve every branch by Gillespie jumps."""
    prior = np.asarray(root_prior, dtype=float)
    if prior.shape != (len(process.states),) or not np.all(np.isfinite(prior)) or np.any(prior < 0):
        raise ValueError("root prior must be finite, nonnegative, and match process states")
    if not np.isclose(float(prior.sum()), 1.0, rtol=0.0, atol=1e-12):
        raise ValueError("root prior must sum to one")
    by_source = {}
    for edge in process.edges:
        if edge.rate_group not in rates:
            raise ValueError(f"missing sampling rate group {edge.rate_group!r}")
        hazard = float(rates[edge.rate_group]) * float(edge.weight)
        if not np.isfinite(hazard) or hazard < 0:
            raise ValueError(f"invalid hazard for rate group {edge.rate_group!r}")
        if hazard:
            by_source.setdefault(edge.source, []).append((edge.target, hazard))

    root_index = int(rng.choice(len(process.states), p=prior))
    states = {tree.root: root_index}
    for parent in tree.preorder():
        source = states[parent]
        for child in tree.children.get(parent, ()):
            branch_length = tree.branch_length(child)
            elapsed, current = 0.0, source
            while elapsed < branch_length:
                active = by_source.get(current, ())
                total = sum(rate for _target, rate in active)
                if not np.isfinite(total):
                    raise FloatingPointError("total outgoing event hazard is non-finite")
                if total <= 0:
                    break
                unit_hazard_time = float(rng.exponential(1.0))
                if not np.isfinite(unit_hazard_time) or unit_hazard_time <= 0:
                    raise FloatingPointError("Gillespie unit-hazard waiting time is invalid")
                remaining = branch_length - elapsed
                if unit_hazard_time >= total * remaining:
                    break
                wait = unit_hazard_time / total
                if not np.isfinite(wait) or wait <= 0 or elapsed + wait <= elapsed:
                    raise FloatingPointError("Gillespie waiting time failed to advance branch time")
                draw = float(rng.random()) * total
                cumulative = 0.0
                current = active[-1][0]
                for target, rate in active:
                    cumulative += rate
                    if draw < cumulative:
                        current = target
                        break
                elapsed += wait
            states[child] = current
    return {node: process.states[index] for node, index in states.items()}


def sample_tip_observation(catalogue, tip_record, latent_state, rng):
    """Sample binary DNA calls while retaining the supplied missingness mask."""
    material = []
    sens = tip_record.get("material_sensitivity", {})
    spec = tip_record.get("material_specificity", {})
    if len(tip_record["material"]) != len(catalogue.material):
        raise ValueError("tip material mask does not match the fixed catalogue")
    for index, (tract, recorded) in enumerate(zip(catalogue.material, tip_record["material"])):
        if recorded is None:
            material.append(None)
            continue
        present = latent_state.material[index] == 1
        sensitivity = float(sens.get(tract.id, 1.0))
        specificity = float(spec.get(tract.id, 1.0))
        probability_positive = sensitivity if present else 1.0 - specificity
        material.append(int(float(rng.random()) < probability_positive))
    return LocusObservation(tuple(material), surveyed_material=frozenset(tip_record.get("surveyed_material", ())))
