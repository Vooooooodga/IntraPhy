"""Read-only compatibility between structural observations and candidate states."""
from __future__ import annotations

from itertools import islice, product
import numpy as np
from .types import Catalogue, ExonConfiguration, ExonSpan, ObservationEvidence
from .material import normalize_exons
from .space import StateSpace


def _outside(exons: tuple[ExonSpan, ...], windows: tuple[ExonSpan, ...]) -> tuple[tuple[int, int], ...]:
    """Compare exon structure outside unknown windows without joining across them."""
    parts = []
    for exon in exons:
        remaining = [(exon.start, exon.end)]
        for w in windows:
            cut = []
            for a, b in remaining:
                if b <= w.start or a >= w.end:
                    cut.append((a, b))
                else:
                    if a < w.start:
                        cut.append((a, w.start))
                    if b > w.end:
                        cut.append((w.end, b))
            remaining = cut
        parts.extend(remaining)
    return tuple(parts)


def compatibility(space: StateSpace, observation: ObservationEvidence,
                  view: str = "evidence", selected: ExonConfiguration | None = None) -> np.ndarray:
    if view not in {"evidence", "annotation"}:
        raise ValueError("Observation view must be evidence or annotation")
    if observation.kind == "excluded":
        return np.ones(len(space.states))
    if observation.kind == "coexisting" and selected is None:
        raise ValueError("Coexisting annotated structures require explicit conditional scenarios")
    natives = (selected,) if selected is not None else observation.configurations
    candidates = list(natives)
    if view == "evidence":
        keys = {v.key for v in natives}
        candidates.extend(a.configuration for a in observation.alternatives if a.replaces_key in keys)
    allowed = np.zeros(len(space.states))
    for i, state in enumerate(space.states):
        if any(p is not None and (state.material[k] == 1) != (p == 1)
               for k, p in enumerate(observation.material_presence)):
            continue
        if observation.kind == "unknown":
            allowed[i] = 1
            continue
        windows = observation.unknown_intervals
        if view == "evidence":
            windows = tuple(sorted(set((*windows, *observation.alternative_exons))))
        for native in candidates:
            # Observed absence never selects unintroduced versus subsequently lost.
            native_exons = normalize_exons(space.catalogue, state.material, native.exons)
            if _outside(state.exons, windows) == _outside(native_exons, windows):
                allowed[i] = 1
                break
    return allowed


class ObservationScenarios:
    """Lazy, re-iterable Cartesian scenarios with shared compatibility vectors."""

    def __init__(self, space: StateSpace, taxa: tuple[str, ...], view: str,
                 max_scenarios: int | None = None):
        if max_scenarios is not None and (not isinstance(max_scenarios, int) or isinstance(max_scenarios, bool)
                                          or max_scenarios <= 0):
            raise ValueError("max_observation_scenarios must be a positive integer or None")
        if view not in {"evidence", "annotation"}:
            raise ValueError("Observation view must be evidence or annotation")
        self._space, self._taxa, self._view = space, taxa, view
        self._by_species = {o.species: o for o in space.catalogue.observations}
        self._coexist = [(s, self._by_species[s].configurations) for s in taxa
                         if s in self._by_species and self._by_species[s].kind == "coexisting"]
        count = 1
        for _, choices in self._coexist:
            count *= len(choices)
        self.theoretical_count = count
        self.limit = max_scenarios
        if max_scenarios is not None and count > max_scenarios:
            raise ValueError("observation_scenarios_incomplete: too many coexisting structures")
        self._vectors = {}
        for species in taxa:
            observation = self._by_species.get(species)
            if observation is None:
                self._vectors[(species, None)] = self._freeze(np.ones(len(space.states)))
            elif observation.kind == "coexisting":
                for selected in observation.configurations:
                    self._vectors[(species, selected.key)] = self._freeze(
                        compatibility(space, observation, view, selected))
            else:
                self._vectors[(species, None)] = self._freeze(compatibility(space, observation, view))
            for (vector_species, vector_key), vector in self._vectors.items():
                if vector_species == species and not vector.any():
                    raise ValueError(f"Observation incompatible with candidate catalogue: {species}")

    @staticmethod
    def _freeze(vector):
        vector.setflags(write=False)
        return vector

    def __len__(self):
        return self.theoretical_count

    def __iter__(self):
        for choices in product(*(cs for _, cs in self._coexist)):
            selected = dict(zip((s for s, _ in self._coexist), choices))
            weights = {}
            for species in self._taxa:
                config = selected.get(species)
                weights[species] = self._vectors[(species, config.key if config else None)]
            label = ";".join(f"{s}={config.key}" for s, config in sorted(selected.items())) or "single_structure"
            yield label, weights

    def __getitem__(self, index):
        if isinstance(index, slice):
            return list(self)[index]
        if index < 0:
            index += len(self)
        if index < 0 or index >= len(self):
            raise IndexError(index)
        return next(islice(self, index, None))


def observation_scenarios(space: StateSpace, taxa: tuple[str, ...], view: str,
                          max_scenarios: int | None = None):
    return ObservationScenarios(space, taxa, view, max_scenarios)


def observation_scenario_count(space: StateSpace, taxa: tuple[str, ...]) -> int:
    count = 1
    by_species = {o.species: o for o in space.catalogue.observations}
    for species in taxa:
        observation = by_species.get(species)
        if observation is not None and observation.kind == "coexisting":
            count *= len(observation.configurations)
    return count
