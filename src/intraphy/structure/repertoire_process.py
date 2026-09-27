"""Bounded event-table processes over declared shared-material repertoires.

This module supplies a finite, conditional event-table CTMC primitive. It does
not establish biological completeness, a general structural turnover model, or
an ML/rate-fitting/CLI interface.
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass
import math
import numpy as np

from .material import valid_configuration
from .repertoire import ExonRepertoire, delete_repertoire_material
from .types import Catalogue, MaterialState


@dataclass(frozen=True)
class RepertoireEvent:
    source: ExonRepertoire
    target: ExonRepertoire
    kind: str
    opportunity: str
    weight: float = 1.0


@dataclass(frozen=True)
class RepertoireProcess:
    """Finite conditional event table over declared repertoire states.

    Omitted opportunities are unavailable under the caller-declared finite
    model. Structural or insertion coupling is not inferred or biologically
    validated by this class.
    """
    states: tuple[ExonRepertoire, ...]
    events: tuple[RepertoireEvent, ...]
    provenance: str
    catalogue: Catalogue

    def __post_init__(self):
        states = tuple(self.states)
        if not states or any(not isinstance(s, ExonRepertoire) for s in states):
            raise ValueError("A repertoire process requires nonempty repertoire states")
        if len(set(states)) != len(states):
            raise ValueError("Duplicate repertoire states are forbidden")
        if len({len(s.material) for s in states}) != 1:
            raise ValueError("All process states require one material dimension")
        if not isinstance(self.provenance, str) or not self.provenance.strip():
            raise ValueError("Process provenance is required")
        if not isinstance(self.catalogue, Catalogue) or self.catalogue.status != "qualified":
            raise ValueError("A qualified catalogue is required")
        if any(len(s.material) != len(self.catalogue.material) or
               any(not valid_configuration(self.catalogue, c) or not set(c.exons) <= set(self.catalogue.spans)
                   for c in s.configurations) for s in states):
            raise ValueError("Process states must be valid configurations in the bound catalogue")
        index = {s: i for i, s in enumerate(states)}
        seen = {}
        totals = {}
        unique = []
        for event in self.events:
            if not isinstance(event, RepertoireEvent) or event.source not in index or event.target not in index:
                raise ValueError("Event endpoints must be declared process states")
            if event.source == event.target or not isinstance(event.kind, str) or not event.kind.strip() or not isinstance(event.opportunity, str) or not event.opportunity.strip():
                raise ValueError("Events require distinct endpoints, kind, and opportunity")
            if isinstance(event.weight, bool) or not isinstance(event.weight, (int, float)) or not math.isfinite(event.weight) or event.weight <= 0:
                raise ValueError("Event weights must be finite and positive")
            changes = [(a, b) for a, b in zip(event.source.material, event.target.material) if a != b]
            if event.kind == "dna_deletion" and (len(changes) != 1 or changes[0] != (1, 2)):
                raise ValueError("dna_deletion must change exactly one material from present to deleted")
            if event.kind == "dna_insertion" and (len(changes) != 1 or changes[0] != (0, 1)):
                raise ValueError("dna_insertion must change exactly one material from unintroduced to present")
            if event.kind not in {"dna_deletion", "dna_insertion"} and changes:
                raise ValueError("Non-material events cannot change lifecycle state")
            if event.kind == "dna_deletion":
                material_id = self.catalogue.material[next(i for i, pair in enumerate(zip(event.source.material, event.target.material)) if pair == (1, 2))].id
                if event.opportunity != material_id or event.target != delete_repertoire_material(self.catalogue, event.source, material_id):
                    raise ValueError("dna_deletion event does not match the bound catalogue operation")
            if event.kind == "dna_insertion":
                changed = next(i for i, pair in enumerate(zip(event.source.material, event.target.material)) if pair == (0, 1))
                if event.opportunity != self.catalogue.material[changed].id:
                    raise ValueError("dna_insertion opportunity must name the changed material")
            key = (event.source, event.kind, event.opportunity, event.target)
            if key in seen:
                if seen[key] != event.weight:
                    raise ValueError("Conflicting weights for duplicate event identity")
                continue
            seen[key] = event.weight
            totals[(event.source, event.kind, event.opportunity)] = totals.get((event.source, event.kind, event.opportunity), 0.) + event.weight
            unique.append(event)
        if any(not math.isclose(v, 1., rel_tol=0., abs_tol=1e-12) for v in totals.values()):
            raise ValueError("Opportunity outcome weights must sum to one")
        object.__setattr__(self, "states", states)
        object.__setattr__(self, "events", tuple(unique))

    @property
    def index(self):
        return {state: i for i, state in enumerate(self.states)}

    def generator(self, rates: dict[str, float], *, scale: float = 1., include_marks: bool = True):
        if type(include_marks) is not bool:
            raise TypeError("include_marks must be Boolean")
        if set(rates) != {event.kind for event in self.events}:
            raise ValueError("Explicit rates are required for every event kind")
        if any(isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) or v < 0 for v in rates.values()):
            raise ValueError("Rates must be finite and nonnegative")
        if isinstance(scale, bool) or not isinstance(scale, (int, float)) or not math.isfinite(scale) or scale <= 0:
            raise ValueError("Scale must be finite and positive")
        n = len(self.states); q = np.zeros((n, n)); index = self.index
        marks = {kind: np.zeros((n, n)) for kind in rates} if include_marks else {}
        for event in self.events:
            value = scale * rates[event.kind] * event.weight
            if not math.isfinite(value): raise ArithmeticError("Generator rate overflow")
            i, j = index[event.source], index[event.target]
            q[i, j] += value
            if include_marks: marks[event.kind][i, j] += value
        np.fill_diagonal(q, -q.sum(axis=1))
        if not np.isfinite(q).all() or any(not np.isfinite(m).all() for m in marks.values()):
            raise ArithmeticError("Generator accumulation overflow")
        return q, marks


def deletion_process(catalogue: Catalogue, seeds, *, provenance: str, max_states: int | None = None):
    """Build a deletion-only conditional closure over shared repertoire states.

    Generic opportunities omitted from this declared event table are
    unavailable under this conditional process; the table does not establish
    biological state coverage.
    """
    if not isinstance(catalogue, Catalogue) or catalogue.status != "qualified":
        raise ValueError("A qualified catalogue is required")
    if not isinstance(provenance, str) or not provenance.strip():
        raise ValueError("Process provenance is required")
    if max_states is not None and (type(max_states) is not int or max_states <= 0):
        raise ValueError("max_states must be a positive integer or None")
    seeds = tuple(seeds)
    if not seeds:
        raise ValueError("At least one seed repertoire is required")
    queue = deque(); queued = set(); states = []; known = set()
    for seed in seeds:
        if not isinstance(seed, ExonRepertoire) or len(seed.material) != len(catalogue.material) or any(not valid_configuration(catalogue, c) or not set(c.exons) <= set(catalogue.spans) for c in seed.configurations):
            raise ValueError("Seed repertoire is invalid for the catalogue")
        if seed in queued:
            continue
        if max_states is not None and len(queued) >= max_states:
            raise ValueError("state_space_incomplete")
        queued.add(seed); queue.append(seed)
    events = []
    while queue:
        source = queue.popleft()
        if source in known: continue
        known.add(source); states.append(source)
        for k, material in enumerate(catalogue.material):
            if source.material[k] != int(MaterialState.PRESENT): continue
            target = delete_repertoire_material(catalogue, source, material.id)
            events.append(RepertoireEvent(source, target, "dna_deletion", material.id, 1.))
            if target not in known and target not in queued:
                if max_states is not None and len(queued) >= max_states:
                    raise ValueError("state_space_incomplete")
                queued.add(target); queue.append(target)
    return RepertoireProcess(tuple(states), tuple(events), provenance, catalogue)
