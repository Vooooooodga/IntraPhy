"""Finite shared-material exon-repertoire observation constraints.

This module provides constraint vectors for declared finite repertoire
candidates only. It does not implement a transition process, root prior,
detection model, rate fitting, or a joint maximum-likelihood model.
"""
from __future__ import annotations

from dataclasses import dataclass
from collections import deque
import numpy as np

from .types import ExonConfiguration, ObservationEvidence, Catalogue, MaterialState
from .material import normalize_exons, valid_configuration
from .observations import compatibility
from .space import StateSpace


@dataclass(frozen=True)
class ExonRepertoire:
    """Canonical nonempty set of ordered exon configurations sharing DNA."""

    configurations: tuple[ExonConfiguration, ...]
    material: tuple[int, ...]

    def __post_init__(self):
        configs = tuple(self.configurations)
        material = tuple(self.material)
        if not configs:
            raise ValueError("An exon repertoire must contain at least one configuration")
        if any(not isinstance(c, ExonConfiguration) for c in configs):
            raise TypeError("Repertoire members must be ExonConfiguration objects")
        if any(type(v) is not int or v not in (0, 1, 2) for v in material):
            raise ValueError("Repertoire material entries must be integer lifecycle states 0, 1, or 2")
        if any(type(v) is not int or v not in (0, 1, 2)
               for c in configs for v in c.material):
            raise ValueError("Member material entries must be integer lifecycle states 0, 1, or 2")
        if any(c.material != material for c in configs):
            raise ValueError("All repertoire configurations must share one material tuple")
        canonical = tuple(sorted(set(configs)))
        object.__setattr__(self, "configurations", canonical)
        object.__setattr__(self, "material", material)


def _matching(repertoire, observation, space, view):
    natives = tuple(sorted(set(observation.configurations)))
    if observation.kind in {"unknown", "excluded"}:
        return True
    candidates = []
    members = tuple(sorted(repertoire.configurations, key=lambda c: c.key))
    for native in natives:
        allowed = compatibility(space, observation, view, native)
        candidates.append(tuple(i for i, member in enumerate(members)
                                if allowed[space.index[member]] > 0))
    if any(not values for values in candidates):
        return False
    left_match, right_match = {}, {}
    for start in range(len(candidates)):
        queue, seen_left, previous = deque([start]), {start}, {}
        free = None
        while queue and free is None:
            current = queue.popleft()
            for member in candidates[current]:
                if member in previous:
                    continue
                previous[member] = current
                owner = right_match.get(member)
                if owner is None:
                    free = member; break
                if owner not in seen_left:
                    seen_left.add(owner); queue.append(owner)
        if free is None:
            return False
        member = free
        while member is not None:
            current = previous[member]
            old = left_match.get(current)
            right_match[member] = current; left_match[current] = member
            member = old
    return True


def repertoire_compatibility(space: StateSpace, repertoires, observation: ObservationEvidence,
                             view: str = "evidence", *, complete_repertoire: bool = False) -> np.ndarray:
    """Return a 0/1 vector over declared finite joint repertoire candidates only.

    The vector has one entry per declared joint repertoire candidate and reuses
    the exact ``space.index`` to validate every member.
    This layer has no transition process, root prior, detection model, rate
    fitting, or CLI and never claims to implement joint ML. The default is a
    lower-bound constraint, allowing unannotated extra members. Complete mode
    requires every member to be matched for resolved observations; partial and
    unknown completeness declarations are rejected.
    """
    if view not in {"evidence", "annotation"}:
        raise ValueError("Observation view must be evidence or annotation")
    if type(complete_repertoire) is not bool:
        raise TypeError("complete_repertoire must be a Boolean")
    repertoires = tuple(repertoires)
    if any(not isinstance(r, ExonRepertoire) for r in repertoires):
        raise TypeError("Repertoire candidates must be ExonRepertoire objects")
    if observation.kind in {"partial", "unknown", "excluded"} and complete_repertoire:
        raise ValueError("Complete repertoire constraints are unsupported for partial, unknown, or excluded observations")
    result = np.zeros(len(repertoires))
    if len(observation.material_presence) != len(space.catalogue.material):
        raise ValueError("Observation material vector does not match catalogue material dimension")
    for _, repertoire in enumerate(repertoires):
        if any(member not in space.index for member in repertoire.configurations):
            raise ValueError("Repertoire member is outside the declared candidate state set")
        if len(observation.material_presence) != len(repertoire.material):
            raise ValueError("Observation material vector does not match repertoire material dimension")
        if observation.kind == "excluded":
            result[_] = 1
            continue
        if observation.material_presence and any(p is not None and (repertoire.material[k] == 1) != (p == 1)
                                                for k, p in enumerate(observation.material_presence)):
            continue
        if observation.kind == "unknown":
            ok = True
        else:
            if observation.kind in {"observed", "partial"}:
                allowed = compatibility(space, observation, view)
                ok = any(allowed[space.index[member]] > 0 for member in repertoire.configurations)
            else:
                ok = _matching(repertoire, observation, space, view)
            if ok and complete_repertoire:
                ok = (len(repertoire.configurations) == (1 if observation.kind in {"observed", "partial"}
                                                         else len(set(observation.configurations))))
        if ok:
            result[_] = 1
    return result


def delete_repertoire_material(catalogue: Catalogue, repertoire: ExonRepertoire,
                               material_id: str) -> ExonRepertoire:
    """Apply one shared PRESENT-to-DELETED material change to every member."""
    ids = tuple(m.id for m in catalogue.material)
    if material_id not in ids or len(repertoire.material) != len(ids):
        raise ValueError("Material identifier or repertoire dimension is invalid")
    index = ids.index(material_id)
    if repertoire.material[index] != int(MaterialState.PRESENT):
        raise ValueError("Only present material can be deleted")
    material = list(repertoire.material); material[index] = int(MaterialState.DELETED)
    members = []
    declared_spans = frozenset(catalogue.spans)
    for member in repertoire.configurations:
        if not valid_configuration(catalogue, member) or not set(member.exons) <= declared_spans:
            raise ValueError("Repertoire member is not a valid catalogue configuration")
        exons = normalize_exons(catalogue, tuple(material), member.exons)
        if not set(exons) <= declared_spans:
            raise ValueError("Shared deletion leaves the declared exon catalogue")
        members.append(ExonConfiguration(exons, tuple(material)))
    return ExonRepertoire(tuple(members), tuple(material))
