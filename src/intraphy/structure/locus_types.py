"""Finite, evidence-conditioned labelled DNA copy states."""
from __future__ import annotations

from dataclasses import dataclass
import math


def _identifier(value: str, name: str) -> None:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty string")


@dataclass(frozen=True, order=True)
class MaterialTract:
    id: str
    start: int
    end: int
    evidence: tuple[str, ...] = ()

    def __post_init__(self):
        _identifier(self.id, "material id")
        if type(self.start) is not int or type(self.end) is not int or not 0 <= self.start < self.end:
            raise ValueError("Material tracts require 0 <= start < end integer coordinates")
        if any(not isinstance(x, str) or not x.strip() for x in self.evidence):
            raise ValueError("Material evidence entries must be nonempty strings")


@dataclass(frozen=True)
class CopySlot:
    """A labelled genomic copy position supported by supplied evidence."""
    id: str
    material_ids: tuple[str, ...]
    homologous_to: tuple[str, ...] = ()
    collinear_with: tuple[str, ...] = ()
    orientation: str = "+"
    evidence: tuple[str, ...] = ()

    def __post_init__(self):
        _identifier(self.id, "copy id")
        if not self.material_ids or len(set(self.material_ids)) != len(self.material_ids):
            raise ValueError("A copy slot requires unique material tract identifiers")
        if any(len(set(values)) != len(values) for values in (self.homologous_to, self.collinear_with, self.evidence)):
            raise ValueError("Copy homology, collinearity, and evidence entries must be unique")
        if self.orientation not in {"+", "-", "unknown"}:
            raise ValueError("Copy orientation must be +, -, or unknown")
        for values in (self.material_ids, self.homologous_to, self.collinear_with, self.evidence):
            for value in values:
                _identifier(value, "copy metadata entry")


@dataclass(frozen=True)
class EventOpportunity:
    """One outcome of a copy duplication or continuous DNA deletion."""
    id: str
    outcome_id: str
    kind: str
    rate_group: str
    weight: float = 1.0
    preconditions: tuple[tuple[str, int], ...] = ()
    material_gains: tuple[str, ...] = ()
    material_deletions: tuple[str, ...] = ()
    source_copy: str | None = None
    target_copy: str | None = None
    interval: tuple[int, int] | None = None

    def __post_init__(self):
        for value, name in ((self.id, "opportunity id"), (self.outcome_id, "opportunity outcome id"),
                            (self.kind, "event kind"), (self.rate_group, "rate group")):
            _identifier(value, name)
        if isinstance(self.weight, bool) or not isinstance(self.weight, (int, float)) or not math.isfinite(self.weight) or self.weight <= 0:
            raise ValueError("Opportunity outcome weight must be finite and positive")
        for value in (self.source_copy, self.target_copy):
            if value is not None:
                _identifier(value, "copy id")
        if len(set(self.material_gains)) != len(self.material_gains) or len(set(self.material_deletions)) != len(self.material_deletions):
            raise ValueError("Material gain/deletion IDs must be unique within an outcome")
        if set(self.material_gains) & set(self.material_deletions):
            raise ValueError("An outcome cannot gain and delete the same material")
        if len({key for key, _ in self.preconditions}) != len(self.preconditions):
            raise ValueError("Precondition material IDs must be unique")
        for material_id, state in self.preconditions:
            _identifier(material_id, "precondition material")
            if type(state) is not int or state not in (0, 1, 2):
                raise ValueError("Precondition states must be 0, 1, or 2")
        for value in (*self.material_gains, *self.material_deletions):
            _identifier(value, "event material")
        if self.interval is not None and (len(self.interval) != 2 or any(type(x) is not int for x in self.interval)
                                          or self.interval[0] >= self.interval[1]):
            raise ValueError("Event interval must be an increasing integer coordinate pair")


@dataclass(frozen=True)
class LocusCatalogue:
    material: tuple[MaterialTract, ...]
    copies: tuple[CopySlot, ...]
    opportunities: tuple[EventOpportunity, ...]
    provenance: str = ""
    state_model: str = "irreversible"

    def __post_init__(self):
        from .locus_validation import validate_locus_catalogue
        validate_locus_catalogue(self)

    @property
    def material_states(self) -> tuple[int, ...]:
        return (0, 1) if self.state_model == "binary" else (0, 1, 2)


@dataclass(frozen=True)
class LocusState:
    material: tuple[int, ...]

    def __post_init__(self):
        object.__setattr__(self, "material", tuple(self.material))
        if any(type(x) is not int or x not in (0, 1, 2) for x in self.material):
            raise ValueError("Material state entries must be 0, 1, or 2")


@dataclass(frozen=True)
class ProcessEdge:
    source: int
    target: int
    opportunity_id: str
    outcome_id: str
    kind: str
    rate_group: str
    weight: float

    def __post_init__(self):
        if type(self.source) is not int or type(self.target) is not int or self.source < 0 or self.target < 0 or self.source == self.target:
            raise ValueError("Process edge endpoints must be distinct nonnegative state indices")
        for value in (self.opportunity_id, self.outcome_id, self.kind, self.rate_group):
            _identifier(value, "process edge identifier")
        if isinstance(self.weight, bool) or not isinstance(self.weight, (int, float)) or not math.isfinite(self.weight) or self.weight <= 0:
            raise ValueError("Process edge weight must be finite and positive")


@dataclass(frozen=True)
class LocusProcess:
    catalogue: LocusCatalogue
    states: tuple[LocusState, ...]
    edges: tuple[ProcessEdge, ...]

    def __post_init__(self):
        from .locus_events import state_is_valid
        if not self.states or len(set(self.states)) != len(self.states):
            raise ValueError("A locus process requires unique reachable states")
        if any(not isinstance(s, LocusState) or not state_is_valid(self.catalogue, s) for s in self.states):
            raise ValueError("Process states must be valid for the catalogue")
        if any(e.source < 0 or e.source >= len(self.states) or e.target < 0 or e.target >= len(self.states) for e in self.edges):
            raise ValueError("Process edges must connect indexed states")

    @property
    def index(self):
        return {state: i for i, state in enumerate(self.states)}


@dataclass(frozen=True)
class LocusObservation:
    """DNA tip calls; None is unobserved, never biological absence."""
    material: tuple[int | None, ...]
    surveyed_material: frozenset[str] = frozenset()

    def __post_init__(self):
        object.__setattr__(self, "material", tuple(self.material))
        if any(x is not None and (type(x) is not int or x not in (0, 1)) for x in self.material):
            raise ValueError("Observed DNA presence is 0, 1, or None")
        object.__setattr__(self, "surveyed_material", frozenset(self.surveyed_material))
        for value in self.surveyed_material:
            _identifier(value, "surveyed material id")
