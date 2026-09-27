"""Evidence-conditioned labelled-copy locus state types.

Coordinates use one ordered, interbase coordinate axis supplied by the caller.
The catalogue records the finite copy/feature support of a model; it does not
infer unlisted homologues or unrestricted copy birth.
"""
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
class SpliceFeature:
    """A selectable structural feature; availability does not imply usage."""
    id: str
    kind: str
    required_material: tuple[str, ...]
    prerequisites: tuple[str, ...] = ()
    copy_id: str | None = None
    start: int | None = None
    end: int | None = None
    donor: int | None = None
    acceptor: int | None = None
    mxe_group: str | None = None
    evidence: tuple[str, ...] = ()

    def __post_init__(self):
        _identifier(self.id, "feature id")
        _identifier(self.kind, "feature kind")
        if self.kind not in {"exon", "splice"}:
            raise ValueError("Feature kind must be exon or splice")
        if len(set(self.required_material)) != len(self.required_material):
            raise ValueError("Feature required_material entries must be unique")
        for value in (*self.required_material, *self.prerequisites):
            _identifier(value, "feature dependency")
        if self.copy_id is not None:
            _identifier(self.copy_id, "feature copy id")
        if self.mxe_group is not None:
            _identifier(self.mxe_group, "MXE group")
        if (self.start is None) != (self.end is None):
            raise ValueError("Feature interval requires both start and end")
        if self.start is not None and (type(self.start) is not int or type(self.end) is not int or not 0 <= self.start < self.end):
            raise ValueError("Feature interval requires 0 <= start < end")
        for value in (self.donor, self.acceptor):
            if value is not None and type(value) is not int:
                raise ValueError("Splice endpoint coordinates must be integers")
        if self.donor is not None and self.acceptor is not None and self.donor >= self.acceptor:
            raise ValueError("A directed exon requires donor < acceptor on the declared axis")
        if self.kind == "exon" and (self.start is None or self.copy_id is None):
            raise ValueError("An exon feature requires an interval and labelled copy")
        if self.kind == "splice" and (self.donor is None or self.acceptor is None or len(self.prerequisites) < 2):
            raise ValueError("A splice feature requires donor/acceptor endpoints and two prerequisite exon features")
        if any(not isinstance(x, str) or not x.strip() for x in self.evidence):
            raise ValueError("Feature evidence entries must be nonempty strings")


@dataclass(frozen=True)
class EventOpportunity:
    """One explicitly supported outcome of a typed biological opportunity.

    Records sharing ``id`` are alternative outcomes and their weights sum to
    one. The same opportunity is gated as a whole by ``preconditions``.
    """
    id: str
    outcome_id: str
    kind: str
    rate_group: str
    weight: float = 1.0
    preconditions: tuple[tuple[str, int], ...] = ()
    required_features: tuple[str, ...] = ()
    forbidden_features: tuple[str, ...] = ()
    material_gains: tuple[str, ...] = ()
    material_deletions: tuple[str, ...] = ()
    feature_on: tuple[str, ...] = ()
    feature_off: tuple[str, ...] = ()
    source_copy: str | None = None
    target_copy: str | None = None
    feature_map: tuple[tuple[str, str], ...] = ()
    context_features: tuple[str, ...] = ()
    graft_features: tuple[str, ...] = ()
    interval: tuple[int, int] | None = None

    def __post_init__(self):
        for value, name in ((self.id, "opportunity id"), (self.kind, "event kind"), (self.rate_group, "rate group")):
            _identifier(value, name)
        _identifier(self.outcome_id, "opportunity outcome id")
        if isinstance(self.weight, bool) or not isinstance(self.weight, (int, float)) or not math.isfinite(self.weight) or self.weight <= 0:
            raise ValueError("Opportunity outcome weight must be finite and positive")
        if self.source_copy is not None:
            _identifier(self.source_copy, "source copy")
        if self.target_copy is not None:
            _identifier(self.target_copy, "target copy")
        if len(set(self.material_gains)) != len(self.material_gains) or len(set(self.material_deletions)) != len(self.material_deletions):
            raise ValueError("Material gain/deletion IDs must be unique within an outcome")
        if any(len(set(values)) != len(values) for values in (self.feature_on, self.feature_off,
                   self.context_features, self.required_features, self.forbidden_features,
                   self.graft_features)):
            raise ValueError("Feature IDs in each opportunity field must be unique")
        if set(self.material_gains) & set(self.material_deletions) or set(self.feature_on) & set(self.feature_off):
            raise ValueError("An event outcome cannot both gain and delete the same object")
        if len({key for key, _ in self.preconditions}) != len(self.preconditions):
            raise ValueError("Precondition material IDs must be unique")
        for material_id, state in self.preconditions:
            _identifier(material_id, "precondition material")
            if type(state) is not int or state not in (0, 1, 2):
                raise ValueError("Precondition states are 0=unintroduced, 1=present, 2=deleted")
        for value in (*self.material_gains, *self.material_deletions, *self.feature_on,
                      *self.feature_off, *self.context_features, *self.required_features,
                      *self.forbidden_features, *self.graft_features):
            _identifier(value, "event target")
        if set(self.required_features) & set(self.forbidden_features):
            raise ValueError("An opportunity cannot require and forbid the same feature")
        if len(set(self.feature_map)) != len(self.feature_map) or any(len(pair) != 2 for pair in self.feature_map):
            raise ValueError("feature_map must contain unique source/target feature pairs")
        if self.interval is not None and (len(self.interval) != 2 or any(type(x) is not int for x in self.interval) or self.interval[0] >= self.interval[1]):
            raise ValueError("Event interval must be an increasing integer coordinate pair")


@dataclass(frozen=True)
class LocusCatalogue:
    material: tuple[MaterialTract, ...]
    copies: tuple[CopySlot, ...]
    features: tuple[SpliceFeature, ...]
    opportunities: tuple[EventOpportunity, ...]
    mxe_groups: tuple[tuple[str, tuple[str, ...]], ...] = ()
    provenance: str = ""

    def __post_init__(self):
        from .locus_validation import validate_locus_catalogue
        validate_locus_catalogue(self)


@dataclass(frozen=True)
class LocusState:
    material: tuple[int, ...]
    active_features: frozenset[str]

    def __post_init__(self):
        object.__setattr__(self, "material", tuple(self.material))
        if any(type(x) is not int or x not in (0, 1, 2) for x in self.material):
            raise ValueError("Material state entries are 0=unintroduced, 1=present, 2=deleted")
        object.__setattr__(self, "active_features", frozenset(self.active_features))
        for feature in self.active_features:
            _identifier(feature, "active feature")


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
        if not self.states or len(set(self.states)) != len(self.states):
            raise ValueError("A locus process requires unique reachable states")
        if any(not isinstance(s, LocusState) for s in self.states):
            raise TypeError("Process states must be LocusState values")
        if any(len(s.material) != len(self.catalogue.material) for s in self.states):
            raise ValueError("Process material dimension disagrees with catalogue")
        if any(not e.source >= 0 or e.source >= len(self.states) or e.target < 0 or e.target >= len(self.states) or e.source == e.target for e in self.edges):
            raise ValueError("Process edges must connect distinct indexed states")
        if any(not math.isfinite(e.weight) or e.weight <= 0 for e in self.edges):
            raise ValueError("Process edge weights must be finite and positive")

    @property
    def index(self):
        return {state: i for i, state in enumerate(self.states)}


@dataclass(frozen=True)
class LocusObservation:
    """Tip constraints; None means unobserved, never biological absence."""
    material: tuple[int | None, ...]
    features: tuple[tuple[str, int | None], ...] = ()
    surveyed_features: frozenset[str] = frozenset()
    observed_paths: tuple[tuple[str, ...], ...] = ()
    surveyed_material: frozenset[str] = frozenset()

    def __post_init__(self):
        if any(x is not None and (type(x) is not int or x not in (0, 1)) for x in self.material):
            raise ValueError("Observed DNA presence is 0, 1, or None; deletion history is latent")
        if len({key for key, _ in self.features}) != len(self.features):
            raise ValueError("Feature observations must have unique IDs")
        if any(not isinstance(key, str) or not key.strip() or
               (value is not None and (type(value) is not int or value not in (0, 1)))
               for key, value in self.features):
            raise ValueError("Feature observations require an ID and value 0, 1, or None")
        object.__setattr__(self, "surveyed_features", frozenset(self.surveyed_features))
        if not self.surveyed_features <= {key for key, _ in self.features}:
            raise ValueError("Surveyed features require an explicit observation entry")
        object.__setattr__(self, "surveyed_material", frozenset(self.surveyed_material))
        if any(not isinstance(value, str) or not value.strip() for value in self.surveyed_material):
            raise ValueError("Surveyed material IDs must be nonempty strings")
        if any(not path or any(not isinstance(value, str) or not value.strip() for value in path)
               or len(set(path)) != len(path) for path in self.observed_paths):
            raise ValueError("Observed splice paths must contain unique exon feature IDs")
