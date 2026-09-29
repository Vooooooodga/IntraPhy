"""Immutable biological objects. Coordinates are interbase, half-open and local."""
from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum
from typing import Iterable


class MaterialState(IntEnum):
    UNINTRODUCED = 0
    PRESENT = 1
    DELETED = 2


@dataclass(frozen=True, order=True)
class ExonSpan:
    start: int
    end: int

    def __post_init__(self) -> None:
        if isinstance(self.start, bool) or isinstance(self.end, bool):
            raise ValueError("Exon coordinates must be integers, not Boolean values")
        if not isinstance(self.start, int) or not isinstance(self.end, int):
            raise ValueError("Exon coordinates must be integers")
        if not 0 <= self.start < self.end:
            raise ValueError("An exon requires 0 <= start < end")

    def overlaps(self, other: ExonSpan) -> bool:
        return self.start < other.end and other.start < self.end


@dataclass(frozen=True, order=True)
class Material:
    id: str
    start: int
    end: int

    def __post_init__(self) -> None:
        ExonSpan(self.start, self.end)
        if not self.id:
            raise ValueError("Material requires an explicit source identifier")


@dataclass(frozen=True, order=True)
class ExonConfiguration:
    exons: tuple[ExonSpan, ...]
    material: tuple[int, ...] = ()

    def __post_init__(self) -> None:
        if tuple(sorted(self.exons)) != self.exons:
            raise ValueError("Exons must be ordered in the common transcriptional axis")
        if any(a.end > b.start for a, b in zip(self.exons, self.exons[1:])):
            raise ValueError("A single configuration cannot contain overlapping exons")
        if any(type(s) is not int or s not in (0, 1, 2) for s in self.material):
            raise ValueError("Material states are 0=unintroduced, 1=present, 2=deleted")

    @property
    def key(self) -> str:
        intervals = ",".join(f"{e.start}:{e.end}" for e in self.exons) or "empty"
        return intervals + "|" + "".join(map(str, self.material))


@dataclass(frozen=True)
class ExonInstance:
    id: str
    family: str
    species: str
    locus: str
    contig: str
    start: int
    end: int
    strand: str
    transcripts: tuple[str, ...]
    cds: tuple[tuple[int, int, int | None], ...] = ()
    source: str = "annotation"
    comparison_start: int | None = None
    comparison_end: int | None = None

    def __post_init__(self) -> None:
        ExonSpan(self.start, self.end)
        if self.strand not in ("+", "-"):
            raise ValueError("An exon requires a resolved strand")


@dataclass(frozen=True)
class ExonCorrespondence:
    id: str
    exon_ids: tuple[str, ...]
    relation: str
    blocks: tuple[tuple[str, int, int, int, int], ...] = ()
    evidence: tuple[str, ...] = ()
    status: str = "candidate"


@dataclass(frozen=True)
class ConfigurationAlternative:
    """Atomic uncertainty about one annotation, never an observed second isoform."""
    configuration: ExonConfiguration
    replaces_key: str
    source_species: str
    source_transcript: str
    evidence: str
    identity: float

    def __post_init__(self) -> None:
        if not all((self.replaces_key, self.source_species, self.source_transcript, self.evidence)):
            raise ValueError("A structure alternative needs the original and an evidence source")
        if not 0 <= self.identity <= 1:
            raise ValueError("Alternative identity must lie in [0, 1]")


@dataclass(frozen=True)
class InsertionPayload:
    """A source-supported exon module carried by one declared inserted tract."""
    material_id: str
    exons: tuple[ExonSpan, ...]
    source_id: str
    evidence: str

    def __post_init__(self) -> None:
        if not all((self.material_id, self.source_id, self.evidence, self.exons)):
            raise ValueError("Insertion payload requires a source, evidence and exon structure")
        ExonConfiguration(self.exons)
        if any(a.end >= b.start for a, b in zip(self.exons, self.exons[1:])):
            raise ValueError("Payload exons require nonempty intervening intervals")


@dataclass(frozen=True)
class ObservationEvidence:
    species: str
    configurations: tuple[ExonConfiguration, ...]
    kind: str = "observed"
    reasons: tuple[str, ...] = ()
    native_exon_ids: tuple[str, ...] = ()
    # Unknown is not absence. None means uninformative material status.
    material_presence: tuple[int | None, ...] = ()
    unknown_intervals: tuple[ExonSpan, ...] = ()
    alternative_exons: tuple[ExonSpan, ...] = ()
    alternatives: tuple[ConfigurationAlternative, ...] = ()

    def __post_init__(self) -> None:
        if self.kind not in {"observed", "partial", "unknown", "coexisting", "excluded"}:
            raise ValueError(f"Unsupported observation kind: {self.kind}")
        if self.kind == "coexisting" and len(set(self.configurations)) < 2:
            raise ValueError("Coexisting structures require at least two distinct configurations")
        if any(x not in (0, 1, None) for x in self.material_presence):
            raise ValueError("Observed material presence is 0, 1, or unknown")


@dataclass(frozen=True)
class ElementaryEdit:
    source: ExonConfiguration
    target: ExonConfiguration
    kind: str
    opportunity: str
    footprint: tuple[int, int]
    affected: tuple[ExonSpan, ...]
    consequences: tuple[str, ...] = ()
    material_id: str | None = None
    weight: float = 1.0

    @property
    def event_key(self) -> str:
        return f"{self.kind}:{self.footprint[0]}:{self.footprint[1]}:{self.material_id or ''}"


@dataclass(frozen=True)
class Catalogue:
    family: str
    unit: str
    length: int
    spans: tuple[ExonSpan, ...]
    junctions: tuple[tuple[int, int], ...]
    material: tuple[Material, ...] = ()
    observations: tuple[ObservationEvidence, ...] = ()
    exon_instances: tuple[ExonInstance, ...] = ()
    alignment_offset: int = 0
    status: str = "qualified"
    reasons: tuple[str, ...] = ()
    discovery: str = "annotation_discovered"
    boundary_candidates: tuple[ExonSpan, ...] = ()
    insertion_payloads: tuple[InsertionPayload, ...] = ()
    observation_unit: str = "transcript_configuration"

    def __post_init__(self) -> None:
        if not self.family or not self.unit:
            raise ValueError("Family and local unit identifiers must be explicit")
        if self.status not in {"qualified", "unresolved"}:
            raise ValueError("Catalogue status must be qualified or unresolved")
        if self.discovery not in {"annotation_discovered", "independent_catalogue"}:
            raise ValueError("An explicit supported discovery scheme is required")
        if self.observation_unit not in {"transcript_configuration", "genomic_exon_spans"}:
            raise ValueError("Unsupported catalogue observation unit")
        if type(self.alignment_offset) is not int or self.alignment_offset < 0:
            raise ValueError("Alignment offset must be a nonnegative integer")
        if type(self.length) is not int or self.length < 1:
            raise ValueError("A catalogue requires a positive coordinate length")
        if len({m.id for m in self.material}) != len(self.material):
            raise ValueError("Material identifiers must be unique")
        ordered = sorted(self.material, key=lambda m: m.start)
        if any(a.end > b.start for a, b in zip(ordered, ordered[1:])):
            raise ValueError("Overlapping material histories require a different model")
        if any(e.end > self.length for e in (*self.spans, *self.boundary_candidates)):
            raise ValueError("Exon outside the declared local coordinate range")
        if any(m.end > self.length for m in self.material):
            raise ValueError("Material outside the declared coordinate range")
        if any(not 0 <= d < a <= self.length for d, a in self.junctions):
            raise ValueError("Junctions require ordered donor and acceptor interbase positions")
        if len({o.species for o in self.observations}) != len(self.observations):
            raise ValueError("Duplicate species observation")
        for o in self.observations:
            if len(o.material_presence) != len(self.material):
                raise ValueError("Observation material vector does not match catalogue")
            if not o.species:
                raise ValueError("Empty observed species identifier")
            if any(a.replaces_key not in {v.key for v in o.configurations} for a in o.alternatives):
                raise ValueError("Alternative must identify one of the native configurations")
            for v in (*o.configurations, *(a.configuration for a in o.alternatives)):
                if len(v.material) != len(self.material) or any(e.end > self.length for e in v.exons):
                    raise ValueError("Observed configuration does not match the coordinate/material catalogue")
            if any(e.end > self.length for e in (*o.unknown_intervals, *o.alternative_exons)):
                raise ValueError("Uncertainty interval lies outside the declared catalogue")
        sources = {m.id: m for m in self.material}
        for payload in self.insertion_payloads:
            source = sources.get(payload.material_id)
            if source is None or any(not source.start <= e.start < e.end <= source.end for e in payload.exons):
                raise ValueError("Insertion payload is outside its declared source tract")
            if not set(payload.exons) <= set(self.spans):
                raise ValueError("Payload exons must be included in the declared catalogue")



def canonical_exons(exons: Iterable[ExonSpan]) -> tuple[ExonSpan, ...]:
    return tuple(sorted(set(exons)))
