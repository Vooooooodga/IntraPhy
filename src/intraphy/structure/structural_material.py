"""Condition away only sequence tracts internal to every known physical exon."""
from __future__ import annotations

from dataclasses import replace

from .types import Catalogue


_REASON = "conditioned_internal_sequence_indels"


def _partial_window_touches_material(locus, species, alignment, absolute_start, absolute_end):
    if not getattr(locus, "partial", False):
        return False
    if not getattr(locus, "partial_location_known", False):
        return True
    intervals = getattr(locus, "partial_intervals", ())
    if not intervals:
        return True
    columns = alignment.columns.get(species)
    if columns is None:
        return True
    for genomic_start, genomic_end in intervals:
        if locus.strand == "-":
            sequence_start = locus.search_end - genomic_end
            sequence_end = locus.search_end - genomic_start
        else:
            sequence_start = genomic_start - (locus.search_start - 1)
            sequence_end = genomic_end - (locus.search_start - 1)
        if sequence_start < 0 or sequence_end > len(columns) or sequence_start >= sequence_end:
            return True
        projected = columns[sequence_start:sequence_end]
        if not projected:
            return True
        left, right = min(projected), max(projected) + 1
        # Touching an endpoint can invalidate strict containment too.
        if left <= absolute_end and absolute_start <= right:
            return True
    return False


def _tract_is_internal(material, *, loci, ids, spans, instances, alignment, component_start):
    absolute_start, absolute_end = component_start + material.start, component_start + material.end
    component_ids = tuple(ids)
    if any(alignment.anomalies.get(eid) for eid in component_ids):
        return False
    for locus in loci:
        if not getattr(locus, "paths", None):
            return False
        local = tuple(sorted({spans[eid] for eid in component_ids
                              if instances[eid].species == locus.species}))
        if not local or any(left.end >= right.start for left, right in zip(local, local[1:])):
            return False
        containers = [exon for exon in local
                      if exon.start < material.start and material.end < exon.end]
        if len(containers) != 1:
            return False
        if _partial_window_touches_material(
                locus, locus.species, alignment, absolute_start, absolute_end):
            return False
    return True


def _without_material_dimension(configuration, retained_indices):
    return replace(configuration,
        material=tuple(configuration.material[i] for i in retained_indices))


def condition_internal_sequence_material(catalogue: Catalogue, *, loci, ids, spans,
                                         instances, alignment, component_start):
    """Project out only tracts strictly internal to one clean exon in every locus.

    Observation kinds, configurations' exon spans, uncertainty intervals, source
    coordinates and material-null status are otherwise left untouched.
    """
    if catalogue.observation_unit != "genomic_exon_spans" or not catalogue.material:
        return catalogue
    conditioned = tuple(i for i, material in enumerate(catalogue.material)
        if _tract_is_internal(material, loci=loci, ids=ids, spans=spans,
                              instances=instances, alignment=alignment,
                              component_start=component_start))
    if not conditioned:
        return catalogue
    retained = tuple(i for i in range(len(catalogue.material)) if i not in conditioned)
    observations = []
    for observation in catalogue.observations:
        observations.append(replace(observation,
            configurations=tuple(_without_material_dimension(config, retained)
                                 for config in observation.configurations),
            material_presence=tuple(observation.material_presence[i] for i in retained)))
    return replace(catalogue,
        material=tuple(catalogue.material[i] for i in retained),
        observations=tuple(observations),
        reasons=tuple(sorted(set(catalogue.reasons) | {_REASON})))
