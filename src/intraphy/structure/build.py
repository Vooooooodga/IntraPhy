"""Sequence/annotation evidence to local exon configurations, before tree inference."""
from __future__ import annotations
from collections import defaultdict
from dataclasses import replace
from pathlib import Path
import re

from .alignment import FamilyAlignment, gap_supported, paired_support
from .types import Catalogue, ExonConfiguration, ExonSpan, Material, ObservationEvidence
from .material import normalize_exons
from .alternatives import structure_alternatives
from .genomic_observations import genomic_span_observation


def _components(alignment: FamilyAlignment):
    exons = [(e.id, alignment.exons[e.id].start, alignment.exons[e.id].end)
             for locus in alignment.loci for e in locus.exons]
    # A continuous indel spanning several exons must keep them in one unit.
    gaps = sorted({(m.start(), m.end()) for row in alignment.rows.values() for m in re.finditer("-+", row)
                   if any(a < m.end() and m.start() < b for _, a, b in exons)})
    intervals = sorted([(a, b, eid) for eid, a, b in exons] + [(a, b, None) for a, b in gaps],
                       key=lambda x: (x[0], x[1], x[2] or ""))
    groups, current = [], []
    left = right = None
    for a, b, eid in intervals:
        if right is None or a >= right:
            if current:
                groups.append((left, right, tuple(sorted(set(current)))))
            left, right, current = a, b, []
        else:
            right = max(right, b)
        if eid:
            current.append(eid)
    if current:
        groups.append((left, right, tuple(sorted(set(current)))))
    return groups, gaps


def _blocks(aligned_a, aligned_b, start, end):
    runs, run = [], None
    for i in range(start, end):
        supported = aligned_a[i] in "ACGT" and aligned_b[i] in "ACGT"
        if supported and run is None:
            run = i
        if not supported and run is not None:
            runs.append((run, i))
            run = None
    if run is not None:
        runs.append((run, end))
    return runs


def _record_physical_exons(family, unit, species, positive, instances, spans, start,
                           kind, local_reasons, fatal, correspondence, coordinates, relation):
    for eid in positive:
        exon, span = instances[eid], spans[eid]
        correspondence.append({"family_id": family, "unit_id": unit, "species": species,
            "exon_id": eid, "contig": exon.contig, "source_start0": exon.start, "source_end0": exon.end,
            "strand": exon.strand, "local_start0": span.start, "local_end0": span.end,
            "local_exon_count": relation, "status": kind,
            "transcripts": ";".join(exon.transcripts),
            "reasons": ";".join(sorted(set(local_reasons) | fatal))})
        coordinates.append({"family_id": family, "unit_id": unit, "species": species, "exon_id": eid,
            "alignment_start0": span.start+start, "alignment_end0": span.end+start,
            "source_contig": exon.contig, "source_start0": exon.start, "source_end0": exon.end,
            "source_strand": exon.strand, "coordinate_system": "0-based-half-open"})


def _material(rows, start, end, gaps, flank, identity):
    local = [(a, b) for a, b in gaps if start <= a < b <= end]
    overlapping = any(a < d and c < b and (a, b) != (c, d)
                      for i, (a, b) in enumerate(local) for c, d in local[i+1:])
    if overlapping:
        return (), {}, ("overlapping_indel_tracts_unresolved",)
    material = tuple(Material(f"tract_{a}_{b}", a-start, b-start) for a, b in local)
    observations = {s: [] for s in rows}
    reasons = set()
    for a, b in local:
        for species, row in rows.items():
            fragment = row[a:b]
            if set(fragment) <= set("ACGT"):
                value = 1
            elif set(fragment) == {"-"} and gap_supported(rows, species, a, b, flank, identity):
                value = 0
            else:
                value = None
                reasons.add("unqualified_gap_or_ambiguous_material")
            observations[species].append(value)
    return material, {s: tuple(v) for s, v in observations.items()}, tuple(sorted(reasons))


def build_catalogues(alignment: FamilyAlignment, *, minimum_identity: float = .7,
                     anchor_bases: int = 12, anchor_identity: float = .8,
                     observation_unit: str = "transcript_configuration"):
    if not 0 <= minimum_identity <= 1 or not 0 <= anchor_identity <= 1 or anchor_bases < 1:
        raise ValueError("Invalid explicitly declared evidence thresholds")
    if observation_unit not in {"transcript_configuration", "genomic_exon_spans"}:
        raise ValueError("Unsupported catalogue observation unit")
    components, gaps = _components(alignment)
    instances = {e.id: e for l in alignment.loci for e in l.exons}
    catalogues, correspondence, candidates, coordinates = [], [], [], []
    for number, (start, end, ids) in enumerate(components, 1):
        family = alignment.loci[0].family
        unit = f"unit_{number:04d}"
        spans = {eid: ExonSpan(alignment.exons[eid].start-start, alignment.exons[eid].end-start) for eid in ids}
        material, material_obs, reasons = _material(alignment.rows, start, end, gaps, anchor_bases, anchor_identity)
        fatal = set(reasons) & {"overlapping_indel_tracts_unresolved"}
        anomalies = {reason for eid in ids for reason in alignment.anomalies.get(eid, ())}
        fatal.update(anomalies)
        boundary_candidates = tuple(sorted(set(spans.values())))
        junctions = set()
        observed_native = {}
        for locus in alignment.loci:
            configs = []
            if observation_unit == "transcript_configuration":
                for tx, path in sorted(locus.paths.items()):
                    local = tuple(spans[eid] for eid in path if eid in spans)
                    if local:
                        if local != tuple(sorted(local)) or any(a.end >= b.start for a, b in zip(local, local[1:])):
                            fatal.add("non_collinear_or_zero_length_spacer_annotation")
                            continue
                        junctions.update((a.end, b.start) for a, b in zip(local, local[1:]))
                    configs.append((tx, local))
            else:
                local = tuple(sorted({spans[eid] for eid in ids if instances[eid].species == locus.species}))
                if not any(a.end >= b.start for a, b in zip(local, local[1:])):
                    junctions.update((a.end, b.start) for a, b in zip(local, local[1:]))
            observed_native[locus.species] = configs or [("unresolved", ())]
        proto = Catalogue(family, unit, end-start, boundary_candidates, tuple(sorted(junctions)), material,
                          exon_instances=tuple(replace(instances[eid], comparison_start=spans[eid].start, comparison_end=spans[eid].end) for eid in ids), alignment_offset=start,
                          boundary_candidates=boundary_candidates, observation_unit=observation_unit)
        observations = []
        for locus in alignment.loci:
            species, row = locus.species, alignment.rows[locus.species]
            presence = material_obs.get(species, (None,)*len(material))
            positive = tuple(eid for eid in ids if instances[eid].species == species)
            if observation_unit == "genomic_exon_spans":
                observation = genomic_span_observation(species=species, proto=proto, locus=locus,
                    ids=ids, spans=spans, instances=instances, alignment=alignment,
                    material_presence=presence, all_material_presence=material_obs, start=start,
                    minimum_identity=minimum_identity, anchor_bases=anchor_bases,
                    anchor_identity=anchor_identity)
                native_configs, kind = observation.configurations, observation.kind
                local_reasons, unknown = list(observation.reasons), observation.unknown_intervals
                if fatal:
                    kind = "excluded"
                    observation = replace(observation, kind=kind,
                        reasons=tuple(sorted(set(observation.reasons) | fatal)))
                observations.append(observation)
                relation = str(len(positive))
                _record_physical_exons(family, unit, species, positive, instances, spans, start,
                    kind, local_reasons, fatal, correspondence, coordinates, relation)
                continue
            else:
                # A native sequence tells us that absent material is absent, but not
                # whether it was never introduced or was subsequently deleted.
                native_material = tuple(1 if p == 1 else 0 for p in presence)
                native_configs = tuple(sorted({ExonConfiguration(normalize_exons(proto, native_material, exons), native_material)
                                               for _, exons in observed_native[species]}))
                unknown = tuple(sorted(set([ExonSpan(m.start(), m.end()) for m in re.finditer("[^ACGT-]+", row[start:end])]
                    + [ExonSpan(m.start, m.end) for m, p in zip(material, presence) if p is None])))
                local_reasons = []
                if not locus.paths:
                    unknown = (ExonSpan(0, end-start),)
                    local_reasons.append("locus_has_no_exon_annotation")
                if locus.partial:
                    unknown = (ExonSpan(0, end-start),)
                    local_reasons.append("partial_annotation_structure_unresolved")
                alternatives = set()
                for eid, candidate in spans.items():
                    source = instances[eid].species
                    if source == species:
                        continue
                    paired, identity = paired_support(row, alignment.rows[source], candidate.start+start, candidate.end+start)
                    target_bases = sum(b in "ACGT" for b in row[candidate.start+start:candidate.end+start])
                    source_bases = sum(b in "ACGT" for b in alignment.rows[source][candidate.start+start:candidate.end+start])
                    if paired >= min(anchor_bases, max(1, source_bases)) and identity >= minimum_identity:
                        local_reasons.append("sequence_supported")
                    elif target_bases > 0 and not positive:
                        local_reasons.append("unresolved_sequence_correspondence")
                atomic, trace = structure_alternatives(proto, alignment, species, native_configs,
                    observed_native, material_obs, minimum_identity=minimum_identity,
                    anchor_bases=anchor_bases, anchor_identity=anchor_identity)
                candidates.extend(trace)
                if atomic:
                    local_reasons.append("sequence_supported_structure_alternative")
                if unknown:
                    local_reasons.append("assembly_ambiguous_bases")
                if any(p is None for p in presence):
                    local_reasons.append("material_state_unresolved")
                if positive:
                    qualified = False
                    for other, target in alignment.rows.items():
                        if other == species:
                            continue
                        paired, identity = paired_support(row, target, start, end)
                        if paired >= min(anchor_bases, max(1, sum(b in "ACGT" for b in row[start:end]))) and identity >= minimum_identity:
                            qualified = True
                        if material_obs.get(other) and all(p == 0 for p in material_obs[other]):
                            qualified = True
                    if not qualified and len(alignment.loci) > 1:
                        local_reasons.append("unresolved_sequence_correspondence")
                kind = "coexisting" if len(native_configs) > 1 else "observed"
                if "unresolved_sequence_correspondence" in local_reasons:
                    kind = "unknown"
                elif unknown or any(p is None for p in presence) or atomic:
                    kind = "coexisting" if len(native_configs) > 1 else "partial"
                if fatal:
                    kind = "excluded"
                observations.append(ObservationEvidence(species, native_configs, kind,
                    tuple(sorted(set(local_reasons))), positive, presence, unknown,
                    tuple(sorted(alternatives)), atomic))
                counts = [len(exons) for _, exons in observed_native[species]]
                relation = "coexisting" if len(set(counts)) > 1 else str(max(counts, default=0))
            _record_physical_exons(family, unit, species, positive, instances, spans, start,
                kind, local_reasons, fatal, correspondence, coordinates, relation)
        status = "unresolved" if fatal else "qualified"
        catalogues.append(replace(proto, observations=tuple(observations), status=status,
                                  reasons=tuple(sorted(set(reasons)|fatal))))
    return tuple(catalogues), correspondence, candidates, coordinates
