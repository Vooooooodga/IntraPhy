"""Construct species observations from physical exon intervals, not isoforms."""
from __future__ import annotations

import re
from .alignment import gap_supported, paired_support
from .types import ExonConfiguration, ExonSpan, ObservationEvidence


def genomic_span_observation(*, species, proto, locus, ids, spans, instances, alignment,
                             material_presence, all_material_presence, start, minimum_identity, anchor_bases,
                             anchor_identity):
    row = alignment.rows[species]
    positive = tuple(eid for eid in ids if instances[eid].species == species)
    local = tuple(sorted({spans[eid] for eid in positive}))
    material = tuple(1 if value == 1 else 0 for value in material_presence)
    reasons = set()
    localized_partial_applied = False
    unknown = [ExonSpan(match.start(), match.end())
               for match in re.finditer("[^ACGT-]+", row[start:start + proto.length])]
    for genomic_start, genomic_end in getattr(locus, "partial_intervals", ()):
        if locus.strand == "-":
            sequence_start = locus.search_end - genomic_end
            sequence_end = locus.search_end - genomic_start
        else:
            sequence_start = genomic_start - (locus.search_start - 1)
            sequence_end = genomic_end - (locus.search_start - 1)
        columns = alignment.columns[species][sequence_start:sequence_end]
        if columns:
            left = max(start, min(columns))
            right = min(start + proto.length, max(columns) + 1)
            if left < right:
                unknown.append(ExonSpan(left - start, right - start))
                localized_partial_applied = True
    unknown.extend(ExonSpan(m.start, m.end) for m, value in
                   zip(proto.material, material_presence) if value is None)

    def material_covers(candidate):
        intervals = sorted((m.start, m.end) for m, value in
                           zip(proto.material, material_presence)
                           if value == 0 and m.start < candidate.end and candidate.start < m.end)
        cursor = candidate.start
        for left, right in intervals:
            if left > cursor:
                return False
            cursor = max(cursor, right)
            if cursor >= candidate.end:
                return True
        return False

    def absent_material_covers_exons(other_presence):
        intervals = sorted((m.start, m.end) for m, value in
                           zip(proto.material, other_presence) if value == 0)
        for exon in local:
            cursor = exon.start
            for left, right in intervals:
                if left > cursor:
                    break
                cursor = max(cursor, right)
                if cursor >= exon.end:
                    break
            if cursor < exon.end:
                return False
        return bool(local)

    def qualified_absence(candidate, source):
        a, b = candidate.start + start, candidate.end + start
        fragment = row[a:b]
        source_bases = sum(base in "ACGT" for base in alignment.rows[source][a:b])
        return (source_bases > 0 and set(fragment) == {"-"} and
                (material_covers(candidate) or gap_supported(alignment.rows, species, a, b,
                                                              anchor_bases, anchor_identity)))

    if not positive:
        reasons.add("annotation_missing_or_no_physical_exon_spans")
        absent = [qualified_absence(spans[eid], instances[eid].species) for eid in ids]
        unlocalized_partial = (locus.partial and
                               not getattr(locus, "partial_location_known", False))
        if absent and all(absent) and not unlocalized_partial and not unknown:
            return ObservationEvidence(species, (ExonConfiguration((), material),), "observed",
                ("all_candidate_spans_qualified_absent",), (), material_presence)
        reasons.add("genomic_absence_unqualified")
        return ObservationEvidence(species, (), "unknown", tuple(sorted(reasons)), (),
                                   material_presence, tuple(sorted(set(unknown))))

    # Keep source rows untouched; overlapping/abutting boundaries make only this
    # tip's span geometry unresolved.
    if any(left.end >= right.start for left, right in zip(local, local[1:])):
        reasons.add("overlapping_or_adjacent_physical_exon_spans")
        return ObservationEvidence(species, (), "unknown", tuple(sorted(reasons)), positive,
                                   material_presence, tuple(sorted(set(unknown))))

    for eid in ids:
        if instances[eid].species == species:
            continue
        candidate = spans[eid]
        if any(candidate.start < exon.end and exon.start < candidate.end for exon in local):
            continue
        source = instances[eid].species
        if qualified_absence(candidate, source):
            continue
        a, b = candidate.start + start, candidate.end + start
        paired, identity = paired_support(row, alignment.rows[source], a, b)
        target_bases = sum(base in "ACGT" for base in row[a:b])
        if target_bases > 0 or paired == 0 or identity < minimum_identity:
            unknown.append(candidate)
            reasons.add("projected_span_structure_unqualified")

    if (locus.partial and not getattr(locus, "partial_location_known", False)) or not locus.paths:
        unknown.append(ExonSpan(0, proto.length))
        reasons.add("partial_or_missing_annotation_structure")
    elif localized_partial_applied:
        reasons.add("localized_partial_annotation_boundary")

    qualified_placement = False
    for other, target in alignment.rows.items():
        if other == species:
            continue
        paired, identity = paired_support(row, target, start, start + proto.length)
        bases = sum(base in "ACGT" for base in row[start:start + proto.length])
        if paired >= min(anchor_bases, max(1, bases)) and identity >= minimum_identity:
            qualified_placement = True
        other_presence = all_material_presence.get(other, ())
        if other_presence and absent_material_covers_exons(other_presence):
            qualified_placement = True
    if len(alignment.loci) > 1 and not qualified_placement:
        unknown.append(ExonSpan(0, proto.length))
        reasons.add("unresolved_sequence_correspondence")

    if unknown:
        reasons.add("local_physical_exons_retained_outside_unknown_windows")
    kind = "partial" if unknown else "observed"
    return ObservationEvidence(species, (ExonConfiguration(local, material),), kind,
        tuple(sorted(reasons)), positive, material_presence, tuple(sorted(set(unknown))))
