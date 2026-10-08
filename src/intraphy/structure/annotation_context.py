"""Conservative annotation context for homologous candidate exon intervals."""
from __future__ import annotations

from .alignment import paired_support
from .types import ExonSpan


def _partial_alignment_spans(*, species, locus, alignment):
    if not locus.partial_intervals:
        return ()
    columns = alignment.columns.get(species)
    if columns is None:
        return None
    spans = []
    for genomic_start, genomic_end in locus.partial_intervals:
        if locus.strand == "-":
            sequence_start = locus.search_end - genomic_end
            sequence_end = locus.search_end - genomic_start
        else:
            sequence_start = genomic_start - (locus.search_start - 1)
            sequence_end = genomic_end - (locus.search_start - 1)
        if not 0 <= sequence_start < sequence_end <= len(columns):
            return None
        partial_columns = columns[sequence_start:sequence_end]
        if not partial_columns:
            return None
        spans.append(ExonSpan(min(partial_columns), max(partial_columns) + 1))
    return tuple(spans)


def _overlaps_any(span, other_spans):
    return any(span.start < other.end and other.start < span.end for other in other_spans)


def annotation_supported_nonexonic_interval(*, species, source_species, candidate, source_id,
        locus, alignment, start, anchor_bases, minimum_identity, anchor_identity):
    """Whether a retained interval is intronic between exons on a provided path.

    This is conditional on the supplied whole-locus annotation. It does not
    assert that the annotation includes every biological transcript.
    """
    required = ("paths", "exons", "partial", "partial_intervals",
                "partial_location_known", "strand", "search_start", "search_end")
    if any(not hasattr(locus, key) for key in required):
        return False
    if (not locus.paths or not locus.exons or locus.strand not in {"+", "-"}
            or locus.partial and not locus.partial_location_known
            or locus.partial and not locus.partial_intervals):
        return False
    if not all(hasattr(alignment, key) for key in ("loci", "rows", "exons", "anomalies", "columns")):
        return False
    row = alignment.rows.get(species)
    source = next((value for value in alignment.loci if value.species == source_species), None)
    if row is None or source is None:
        return False
    if not hasattr(source, "exons") or source_id not in {exon.id for exon in source.exons}:
        return False

    whole = ExonSpan(start + candidate.start, start + candidate.end)
    if whole.end > len(row) or not all(base in "ACGT" for base in row[whole.start:whole.end]):
        return False
    context_start, context_end = whole.start - anchor_bases, whole.end + anchor_bases
    if context_start < 0 or context_end > len(row):
        return False

    exon_ids = {exon.id for exon in locus.exons}
    if source_id not in alignment.exons:
        return False
    candidate_projection = alignment.exons[source_id]
    if candidate_projection != whole:
        return False
    if any(exon_id not in alignment.exons for exon_id in exon_ids):
        return False
    for exon_id in exon_ids:
        if alignment.exons[exon_id].overlaps(whole):
            return False

    path_pairs = []
    for path in locus.paths.values():
        if not path or any(exon_id not in exon_ids for exon_id in path):
            continue
        for left_id, right_id in zip(path, path[1:]):
            left, right = alignment.exons[left_id], alignment.exons[right_id]
            if left.end <= whole.start and whole.end <= right.start:
                path_pairs.append((left_id, right_id))
    if not path_pairs:
        return False

    # Assess candidate and path evidence separately; one ambiguous transcript
    # path does not invalidate another fully supported path.
    if source_id in alignment.anomalies:
        return False

    source_row = alignment.rows.get(source.species)
    source_whole = alignment.exons[source_id]
    if source_row is None or source_whole.end > len(source_row):
        return False
    paired, identity = paired_support(row, source_row, whole.start, whole.end)
    if paired != whole.end - whole.start or identity < minimum_identity:
        return False
    for left, right in ((context_start, whole.start), (whole.end, context_end)):
        paired, identity = paired_support(row, source_row, left, right)
        if paired != anchor_bases or identity < anchor_identity:
            return False

    partial_spans = _partial_alignment_spans(species=species, locus=locus,
                                             alignment=alignment)
    if partial_spans is None or _overlaps_any(ExonSpan(context_start, context_end), partial_spans):
        return False

    for left_id, right_id in path_pairs:
        if left_id in alignment.anomalies or right_id in alignment.anomalies:
            continue
        left, right = alignment.exons[left_id], alignment.exons[right_id]
        if (left.end - left.start < anchor_bases
                or right.end - right.start < anchor_bases):
            continue
        left_boundary = ExonSpan(left.end - anchor_bases, left.end)
        right_boundary = ExonSpan(right.start, right.start + anchor_bases)
        if (left_boundary.start < left.start or right_boundary.end > right.end
                or right_boundary.end > len(row)):
            continue
        if (not all(base in "ACGT" for base in row[left_boundary.start:left_boundary.end])
                or not all(base in "ACGT" for base in row[right_boundary.start:right_boundary.end])):
            continue
        if _overlaps_any(left_boundary, partial_spans) or _overlaps_any(right_boundary, partial_spans):
            continue
        if any(partial.start <= left.end <= partial.end
               or partial.start <= right.start <= partial.end
               for partial in partial_spans):
            continue
        return True
    return False
