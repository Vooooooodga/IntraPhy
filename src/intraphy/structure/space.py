"""Complete enumeration of a declared finite exon/boundary/material catalogue."""
from __future__ import annotations

from dataclasses import dataclass, replace, field
from functools import cached_property
from bisect import bisect_right
from itertools import product
import math

from .types import Catalogue, ElementaryEdit, ExonConfiguration, ExonSpan
from .material import normalize_exons, valid_configuration
from .edits import elementary_edits


@dataclass(frozen=True)
class StateSpace:
    catalogue: Catalogue
    states: tuple[ExonConfiguration, ...]
    edits: tuple[ElementaryEdit, ...]
    complete: bool
    reason: str
    limit: int | None
    diagnostics: dict = field(default_factory=dict, compare=False)

    @cached_property
    def index(self) -> dict[ExonConfiguration, int]:
        return {state: i for i, state in enumerate(self.states)}

    @cached_property
    def indexed_edits(self):
        index = self.index
        return tuple((index[e.source], index[e.target], e) for e in self.edits)


def close_span_catalogue(c: Catalogue, max_spans: int | None = None) -> Catalogue:
    """Close *coordinates*, not an unconstrained universe of arbitrary sequence.

    Include fused and shifted ancestral intermediates and material-trimmed ends.
    Candidate endpoints are declared by observed exon ends and material tracts.
    Completion here never claims to discover all historical exons.
    """
    starts = {e.start for e in c.spans}
    ends = {e.end for e in c.spans}
    for d, a in c.junctions:
        ends.add(d)
        starts.add(a)
    for m in c.material:
        starts.add(m.end)
        ends.add(m.start)
    spans = tuple(ExonSpan(a, b) for a in sorted(starts) for b in sorted(ends) if a < b)
    if max_spans is not None and (type(max_spans) is not int or max_spans < 1):
        raise ValueError("max_spans must be a positive integer or None")
    if max_spans is not None and len(spans) > max_spans:
        raise ValueError("state_space_incomplete: exon-boundary catalogue exceeds max_spans")
    return replace(c, spans=spans, boundary_candidates=c.boundary_candidates or c.spans)


def estimate_geometry_count(spans):
    """Exact interval-scheduling DP count, including the empty configuration."""
    ordered = tuple(sorted(spans))
    starts = [e.start for e in ordered]
    counts = [1]*(len(ordered)+1)
    for i in range(len(ordered)-1, -1, -1):
        following = bisect_right(starts, ordered[i].end)
        counts[i] = counts[i+1]+counts[following]
    return counts[0]


def enumerate_space(catalogue: Catalogue, max_states: int | None = None, *, max_spans: int | None = None) -> StateSpace:
    if max_states is not None and (type(max_states) is not int or max_states < 1):
        raise ValueError("max_states must be a positive integer or None")
    diagnostics = {"state_limit": max_states, "finite_declared_catalogue_only": True}
    def failure(c, reason):
        return StateSpace(c, (), (), False, reason, max_states, diagnostics)
    try:
        c = close_span_catalogue(catalogue, max_spans)
    except ValueError as exc:
        return failure(catalogue, str(exc))
    geometries_count = estimate_geometry_count(c.spans)
    diagnostics.update(candidate_span_count=len(c.spans),
        geometry_count_exact=geometries_count, material_combinations=3**len(c.material),
        state_count_upper_bound=geometries_count*3**len(c.material),
        dense_matrix_bytes_at_limit=None if max_states is None else 8*max_states**2)
    if max_states is not None and geometries_count > max_states:
        return failure(c, "state_space_incomplete: exact geometry count exceeds state limit")
    geometries = []
    # Complete enumeration, accelerated by the preflight count and sorted starts.
    starts = [e.start for e in c.spans]
    def visit(current, begin):
        geometries.append(current)
        for i in range(begin, len(c.spans)):
            e = c.spans[i]
            visit(current+(e,), bisect_right(starts, e.end))
    visit((), 0)
    states = set()
    span_set = frozenset(c.spans)
    for material in product((0, 1, 2), repeat=len(c.material)):
        for exons in geometries:
            normalized = normalize_exons(c, tuple(material), exons)
            if not set(normalized) <= span_set:
                return failure(c, "state_space_incomplete: unrepresented material boundary")
            states.add(ExonConfiguration(normalized, tuple(material)))
            if max_states is not None and len(states) > max_states:
                return failure(c, "state_space_incomplete: material/configuration state limit")
    ordered = tuple(sorted(states))
    edits = tuple(e for state in ordered for e in elementary_edits(c, state))
    if any(e.target not in states for e in edits):
        return failure(c, "state_space_incomplete: transition leaves catalogue")
    diagnostics.update(actual_state_count=len(ordered), elementary_edges=len(edits),
                       dense_matrix_bytes=8*len(ordered)**2)
    return StateSpace(c, ordered, edits, True, "candidate_catalogue_closed", max_states, diagnostics)
