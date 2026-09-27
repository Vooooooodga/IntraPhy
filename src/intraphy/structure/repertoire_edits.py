"""Experimental shared-opportunity repertoire event generation.

The named policy is a conditional finite catalogue model. It does not claim
biological completeness, transcript usage, copy lineage, or general turnover.
"""
from __future__ import annotations

from collections import defaultdict, deque
from .edits import elementary_edits
from .material import valid_configuration
from .repertoire import ExonRepertoire, delete_repertoire_material
from .repertoire_process import RepertoireEvent, RepertoireProcess
from .types import Catalogue, ExonConfiguration, ExonSpan

POLICY = "shared_coordinate_opportunities_v1"

def _rkey(r):
    return (r.material, tuple(c.key for c in r.configurations))


def _structural_events(catalogue, source):
    by = defaultdict(list)
    for member in source.configurations:
        for edit in elementary_edits(catalogue, member):
            if edit.kind in {"dna_deletion", "dna_insertion"}:
                continue
            by[(edit.kind, edit.opportunity, edit.footprint)].append((member, edit))
    outcomes = defaultdict(dict)
    for key, pairs in by.items():
        kind, opportunity, footprint = key
        affected = []
        if opportunity.startswith("exon:"):
            bits = tuple(map(int, opportunity.split(":")[1:]))
            affected = [m for m in source.configurations if any((e.start, e.end) == bits for e in m.exons)]
        elif opportunity.startswith("pair:"):
            bits = tuple(map(int, opportunity.split(":")[1:]))
            affected = [m for m in source.configurations if any((a.start,a.end,b.start,b.end)==bits for a,b in zip(m.exons,m.exons[1:]))]
        elif kind == "acceptor_shift":
            affected = [m for m in source.configurations if any(e.start == footprint[0] for e in m.exons)]
        elif kind == "donor_shift":
            affected = [m for m in source.configurations if any(e.end == footprint[0] for e in m.exons)]
        else:
            affected = [m for m in source.configurations if any(footprint[0] in (e.start,e.end) for e in m.exons)]
        mapping = {}
        for m, e in pairs:
            if m in affected:
                if m in mapping: raise ValueError("Ambiguous shared opportunity")
                mapping[m] = e
        if kind == "exonization":
            candidate = ExonSpan(*footprint)
            affected = [m for m in source.configurations if candidate not in m.exons and all(not candidate.overlaps(x) and candidate.end != x.start and candidate.start != x.end for x in m.exons)]
            mapping = {}
            for m, e in pairs:
                if m in affected:
                    if m in mapping: raise ValueError("Ambiguous shared opportunity")
                    mapping[m] = e
        if len(mapping) != len(affected):
            continue
        targets = tuple(mapping.get(m).target if m in mapping else m for m in source.configurations)
        target = ExonRepertoire(targets, source.material)
        if target != source: outcomes[(kind, opportunity)][target] = True
    return outcomes


def coupled_process(catalogue: Catalogue, seeds, *, policy: str, option_bank,
                    insertion_outcomes, provenance: str, max_states: int | None = None):
    if policy != POLICY:
        raise ValueError("An explicit supported coupling policy is required")
    if not isinstance(catalogue, Catalogue) or catalogue.status != "qualified" or not isinstance(provenance, str) or not provenance.strip():
        raise ValueError("Qualified catalogue and provenance are required")
    if max_states is not None and (type(max_states) is not int or max_states <= 0):
        raise ValueError("max_states must be positive or None")
    seeds = tuple(seeds); option_bank = tuple(option_bank); insertion_outcomes = tuple(insertion_outcomes)
    if not seeds: raise ValueError("At least one seed is required")
    for config in option_bank:
        if not isinstance(config, ExonConfiguration) or len(config.material) != len(catalogue.material) or not set(config.exons) <= set(catalogue.spans) or not valid_configuration(catalogue, config):
            raise ValueError("Invalid option-bank configuration")
    insertion_weights = {}
    for event in insertion_outcomes:
        if not isinstance(event, RepertoireEvent) or event.kind != "dna_insertion": raise ValueError("Insertion table must contain dna_insertion events")
        key = (event.source, event.target, event.kind, event.opportunity)
        if key in insertion_weights and insertion_weights[key] != event.weight: raise ValueError("Conflicting insertion weights")
        insertion_weights[key] = event.weight
    insertion_totals = defaultdict(float)
    for (source, _target, kind, opportunity), weight in insertion_weights.items():
        insertion_totals[(source, kind, opportunity)] += weight
    if any(abs(value - 1.) > 1e-12 for value in insertion_totals.values()):
        raise ValueError("Insertion opportunity weights must sum to one")
    if insertion_outcomes:
        endpoints = tuple(dict.fromkeys([e.source for e in insertion_outcomes] + [e.target for e in insertion_outcomes]))
        validated_insertions = RepertoireProcess(endpoints, insertion_outcomes, provenance, catalogue).events
    else:
        validated_insertions = ()
    queue, queued, states, events = deque(), set(), [], []
    for seed in seeds:
        if not isinstance(seed, ExonRepertoire) or len(seed.material) != len(catalogue.material) or any(not valid_configuration(catalogue,c) or not set(c.exons) <= set(catalogue.spans) for c in seed.configurations): raise ValueError("Invalid seed")
        if seed not in queued:
            if max_states is not None and len(queued) >= max_states: raise ValueError("state_space_incomplete")
            queued.add(seed); queue.append(seed)
    while queue:
        source = queue.popleft(); states.append(source)
        if max_states is not None and len(queued) > max_states: raise ValueError("state_space_incomplete")
        choices = defaultdict(set)
        for key, targets in _structural_events(catalogue, source).items():
            for target in targets: choices[key].add(target)
        for k, material in enumerate(catalogue.material):
            if source.material[k] == 1:
                target = delete_repertoire_material(catalogue, source, material.id)
                choices[("dna_deletion", material.id)].add(target)
        for option in option_bank:
            if option.material == source.material and option not in source.configurations:
                choices[("option_gain", option.key)].add(ExonRepertoire((*source.configurations, option), source.material))
        if len(source.configurations) > 1:
            for option in source.configurations:
                target = ExonRepertoire(tuple(c for c in source.configurations if c != option), source.material)
                choices[("option_loss", option.key)].add(target)
        for event in validated_insertions:
            if event.source == source:
                events.append(event)
                if event.target not in queued:
                    if max_states is not None and len(queued) >= max_states: raise ValueError("state_space_incomplete")
                    queued.add(event.target); queue.append(event.target)
        for (kind, opportunity), targets in sorted(choices.items()):
            targets = tuple(sorted((t for t in targets if t != source), key=_rkey))
            for target in sorted(targets, key=_rkey):
                events.append(RepertoireEvent(source, target, kind, opportunity, 1/len(targets)))
                if target not in queued:
                    if max_states is not None and len(queued) >= max_states: raise ValueError("state_space_incomplete")
                    queued.add(target); queue.append(target)
    process = RepertoireProcess(tuple(states), tuple(events), provenance, catalogue)
    if any(event.source not in process.index for event in validated_insertions):
        raise ValueError("Insertion table source is unreachable")
    return process
