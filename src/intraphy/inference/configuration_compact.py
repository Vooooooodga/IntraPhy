"""Compact floating-point posteriors for complete configuration CTMCs."""
from __future__ import annotations

import numpy as np
from dataclasses import replace

from ..structure.tree_context import canonical_tree
from .configuration_compact_modes import action_columns as _action_columns
from .configuration_compact_modes import masked_columns as _masked_columns
from .configuration_origin_posterior import evaluate_origin_posterior
from .kernel_cache import KernelCache


def _observable(state):
    exons = tuple((int(exon.start), int(exon.end)) for exon in state.exons)
    dna = tuple(int(value == 1) for value in state.material)
    return exons, dna


def _group(values):
    unique = list(dict.fromkeys(values))
    lookup = {value: i for i, value in enumerate(unique)}
    return unique, np.asarray([lookup[value] for value in values], dtype=np.int64)


def evaluate_compact(space, tree, tips, model, *, max_origins=None,
                     branch_length_mode="supplied", counts=False,
                     block_size=32, cache_bytes=128 * 1024 * 1024,
                     tie_atol=1e-12, tie_rtol=1e-12):
    """Evaluate compact posteriors using exact origin-subset messages.

    Inside and outside message tables retain O(nodes * 2**tracts * states)
    workspace. Joint origin weights still require one sparse scenario sweep.
    Expected edit counts remain unsupported by this compact path.
    """
    if counts:
        raise ValueError("Compact posterior currently supports counts=False only")
    if not space.complete:
        raise ValueError("state_space_incomplete")
    if branch_length_mode not in {"supplied", "unit"}:
        raise ValueError("Unknown branch length mode")
    if type(block_size) is not int or block_size < 1:
        raise ValueError("block_size must be a positive integer")
    if type(cache_bytes) is not int or cache_bytes < 0:
        raise ValueError("cache_bytes must be a nonnegative integer")
    if (not np.isfinite(tie_atol) or not np.isfinite(tie_rtol)
            or tie_atol < 0 or tie_rtol < 0):
        raise ValueError("Tie tolerances must be finite and nonnegative")
    if set(tips) != set(tree.leaf_by_label):
        raise ValueError("Exactly one likelihood vector for every tree tip is required")
    for label, values in tips.items():
        values = np.asarray(values, dtype=float)
        if (values.shape != (len(space.states),) or not np.isfinite(values).all()
                or (values < 0).any()):
            raise ValueError(f"Invalid observation likelihood: {label}")
    if not model.foreground <= set(tree.parent) - {tree.root}:
        raise ValueError("Foreground contains unknown/root nodes")
    context = canonical_tree(tree, model.foreground)
    tree = context.tree
    if any(node != tree.root and len(tree.children.get(node, ())) == 1
           for node in tree.parent):
        raise ValueError("Mixed rate regimes inside a subdivided branch are unsupported")
    model = replace(model, foreground=context.foreground)
    lengths = {child: 1. if branch_length_mode == "unit" else tree.branch_length(child)
               for _, child in tree.edges()}
    configurations, state_groups = _group([_observable(state) for state in space.states])
    exon_configurations, exon_groups = _group([value[0] for value in configurations])
    dna_configurations, dna_groups = _group([value[1] for value in configurations])
    state_exon_groups = np.asarray([exon_groups[group] for group in state_groups])
    state_dna_groups = np.asarray([dna_groups[group] for group in state_groups])
    state_exon_counts = np.asarray([len(state.exons) for state in space.states],
                                   dtype=np.int64)
    generator_cache = KernelCache(maximum_bytes=cache_bytes)
    posterior = evaluate_origin_posterior(space, tree, tips, model, lengths,
        max_origins=max_origins, block_size=block_size, tie_atol=tie_atol,
        tie_rtol=tie_rtol, configurations=configurations,
        state_groups=state_groups, exon_groups=state_exon_groups,
        exon_count=len(exon_configurations), dna_groups=state_dna_groups,
        dna_count=len(dna_configurations),
        state_exon_counts=state_exon_counts,
        generator_cache=generator_cache)
    total = posterior["log_likelihood"]
    base = {"log_likelihood": total, "nodes": posterior["nodes"], "branches": [],
            "origins": posterior["origins"], "posterior_kind": "compact_observable_pairs",
            "root_prior": "uniform_valid_exon_geometries_given_origin_opportunities",
            "origin_prior": "declared_root_weight_and_unit_branch_opportunities_on_canonical_tree",
            "origin_root_weight": model.origin_root_weight,
            "tree_normalization": context.diagnostics(),
            "joint_mode_tolerance": {"absolute": float(tie_atol),
                                      "relative": float(tie_rtol)}}
    if not np.isfinite(total):
        base["nodes"] = {}
        return base
    for edge in tree.edges():
        values = posterior["branches"][edge]
        count_pairs = [
            {"parent_count": parent_count, "child_count": child_count,
             "probability": float(values["exon_count_pair_probabilities"][parent_count,
                                                                           child_count])}
            for parent_count in range(values["exon_count_pair_probabilities"].shape[0])
            for child_count in range(values["exon_count_pair_probabilities"].shape[1])]
        base["branches"].append({"parent": edge[0], "child": edge[1],
            "probability_at_least_one_edit": values["probability_at_least_one_edit"],
            "observable_endpoint_summary": {
                "probability_exon_structure_change": values["probability_exon_structure_change"],
                "probability_dna_presence_change": values["probability_dna_presence_change"],
                "exon_count_pair_probabilities": count_pairs,
                "probability_exon_count_increase": values["probability_exon_count_increase"],
                "probability_exon_count_decrease": values["probability_exon_count_decrease"],
                "probability_exon_count_unchanged": values["probability_exon_count_unchanged"],
                **posterior["modal"][edge]}})
    # Cache counters and byte fields now describe the shared generator cache.
    base["kernel_cache"] = {"hits": generator_cache.hits,
        "misses": generator_cache.misses,
        "retained_bytes": generator_cache.bytes,
        "limit_bytes": generator_cache.maximum_bytes,
        "generator_retained_bytes": generator_cache.bytes,
        "generator_limit_bytes": generator_cache.maximum_bytes}
    return base
