"""Compact floating-point posteriors for complete configuration CTMCs.

This backend retains state and origin priors while replacing dense endpoint
matrices with exponential actions and observable-group contractions. Its model
matches the dense backend up to numerical floating-point tolerance.
"""
from __future__ import annotations

import math
from dataclasses import replace
import numpy as np
from scipy.special import logsumexp

from ..structure.origins import origin_scenarios
from ..structure.tree_context import canonical_tree
from .configuration_sparse import _possible, log_action, sparse_generator
from .configuration_compact_modes import branch_scalars, modal_rows
from .kernel_cache import KernelCache


def _log(values):
    with np.errstate(divide="ignore"):
        return np.log(values)


def _observable(state):
    exons = tuple((int(exon.start), int(exon.end)) for exon in state.exons)
    dna = tuple(int(value == 1) for value in state.material)
    return exons, dna


def _group(values):
    unique = list(dict.fromkeys(values))
    lookup = {value: i for i, value in enumerate(unique)}
    return unique, np.asarray([lookup[value] for value in values], dtype=np.int64)


def _action_columns(q, length, log_values):
    """Apply log_action while independently scaling each finite RHS column."""
    shifts = np.max(log_values, axis=0)
    valid = np.isfinite(shifts)
    result = np.full(log_values.shape, -np.inf, dtype=float)
    if valid.any():
        normalized = log_values[:, valid] - shifts[valid][None, :]
        result[:, valid] = log_action(q, length, normalized) + shifts[valid][None, :]
    return result


def _masked_columns(log_vector, groups, group_ids):
    values = np.full((len(log_vector), len(group_ids)), -np.inf, dtype=float)
    for column, group_id in enumerate(group_ids):
        mask = groups == group_id
        values[mask, column] = log_vector[mask]
    return values


def _scenario_payload(space, tree, tips, model, origins, root, lengths, generator_cache):
    generators = {}
    for _, child in tree.edges():
        allowed = tuple(sorted(material for material, origin in origins.items()
                               if origin == child))
        signature = ("Q", allowed, child in model.foreground)
        generators[child] = generator_cache.get_or_compute(signature, lambda child=child:
            sparse_generator(space, model, origins, child))
    inside, messages, siblings = {}, {}, {}
    for node in tree.postorder():
        children = tuple(tree.children.get(node, ()))
        if not children:
            inside[node] = _log(np.asarray(tips[tree.label[node]], dtype=float))
            continue
        terms = []
        for child in children:
            messages[child] = log_action(generators[child], lengths[child], inside[child])
            terms.append(messages[child])
        inside[node] = sum(terms, start=np.zeros(len(root)))
        for index, child in enumerate(children):
            siblings[child] = sum((value for j, value in enumerate(terms) if j != index),
                                  start=np.zeros(len(root)))
    log_root = _log(root.astype(float) / root.sum())
    ll = float(logsumexp(log_root + inside[tree.root]))
    if not np.isfinite(ll):
        if np.isnan(ll) or ll == np.inf:
            raise ArithmeticError("Compact pruning produced a nonfinite log likelihood")
        if _possible(tree, tips, generators, lengths, root.astype(float) / root.sum()):
            raise ArithmeticError("Compact pruning collapsed a structurally possible likelihood")
        return {"log_likelihood": ll, "generators": generators,
                "inside": inside, "outside": {}}
    outside = {tree.root: log_root}
    for parent in tree.preorder():
        for child in tree.children.get(parent, ()):
            base = outside[parent] + siblings[child]
            outside[child] = log_action(generators[child].T, lengths[child], base)
    return {"log_likelihood": ll, "generators": generators,
            "inside": inside, "outside": outside, "siblings": siblings}


def _node_marginals(payload, tree):
    ll = payload["log_likelihood"]
    return {node: np.exp(payload["outside"][node] + payload["inside"][node] - ll)
            for node in tree.preorder()}


def _scenario_key(origins):
    return ("scenario", tuple(sorted(origins.items())))


def _get_payload(cache, generator_cache, key, space, tree, tips, model, origins, root, lengths):
    return cache.get_or_compute(key, lambda: _scenario_payload(
        space, tree, tips, model, origins, root, lengths, generator_cache))


def evaluate_compact(space, tree, tips, model, *, max_origins=None,
                     branch_length_mode="supplied", counts=False,
                     block_size=32, cache_bytes=128 * 1024 * 1024,
                     tie_atol=1e-12, tie_rtol=1e-12):
    """Evaluate full-model state/node/branch posteriors without endpoint matrices.

    Expected edit counts are intentionally unsupported by this compact path.
    Floating-point tie tolerances affect which observable pairs are reported.
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
    generator_cache = KernelCache(maximum_bytes=cache_bytes // 2)
    cache = KernelCache(maximum_bytes=cache_bytes - cache_bytes // 2)
    total = -math.inf
    nodes, branches, origin_logweights = {}, {}, []
    for origins, root, log_prior in origin_scenarios(
            space, tree, max_origins, tips=tips,
            root_weight=model.origin_root_weight):
        key = _scenario_key(origins)
        payload = _get_payload(cache, generator_cache, key, space, tree, tips,
                               model, origins, root, lengths)
        ll = payload["log_likelihood"]
        if not np.isfinite(ll):
            continue
        log_weight = ll + log_prior
        next_total = float(np.logaddexp(total, log_weight))
        old_weight = math.exp(total - next_total) if np.isfinite(total) else 0.
        new_weight = math.exp(log_weight - next_total)
        total = next_total
        node_values = _node_marginals(payload, tree)
        scalar_values = branch_scalars(payload, tree, lengths,
            state_exon_groups, len(exon_configurations),
            state_dna_groups, len(dna_configurations), block_size,
            _masked_columns, _action_columns)
        for node, values in node_values.items():
            if node not in nodes:
                nodes[node] = np.zeros_like(values)
            nodes[node] *= old_weight
            nodes[node] += new_weight * values
        for edge, values in scalar_values.items():
            if edge not in branches:
                branches[edge] = {name: 0. for name in values}
            for name, value in values.items():
                branches[edge][name] *= old_weight
                branches[edge][name] += new_weight * value
        origin_logweights.append((origins, log_weight))

    base = {"log_likelihood": total, "nodes": nodes, "branches": [],
            "origins": [], "posterior_kind": "compact_observable_pairs",
            "root_prior": "uniform_valid_exon_geometries_given_origin_opportunities",
            "origin_prior": "declared_root_weight_and_unit_branch_opportunities_on_canonical_tree",
            "origin_root_weight": model.origin_root_weight,
            "tree_normalization": context.diagnostics(),
            "joint_mode_tolerance": {"absolute": float(tie_atol),
                                      "relative": float(tie_rtol)}}
    if not np.isfinite(total):
        base["nodes"] = {}
        return base
    def get_payload(origins, root):
        return _get_payload(cache, generator_cache, _scenario_key(origins),
            space, tree, tips, model, origins, root, lengths)
    modal = modal_rows(space, tree, tips, lengths, max_origins, total,
        model.origin_root_weight, nodes, configurations, state_groups,
        block_size, tie_atol, tie_rtol, get_payload, _masked_columns,
        _action_columns)
    base["nodes"] = nodes
    base["origins"] = [{"origins": origins,
                         "posterior_weight": float(math.exp(log_weight - total))}
                        for origins, log_weight in origin_logweights]
    for edge in tree.edges():
        values = branches[edge]
        base["branches"].append({"parent": edge[0], "child": edge[1],
            "probability_at_least_one_edit": values["probability_at_least_one_edit"],
            "observable_endpoint_summary": {
                "probability_exon_structure_change": values["probability_exon_structure_change"],
                "probability_dna_presence_change": values["probability_dna_presence_change"],
                **modal[edge]}})
    base["kernel_cache"] = {"hits": cache.hits, "misses": cache.misses,
                            "retained_bytes": cache.bytes,
                            "limit_bytes": cache.maximum_bytes,
                            "generator_retained_bytes": generator_cache.bytes,
                            "generator_limit_bytes": generator_cache.maximum_bytes}
    return base
