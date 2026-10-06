"""Exact likelihood-only pruning over subsets of origin opportunities."""
from __future__ import annotations

import math
import numpy as np
from scipy.special import logsumexp

from ..structure.origins import origin_candidates
from .configuration_sparse import log_action, sparse_generator, _support_after_edge


def _subsets(mask):
    sub = mask
    while True:
        yield sub
        if sub == 0:
            break
        sub = (sub - 1) & mask


def _checked_action(q, length, child_value):
    value = log_action(q, length, child_value)
    if np.isnan(value).any() or np.isposinf(value).any():
        raise ArithmeticError("Sparse origin-subset action produced a nonfinite log likelihood")
    if np.isneginf(value).any():
        reachable = _support_after_edge(q, length, np.isfinite(child_value))
        if np.any(reachable & ~np.isfinite(value)):
            raise ArithmeticError(
                "Sparse origin-subset action collapsed a structurally possible likelihood")
    return value


def _checked_product(left, right):
    value = left + right
    if (np.isnan(value).any() or np.isposinf(value).any()
            or np.any(np.isfinite(left) & np.isfinite(right) & ~np.isfinite(value))):
        raise ArithmeticError("Sparse origin-subset pruning overflowed a log likelihood")
    return value


def _checked_sum(left, right):
    value = np.logaddexp(left, right)
    if (np.isnan(value).any() or np.isposinf(value).any()
            or np.any(np.isfinite(left) & np.isfinite(right) & ~np.isfinite(value))):
        raise ArithmeticError("Sparse origin-subset pruning collapsed a likelihood sum")
    return value


def likelihood_origin_dp(space, tree, tips, model, lengths, cache, *,
                         max_origins=None, sparse_templates=None):
    """Marginalize origin assignments exactly without enumerating scenarios.

    Each inside table maps a subset of tracts introduced strictly below its
    node to the statewise log likelihood. Opportunity assignments on the edge
    to a child are combined with the child's disjoint subset before that
    edge's origin-conditional CTMC action.
    """
    if max_origins is not None:
        if type(max_origins) is not int or max_origins < 1:
            raise ValueError("maximum must be a positive integer or None")
    if set(tips) != set(tree.leaf_by_label):
        raise ValueError("Exactly one likelihood vector for every tree tip is required")
    candidates = origin_candidates(space, tree, tips=tips)
    scenario_count = math.prod(map(len, candidates))
    if max_origins is not None:
        if scenario_count > max_origins:
            raise ValueError("origin_scenarios_incomplete: increase the explicit scenario limit")

    k = len(candidates)
    all_mask = (1 << k) - 1
    material_ids = tuple(m.id for m in space.catalogue.material)
    node_masks = {node: sum(1 << i for i, eligible in enumerate(candidates)
                            if node in eligible) for node in tree.parent}
    tip_logs = {}
    for label, values in tips.items():
        values = np.asarray(values, dtype=float)
        if (values.shape != (len(space.states),) or not np.isfinite(values).all()
                or (values < 0).any()):
            raise ValueError(f"Invalid observation likelihood: {label}")
        with np.errstate(divide="ignore"):
            tip_logs[label] = np.log(values)

    def edge_generator(child, edge_mask):
        allowed = tuple(material_ids[i] for i in range(k) if edge_mask & (1 << i))
        origins = {material_ids[i]: child for i in range(k) if edge_mask & (1 << i)}
        signature = (allowed, child in model.foreground)
        return cache.get_or_compute(signature, lambda: (
            sparse_templates.generator(space, model, origins, child)
            if sparse_templates is not None else sparse_generator(space, model, origins, child)))

    inside = {}
    for node in tree.postorder():
        children = tuple(tree.children.get(node, ()))
        if not children:
            inside[node] = {0: tip_logs[tree.label[node]]}
            continue
        combined = {0: np.zeros(len(space.states))}
        for child in children:
            edge_table = {}
            eligible_here = node_masks[child]
            for assigned_edge in _subsets(eligible_here):
                q = edge_generator(child, assigned_edge)
                for assigned_subtree, child_value in inside[child].items():
                    if assigned_edge & assigned_subtree:
                        continue
                    total_subset = assigned_edge | assigned_subtree
                    value = _checked_action(q, lengths[child], child_value)
                    if total_subset in edge_table:
                        edge_table[total_subset] = _checked_sum(
                            edge_table[total_subset], value)
                    else:
                        edge_table[total_subset] = value
            inside.pop(child)
            next_combined = {}
            for left_mask, left_value in combined.items():
                for right_mask, right_value in edge_table.items():
                    if left_mask & right_mask:
                        continue
                    mask = left_mask | right_mask
                    value = _checked_product(left_value, right_value)
                    if mask in next_combined:
                        next_combined[mask] = _checked_sum(next_combined[mask], value)
                    else:
                        next_combined[mask] = value
            combined = next_combined
        inside[node] = combined

    root_mask = node_masks[tree.root]
    denominator = model.origin_root_weight + len(tree.parent) - 1
    terms = []
    for assigned, values in inside[tree.root].items():
        inherited = all_mask ^ assigned
        if inherited & ~root_mask:
            continue
        required = tuple(1 if inherited & (1 << i) else 0 for i in range(k))
        valid = np.fromiter((state.material == required for state in space.states),
                            dtype=bool, count=len(space.states))
        count = int(valid.sum())
        if count == 0:
            continue
        with np.errstate(divide="ignore"):
            root_ll = float(logsumexp(values[valid]) - math.log(count))
        log_prior = (inherited.bit_count() * math.log(model.origin_root_weight)
                     - k * math.log(denominator))
        terms.append(root_ll + log_prior)
    total = float(logsumexp(terms)) if terms else -math.inf
    if math.isnan(total) or total == math.inf:
        raise ArithmeticError("Sparse origin-subset pruning produced a nonfinite log likelihood")
    return total
