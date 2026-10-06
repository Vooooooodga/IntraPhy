"""Subset outside messages and the exhaustive joint-origin likelihood pass."""
from __future__ import annotations

import math
import numpy as np
from scipy.special import logsumexp

from ..structure.origins import origin_scenarios
from .configuration_origin_dp import (
    _checked_action, _checked_product, _checked_sum, _subsets,
    origin_subset_tables,
)
from .configuration_sparse import likelihood_only


def _convolve(left, right):
    result = {}
    for left_mask, left_value in left.items():
        for right_mask, right_value in right.items():
            if left_mask & right_mask:
                continue
            mask = left_mask | right_mask
            value = _checked_product(left_value, right_value)
            result[mask] = (_checked_sum(result[mask], value)
                            if mask in result else value)
    return result


def _edge_table(child, inside, node_masks, edge_generator, lengths):
    table = {}
    for edge_mask in _subsets(node_masks[child]):
        q = edge_generator(child, edge_mask)
        for below, child_value in inside[child].items():
            if edge_mask & below:
                continue
            mask = edge_mask | below
            acted = _checked_action(q, lengths[child], child_value)
            table[mask] = (_checked_sum(table[mask], acted)
                           if mask in table else acted)
    return table


def _outside_messages(tree, inside, root_outside, node_masks,
                       edge_generator, lengths, nstates):
    outside = {tree.root: root_outside}
    edge_messages = {}
    identity = {0: np.zeros(nstates)}
    for parent in tree.preorder():
        children = tuple(tree.children.get(parent, ()))
        if not children:
            continue
        edge_tables = {child: _edge_table(child, inside, node_masks,
                         edge_generator, lengths) for child in children}
        prefix = [identity]
        for child in children:
            prefix.append(_convolve(prefix[-1], edge_tables[child]))
        suffix = [None] * (len(children) + 1)
        suffix[-1] = identity
        for index in range(len(children) - 1, -1, -1):
            suffix[index] = _convolve(edge_tables[children[index]], suffix[index + 1])
        for index, child in enumerate(children):
            siblings = _convolve(prefix[index], suffix[index + 1])
            edge_outside = {}
            for child_mask in edge_tables[child]:
                value = None
                for sibling_mask, sibling_value in siblings.items():
                    if child_mask & sibling_mask:
                        continue
                    parent_mask = child_mask | sibling_mask
                    if parent_mask not in outside[parent]:
                        continue
                    term = _checked_product(outside[parent][parent_mask], sibling_value)
                    value = _checked_sum(value, term) if value is not None else term
                if value is not None:
                    edge_outside[child_mask] = value
            edge_messages[child] = edge_outside
            child_outside = {}
            for below in inside[child]:
                value = None
                for edge_mask in _subsets(node_masks[child]):
                    if edge_mask & below:
                        continue
                    combined = edge_mask | below
                    if combined not in edge_outside:
                        continue
                    q = edge_generator(child, edge_mask)
                    acted = _checked_action(q.T, lengths[child], edge_outside[combined])
                    value = _checked_sum(value, acted) if value is not None else acted
                if value is not None:
                    child_outside[below] = value
            outside[child] = child_outside
    return outside, edge_messages


def _node_marginals(tree, inside, outside, total, nstates):
    nodes = {}
    for node in tree.preorder():
        log_mass = np.full(nstates, -np.inf)
        for subset, below in inside[node].items():
            above = outside[node].get(subset)
            if above is not None:
                log_mass = _checked_sum(log_mass, _checked_product(above, below))
        nodes[node] = np.exp(log_mass - total)
    return nodes


def _joint_origins(space, tree, tips, model, lengths, max_origins,
                   edge_generator, material_ids, total):
    rows = []
    for origins, root, log_prior in origin_scenarios(
            space, tree, max_origins, tips=tips,
            root_weight=model.origin_root_weight):
        generators = {}
        for _, child in tree.edges():
            edge_mask = sum(1 << i for i, material in enumerate(material_ids)
                            if origins[material] == child)
            generators[child] = edge_generator(child, edge_mask)
        root_prior = root.astype(float)
        root_prior /= root_prior.sum()
        conditional = likelihood_only(tree, tips, generators, lengths, root_prior)
        log_weight = conditional + log_prior
        if np.isfinite(log_weight):
            rows.append({"origins": origins,
                         "posterior_weight": float(math.exp(log_weight - total))})
    if not rows:
        raise ArithmeticError(
            "Origin-subset posterior has finite likelihood but no finite origin scenario")
    return rows


def evaluate_origin_posterior(space, tree, tips, model, lengths, *,
                              max_origins, block_size, tie_atol, tie_rtol,
                              configurations, state_groups, exon_groups,
                              exon_count, dna_groups, dna_count,
                              state_exon_counts, generator_cache):
    """Compute compact posteriors from one subset outside pass and one origin pass."""
    from .configuration_compact_modes import branch_scalars, modal_rows

    tables = origin_subset_tables(space, tree, tips, model, lengths,
        generator_cache, max_origins=max_origins, retain_tables=True)
    total = tables["log_likelihood"]
    if not np.isfinite(total):
        return {"log_likelihood": total, "nodes": {}, "branches": {},
                "modal": {}, "origins": []}
    inside = tables["inside"]
    node_masks = tables["node_masks"]
    edge_generator = tables["edge_generator"]
    outside, edge_messages = _outside_messages(tree, inside,
        tables["root_outside"], node_masks, edge_generator, lengths,
        len(space.states))
    nodes = _node_marginals(tree, inside, outside, total, len(space.states))
    branches = branch_scalars(tree, lengths, inside, edge_messages,
        node_masks, edge_generator, exon_groups, exon_count,
        dna_groups, dna_count, state_exon_counts, block_size, total)
    modal = modal_rows(tree, lengths, inside, edge_messages, nodes,
        node_masks, edge_generator, configurations, state_groups,
        block_size, tie_atol, tie_rtol, total)
    origins = _joint_origins(space, tree, tips, model, lengths, max_origins,
        edge_generator, tables["material_ids"], total)
    return {"log_likelihood": total, "nodes": nodes, "branches": branches,
            "modal": modal, "origins": origins}
