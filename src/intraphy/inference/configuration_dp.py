"""Generalized Sankoff recursion for finite exon configurations.

One numerical owner for arbitrary state count and branch-specific edit costs.
All optimal marginal placements and one compatible joint witness are returned.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np


class SankoffWorkspace:
    """Bounded reuse of inside messages for one fixed tree/state space."""
    def __init__(self):
        self._messages = {}
        self.cache_hits = 0
        self.cache_computes = 0

    def message(self, node, cost, inside):
        signature = (np.asarray(cost).shape, np.asarray(cost).dtype.str, np.asarray(cost).tobytes(),
                     np.asarray(inside).dtype.str, np.asarray(inside).tobytes())
        item = self._messages.get(node)
        if item is not None and item[0] == signature:
            self.cache_hits += 1
            return item[1].copy()
        self.cache_computes += 1
        return None

    def save_message(self, node, cost, inside, value):
        signature = (np.asarray(cost).shape, np.asarray(cost).dtype.str, np.asarray(cost).tobytes(),
                     np.asarray(inside).dtype.str, np.asarray(inside).tobytes())
        self._messages[node] = (signature, np.asarray(value).copy())


@dataclass(frozen=True)
class ConfigurationParsimony:
    cost: float
    nodes: dict[str, tuple[int, ...]]
    pairs: dict[tuple[str, str], tuple[tuple[int, int], ...]]
    witness: dict[str, int]


def sankoff(tree, tips: dict[str, np.ndarray], costs: dict[str, np.ndarray],
            root_allowed: np.ndarray | None = None, tol: float = 1e-9,
            workspace: SankoffWorkspace | None = None) -> ConfigurationParsimony:
    if not tips:
        raise ValueError("At least one tip weight vector is required")
    n = len(next(iter(tips.values())))
    if not n or set(tips) != set(tree.leaf_by_label):
        raise ValueError("Exactly one state vector per tree tip is required")
    for label, vector in tips.items():
        if len(vector) != n or not np.all(np.isfinite(vector)) or np.any(vector < 0):
            raise ValueError(f"Invalid tip likelihood vector: {label}")
    for _, child in tree.edges():
        c = np.asarray(costs[child])
        if c.shape != (n, n) or np.isnan(c).any() or (c < 0).any():
            raise ValueError("Edit costs must be square, nonnegative and without NaN")
    if workspace is not None:
        edge_children = {child for _, child in tree.edges()}
        workspace._messages = {child: item for child, item in workspace._messages.items()
                               if child in edge_children}
    if root_allowed is not None and np.asarray(root_allowed).shape != (n,):
        raise ValueError("Root constraint dimension differs from state space")
    root_cost = np.zeros(n) if root_allowed is None else np.where(root_allowed, 0., np.inf)
    inside, messages, siblings = {}, {}, {}
    for node in tree.postorder():
        children = tuple(tree.children.get(node, ()))
        if not children:
            inside[node] = np.where(tips[tree.label[node]] > 0, 0., np.inf)
        else:
            child_values = []
            for child in children:
                messages[child] = (workspace.message(child, costs[child], inside[child])
                                   if workspace else None)
                if messages[child] is None:
                    messages[child] = np.min(costs[child] + inside[child][None, :], axis=1)
                    if workspace:
                        workspace.save_message(child, costs[child], inside[child], messages[child])
                child_values.append(messages[child])
            inside[node] = sum(child_values, start=np.zeros(n))
            for j, child in enumerate(children):
                # Avoid inf-inf; missing/impossible states remain impossible.
                siblings[child] = sum((v for k, v in enumerate(child_values) if k != j), start=np.zeros(n))
    total = inside[tree.root] + root_cost
    minimum = float(np.min(total))
    if not np.isfinite(minimum):
        return ConfigurationParsimony(minimum, {}, {}, {})
    outside = {tree.root: root_cost}
    pairs = {}
    for parent in tree.preorder():
        for child in tree.children.get(parent, ()):
            base = outside[parent] + siblings[child]
            outside[child] = np.min(base[:, None] + costs[child], axis=0)
            score = base[:, None] + costs[child] + inside[child][None, :]
            pair_indices = np.argwhere(np.isclose(score, minimum, rtol=0, atol=tol))
            pairs[(parent, child)] = tuple(map(tuple, pair_indices.tolist()))
    nodes = {v: tuple(np.flatnonzero(np.isclose(inside[v] + outside[v], minimum, rtol=0, atol=tol)).tolist())
             for v in tree.preorder()}
    witness = {tree.root: nodes[tree.root][0]}
    for parent in tree.preorder():
        for child in tree.children.get(parent, ()):
            value = costs[child][witness[parent], :] + inside[child]
            best = np.min(value)
            witness[child] = int(np.flatnonzero(np.isclose(value, best, rtol=0, atol=tol))[0])
    return ConfigurationParsimony(minimum, nodes, pairs, witness)
