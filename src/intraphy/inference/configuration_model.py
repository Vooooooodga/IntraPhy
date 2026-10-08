"""Operation generators, explicit rate units and origin-conditional CTMCs."""
from __future__ import annotations
from dataclasses import dataclass, replace
import math
import numpy as np
from scipy.special import logsumexp
from ..structure.edits import EDIT_KINDS
from ..structure.space import StateSpace
from ..structure.origins import origin_scenarios, permitted
from ..structure.tree_context import canonical_tree
from .configuration_ctmc import likelihood, transition_matrix, expected_count, probability_any_change
from .configuration_origin_dp import likelihood_origin_dp


@dataclass(frozen=True)
class RateModel:
    rates: dict[str, float]
    scale: float = 1.
    foreground_multiplier: float = 1.
    foreground: frozenset[str] = frozenset()
    origin_root_weight: float = 1.

    def __post_init__(self):
        if not math.isfinite(self.origin_root_weight) or self.origin_root_weight <= 0:
            raise ValueError("Root opportunity weight must be finite and positive")
        if set(self.rates) != set(EDIT_KINDS):
            raise ValueError("Specify every elementary edit rate, including explicit zeros")
        if any(not np.isfinite(r) or r < 0 for r in self.rates.values()):
            raise ValueError("Edit rates must be finite and nonnegative")
        if not np.isfinite(self.scale) or self.scale <= 0:
            raise ValueError("Rate scale must be finite and positive")
        if not np.isfinite(self.foreground_multiplier) or self.foreground_multiplier < 0:
            raise ValueError("Foreground multiplier must be finite and nonnegative")


def generator(space: StateSpace, model: RateModel, origins: dict[str, str], child: str, *, include_marks=True):
    if not space.complete:
        raise ValueError("state_space_incomplete: refusing a renormalized truncated generator")
    index = space.index
    q = np.zeros((len(index), len(index)))
    marked = {kind: np.zeros_like(q) for kind in EDIT_KINDS} if include_marks else {}
    multiplier = model.scale * (model.foreground_multiplier if child in model.foreground else 1.)
    for i, j, edit in space.indexed_edits:
        if permitted(edit, origins, child):
            value = model.rates[edit.kind] * edit.weight * multiplier
            q[i, j] += value
            if include_marks:
                marked[edit.kind][i, j] += value
    np.fill_diagonal(q, -q.sum(axis=1))
    return q, marked


def evaluate_model(space: StateSpace, tree, tips, model: RateModel, *,
                   max_origins: int | None = None, posterior: bool = True,
                   counts: bool = True, branch_length_mode: str = "supplied",
                   backend: str = "dense", sparse_templates=None):
    if not space.complete:
        raise ValueError("state_space_incomplete")
    if not model.foreground <= set(tree.parent)-{tree.root}:
        raise ValueError("Foreground contains unknown/root nodes")
    if branch_length_mode not in {"supplied", "unit"}:
        raise ValueError("Unknown branch length mode")
    if backend not in {"dense", "sparse"}:
        raise ValueError("Unknown likelihood backend")
    if backend == "sparse" and (posterior or counts):
        raise ValueError("Sparse backend currently supports likelihood only; set posterior=False, counts=False")
    if sparse_templates is not None and backend != "sparse":
        raise ValueError("Sparse generator templates require the sparse backend")
    context = canonical_tree(tree, model.foreground)
    tree = context.tree
    if any(n != tree.root and len(tree.children.get(n, ())) == 1 for n in tree.parent):
        raise ValueError("Mixed rate regimes inside a subdivided branch are unsupported; "
                         "use canonical biological branches rather than changing origin opportunities")
    model = replace(model, foreground=context.foreground)
    lengths = {c: 1. if branch_length_mode == "unit" else tree.branch_length(c) for _, c in tree.edges()}
    from .kernel_cache import KernelCache
    cache = KernelCache()
    total = -math.inf
    nodes = {n: np.zeros(len(space.states)) for n in tree.parent} if posterior else {}
    endpoints = {e: np.zeros((len(space.states), len(space.states))) for e in tree.edges()} if posterior else {}
    branch_values = {e: {"probability_at_least_one_edit": 0.,
        **({"expected_edits": 0., **{f"expected_{k}": 0. for k in EDIT_KINDS}} if counts else {})}
        for e in tree.edges()} if posterior else {}
    origin_weights = []
    def base_output(log_likelihood):
        return {"log_likelihood": log_likelihood, "nodes": {}, "branches": [], "origins": [],
                "root_prior": "uniform_valid_exon_geometries_given_origin_opportunities",
                "origin_prior": "declared_root_weight_and_unit_branch_opportunities_on_canonical_tree",
                "origin_root_weight": model.origin_root_weight,
                "tree_normalization": context.diagnostics(),
                "kernel_cache": {"hits": cache.hits, "misses": cache.misses,
                                 "retained_bytes": cache.bytes, "limit_bytes": cache.maximum_bytes}}
    if backend == "sparse":
        total = likelihood_origin_dp(space, tree, tips, model, lengths, cache,
                                     max_origins=max_origins,
                                     sparse_templates=sparse_templates)
        return base_output(total)
    for origins, root, log_prior in origin_scenarios(space, tree, max_origins, tips=tips, root_weight=model.origin_root_weight):
        matrices, generators, marked = {}, {}, {}
        for _, child in tree.edges():
            signature = (tuple(k for k, v in origins.items() if v == child), lengths[child], child in model.foreground)
            def compute(child=child):
                q, b = generator(space, model, origins, child, include_marks=counts and posterior)
                return q, b, transition_matrix(q, lengths[child])
            generators[child], marked[child], matrices[child] = cache.get_or_compute(signature, compute)
        result = likelihood(tree, tips, matrices, root.astype(float)/root.sum(), posterior)
        ll = result.log_likelihood+log_prior
        origin_weights.append((origins, ll))
        if not np.isfinite(ll):
            continue
        next_total = float(np.logaddexp(total, ll))
        old_weight = float(np.exp(total-next_total))
        new_weight = float(np.exp(ll-next_total))
        total = next_total
        if not posterior:
            continue
        # Online mixture avoids retaining a dense tree likelihood per origin scenario.
        for node, values in result.nodes.items():
            nodes[node] *= old_weight
            nodes[node] += new_weight*values
        for edge, endpoint in result.endpoints.items():
            _, child = edge
            endpoints[edge] *= old_weight
            endpoints[edge] += new_weight*endpoint
            values = branch_values[edge]
            for key in values:
                values[key] *= old_weight
            values["probability_at_least_one_edit"] += new_weight*probability_any_change(
                generators[child], matrices[child], endpoint, lengths[child])
            if counts:
                values["expected_edits"] += new_weight*expected_count(
                    generators[child], matrices[child], endpoint, lengths[child])
                for kind, b in marked[child].items():
                    if b.any():
                        values[f"expected_{kind}"] += new_weight*expected_count(
                            generators[child], matrices[child], endpoint, lengths[child], b)
    out = base_output(total)
    if not posterior or not np.isfinite(total):
        return out
    out["origins"] = [{"origins": origins, "posterior_weight": float(np.exp(ll-total))}
                      for origins, ll in origin_weights if np.isfinite(ll)]
    out["nodes"] = nodes
    for edge, endpoint in endpoints.items():
        out["branches"].append({"parent": edge[0], "child": edge[1],
            "probability_different_endpoints": 1.-float(np.trace(endpoint)),
            "endpoint_probabilities": endpoint, **branch_values[edge]})
    return out
