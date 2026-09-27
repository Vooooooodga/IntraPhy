"""Sparse CTMC likelihoods for evidence-conditioned locus processes."""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Mapping

import numpy as np
from scipy.sparse import coo_matrix, bmat, csr_matrix
from scipy.sparse.linalg import expm_multiply
from intraphy.structure.locus_observations import observation_emission


@dataclass(frozen=True)
class LocusLikelihoodResult:
    log_likelihood: float
    node_posteriors: Mapping[str, np.ndarray] | None = None
    branch_event_counts: Mapping[str, Mapping[str, float]] | None = None


def process_generator(process, rates: Mapping[str, float]) -> csr_matrix:
    """Build a row generator; distinct edge/opportunity intensities add."""
    n = len(process.states)
    groups = {event.rate_group for event in process.catalogue.opportunities}
    groups.update(edge.rate_group for edge in process.edges)
    missing = groups - set(rates)
    if missing:
        raise ValueError(f"Missing rate groups: {sorted(missing)}")
    for group in rates:
        raw = rates[group]
        if isinstance(raw, (bool, np.bool_)) or not isinstance(raw, (int, float, np.number)):
            raise ValueError(f"Rate {group!r} must be numeric, not boolean or text")
        value = float(raw)
        if not math.isfinite(value) or value < 0:
            raise ValueError(f"Rate {group!r} must be finite and nonnegative")
    rows, cols, data = [], [], []
    outgoing = np.zeros(n, dtype=float)
    for edge in process.edges:
        if edge.rate_group not in rates:
            raise ValueError(f"Missing rate for group {edge.rate_group!r}")
        rate = float(rates[edge.rate_group])
        if not math.isfinite(rate) or rate < 0:
            raise ValueError(f"Rate {edge.rate_group!r} must be finite and nonnegative")
        value = rate * float(edge.weight)
        if not math.isfinite(value):
            raise ValueError("Transition intensity is non-finite")
        if value:
            rows.append(edge.source); cols.append(edge.target); data.append(value)
            outgoing[edge.source] += value
    rows.extend(range(n)); cols.extend(range(n)); data.extend(-outgoing)
    q = coo_matrix((data, (rows, cols)), shape=(n, n), dtype=float).tocsr()
    q.sum_duplicates()
    if not np.all(np.isfinite(q.data)):
        raise ValueError("Generator contains non-finite values")
    return q


def _vector(process, value, name, sensitivity=None, specificity=None, allow_zero=False):
    if hasattr(value, "material") and hasattr(value, "features"):
        result = observation_emission(process.catalogue, process.states, value,
                                      sensitivity=sensitivity, specificity=specificity)
    else:
        result = np.asarray(value, dtype=float)
    if result.shape != (len(process.states),) or not np.all(np.isfinite(result)) or np.any(result < 0):
        raise ValueError(f"{name} must be a finite nonnegative vector with one entry per state")
    if np.any(result > 1.0 + 1e-12):
        raise ValueError(f"{name} emission probabilities must not exceed one")
    scale = float(np.max(result))
    if scale <= 0 and not allow_zero:
        raise ValueError(f"{name} has zero total support")
    if scale <= 0:
        return result, -math.inf
    return result / scale, math.log(scale)


def _action(matrix, vector, branch_length):
    if not math.isfinite(branch_length) or branch_length < 0:
        raise ValueError("Tree branch lengths must be finite and nonnegative")
    if branch_length == 0:
        return np.asarray(vector, dtype=float).copy()
    result = np.asarray(expm_multiply(matrix * branch_length, vector), dtype=float)
    if not np.all(np.isfinite(result)):
        raise FloatingPointError("Sparse exponential action returned non-finite values")
    # Sparse exponential actions can produce tiny negative round-off only.
    tolerance = 1e-12 * max(1.0, float(np.max(np.abs(result))))
    if np.min(result) < -tolerance:
        raise FloatingPointError("Sparse exponential action produced a negative likelihood")
    result[(result < 0) & (result >= -tolerance)] = 0.0
    return result


def evaluate_locus(process, tree, tips: Mapping[str, object], root_prior,
                   rates: Mapping[str, float], *, posterior=True, counts=False,
                   feature_sensitivity=None, feature_specificity=None):
    """Evaluate one connected locus by sparse-vector pruning.

    Tip values may be LocusObservation objects or explicit emission vectors.
    `root_prior` is an explicit probability vector; no rate-derived root model
    is imposed. The returned node marginals are conditional on all tips.
    """
    q = process_generator(process, rates)
    n = q.shape[0]
    raw_prior = np.asarray(root_prior, dtype=float)
    if (np.asarray(root_prior).dtype.kind == "b" or raw_prior.shape != (n,) or not np.all(np.isfinite(raw_prior))
            or np.any(raw_prior < 0) or not math.isclose(float(raw_prior.sum()), 1.0,
                                                          rel_tol=0.0, abs_tol=1e-12)):
        raise ValueError("root_prior must be a finite probability vector that sums to one")
    prior = raw_prior.copy()
    emissions, emission_scales = {}, {}
    expected_labels = {tree.label[node] for node in tree.leaves}
    if set(tips) != expected_labels:
        raise ValueError(f"Tip labels must exactly match tree leaves; expected {sorted(expected_labels)}")
    for node in tree.postorder():
        if tree.children.get(node):
            continue
        label = tree.label[node]
        raw = tips[label]
        emissions[node], emission_scales[node] = _vector(
            process, raw, f"tip {label!r}", feature_sensitivity, feature_specificity,
            allow_zero=True)
        if not np.any(emissions[node]):
            return LocusLikelihoodResult(-math.inf, None, None)
    inside, log_scales, upward = {}, {}, {}
    for node in tree.postorder():
        children = tree.children.get(node, ())
        if not children:
            inside[node] = emissions[node]
            log_scales[node] = emission_scales[node]
            continue
        value = np.ones(n, dtype=float)
        log_scale = 0.0
        for child in children:
            child_value = _action(q, inside[child], tree.branch_length(child))
            upward[child] = child_value
            factor = float(np.max(child_value))
            if factor == 0:
                return LocusLikelihoodResult(-math.inf, None, None)
            if factor < 0 or not math.isfinite(factor):
                raise FloatingPointError(f"Invalid pruning likelihood on branch to {child!r}")
            value *= child_value / factor
            log_scale += log_scales[child] + math.log(factor)
            factor = float(np.max(value))
            if factor == 0:
                return LocusLikelihoodResult(-math.inf, None, None)
            if factor < 0 or not math.isfinite(factor):
                raise FloatingPointError("Pruning product has negative or non-finite mass")
            value /= factor
            log_scale += math.log(factor)
        inside[node], log_scales[node] = value, log_scale
    root = tree.root
    likelihood = float(prior @ inside[root])
    if likelihood == 0:
        return LocusLikelihoodResult(-math.inf, None, None)
    if likelihood < 0 or not math.isfinite(likelihood):
        raise FloatingPointError("Data likelihood is negative or non-finite")
    log_likelihood = math.log(likelihood) + log_scales[root]
    if not math.isfinite(log_likelihood):
        raise FloatingPointError("Log likelihood is non-finite")
    if posterior or counts:
        node_posteriors, outside = _posterior_messages(tree, q, prior, inside, upward)
    else:
        node_posteriors, outside = None, None
    expected = _branch_counts(tree, process, q, rates, inside, outside, upward) if counts else None
    return LocusLikelihoodResult(log_likelihood, node_posteriors if posterior else None, expected)


def _posterior_messages(tree, q, prior, inside, upward):
    outside = {tree.root: prior.copy()}
    result = {}
    q_transpose = q.T.tocsr()
    for node in tree.preorder():
        marginal = outside[node] * inside[node]
        total = float(marginal.sum())
        if total <= 0 or not math.isfinite(total):
            raise FloatingPointError(f"Node posterior is undefined at {node!r}")
        result[node] = marginal / total
        children = tuple(tree.children.get(node, ()))
        for child in children:
            context = outside[node].copy()
            for sibling in children:
                if sibling != child:
                    context *= upward[sibling]
                    scale = float(np.max(context))
                    if scale <= 0 or not math.isfinite(scale):
                        raise FloatingPointError(f"Outside context is undefined above {child!r}")
                    context /= scale
            propagated = _action(q_transpose, context, tree.branch_length(child))
            scale = float(propagated.sum())
            if scale <= 0 or not math.isfinite(scale):
                raise FloatingPointError(f"Outside likelihood is undefined below {child!r}")
            outside[child] = propagated / scale
    return result, outside


def _branch_counts(tree, process, q, rates, inside, outside, upward):
    """Expected per-edge-branch jumps, grouped by rate class."""
    n = q.shape[0]
    mark_keys = {(event.id, event.rate_group) for event in process.catalogue.opportunities}
    mark_keys.update((edge.opportunity_id, edge.rate_group) for edge in process.edges)
    opportunity_groups = {}
    group_opportunities = {}
    for opportunity_id, rate_group in mark_keys:
        mark_name = f"opportunity:{opportunity_id}:{rate_group}"
        opportunity_groups.setdefault(opportunity_id, []).append(mark_name)
        group_opportunities.setdefault(rate_group, []).append(
            mark_name)
    contexts = []
    for parent, child in tree.edges():
        siblings = tree.children[parent]
        context = outside[parent].copy()
        for sibling in siblings:
            if sibling != child:
                factor = upward[sibling]
                context *= factor
                scale = float(np.max(context))
                if scale <= 0 or not math.isfinite(scale):
                    raise FloatingPointError(f"Outside context is undefined on branch {child!r}")
                context /= scale
        child_inside = inside[child]
        branch_length = tree.branch_length(child)
        denominator = float(context @ upward[child])
        if denominator <= 0 or not math.isfinite(denominator):
            raise FloatingPointError(f"Cannot normalize event counts on branch {child!r}")
        contexts.append((child, context, inside[child], branch_length, denominator))
    result = {child: {"opportunity:" + opportunity_id: 0.0
                      for opportunity_id in opportunity_groups}
              for _parent, child in tree.edges()}
    for group in group_opportunities:
        for branch in result.values():
            branch["rate_group:" + group] = 0.0
    zero = csr_matrix((n, n), dtype=float)
    for opportunity_id, rate_group in sorted(mark_keys):
        rows, cols, data = [], [], []
        for edge in process.edges:
            if edge.opportunity_id == opportunity_id and edge.rate_group == rate_group:
                value = float(rates[edge.rate_group]) * float(edge.weight)
                if value:
                    rows.append(edge.source); cols.append(edge.target); data.append(value)
        mark = coo_matrix((data, (rows, cols)), shape=(n, n), dtype=float).tocsr()
        if not mark.nnz:
            continue
        block = bmat([[q, mark], [zero, q]], format="csr")
        mark_name = f"opportunity:{opportunity_id}:{rate_group}"
        for child, context, child_inside, branch_length, denominator in contexts:
            action = _action(block, np.concatenate((np.zeros(n), child_inside)), branch_length)
            value = float(context @ action[:n]) / denominator
            if value < 0 and value > -1e-10:
                value = 0.0
            if value < 0 or not math.isfinite(value):
                raise FloatingPointError("Expected event count is negative or non-finite")
            result[child]["opportunity:" + opportunity_id] += value
            result[child]["rate_group:" + rate_group] += value
    return result
