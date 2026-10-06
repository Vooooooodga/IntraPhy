"""Exact sparse-generator likelihood evaluation without dense transition matrices."""
from __future__ import annotations

import math
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply
from scipy.special import logsumexp

from ..structure.origins import permitted
from .configuration_uniformization import log_uniformization_action


def sparse_generator(space, model, origins, child):
    """Build the same complete-state row generator as configuration_model.generator."""
    if not space.complete:
        raise ValueError("state_space_incomplete: refusing a renormalized truncated generator")
    rows, cols, values = [], [], []
    multiplier = model.scale * (model.foreground_multiplier
                                if child in model.foreground else 1.)
    for i, j, edit in space.indexed_edits:
        if permitted(edit, origins, child):
            value = model.rates[edit.kind] * edit.weight * multiplier
            if not np.isfinite(value) or value < 0:
                raise ArithmeticError("Sparse generator contains an invalid transition rate")
            if value:
                rows.append(i)
                cols.append(j)
                values.append(value)
    n = len(space.states)
    q = csr_matrix((values, (rows, cols)), shape=(n, n), dtype=float)
    q.sum_duplicates()
    row_sums = np.asarray(q.sum(axis=1)).ravel()
    if not np.isfinite(row_sums).all():
        raise ArithmeticError("Sparse generator row sum overflowed")
    q.setdiag(-row_sums)
    q.eliminate_zeros()
    if not np.isfinite(q.data).all():
        raise ArithmeticError("Sparse generator contains nonfinite values")
    return q


def log_action(q, length, log_values):
    """Return log(exp(tQ) exp(log_values)) for vector or matrix right sides."""
    log_values = np.asarray(log_values, dtype=float)
    if not np.isfinite(length) or length < 0:
        raise ValueError("Branch length must be finite and nonnegative")
    if (q.ndim != 2 or q.shape[0] != q.shape[1]
            or log_values.ndim not in (1, 2) or log_values.shape[0] != q.shape[0]):
        raise ValueError("Sparse generator and right-hand side dimensions do not agree")
    if not np.isfinite(q.data).all():
        raise ArithmeticError("Sparse generator contains nonfinite values")
    if np.isnan(log_values).any() or np.isposinf(log_values).any():
        raise ArithmeticError("Log likelihood input contains NaN or positive infinity")
    if length == 0 or q.nnz == 0:
        return log_values.copy()
    values = log_values[:, None] if log_values.ndim == 1 else log_values
    finite = np.isfinite(values)
    if not finite.any():
        return np.full(log_values.shape, -np.inf)
    shift = float(np.max(values[finite]))
    shifts = np.full(values.shape[1], shift)
    valid = np.any(finite, axis=0)
    normalized = np.full(values.shape, -np.inf)
    with np.errstate(over="ignore"):
        normalized[:, valid] = values[:, valid] - shift
    with np.errstate(under="ignore"):
        vector = np.exp(normalized)
    lost_input = finite & (~np.isfinite(normalized) | (vector == 0))
    fallback = valid & np.any(lost_input, axis=0)
    fast_columns = np.flatnonzero(valid & ~fallback)
    a = q * length
    if not np.isfinite(a.data).all():
        raise ArithmeticError("Scaled sparse generator contains nonfinite values")
    trace = float(a.diagonal().sum())
    if not np.isfinite(trace):
        raise ArithmeticError("Scaled sparse generator trace overflowed")
    output = np.full(values.shape, -np.inf)
    if fast_columns.size:
        result = expm_multiply(a, vector[:, fast_columns], traceA=trace)
        if not np.isfinite(result).all() or np.any(result < 0):
            raise ArithmeticError("Sparse CTMC action produced nonfinite or negative values")
        for index, column in enumerate(fast_columns):
            positive = result[:, index] > 0
            reachable = _support_after_edge(
                q, length, np.isfinite(normalized[:, column]))
            if np.any(reachable & ~positive):
                fallback[column] = True
                continue
            output[positive, column] = (
                np.log(result[positive, index]) + shifts[column])
    fallback_columns = np.flatnonzero(fallback)
    if fallback_columns.size:
        output[:, fallback_columns] = log_uniformization_action(
            q, length, values[:, fallback_columns], _support_after_edge)
    return output[:, 0] if log_values.ndim == 1 else output


def _support_after_edge(q, length, child_support):
    """Find parent states with positive CTMC probability of reaching support."""
    if length == 0:
        return child_support.copy()
    reverse = q.copy().tocsr()
    reverse.setdiag(0)
    reverse.eliminate_zeros()
    reverse = reverse.transpose().tocsr()
    supported = child_support.copy()
    stack = list(np.flatnonzero(supported))
    while stack:
        state = stack.pop()
        for parent in reverse.indices[reverse.indptr[state]:reverse.indptr[state + 1]]:
            if not supported[parent]:
                supported[parent] = True
                stack.append(int(parent))
    return supported


def _possible(tree, tips, generators, lengths, root_prior):
    """Boolean pruning distinguishes structural impossibility from numeric collapse."""
    inside = {}
    for node in tree.postorder():
        children = tuple(tree.children.get(node, ()))
        if not children:
            inside[node] = np.asarray(tips[tree.label[node]]) > 0
            continue
        support = np.ones(len(root_prior), dtype=bool)
        for child in children:
            child_support = _support_after_edge(
                generators[child], lengths[child], inside[child])
            support &= child_support
        inside[node] = support
    return bool(np.any((root_prior > 0) & inside[tree.root]))


def likelihood_only(tree, tips, generators, lengths, root_prior):
    """Exact log likelihood using sparse exponential actions and log scaling."""
    root_prior = np.asarray(root_prior, dtype=float)
    n = len(root_prior)
    if (n == 0 or root_prior.shape != (n,) or not np.isfinite(root_prior).all()
            or (root_prior < 0).any() or not np.isclose(root_prior.sum(), 1)):
        raise ValueError("Root prior must be an explicit normalized nonnegative vector")
    if set(tips) != set(tree.leaf_by_label):
        raise ValueError("Exactly one likelihood vector for every tree tip is required")
    for label, values in tips.items():
        values = np.asarray(values, dtype=float)
        if (values.shape != (n,) or not np.isfinite(values).all()
                or (values < 0).any()):
            raise ValueError(f"Invalid observation likelihood: {label}")
    log_tips = {}
    for label, values in tips.items():
        with np.errstate(divide="ignore"):
            log_tips[label] = np.log(np.asarray(values, dtype=float))
    inside = {}
    for node in tree.postorder():
        children = tuple(tree.children.get(node, ()))
        if not children:
            inside[node] = log_tips[tree.label[node]]
            continue
        value = np.zeros(n)
        for child in children:
            value += log_action(generators[child], lengths[child], inside[child])
        if np.isposinf(value).any():
            raise ArithmeticError("Sparse pruning overflowed its log likelihood")
        inside[node] = value
    with np.errstate(divide="ignore"):
        log_root = np.log(root_prior)
    ll = float(logsumexp(log_root + inside[tree.root]))
    if ll == -math.inf and _possible(tree, tips, generators, lengths, root_prior):
        raise ArithmeticError("Sparse pruning collapsed a structurally possible likelihood")
    if math.isnan(ll) or ll == math.inf:
        raise ArithmeticError("Sparse pruning produced a nonfinite log likelihood")
    return ll
