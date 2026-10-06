"""Log-domain sparse uniformization for rare CTMC action underflow."""
from __future__ import annotations

import math

import numpy as np
from scipy.special import logsumexp
from scipy.sparse import csr_matrix


def _log_uniformized_matrix(action, rate):
    """Return log(I + action/rate), preserving every positive sparse entry."""
    matrix = action.copy().tocsr()
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    if not np.isfinite(matrix.data).all():
        raise ArithmeticError("Sparse uniformization matrix contains nonfinite rates")
    diagonal = matrix.diagonal()
    if (diagonal > 0).any():
        raise ArithmeticError("Sparse uniformization requires a nonpositive diagonal")
    log_rate = math.log(rate)
    data, indices, indptr = [], [], [0]
    for row in range(matrix.shape[0]):
        start, end = matrix.indptr[row:row + 2]
        for offset in range(start, end):
            column = int(matrix.indices[offset])
            value = float(matrix.data[offset])
            if column == row:
                continue
            if value < 0:
                raise ArithmeticError("Sparse uniformization requires nonnegative off-diagonal rates")
            if value > 0:
                data.append(math.log(value) - log_rate)
                indices.append(column)
        stay = rate + float(diagonal[row])
        if stay < 0:
            raise ArithmeticError("Sparse uniformization rate is below an exit rate")
        if stay > 0:
            data.append(math.log(stay) - log_rate)
            indices.append(row)
        indptr.append(len(data))
    log_matrix = csr_matrix((np.asarray(data), np.asarray(indices, dtype=np.int32),
                             np.asarray(indptr, dtype=np.int64)),
                            shape=matrix.shape)
    row_logs = np.full(matrix.shape[0], -np.inf)
    for row in range(matrix.shape[0]):
        start, end = indptr[row], indptr[row + 1]
        if end > start:
            row_logs[row] = logsumexp(data[start:end])
    log_beta = float(np.max(row_logs, initial=-np.inf))
    if not np.isfinite(log_beta):
        raise ArithmeticError("Sparse uniformization has no positive transition rows")
    return log_matrix, log_beta


def _log_matvec(matrix: csr_matrix, vector: np.ndarray) -> np.ndarray:
    """Apply a nonnegative CSR matrix in the log semiring."""
    output = np.full(matrix.shape[0], -np.inf)
    counts = np.diff(matrix.indptr)
    rows = np.flatnonzero(counts)
    if rows.size:
        starts = matrix.indptr[rows]
        terms = matrix.data + vector[matrix.indices]
        output[rows] = np.logaddexp.reduceat(terms, starts)
    return output


def _one_column(action, length, log_values, support_after_edge):
    reachable = np.asarray(support_after_edge(
        action, length, np.isfinite(log_values)), dtype=bool)
    if reachable.shape != log_values.shape:
        raise ArithmeticError("Sparse action support has an invalid shape")
    if not reachable.any():
        return np.full(log_values.shape, -np.inf)
    if not np.isfinite(log_values[reachable]).any():
        raise ArithmeticError("Sparse action support conflicts with its finite input")

    diagonal = action.diagonal()
    rate = float(np.max(-diagonal, initial=0.))
    if not np.isfinite(rate) or rate <= 0:
        raise ArithmeticError("Sparse action cannot be uniformly scaled")
    log_matrix, log_beta = _log_uniformized_matrix(action, rate)
    log_theta = math.log(rate) + math.log(length)
    theta = rate * length
    if not np.isfinite(theta):
        raise ArithmeticError("Sparse uniformization parameter overflowed")
    log_weight = -theta
    current = log_values.copy()
    total = np.full(log_values.shape, -np.inf)
    log_epsilon = math.log(np.finfo(float).eps)
    log_z = log_theta + log_beta
    maximum_input = float(np.max(log_values[np.isfinite(log_values)]))
    k = 0

    while True:
        total = np.logaddexp(total, log_weight + current)
        if np.isfinite(total[reachable]).all():
            min_output = float(np.min(total[reachable]))
            log_k2 = math.log(k + 2)
            if log_k2 > log_z:
                log_ratio = log_z - log_k2
                ratio = math.exp(log_ratio)
                if ratio < 1.:
                    log_tail = (-theta + (k + 1) * log_z
                                - math.lgamma(k + 2)
                                - math.log1p(-ratio) + maximum_input)
                    if (not np.isnan(log_tail)
                            and log_tail <= log_epsilon + min_output):
                        return total
        current = _log_matvec(log_matrix, current)
        k += 1
        log_weight += log_theta - math.log(k)
        if np.isnan(current).any() or np.isposinf(current).any():
            raise ArithmeticError("Log-domain sparse uniformization became nonfinite")


def log_uniformization_action(action, length, log_values, support_after_edge):
    """Recover affected vector/matrix columns with a relative tail bound."""
    values = np.asarray(log_values, dtype=float)
    matrix = values[:, None] if values.ndim == 1 else values
    if matrix.ndim != 2 or matrix.shape[0] != action.shape[0]:
        raise ValueError("Sparse uniformization dimensions do not agree")
    result = np.column_stack([
        _one_column(action, length, matrix[:, column], support_after_edge)
        for column in range(matrix.shape[1])])
    if np.isnan(result).any() or np.isposinf(result).any():
        raise ArithmeticError("Log-domain sparse uniformization produced invalid output")
    return result[:, 0] if values.ndim == 1 else result
