"""Joint finite-repertoire CTMC evaluation.

This module evaluates an explicitly declared repertoire process.  It does not
invent event rates, detection probabilities, root priors, or branch units.
"""
from __future__ import annotations

import numpy as np

from .configuration_ctmc import (
    likelihood, transition_matrix, expected_count, probability_any_change,
)


def _provenance(name: str, value: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty provenance string")
    return value.strip()


def evaluate_repertoire_model(process, tree, tips, *, rates,
                              root_prior, root_prior_provenance,
                              observation_provenance, rate_provenance,
                              branch_length_unit,
                              scale=1., posterior=True, counts=True):
    """Evaluate one joint repertoire process on a supplied tree.

    ``process.generator`` owns the finite declared event table and returns
    ``(Q, marks)``.  ``tips`` contains one emission vector per leaf label.
    Branch lengths are always read from the supplied tree; no unit fallback or
    tree canonicalization is performed here.
    """
    _provenance("root_prior_provenance", root_prior_provenance)
    _provenance("observation_provenance", observation_provenance)
    _provenance("rate_provenance", rate_provenance)
    _provenance("branch_length_unit", branch_length_unit)
    from ..structure.repertoire_process import RepertoireProcess
    if not isinstance(process, RepertoireProcess):
        raise TypeError("process must be a validated RepertoireProcess")
    if not isinstance(posterior, bool) or not isinstance(counts, bool):
        raise TypeError("posterior and counts must be Boolean")
    states = tuple(process.states)
    n = len(states)
    if n == 0:
        raise ValueError("Repertoire process must declare at least one state")
    prior = np.asarray(root_prior, dtype=float)
    if prior.shape != (n,) or not np.isfinite(prior).all() or (prior < 0).any() or not np.isclose(prior.sum(), 1., rtol=0., atol=1e-12):
        raise ValueError("root_prior must be a finite normalized vector matching process.states")
    if isinstance(scale, (bool, np.bool_)) or not np.isfinite(scale) or scale <= 0:
        raise ValueError("scale must be finite and positive")
    if not isinstance(tips, dict) or set(tips) != set(tree.leaf_by_label):
        raise ValueError("tips must provide exactly one emission vector per tree leaf label")
    emissions = {}
    for label, vector in tips.items():
        values = np.asarray(vector, dtype=float)
        if values.shape != (n,) or not np.isfinite(values).all() or (values < 0).any() or (values > 1).any():
            raise ValueError(f"Invalid tip emission vector: {label}")
        emissions[label] = values
    lengths = {}
    for _, child in tree.edges():
        length = tree.branch_length(child)
        if not np.isfinite(length) or length < 0:
            raise ValueError(f"Invalid branch length for {child}")
        lengths[child] = float(length)
    declared_rates = dict(rates)
    q, marks = process.generator(declared_rates, scale=scale, include_marks=counts and posterior)
    q = np.asarray(q, dtype=float)
    if q.shape != (n, n):
        raise ValueError("Process generator dimension does not match process.states")
    transitions = {child: transition_matrix(q, length) for child, length in lengths.items()}
    result = likelihood(tree, emissions, transitions, prior, posterior=posterior)
    output = {
        "log_likelihood": result.log_likelihood,
        "nodes": result.nodes if posterior else {},
        "branches": {},
        "state_count": n,
        "states": states,
        "catalogue": {"family": process.catalogue.family, "unit": process.catalogue.unit},
        "event_table": tuple(getattr(process, "events", ())),
        "finite_declared_event_table_only": True,
        "root_prior_provenance": root_prior_provenance.strip(),
        "observation_provenance": observation_provenance.strip(),
        "branch_length_unit": branch_length_unit.strip(),
        "rate_provenance": rate_provenance.strip(),
        "process_provenance": process.provenance,
        "rates": dict(declared_rates),
        "scale": float(scale),
        "root_prior": prior.copy(),
    }
    if not posterior:
        return output
    for edge, endpoint in result.endpoints.items():
        _, child = edge
        item = {"endpoint_matrix": endpoint,
                "probability_at_least_one_edit": probability_any_change(q, transitions[child], endpoint, lengths[child]),
                "probability_different_endpoints": 0.}
        different = 1. - float(np.trace(endpoint))
        if not np.isfinite(different) or not -1e-8 <= different <= 1. + 1e-8:
            raise ArithmeticError("Invalid probability of different endpoints")
        item["probability_different_endpoints"] = float(np.clip(different, 0., 1.))
        if counts:
            item["expected_edits"] = expected_count(q, transitions[child], endpoint, lengths[child])
            for mark, matrix in (marks or {}).items():
                item[f"expected_{mark}"] = expected_count(q, transitions[child], endpoint, lengths[child], np.asarray(matrix, dtype=float))
        output["branches"][edge] = item
    return output
