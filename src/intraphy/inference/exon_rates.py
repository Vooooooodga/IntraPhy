"""Pooled rate-scale estimation; no unconstrained per-exon rate estimates.

Relative operation rates and root/origin priors are fixed. A small number of
shared multipliers is estimated across a declared collection of independent
genes. Annotation-discovered catalogues give conditional exploratory fits only.
"""
from __future__ import annotations
from collections import defaultdict
from dataclasses import dataclass, replace
import math
from pathlib import Path
import numpy as np
from scipy.optimize import minimize, minimize_scalar

from ..structure.serialization import read_catalogues
from ..structure.validation import validate_collection
from ..structure.tree_context import canonical_tree
from ..structure.space import StateSpace, enumerate_space
from ..structure.observations import observation_scenarios
from .configuration_model import RateModel, evaluate_model


@dataclass(frozen=True)
class InferenceUnit:
    family: str
    unit: str
    space: StateSpace
    tree: object
    tips: dict[str, np.ndarray]


def load_units(paths, tree, *, max_states=None, observation_view="evidence"):
    units, excluded = [], []
    catalogues = validate_collection(tuple(c for path in paths for c in read_catalogues(path)))
    for c in catalogues:
        space = enumerate_space(c, max_states)
        reason = None
        if c.status != "qualified":
            reason = "unqualified_correspondence"
        elif not space.complete:
            reason = "state_space_incomplete"
        elif any(o.kind == "coexisting" for o in c.observations):
            reason = "coexisting_structures"
        if reason:
            excluded.append({"family": c.family, "unit": c.unit, "reason": reason})
            continue
        _, tips = observation_scenarios(space, tuple(sorted(tree.leaf_by_label)), observation_view)[0]
        if sum(not np.all(v == v[0]) for v in tips.values()) < 2:
            excluded.append({"family": c.family, "unit": c.unit, "reason": "insufficient_observed_taxa"})
            continue
        units.append(InferenceUnit(c.family, c.unit, space, tree, tips))
    if not units:
        raise ValueError("No estimable input units; inspect the configuration diagnostics")
    return tuple(units), excluded


def collection_log_likelihood(units, model, weights=None):
    weights = {u.family: 1 for u in units} if weights is None else weights
    if any(not np.isfinite(v) or v < 0 for v in weights.values()):
        raise ValueError("Gene multiplicities must be finite and nonnegative")
    total = 0.
    for u in units:
        count = weights.get(u.family, 0)
        if count:
            value = evaluate_model(u.space, u.tree, u.tips, model, posterior=False, counts=False)["log_likelihood"]
            if not np.isfinite(value):
                return -math.inf
            total += count*value
    return float(total)


def _curvature(fun, x, h=1e-3):
    x = np.asarray(x, float)
    n, hessian = len(x), np.zeros((len(x), len(x)))
    center = fun(x)
    for i in range(n):
        ei = np.eye(n)[i]*h
        hessian[i, i] = (fun(x+ei)-2*center+fun(x-ei))/h**2
        for j in range(i):
            ej = np.eye(n)[j]*h
            hessian[i, j] = hessian[j, i] = (fun(x+ei+ej)-fun(x+ei-ej)-fun(x-ei+ej)+fun(x-ei-ej))/(4*h**2)
    return hessian


def fit_scale(units, base: RateModel, *, foreground=frozenset(), weights=None,
              bounds=(-9.210340371976184, 9.210340371976184), nested_null=None):
    if not np.isfinite(bounds).all() or bounds[0] >= bounds[1]:
        raise ValueError("Log-scale search bounds must be finite and ordered")
    genes = {u.family for u in units}
    if weights is not None:
        if set(weights)-genes or any(not np.isfinite(v) or v < 0 for v in weights.values()):
            raise ValueError("Gene weights must name existing genes and be finite and nonnegative")
        if not any(weights.values()):
            return {"status": "no_positive_gene_weight", "converged": False, "log_likelihood": None}
    active = genes if weights is None else {g for g in genes if weights.get(g, 0) > 0}
    if not units or len(genes) < 2:
        return {"status": "fewer_than_two_genes", "converged": False, "log_likelihood": None}
    if foreground:
        tree = units[0].tree
        all_children = set(tree.parent)-{tree.root}
        if not foreground < all_children:
            raise ValueError("Foreground must be a nonempty proper subset of non-root branches")
    dimension = 2 if foreground else 1
    def model_at(x):
        return replace(base, scale=math.exp(float(x[0])), foreground=foreground,
                       foreground_multiplier=math.exp(float(x[1])) if foreground else 1.)
    def objective(x):
        return -collection_log_likelihood(units, model_at(x), weights)
    grid = np.linspace(bounds[0], bounds[1], 13)
    starts = []
    if dimension == 1:
        for a, b in zip(grid, grid[1:]):
            fit = minimize_scalar(lambda t: objective([t]), method="bounded", bounds=(a, b), options={"xatol": 1e-6})
            starts.append((float(fit.fun), np.array([fit.x]), bool(fit.success)))
    else:
        initial = [(0, 0), (-2, 0), (2, 0), (0, -2), (0, 2)]
        if nested_null is not None and nested_null.get("scale", 0) > 0:
            initial.insert(0, (math.log(nested_null["scale"]), 0))
        for a, b in initial:
            x = np.clip([a, b], bounds[0], bounds[1])
            fit = minimize(objective, x, method="L-BFGS-B", bounds=[bounds]*2)
            starts.append((float(fit.fun), np.asarray(fit.x), bool(fit.success)))
    finite = [r for r in starts if np.isfinite(r[0]) and r[2]]
    if not finite:
        return {"status": "optimizer_failed_or_impossible_observation", "converged": False, "log_likelihood": None}
    score, theta, _ = min(finite, key=lambda r: r[0])
    if nested_null is not None and -score < nested_null["log_likelihood"]-1e-7:
        return {"status": "nested_optimization_failure", "converged": False,
                "log_likelihood": -score, "null_log_likelihood": nested_null["log_likelihood"]}
    boundary = any(abs(x-bounds[0]) < 1e-3 or abs(x-bounds[1]) < 1e-3 for x in theta)
    hessian = _curvature(objective, theta)
    eigen = np.linalg.eigvalsh(hessian) if np.isfinite(hessian).all() else np.array([math.nan])
    identifiable = bool(np.isfinite(eigen).all() and eigen.min() > 1e-6)
    fitted = model_at(theta)
    status = "at_numerical_boundary" if boundary else "nonidentifiable_curvature" if not identifiable else "estimated_conditional_scale"
    return {"status": status, "converged": True, "valid_for_resampling": not boundary and identifiable,
            "scale": fitted.scale, "foreground_multiplier": fitted.foreground_multiplier,
            "foreground_children": sorted(foreground), "log_likelihood": -score,
            "estimated_parameter_count": dimension, "AIC_conditional": 2*dimension+2*score,
            "curvature_eigenvalues": eigen.tolist(), "log_search_bounds": list(bounds),
            "gene_count": len(genes), "unit_count": len(units),
            "positive_weight_gene_count": len(active),
            "weighted_gene_count": sum(weights.values()) if weights is not None else len(genes),
            "interpretation": "conditional_on_fixed_relative_rates_root_prior_and_candidate_catalogue",
            "formal_asymptotic_P": None}


def gene_bootstrap(units, base, replicates, seed, *, foreground=frozenset()):
    if replicates < 0:
        raise ValueError("Bootstrap replicates cannot be negative")
    genes = tuple(sorted({u.family for u in units}))
    if len(genes) < 2:
        raise ValueError("The original gene collection must contain at least two genes")
    rng = np.random.default_rng(seed)
    draws = []
    for draw in range(replicates):
        sample = rng.choice(genes, size=len(genes), replace=True)
        weights = {gene: int(sum(sample == gene)) for gene in genes}
        result = fit_scale(units, base, foreground=foreground, weights=weights)
        draws.append({"replicate": draw+1, "gene_multiplicities": weights, **result})
    all_valid = bool(draws) and all(d.get("valid_for_resampling") for d in draws)
    interval = np.quantile([d["scale"] for d in draws], [.025, .975]).tolist() if all_valid else None
    return {"replicates": draws, "percentile_interval": interval,
            "resampling_unit": "gene_with_all_its_local_configurations", "seed": seed,
            "interval_status": "conditional_gene_bootstrap" if all_valid else "not_reported_due_to_failed_or_boundary_draws"}
