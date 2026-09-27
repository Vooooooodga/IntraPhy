"""Conditional pooled scalar fitting for explicitly declared repertoire processes."""
from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np
from scipy.optimize import minimize_scalar

from ..structure.validation import validate_collection
from .repertoire_model import evaluate_repertoire_model
from ..structure.repertoire_process import RepertoireProcess


@dataclass(frozen=True)
class RepertoireInferenceUnit:
    process: object
    tree: object
    tips: dict
    root_prior: np.ndarray
    root_prior_provenance: str
    observation_provenance: str

    def __post_init__(self):
        if not isinstance(self.process, RepertoireProcess):
            raise TypeError("process must be a RepertoireProcess")
        _prov("root_prior_provenance", self.root_prior_provenance)
        _prov("observation_provenance", self.observation_provenance)


def _prov(name, value):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty provenance string")
    return value.strip()


def fit_repertoire_scale(units, *, rates, rate_provenance, branch_length_unit,
                         collection_provenance, log_bounds):
    units = tuple(units)
    _prov("rate_provenance", rate_provenance); _prov("branch_length_unit", branch_length_unit)
    _prov("collection_provenance", collection_provenance)
    if len(log_bounds) != 2 or not all(isinstance(x, (int, float)) and not isinstance(x, bool) and np.isfinite(x) for x in log_bounds) or log_bounds[0] >= log_bounds[1]:
        raise ValueError("log_bounds must be two finite ordered values")
    try:
        endpoint_scales = [math.exp(float(x)) for x in log_bounds]
    except (OverflowError, ValueError):
        raise ValueError("log_bounds produce an invalid positive scale")
    if not all(np.isfinite(x) and x > 0 for x in endpoint_scales):
        raise ValueError("log_bounds produce an invalid positive scale")
    if not units:
        return {"status": "fewer_than_two_genes", "converged": False, "log_likelihood": None}
    if any(not isinstance(u, RepertoireInferenceUnit) for u in units):
        raise TypeError("units must contain RepertoireInferenceUnit objects")
    catalogues = tuple(u.process.catalogue for u in units)
    validate_collection(catalogues)
    genes = {c.family for c in catalogues}
    if len(genes) < 2:
        return {"status": "fewer_than_two_genes", "converged": False, "log_likelihood": None,
                "gene_count": len(genes), "unit_count": len(units)}
    if not isinstance(rates, dict) or any(not isinstance(k, str) or not isinstance(v, (int, float)) or isinstance(v, bool) or not np.isfinite(v) or v < 0 for k, v in rates.items()):
        raise ValueError("rates must be finite nonnegative numeric values")
    kinds = set().union(*(set(e.kind for e in u.process.events) for u in units))
    if set(rates) != kinds:
        raise ValueError("rates must match the union of process event kinds, including explicit zeros")
    def score(x):
        scale = math.exp(float(x))
        total = 0.
        impossible = False
        for unit in units:
            local = {k: rates[k] for k in {e.kind for e in unit.process.events}}
            out = evaluate_repertoire_model(unit.process, unit.tree, unit.tips,
                rates=local, root_prior=unit.root_prior,
                root_prior_provenance=unit.root_prior_provenance,
                observation_provenance=unit.observation_provenance,
                rate_provenance=rate_provenance, branch_length_unit=branch_length_unit,
                scale=scale, posterior=False, counts=False)
            value = float(out["log_likelihood"])
            if not np.isfinite(value):
                if np.isnan(value) or value > 0: raise ArithmeticError("invalid unit likelihood")
                impossible = True
            else: total += value
        return -math.inf if impossible else total
    lo, hi = map(float, log_bounds); grid = np.linspace(lo, hi, 13); candidates = []
    metadata = {"log_search_bounds": [lo, hi], "gene_count": len(genes), "unit_count": len(units),
                "rates": dict(rates), "rate_provenance": rate_provenance.strip(),
                "branch_length_unit": branch_length_unit.strip(), "collection_provenance": collection_provenance.strip(),
                "process_provenance": {f"{u.process.catalogue.family}/{u.process.catalogue.unit}": u.process.provenance for u in units},
                "unit_assumptions": [{"family": u.process.catalogue.family, "unit": u.process.catalogue.unit,
                                      "root_prior": np.asarray(u.root_prior, float).tolist(),
                                      "root_prior_provenance": u.root_prior_provenance,
                                      "observation_provenance": u.observation_provenance} for u in units]}
    initial = score((lo + hi) / 2.)
    scale_independent = (all(v == 0 for v in rates.values()) or
                         all(float(u.tree.branch_length(child)) == 0 for u in units for _, child in u.tree.edges()) or
                         all(np.all(np.asarray(v, float) == np.asarray(v, float)[0]) for u in units for v in u.tips.values()))
    if scale_independent:
        if not np.isfinite(initial):
            return {**metadata, "status": "optimizer_failed_or_impossible_observation", "converged": False,
                    "valid_fit": False, "log_likelihood": None, "scale": None, "candidate_scale": None}
        return {**metadata, "status": "nonidentifiable_curvature", "converged": False,
                "valid_fit": False, "log_likelihood": initial, "scale": None,
                "candidate_scale": math.exp((lo + hi) / 2.)}
    failed_intervals = []
    for a, b in zip(grid[:-1], grid[1:]):
        fit = minimize_scalar(lambda x: -score(x), method="bounded", bounds=(a, b), options={"xatol": 1e-8})
        if fit.success and np.isfinite(fit.fun) and np.isfinite(fit.x): candidates.append((float(fit.fun), float(fit.x)))
        else: failed_intervals.append([float(a), float(b)])
    for endpoint in (lo, hi):
        value = score(endpoint)
        if np.isfinite(value): candidates.append((-value, endpoint))
    if not candidates:
        return {**metadata, "status": "optimizer_failed_or_impossible_observation", "converged": False,
                "valid_fit": False, "log_likelihood": None, "scale": None, "candidate_scale": None,
                "failed_intervals": failed_intervals}
    objective, x = min(candidates, key=lambda z: z[0]); h = min(1e-3, (x-lo)/2, (hi-x)/2)
    boundary = abs(x-lo) <= 1e-6 or abs(x-hi) <= 1e-6
    curvature = None if boundary or h <= 0 else (-score(x+h) + 2*score(x) - score(x-h)) / h**2
    status = "optimizer_partial_failure" if failed_intervals else "at_numerical_boundary" if boundary else "nonidentifiable_curvature" if curvature is None or not np.isfinite(curvature) or curvature <= 1e-8 else "estimated_conditional_scale"
    return {**metadata, "status": status, "converged": not failed_intervals, "valid_fit": status == "estimated_conditional_scale", "scale": math.exp(x) if status == "estimated_conditional_scale" else None, "candidate_scale": math.exp(x), "failed_intervals": failed_intervals,
            "log_likelihood": -objective, "log_scale": x, "curvature": None if curvature is None else float(curvature),
            "root_prior_provenance": [u.root_prior_provenance for u in units],
            "observation_provenance": [u.observation_provenance for u in units],
            "unit_assumptions": [{"family": u.process.catalogue.family, "unit": u.process.catalogue.unit,
                                  "root_prior": np.asarray(u.root_prior, float).tolist(),
                                  "root_prior_provenance": u.root_prior_provenance,
                                  "observation_provenance": u.observation_provenance,
                                  "process_provenance": u.process.provenance} for u in units],
            "interpretation": "conditional_on_fixed_relative_rates_root_priors_and_declared_processes"}
