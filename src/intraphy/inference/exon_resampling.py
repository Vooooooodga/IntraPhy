"""Conditional finite-catalogue simulation and foreground model comparison.

This does NOT calibrate the FASTA/GFF discovery pipeline. A Monte Carlo P value
is withheld unless the collection is explicitly independently specified and
all fitted bootstrap draws are valid. Failed replicates are never discarded.
"""
from __future__ import annotations
from dataclasses import replace
from itertools import product
import math
import numpy as np
from .configuration_model import generator, RateModel
from .configuration_ctmc import transition_matrix
from .exon_rates import fit_scale, InferenceUnit
from ..structure.tree_context import canonical_tree


def sample_unit_history(unit: InferenceUnit, model: RateModel, rng):
    """Return simulated observations together with latent assignments and origins.

    The draw order is shared with :func:`simulate_unit`; callers that only need
    tip observations should keep using that stable API.
    """
    context = canonical_tree(unit.tree, model.foreground)
    space, tree = unit.space, context.tree
    if any(n != tree.root and len(tree.children.get(n, ())) == 1 for n in tree.parent):
        raise ValueError("Mixed rate regimes inside a subdivided branch are unsupported")
    model = replace(model, foreground=context.foreground)
    nodes = tuple(sorted(tree.parent))
    opportunity = np.array([model.origin_root_weight if n == tree.root else 1. for n in nodes])
    opportunity /= opportunity.sum()
    origins = {m.id: str(rng.choice(nodes, p=opportunity)) for m in space.catalogue.material}
    required = tuple(1 if origins[m.id] == tree.root else 0 for m in space.catalogue.material)
    valid_root = [i for i, state in enumerate(space.states) if state.material == required]
    assigned = {tree.root: int(rng.choice(valid_root))}
    for parent in tree.preorder():
        for child in tree.children.get(parent, ()):
            q, _ = generator(space, model, origins, child)
            p = transition_matrix(q, tree.branch_length(child))
            # Generator is finite/closed. Adjustment is roundoff only for RNG.
            values = p[assigned[parent]].copy()
            values[np.argmax(values)] += 1-values.sum()
            assigned[child] = int(rng.choice(len(space.states), p=values))
    tips = {}
    for species, node in tree.leaf_by_label.items():
        if np.all(unit.tips[species] == unit.tips[species][0]):
            tips[species] = np.ones(len(space.states))
            continue
        simulated = space.states[assigned[node]]
        # Observations cannot distinguish unintroduced material from material
        # that was deleted; do not leak latent origin labels into the data.
        tips[species] = np.array([float(s.exons == simulated.exons and
            all((a == 1) == (b == 1) for a, b in zip(s.material, simulated.material))) for s in space.states])
    return replace(unit, tree=tree, tips=tips), assigned, origins


def simulate_unit(unit: InferenceUnit, model: RateModel, rng) -> InferenceUnit:
    """Simulate tip observations while preserving the historical return API."""
    simulated, _, _ = sample_unit_history(unit, model, rng)
    return simulated


def foreground_bootstrap(units, base: RateModel, foreground, replicates: int, seed: int):
    foreground = frozenset(foreground)
    if not foreground:
        raise ValueError("At least one foreground child node is required")
    if replicates < 0:
        raise ValueError("Replicate count cannot be negative")
    # Null and alternative must use the same opportunity catalogue. A rate
    # boundary inserted inside a unary edge would otherwise change its prior.
    for unit in units:
        if set(canonical_tree(unit.tree).tree.parent) != set(canonical_tree(unit.tree, foreground).tree.parent):
            raise ValueError("Foreground on a subdivided edge is not supported for rate comparison; "
                             "declare the contrast on canonical biological branches")
    null = fit_scale(units, base)
    alternative = fit_scale(units, base, foreground=foreground, nested_null=null)
    result = {"null": null, "alternative": alternative, "parametric_replicates": [],
              "conditional_monte_carlo_P": None, "discovery_pipeline_calibrated": False,
              "test_scope": "one_prespecified_foreground_rate_contrast_over_the_entire_gene_collection",
              "multiple_testing": "one_test; no hidden_per_gene_tests", "seed": seed}
    if not null.get("valid_for_resampling") or not alternative.get("valid_for_resampling"):
        result["status"] = "invalid_observed_fit"
        return result
    observed = nested_statistic(null, alternative)
    if observed is None:
        result["status"] = "nested_optimization_failure"
        return result
    result["likelihood_ratio"] = observed
    if any(u.space.catalogue.discovery != "independent_catalogue" for u in units):
        result["status"] = "annotation_discovered_catalogue_not_calibrated"
        return result
    # Partial/candidate coarsening needs its own observation simulator. The
    # implemented calibration accepts exact visible configurations or missing tips.
    for u in units:
        for values in u.tips.values():
            allowed = [s for s, v in zip(u.space.states, values) if v > 0]
            observable = {(s.exons, tuple(x == 1 for x in s.material)) for s in allowed}
            if not np.all(values == values[0]) and len(observable) != 1:
                result["status"] = "partial_observation_process_not_calibrated"
                return result
    rng = np.random.default_rng(seed)
    null_model = replace(base, scale=null["scale"], foreground=frozenset(), foreground_multiplier=1.)
    for i in range(replicates):
        simulation = tuple(simulate_unit(u, null_model, rng) for u in units)
        nfit = fit_scale(simulation, base)
        afit = fit_scale(simulation, base, foreground=foreground, nested_null=nfit)
        valid = bool(nfit.get("valid_for_resampling") and afit.get("valid_for_resampling"))
        statistic = nested_statistic(nfit, afit) if valid else None
        valid = valid and statistic is not None
        result["parametric_replicates"].append({"replicate": i+1, "valid": valid,
            "likelihood_ratio": statistic, "null_status": nfit["status"], "alternative_status": afit["status"]})
    draws = result["parametric_replicates"]
    if draws and all(d["valid"] for d in draws):
        result["conditional_monte_carlo_P"] = (1+sum(d["likelihood_ratio"] >= observed for d in draws))/(replicates+1)
        result["monte_carlo_resolution"] = 1/(replicates+1)
        result["status"] = "conditional_finite_catalogue_test_only"
    else:
        result["status"] = "bootstrap_not_run_or_contains_failed_draws"
    return result


def nested_statistic(null, alternative, tolerance=1e-7):
    """Clamp numerical roundoff only; a genuinely worse nested fit is a failure."""
    a, n = alternative.get("log_likelihood"), null.get("log_likelihood")
    if a is None or n is None or not np.isfinite([a, n]).all() or a < n-tolerance:
        return None
    return max(0., 2*(a-n))
