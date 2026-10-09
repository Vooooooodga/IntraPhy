"""Independent binary-character reference for genomic-exon verification."""
from __future__ import annotations

import math

import numpy as np
from scipy.optimize import minimize

from ..inference.ctmc import _pattern_log_likelihood
from ..inference.posterior import _posterior_messages
from ..structure.origins import origin_scenarios
from ..topology import SpeciesTree


_PROBABILITY_TOLERANCE = 1e-10


def _clean_probability_vector(values, context):
    probabilities = np.asarray(values, dtype=float)
    if probabilities.ndim != 1 or probabilities.size == 0:
        raise ValueError(f"{context} must be a nonempty probability vector")
    if np.any(~np.isfinite(probabilities)):
        raise ValueError(f"{context} contains a non-finite probability")
    if (np.any(probabilities < -_PROBABILITY_TOLERANCE) or
            np.any(probabilities > 1.0 + _PROBABILITY_TOLERANCE)):
        raise ValueError(f"{context} contains a probability outside [0, 1]")
    probabilities = np.clip(probabilities, 0.0, 1.0)
    total = math.fsum(float(value) for value in probabilities)
    if not math.isclose(total, 1.0, rel_tol=0.0, abs_tol=1e-8):
        raise ValueError(f"{context} probabilities do not sum to one")
    return probabilities / total


def _features(space):
    spans = sorted({span for state in space.states for span in state.exons})
    features = [{"feature_id": f"exon:{span.start}:{span.end}", "kind": "exon_span",
                 "span": [span.start, span.end]} for span in spans]
    features.extend({"feature_id": f"material:{material.id}", "kind": "material_presence",
                     "material_id": material.id}
                    for material in space.catalogue.material)
    bits = np.array([
        [span in state.exons for span in spans] +
        [state.material[i] == 1 for i in range(len(space.catalogue.material))]
        for state in space.states
    ], dtype=bool)
    return features, bits


def _exact_binary_observations(unit, features, bits):
    space, tree, tips = unit["space"], unit["tree"], unit["tips"]
    if not space.complete or not space.states:
        raise ValueError("independent baseline requires a complete native state space")
    if set(tips) != set(tree.leaf_by_label):
        raise ValueError("tip compatibility vectors must exactly match tree leaf labels")
    by_species = {item.species: item for item in space.catalogue.observations}
    for taxon in tree.leaf_by_label:
        item = by_species.get(taxon)
        if item is not None and item.kind in {"partial", "coexisting"}:
            raise ValueError(f"unsupported {item.kind} observation for independent binary baseline: {taxon}")
    result = {feature["feature_id"]: {} for feature in features}
    for taxon in sorted(tree.leaf_by_label):
        vector = np.asarray(tips[taxon], dtype=float)
        if vector.shape != (len(space.states),) or np.any(~np.isfinite(vector)) or np.any(vector < 0):
            raise ValueError(f"invalid tip compatibility vector: {taxon}")
        if np.any((vector != 0.0) & (vector != 1.0)):
            raise ValueError(f"independent baseline requires exact binary compatibility: {taxon}")
        allowed = vector > 0
        if not allowed.any():
            raise ValueError(f"empty tip compatibility vector: {taxon}")
        if np.all(allowed):
            values = [None] * len(features)
        else:
            values = []
            for column in range(len(features)):
                possible = np.unique(bits[allowed, column])
                if len(possible) != 1:
                    raise ValueError(f"tip observation is partial for visible feature at {taxon}")
                values.append(int(possible[0]))
        for feature, value in zip(features, values):
            result[feature["feature_id"]][taxon] = value
    return result


def _root_feature_marginals(space, tree, features, bits):
    all_states = np.ones(len(space.states), dtype=float)
    unrestricted_tips = {label: all_states for label in tree.leaf_by_label}
    feature_index = {feature["feature_id"]: i for i, feature in enumerate(features)}
    totals = np.zeros(len(features), dtype=float)
    weight_total = 0.0
    scenario_count = 0
    for _origins, root_mask, log_prior in origin_scenarios(
            space, tree, tips=unrestricted_tips, root_weight=1.0):
        legal = np.flatnonzero(root_mask)
        if not len(legal):
            continue
        weight = math.exp(log_prior)
        totals += weight * np.mean(bits[legal], axis=0)
        weight_total += weight
        scenario_count += 1
    if scenario_count == 0 or weight_total <= 0:
        raise ValueError("native root prior has no legal origin scenarios")
    if not math.isclose(weight_total, 1.0, rel_tol=0.0, abs_tol=1e-10):
        raise ValueError(f"full material-origin prior did not resolve unit mass: {weight_total}")
    return ({key: float(totals[index] / weight_total)
             for key, index in feature_index.items()},
            {"method": "full_material_origin_prior_uniform_legal_root_states",
             "root_weight": 1.0, "branch_opportunity_weight": 1.0,
             "origin_scenarios": scenario_count,
             "retained_prior_mass": weight_total,
             "independent_root_marginals": True,
             "joint_root_distribution_matches_native": False})


def _branch_lengths(tree, branch_length_mode):
    if branch_length_mode not in {"supplied", "unit"}:
        raise ValueError("branch_length_mode must be supplied or unit")
    return {child: (tree.branch_length(child) if branch_length_mode == "supplied" else 1.0)
            for _parent, child in tree.edges()}


def _tree_with_lengths(tree, lengths):
    rows = []
    for node, parent in tree.parent.items():
        rows.append({"node_id": node, "parent_id": parent, "label": tree.label[node],
                     "branch_length": "" if not parent else lengths[node]})
    return SpeciesTree(rows)


def _validate_foreground(tree, foreground):
    children = frozenset(foreground)
    edges = {child for _parent, child in tree.edges()}
    if not children or not children <= edges or children == edges:
        raise ValueError("foreground must select canonical edges and leave positive background exposure")
    return children


def _evaluate(space, tree, observations, root_probabilities, rate, rho, foreground, lengths):
    values = []
    for feature_id, pattern in observations.items():
        root_p = root_probabilities[feature_id]
        value = _pattern_log_likelihood(tree, pattern, rate, rate, rho, foreground, root_p)
        if math.isnan(value) or value == math.inf:
            raise ArithmeticError("non-finite independent-character likelihood")
        if value == -math.inf:
            return -math.inf
        values.append(value)
    return float(math.fsum(values))


def _fit_rate(space, tree, observations, root_probabilities, rho, foreground, lengths):
    total_exposure = math.fsum(lengths.values())
    if not math.isfinite(total_exposure) or total_exposure <= 0:
        raise ValueError("positive finite branch exposure is required")
    scale = 1.0 / total_exposure
    starts = [0.0, *(scale * x for x in (0.01, 0.1, 1.0, 10.0, 100.0))]
    probe_values = [_evaluate(space, tree, observations, root_probabilities,
                              rate, rho, foreground, lengths)
                    for rate in starts]
    finite_probe = [value for value in probe_values if math.isfinite(value)]
    if len(finite_probe) == len(probe_values) and max(finite_probe) - min(finite_probe) <= 1e-12:
        return {"status": "no_rate_information", "rate": None,
                "log_likelihood": finite_probe[0],
                "diagnostic": "likelihood is flat over the prespecified rate probes"}
    successful_candidates = []
    failed_candidates = []
    unresolved_failures = []
    for start in starts:
        if start == 0:
            try:
                successful_candidates.append((0.0, _evaluate(
                    space, tree, observations, root_probabilities, 0.0, rho, foreground, lengths),
                    "explicit zero-rate boundary"))
            except (ArithmeticError, ValueError) as exc:
                unresolved_failures.append(str(exc))
            continue
        try:
            result = minimize(
                lambda x: -_evaluate(space, tree, observations, root_probabilities,
                                      float(x[0]), rho, foreground, lengths),
                x0=np.array([start]), method="L-BFGS-B", bounds=[(0.0, None)],
                options={"maxiter": 1000, "ftol": 1e-12, "gtol": 1e-8},
            )
            rate = max(0.0, float(result.x[0]))
            ll = _evaluate(space, tree, observations, root_probabilities,
                           rate, rho, foreground, lengths)
            if math.isfinite(ll):
                candidate = (rate, ll, str(result.message))
                (successful_candidates if result.success else failed_candidates).append(candidate)
            elif ll == -math.inf:
                (successful_candidates if result.success else failed_candidates).append(
                    (rate, ll, str(result.message)))
            else:
                unresolved_failures.append("optimizer produced non-finite likelihood")
        except (ArithmeticError, ValueError, FloatingPointError) as exc:
            unresolved_failures.append(str(exc))
    finite = [(rate, ll) for rate, ll, _message in successful_candidates if math.isfinite(ll)]
    if not finite:
        if failed_candidates or unresolved_failures:
            details = [message for _rate, _ll, message in failed_candidates if message]
            details.extend(unresolved_failures)
            return {"status": "optimization_unresolved", "rate": None,
                    "log_likelihood": None,
                    "diagnostic": "; ".join(details) or "no successful finite rate candidate"}
        if any(ll == -math.inf for _rate, ll, _message in successful_candidates):
            return {"status": "impossible_observation_under_model", "rate": None,
                    "log_likelihood": None, "diagnostic": "all evaluated rates have zero likelihood"}
        return {"status": "optimizer_failed", "rate": None, "log_likelihood": None,
                "diagnostic": "no successful finite rate candidate"}
    rate, ll = max(finite, key=lambda pair: (pair[1], -pair[0]))
    failed_above = [(candidate_rate, candidate_ll, message)
                    for candidate_rate, candidate_ll, message in failed_candidates
                    if math.isfinite(candidate_ll) and candidate_ll > ll + 1e-8]
    if unresolved_failures or failed_above:
        details = [message for _candidate_rate, _candidate_ll, message in failed_above if message]
        details.extend(unresolved_failures)
        return {"status": "optimization_unresolved", "rate": rate, "log_likelihood": ll,
                "diagnostic": "; ".join(details) or "failed candidate exceeds best successful likelihood"}
    lower_failed = [message for _candidate_rate, candidate_ll, message in failed_candidates
                    if candidate_ll <= ll + 1e-8 and message]
    if rate == 0.0:
        return {"status": "zero_boundary", "rate": rate, "log_likelihood": ll,
                "diagnostic": "maximum likelihood occurs at zero rate" +
                ("; lower failed candidates: " + "; ".join(lower_failed) if lower_failed else "")}
    probe = max(rate * 2.0, scale * 200.0)
    probe_ll = _evaluate(space, tree, observations, root_probabilities,
                         probe, rho, foreground, lengths)
    if probe_ll >= ll - 1e-8:
        return {"status": "upper_tail_unresolved", "rate": rate, "log_likelihood": ll,
                "diagnostic": "a larger-rate probe is not lower in likelihood"}
    return {"status": "estimated", "rate": rate, "log_likelihood": ll,
            "diagnostic": "finite interior maximum found" +
            ("; lower failed candidates: " + "; ".join(lower_failed) if lower_failed else "")}


def _illegal_node_mass(space, tree, features, node_marginals):
    feature_bits = []
    for state in space.states:
        state_exons = {(e.start, e.end) for e in state.exons}
        feature_bits.append(tuple(
            (tuple(feature["span"]) in state_exons if feature["kind"] == "exon_span"
             else state.material[next(i for i, m in enumerate(space.catalogue.material)
                                      if m.id == feature["material_id"])] == 1)
            for feature in features))
    state_signatures = set(feature_bits)
    result = {}
    for node in tree.preorder():
        legal_mass = 0.0
        clean_marginals = {}
        for feature in features:
            key = feature["feature_id"]
            clean_marginals[key] = _clean_probability_vector(
                node_marginals[key][node], f"node {node} feature {key}")
        for signature in state_signatures:
            probability = 1.0
            for feature, bit in zip(features, signature):
                p = clean_marginals[feature["feature_id"]][1]
                probability *= p if bit else 1.0 - p
            legal_mass += probability
        illegal = 1.0 - legal_mass
        if illegal < -1e-10 or illegal > 1.0 + 1e-10:
            raise ArithmeticError("independent node illegal mass is outside [0, 1]")
        result[node] = min(1.0, max(0.0, illegal))
    return result


def fit_independent_exon_characters(unit, foreground, *, foreground_multiplier=1.0,
                                    branch_length_mode="supplied"):
    """Fit one independent binary rate to exact visible features in one unit.

    No generating-rate or truth argument is accepted. The rate is per binary
    character and branch-length unit; it is distinct from native per-edit mu.
    """
    rho = float(foreground_multiplier)
    if not math.isfinite(rho) or rho < 0:
        raise ValueError("foreground_multiplier must be finite and nonnegative")
    space, tree = unit["space"], unit["tree"]
    foreground = _validate_foreground(tree, foreground)
    lengths = _branch_lengths(tree, branch_length_mode)
    model_tree = _tree_with_lengths(tree, lengths)
    if not math.fsum(lengths[c] for c in foreground) > 0 or not math.fsum(
            length for child, length in lengths.items() if child not in foreground) > 0:
        raise ValueError("foreground and background must both have positive branch exposure")
    features, bits = _features(space)
    observations = _exact_binary_observations(unit, features, bits)
    root_probabilities, root_prior = _root_feature_marginals(space, model_tree, features, bits)
    fit = _fit_rate(space, model_tree, observations, root_probabilities, rho, foreground, lengths)
    output = {
        "status": fit["status"], "rate": fit["rate"], "log_likelihood": fit["log_likelihood"],
        "diagnostic": fit["diagnostic"], "rate_units": "per_binary_character_per_branch_length_unit",
        "foreground_multiplier": rho, "branch_length_mode": branch_length_mode,
        "features": features, "tip_observations": observations,
        "root_feature_probabilities": root_probabilities, "root_prior": root_prior,
        "native_rate_comparable_to_binary_rate": False,
        "foreground_children": sorted(foreground), "node_illegal_mass": None,
        "branch_posteriors": None,
    }
    if fit["rate"] is not None and fit["status"] in {"estimated", "zero_boundary"}:
        node_marginals, edge_posteriors = {}, {}
        try:
            for feature in features:
                key = feature["feature_id"]
                node, edge = _posterior_messages(model_tree, observations[key], fit["rate"], fit["rate"],
                    rho, foreground, root_probabilities[key])
                node_marginals[key] = node
                edge_posteriors[key] = edge
        except (ValueError, FloatingPointError) as exc:
            output["status"] = "posterior_failed"
            output["diagnostic"] = str(exc)
            return output
        rows = []
        for parent, child in tree.edges():
            feature_rows = []
            for feature in features:
                joint = edge_posteriors[feature["feature_id"]][(parent, child)]
                feature_rows.append({"feature_id": feature["feature_id"], "kind": feature["kind"],
                    "p00": float(joint[0, 0]), "p01": float(joint[0, 1]),
                    "p10": float(joint[1, 0]), "p11": float(joint[1, 1])})
            rows.append({"parent_node": parent, "child_node": child, "features": feature_rows})
        output["branch_posteriors"] = rows
        output["node_illegal_mass"] = _illegal_node_mass(space, tree, features, node_marginals)
    return output


def score_independent_endpoint_posteriors(fit, truth_feature_pairs):
    """Score complete true visible endpoint pairs without refitting the model.

    ``truth_feature_pairs`` maps ``(parent, child)`` to feature IDs and their
    ``(parent_bit, child_bit)`` truth pair. Truth is consumed only by scoring.
    """
    if (fit.get("status") not in {"estimated", "zero_boundary"} or
            fit.get("branch_posteriors") is None):
        return {"status": "no_call", "edges": []}
    scored = []
    for row in fit["branch_posteriors"]:
        edge = (row["parent_node"], row["child_node"])
        truth = truth_feature_pairs.get(edge)
        if truth is None:
            raise ValueError(f"missing truth feature pairs for edge {edge}")
        posterior = {}
        for item in row["features"]:
            probabilities = _clean_probability_vector(
                [item[name] for name in ("p00", "p01", "p10", "p11")],
                f"edge {edge} feature {item.get('feature_id')}")
            posterior[item["feature_id"]] = {
                **item, **{name: float(probabilities[index])
                           for index, name in enumerate(("p00", "p01", "p10", "p11"))}}
        if set(truth) != set(posterior):
            raise ValueError(f"truth feature roster differs from fitted roster for edge {edge}")
        log_probability = 0.0
        any_change_same = 1.0
        for feature in fit["features"]:
            key = feature["feature_id"]
            start, end = truth[key]
            if start not in (0, 1) or end not in (0, 1):
                raise ValueError("truth feature endpoint bits must be 0 or 1")
            p = posterior[key][f"p{start}{end}"]
            log_probability = -math.inf if p <= 0 else log_probability + math.log(p)
            if feature["kind"] == "exon_span":
                item = posterior[key]
                any_change_same *= item["p00"] + item["p11"]
        scored.append({"parent_node": edge[0], "child_node": edge[1],
                       "true_full_endpoint_log_probability": log_probability,
                       "p_any_exon_geometry_endpoint_change": 1.0 - any_change_same})
    return {"status": "scored", "edges": scored}


__all__ = ["fit_independent_exon_characters", "score_independent_endpoint_posteriors"]
