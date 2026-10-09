"""Truth-only scores for native and factorized exon endpoint histories."""
from __future__ import annotations

import math

import numpy as np

from ..inference.configuration_model import RateModel, evaluate_model
from ..structure.edits import EDIT_KINDS


def native_endpoint_scores(unit, assigned, mu, foreground, multiplier):
    """Score latent visible endpoint pairs using dense full-origin posteriors."""
    model = RateModel({kind: float(mu) for kind in EDIT_KINDS},
        foreground=frozenset(foreground), foreground_multiplier=float(multiplier))
    result = evaluate_model(unit.space, unit.tree, unit.tips, model,
        posterior=True, counts=False, backend="dense")
    if not math.isfinite(float(result.get("log_likelihood", math.nan))):
        raise ArithmeticError("Native history scoring likelihood is nonfinite")
    states = unit.space.states
    visible = [_visible(state) for state in states]
    branches = {(branch["parent"], branch["child"]): branch
                for branch in result.get("branches", ())}
    expected = set(unit.tree.edges())
    if set(branches) != expected:
        raise ArithmeticError("Native history scoring returned incomplete branch endpoints")
    changed = np.fromiter((left.exons != right.exons for left in states for right in states),
                          dtype=bool, count=len(states) ** 2).reshape(
                              (len(states), len(states)))
    scores = []
    for parent, child in unit.tree.edges():
        endpoint = np.asarray(branches[(parent, child)].get("endpoint_probabilities"), dtype=float)
        if endpoint.shape != (len(states), len(states)) or not np.isfinite(endpoint).all():
            raise ArithmeticError(f"Invalid endpoint matrix on {parent}->{child}")
        tolerance = 1e-10
        if (endpoint < -tolerance).any() or (endpoint > 1. + tolerance).any():
            raise ArithmeticError(f"Endpoint probabilities outside [0,1] on {parent}->{child}")
        endpoint = np.clip(endpoint, 0., 1.)
        total = float(endpoint.sum())
        if not math.isfinite(total) or abs(total - 1.) > 1e-7:
            raise ArithmeticError(f"Endpoint probabilities do not sum to one on {parent}->{child}")
        before, after = states[assigned[parent]], states[assigned[child]]
        truth_parent, truth_child = _visible(before), _visible(after)
        parent_indices = [index for index, value in enumerate(visible) if value == truth_parent]
        child_indices = [index for index, value in enumerate(visible) if value == truth_child]
        pair_probability = float(endpoint[np.ix_(parent_indices, child_indices)].sum())
        change_probability = float(endpoint[changed].sum())
        if (not math.isfinite(pair_probability) or not math.isfinite(change_probability)
                or pair_probability < -tolerance or pair_probability > 1. + tolerance
                or change_probability < -tolerance or change_probability > 1. + tolerance):
            raise ArithmeticError(f"Invalid endpoint summary on {parent}->{child}")
        pair_probability = min(1., max(0., pair_probability))
        change_probability = min(1., max(0., change_probability))
        truth_changed = before.exons != after.exons
        scores.append({"parent_node": parent, "child_node": child,
            "probability_exon_geometry_change": change_probability,
            "truth_exon_geometry_changed": truth_changed,
            "brier": (change_probability - float(truth_changed)) ** 2,
            "true_visible_pair_probability": pair_probability,
            "true_visible_pair_log_score": (math.log(pair_probability)
                if pair_probability > 0 else None),
            "zero_probability_true_pair": pair_probability <= 0})
    return scores


def independent_truth_pairs(fit, sampled, assigned, tree):
    """Return truth feature pairs without passing them into the independent fit."""
    pairs = {}
    for parent, child in tree.edges():
        before = sampled.space.states[assigned[parent]]
        after = sampled.space.states[assigned[child]]
        before_spans = {(span.start, span.end) for span in before.exons}
        after_spans = {(span.start, span.end) for span in after.exons}
        before_material = {item.id: int(before.material[index] == 1)
                           for index, item in enumerate(sampled.space.catalogue.material)}
        after_material = {item.id: int(after.material[index] == 1)
                          for index, item in enumerate(sampled.space.catalogue.material)}
        calls = {}
        for feature in fit.get("features", ()):
            if feature["kind"] in {"span", "exon_span"}:
                span = tuple(feature["span"])
                pair = (int(span in before_spans), int(span in after_spans))
            elif feature["kind"] in {"material", "material_presence"}:
                material_id = feature["material_id"]
                pair = (before_material[material_id], after_material[material_id])
            else:
                raise ValueError(f"Unsupported independent feature kind: {feature['kind']}")
            calls[feature["feature_id"]] = pair
        pairs[(parent, child)] = calls
    return pairs


def _mean(values):
    return float(np.mean(values)) if values else None


def panel_summary(genes, arm, expected_branches):
    entries = [gene.get(arm, {}) for gene in genes]
    complete_gene = [entry.get("status") == "scored"
                     and entry.get("branch_count") == expected_branches for entry in entries]
    complete = len(genes) == 4 and all(complete_gene)
    failed = sum(entry.get("status") == "failed" for entry in entries)
    count = sum(complete_gene)
    result = {"status": "scored" if complete else "incomplete", "requested_genes": 4,
        "complete_genes": count, "failed_genes": failed,
        "unresolved_genes": max(0, len(genes) - count - failed),
        "mean_endpoint_brier": None, "mean_true_visible_pair_log_score": None,
        "log_score_status": "incomplete"}
    if complete:
        result["mean_endpoint_brier"] = _mean([entry["mean_endpoint_brier"] for entry in entries])
        zeros = sum(entry.get("zero_probability_true_pairs", 0) for entry in entries)
        result["mean_true_visible_pair_log_score"] = (None if zeros else _mean(
            [entry["mean_true_visible_pair_log_score"] for entry in entries]))
        result["log_score_status"] = "negative_infinity" if zeros else "finite"
        result["zero_probability_true_pairs"] = zeros
    return result


def paired_panel_summary(native, independent):
    complete = native.get("status") == independent.get("status") == "scored"
    out = {"status": "scored" if complete else "incomplete",
        "n_genes": 4 if complete else 0, "log_score_status": "incomplete"}
    if not complete:
        return out
    out.update({"native_mean_endpoint_brier": native["mean_endpoint_brier"],
        "independent_mean_endpoint_brier": independent["mean_endpoint_brier"],
        "independent_minus_native_brier": independent["mean_endpoint_brier"]
            - native["mean_endpoint_brier"]})
    zero_logs = (native.get("log_score_status") == "negative_infinity"
                 or independent.get("log_score_status") == "negative_infinity")
    out["native_mean_log_score"] = (None if zero_logs else
        native["mean_true_visible_pair_log_score"])
    out["independent_mean_log_score"] = (None if zero_logs else
        independent["mean_true_visible_pair_log_score"])
    out["independent_minus_native_log_score"] = (None if zero_logs else
        independent["mean_true_visible_pair_log_score"]
        - native["mean_true_visible_pair_log_score"])
    out["log_score_status"] = "negative_infinity" if zero_logs else "finite"
    return out


def summarize_calibration(records, requested):
    out = {"schema": "intraphy.exon-comparison-calibration-summary/1",
        "requested_replicates": requested, "records_written": len(records),
        "simulation_incomplete_replicates": sum(record.get("simulation_status") != "complete"
                                                  for record in records),
        "records_with_foreground_fit": sum((record.get("foreground_comparison") or {}).get("status")
            not in {"not_run", "failed", None} for record in records),
        "p_values_reported": False, "confidence_intervals_reported": False}
    for key, panel_key in (("native_history", "joint_history_panel"),
                           ("independent", "independent_history_panel")):
        panels = [record.get(panel_key, {}) for record in records]
        complete = [panel for panel in panels if panel.get("status") == "scored"]
        zeros = sum(panel.get("log_score_status") == "negative_infinity" for panel in complete)
        out[key] = {"requested_replicates": requested, "scored_replicates": len(complete),
            "incomplete_replicates": requested - len(complete),
            "mean_endpoint_brier": _mean([panel["mean_endpoint_brier"] for panel in complete]),
            "mean_true_visible_pair_log_score": None if zeros else _mean(
                [panel["mean_true_visible_pair_log_score"] for panel in complete]),
            "negative_infinity_log_score_replicates": zeros}
    paired = [record["paired_panel"] for record in records
              if record.get("paired_panel", {}).get("status") == "scored"]
    zeros = sum(panel.get("log_score_status") == "negative_infinity" for panel in paired)
    out["paired_panel"] = {"requested_replicates": requested, "scored_replicates": len(paired),
        "incomplete_replicates": requested - len(paired),
        "native_mean_endpoint_brier": _mean([panel["native_mean_endpoint_brier"] for panel in paired]),
        "independent_mean_endpoint_brier": _mean(
            [panel["independent_mean_endpoint_brier"] for panel in paired]),
        "independent_minus_native_brier": _mean(
            [panel["independent_minus_native_brier"] for panel in paired]),
        "negative_infinity_log_score_replicates": zeros,
        "log_score_status": "negative_infinity" if zeros else "finite",
        "native_mean_log_score": None if zeros else _mean(
            [panel.get("native_mean_log_score") for panel in paired]),
        "independent_mean_log_score": None if zeros else _mean(
            [panel.get("independent_mean_log_score") for panel in paired]),
        "independent_minus_native_log_score": None if zeros else _mean(
            [panel.get("independent_minus_native_log_score") for panel in paired])}
    out["paired_comparison_scope"] = "same-four-gene-complete-panels-only"
    out["foreground_fit_status_counts"] = {}
    for record in records:
        status = (record.get("foreground_comparison") or {}).get("status") or "unavailable"
        counts = out["foreground_fit_status_counts"]
        counts[status] = counts.get(status, 0) + 1
    return out


def _visible(state):
    return (tuple((span.start, span.end) for span in state.exons),
            tuple(int(value == 1) for value in state.material))
