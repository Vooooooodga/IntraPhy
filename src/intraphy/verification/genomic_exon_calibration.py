"""Scoring summaries for fixed-catalogue genomic exon calibration."""
from __future__ import annotations

import numpy as np


def score_replicate(arm, replicate, units, latent, model):
    from ..inference.configuration_model import evaluate_model
    from ..inference.genomic_exon_branches import branch_change_rows
    from ..inference.genomic_exon_output import state_rows

    output = []
    for unit, (observed, assigned, _) in zip(units, latent):
        posterior = evaluate_model(observed.space, observed.tree, observed.tips,
                                   model, posterior=True, counts=False)
        result = {"family_id": observed.family, "unit_id": observed.unit,
                  "states": state_rows(observed.space), "ctmc": posterior}
        rows = branch_change_rows([{"family_id": observed.family, "units": [result]}],
                                  {}, observed.tree)
        by_branch = {}
        for row in rows:
            by_branch.setdefault((row["parent_node_id"], row["child_node_id"]), []).append(row)
        for (parent, child), tied_rows in by_branch.items():
            before = observed.space.states[assigned[parent]]
            after = observed.space.states[assigned[child]]
            truth_parent = (tuple((e.start, e.end) for e in before.exons),
                            tuple(int(x == 1) for x in before.material))
            truth_child = (tuple((e.start, e.end) for e in after.exons),
                           tuple(int(x == 1) for x in after.material))
            modal_pairs = {(tuple(tuple(span) for span in row["parent_exons"]),
                             tuple(row["parent_dna_presence"]),
                             tuple(tuple(span) for span in row["child_exons"]),
                             tuple(row["child_dna_presence"])) for row in tied_rows}
            hit = (truth_parent[0], truth_parent[1], truth_child[0], truth_child[1]) in modal_pairs
            p = float(tied_rows[0]["probability_exon_structure_change"])
            if not np.isfinite(p) or p < -1e-12 or p > 1 + 1e-12:
                raise FloatingPointError(f"Invalid endpoint-change probability on {parent}->{child}: {p}")
            p = min(1., max(0., p))
            changed = truth_parent[0] != truth_child[0]
            output.append({"replicate": replicate, "arm": arm,
                "unit_id": observed.unit, "parent_node_id": parent,
                "child_node_id": child, "truth_changed": changed,
                "probability_change": p, "modal_tie_count": len(modal_pairs),
                "modal_hit": hit, "modal_credit": 1 / len(modal_pairs) if hit else 0.,
                "truth_parent_exons": truth_parent[0], "truth_child_exons": truth_child[0],
                "modal_pairs": sorted(modal_pairs), "brier": (p-float(changed))**2,
                "accuracy": (p >= .5) == changed,
                "truth_parent_dna_presence": truth_parent[1],
                "truth_child_dna_presence": truth_child[1]})
    expected = len(units) * sum(1 for _ in latent[0][0].tree.edges())
    if len(output) != expected:
        raise ValueError(f"Incomplete branch scores: expected {expected}, received {len(output)}")
    return output


def summarize_arm(arm, replicate_rows, branch_rows, requested, units_per_replicate, branches_per_unit):
    status_column = f"{arm}_status"
    successful_ids = {row["replicate"] for row in replicate_rows if row[status_column] == "scored"}
    scored = [row for row in branch_rows if row["replicate"] in successful_ids]
    replicate_means = [float(np.mean([row["brier"] for row in scored
                                      if row["replicate"] == rep]))
                       for rep in sorted(successful_ids)]
    brier = float(np.mean(replicate_means)) if replicate_means else None
    se = (float(np.std(replicate_means, ddof=1) / np.sqrt(len(replicate_means)))
          if len(replicate_means) > 1 else None)
    bins = []
    for index in range(10):
        lo, hi = index / 10, (index + 1) / 10
        values = [row for row in scored if lo <= row["probability_change"] < hi or
                  index == 9 and row["probability_change"] == 1.]
        bins.append({"bin": index + 1, "lower": lo, "upper": hi, "count": len(values),
            "mean_predicted": float(np.mean([r["probability_change"] for r in values])) if values else None,
            "actual_frequency": float(np.mean([r["truth_changed"] for r in values])) if values else None})
    truths = [row["truth_changed"] for row in scored]
    predictions = [row["probability_change"] for row in scored]
    return {"scope": "fixed_catalogue_conditional_calibration",
        "discovery_pipeline_calibrated": False,
        "conditional_on_scored_replicates": True,
        "requested_replicates": requested, "scored_replicates": len(successful_ids),
        "failed_replicates": sum(r[status_column] == "failed" for r in replicate_rows),
        "nonidentified_replicates": sum(r[status_column] == "nonidentified" for r in replicate_rows),
        "requested_branch_rows": requested * units_per_replicate * branches_per_unit,
        "scored_branch_rows": len(scored),
        "failed_branch_rows": requested * units_per_replicate * branches_per_unit - len(scored),
        "brier": brier, "brier_replicate_mean_monte_carlo_se": se,
        "accuracy_threshold": .5,
        "accuracy_tie_rule": "probability_equal_to_0.5_is_predicted_change",
        "accuracy": (float(np.mean([(p >= .5) == truth for p, truth in zip(predictions, truths)]))
                     if scored else None),
        "true_change_prevalence": float(np.mean(truths)) if truths else None,
        "joint_modal_tie_hit_rate": float(np.mean([r["modal_hit"] for r in scored])) if scored else None,
        "joint_modal_credit_mean": float(np.mean([r["modal_credit"] for r in scored])) if scored else None,
        "reliability_bins": bins}
