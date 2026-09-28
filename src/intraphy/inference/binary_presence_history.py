"""Domain-labelled posterior exports for binary presence sites."""
from __future__ import annotations

from intraphy.inference.posterior import _expected_transition_count, _posterior_messages


def _site_status(site):
    if not site["eligible"]:
        return "excluded"
    if not site["known_states"]:
        return "all_unknown"
    if site["known_states"] == {1}:
        return "no_absence_observed"
    if site["known_states"] == {0}:
        return "no_presence_observed"
    return "observed_contrast"


def _history_for_family(family, eligible, tree, fit, status, root_frequency, root_presence,
                        expected_counts, domain):
    record = {"family_id": family, "status": status, "model": domain["model"], "sites": []}
    if fit is None or not fit.get("optimizer_success") or not fit.get("rates"):
        record["posterior_available"] = False
        record["posterior_unavailable_reason"] = "no_supported_rate_fit_or_fixed_rates"
        return record
    gain, loss = fit["rates"]["gain"], fit["rates"]["loss"]
    rho = None if root_frequency == "stationary" else float(root_presence)
    for site_id, site in sorted(eligible.items()):
        record_site = {"site_id": site_id, "status": _site_status(site)}
        if not site["known_states"] or 1 not in site["known_states"]:
            record_site["posterior_available"] = False
            record_site["posterior_unavailable_reason"] = "site_not_selected_by_observed_at_least_one_ascertainment"
            record["sites"].append(record_site)
            continue
        node, edge = _posterior_messages(tree, site["observations"], gain, loss, 1.0, frozenset(), rho)
        record_site.update({"posterior_available": True, "nodes": [
            {"node_id": node_id, "label": tree.label[node_id], domain["node_p0"]: float(prob[0]),
             domain["node_p1"]: float(prob[1])}
            for node_id, prob in sorted(node.items())
        ]})
        edges = []
        for (parent, child), joint in sorted(edge.items()):
            item = {"parent_node": parent, "child_node": child,
                    domain["edge_p01"]: float(joint[0, 1]),
                    domain["edge_p10"]: float(joint[1, 0]),
                    "p_any_endpoint_change": float(joint[0, 1] + joint[1, 0])}
            if expected_counts:
                item[domain["expected_01"]] = _expected_transition_count(
                    joint, gain, loss, tree.branch_length(child), 1.0, 0, 1
                )
                item[domain["expected_10"]] = _expected_transition_count(
                    joint, gain, loss, tree.branch_length(child), 1.0, 1, 0
                )
            edges.append(item)
        record_site["edges"] = edges
        record["sites"].append(record_site)
    record["posterior_available"] = any(site.get("posterior_available") for site in record["sites"])
    if not record["posterior_available"]:
        record["posterior_unavailable_reason"] = "no_selected_eligible_sites"
    return record
