"""Conditional phylogenetic inference of homologous DNA presence states."""
from __future__ import annotations

import math
from numbers import Real

from intraphy.inference.ctmc import _compress_patterns
from intraphy.inference.dna_fit import fit_dna_family
from intraphy.inference.dna_output import write_dna_outputs
from intraphy.inference.posterior import _expected_transition_count, _posterior_messages


def analyze_dna_presence(
    rows, tree, output_dir, *, root_frequency="stationary", root_presence=0.5,
    fixed_rates=None, branch_length_unit="supplied_tree_units", expected_counts=False,
    workers=1,
):
    """Fit family-shared DNA presence gain/loss rates and export conditional histories.

    Each row supplies ``family_id``, ``site_id``, ``species``, ``state`` (0, 1,
    or ``unknown``), ``eligible`` (Boolean), ``reason``, and ``evidence``.
    Sites are discovered conditional on at least one observed present tip.
    """
    rows = tuple(rows)
    _validate_options(root_frequency, root_presence, fixed_rates, branch_length_unit, expected_counts, workers)
    _validate_tree(tree)
    families = _group_rows(rows, tree)
    fit_rows = []
    history_rows = []
    fixed = dict(fixed_rates or {})
    for family, sites in sorted(families.items()):
        eligible = {site_id: site for site_id, site in sites.items() if site["eligible"]}
        excluded = [site for site_id, site in sorted(sites.items()) if not site["eligible"]]
        informative_sites = {
            site_id: site["observations"] for site_id, site in eligible.items()
            if 0 in site["known_states"] and 1 in site["known_states"]
        }
        selected_sites = {
            site_id: site for site_id, site in eligible.items() if 1 in site["known_states"]
        }
        any_contrast = bool(informative_sites)
        known_states = {state for site in eligible.values() for state in site["known_states"]}
        reasons = []
        site_records = []
        for site_id, site in sorted(sites.items()):
            site_status = _site_status(site)
            site_records.append({
                "site_id": site_id, "eligible": site["eligible"], "status": site_status,
                "known_tip_count": len(site["known_states_by_tip"]),
                "observed_states": sorted(site["known_states"]),
                "reason": site["reason"], "evidence": site["evidence"],
            })
        fixed_complete = set(fixed) == {"gain", "loss"}
        if not eligible:
            status = "all_excluded"
            reasons.append("no_eligible_sites")
            fit = None
        elif all(not site["known_states"] for site in eligible.values()):
            status = "all_unknown"
            reasons.append("eligible_sites_have_no_known_tip_states")
            fit = None
        elif not known_states or known_states == {0}:
            status = "no_observed_presence"
            reasons.append("no_eligible_site_has_an_observed_present_tip")
            fit = None
        elif not any_contrast and not fixed_complete:
            status = "no_observed_contrast"
            reasons.append("no_site_has_both_observed_absence_and_presence")
            fit = None
        else:
            # Retain eligible discovered sites, including constant-present
            # patterns, under each site's fixed observed-tip mask.
            fit = fit_dna_family(
                tree, {site_id: site["observations"] for site_id, site in selected_sites.items()},
                root_frequency=root_frequency, root_presence=root_presence,
                fixed_rates=fixed, workers=workers,
            )
            if not fit.get("optimizer_success"):
                status = "optimizer_unresolved"
                reasons.append("rate_fit_did_not_converge_or_was_non_finite")
            else:
                status = "estimated" if not fixed_complete else "fixed_rates"
                if fit.get("boundary_rates"):
                    reasons.append("one_or_more_rates_at_zero_boundary")
                if fit.get("locally_flat"):
                    reasons.append("local_curvature_is_flat")
        fit_row = {
            "family_id": family, "status": status, "reasons": reasons,
            "model": "dna-presence-ctmc", "root_frequency": root_frequency,
            "root_presence": root_presence if root_frequency == "fixed" else None,
            "branch_length_unit": branch_length_unit,
            "ascertainment": "observed-at-least-one;conditioned on each site's fixed observation mask",
            "n_sites": len(sites), "n_eligible_sites": len(eligible),
            "n_excluded_sites": len(excluded), "n_all_unknown_sites": sum(not s["known_states"] for s in eligible.values()),
            "n_variable_sites": len(informative_sites),
            "n_constant_present_sites": sum(s["known_states"] == {1} for s in eligible.values()),
            "n_single_known_present_sites": sum(
                s["known_states"] == {1} and len(s["known_states_by_tip"]) == 1
                for s in eligible.values()
            ),
            "n_constant_absent_sites": sum(s["known_states"] == {0} for s in eligible.values()),
            "n_compressed_patterns": len(_compress_patterns(
                [s["observations"] for s in selected_sites.values()], tuple(sorted(tree.leaf_by_label))
            )),
            "workers_requested": workers,
            "independent_site_assumption": "site likelihoods are summed; linked sites make this a composite likelihood",
            "gain_rate": fit["rates"].get("gain") if fit else None,
            "loss_rate": fit["rates"].get("loss") if fit else None,
            "log_likelihood": _finite_or_none(fit.get("log_likelihood")) if fit else None,
            "optimizer_success": fit.get("optimizer_success", False) if fit else False,
            "rate_information_status": _rate_information_status(fit, status),
            "optimizer_status": fit.get("optimizer_status") if fit else None,
            "optimizer_message": fit.get("optimizer_message") if fit else None,
            "iterations": fit.get("iterations") if fit else None,
            "boundary_rates": fit.get("boundary_rates", []) if fit else [],
            "locally_flat": fit.get("locally_flat") if fit else None,
            "curvature_eigenvalues": fit.get("curvature_eigenvalues", []) if fit else [],
            "selected_start": fit.get("selected_start") if fit else None,
            "starts": fit.get("starts", []) if fit else [],
            "excluded_sites": [{"site_id": site_id, "reason": site["reason"], "evidence": site["evidence"]}
                               for site_id, site in sorted(sites.items()) if not site["eligible"]],
            "sites": site_records,
        }
        fit_rows.append(fit_row)
        history_rows.append(_history_for_family(
            family, eligible, tree, fit, status, root_frequency, root_presence,
            expected_counts,
        ))
    write_dna_outputs(output_dir, tree, fit_rows, history_rows)
    return {"model": "dna-presence-ctmc", "families": fit_rows, "history": history_rows}


def _validate_options(root_frequency, root_presence, fixed_rates, branch_length_unit, expected_counts, workers):
    if root_frequency not in ("stationary", "fixed"):
        raise ValueError("root_frequency must be 'stationary' or 'fixed'")
    if isinstance(root_presence, bool) or not isinstance(root_presence, Real) or not math.isfinite(float(root_presence)) or not 0 <= float(root_presence) <= 1:
        raise ValueError("root_presence must be finite and between 0 and 1")
    if branch_length_unit != "supplied_tree_units":
        raise ValueError("branch_length_unit must be 'supplied_tree_units'; units are not inferred")
    if type(expected_counts) is not bool:
        raise ValueError("expected_counts must be Boolean")
    if type(workers) is not int or workers < 1:
        raise ValueError("workers must be a positive integer")
    fixed = dict(fixed_rates or {})
    if set(fixed) - {"gain", "loss"}:
        raise ValueError("fixed_rates may contain only 'gain' and 'loss'")
    for rate, value in fixed.items():
        if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(float(value)) or float(value) < 0:
            raise ValueError(f"fixed rate {rate!r} must be finite and nonnegative")
    if fixed and set(fixed) != {"gain", "loss"}:
        raise ValueError("fixed_rates must specify both gain and loss together")


def _validate_tree(tree):
    labels = set(tree.leaf_by_label)
    if not labels:
        raise ValueError("tree must have at least one tip")
    for _parent, child in tree.edges():
        length = tree.branch_length(child)
        if not math.isfinite(length) or length < 0:
            raise ValueError("all supplied non-root branch lengths must be finite and nonnegative")


def _group_rows(rows, tree):
    labels = set(tree.leaf_by_label)
    grouped = {}
    seen = set()
    for index, row in enumerate(rows):
        family = str(row.get("family_id", "")).strip()
        site = str(row.get("site_id", "")).strip()
        species = str(row.get("species", "")).strip()
        if not family or not site or not species:
            raise ValueError(f"row {index + 1} has an empty family_id, site_id, or species")
        if species not in labels:
            raise ValueError(f"row {index + 1} has species absent from tree: {species!r}")
        key = (family, site, species)
        if key in seen:
            raise ValueError(f"duplicate family/site/species row: {key!r}")
        seen.add(key)
        state = row.get("state")
        if isinstance(state, bool) or state not in (0, 1, "0", "1", "unknown"):
            raise ValueError(f"invalid DNA presence state in row {index + 1}: {state!r}")
        state = int(state) if state in {0, 1, "0", "1"} else "unknown"
        eligible = _parse_bool(row.get("eligible"), index)
        bucket = grouped.setdefault((family, site), {"rows": {}, "eligibility": set()})
        bucket["rows"][species] = dict(row, state=state)
        bucket["eligibility"].add(eligible)
    by_family = {}
    for (family, site), bucket in grouped.items():
        missing = labels - set(bucket["rows"])
        extra = set(bucket["rows"]) - labels
        if missing or extra:
            raise ValueError(f"family/site {family!r}/{site!r} tip mismatch; missing={sorted(missing)}, extra={sorted(extra)}")
        if len(bucket["eligibility"]) != 1:
            raise ValueError(f"conflicting eligible flags among species for family/site {family!r}/{site!r}")
        observations = {label: bucket["rows"][label]["state"] for label in labels}
        known = {value for value in observations.values() if value in {0, 1}}
        reasons = sorted({str(item.get("reason", "") or "") for item in bucket["rows"].values()} - {""})
        evidence = sorted({str(item.get("evidence", "") or "") for item in bucket["rows"].values()} - {""})
        by_family.setdefault(family, {})[site] = {
            "eligible": next(iter(bucket["eligibility"])), "observations": observations,
            "known_states": known, "known_states_by_tip": {k: v for k, v in observations.items() if v in {0, 1}},
            "reason": ";".join(reasons), "evidence": ";".join(evidence),
        }
    return by_family


def _parse_bool(value, index):
    if type(value) is bool:
        return value
    if value in ("true", "True", "1", 1):
        return True
    if value in ("false", "False", "0", 0):
        return False
    raise ValueError(f"eligible must be Boolean in row {index + 1}")


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


def _finite_or_none(value):
    return float(value) if isinstance(value, Real) and math.isfinite(float(value)) else None


def _rate_information_status(fit, status):
    if fit is None:
        return status
    if not fit.get("optimizer_success"):
        return "optimizer_unresolved"
    if status == "fixed_rates":
        return "externally_fixed_rates"
    if fit.get("locally_flat"):
        return "information_poor_local_curvature_flat"
    if fit.get("boundary_rates"):
        return "conditional_estimate_at_rate_boundary"
    return "conditional_optimizer_estimate; global optimum not certified"


def _history_for_family(family, eligible, tree, fit, status, root_frequency, root_presence, expected_counts):
    record = {"family_id": family, "status": status, "sites": []}
    if fit is None or not fit.get("optimizer_success") or not fit.get("rates"):
        record["posterior_available"] = False
        record["posterior_unavailable_reason"] = "no_supported_rate_fit_or_fixed_rates"
        return record
    gain, loss = fit["rates"]["gain"], fit["rates"]["loss"]
    rho = None if root_frequency == "stationary" else float(root_presence)
    for site_id, site in sorted(eligible.items()):
        observations = site["observations"]
        record_site = {"site_id": site_id, "status": _site_status(site)}
        if not site["known_states"] or 1 not in site["known_states"]:
            record_site["posterior_available"] = False
            record_site["posterior_unavailable_reason"] = "site_not_selected_by_observed_at_least_one_ascertainment"
            record["sites"].append(record_site)
            continue
        node, edge = _posterior_messages(tree, observations, gain, loss, 1.0, frozenset(), rho)
        record_site["posterior_available"] = True
        record_site["nodes"] = [
            {"node_id": node_id, "label": tree.label[node_id], "p_absent": float(prob[0]), "p_present": float(prob[1])}
            for node_id, prob in sorted(node.items())
        ]
        edges = []
        for (parent, child), joint in sorted(edge.items()):
            item = {
                "parent_node": parent, "child_node": child,
                "p_absent_to_present_endpoint": float(joint[0, 1]),
                "p_present_to_absent_endpoint": float(joint[1, 0]),
                "p_any_endpoint_change": float(joint[0, 1] + joint[1, 0]),
            }
            if expected_counts:
                item["expected_gain_site_state_change_count"] = _expected_transition_count(
                    joint, gain, loss, tree.branch_length(child), 1.0, 0, 1
                )
                item["expected_loss_site_state_change_count"] = _expected_transition_count(
                    joint, gain, loss, tree.branch_length(child), 1.0, 1, 0
                )
            edges.append(item)
        record_site["edges"] = edges
        record["sites"].append(record_site)
    record["posterior_available"] = any(site.get("posterior_available") for site in record["sites"])
    if not record["posterior_available"]:
        record["posterior_unavailable_reason"] = "no_selected_eligible_sites"
    return record
