"""Shared fitting and history orchestration for separate binary-presence domains."""
from __future__ import annotations

import math
from numbers import Real

from intraphy.inference.ctmc import _compress_patterns
from intraphy.inference.dna_fit import fit_binary_family
from intraphy.inference.dna_output import write_binary_presence_outputs
from intraphy.inference.binary_presence_history import _history_for_family, _site_status


DNA_DOMAIN = {
    "model": "dna-presence-ctmc", "observation_type": "homologous_dna_presence",
    "analysis_scope": "single-copy-dna-presence", "file_prefix": "dna",
    "state_0": "DNA absent at the homologous position",
    "state_1": "DNA present at the homologous position",
    "transition_01": "structural material gain; source and mechanism are unresolved",
    "transition_10": "structural material loss; mechanism is unresolved",
    "event_interpretation": "expected CTMC site-state changes are not counts of molecular lesions",
    "node_p0": "p_absent", "node_p1": "p_present",
    "edge_p01": "p_absent_to_present_endpoint", "edge_p10": "p_present_to_absent_endpoint",
    "expected_01": "expected_gain_site_state_change_count",
    "expected_10": "expected_loss_site_state_change_count",
    "history_interpretation": "node and edge quantities describe DNA presence-state histories conditional on fitted rates and observed tip states",
}

INTRON_DOMAIN = {
    "model": "intron-position-ctmc", "observation_type": "intron_position_presence",
    "analysis_scope": "aligned-genomic-intron-position-presence", "file_prefix": "intron",
    "state_0": "coding sequence continuous at the homologous coding position",
    "state_1": "intron interrupts coding sequence at the homologous coding position",
    "transition_01": "intron gain at a homologous coding position; mechanism is unresolved",
    "transition_10": "intron loss at a homologous coding position; mechanism is unresolved",
    "event_interpretation": "expected CTMC intron-position site-state changes are not counts of physical lesions",
    "node_p0": "p_no_intron", "node_p1": "p_intron_present",
    "edge_p01": "p_no_intron_to_intron_present_endpoint",
    "edge_p10": "p_intron_present_to_no_intron_endpoint",
    "expected_01": "expected_intron_position_gain_site_state_change_count",
    "expected_10": "expected_intron_position_loss_site_state_change_count",
    "history_interpretation": "node and edge quantities describe intron-position state histories conditional on fitted rates and observed tip states",
}


def analyze_binary_presence(rows, tree, output_dir, *, domain, root_frequency="stationary",
                            root_presence=0.5, fixed_rates=None,
                            branch_length_unit="supplied_tree_units", expected_counts=False,
                            workers=1):
    """Analyze one binary domain per call; domains are never combined in a fit."""
    rows = tuple(rows)
    _validate_options(root_frequency, root_presence, fixed_rates, branch_length_unit, expected_counts, workers)
    _validate_tree(tree)
    families = _group_rows(rows, tree, domain)
    fit_rows, history_rows = [], []
    fixed = dict(fixed_rates or {})
    for family, sites in sorted(families.items()):
        eligible = {site_id: site for site_id, site in sites.items() if site["eligible"]}
        excluded = [site for site_id, site in sorted(sites.items()) if not site["eligible"]]
        informative_sites = {
            site_id: site["observations"] for site_id, site in eligible.items()
            if 0 in site["known_states"] and 1 in site["known_states"]
        }
        selected_sites = {site_id: site for site_id, site in eligible.items() if 1 in site["known_states"]}
        any_contrast = bool(informative_sites)
        known_states = {state for site in eligible.values() for state in site["known_states"]}
        reasons = []
        site_records = []
        for site_id, site in sorted(sites.items()):
            site_records.append({
                "site_id": site_id, "eligible": site["eligible"], "status": _site_status(site),
                "known_tip_count": len(site["known_states_by_tip"]),
                "observed_states": sorted(site["known_states"]),
                "reason": site["reason"], "evidence": site["evidence"],
            })
        fixed_complete = set(fixed) == {"gain", "loss"}
        if not eligible:
            status, fit = "all_excluded", None
            reasons.append("no_eligible_sites")
        elif all(not site["known_states"] for site in eligible.values()):
            status, fit = "all_unknown", None
            reasons.append("eligible_sites_have_no_known_tip_states")
        elif not known_states or known_states == {0}:
            status, fit = "no_observed_presence", None
            reasons.append("no_eligible_site_has_an_observed_present_tip")
        elif not any_contrast and not fixed_complete:
            status, fit = "no_observed_contrast", None
            reasons.append("no_site_has_both_observed_absence_and_presence")
        else:
            fit = fit_binary_family(
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
        fit_rows.append({
            "family_id": family, "status": status, "reasons": reasons,
            "model": domain["model"], "observation_type": domain["observation_type"],
            "analysis_scope": domain["analysis_scope"],
            "root_frequency": root_frequency,
            "root_presence": root_presence if root_frequency == "fixed" else None,
            "branch_length_unit": branch_length_unit,
            "ascertainment": "observed-at-least-one;conditioned on each site's fixed observation mask",
            "n_sites": len(sites), "n_eligible_sites": len(eligible),
            "n_excluded_sites": len(excluded),
            "n_all_unknown_sites": sum(not s["known_states"] for s in eligible.values()),
            "n_variable_sites": len(informative_sites),
            "n_constant_present_sites": sum(s["known_states"] == {1} for s in eligible.values()),
            "n_single_known_present_sites": sum(
                s["known_states"] == {1} and len(s["known_states_by_tip"]) == 1 for s in eligible.values()
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
        })
        history_rows.append(_history_for_family(
            family, eligible, tree, fit, status, root_frequency, root_presence, expected_counts, domain,
        ))
    write_binary_presence_outputs(output_dir, tree, fit_rows, history_rows, domain)
    return {"model": domain["model"], "observation_type": domain["observation_type"],
            "analysis_scope": domain["analysis_scope"], "families": fit_rows, "history": history_rows}


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
    if not tree.leaf_by_label:
        raise ValueError("tree must have at least one tip")
    for _parent, child in tree.edges():
        length = tree.branch_length(child)
        if not math.isfinite(length) or length < 0:
            raise ValueError("all supplied non-root branch lengths must be finite and nonnegative")


def _group_rows(rows, tree, domain):
    labels, grouped, seen = set(tree.leaf_by_label), {}, set()
    for index, row in enumerate(rows):
        observation_type = row.get("observation_type")
        expected_type = domain["observation_type"]
        if observation_type != expected_type and not (domain is DNA_DOMAIN and observation_type in (None, "")):
            raise ValueError(f"row {index + 1} observation_type must be {expected_type!r}")
        layer = row.get("layer")
        expected_layer = "intron_position" if domain is INTRON_DOMAIN else "exon_presence"
        if layer and layer != expected_layer:
            raise ValueError(f"row {index + 1} layer conflicts with observation_type {expected_type!r}")
        family, site, species = (str(row.get(key, "")).strip() for key in ("family_id", "site_id", "species"))
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
            raise ValueError(f"invalid binary presence state in row {index + 1}: {state!r}")
        state = int(state) if state in {0, 1, "0", "1"} else "unknown"
        eligible = _parse_bool(row.get("eligible"), index)
        bucket = grouped.setdefault((family, site), {"rows": {}, "eligibility": set()})
        bucket["rows"][species] = dict(row, state=state)
        bucket["eligibility"].add(eligible)
    by_family = {}
    for (family, site), bucket in grouped.items():
        missing = labels - set(bucket["rows"])
        if missing:
            raise ValueError(f"family/site {family!r}/{site!r} tip mismatch; missing={sorted(missing)}, extra=[]")
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
