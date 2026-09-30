"""Bounded exact search for modal observable endpoint pairs."""
from __future__ import annotations

import math
import numpy as np
from scipy.special import logsumexp

from ..structure.origins import origin_scenarios


def _probability(value, label):
    if not np.isfinite(value) or value < -1e-10 or value > 1. + 1e-10:
        raise ArithmeticError(f"Invalid compact posterior probability: {label}")
    return min(1., max(0., float(value)))


def branch_scalars(payload, tree, lengths, exon_groups, exon_count,
                   dna_groups, dna_count, block_size, masked_columns,
                   action_columns, state_exon_counts):
    """Compute branch change probabilities from grouped exponential actions."""
    ll = payload["log_likelihood"]
    output = {}
    for parent, child in tree.edges():
        q = payload["generators"][child]
        loga = payload["outside"][parent] + payload["siblings"][child]
        logb = payload["inside"][child]
        same = []
        for groups, count in ((exon_groups, exon_count), (dna_groups, dna_count)):
            equal = 0.
            for start in range(0, count, block_size):
                ids = range(start, min(start + block_size, count))
                acted = action_columns(q, lengths[child],
                                       masked_columns(logb, groups, ids))
                parent_values = masked_columns(loga, groups, ids)
                masses = logsumexp(parent_values + acted, axis=0) - ll
                equal += float(np.exp(masses[np.isfinite(masses)]).sum())
            same.append(_probability(equal, "same observable group"))
        no_edit_log = float(logsumexp(
            loga + q.diagonal() * lengths[child] + logb) - ll)
        no_edit = float(np.exp(no_edit_log)) if np.isfinite(no_edit_log) else 0.
        no_edit = _probability(no_edit, "no edit")
        maximum_count = int(np.max(state_exon_counts, initial=0))
        count_joint = np.zeros((maximum_count + 1, maximum_count + 1), dtype=float)
        for start in range(0, maximum_count + 1, block_size):
            child_counts = range(start, min(start + block_size, maximum_count + 1))
            acted = action_columns(q, lengths[child],
                                   masked_columns(logb, state_exon_counts, child_counts))
            for parent_count in range(maximum_count + 1):
                parent_values = np.full_like(loga, -np.inf)
                parent_mask = state_exon_counts == parent_count
                parent_values[parent_mask] = loga[parent_mask]
                log_mass = logsumexp(parent_values[:, None] + acted, axis=0) - ll
                if np.isnan(log_mass).any() or np.isposinf(log_mass).any():
                    raise ArithmeticError("Invalid compact exon-count endpoint mass")
                masses = np.zeros(len(tuple(child_counts)), dtype=float)
                finite = np.isfinite(log_mass)
                masses[finite] = np.exp(log_mass[finite])
                count_joint[parent_count, start:start + len(masses)] += masses
        count_up = float(np.triu(count_joint, k=1).sum())
        count_down = float(np.tril(count_joint, k=-1).sum())
        count_same = float(np.trace(count_joint))
        output[(parent, child)] = {
            "probability_exon_structure_change": _probability(1. - same[0], "exon change"),
            "probability_dna_presence_change": _probability(1. - same[1], "DNA change"),
            "probability_at_least_one_edit": _probability(1. - no_edit, "any edit"),
            "exon_count_pair_probabilities": count_joint,
            "probability_exon_count_increase": _probability(count_up, "exon-count increase"),
            "probability_exon_count_decrease": _probability(count_down, "exon-count decrease"),
            "probability_exon_count_unchanged": _probability(count_same, "unchanged exon count"),
        }
    return output


def modal_rows(space, tree, tips, lengths, max_origins, total_ll,
               origin_root_weight, node_marginals, configurations,
               state_groups, block_size, tie_atol, tie_rtol,
               get_payload, masked_columns, action_columns):
    """Find all numerically tied modes using posterior parent-mass bounds."""
    group_count = len(configurations)
    output = {}
    for parent, child in tree.edges():
        parent_mass = np.bincount(state_groups,
            weights=node_marginals[parent], minlength=group_count)
        order = np.argsort(-parent_mass, kind="stable")
        incumbent, candidates = -1., []
        finished, evaluated = False, 0
        for start in range(0, group_count, block_size):
            selected = order[start:min(start + block_size, group_count)]
            joint = np.zeros((len(selected), group_count), dtype=float)
            for origins, root, log_prior in origin_scenarios(
                    space, tree, max_origins, tips=tips,
                    root_weight=origin_root_weight):
                payload = get_payload(origins, root)
                ll = payload["log_likelihood"]
                if not np.isfinite(ll):
                    continue
                weight = math.exp(ll + log_prior - total_ll)
                loga = payload["outside"][parent] + payload["siblings"][child]
                masked = masked_columns(loga, state_groups, selected)
                acted = action_columns(payload["generators"][child].T,
                                       lengths[child], masked)
                log_joint = acted + payload["inside"][child][:, None] - ll
                state_mass = np.exp(log_joint) * weight
                for column in range(len(selected)):
                    joint[column] += np.bincount(state_groups,
                        weights=state_mass[:, column], minlength=group_count)
            for local, group_id in enumerate(selected):
                evaluated += 1
                row = joint[local]
                row_max = float(row.max())
                if row_max > incumbent:
                    incumbent = row_max
                    tolerance = tie_atol + tie_rtol * abs(incumbent)
                    candidates = [(p, c, value) for p, c, value in candidates
                                  if abs(value - incumbent) <= tolerance]
                tolerance = tie_atol + tie_rtol * abs(incumbent)
                for child_group in np.flatnonzero(np.abs(row - incumbent) <= tolerance):
                    candidates.append((int(group_id), int(child_group),
                                       float(row[child_group])))
                next_row = start + local + 1
                if next_row < group_count:
                    next_bound = float(parent_mass[order[next_row]])
                    if next_bound + tolerance < incumbent:
                        finished = True
                        break
            if finished:
                break
        modes = [{
            "parent_exons": configurations[parent_group][0],
            "child_exons": configurations[child_group][0],
            "parent_dna_presence": configurations[parent_group][1],
            "child_dna_presence": configurations[child_group][1],
            "joint_configuration_probability": probability,
        } for parent_group, child_group, probability in candidates
          if abs(probability - incumbent) <= tie_atol + tie_rtol * abs(incumbent)]
        output[(parent, child)] = {"joint_modes": modes,
            "mode_tolerance": {"absolute": float(tie_atol), "relative": float(tie_rtol)},
            "parent_group_rows_evaluated": evaluated,
            "parent_group_count": group_count}
    return output
