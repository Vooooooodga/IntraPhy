"""Blockwise compact branch summaries from origin-subset edge messages."""
from __future__ import annotations

import math
import numpy as np
from scipy.special import logsumexp

from .configuration_origin_dp import _checked_product, _subsets
from .configuration_sparse import _support_after_edge, log_action

ENDPOINT_CATEGORY_MASS_TOLERANCE = 1e-8


def _probability(value, label):
    if not np.isfinite(value) or value < -1e-10 or value > 1. + 1e-10:
        raise ArithmeticError(f"Invalid compact posterior probability: {label}")
    return min(1., max(0., float(value)))


def action_columns(q, length, log_values):
    """Apply a sparse exponential with independent column scaling and checks."""
    if np.isnan(log_values).any() or np.isposinf(log_values).any():
        raise ArithmeticError("Log likelihood input contains NaN or positive infinity")
    shifts = np.max(log_values, axis=0)
    valid = np.isfinite(shifts)
    result = np.full(log_values.shape, -np.inf, dtype=float)
    if valid.any():
        normalized = log_values[:, valid] - shifts[valid][None, :]
        acted = log_action(q, length, normalized)
        for column in range(acted.shape[1]):
            if np.isneginf(acted[:, column]).any():
                reachable = _support_after_edge(
                    q, length, np.isfinite(normalized[:, column]))
                if np.any(reachable & ~np.isfinite(acted[:, column])):
                    raise ArithmeticError(
                        "Sparse origin-subset action collapsed a reachable likelihood")
        result[:, valid] = acted + shifts[valid][None, :]
    if np.isnan(result).any() or np.isposinf(result).any():
        raise ArithmeticError("Sparse origin-subset action produced a nonfinite log likelihood")
    return result


def masked_columns(log_vector, groups, group_ids):
    values = np.full((len(log_vector), len(group_ids)), -np.inf, dtype=float)
    for column, group_id in enumerate(group_ids):
        mask = groups == group_id
        values[mask, column] = log_vector[mask]
    return values


def _components(child, inside, edge_outside, node_masks, edge_generator):
    for below, child_value in inside[child].items():
        for edge_mask in _subsets(node_masks[child]):
            if edge_mask & below:
                continue
            combined = edge_mask | below
            if combined in edge_outside[child]:
                yield (edge_generator(child, edge_mask),
                       edge_outside[child][combined], child_value)


def branch_scalars(tree, lengths, inside, edge_outside, node_masks,
                   edge_generator, exon_groups, exon_count,
                   dna_groups, dna_count, joint_groups, joint_count,
                   state_exon_counts, block_size, total):
    """Marginalize observable branch statistics across origin subsets."""
    output = {}
    maximum_count = int(np.max(state_exon_counts, initial=0))
    for parent, child in tree.edges():
        equal = [0., 0.]
        equal_joint = 0.
        no_edit = 0.
        can_edit = False
        count_joint = np.zeros((maximum_count + 1, maximum_count + 1))
        for q, loga, logb in _components(child, inside, edge_outside,
                                          node_masks, edge_generator):
            can_edit |= lengths[child] > 0 and q.nnz > 0
            no_edit_log = float(logsumexp(
                loga + q.diagonal() * lengths[child] + logb) - total)
            if np.isfinite(no_edit_log):
                no_edit += math.exp(no_edit_log)
            groupings = [(0, exon_groups, exon_count),
                         (1, dna_groups, dna_count)]
            if exon_count > 1 and dna_count > 1:
                groupings.append((2, joint_groups, joint_count))
            for slot, groups, count in groupings:
                for start in range(0, count, block_size):
                    ids = range(start, min(start + block_size, count))
                    acted = action_columns(q, lengths[child],
                                           masked_columns(logb, groups, ids))
                    parent_values = masked_columns(loga, groups, ids)
                    masses = logsumexp(parent_values + acted, axis=0) - total
                    if np.isnan(masses).any() or np.isposinf(masses).any():
                        raise ArithmeticError(
                            "Invalid compact observable endpoint equality mass")
                    mass = float(np.exp(masses[np.isfinite(masses)]).sum())
                    if slot == 2:
                        equal_joint += mass
                    else:
                        equal[slot] += mass
            for start in range(0, maximum_count + 1, block_size):
                child_counts = range(start, min(start + block_size, maximum_count + 1))
                acted = action_columns(q, lengths[child],
                    masked_columns(logb, state_exon_counts, child_counts))
                for parent_count in range(maximum_count + 1):
                    parent_mask = state_exon_counts == parent_count
                    log_mass = logsumexp(
                        np.where(parent_mask, loga, -np.inf)[:, None] + acted,
                        axis=0) - total
                    if np.isnan(log_mass).any() or np.isposinf(log_mass).any():
                        raise ArithmeticError("Invalid compact exon-count endpoint mass")
                    masses = np.zeros(len(tuple(child_counts)), dtype=float)
                    finite = np.isfinite(log_mass)
                    masses[finite] = np.exp(log_mass[finite])
                    count_joint[parent_count, start:start + len(masses)] += masses
        count_up = float(np.triu(count_joint, k=1).sum())
        count_down = float(np.tril(count_joint, k=-1).sum())
        count_same = float(np.trace(count_joint))
        any_edit = 1. - no_edit if can_edit else 0.
        if dna_count == 1:
            equal_joint = equal[0]
        elif exon_count == 1:
            equal_joint = equal[1]
        categories = tuple(_probability(value, label) for value, label in zip((
            equal_joint, equal[1] - equal_joint, equal[0] - equal_joint,
            1. - equal[0] - equal[1] + equal_joint), (
            "neither changed", "exon-only change", "material-only change",
            "joint change")))
        category_mass = sum(categories)
        if abs(category_mass - 1.) > ENDPOINT_CATEGORY_MASS_TOLERANCE:
            raise ArithmeticError(
                "Joint endpoint category probabilities do not sum to one "
                f"within {ENDPOINT_CATEGORY_MASS_TOLERANCE:g}: {category_mass}")
        output[(parent, child)] = {
            "probability_exon_structure_change": _probability(1. - equal[0], "exon change"),
            "probability_dna_presence_change": _probability(1. - equal[1], "DNA change"),
            "probability_neither_changed": categories[0],
            "probability_exon_only_changed": categories[1],
            "probability_material_only_changed": categories[2],
            "probability_both_changed": categories[3],
            "probability_at_least_one_edit": _probability(any_edit, "any edit"),
            "exon_count_pair_probabilities": count_joint,
            "probability_exon_count_increase": _probability(count_up, "exon-count increase"),
            "probability_exon_count_decrease": _probability(count_down, "exon-count decrease"),
            "probability_exon_count_unchanged": _probability(count_same, "unchanged exon count"),
        }
    return output


def modal_rows(tree, lengths, inside, edge_outside, node_marginals,
               node_masks, edge_generator, configurations, state_groups,
               block_size, tie_atol, tie_rtol, total):
    """Find tied observable endpoint modes with the existing parent-mass bound."""
    group_count = len(configurations)
    output = {}
    for parent, child in tree.edges():
        parent_mass = np.bincount(state_groups, weights=node_marginals[parent],
                                  minlength=group_count)
        order = np.argsort(-parent_mass, kind="stable")
        incumbent, candidates = -1., []
        finished, evaluated = False, 0
        for start in range(0, group_count, block_size):
            selected = order[start:min(start + block_size, group_count)]
            joint = np.zeros((len(selected), group_count))
            for q, loga, logb in _components(child, inside, edge_outside,
                                              node_masks, edge_generator):
                masked = masked_columns(loga, state_groups, selected)
                acted = action_columns(q.T, lengths[child], masked)
                state_log_mass = _checked_product(acted, logb[:, None]) - total
                if np.isnan(state_log_mass).any() or np.isposinf(state_log_mass).any():
                    raise ArithmeticError("Invalid compact joint endpoint mass")
                state_mass = np.exp(state_log_mass)
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
                    if parent_mass[order[next_row]] + tolerance < incumbent:
                        finished = True
                        break
            if finished:
                break
        modes = [{
            "parent_exons": configurations[p][0],
            "child_exons": configurations[c][0],
            "parent_dna_presence": configurations[p][1],
            "child_dna_presence": configurations[c][1],
            "joint_configuration_probability": probability,
        } for p, c, probability in candidates
          if abs(probability - incumbent) <= tie_atol + tie_rtol * abs(incumbent)]
        output[(parent, child)] = {"joint_modes": modes,
            "mode_tolerance": {"absolute": float(tie_atol), "relative": float(tie_rtol)},
            "parent_group_rows_evaluated": evaluated,
            "parent_group_count": group_count}
    return output
