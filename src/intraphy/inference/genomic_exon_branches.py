"""Summaries of net observable endpoint differences on genomic-exon branches."""
from __future__ import annotations

import numpy as np


COORDINATE_SYSTEM = "local_alignment_interbase_0based_halfopen"
INTERPRETATION_SCOPE = (
    "endpoint net differences; does not estimate event counts or molecular mechanisms"
)
JOINT_MODE_ATOL = 1e-12
JOINT_MODE_RTOL = 1e-12


def observable_configuration(state):
    """Merge hidden unintroduced/deleted material states into DNA absence."""
    return (tuple(tuple(span) for span in state["exons"]),
            tuple(int(value == 1) for value in state.get("material", ())))


def _overlap_groups(parent, child):
    nodes = [("p", span) for span in parent] + [("c", span) for span in child]
    adjacency = {i: set() for i in range(len(nodes))}
    for i, (side_a, a) in enumerate(nodes):
        for j in range(i + 1, len(nodes)):
            side_b, b = nodes[j]
            if side_a != side_b and a[0] < b[1] and b[0] < a[1]:
                adjacency[i].add(j)
                adjacency[j].add(i)
    groups, seen = [], set()
    for start in range(len(nodes)):
        if start in seen:
            continue
        stack, component = [start], []
        seen.add(start)
        while stack:
            node = stack.pop()
            component.append(nodes[node])
            for neighbor in adjacency[node] - seen:
                seen.add(neighbor)
                stack.append(neighbor)
        groups.append(component)
    return groups


def classify_exon_change(parent_exons, child_exons, *, dna_changed=False,
                         changed_material_intervals=None):
    """Classify endpoint span geometry by overlap-connected components."""
    if tuple(parent_exons) == tuple(child_exons):
        return ("no_exon_structure_change",)
    classes = []
    for group in _overlap_groups(parent_exons, child_exons):
        old = [span for side, span in group if side == "p"]
        new = [span for side, span in group if side == "c"]
        old.sort()
        new.sort()
        linked_dna_change = dna_changed
        if changed_material_intervals is not None and old + new:
            left, right = min(span[0] for span in old + new), max(span[1] for span in old + new)
            linked_dna_change = any(start < right and end > left
                                    for start, end in changed_material_intervals)
        if not old:
            classes.append("span_gain")
        elif not new:
            classes.append("span_loss")
        elif len(old) == len(new) == 1:
            if old[0] != new[0]:
                classes.append("boundary_change")
            if linked_dna_change:
                classes.append("dna_coupled")
        elif len(old) == 1 and len(new) > 1:
            split = old[0][0] == new[0][0] and old[0][1] == new[-1][1]
            classes.append("complex_change" if linked_dna_change else
                           "net_split" if split else "complex_change")
            if linked_dna_change:
                classes.append("dna_coupled")
        elif len(old) > 1 and len(new) == 1:
            fusion = old[0][0] == new[0][0] and old[-1][1] == new[0][1]
            classes.append("complex_change" if linked_dna_change else
                           "net_fusion" if fusion else "complex_change")
            if linked_dna_change:
                classes.append("dna_coupled")
        else:
            classes.append("complex_change")
            if linked_dna_change:
                classes.append("dna_coupled")
    if dna_changed and changed_material_intervals is None:
        classes = ["complex_change" if value in {"net_split", "net_fusion"} else value
                   for value in classes]
        classes.append("dna_coupled")
    return tuple(dict.fromkeys(classes))


def _changed_material_intervals(parent_mask, child_mask, materials):
    return tuple((item["start"], item["end"]) for before, after, item in
                 zip(parent_mask, child_mask, materials) if before != after)


def _observable_pair(mode):
    return ((tuple(tuple(int(value) for value in span) for span in mode["parent_exons"]),
             tuple(int(value) for value in mode["parent_dna_presence"])),
            (tuple(tuple(int(value) for value in span) for span in mode["child_exons"]),
             tuple(int(value) for value in mode["child_dna_presence"])))


def _descendant_map(tree):
    descendants = {}
    for node in tree.postorder():
        children = tree.children.get(node, ())
        descendants[node] = (sorted(label for child in children
                                    for label in descendants[child]) if children
                             else [tree.label[node]])
    return descendants


def _row(family, result, parent_id, child_id, descendant_species, parent, child,
         probability_exon_change, probability_dna_change, probability, materials,
         offset, count_summary, mode_tolerance=None, evaluated=None, group_count=None):
    dna_changed = parent[1] != child[1]
    changed_intervals = _changed_material_intervals(parent[1], child[1], materials)
    row = {
        "family_id": family, "unit_id": result["unit_id"],
        "parent_node_id": parent_id, "child_node_id": child_id,
        "descendant_species": descendant_species,
        "probability_exon_structure_change": probability_exon_change,
        "probability_dna_presence_change": probability_dna_change,
        **count_summary,
        "joint_configuration_probability": probability,
        "parent_exons": parent[0], "child_exons": child[0],
        "parent_dna_presence": parent[1], "child_dna_presence": child[1],
        "material_ids": [item["id"] for item in materials],
        "changed_material_tracts": [
            {"material_id": item["id"], "start": item["start"],
             "end": item["end"], "parent_presence": before,
             "child_presence": after}
            for before, after, item in zip(parent[1], child[1], materials)
            if before != after],
        "change_classification": classify_exon_change(
            parent[0], child[0], dna_changed=dna_changed,
            changed_material_intervals=changed_intervals),
        "dna_presence_scope": "declared_material_tracts_only",
        "coordinate_system": COORDINATE_SYSTEM,
        "alignment_offset": offset,
        "interpretation_scope": INTERPRETATION_SCOPE,
    }
    if mode_tolerance is not None:
        row.update(joint_mode_tolerance_absolute=mode_tolerance["absolute"],
                   joint_mode_tolerance_relative=mode_tolerance["relative"],
                   parent_group_rows_evaluated=evaluated,
                   parent_group_count=group_count)
    return row


def branch_change_rows(family_results, unit_details, tree):
    """Aggregate CTMC state-pair masses and report tied modal observable pairs."""
    rows, descendant_species = [], _descendant_map(tree)
    for family_result in family_results:
        family = family_result["family_id"]
        for result in family_result["units"]:
            ctmc = result.get("ctmc")
            if not ctmc or not ctmc.get("branches") or not result.get("states"):
                continue
            detail = unit_details.get((family, result["unit_id"]), {})
            catalogue = detail.get("catalogue", {})
            offset = catalogue.get("alignment_offset", 0)
            materials = catalogue.get("material", ())
            state_configs = [observable_configuration(state) for state in result["states"]]
            configurations = list(dict.fromkeys(state_configs))
            group_index = {configuration: index
                           for index, configuration in enumerate(configurations)}
            state_groups = np.asarray([group_index[value] for value in state_configs])
            exon_ids = {value: index for index, value in enumerate(
                dict.fromkeys(configuration[0] for configuration in configurations))}
            exon_groups = np.asarray([exon_ids[configuration[0]]
                                      for configuration in configurations])
            dna_ids = {value: index for index, value in enumerate(
                dict.fromkeys(configuration[1] for configuration in configurations))}
            dna_groups = np.asarray([dna_ids[configuration[1]]
                                     for configuration in configurations])
            for branch in ctmc["branches"]:
                summary = branch.get("observable_endpoint_summary")
                if summary is not None:
                    tolerance = summary["mode_tolerance"]
                    parent_id, child_id = branch["parent"], branch["child"]
                    for mode in summary["joint_modes"]:
                        parent, child = _observable_pair(mode)
                        rows.append(_row(family, result, parent_id, child_id,
                            descendant_species[child_id], parent, child,
                            summary["probability_exon_structure_change"],
                            summary["probability_dna_presence_change"],
                            mode["joint_configuration_probability"], materials, offset,
                            {key: summary[key] for key in (
                                "exon_count_pair_probabilities",
                                "probability_exon_count_increase",
                                "probability_exon_count_decrease",
                                "probability_exon_count_unchanged")},
                            tolerance, summary["parent_group_rows_evaluated"],
                            summary["parent_group_count"]))
                    continue
                endpoint = branch.get("endpoint_probabilities")
                if endpoint is None:
                    continue
                exon_change = exon_groups[:, None] != exon_groups[None, :]
                dna_change = dna_groups[:, None] != dna_groups[None, :]
                endpoint = np.asarray(endpoint, dtype=float)
                if not np.isfinite(endpoint).all() or (endpoint < 0).any():
                    raise ArithmeticError("Invalid dense exon-count endpoint probabilities")
                grouped = np.zeros((len(configurations), len(configurations)), dtype=float)
                for state_index, group in enumerate(state_groups):
                    grouped[group] += np.bincount(state_groups, weights=endpoint[state_index],
                                                  minlength=len(configurations))
                if grouped.size == 0:
                    continue
                max_probability = float(grouped.max())
                tie_tolerance = JOINT_MODE_ATOL + JOINT_MODE_RTOL * abs(max_probability)
                tied_modes = np.argwhere(np.abs(grouped - max_probability) <= tie_tolerance)
                probability_exon_change = float(grouped[exon_change].sum())
                probability_dna_change = float(grouped[dna_change].sum())
                exon_counts = np.asarray([len(configuration[0])
                                          for configuration in state_configs], dtype=np.int64)
                maximum_count = int(exon_counts.max(initial=0))
                count_joint = np.zeros((maximum_count + 1, maximum_count + 1), dtype=float)
                for state_index, parent_count in enumerate(exon_counts):
                    count_joint[parent_count] += np.bincount(
                        exon_counts, weights=endpoint[state_index],
                        minlength=maximum_count + 1)
                count_pairs = [
                    {"parent_count": parent_count, "child_count": child_count,
                     "probability": float(count_joint[parent_count, child_count])}
                    for parent_count in range(maximum_count + 1)
                    for child_count in range(maximum_count + 1)]
                count_increase = float(np.triu(count_joint, k=1).sum())
                count_decrease = float(np.tril(count_joint, k=-1).sum())
                count_unchanged = float(np.trace(count_joint))
                parent_id, child_id = branch["parent"], branch["child"]
                for parent_group, child_group in tied_modes:
                    parent, child = configurations[parent_group], configurations[child_group]
                    probability = float(grouped[parent_group, child_group])
                    rows.append(_row(family, result, parent_id, child_id,
                        descendant_species[child_id], parent, child,
                        probability_exon_change, probability_dna_change,
                        probability, materials, offset,
                        {"exon_count_pair_probabilities": count_pairs,
                         "probability_exon_count_increase": count_increase,
                         "probability_exon_count_decrease": count_decrease,
                         "probability_exon_count_unchanged": count_unchanged},
                        {"absolute": JOINT_MODE_ATOL, "relative": JOINT_MODE_RTOL},
                        len(configurations), len(configurations)))
    return rows
