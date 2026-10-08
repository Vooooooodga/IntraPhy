"""Strict parsing of canonical foreground branch selections."""
from __future__ import annotations

import math

from ..storage.tabular import read_tsv


def read_canonical_foreground(path, context, branch_length_mode="supplied"):
    """Resolve a foreground TSV against a tree normalized without foreground flags."""
    if branch_length_mode not in {"supplied", "unit"}:
        raise ValueError("branch_length_mode must be 'supplied' or 'unit'")
    rows = read_tsv(path)
    if not rows:
        raise ValueError("Foreground branch table has no selected edges")
    tree = context.tree
    edge_set = set(tree.edges())
    labels = {}
    for node, label in tree.label.items():
        labels.setdefault(str(label), set()).add(node)

    def endpoint(token, row_number):
        token = str(token).strip()
        candidates = set()
        if token in tree.parent:
            candidates.add(token)
        candidates.update(labels.get(token, ()))
        if len(candidates) != 1:
            reason = "ambiguous" if candidates else "unknown"
            raise ValueError(f"Foreground row {row_number} has {reason} node reference: {token!r}")
        return next(iter(candidates))

    selected, provenance = [], []
    for row_number, row in enumerate(rows, 2):
        raw_parent = str(row.get("parent_node") or row.get("parent_id") or "").strip()
        raw_child = str(row.get("child_node") or row.get("child_id") or "").strip()
        scope = str(row.get("branch_scope") or row.get("branch") or "").strip()
        by_ids = None
        by_scope = None
        if raw_parent or raw_child:
            if not raw_parent or not raw_child:
                raise ValueError(f"Foreground row {row_number} must provide both parent and child")
            by_ids = (endpoint(raw_parent, row_number), endpoint(raw_child, row_number))
        if scope:
            if scope.count("->") != 1:
                raise ValueError(f"Foreground row {row_number} must use one parent->child branch_scope")
            parent_text, child_text = (part.strip() for part in scope.split("->", 1))
            by_scope = (endpoint(parent_text, row_number), endpoint(child_text, row_number))
        if by_ids is None and by_scope is None:
            raise ValueError(f"Foreground row {row_number} needs parent_node/child_node or branch_scope")
        if by_ids is not None and by_scope is not None and by_ids != by_scope:
            raise ValueError(f"Foreground row {row_number} has conflicting branch identifiers")
        edge = by_ids if by_ids is not None else by_scope
        if edge not in edge_set:
            raise ValueError(
                f"Foreground row {row_number} does not identify a canonical branch: {scope or edge}"
            )
        selected.append(edge)
        provenance.append({"row": row_number, "branch_scope": scope or "NA",
                           "parent_node": edge[0], "child_node": edge[1],
                           "input_children_on_canonical_edge": list(context.edge_sources[edge[1]]),
                           "input_row": dict(row)})
    if len(selected) != len(set(selected)):
        raise ValueError("Foreground branch table contains duplicate canonical edges")
    children = frozenset(child for _parent, child in selected)
    if len(children) == len(edge_set):
        raise ValueError("Foreground comparison requires at least one background canonical edge")
    if branch_length_mode == "supplied":
        exposures = {child: tree.branch_length(child) for _parent, child in tree.edges()}
    else:
        exposures = {child: 1.0 for _parent, child in tree.edges()}
    foreground_exposure = math.fsum(exposures[child] for child in children)
    background_exposure = math.fsum(value for child, value in exposures.items()
                                    if child not in children)
    if foreground_exposure <= 0 or background_exposure <= 0:
        raise ValueError("Foreground and background must each have positive branch exposure")
    return children, provenance, foreground_exposure, background_exposure
