"""Read and explicitly prune rooted species trees without changing branch paths."""
from __future__ import annotations

from pathlib import Path

from ..storage.tabular import read_tsv, write_tsv
from ..topology import SpeciesTree


TREE_FIELDS = ["node_id", "parent_id", "label", "branch_length"]


def read_species_tree_rows(path):
    """Return validated rooted-tree rows from TSV or Newick."""
    source = Path(path)
    if source.suffix.lower() == ".tsv":
        rows = read_tsv(source)
        tree = SpeciesTree(rows)
        return _rows(tree)

    from Bio import Phylo

    phylogeny = Phylo.read(str(source), "newick")
    clades = list(phylogeny.find_clades(order="preorder"))
    parents = {id(child): parent for parent in clades for child in parent.clades}
    ids = {id(clade): f"node_{index}" for index, clade in enumerate(clades, 1)}
    rows = [{
        "node_id": ids[id(clade)],
        "parent_id": ids[id(parents[id(clade)])] if id(clade) in parents else "",
        "label": clade.name or ids[id(clade)],
        "branch_length": (0.0 if id(clade) not in parents else
                           clade.branch_length if clade.branch_length is not None else "NA"),
    } for clade in clades]
    tree = SpeciesTree(rows)
    return _rows(tree)


def select_species_tree_rows(path, species, *, prune=False):
    """Select an exact tip panel, preserving the root and every retained path node."""
    rows = read_species_tree_rows(path)
    tree = SpeciesTree(rows)
    requested = list(species)
    if not requested:
        raise ValueError("Requested species panel must not be empty")
    if len(requested) != len(set(requested)):
        raise ValueError("Requested species panel contains duplicate names")
    requested_set = set(requested)
    current_set = set(tree.leaf_by_label)
    missing = sorted(requested_set - current_set)
    extra = sorted(current_set - requested_set)
    if missing:
        raise ValueError("Species tree is missing requested tips: " + ", ".join(missing))
    if extra and not prune:
        raise ValueError("Species tree has extra tips; explicit pruning is required: " + ", ".join(extra))
    if not extra:
        return rows

    retained = set()
    for label in requested_set:
        node = tree.leaf_by_label[label]
        while node:
            retained.add(node)
            node = tree.parent[node]
    return [row for row in rows if row["node_id"] in retained]


def prune_species_tree(tree_path, output_path, retained_species, *, prune_extra_tips: bool):
    """Write a selected tree as native species_tree.tsv rows."""
    rows = select_species_tree_rows(tree_path, retained_species, prune=prune_extra_tips)
    write_tsv(output_path, rows, TREE_FIELDS)


def _rows(tree):
    return [{
        "node_id": node,
        "parent_id": tree.parent[node],
        "label": tree.label[node],
        "branch_length": "" if tree.length[node] is None else tree.length[node],
    } for node in tree.parent]
