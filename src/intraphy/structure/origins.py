"""Single introduction opportunities for each identified variable DNA tract.

A tract may be inherited from the root, or introduced on one edge. After a
physical loss it cannot reappear. The origin opportunity is a discrete latent
scenario, not inferred from the choice of alignment reference.
"""
from __future__ import annotations
from itertools import product
import math
import numpy as np
from .space import StateSpace


def origin_scenarios(space: StateSpace, tree, maximum: int | None = None, *, tips=None, root_weight: float = 1.):
    if maximum is not None and (type(maximum) is not int or maximum < 1):
        raise ValueError("maximum must be a positive integer or None")
    if not math.isfinite(root_weight) or root_weight <= 0:
        raise ValueError("Root opportunity weight must be finite and positive")
    candidates = origin_candidates(space, tree, tips=tips)
    n = math.prod(map(len, candidates))
    if maximum is not None and n > maximum:
        raise ValueError("origin_scenarios_incomplete: increase the explicit scenario limit")
    # Weights are declared before observations. Removed impossible opportunities
    # keep their prior mass; weights are not renormalized after seeing retention.
    denominator = root_weight + len(tree.parent)-1
    root_masks = {}
    for choices in product(*candidates):
        log_prior = sum(math.log(root_weight) if n == tree.root else 0. for n in choices)
        log_prior -= len(candidates)*math.log(denominator)
        origins = dict(zip((m.id for m in space.catalogue.material), choices))
        required = tuple(1 if origins[m.id] == tree.root else 0 for m in space.catalogue.material)
        if required not in root_masks:
            root_masks[required] = np.array(
                [s.material == required for s in space.states], dtype=bool)
        root = root_masks[required].copy()
        if root.any():
            yield origins, root, log_prior


def origin_candidates(space: StateSpace, tree, *, tips=None):
    """Return the emissions-based opportunity nodes for each material tract."""
    candidates = []
    for k, material in enumerate(space.catalogue.material):
        present = [o.species for o in space.catalogue.observations
                   if o.material_presence[k] == 1 and o.species in tree.leaf_by_label]
        if tips is not None:
            present = [species for species, values in tips.items()
                       if any(values > 0) and all(space.states[i].material[k] == 1
                          for i in np.flatnonzero(values > 0))]
        eligible = []
        for node in sorted(tree.parent):
            descendants = set()
            stack = [node]
            while stack:
                v = stack.pop()
                if v in tree.leaves:
                    descendants.add(tree.label[v])
                stack.extend(tree.children.get(v, ()))
            if set(present) <= descendants:
                eligible.append(node)
        candidates.append(tuple(eligible))
    return tuple(candidates)


def permitted(edit, origins, child: str) -> bool:
    return edit.kind != "dna_insertion" or origins.get(edit.material_id) == child
