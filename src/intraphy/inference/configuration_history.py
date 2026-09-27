"""All optimal origin/structure histories and non-additive elementary edits."""
from __future__ import annotations
from collections import defaultdict
from types import MappingProxyType
import math
import numpy as np

from ..structure.space import StateSpace
from ..structure.origins import origin_scenarios
from ..structure.paths import make_graph
from .configuration_dp import sankoff, SankoffWorkspace
from ..structure.tree_context import canonical_tree


class ReconstructionWorkspace:
    """Per-unit bounded caches for repeated observation reconstructions."""
    def __init__(self, space, tree, costs=None):
        self.space = space
        self.tree = tree
        self.costs_source = tuple(sorted((k, float(v)) for k, v in (costs or {}).items()))
        self.costs = MappingProxyType(dict(self.costs_source)) if costs is not None else None
        self.graphs = {}
        self.dp = SankoffWorkspace()
        self.graph_cache_hits = 0
        self.graph_cache_computes = 0

    def reconstruct(self, tips, *, max_origins=None):
        return _reconstruct(self.space, self.tree, tips, max_origins=max_origins,
                            costs=self.costs, workspace=self)

    def diagnostics(self):
        return {"graph_cache_hits": self.graph_cache_hits,
                "graph_cache_computes": self.graph_cache_computes,
                "message_cache_hits": self.dp.cache_hits,
                "message_cache_computes": self.dp.cache_computes,
                "retained_graph_entries": len(self.graphs),
                "retained_message_entries": len(self.dp._messages)}


def reconstruct(space: StateSpace, tree, tips, *, max_origins=None, costs=None,
                normalize_tree=True, workspace=None):
    if normalize_tree:
        tree = canonical_tree(tree).tree
    if workspace is None:
        return _reconstruct(space, tree, tips, max_origins=max_origins, costs=costs,
                            workspace=None)
    elif (workspace.space is not space or workspace.tree is not tree or
          workspace.costs_source != tuple(sorted((k, float(v)) for k, v in (costs or {}).items()))):
        raise ValueError("ReconstructionWorkspace is bound to a different space, tree, or costs")
    return _reconstruct(space, tree, tips, max_origins=max_origins, costs=costs,
                        workspace=workspace)


def _reconstruct(space: StateSpace, tree, tips, *, max_origins=None, costs=None,
                 workspace=None):
    scenarios, optimum = [], math.inf
    graph_cache = workspace.graphs if workspace else {}
    for origins, root, _ in origin_scenarios(space, tree, max_origins, tips=tips):
        # All edges without an introduction share a graph.
        graphs = {}
        for _, child in tree.edges():
            signature = tuple(sorted(k for k, v in origins.items() if v == child))
            key = signature
            if key not in graph_cache:
                graph_cache[key] = make_graph(space, origins, child, costs)
                if workspace: workspace.graph_cache_computes += 1
            else:
                if workspace: workspace.graph_cache_hits += 1
            graphs[child] = graph_cache[key]
        result = sankoff(tree, tips, {child: graph.distance for child, graph in graphs.items()}, root,
                         workspace=workspace.dp if workspace else None)
        if result.cost < optimum-1e-9:
            optimum, scenarios = result.cost, []
        if math.isfinite(result.cost) and abs(result.cost-optimum) <= 1e-9:
            scenarios.append((origins, result, graphs))
    if not scenarios:
        return {"minimum_cost": None, "status": "no_compatible_history", "events": [], "witness": [], "nodes": {}}
    global_nodes = defaultdict(set)
    global_events = {}
    branch_pairs = defaultdict(set)
    for parent, child in tree.edges():
        alternatives = []
        details = {}
        path_cache = {}
        for origins, result, graphs in scenarios:
            for node, values in result.nodes.items():
                global_nodes[node].update(values)
            graph = graphs[child]
            for source, target in result.pairs[(parent, child)]:
                branch_pairs[(parent, child)].add((source, target))
                key = (id(graph), source, target)
                if key not in path_cache:
                    path_cache[key] = graph.path_event_bounds(source, target)
                alternatives.append(path_cache[key])
                best = graph.distance[source, target]
                for i, j, edit, weight in graph.edges:
                    if np.isclose(graph.distance[source, i]+weight+graph.distance[j, target], best, rtol=0, atol=1e-9):
                        details.setdefault(edit.event_key, []).append(edit)
        event_keys = sorted(set().union(*(a.keys() for a in alternatives)))
        for key in event_keys:
            lower = min(a.get(key, (0, 0))[0] for a in alternatives)
            upper = max(a.get(key, (0, 0))[1] for a in alternatives)
            edits = details[key]
            one = edits[0]
            global_events[(parent, child, key)] = {
                "parent": parent, "child": child, "edit_key": key, "operation": one.kind,
                "start": one.footprint[0], "end": one.footprint[1],
                "support": "required" if lower else "possible", "minimum_count": lower,
                "maximum_count": upper, "affected_spans": sorted({(e.start, e.end) for x in edits for e in x.affected}),
                "consequences": sorted({c for e in edits for c in e.consequences}),
                "molecular_mutation_count": "not_identified"}
    origins, representative, graphs = scenarios[0]
    witness = []
    for parent, child in tree.edges():
        for order, edit in enumerate(graphs[child].witness_path(representative.witness[parent], representative.witness[child]), 1):
            witness.append({"parent": parent, "child": child, "order": order, "edit_key": edit.event_key,
                "operation": edit.kind, "source": space.index[edit.source], "target": space.index[edit.target],
                "consequences": list(edit.consequences), "molecular_mutation_count": "not_identified"})
    unit_costs = not costs or all(float(v) == 1 for v in costs.values())
    if unit_costs and len(witness) != round(optimum):
        raise ArithmeticError("Representative edit history does not attain the global minimum")
    return {"status": "conditional_on_catalogue_and_observation", "minimum_cost": optimum,
            "minimum_structural_edits": int(round(optimum)) if unit_costs else None,
            "representative_edit_count": len(witness), "optimal_origin_scenarios": len(scenarios),
            "representative_origins": origins, "representative_states": representative.witness,
            "events": list(global_events.values()), "witness": witness,
            "nodes": {node: sorted(values) for node, values in global_nodes.items()},
            "pairs": [{"parent": p, "child": c, "pairs": sorted(values)} for (p, c), values in branch_pairs.items()]}
