"""Regression tests for exact parsimony workspace reuse (V19)."""
import unittest
import numpy as np

from intraphy.inference.configuration_dp import SankoffWorkspace, sankoff
from intraphy.inference.configuration_history import ReconstructionWorkspace, reconstruct
from intraphy.structure.edits import EDIT_KINDS
from intraphy.structure.paths import make_graph
from intraphy.structure.space import enumerate_space
from intraphy.structure.types import Catalogue, ExonConfiguration as C, ExonSpan as E, Material
from intraphy.topology import SpeciesTree


def _tree():
    return SpeciesTree([
        {"node_id": "r", "parent_id": "", "label": "r"},
        {"node_id": "x", "parent_id": "r", "label": "x", "branch_length": .2},
        {"node_id": "A", "parent_id": "x", "label": "A", "branch_length": .3},
        {"node_id": "B", "parent_id": "x", "label": "B", "branch_length": .3},
        {"node_id": "C", "parent_id": "r", "label": "C", "branch_length": .4},
    ])


def _space():
    return enumerate_space(Catalogue(
        "g", "u", 20, (E(0, 20), E(0, 8), E(12, 20)), ((8, 12),),
        boundary_candidates=(E(0, 20), E(0, 8), E(12, 20))))


class ConfigurationReuseTests(unittest.TestCase):
    def test_repeated_workspace_reconstruction_is_full_result_equal(self):
        space, tree = _space(), _tree()
        tips = {"A": np.array([1., 0., 0., 0.]), "B": np.array([0., 1., 0., 0.]),
                "C": np.ones(len(space.states))}
        tips = {k: np.eye(len(space.states))[i] if i < len(space.states) else np.ones(len(space.states))
                for i, k in enumerate(("A", "B", "C"))}
        workspace = ReconstructionWorkspace(space, tree)
        first = workspace.reconstruct(tips)
        second = workspace.reconstruct({k: v.copy() for k, v in tips.items()})
        uncached = reconstruct(space, tree, tips)
        self.assertEqual(first, second)
        self.assertEqual(first, uncached)
        self.assertGreater(workspace.graph_cache_hits, 0)
        self.assertGreater(workspace.dp.cache_hits, 0)

    def test_material_origin_graphs_reuse_only_identical_origin_signatures(self):
        catalogue = Catalogue("g", "u", 20, (E(0, 20), E(0, 8), E(12, 20)), ((8, 12),),
                              material=(Material("m", 8, 12),),
                              boundary_candidates=(E(0, 20), E(0, 8), E(12, 20)))
        space, tree = enumerate_space(catalogue), _tree()
        tips = {s: np.ones(len(space.states)) for s in "ABC"}
        workspace = ReconstructionWorkspace(space, tree)
        result = workspace.reconstruct(tips)
        self.assertEqual(result, reconstruct(space, tree, tips))
        constrained = {s: v.copy() for s, v in tips.items()}
        constrained["A"] = np.array([float(state.material == (1,)) for state in space.states])
        self.assertEqual(workspace.reconstruct(constrained), reconstruct(space, tree, constrained))
        self.assertEqual(workspace.reconstruct(tips), reconstruct(space, tree, tips))
        self.assertTrue(result["optimal_origin_scenarios"] >= 1)
        self.assertGreaterEqual(workspace.graph_cache_computes, 1)
        self.assertGreater(workspace.graph_cache_hits, 0)
        self.assertEqual(workspace.diagnostics()["retained_graph_entries"], workspace.graph_cache_computes)

    def test_sankoff_workspace_isolated_by_float64_int64_cost_bytes(self):
        tree = _tree()
        tips = {s: np.array([1., 0., 0.]) for s in "ABC"}
        costs = {c: np.array([[0, 1, 4], [2, 0, 1], [3, 2, 0]], dtype=float) for _, c in tree.edges()}
        workspace = SankoffWorkspace()
        a = sankoff(tree, tips, costs, workspace=workspace)
        integer = {c: value.view(np.int64) for c, value in costs.items()}
        b = sankoff(tree, tips, integer, workspace=workspace)
        uncached = sankoff(tree, tips, integer)
        self.assertEqual(b, uncached)
        self.assertEqual(workspace.cache_hits, 0)

    def test_workspace_does_not_reuse_changed_tips_root_mask_or_topology(self):
        tree, space = _tree(), _space()
        costs = {c: make_graph(space, {}, c).distance for _, c in tree.edges()}
        tips = {s: np.eye(len(space.states))[0] for s in "ABC"}
        workspace = SankoffWorkspace()
        first = sankoff(tree, tips, costs, workspace=workspace)
        changed = dict(tips); changed["B"] = np.eye(len(space.states))[-1]
        second = sankoff(tree, changed, costs, workspace=workspace)
        masked = sankoff(tree, tips, costs, root_allowed=np.array([True] + [False] * (len(space.states)-1)), workspace=workspace)
        self.assertEqual(second, sankoff(tree, changed, costs))
        self.assertEqual(masked, sankoff(tree, tips, costs,
                                         root_allowed=np.array([True] + [False] * (len(space.states)-1))))
        self.assertGreaterEqual(workspace.cache_computes, 1)

    def test_inplace_cost_mutation_matches_uncached(self):
        tree, space = _tree(), _space()
        costs = {c: make_graph(space, {}, c).distance.copy() for _, c in tree.edges()}
        tips = {s: np.eye(len(space.states))[i] for i, s in enumerate("ABC")}
        workspace = SankoffWorkspace()
        sankoff(tree, tips, costs, workspace=workspace)
        child = next(iter(costs)); costs[child][0, 1] += 7
        self.assertEqual(sankoff(tree, tips, costs, workspace=workspace),
                         sankoff(tree, tips, costs))

    def test_changed_topology_prunes_message_cache_and_matches_uncached(self):
        space = _space(); old = _tree()
        tips = {s: np.eye(len(space.states))[i] for i, s in enumerate("ABC")}
        costs = {c: make_graph(space, {}, c).distance for _, c in old.edges()}
        workspace = SankoffWorkspace(); sankoff(old, tips, costs, workspace=workspace)
        star = SpeciesTree([{"node_id":"r","parent_id":"","label":"r"},
            *[{"node_id":s,"parent_id":"r","label":s,"branch_length":.3} for s in "ABC"]])
        star_costs = {c: make_graph(space, {}, c).distance for _, c in star.edges()}
        got = sankoff(star, tips, star_costs, workspace=workspace)
        self.assertEqual(got, sankoff(star, tips, star_costs))
        self.assertLessEqual(len(workspace._messages), len(tuple(star.edges())))

    def test_frozen_cost_copy_rejects_mutation_and_changed_wrapper_costs(self):
        space, tree = _space(), _tree()
        original = {"split": 1.0}
        workspace = ReconstructionWorkspace(space, tree, original)
        with self.assertRaises(TypeError): workspace.costs["split"] = 2.0
        original["split"] = 2.0
        tips = {s: np.ones(len(space.states)) for s in "ABC"}
        self.assertEqual(workspace.reconstruct(tips), reconstruct(space, tree, tips, costs={"split": 1.0}))
        with self.assertRaises(ValueError):
            reconstruct(space, tree, tips, costs=original, workspace=workspace, normalize_tree=False)

    def test_no_compatible_and_cross_space_cases(self):
        space, tree = _space(), _tree()
        tips = {s: np.zeros(len(space.states)) for s in "ABC"}
        workspace = ReconstructionWorkspace(space, tree)
        self.assertEqual(workspace.reconstruct(tips), reconstruct(space, tree, tips))
        other = _space()
        with self.assertRaises(ValueError):
            reconstruct(other, tree, tips, workspace=workspace, normalize_tree=False)

    def test_graph_cache_key_reuses_cross_edge_unintroduced_graph(self):
        catalogue = Catalogue("g", "u", 20, (E(0, 20),), (), material=(Material("m", 0, 20),),
                              boundary_candidates=(E(0, 20),))
        space, tree = enumerate_space(catalogue), _tree()
        workspace = ReconstructionWorkspace(space, tree)
        tips = {s: np.ones(len(space.states)) for s in "ABC"}
        workspace.reconstruct(tips)
        self.assertLess(workspace.graph_cache_computes, len(tuple(tree.edges())))


if __name__ == "__main__":
    unittest.main()
