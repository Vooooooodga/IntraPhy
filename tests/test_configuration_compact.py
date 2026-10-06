"""Dense-oracle tests for exact compact configuration posteriors."""
import unittest

import numpy as np
from scipy.sparse import csr_matrix

from intraphy.inference.configuration_compact import _action_columns, evaluate_compact
from intraphy.inference.configuration_model import RateModel, evaluate_model
from intraphy.structure.edits import EDIT_KINDS
from intraphy.structure.space import enumerate_space
from intraphy.structure.types import Catalogue, ExonSpan, Material
from intraphy.topology import SpeciesTree


def _tree(length_a=0.2, length_b=0.3):
    return SpeciesTree([
        {"node_id": "r", "parent_id": "", "label": "r"},
        {"node_id": "a", "parent_id": "r", "label": "A", "branch_length": length_a},
        {"node_id": "b", "parent_id": "r", "label": "B", "branch_length": length_b},
    ])


def _space(exons=1, materials=0):
    spans = tuple(ExonSpan(3*i, 3*i+2) for i in range(exons))
    tracts = tuple(Material(f"m{i}", 100+10*i, 101+10*i) for i in range(materials))
    return enumerate_space(Catalogue("f", "u", max(3*exons, 101+10*materials),
        spans, (), material=tracts, boundary_candidates=spans))


def _model(*, zero=False, origin_root_weight=1.):
    return RateModel({kind: (0. if zero else 0.04) for kind in EDIT_KINDS},
                     origin_root_weight=origin_root_weight)


def _observable(state):
    return (tuple((e.start, e.end) for e in state.exons),
            tuple(int(value == 1) for value in state.material))


def _dense_grouped(space, endpoint):
    configs = list(dict.fromkeys(_observable(state) for state in space.states))
    ids = {value: i for i, value in enumerate(configs)}
    state_ids = [ids[_observable(state)] for state in space.states]
    grouped = np.zeros((len(configs), len(configs)))
    for i, parent in enumerate(state_ids):
        for j, child in enumerate(state_ids):
            grouped[parent, child] += endpoint[i, j]
    return configs, grouped


class ConfigurationCompactTests(unittest.TestCase):
    def test_action_columns_rejects_nonfinite_positive_mass_and_preserves_zero_columns(self):
        q = csr_matrix((2, 2), dtype=float)
        for invalid in (np.nan, np.inf):
            with self.subTest(invalid=invalid):
                values = np.asarray([[0., invalid], [-1., -np.inf]])
                with self.assertRaisesRegex(ArithmeticError, "NaN or positive infinity"):
                    _action_columns(q, 0.5, values)
        values = np.asarray([[0., -np.inf], [-2., -np.inf]])
        result = _action_columns(q, 0.5, values)
        np.testing.assert_array_equal(result, values)
        all_zero = _action_columns(q, 0., np.full((2, 2), -np.inf))
        self.assertTrue(np.isneginf(all_zero).all())

    def _compare_dense(self, space, tree, tips, model, **compact_options):
        dense = evaluate_model(space, tree, tips, model, counts=False)
        compact = evaluate_compact(space, tree, tips, model, **compact_options)
        self.assertAlmostEqual(compact["log_likelihood"], dense["log_likelihood"], places=10)
        self.assertEqual(compact["posterior_kind"], "compact_observable_pairs")
        self.assertAlmostEqual(sum(x["posterior_weight"] for x in compact["origins"]), 1.)
        for node, probabilities in dense["nodes"].items():
            np.testing.assert_allclose(compact["nodes"][node], probabilities, atol=2e-10)
        dense_branches = {(b["parent"], b["child"]): b for b in dense["branches"]}
        compact_branches = {(b["parent"], b["child"]): b for b in compact["branches"]}
        self.assertEqual(set(dense_branches), set(compact_branches))
        for edge, expected in dense_branches.items():
            branch = compact_branches[edge]
            self.assertNotIn("endpoint_probabilities", branch)
            groups, grouped = _dense_grouped(space, expected["endpoint_probabilities"])
            same_exons = sum(grouped[i, j] for i, a in enumerate(groups)
                             for j, b in enumerate(groups) if a[0] == b[0])
            same_dna = sum(grouped[i, j] for i, a in enumerate(groups)
                           for j, b in enumerate(groups) if a[1] == b[1])
            summary = branch["observable_endpoint_summary"]
            exon_counts = np.asarray([len(state.exons) for state in space.states])
            maximum_count = int(exon_counts.max(initial=0))
            expected_count_pairs = np.zeros((maximum_count + 1, maximum_count + 1))
            for parent_state, parent_count in enumerate(exon_counts):
                for child_state, child_count in enumerate(exon_counts):
                    expected_count_pairs[parent_count, child_count] += expected[
                        "endpoint_probabilities"][parent_state, child_state]
            reported_pairs = summary["exon_count_pair_probabilities"]
            self.assertEqual(len(reported_pairs), (maximum_count + 1) ** 2)
            observed_count_pairs = np.zeros_like(expected_count_pairs)
            for pair in reported_pairs:
                observed_count_pairs[pair["parent_count"], pair["child_count"]] = pair["probability"]
            np.testing.assert_allclose(observed_count_pairs, expected_count_pairs, atol=2e-10)
            self.assertAlmostEqual(float(observed_count_pairs.sum()), 1., places=9)
            self.assertAlmostEqual(summary["probability_exon_count_increase"],
                                   float(np.triu(expected_count_pairs, 1).sum()), places=9)
            self.assertAlmostEqual(summary["probability_exon_count_decrease"],
                                   float(np.tril(expected_count_pairs, -1).sum()), places=9)
            self.assertAlmostEqual(summary["probability_exon_count_unchanged"],
                                   float(np.trace(expected_count_pairs)), places=9)
            for node, axis in ((edge[0], 0), (edge[1], 1)):
                expected_node_counts = np.bincount(exon_counts,
                    weights=dense["nodes"][node], minlength=maximum_count + 1)
                np.testing.assert_allclose(observed_count_pairs.sum(axis=1-axis),
                                           expected_node_counts, atol=2e-10)
            self.assertAlmostEqual(summary["probability_exon_structure_change"],
                                   1. - same_exons, places=9)
            self.assertAlmostEqual(summary["probability_dna_presence_change"],
                                   1. - same_dna, places=9)
            self.assertAlmostEqual(branch["probability_at_least_one_edit"],
                                   expected["probability_at_least_one_edit"], places=9)
            maximum = float(grouped.max())
            modes = {(groups[i], groups[j]) for i, j in np.argwhere(
                np.isclose(grouped, maximum, atol=1e-12, rtol=1e-12))}
            observed = {(tuple(pair["parent_exons"]), tuple(pair["parent_dna_presence"]),
                         tuple(pair["child_exons"]), tuple(pair["child_dna_presence"]))
                        for pair in summary["joint_modes"]}
            expected_modes = {(p[0], p[1], c[0], c[1]) for p, c in modes}
            self.assertEqual(observed, expected_modes)
        return compact

    def test_nonreversible_origin_mixture_and_hidden_material_collapse(self):
        space, tree = _space(exons=1, materials=1), _tree()
        tips = {
            "A": np.asarray([float(state.material == (1,)) for state in space.states]),
            "B": np.ones(len(space.states)),
        }
        rates = {kind: 0.025 for kind in EDIT_KINDS}
        rates["dna_insertion"], rates["dna_deletion"] = 0.11, 0.015
        model = RateModel(rates, origin_root_weight=3.)
        result = self._compare_dense(space, tree, tips, model, block_size=1,
                                     cache_bytes=0)
        self.assertGreater(len(result["origins"]), 1)
        self.assertTrue(any(state.material == (0,) for state in space.states))
        self.assertTrue(any(state.material == (2,) for state in space.states))

    def test_informative_mixed_exon_material_and_foreground_rates(self):
        space, tree = _space(exons=2, materials=1), _tree(0.15, 0.4)
        tips = {
            "A": np.asarray([float(bool(state.exons) and state.material == (1,))
                             for state in space.states]),
            "B": np.asarray([1. if not state.exons else 0.2
                             for state in space.states]),
        }
        model = RateModel({kind: 0.035 for kind in EDIT_KINDS},
            foreground=frozenset({"a"}), foreground_multiplier=2.7,
            origin_root_weight=1.8)
        result = self._compare_dense(space, tree, tips, model, block_size=2)
        self.assertTrue(any(branch["observable_endpoint_summary"]
                            ["probability_exon_structure_change"] > 0
                            for branch in result["branches"]))
        self.assertTrue(any(branch["observable_endpoint_summary"]
                            ["probability_dna_presence_change"] > 0
                            for branch in result["branches"]))

    def test_symmetric_ties_and_near_distinct_mode(self):
        space, tree = _space(exons=1), _tree()
        tied_tips = {label: np.ones(len(space.states)) for label in ("A", "B")}
        tied = evaluate_compact(space, tree, tied_tips, _model(zero=True))
        self.assertTrue(all(len(branch["observable_endpoint_summary"]["joint_modes"]) ==
                            len(space.states) for branch in tied["branches"]))
        tips = {"A": np.ones(len(space.states)),
                "B": np.asarray([1., 1. - 1e-9])}
        distinct = evaluate_compact(space, tree, tips, _model(zero=True))
        self.assertTrue(all(len(branch["observable_endpoint_summary"]["joint_modes"]) == 1
                            for branch in distinct["branches"]))

    def test_zero_length_rates_unknown_tips_and_large_complete_space(self):
        small, tree = _space(exons=1, materials=1), _tree(0., 0.)
        unknown = {label: np.ones(len(small.states)) for label in ("A", "B")}
        result = self._compare_dense(small, tree, unknown, _model(zero=True),
                                     block_size=1, cache_bytes=0)
        self.assertTrue(all(branch["probability_at_least_one_edit"] == 0.
                            for branch in result["branches"]))
        large = _space(exons=5)
        self.assertGreater(len(large.states), 64)
        tips = {label: np.ones(len(large.states)) for label in ("A", "B")}
        large_result = self._compare_dense(large, _tree(), tips, _model(zero=True),
                                            block_size=16, cache_bytes=0)
        self.assertTrue(all("endpoint_probabilities" not in branch
                            for branch in large_result["branches"]))
        self.assertTrue(all(branch["probability_at_least_one_edit"] == 0.
                            for branch in large_result["branches"]))

    def test_counts_are_rejected_in_compact_backend(self):
        space, tree = _space(), _tree()
        tips = {label: np.ones(len(space.states)) for label in ("A", "B")}
        with self.assertRaisesRegex(ValueError, "counts=False"):
            evaluate_compact(space, tree, tips, _model(), counts=True)
