"""Dense-oracle contracts for exact origin-subset compact posteriors."""
import unittest
from unittest.mock import patch

import numpy as np

from intraphy.inference.configuration_compact import evaluate_compact
from intraphy.inference.configuration_model import RateModel, evaluate_model
from intraphy.inference.configuration_origin_posterior import origin_scenarios
from intraphy.structure.edits import EDIT_KINDS
from intraphy.structure.origins import origin_scenarios as enumerate_origins
from intraphy.structure.space import enumerate_space
from intraphy.structure.types import Catalogue, ExonSpan, Material
from intraphy.topology import SpeciesTree


def _space(material_count=1):
    exon = ExonSpan(0, 2)
    materials = tuple(Material(f"m{i}", 10 + 4*i, 11 + 4*i)
                      for i in range(material_count))
    return enumerate_space(Catalogue("f", "posterior", 16, (exon,), (),
        material=materials, boundary_candidates=(exon,)))


def _polytomy_tree():
    return SpeciesTree([
        {"node_id": "r", "parent_id": "", "label": "r"},
        {"node_id": "a", "parent_id": "r", "label": "A", "branch_length": 0.},
        {"node_id": "x", "parent_id": "r", "label": "x", "branch_length": 0.3},
        {"node_id": "d", "parent_id": "r", "label": "D", "branch_length": 0.4},
        {"node_id": "b", "parent_id": "x", "label": "B", "branch_length": 0.2},
        {"node_id": "c", "parent_id": "x", "label": "C", "branch_length": 0.5},
    ])


def _model(*, foreground=(), root_weight=2.3, zero=False):
    rates = {kind: (0. if zero else 0.035) for kind in EDIT_KINDS}
    if not zero:
        rates["dna_insertion"] = 0.13
        rates["dna_deletion"] = 0.009
    return RateModel(rates, foreground=frozenset(foreground),
        foreground_multiplier=2.4, origin_root_weight=root_weight)


def _observable(state):
    return (tuple((exon.start, exon.end) for exon in state.exons),
            tuple(int(value == 1) for value in state.material))


def _assert_posterior_matches_dense(test, space, tree, tips, model):
    dense = evaluate_model(space, tree, tips, model, counts=False)
    compact = evaluate_compact(space, tree, tips, model, cache_bytes=0)
    test.assertAlmostEqual(compact["log_likelihood"], dense["log_likelihood"], delta=2e-10)

    def weights(result):
        return {tuple(sorted(row["origins"].items())): row["posterior_weight"]
                for row in result["origins"]}

    dense_weights, compact_weights = weights(dense), weights(compact)
    test.assertEqual(set(compact_weights), set(dense_weights))
    test.assertAlmostEqual(sum(compact_weights.values()), 1., places=10)
    for scenario, expected in dense_weights.items():
        test.assertAlmostEqual(compact_weights[scenario], expected, delta=2e-9)
    for node, expected in dense["nodes"].items():
        np.testing.assert_allclose(compact["nodes"][node], expected, atol=2e-9, rtol=0.)

    unique_groups = list(dict.fromkeys(_observable(state) for state in space.states))
    group_ids = {value: i for i, value in enumerate(unique_groups)}
    grouped_ids = np.asarray([group_ids[_observable(state)] for state in space.states])
    dense_branches = {(row["parent"], row["child"]): row for row in dense["branches"]}
    compact_branches = {(row["parent"], row["child"]): row for row in compact["branches"]}
    test.assertEqual(set(compact_branches), set(dense_branches))
    exon_counts = np.asarray([len(state.exons) for state in space.states])
    max_count = int(exon_counts.max(initial=0))
    for edge, expected in dense_branches.items():
        observed = compact_branches[edge]
        summary = observed["observable_endpoint_summary"]
        endpoint = expected["endpoint_probabilities"]
        count_pairs = np.zeros((max_count + 1, max_count + 1))
        grouped_pairs = np.zeros((len(unique_groups), len(unique_groups)))
        for i, parent_count in enumerate(exon_counts):
            for j, child_count in enumerate(exon_counts):
                count_pairs[parent_count, child_count] += endpoint[i, j]
                grouped_pairs[grouped_ids[i], grouped_ids[j]] += endpoint[i, j]
        reported_counts = np.zeros_like(count_pairs)
        for row in summary["exon_count_pair_probabilities"]:
            reported_counts[row["parent_count"], row["child_count"]] = row["probability"]
        np.testing.assert_allclose(reported_counts, count_pairs, atol=2e-9, rtol=0.)
        test.assertAlmostEqual(summary["probability_exon_structure_change"],
            1. - sum(grouped_pairs[i, j] for i, left in enumerate(unique_groups)
                     for j, right in enumerate(unique_groups) if left[0] == right[0]),
            delta=2e-9)
        test.assertAlmostEqual(summary["probability_dna_presence_change"],
            1. - sum(grouped_pairs[i, j] for i, left in enumerate(unique_groups)
                     for j, right in enumerate(unique_groups) if left[1] == right[1]),
            delta=2e-9)
        test.assertAlmostEqual(observed["probability_at_least_one_edit"],
                               expected["probability_at_least_one_edit"], delta=2e-9)
        maximum = grouped_pairs.max()
        expected_modes = {(unique_groups[i], unique_groups[j])
            for i, j in np.argwhere(np.isclose(grouped_pairs, maximum,
                                               atol=1e-12, rtol=1e-12))}
        # Compare the public observable tuples directly, preserving tied pairs.
        normalized_observed = {
            ((tuple(tuple(span) for span in pair["parent_exons"]),
              tuple(pair["parent_dna_presence"])),
             (tuple(tuple(span) for span in pair["child_exons"]),
              tuple(pair["child_dna_presence"])))
            for pair in summary["joint_modes"]}
        test.assertEqual(normalized_observed, expected_modes)
    return dense, compact


class OriginSubsetPosteriorTests(unittest.TestCase):
    def test_adjoint_posterior_on_internal_node_polytomy_and_foreground_zero_edge(self):
        space, tree = _space(), _polytomy_tree()
        tips = {
            label: np.asarray([0.07 + 0.013 * ((i + offset) % 7)
                               for i in range(len(space.states))])
            for label, offset in (("A", 0), ("B", 1), ("C", 2), ("D", 3))
        }
        dense, compact = _assert_posterior_matches_dense(self, space, tree, tips,
            _model(foreground={"x", "b"}))
        self.assertIn("x", compact["nodes"])
        self.assertGreater(len(dense["origins"]), 1)
        zero_edge = next(branch for branch in compact["branches"]
                         if branch["child"] == "a")
        self.assertEqual(zero_edge["probability_at_least_one_edit"], 0.)

    def test_two_material_origin_weights_keep_impossible_scenario_prior_mass(self):
        space = _space(material_count=2)
        tree = SpeciesTree([
            {"node_id": "r", "parent_id": "", "label": "r"},
            {"node_id": "a", "parent_id": "r", "label": "A", "branch_length": 0.3},
            {"node_id": "b", "parent_id": "r", "label": "B", "branch_length": 0.6},
        ])
        tips = {
            "A": np.asarray([float(state.material == (1, 0))
                             for state in space.states]),
            "B": np.asarray([float(state.material == (0, 1))
                             for state in space.states]),
        }
        base_model = _model(root_weight=3.7)
        rates = dict(base_model.rates)
        rates["dna_deletion"] = 0.
        model = RateModel(rates, origin_root_weight=3.7)
        scenarios = list(enumerate_origins(space, tree, tips=tips,
                                           root_weight=model.origin_root_weight))
        self.assertGreater(len({int(root.sum()) for _, root, _ in scenarios}), 1)
        dense, compact = _assert_posterior_matches_dense(self, space, tree, tips, model)
        dense_weights = {tuple(sorted(row["origins"].items())): row["posterior_weight"]
                         for row in dense["origins"]}
        self.assertLess(len(dense_weights), len(scenarios))
        self.assertAlmostEqual(sum(row["posterior_weight"] for row in compact["origins"]),
                               1., places=10)

    def test_compact_posterior_uses_one_explicit_origin_weight_sweep(self):
        space, tree = _space(), _polytomy_tree()
        tips = {label: np.ones(len(space.states)) for label in ("A", "B", "C", "D")}
        original = origin_scenarios
        calls = []

        def counted(*args, **kwargs):
            calls.append(1)
            yield from original(*args, **kwargs)

        with patch("intraphy.inference.configuration_origin_posterior.origin_scenarios",
                   side_effect=counted):
            result = evaluate_compact(space, tree, tips, _model(), cache_bytes=0)
        self.assertEqual(len(calls), 1)
        self.assertEqual(len(result["origins"]), len(list(enumerate_origins(
            space, tree, tips=tips, root_weight=2.3))))


if __name__ == "__main__":
    unittest.main()
