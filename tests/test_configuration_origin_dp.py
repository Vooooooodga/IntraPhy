"""Regression contracts for likelihood-only origin-subset dynamic programming.

The reference likelihood below deliberately enumerates every origin scenario and
calls the existing single-scenario sparse pruning routine. It remains independent
of evaluate_model's likelihood-only implementation.
"""
import math
import unittest
from dataclasses import replace
from unittest.mock import patch

import numpy as np

from intraphy.inference.configuration_ctmc import likelihood, transition_matrix
from intraphy.inference.configuration_model import RateModel, evaluate_model, generator
from intraphy.inference import configuration_sparse
from intraphy.inference.configuration_sparse import likelihood_only, sparse_generator
from intraphy.structure.edits import EDIT_KINDS
from intraphy.structure.origins import origin_scenarios
from intraphy.structure.space import enumerate_space
from intraphy.structure.tree_context import canonical_tree
from intraphy.structure.types import Catalogue, ExonSpan, Material
from intraphy.topology import SpeciesTree


def tree(lengths=(0.2, 0.3)):
    rows = [{"node_id": "r", "parent_id": "", "label": "r"}]
    for node, label, length in zip(("a", "b"), ("A", "B"), lengths):
        rows.append({"node_id": node, "parent_id": "r", "label": label,
                     "branch_length": length})
    return SpeciesTree(rows)


def space(materials=0, exons=1):
    spans = tuple(ExonSpan(4*i, 4*i+2) for i in range(exons))
    tracts = tuple(Material(f"m{i}", 100+4*i, 101+4*i)
                   for i in range(materials))
    return enumerate_space(Catalogue("f", "u", max(3, 101+4*materials),
        spans, (), material=tracts, boundary_candidates=spans))


def coupled_space():
    left, right = ExonSpan(1, 3), ExonSpan(5, 7)
    catalogue = Catalogue("f", "coupled", 8, (left, right), ((3, 5),),
        material=(Material("m0", 2, 6),), boundary_candidates=(left, right))
    return enumerate_space(catalogue)


def rates(value=0.04, **changes):
    values = {kind: value for kind in EDIT_KINDS}
    values.update(changes)
    return RateModel(values)


def legacy_origin_enumeration(state_space, species_tree, tips, model, *,
                              max_origins=None, branch_length_mode="supplied"):
    """Literal origin-scenario sum with one sparse pruning per scenario."""
    context = canonical_tree(species_tree, model.foreground)
    canonical = context.tree
    normalized_model = RateModel(model.rates, model.scale,
        model.foreground_multiplier, context.foreground, model.origin_root_weight)
    lengths = {child: (1. if branch_length_mode == "unit"
                       else canonical.branch_length(child))
               for _, child in canonical.edges()}
    total = -math.inf
    for origins, root_mask, log_prior in origin_scenarios(
            state_space, canonical, max_origins, tips=tips,
            root_weight=normalized_model.origin_root_weight):
        root = root_mask.astype(float)
        root /= root.sum()
        generators = {child: sparse_generator(
            state_space, normalized_model, origins, child)
            for _, child in canonical.edges()}
        conditional = likelihood_only(canonical, tips, generators, lengths, root)
        total = float(np.logaddexp(total, conditional + log_prior))
    return total


def dense_origin_enumeration(state_space, species_tree, tips, model, *,
                             max_origins=None, branch_length_mode="supplied"):
    """Dense CTMC oracle, independently evaluating every origin scenario."""
    context = canonical_tree(species_tree, model.foreground)
    canonical = context.tree
    normalized_model = RateModel(model.rates, model.scale,
        model.foreground_multiplier, context.foreground, model.origin_root_weight)
    lengths = {child: (1. if branch_length_mode == "unit"
                       else canonical.branch_length(child))
               for _, child in canonical.edges()}
    total = -math.inf
    for origins, root_mask, log_prior in origin_scenarios(
            state_space, canonical, max_origins, tips=tips,
            root_weight=normalized_model.origin_root_weight):
        root = root_mask.astype(float)
        root /= root.sum()
        matrices = {}
        for _, child in canonical.edges():
            q, _ = generator(state_space, normalized_model, origins, child,
                             include_marks=False)
            matrices[child] = transition_matrix(q, lengths[child])
        conditional = likelihood(canonical, tips, matrices, root,
                                 posterior=False).log_likelihood
        total = float(np.logaddexp(total, conditional + log_prior))
    return total


def assert_likelihood_matches_oracle(test, state_space, species_tree, tips,
                                     model, **kwargs):
    expected = legacy_origin_enumeration(state_space, species_tree, tips, model,
                                         **kwargs)
    dense_expected = dense_origin_enumeration(state_space, species_tree, tips,
                                              model, **kwargs)
    test.assertAlmostEqual(expected, dense_expected, delta=2e-11)
    result = evaluate_model(state_space, species_tree, tips, model,
        posterior=False, counts=False, backend="sparse", **kwargs)
    test.assertAlmostEqual(result["log_likelihood"], expected, delta=2e-11)
    return result


class OriginSubsetDPTests(unittest.TestCase):
    def test_k_zero_through_three_with_unknown_and_weighted_emissions(self):
        for k in range(4):
            with self.subTest(materials=k):
                state_space, species_tree = space(k, exons=1), tree()
                unknown = {label: np.ones(len(state_space.states))
                           for label in ("A", "B")}
                assert_likelihood_matches_oracle(self, state_space, species_tree,
                    unknown, rates())
                tips = {
                    "A": np.asarray([1. if i % 3 else 0.17
                                      for i in range(len(state_space.states))], dtype=float),
                    "B": np.asarray([0.11 if i % 2 else 1.
                                      for i in range(len(state_space.states))], dtype=float),
                }
                fit_model = RateModel({**rates().rates, "dna_insertion": 0.13,
                    "dna_deletion": 0.021}, origin_root_weight=2.4)
                assert_likelihood_matches_oracle(self, state_space, species_tree,
                    tips, fit_model)

    def test_multifurcation_foreground_zero_length_and_zero_rate(self):
        state_space = space(2, exons=1)
        species_tree = SpeciesTree([
            {"node_id": "r", "parent_id": "", "label": "r"},
            {"node_id": "a", "parent_id": "r", "label": "A", "branch_length": 0.},
            {"node_id": "b", "parent_id": "r", "label": "B", "branch_length": 0.3},
            {"node_id": "c", "parent_id": "r", "label": "C", "branch_length": 0.5},
        ])
        tips = {label: np.asarray([float((i + offset) % 3 != 0) + 0.08
                                    for i in range(len(state_space.states))])
                for label, offset in (("A", 0), ("B", 1), ("C", 2))}
        foreground_model = RateModel(rates(0.06).rates,
            foreground_multiplier=2.8, foreground=frozenset({"b"}),
            origin_root_weight=0.7)
        assert_likelihood_matches_oracle(self, state_space, species_tree, tips,
                                         foreground_model)
        zero_model = RateModel({kind: 0. for kind in EDIT_KINDS})
        assert_likelihood_matches_oracle(self, state_space, species_tree, tips,
                                         zero_model)

    def test_impossible_support_and_deep_finite_emissions(self):
        state_space = space(exons=1)
        species_tree = tree()
        absent = np.asarray([float(not state.exons) for state in state_space.states])
        present = 1. - absent
        zero_model = RateModel({kind: 0. for kind in EDIT_KINDS})
        impossible = {"A": absent, "B": present}
        self.assertEqual(assert_likelihood_matches_oracle(self, state_space,
            species_tree, impossible, zero_model)["log_likelihood"], -math.inf)
        tiny = 1e-280
        finite = {"A": np.full(len(state_space.states), tiny),
                  "B": np.full(len(state_space.states), tiny)}
        finite_result = assert_likelihood_matches_oracle(self, state_space,
            tree((0., 0.)), finite, zero_model)
        self.assertAlmostEqual(finite_result["log_likelihood"], 2*math.log(tiny),
                               delta=2e-11)

    def test_root_opportunity_weight_and_scenario_root_normalization(self):
        state_space, species_tree = space(2), tree()
        tips = {label: np.ones(len(state_space.states)) for label in ("A", "B")}
        weighted_model = RateModel(rates(0.03).rates, origin_root_weight=3.5)
        scenarios = list(origin_scenarios(state_space, species_tree, tips=tips,
                                          root_weight=weighted_model.origin_root_weight))
        self.assertGreater(len(scenarios), 1)
        self.assertTrue(all(root.any() for _, root, _ in scenarios))
        # Each scenario root vector is normalized separately by the oracle.
        assert_likelihood_matches_oracle(self, state_space, species_tree, tips,
                                         weighted_model)

    def test_impossible_root_origin_keeps_its_declared_prior_mass(self):
        state_space, species_tree = space(1, exons=0), tree()
        absent = np.asarray([float(state.material == (0,))
                             for state in state_space.states])
        present = np.asarray([float(state.material == (1,))
                              for state in state_space.states])
        tips = {"A": absent, "B": present}
        edit_rates = {kind: 0. for kind in EDIT_KINDS}
        edit_rates["dna_insertion"] = 0.2
        fit_model = RateModel(edit_rates, origin_root_weight=2.)
        scenarios = list(origin_scenarios(state_space, species_tree, tips=tips,
                                          root_weight=fit_model.origin_root_weight))
        self.assertEqual({origins["m0"] for origins, _, _ in scenarios}, {"r", "b"})
        root_scenario = next(row for row in scenarios if row[0]["m0"] == "r")
        self.assertAlmostEqual(math.exp(root_scenario[2]), 0.5)
        result = assert_likelihood_matches_oracle(self, state_space, species_tree,
                                                  tips, fit_model)
        # The root-origin scenario contributes zero likelihood; surviving edge
        # scenarios keep their original prior denominator.
        expected_b = next(row for row in scenarios if row[0]["m0"] == "b")
        self.assertAlmostEqual(math.exp(expected_b[2]), 0.25)
        self.assertTrue(np.isfinite(result["log_likelihood"]))
        # Only the b-origin scenario has positive likelihood here: insertion
        # on its 0.3-length branch has probability 1-exp(-0.2*0.3).
        self.assertAlmostEqual(result["log_likelihood"],
            math.log(0.25 * -math.expm1(-0.2 * 0.3)), delta=2e-11)

    def test_canonical_unary_collapse_preserves_likelihood_and_opportunity_count(self):
        state_space = space(1, exons=1)
        subdivided = SpeciesTree([
            {"node_id": "r", "parent_id": "", "label": "r"},
            {"node_id": "u", "parent_id": "r", "label": "u", "branch_length": 0.1},
            {"node_id": "a", "parent_id": "u", "label": "A", "branch_length": 0.2},
            {"node_id": "b", "parent_id": "r", "label": "B", "branch_length": 0.4},
        ])
        tips = {"A": np.asarray([float(state.material == (1,))
                                  for state in state_space.states]),
                "B": np.ones(len(state_space.states))}
        fit_model = RateModel(rates(0.045).rates, origin_root_weight=1.7)
        actual = assert_likelihood_matches_oracle(self, state_space, subdivided,
                                                  tips, fit_model)
        collapsed = canonical_tree(subdivided).tree
        expected = assert_likelihood_matches_oracle(self, state_space, collapsed,
                                                    tips, fit_model)
        self.assertAlmostEqual(actual["log_likelihood"],
                               expected["log_likelihood"], delta=2e-11)
        mixed_regime = RateModel(fit_model.rates, foreground=frozenset({"u"}))
        with self.assertRaisesRegex(ValueError, "Mixed rate regimes"):
            evaluate_model(state_space, subdivided, tips, mixed_regime,
                posterior=False, counts=False, backend="sparse")

    def test_explicit_limit_precedes_root_filtering_and_default_exceeds_64(self):
        state_space, species_tree = space(4, exons=0), tree()
        unknown = {label: np.ones(len(state_space.states)) for label in ("A", "B")}
        scenarios = list(origin_scenarios(state_space, species_tree, tips=unknown))
        self.assertGreater(len(scenarios), 64)
        assert_likelihood_matches_oracle(self, state_space, species_tree, unknown, rates())
        with self.assertRaisesRegex(ValueError, "origin_scenarios_incomplete"):
            evaluate_model(state_space, species_tree, unknown, rates(), max_origins=1,
                           posterior=False, counts=False, backend="sparse")
        # An unbounded call retains the full opportunity mixture by default.
        self.assertAlmostEqual(evaluate_model(state_space, species_tree, unknown,
            rates(), posterior=False, counts=False, backend="sparse")["log_likelihood"],
            0., places=11)

    def test_explicit_limit_counts_assignments_with_impossible_root_masks(self):
        complete = space(1, exons=0)
        state_space = replace(complete,
            states=tuple(state for state in complete.states if state.material != (1,)),
            edits=())
        tips = {label: np.ones(len(state_space.states)) for label in ("A", "B")}
        surviving = list(origin_scenarios(state_space, tree(), tips=tips))
        self.assertEqual(len(surviving), 2)
        with self.assertRaisesRegex(ValueError, "origin_scenarios_incomplete"):
            evaluate_model(state_space, tree(), tips, rates(), max_origins=2,
                posterior=False, counts=False, backend="sparse")

    def test_tip_roster_and_explicit_origin_limit_validation(self):
        state_space, species_tree = space(1), tree()
        tips = {label: np.ones(len(state_space.states)) for label in ("A", "B")}
        for invalid_tips in ({"A": tips["A"]},
                             {**tips, "C": np.ones(len(state_space.states))}):
            with self.subTest(roster=tuple(invalid_tips)):
                with self.assertRaises(ValueError):
                    evaluate_model(state_space, species_tree, invalid_tips, rates(),
                        posterior=False, counts=False, backend="sparse")
        for invalid_limit in (0, -1, 1.5, True):
            with self.subTest(max_origins=invalid_limit):
                with self.assertRaises(ValueError):
                    evaluate_model(state_space, species_tree, tips, rates(),
                        max_origins=invalid_limit, posterior=False, counts=False,
                        backend="sparse")

    def test_origin_dp_reduces_exponential_actions_without_timing_threshold(self):
        state_space, species_tree = space(4, exons=0), tree()
        tips = {label: np.ones(len(state_space.states)) for label in ("A", "B")}
        fit_model = rates(0.05)
        original_action = configuration_sparse.expm_multiply

        def count_actions(run):
            calls = []

            def counted(*args, **kwargs):
                calls.append(1)
                return original_action(*args, **kwargs)

            with patch("intraphy.inference.configuration_sparse.expm_multiply",
                       side_effect=counted):
                value = run()
            return value, len(calls)

        legacy, legacy_actions = count_actions(lambda: legacy_origin_enumeration(
            state_space, species_tree, tips, fit_model))
        optimized, optimized_actions = count_actions(lambda: evaluate_model(
            state_space, species_tree, tips, fit_model, posterior=False,
            counts=False, backend="sparse")["log_likelihood"])
        self.assertAlmostEqual(optimized, legacy, delta=2e-11)
        self.assertGreater(legacy_actions, 64)
        self.assertLess(optimized_actions, legacy_actions)

    def test_nonempty_tract_origins_couple_to_shared_deletion_process(self):
        state_space, species_tree = coupled_space(), tree((0.15, 0.4))
        self.assertTrue(any(edit.kind == "dna_deletion" and len(edit.affected) == 2
                            and edit.source.exons != edit.target.exons
                            for edit in state_space.edits))
        tips = {
            "A": np.asarray([float(state.material == (1,)) for state in state_space.states]),
            "B": np.asarray([1. if state.material == (1,) else 0.23
                             for state in state_space.states]),
        }
        scenarios = list(origin_scenarios(state_space, species_tree, tips=tips))
        self.assertGreater(len({int(root.sum()) for _, root, _ in scenarios}), 1)
        fit_model = RateModel({**rates(0.035).rates,
            "dna_insertion": 0.17, "dna_deletion": 0.012}, origin_root_weight=2.)
        assert_likelihood_matches_oracle(self, state_space, species_tree, tips,
                                         fit_model)

    def test_reachable_subset_action_collapse_is_not_hidden_by_finite_mixture(self):
        state_space, species_tree = space(1, exons=0), tree()
        tips = {label: np.ones(len(state_space.states)) for label in ("A", "B")}
        original = configuration_sparse.log_action
        calls = 0

        def collapse_one_coordinate(q, length, log_values):
            nonlocal calls
            result = original(q, length, log_values)
            if calls == 0:
                result = result.copy()
                result[0] = -np.inf
            calls += 1
            return result

        with patch("intraphy.inference.configuration_origin_dp.log_action",
                   side_effect=collapse_one_coordinate):
            with self.assertRaises(ArithmeticError):
                evaluate_model(state_space, species_tree, tips, rates(0.05),
                    posterior=False, counts=False, backend="sparse")
        self.assertGreater(calls, 0)


if __name__ == "__main__":
    unittest.main()
