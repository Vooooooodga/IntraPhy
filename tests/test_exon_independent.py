import inspect
import math
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from intraphy.inference.ctmc import _transition_matrix
from intraphy.structure.types import (Catalogue, ExonConfiguration, ExonSpan,
                                      Material, ObservationEvidence)
from intraphy.structure.space import StateSpace
from intraphy.topology import SpeciesTree
from intraphy.verification.exon_independent import (
    _evaluate, _features, _fit_rate, _illegal_node_mass, _root_feature_marginals,
    fit_independent_exon_characters, score_independent_endpoint_posteriors,
)


def _tree():
    return SpeciesTree([
        {"node_id": "R", "parent_id": "", "branch_length": ""},
        {"node_id": "A", "parent_id": "R", "label": "A", "branch_length": "0.7"},
        {"node_id": "B", "parent_id": "R", "label": "B", "branch_length": "1.2"},
    ])


def _space(states, *, materials=(), observations=()):
    spans = tuple(sorted({e for state in states for e in state.exons}))
    catalogue = Catalogue(family="g", unit="u", length=20, spans=spans,
        junctions=(), material=tuple(materials), observations=tuple(observations),
        discovery="independent_catalogue", observation_unit="genomic_exon_spans")
    return StateSpace(catalogue, tuple(states), (), True, "test", None)


def _unit(space, tips):
    return {"space": space, "tree": _tree(), "tips": tips}


class IndependentExonBaselineTests(unittest.TestCase):
    def test_dense_two_state_likelihood_matches_pruning_and_zero_rate(self):
        span = ExonSpan(2, 5)
        states = (ExonConfiguration((), ()), ExonConfiguration((span,), ()))
        space = _space(states)
        tree = _tree()
        features, _bits = _features(space)
        feature_id = features[0]["feature_id"]
        pattern = {"A": 1, "B": 0}
        observations = {feature_id: pattern}
        root = {feature_id: 0.5}
        foreground = frozenset({"A"})
        actual = _evaluate(space, tree, observations, root, 0.31, 1.7,
                           foreground, {"A": 0.7, "B": 1.2})
        matrix_a = _transition_matrix(0.31, 0.31, 0.7, 1.7)
        matrix_b = _transition_matrix(0.31, 0.31, 1.2, 1.0)
        expected = math.log(0.5 * matrix_a[0, 1] * matrix_b[0, 0] +
                            0.5 * matrix_a[1, 1] * matrix_b[1, 0])
        self.assertAlmostEqual(actual, expected, places=12)
        zero = _evaluate(space, tree, observations, root, 0.0, 1.7,
                         foreground, {"A": 0.7, "B": 1.2})
        self.assertEqual(zero, -math.inf)

    def test_zero_rate_boundary_and_complete_missing_tip(self):
        span = ExonSpan(2, 5)
        states = (ExonConfiguration((), ()), ExonConfiguration((span,), ()))
        space = _space(states)
        tips = {"A": np.array([1.0, 0.0]), "B": np.array([1.0, 0.0])}
        fit = fit_independent_exon_characters(_unit(space, tips), {"A"})
        self.assertEqual(fit["status"], "zero_boundary")
        self.assertEqual(fit["rate"], 0.0)
        self.assertEqual(fit["tip_observations"]["exon:2:5"]["B"], 0)
        unknown_tips = {"A": np.array([1.0, 0.0]), "B": np.ones(2)}
        unknown = fit_independent_exon_characters(_unit(space, unknown_tips), {"A"})
        self.assertIsNone(unknown["tip_observations"]["exon:2:5"]["B"])

    def test_full_zero_feature_is_retained(self):
        span = ExonSpan(2, 5)
        states = (ExonConfiguration((), ()), ExonConfiguration((span,), ()))
        space = _space(states)
        features, bits = _features(space)
        self.assertEqual([f["feature_id"] for f in features], ["exon:2:5"])
        self.assertEqual(bits[:, 0].tolist(), [False, True])
        tips = {"A": np.array([1.0, 0.0]), "B": np.array([1.0, 0.0])}
        fit = fit_independent_exon_characters(_unit(space, tips), {"A"})
        self.assertEqual(fit["tip_observations"]["exon:2:5"], {"A": 0, "B": 0})

    def test_root_marginal_uses_full_origin_prior_not_observed_presence_filter(self):
        material = Material("m", 4, 8)
        states = tuple(ExonConfiguration((), (value,)) for value in (0, 1, 2))
        observation = ObservationEvidence("A", (states[1],), "observed", (), (), (1,))
        space = _space(states, materials=(material,), observations=(observation,))
        features, bits = _features(space)
        probabilities, prior = _root_feature_marginals(space, _tree(), features, bits)
        self.assertAlmostEqual(probabilities["material:m"], 1.0 / 3.0, places=12)
        self.assertEqual(prior["origin_scenarios"], 3)
        self.assertFalse(prior["joint_root_distribution_matches_native"])

    def test_partial_observations_are_rejected_and_whole_tip_unknown_is_preserved(self):
        span = ExonSpan(2, 5)
        states = (ExonConfiguration((), ()), ExonConfiguration((span,), ()))
        partial = ObservationEvidence("A", (states[0],), "partial", (), (), ())
        space = _space(states, observations=(partial,))
        with self.assertRaisesRegex(ValueError, "partial"):
            fit_independent_exon_characters(_unit(space, {"A": [1, 0], "B": [1, 0]}), {"A"})
        clean = _space(states)
        fit = fit_independent_exon_characters(
            _unit(clean, {"A": [1, 0], "B": [1, 1]}), {"A"})
        self.assertIsNone(fit["tip_observations"]["exon:2:5"]["B"])

    def test_all_unknown_tips_are_flat_and_rate_is_unidentified(self):
        span = ExonSpan(2, 5)
        space = _space((ExonConfiguration((), ()), ExonConfiguration((span,), ())))
        fit = fit_independent_exon_characters(
            _unit(space, {"A": np.ones(2), "B": np.ones(2)}), {"A"})
        self.assertEqual(fit["status"], "no_rate_information")
        self.assertIsNone(fit["rate"])
        self.assertIsNone(fit["branch_posteriors"])

    def test_unsuccessful_optimizer_is_not_scored(self):
        span = ExonSpan(2, 5)
        space = _space((ExonConfiguration((), ()), ExonConfiguration((span,), ())))
        failed = SimpleNamespace(x=np.array([0.2]), success=False, message="iteration limit")
        with patch("intraphy.verification.exon_independent.minimize", return_value=failed):
            fit = fit_independent_exon_characters(
                _unit(space, {"A": [1, 0], "B": [0, 1]}), {"A"})
        self.assertEqual(fit["status"], "optimization_unresolved")
        self.assertIsNone(fit["branch_posteriors"])
        self.assertEqual(score_independent_endpoint_posteriors(fit, {}),
                         {"status": "no_call", "edges": []})

    def test_failed_candidate_above_successful_zero_is_unresolved(self):
        failed = SimpleNamespace(x=np.array([2.0]), success=False, message="did not converge")
        objective = lambda _space, _tree, _observations, _root, rate, _rho, _fg, _lengths: -(rate - 2.0) ** 2
        with patch("intraphy.verification.exon_independent.minimize", return_value=failed), \
                patch("intraphy.verification.exon_independent._evaluate", side_effect=objective):
            fit = _fit_rate(None, None, {}, {}, 1.0, {"A"}, {"A": 1.0, "B": 1.0})
        self.assertEqual(fit["status"], "optimization_unresolved")

    def test_lower_failed_candidates_do_not_block_successful_best_fit(self):
        results = [SimpleNamespace(x=np.array([2.0]), success=True, message="converged")]
        results.extend(SimpleNamespace(x=np.array([10.0]), success=False, message="lower failed point")
                       for _ in range(4))
        objective = lambda _space, _tree, _observations, _root, rate, _rho, _fg, _lengths: -(rate - 2.0) ** 2
        with patch("intraphy.verification.exon_independent.minimize", side_effect=results), \
                patch("intraphy.verification.exon_independent._evaluate", side_effect=objective):
            fit = _fit_rate(None, None, {}, {}, 1.0, {"A"}, {"A": 1.0, "B": 1.0})
        self.assertEqual(fit["status"], "estimated")
        self.assertEqual(fit["rate"], 2.0)
        self.assertIn("lower failed candidates", fit["diagnostic"])

    def test_illegal_node_mass_is_not_renormalized(self):
        left, right = ExonSpan(2, 7), ExonSpan(4, 9)
        states = (ExonConfiguration((), ()), ExonConfiguration((left,), ()),
                  ExonConfiguration((right,), ()))
        space = _space(states)
        features, _ = _features(space)
        marginals = {feature["feature_id"]:
                     {node: np.array([2 / 3, 1 / 3]) for node in _tree().preorder()}
                     for feature in features}
        illegal = _illegal_node_mass(space, _tree(), features, marginals)
        self.assertAlmostEqual(illegal["R"], 1 / 9, places=12)

    def test_illegal_mass_includes_material_and_exon_joint_constraints(self):
        span = ExonSpan(4, 7)
        material = Material("m", 3, 8)
        states = (ExonConfiguration((), (0,)), ExonConfiguration((), (1,)),
                  ExonConfiguration((span,), (1,)), ExonConfiguration((), (2,)))
        space = _space(states, materials=(material,))
        features, _ = _features(space)
        marginals = {feature["feature_id"]:
                     {node: np.array([0.5, 0.5]) for node in _tree().preorder()}
                     for feature in features}
        illegal = _illegal_node_mass(space, _tree(), features, marginals)
        self.assertAlmostEqual(illegal["R"], 0.25, places=12)

    def test_node_marginal_nan_is_rejected(self):
        span = ExonSpan(2, 5)
        space = _space((ExonConfiguration((), ()), ExonConfiguration((span,), ())))
        features, _ = _features(space)
        marginals = {features[0]["feature_id"]:
                     {node: np.array([0.5, 0.5]) for node in _tree().preorder()}}
        marginals[features[0]["feature_id"]]["R"] = np.array([math.nan, 0.5])
        with self.assertRaisesRegex(ValueError, "non-finite"):
            _illegal_node_mass(space, _tree(), features, marginals)

    def test_endpoint_scoring_and_fit_has_no_truth_input(self):
        self.assertNotIn("truth", inspect.signature(fit_independent_exon_characters).parameters)
        fit = {"status": "estimated", "features": [
            {"feature_id": "exon:2:5", "kind": "exon_span"},
            {"feature_id": "material:m", "kind": "material_presence"}],
            "branch_posteriors": [
                {"parent_node": "R", "child_node": "A", "features": [
                    {"feature_id": "exon:2:5", "kind": "exon_span",
                     "p00": 0.4, "p01": 0.3, "p10": 0.2, "p11": 0.1},
                    {"feature_id": "material:m", "kind": "material_presence",
                     "p00": 0.5, "p01": 0.2, "p10": 0.2, "p11": 0.1}]},
                {"parent_node": "R", "child_node": "B", "features": [
                    {"feature_id": "exon:2:5", "kind": "exon_span",
                     "p00": 0.5, "p01": 0.1, "p10": 0.1, "p11": 0.3},
                    {"feature_id": "material:m", "kind": "material_presence",
                     "p00": 0.4, "p01": 0.2, "p10": 0.2, "p11": 0.2}]}]}
        truth = {("R", "A"): {"exon:2:5": (0, 1)},
                 ("R", "B"): {"exon:2:5": (0, 0), "material:m": (0, 0)}}
        truth[("R", "A")]["material:m"] = (1, 1)
        score = score_independent_endpoint_posteriors(fit, truth)
        self.assertEqual(score["status"], "scored")
        self.assertEqual(len(score["edges"]), 2)
        self.assertAlmostEqual(score["edges"][0]["p_any_exon_geometry_endpoint_change"], 0.5)
        self.assertAlmostEqual(score["edges"][0]["true_full_endpoint_log_probability"],
                               math.log(0.3) + math.log(0.1))
        self.assertAlmostEqual(score["edges"][1]["p_any_exon_geometry_endpoint_change"], 0.2)

    def test_zero_endpoint_probability_scores_as_negative_infinity(self):
        fit = {"status": "estimated", "features": [
            {"feature_id": "exon:2:5", "kind": "exon_span"}],
            "branch_posteriors": [{"parent_node": "R", "child_node": "A", "features": [
                {"feature_id": "exon:2:5", "kind": "exon_span",
                 "p00": 0.0, "p01": 0.5, "p10": 0.25, "p11": 0.25}]}]}
        truth = {("R", "A"): {"exon:2:5": (0, 0)}}
        scored = score_independent_endpoint_posteriors(fit, truth)
        self.assertEqual(scored["edges"][0]["true_full_endpoint_log_probability"], -math.inf)

    def test_endpoint_scorer_rejects_nan_and_negative_probabilities(self):
        base = {"status": "estimated", "features": [
            {"feature_id": "exon:2:5", "kind": "exon_span"}],
            "branch_posteriors": [{"parent_node": "R", "child_node": "A", "features": [
                {"feature_id": "exon:2:5", "kind": "exon_span",
                 "p00": 0.25, "p01": 0.25, "p10": 0.25, "p11": 0.25}]}]}
        truth = {("R", "A"): {"exon:2:5": (0, 0)}}
        for bad in (math.nan, -0.01):
            fit = {**base, "branch_posteriors": [
                {**base["branch_posteriors"][0], "features": [
                    {**base["branch_posteriors"][0]["features"][0], "p00": bad}]}]}
            with self.assertRaisesRegex(ValueError, "probability"):
                score_independent_endpoint_posteriors(fit, truth)


if __name__ == "__main__":
    unittest.main()
