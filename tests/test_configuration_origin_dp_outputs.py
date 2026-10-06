"""Fit agreement and retained posterior-route outputs for origin-subset DP."""
import unittest
from unittest.mock import patch

import numpy as np

from intraphy.inference.configuration_compact import evaluate_compact
from intraphy.inference.configuration_model import RateModel, evaluate_model
from intraphy.inference.genomic_exon_rates import fit_family_rate
from intraphy.structure.edits import EDIT_KINDS
from intraphy.structure.space import enumerate_space
from intraphy.structure.types import Catalogue, Material
from intraphy.topology import SpeciesTree


def _tree():
    return SpeciesTree([
        {"node_id": "r", "parent_id": "", "label": "r"},
        {"node_id": "a", "parent_id": "r", "label": "A", "branch_length": 0.35},
        {"node_id": "b", "parent_id": "r", "label": "B", "branch_length": 0.65},
    ])


def _space():
    material = Material("m0", 2, 3)
    return enumerate_space(Catalogue("f", "fit", 4, (), (),
        material=(material,), boundary_candidates=()))


def _origin_weights(result):
    return {tuple(sorted(row["origins"].items())): row["posterior_weight"]
            for row in result["origins"]}


class OriginSubsetOutputTests(unittest.TestCase):
    def test_fit_agrees_with_dense_full_enumeration_and_keeps_posterior_outputs(self):
        state_space, species_tree = _space(), _tree()
        present = np.asarray([float(state.material == (1,))
                              for state in state_space.states])
        absent = 1. - present
        units = [{"space": state_space, "tree": species_tree,
                  "tips": {"A": present, "B": present}} for _ in range(3)]
        units.append({"space": state_space, "tree": species_tree,
                      "tips": {"A": present, "B": absent}})

        sparse_fit = fit_family_rate(units)
        original = evaluate_model

        def dense_reference(*args, **kwargs):
            kwargs.pop("sparse_templates", None)
            kwargs["backend"] = "dense"
            return original(*args, **kwargs)

        with patch("intraphy.inference.genomic_exon_rates.evaluate_model",
                   side_effect=dense_reference):
            dense_fit = fit_family_rate(units)

        self.assertEqual(sparse_fit["status"], dense_fit["status"])
        self.assertEqual(sparse_fit["converged"], dense_fit["converged"])
        self.assertIsNotNone(sparse_fit["mu"])
        self.assertAlmostEqual(sparse_fit["mu"], dense_fit["mu"], delta=2e-7)
        self.assertAlmostEqual(sparse_fit["log_likelihood"],
                               dense_fit["log_likelihood"], delta=2e-9)
        self.assertEqual([x for x, _ in sparse_fit["starting_log_likelihoods"]],
                         [x for x, _ in dense_fit["starting_log_likelihoods"]])

        fitted_posteriors = []
        for fitted_mu in (sparse_fit["mu"], dense_fit["mu"]):
            fitted_model = RateModel({kind: fitted_mu for kind in EDIT_KINDS})
            dense = evaluate_model(state_space, species_tree,
                                   units[0]["tips"], fitted_model)
            compact = evaluate_compact(state_space, species_tree,
                                       units[0]["tips"], fitted_model)
            self._assert_posterior_routes_agree(dense, compact)
            fitted_posteriors.append(dense)

        sparse_mu_posterior, dense_mu_posterior = fitted_posteriors
        sparse_weights = _origin_weights(sparse_mu_posterior)
        dense_weights = _origin_weights(dense_mu_posterior)
        self.assertEqual(set(sparse_weights), set(dense_weights))
        for scenario in sparse_weights:
            self.assertAlmostEqual(sparse_weights[scenario], dense_weights[scenario],
                                   delta=2e-6)
        self.assertEqual(len(sparse_mu_posterior["branches"]),
                         len(dense_mu_posterior["branches"]))
        for sparse_branch, dense_branch in zip(sparse_mu_posterior["branches"],
                                               dense_mu_posterior["branches"]):
            self.assertEqual((sparse_branch["parent"], sparse_branch["child"]),
                             (dense_branch["parent"], dense_branch["child"]))
            np.testing.assert_allclose(sparse_branch["endpoint_probabilities"],
                                       dense_branch["endpoint_probabilities"],
                                       atol=2e-6, rtol=0.)

    def test_dense_and_compact_retain_joint_origin_posterior(self):
        state_space, species_tree = _space(), _tree()
        present = np.asarray([float(state.material == (1,))
                              for state in state_space.states])
        tips = {"A": present, "B": np.ones(len(state_space.states))}
        fit_model = RateModel({kind: 0.04 for kind in EDIT_KINDS},
                              origin_root_weight=2.)
        dense = evaluate_model(state_space, species_tree, tips, fit_model)
        compact = evaluate_compact(state_space, species_tree, tips, fit_model)
        self._assert_posterior_routes_agree(dense, compact)

    def _assert_posterior_routes_agree(self, dense, compact):
        self.assertAlmostEqual(dense["log_likelihood"], compact["log_likelihood"],
                               delta=2e-10)
        for result in (dense, compact):
            self.assertTrue(result["origins"])
            self.assertAlmostEqual(sum(row["posterior_weight"]
                                       for row in result["origins"]), 1., places=10)
            self.assertTrue(result["branches"])
        self.assertTrue(all("endpoint_probabilities" in branch
                            for branch in dense["branches"]))
        self.assertTrue(all("joint_modes" in branch["observable_endpoint_summary"]
                            for branch in compact["branches"]))
        dense_weights, compact_weights = _origin_weights(dense), _origin_weights(compact)
        self.assertEqual(set(dense_weights), set(compact_weights))
        for scenario in dense_weights:
            self.assertAlmostEqual(dense_weights[scenario], compact_weights[scenario],
                                   delta=2e-6)


if __name__ == "__main__":
    unittest.main()
