"""Contract fixtures for pooled repertoire scale fitting; intentionally not run here."""
import unittest
from dataclasses import replace
import numpy as np
from unittest.mock import patch

from intraphy.inference.repertoire_rates import RepertoireInferenceUnit, fit_repertoire_scale
import intraphy.inference.repertoire_rates as rates_module
from intraphy.structure.repertoire_process import deletion_process
from intraphy.structure.repertoire import ExonRepertoire
from intraphy.structure.types import Catalogue, ExonConfiguration as C, ExonSpan as E, Material
from intraphy.topology import SpeciesTree


def _tree():
    return SpeciesTree([{"node_id":"r","parent_id":"","label":"r"},
                        {"node_id":"A","parent_id":"r","label":"A","branch_length":1.},
                        {"node_id":"B","parent_id":"r","label":"B","branch_length":1.}])


def _unit(family, tips=None):
    c = Catalogue(family, "u", 10, (E(0, 10),), (), material=(Material("m", 0, 10),),
                  boundary_candidates=(E(0, 10),))
    seed = ExonRepertoire((C((E(0, 10),), (1,)),), (1,))
    p = deletion_process(c, (seed,), provenance="synthetic deletion process")
    emissions = tips or {"A": np.array([1., 0.]), "B": np.array([0., 1.])}
    return RepertoireInferenceUnit(p, _tree(), emissions, np.array([1., 0.]),
                                   "root present prior", "synthetic two-tip observation")


class RepertoireRateContractTests(unittest.TestCase):
    def test_two_independent_genes_have_interior_log_two_scale(self):
        result = fit_repertoire_scale((_unit("g1"), _unit("g2")),
            rates={"dna_deletion": 1.}, rate_provenance="fixed synthetic base rate",
            branch_length_unit="unit branches", collection_provenance="independent genes",
            log_bounds=(-4., 4.))
        self.assertEqual(result["status"], "estimated_conditional_scale")
        self.assertAlmostEqual(result["scale"], np.log(2.), places=3)

    def test_single_gene_is_rejected(self):
        self.assertEqual(fit_repertoire_scale((_unit("g1"),), rates={"dna_deletion": 1.},
            rate_provenance="r", branch_length_unit="u", collection_provenance="c",
            log_bounds=(-2., 2.))["status"], "fewer_than_two_genes")

    def test_unknown_is_nonidentifiable_and_all_present_is_boundary(self):
        unknown = {"A": np.ones(2), "B": np.ones(2)}
        self.assertEqual(fit_repertoire_scale((_unit("g1", unknown), _unit("g2", unknown)),
            rates={"dna_deletion": 1.}, rate_provenance="r", branch_length_unit="u", collection_provenance="c",
            log_bounds=(-4., 4.))["status"], "nonidentifiable_curvature")
        present = {"A": np.array([1., 0.]), "B": np.array([1., 0.])}
        self.assertEqual(fit_repertoire_scale((_unit("g1", present), _unit("g2", present)),
            rates={"dna_deletion": 1.}, rate_provenance="r", branch_length_unit="u", collection_provenance="c",
            log_bounds=(-4., 4.))["status"], "at_numerical_boundary")

    def test_duplicate_collection_and_invalid_contracts_fail(self):
        kwargs = dict(rates={"dna_deletion": 1.}, rate_provenance="r", branch_length_unit="u",
                      collection_provenance="c", log_bounds=(-2., 2.))
        with self.assertRaises(ValueError): fit_repertoire_scale((_unit("g1"), _unit("g1")), **kwargs)
        with self.assertRaises(ValueError): fit_repertoire_scale((_unit("g1"), _unit("g2")), rates={"dna_deletion": -1.}, **{k:v for k,v in kwargs.items() if k != "rates"})
        with self.assertRaises(ValueError): fit_repertoire_scale((_unit("g1"), _unit("g2")), **{**kwargs, "log_bounds": (1., 1.)})

    def test_impossible_later_unit_is_not_dropped(self):
        impossible = _unit("g2", {"A": np.zeros(2), "B": np.zeros(2)})
        result = fit_repertoire_scale((_unit("g1"), impossible), rates={"dna_deletion": 1.},
            rate_provenance="r", branch_length_unit="u", collection_provenance="c", log_bounds=(-4., 4.))
        self.assertFalse(result.get("valid_fit", False))
        self.assertIn(result["status"], {"optimizer_failed_or_impossible_observation", "optimizer_partial_failure"})

    def test_overlap_provenance_and_failed_optimizer_contracts(self):
        a, b = _unit("g1"), _unit("g1")
        b = replace(b, process=replace(b.process, catalogue=replace(b.process.catalogue, unit="v", alignment_offset=5)))
        with self.assertRaises(ValueError): fit_repertoire_scale((a, b, _unit("g2")), rates={"dna_deletion": 1.}, rate_provenance="r", branch_length_unit="u", collection_provenance="c", log_bounds=(-2., 2.))
        with self.assertRaises(ValueError): fit_repertoire_scale((a, _unit("g2")), rates={"dna_deletion": 1.}, rate_provenance="", branch_length_unit="u", collection_provenance="c", log_bounds=(-2., 2.))
        class Failed:
            success = False; fun = np.inf; x = np.nan
        with patch("intraphy.inference.repertoire_rates.minimize_scalar", return_value=Failed()):
            result = fit_repertoire_scale((a, _unit("g2")), rates={"dna_deletion": 1.}, rate_provenance="r", branch_length_unit="u", collection_provenance="c", log_bounds=(-2., 2.))
        self.assertFalse(result.get("valid_fit", False))

    def test_boundary_curvature_is_none_and_narrow_bounds_are_valid_contract(self):
        present = {"A": np.array([1., 0.]), "B": np.array([1., 0.])}
        result = fit_repertoire_scale((_unit("g1", present), _unit("g2", present)), rates={"dna_deletion": 1.}, rate_provenance="r", branch_length_unit="u", collection_provenance="c", log_bounds=(-.01, .01))
        self.assertIsNone(result["curvature"])

    def test_zero_rate_unknown_data_is_nonidentifiable(self):
        unknown = {"A": np.ones(2), "B": np.ones(2)}
        result = fit_repertoire_scale((_unit("g1", unknown), _unit("g2", unknown)), rates={"dna_deletion": 0.}, rate_provenance="r", branch_length_unit="u", collection_provenance="c", log_bounds=(-2., 2.))
        self.assertEqual(result["status"], "nonidentifiable_curvature")

    def test_zero_rate_validation_and_impossible_observation(self):
        unknown = {"A": np.ones(2), "B": np.ones(2)}
        with self.assertRaises(ValueError):
            fit_repertoire_scale((replace(_unit("g1", unknown), root_prior=np.array([1., 1.])), _unit("g2", unknown)), rates={"dna_deletion": 0.}, rate_provenance="r", branch_length_unit="u", collection_provenance="c", log_bounds=(-2., 2.))
        impossible = fit_repertoire_scale((_unit("g1"), _unit("g2")), rates={"dna_deletion": 0.}, rate_provenance="r", branch_length_unit="u", collection_provenance="c", log_bounds=(-2., 2.))
        self.assertFalse(impossible["valid_fit"])
        malformed = _unit("g1", {"A": np.array([-1., 0.]), "B": np.array([1., 1.])})
        with self.assertRaises(ValueError):
            fit_repertoire_scale((malformed, _unit("g2", unknown)), rates={"dna_deletion": 0.}, rate_provenance="r", branch_length_unit="u", collection_provenance="c", log_bounds=(-2., 2.))

    def test_nonfinite_midpoint_does_not_abort_search(self):
        calls = []
        def fake(process, tree, tips, **kwargs):
            calls.append(kwargs["scale"])
            x = np.log(kwargs["scale"])
            return {"log_likelihood": -np.inf if abs(x) < 1e-12 else -(x - .2) ** 2}
        with patch.object(rates_module, "evaluate_repertoire_model", side_effect=fake):
            result = fit_repertoire_scale((_unit("g1"), _unit("g2")), rates={"dna_deletion": 1.}, rate_provenance="r", branch_length_unit="u", collection_provenance="c", log_bounds=(-2., 2.))
        self.assertTrue(calls)
        self.assertGreater(len(calls), 2)
        self.assertTrue(np.isfinite(result["candidate_scale"]))
        self.assertGreater(result["log_likelihood"], -0.1)


if __name__ == "__main__":
    unittest.main()
