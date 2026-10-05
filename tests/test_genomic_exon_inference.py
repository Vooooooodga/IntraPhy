"""Regression contracts for the genomic-span fit and existing CTMC kernel."""
from types import SimpleNamespace
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from intraphy.inference.configuration_model import RateModel, evaluate_model, generator
from intraphy.inference.genomic_exon_rates import fit_family_rate, fixed_family_fit
from intraphy.inference.genomic_exon_run import infer_genomic_exons
from intraphy.structure.edits import EDIT_KINDS
from intraphy.structure.space import enumerate_space
from intraphy.structure.types import Catalogue, ExonConfiguration, ExonSpan, Material, ObservationEvidence
from intraphy.structure.serialization import write_catalogues
from intraphy.topology import SpeciesTree


class _Tree:
    def __init__(self, length):
        self.length = length

    def edges(self):
        yield "root", "tip"

    def branch_length(self, child):
        return self.length


def _unit(length, pattern):
    return {"space": object(), "tree": _Tree(length), "tips": {"pattern": pattern}}


def _analytic_log_likelihood(space, tree, tips, model, **kwargs):
    exposure = model.rates[EDIT_KINDS[0]] * tree.branch_length("tip")
    pattern = tips["pattern"]
    if pattern == "constant":
        return {"log_likelihood": 0.}
    if pattern == "concordant":
        return {"log_likelihood": -exposure}
    if pattern == "discordant":
        return {"log_likelihood": -np.exp(-exposure)}
    if pattern == "finite":
        return {"log_likelihood": -(exposure-2.)**2}
    raise AssertionError(pattern)


class GenomicExonFitTests(unittest.TestCase):
    def test_state_space_only_enumerates_complete_units_without_inference_or_state_rows(self):
        span = ExonSpan(0, 2)
        qualified = Catalogue("f", "u", 15,
            (span, ExonSpan(3, 5), ExonSpan(6, 8), ExonSpan(9, 11),
             ExonSpan(12, 14)), (),
            observations=(ObservationEvidence("A", (ExonConfiguration(()),)),
                          ObservationEvidence("B", (ExonConfiguration((span,)),)),),
            observation_unit="genomic_exon_spans")
        unqualified = Catalogue("g", "v", 2, (), (), status="unresolved",
            reasons=("no_annotated_exons",), observation_unit="genomic_exon_spans")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "input"
            root.mkdir()
            (root / "species_tree.tsv").write_text(
                "node_id\tparent_id\tlabel\tbranch_length\n"
                "root\t\troot\t\nA\troot\tA\t0.1\nB\troot\tB\t0.2\n",
                encoding="utf-8")
            catalogue_path = Path(directory) / "catalogues.jsonl"
            write_catalogues(catalogue_path, (qualified, unqualified))
            output = Path(directory) / "output"
            with patch("intraphy.inference.genomic_exon_family.fit_family_rate",
                       side_effect=AssertionError("rate fit called")), \
                    patch("intraphy.inference.genomic_exon_family.fixed_family_fit",
                          side_effect=AssertionError("fixed fit called")), \
                    patch("intraphy.inference.genomic_exon_family.evaluate_model",
                          side_effect=AssertionError("CTMC called")), \
                    patch("intraphy.inference.genomic_exon_family.evaluate_compact",
                          side_effect=AssertionError("compact CTMC called")), \
                    patch("intraphy.inference.genomic_exon_run.state_rows",
                          side_effect=AssertionError("state rows materialized")), \
                    patch("intraphy.inference.genomic_exon_family.state_rows",
                          side_effect=AssertionError("state rows materialized")):
                diagnostics = infer_genomic_exons(root, output,
                    configurations=catalogue_path, state_space_only=True)
            qualified_result = next(row for row in diagnostics if row["family_id"] == "f")
            self.assertTrue(qualified_result["state_space_complete"])
            self.assertGreater(qualified_result["state_count"], 64)
            self.assertEqual(qualified_result["status"], "eligible_conditional_unit")
            self.assertEqual(qualified_result["probability_status"], "not_computed")
            unqualified_result = next(row for row in diagnostics if row["family_id"] == "g")
            self.assertEqual(unqualified_result["status"], "unqualified_catalogue")
            self.assertIsNone(unqualified_result["state_count"])
            manifest = json.loads((output / "run_result.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["status"], "state_space_only")
            self.assertEqual(set(manifest["artifacts"]), {
                "exon_configurations.jsonl", "species_tree.tsv", "state_space_diagnostics.jsonl"})
            self.assertFalse((output / "exon_structure_fit.json").exists())
            self.assertFalse((output / "exon_history.json").exists())
            progress = [json.loads(line) for line in
                        (output / "state_space_diagnostics.jsonl").read_text(encoding="utf-8").splitlines()]
            self.assertEqual([row["stage"] for row in progress], [
                "enumeration_started", "enumeration_and_eligibility_complete",
                "unit_started", "unit_assessment_complete"])

    def test_concordant_one_edge_optimizer_harness_selects_zero_boundary(self):
        with patch("intraphy.inference.genomic_exon_rates.evaluate_model", _analytic_log_likelihood):
            fit = fit_family_rate([_unit(1., "concordant")])
        self.assertEqual(fit["status"], "zero_boundary")
        self.assertEqual(fit["mu"], 0.)

    def test_rate_fit_uses_sparse_likelihood_backend(self):
        calls = []
        def likelihood(*args, **kwargs):
            calls.append(kwargs.copy())
            return _analytic_log_likelihood(*args, **kwargs)
        with patch("intraphy.inference.genomic_exon_rates.evaluate_model", likelihood):
            fixed_family_fit([_unit(1., "concordant")],
                             RateModel({kind: .1 for kind in EDIT_KINDS}))
        self.assertTrue(calls)
        self.assertTrue(all(call["backend"] == "sparse" for call in calls))

    def test_discordant_single_unit_has_no_finite_rate_identification(self):
        with patch("intraphy.inference.genomic_exon_rates.evaluate_model", _analytic_log_likelihood):
            fit = fit_family_rate([_unit(1., "discordant")])
        self.assertEqual(fit["status"], "upper_tail_unresolved")
        self.assertIsNone(fit.get("mu"))

    def test_one_variable_and_constant_units_identify_finite_positive_rate(self):
        units = [_unit(1., "finite"), _unit(1., "constant"), _unit(1., "constant")]
        with patch("intraphy.inference.genomic_exon_rates.evaluate_model", _analytic_log_likelihood):
            fit = fit_family_rate(units)
        self.assertEqual(fit["status"], "estimated_conditional_composite_rate")
        self.assertAlmostEqual(fit["mu"], 2., delta=.02)

    def test_proportional_tree_length_scaling_preserves_likelihood_and_scales_rate(self):
        units_a = [_unit(1., "finite"), _unit(1., "constant"), _unit(1., "constant")]
        units_b = [_unit(4., "finite"), _unit(4., "constant"), _unit(4., "constant")]
        with patch("intraphy.inference.genomic_exon_rates.evaluate_model", _analytic_log_likelihood):
            first, second = fit_family_rate(units_a), fit_family_rate(units_b)
        self.assertAlmostEqual(first["log_likelihood"], second["log_likelihood"], places=5)
        self.assertAlmostEqual(first["mu"], 4*second["mu"], delta=.02)

    def test_fixed_family_fit_uses_existing_configuration_ctmc_kernel(self):
        catalogue = Catalogue("f", "u", 4, (ExonSpan(0, 2),), (),
            observations=(ObservationEvidence("A", (ExonConfiguration((ExonSpan(0, 2),)),)),
                          ObservationEvidence("B", (ExonConfiguration(()),))),
            boundary_candidates=(ExonSpan(0, 2),), observation_unit="genomic_exon_spans")
        space = enumerate_space(catalogue)
        tree = SpeciesTree([{"node_id": "r", "parent_id": "", "label": "r"},
            {"node_id": "a", "parent_id": "r", "label": "A", "branch_length": .2},
            {"node_id": "b", "parent_id": "r", "label": "B", "branch_length": .3}])
        tips = {"A": np.array([0., 1.]), "B": np.array([1., 0.])}
        model = RateModel({kind: .04 for kind in EDIT_KINDS}, scale=1.)
        expected = evaluate_model(space, tree, tips, model, posterior=False, counts=False)["log_likelihood"]
        fitted = fixed_family_fit([{"space": space, "tree": tree, "tips": tips}], model)
        self.assertAlmostEqual(fitted["log_likelihood"], expected)

    def test_real_two_state_kernel_fits_one_variable_with_three_constant_units(self):
        exon = ExonSpan(0, 2)
        catalogue = Catalogue("f", "u", 2, (exon,), (), boundary_candidates=(exon,),
                              observation_unit="genomic_exon_spans")
        space = enumerate_space(catalogue)
        self.assertEqual(len(space.states), 2)
        tree = SpeciesTree([{"node_id": "r", "parent_id": "", "label": "r"},
            {"node_id": "a", "parent_id": "r", "label": "A", "branch_length": .5},
            {"node_id": "b", "parent_id": "r", "label": "B", "branch_length": .5}])
        present = np.array([0., 1.])
        absent = np.array([1., 0.])
        same = {"A": present, "B": present}
        different = {"A": present, "B": absent}
        units = [{"unit_id": str(i), "space": space, "tree": tree,
                  "tips": same if i < 3 else different} for i in range(4)]
        fit = fit_family_rate(units)
        self.assertEqual(fit["status"], "estimated_conditional_composite_rate")
        self.assertAlmostEqual(fit["mu"], np.log(2.)/2., delta=.02)

    def test_real_two_state_same_tip_states_have_zero_boundary_fit(self):
        exon = ExonSpan(0, 2)
        space = enumerate_space(Catalogue("f", "u", 2, (exon,), (),
            boundary_candidates=(exon,), observation_unit="genomic_exon_spans"))
        tree = SpeciesTree([{"node_id": "r", "parent_id": "", "label": "r"},
            {"node_id": "a", "parent_id": "r", "label": "A", "branch_length": .5},
            {"node_id": "b", "parent_id": "r", "label": "B", "branch_length": .5}])
        present = np.array([0., 1.])
        fit = fit_family_rate([{"unit_id": "u", "space": space, "tree": tree,
                                "tips": {"A": present, "B": present}}])
        self.assertEqual(fit["status"], "zero_boundary")
        self.assertEqual(fit["mu"], 0.)

    def test_real_kernel_rate_is_inverse_to_proportional_tree_length(self):
        exon = ExonSpan(0, 2)
        space = enumerate_space(Catalogue("f", "u", 2, (exon,), (),
            boundary_candidates=(exon,), observation_unit="genomic_exon_spans"))
        present, absent = np.array([0., 1.]), np.array([1., 0.])
        results = []
        for branch in (.5, 2.):
            tree = SpeciesTree([{"node_id": "r", "parent_id": "", "label": "r"},
                {"node_id": "a", "parent_id": "r", "label": "A", "branch_length": branch},
                {"node_id": "b", "parent_id": "r", "label": "B", "branch_length": branch}])
            same, different = {"A": present, "B": present}, {"A": present, "B": absent}
            units = [{"unit_id": str(i), "space": space, "tree": tree,
                      "tips": same if i < 3 else different} for i in range(4)]
            results.append(fit_family_rate(units))
        self.assertAlmostEqual(results[0]["log_likelihood"], results[1]["log_likelihood"], places=5)
        self.assertAlmostEqual(results[0]["mu"], 4*results[1]["mu"], delta=.02)

    def test_shared_dna_deletion_uses_original_one_way_kernel(self):
        left, right = ExonSpan(0, 4), ExonSpan(6, 10)
        catalogue = Catalogue("f", "u", 10, (left, right), ((4, 6),),
            (Material("tract", 4, 6),), boundary_candidates=(left, right),
            observation_unit="genomic_exon_spans")
        space = enumerate_space(catalogue)
        source = ExonConfiguration((left, right), (1,))
        target = ExonConfiguration((ExonSpan(0, 10),), (2,))
        model = RateModel({kind: .2 for kind in EDIT_KINDS})
        q, _ = generator(space, model, {"tract": "r"}, "tip")
        edit = next(edit for edit in space.edits
                    if edit.source == source and edit.target == target and edit.kind == "dna_deletion")
        self.assertAlmostEqual(q[space.index[source], space.index[target]], .2*edit.weight)
        self.assertEqual(q[space.index[target], space.index[source]], 0.)


if __name__ == "__main__":
    unittest.main()
