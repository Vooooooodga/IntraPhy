from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

from intraphy.inference.locus_bootstrap import (
    _simulated_units, bootstrap_locus_comparison, validate_dna_bootstrap_design,
)
from intraphy.inference.locus_rates import LocusFitUnit
from intraphy.structure.locus_sampling import sample_history, sample_tip_observation
from intraphy.structure.locus_types import (
    CopySlot, EventOpportunity, LocusCatalogue, LocusProcess, LocusState,
    MaterialTract, ProcessEdge, SpliceFeature,
)
from intraphy.topology import SpeciesTree


def _fixture(material=(0,), features=None, sensitivity=None, specificity=None,
             latent_features=()):
    feature_rows = [
        {"id": f.id, "kind": f.kind, "required_material": list(f.required_material),
         "prerequisites": list(f.prerequisites), "copy_id": f.copy_id,
         "start": f.start, "end": f.end, "donor": f.donor, "acceptor": f.acceptor,
         "mxe_group": f.mxe_group, "evidence": list(f.evidence)}
        for f in latent_features
    ]
    catalogue = LocusCatalogue(
        (MaterialTract("A", 0, 10),), (CopySlot("cA", ("A",)),), tuple(latent_features),
        (EventOpportunity("del", "delA", "dna_deletion", "loss",
                          material_deletions=("A",), interval=(0, 10)),),
    )
    present, absent = LocusState((1,), frozenset()), LocusState((2,), frozenset())
    process = LocusProcess(catalogue, (present, absent),
        (ProcessEdge(0, 1, "del", "delA", "dna_deletion", "loss", 1.0),))
    tree = SpeciesTree([
        {"node_id": "root", "parent_id": "", "label": "root"},
        {"node_id": "tip", "parent_id": "root", "label": "tip", "branch_length": 1.0},
    ])
    tip = {"material": list(material), "features": features or {},
           "surveyed_features": [], "sensitivity": {}, "specificity": {},
           "observed_paths": [], "surveyed_material": ["A"],
           "material_sensitivity": sensitivity or {}, "material_specificity": specificity or {},
           "material_evidence": {"A": "fixed test evidence"}}
    raw_unit = {"family": "f", "unit": "u",
                "catalogue": {"material": [{"id": "A"}], "features": feature_rows},
                "observations": {"tips": {"tip": tip}}}
    unit = SimpleNamespace(family="f", unit="u", process=process, tree=tree,
                           root_prior=np.array([1.0, 0.0]), tips={"tip": np.array([1.0, 0.0])})
    model_record = {"units": [raw_unit]}
    bundle = SimpleNamespace(units=(unit,), model_record=model_record,
        rates={"loss": {"mode": "fit", "value": 0.7}},
        fit_units=lambda: (LocusFitUnit(process, tree, {"tip": np.array([1.0, 0.0])},
                                        np.array([1.0, 0.0]), name="f/u"),))
    return bundle, process, tree, tip


class FixedRng:
    def choice(self, _size, p):
        return int(np.flatnonzero(np.asarray(p) > 0)[0])

    def exponential(self, _scale):
        return 0.25

    def random(self):
        return 0.5


class LocusBootstrapTests(unittest.TestCase):
    def test_gillespie_sampler_follows_one_explicit_edge(self):
        _bundle, process, tree, _tip = _fixture()
        states = sample_history(process, tree, np.array([1.0, 0.0]), {"loss": 1.0}, FixedRng())
        self.assertEqual(states["root"], process.states[0])
        self.assertEqual(states["tip"], process.states[1])

    def test_missing_material_mask_is_retained(self):
        _bundle, _process, _tree, tip = _fixture(material=(None,))
        observation = sample_tip_observation(
            _bundle.units[0].process.catalogue, tip,
            _bundle.units[0].process.states[0], FixedRng())
        self.assertEqual(observation.material, (None,))

    def test_material_detection_sensitivity_zero_makes_present_call_negative(self):
        bundle, _process, _tree, tip = _fixture(material=(1,), sensitivity={"A": 0.0},
                                                 specificity={"A": 1.0})
        observation = sample_tip_observation(bundle.units[0].process.catalogue, tip,
            bundle.units[0].process.states[0], FixedRng())
        self.assertEqual(observation.material, (0,))

    def test_same_seed_reproduces_simulated_material_calls(self):
        bundle, *_ = _fixture(material=(0,))
        def simulate():
            rng = np.random.default_rng(np.random.SeedSequence([83, 4]))
            return _simulated_units(bundle, {"loss": 0.7}, rng)[0].tips["tip"]
        np.testing.assert_array_equal(simulate(), simulate())

    def test_feature_metadata_is_rejected_before_comparison(self):
        bundle, *_ = _fixture(features={"exon": 1})
        with self.assertRaisesRegex(ValueError, "features"):
            bootstrap_locus_comparison(bundle, ("loss",), 2, 4,
                                       sampling_design="fixed_catalogue")

    def test_latent_catalogue_features_are_allowed_without_tip_observations(self):
        feature = SpliceFeature("exA", "exon", ("A",), copy_id="cA", start=0, end=10)
        bundle, *_ = _fixture(latent_features=(feature,))
        validate_dna_bootstrap_design(bundle)

    def test_misaligned_material_order_is_rejected(self):
        bundle, *_ = _fixture()
        bundle.model_record["units"][0]["catalogue"]["material"] = [
            {"id": "wrong"}]
        with self.assertRaisesRegex(ValueError, "material ID order"):
            validate_dna_bootstrap_design(bundle)

    def test_failed_replicate_is_retained_and_invalidates_p_value(self):
        bundle, *_ = _fixture()
        good_fit = SimpleNamespace(optimizer_success=True, rates={"loss": 0.4},
                                   fitted_log_likelihood=-1.0)
        bad_fit = SimpleNamespace(optimizer_success=False, rates={"loss": 0.4},
                                  fitted_log_likelihood=-1.1)
        good = SimpleNamespace(null_fit=good_fit, full_fit=good_fit,
                               statistic=0.2, status="evaluated_within_numerical_tolerance")
        bad = SimpleNamespace(null_fit=bad_fit, full_fit=good_fit,
                              statistic=None, status="one_or_both_fits_unsuccessful")
        with patch("intraphy.inference.locus_bootstrap.compare_locus_rates",
                   side_effect=(good, good, bad)):
            result = bootstrap_locus_comparison(bundle, ("loss",), 2, 17,
                sampling_design="fixed_catalogue")
        self.assertIsNone(result.p_value)
        self.assertEqual(result.status, "incomplete")
        self.assertEqual(result.completed_replicates, 1)
        self.assertEqual([row["status"] for row in result.replicates], ["success", "failed"])

    def test_failed_observed_comparison_is_retained_without_sampling(self):
        bundle, *_ = _fixture()
        failed_fit = SimpleNamespace(optimizer_success=False, rates={"loss": 0.7},
                                     fitted_log_likelihood=-1.0)
        failed = SimpleNamespace(null_fit=failed_fit, full_fit=failed_fit,
                                 statistic=None, status="one_or_both_fits_unsuccessful")
        with patch("intraphy.inference.locus_bootstrap.compare_locus_rates", return_value=failed):
            with patch("intraphy.inference.locus_bootstrap.sample_history") as sampler:
                result = bootstrap_locus_comparison(bundle, ("loss",), 3, 9,
                    sampling_design="fixed_catalogue")
        sampler.assert_not_called()
        self.assertIs(result.observed_comparison, failed)
        self.assertEqual(result.requested_replicates, 3)
        self.assertEqual(result.completed_replicates, 0)
        self.assertEqual(result.workers, 0)
        self.assertEqual(result.replicates, ())
        self.assertIsNone(result.p_value)


if __name__ == "__main__":
    unittest.main()
