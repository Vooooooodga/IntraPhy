import json
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from unittest import mock

import numpy as np

from intraphy.commands.parser import build_parser
from intraphy.inference import configuration_compact
from intraphy.inference import genomic_exon_calibration as calibration
from intraphy.inference.configuration_model import RateModel, evaluate_model
from intraphy.inference.exon_rates import InferenceUnit
from intraphy.inference.exon_resampling import sample_unit_history, simulate_unit
from intraphy.inference.genomic_exon_rates import fit_family_rate
from intraphy.inference.genomic_exon_branches import branch_change_rows
from intraphy.inference.genomic_exon_output import state_rows
from intraphy.verification.genomic_exon_calibration import score_replicate, summarize_arm
from intraphy.structure.edits import EDIT_KINDS
from intraphy.structure.space import enumerate_space
from intraphy.structure.types import ExonConfiguration


class GenomicExonCalibrationTests(unittest.TestCase):
    @staticmethod
    def template(scenario="geometry"):
        catalogue = calibration._catalogues(scenario)[0]
        space = enumerate_space(catalogue)
        return calibration._unit(catalogue, space, catalogue.unit)

    def test_sampler_history_preserves_legacy_api_and_seeded_draws(self):
        for scenario in ("geometry", "shared-deletion"):
            with self.subTest(scenario=scenario):
                unit = self.template(scenario)
                model = RateModel({kind: .2 for kind in EDIT_KINDS})
                legacy = simulate_unit(unit, model, np.random.default_rng(73))
                sampled, assigned, origins = sample_unit_history(
                    unit, model, np.random.default_rng(73))
                self.assertIsInstance(legacy, InferenceUnit)
                self.assertIsInstance(sampled, InferenceUnit)
                self.assertEqual(legacy.tips.keys(), sampled.tips.keys())
                for taxon in legacy.tips:
                    np.testing.assert_array_equal(legacy.tips[taxon], sampled.tips[taxon])
                self.assertEqual(set(assigned), set(unit.tree.parent))
                self.assertEqual(set(origins), {item.id for item in unit.space.catalogue.material})

    def test_zero_and_deleted_material_states_are_observationally_identical(self):
        unit = self.template("shared-deletion")
        candidates = {}
        for state in unit.space.states:
            if state.material:
                candidates.setdefault(state.exons, {})[state.material[0]] = state
        pair = next(states for states in candidates.values() if 0 in states and 2 in states)
        state0, state2 = pair[0], pair[2]
        self.assertEqual(calibration._visible(state0), calibration._visible(state2))
        config0 = ExonConfiguration(state0.exons, state0.material)
        config2 = ExonConfiguration(state2.exons, state2.material)
        np.testing.assert_array_equal(calibration._exact_tip(unit.space, config0),
                                      calibration._exact_tip(unit.space, config2))

    def test_fixed_catalogues_declare_genomic_span_observation(self):
        for scenario in ("geometry", "shared-deletion"):
            with self.subTest(scenario=scenario):
                catalogues = calibration._catalogues(scenario)
                self.assertEqual(len(catalogues), 2)
                self.assertTrue(all(item.discovery == "independent_catalogue" for item in catalogues))
                self.assertTrue(all(item.observation_unit == "genomic_exon_spans" for item in catalogues))
                self.assertTrue(all(enumerate_space(item).complete for item in catalogues))
        geometry = calibration._catalogues("geometry")[0]
        self.assertEqual(len(geometry.spans), 3)
        self.assertTrue(geometry.junctions)
        deletion = calibration._catalogues("shared-deletion")[0]
        self.assertEqual(len(deletion.material), 1)
        self.assertEqual((deletion.material[0].start, deletion.material[0].end), (0, deletion.length))

    def test_scoring_uses_compact_posterior_and_matches_dense_reference(self):
        model = RateModel({kind: .2 for kind in EDIT_KINDS})
        for index, scenario in enumerate(("geometry", "shared-deletion")):
            with self.subTest(scenario=scenario):
                unit = self.template(scenario)
                observed, assigned, origins = sample_unit_history(
                    unit, model, np.random.default_rng(311 + index))
                latent = [(observed, assigned, origins)]
                with mock.patch.object(configuration_compact, "evaluate_compact",
                                       wraps=configuration_compact.evaluate_compact) as compact:
                    scored = score_replicate("fixed", index + 1, [unit], latent, model)
                self.assertEqual(compact.call_count, 1)

                dense = evaluate_model(observed.space, observed.tree, observed.tips,
                                       model, posterior=True, counts=False)
                dense_result = {"family_id": observed.family, "unit_id": observed.unit,
                                "states": state_rows(observed.space), "ctmc": dense}
                reference_rows = branch_change_rows(
                    [{"family_id": observed.family, "units": [dense_result]}],
                    {}, observed.tree)
                reference = {}
                for row in reference_rows:
                    key = (row["parent_node_id"], row["child_node_id"])
                    reference.setdefault(key, {"probability": row["probability_exon_structure_change"],
                        "modes": set()})["modes"].add((
                            tuple(tuple(span) for span in row["parent_exons"]),
                            tuple(row["parent_dna_presence"]),
                            tuple(tuple(span) for span in row["child_exons"]),
                            tuple(row["child_dna_presence"])))

                self.assertEqual(len(scored), len(reference))
                for row in scored:
                    key = (row["parent_node_id"], row["child_node_id"])
                    self.assertAlmostEqual(row["probability_change"],
                                           reference[key]["probability"], places=12)
                    self.assertEqual(set(row["modal_pairs"]), reference[key]["modes"])

    def test_zero_rate_all_unknown_tips_preserve_dense_modal_ties(self):
        unit = self.template("geometry")
        truth_model = RateModel({kind: .2 for kind in EDIT_KINDS})
        zero_model = RateModel({kind: 0. for kind in EDIT_KINDS})
        sampled, assigned, origins = sample_unit_history(
            unit, truth_model, np.random.default_rng(417))
        observed = replace(sampled, tips={
            species: np.ones(len(sampled.space.states)) for species in sampled.tips
        })

        with mock.patch.object(configuration_compact, "evaluate_compact",
                               wraps=configuration_compact.evaluate_compact) as compact:
            scored = score_replicate(
                "fixed", 1, [unit], [(observed, assigned, origins)], zero_model)
        self.assertEqual(compact.call_count, 1)

        dense = evaluate_model(observed.space, observed.tree, observed.tips,
                               zero_model, posterior=True, counts=False)
        dense_result = {"family_id": observed.family, "unit_id": observed.unit,
                        "states": state_rows(observed.space), "ctmc": dense}
        reference_rows = branch_change_rows(
            [{"family_id": observed.family, "units": [dense_result]}],
            {}, observed.tree)
        reference = {}
        for row in reference_rows:
            key = (row["parent_node_id"], row["child_node_id"])
            reference.setdefault(key, {"probability": row["probability_exon_structure_change"],
                "modes": set()})["modes"].add((
                    tuple(tuple(span) for span in row["parent_exons"]),
                    tuple(row["parent_dna_presence"]),
                    tuple(tuple(span) for span in row["child_exons"]),
                    tuple(row["child_dna_presence"])))

        self.assertEqual(len(scored), len(reference))
        for row in scored:
            key = (row["parent_node_id"], row["child_node_id"])
            self.assertAlmostEqual(row["probability_change"],
                                   reference[key]["probability"], places=12)
            self.assertEqual(set(row["modal_pairs"]), reference[key]["modes"])
        self.assertTrue(any(row["modal_tie_count"] > 1 for row in scored))

    def test_fit_receives_only_observed_units_and_outputs_scores(self):
        captured = []

        def fitted(units):
            captured.extend(units)
            return {"converged": True, "mu": .2,
                    "status": "estimated_conditional_composite_rate",
                    "valid_for_resampling": True}

        with tempfile.TemporaryDirectory() as temporary, \
                mock.patch.object(calibration, "fit_family_rate", side_effect=fitted):
            result = calibration.calibrate_genomic_exons(temporary, "geometry", .2, 2, 1, 19)
            self.assertEqual(len(captured), 2)
            self.assertTrue(all(set(unit) == {"family", "unit", "space", "tree", "tips"}
                                for unit in captured))
            self.assertTrue(all(unit["unit"].endswith(("draw001", "draw002")) for unit in captured))
            self.assertEqual(result["estimated"]["scored_replicates"], 1)
            self.assertIsNotNone(result["estimated"]["brier"])
            self.assertFalse(result["estimated"]["discovery_pipeline_calibrated"])
            self.assertEqual(result["estimated"]["scored_branch_rows"], 12)
            for filename in ("branches.tsv", "reliability_fixed.tsv", "reliability_estimated.tsv"):
                self.assertTrue((Path(temporary) / filename).exists())
            metadata = json.loads((Path(temporary) / "calibration_metadata.json").read_text())
            self.assertFalse(metadata["discovery_pipeline_calibrated"])
            self.assertEqual(metadata["posterior_backend"], "compact_origin_subset")
            self.assertEqual(metadata["generating_model"]["rates"]["split"], .2)
            self.assertIn("each other canonical node", metadata["material_origin_sampling"])

    def test_all_empty_shared_deletion_tips_withhold_rate_at_finite_tail_probe(self):
        units = []
        for index in range(8):
            catalogue = calibration._catalogues("shared-deletion")[index % 2]
            space = enumerate_space(catalogue)
            template = calibration._unit(catalogue, space, catalogue.unit)
            empty = np.asarray([float(calibration._visible(state) == ((), (0,)))
                                for state in space.states])
            self.assertTrue(np.any(empty))
            units.append({"space": space, "tree": template.tree,
                "tips": {species: empty.copy() for species in template.tips}})

        fit = fit_family_rate(units)
        self.assertEqual(fit["status"], "upper_tail_unresolved")
        self.assertFalse(fit["converged"])
        self.assertIsNone(fit["mu"])
        diagnostic = fit["upper_tail_diagnostic"]
        self.assertEqual(diagnostic["x_probe"],
                         2. * max(fit["selected_candidate"]["x"], 100.))
        self.assertIsNotNone(diagnostic["log_likelihood"])
        self.assertIn("finite_probe_only", diagnostic["scope"])

    def test_finite_high_rate_peak_passes_tail_probe(self):
        unit = self.template()
        fit_unit = {"space": unit.space, "tree": unit.tree, "tips": unit.tips}
        exposure = sum(unit.tree.branch_length(child) for _, child in unit.tree.edges())

        def peaked_likelihood(space, tree, tips, model, **kwargs):
            x = next(iter(model.rates.values())) * exposure
            return {"log_likelihood": -(x - 250.) ** 2}

        with mock.patch("intraphy.inference.genomic_exon_rates.evaluate_model",
                        side_effect=peaked_likelihood):
            fit = fit_family_rate([fit_unit])
        self.assertEqual(fit["status"], "estimated_conditional_composite_rate")
        self.assertTrue(fit["converged"])
        self.assertAlmostEqual(fit["dimensionless_rate_x"], 250., places=3)
        self.assertGreater(fit["dimensionless_rate_x"], 200.)
        self.assertGreater(fit["upper_tail_diagnostic"]["x_probe"],
                           fit["dimensionless_rate_x"])
        self.assertLess(fit["upper_tail_diagnostic"]["log_likelihood"],
                        fit["log_likelihood"])
        self.assertEqual(fit["upper_tail_diagnostic"]["reason"], "lower_probe")
        self.assertIn("finite_probe_only", fit["upper_tail_diagnostic"]["scope"])

    def test_zero_rate_nonidentified_replicates_remain_and_fixed_arm_scores(self):
        with tempfile.TemporaryDirectory() as temporary, \
                mock.patch.object(calibration, "fit_family_rate", return_value={
                    "converged": False, "mu": None, "status": "no_rate_information"}):
            result = calibration.calibrate_genomic_exons(temporary, "geometry", 0., 1, 2, 7)
            replicates = json.loads((Path(temporary) / "replicates.json").read_text())
            self.assertEqual(len(replicates), 2)
            self.assertEqual([row["estimated_status"] for row in replicates], ["nonidentified"] * 2)
            self.assertEqual([row["fixed_status"] for row in replicates], ["scored"] * 2)
            self.assertEqual(result["estimated"]["requested_replicates"], 2)
            self.assertEqual(result["estimated"]["nonidentified_replicates"], 2)
            self.assertEqual(result["estimated"]["scored_replicates"], 0)
            self.assertEqual(result["fixed"]["scored_replicates"], 2)
            metadata = json.loads((Path(temporary) / "calibration_metadata.json").read_text())
            self.assertEqual(metadata["posterior_backend"], "compact_origin_subset")
            self.assertEqual(metadata["scenario"], "geometry")

    def test_numerical_fit_failure_is_recorded_without_suppressing_fixed_scores(self):
        with tempfile.TemporaryDirectory() as temporary, \
                mock.patch.object(calibration, "fit_family_rate", side_effect=FloatingPointError("fit overflow")):
            result = calibration.calibrate_genomic_exons(temporary, "geometry", .2, 1, 1, 8)
            row = json.loads((Path(temporary) / "replicates.json").read_text())[0]
            self.assertEqual(row["fit_error_type"], "FloatingPointError")
            self.assertEqual(row["fit_error_message"], "fit overflow")
            self.assertEqual(row["fixed_status"], "scored")
            self.assertEqual(row["estimated_status"], "failed")
            self.assertEqual(result["fixed"]["scored_replicates"], 1)
            self.assertEqual(result["estimated"]["failed_replicates"], 1)

    def test_nonconverged_rate_with_numeric_value_is_not_scored(self):
        with tempfile.TemporaryDirectory() as temporary, \
                mock.patch.object(calibration, "fit_family_rate", return_value={
                    "converged": False, "mu": .2, "status": "optimizer_failed"}):
            result = calibration.calibrate_genomic_exons(temporary, "geometry", .2, 1, 1, 15)
            row = json.loads((Path(temporary) / "replicates.json").read_text())[0]
            self.assertEqual(row["estimated_mu"], None)
            self.assertEqual(row["estimated_status"], "failed")
            self.assertEqual(row["fixed_status"], "scored")
            self.assertEqual(result["estimated"]["scored_replicates"], 0)
            self.assertEqual(result["fixed"]["scored_replicates"], 1)

    def test_incomplete_and_numeric_scoring_errors_are_recorded(self):
        unit = self.template()
        model = RateModel({kind: .2 for kind in EDIT_KINDS})
        observed, assigned, origins = sample_unit_history(unit, model, np.random.default_rng(31))
        with mock.patch("intraphy.inference.genomic_exon_branches.branch_change_rows", return_value=[]):
            with self.assertRaisesRegex(ValueError, "Incomplete branch scores"):
                score_replicate("fixed", 1, [unit], [(observed, assigned, origins)], model)
        with tempfile.TemporaryDirectory() as temporary, \
                mock.patch.object(calibration, "fit_family_rate", return_value={
                    "converged": True, "mu": .2, "status": "estimated_conditional_composite_rate"}), \
                mock.patch("intraphy.inference.genomic_exon_branches.branch_change_rows",
                           side_effect=FloatingPointError("endpoint overflow")):
            calibration.calibrate_genomic_exons(temporary, "geometry", .2, 1, 1, 15)
            row = json.loads((Path(temporary) / "replicates.json").read_text())[0]
            self.assertEqual(row["fixed_error_type"], "FloatingPointError")
            self.assertEqual(row["fixed_error_message"], "endpoint overflow")
            self.assertEqual(row["estimated_error_type"], "FloatingPointError")
            self.assertEqual(row["fixed_status"], "failed")

    def test_summary_mcse_uses_replicate_means_and_places_one_in_final_bin(self):
        replicates = [{"replicate": 1, "fixed_status": "scored"},
                      {"replicate": 2, "fixed_status": "scored"},
                      {"replicate": 3, "fixed_status": "failed"}]
        branches = [
            {"replicate": 1, "probability_change": 0., "truth_changed": False, "brier": 0., "modal_hit": True, "modal_credit": 1.},
            {"replicate": 1, "probability_change": 1., "truth_changed": True, "brier": 0., "modal_hit": True, "modal_credit": .5},
            {"replicate": 2, "probability_change": 0., "truth_changed": True, "brier": 1., "modal_hit": False, "modal_credit": 0.},
            {"replicate": 2, "probability_change": 1., "truth_changed": False, "brier": 1., "modal_hit": False, "modal_credit": 0.},
        ]
        result = summarize_arm("fixed", replicates, branches, 3, 1, 2)
        self.assertEqual(result["brier"], .5)
        self.assertEqual(result["brier_replicate_mean_monte_carlo_se"], .5)
        self.assertEqual(result["failed_replicates"], 1)
        self.assertEqual(result["scored_replicates"], 2)
        self.assertEqual(result["requested_branch_rows"], 6)
        self.assertEqual(result["scored_branch_rows"], 4)
        self.assertEqual(result["reliability_bins"][-1]["count"], 2)

    def test_cli_runs_small_real_model_and_writes_output(self):
        from intraphy.cli import main
        with tempfile.TemporaryDirectory() as temporary:
            status = main(["calibrate-exons", "--scenario", "geometry", "--rate", ".2",
                "--units", "2", "--replicates", "1", "--seed", "12",
                "--output-dir", temporary])
            self.assertEqual(status, 0)
            summary = json.loads((Path(temporary) / "summary.json").read_text())
            self.assertEqual(summary["fixed"]["requested_replicates"], 1)
            self.assertTrue((Path(temporary) / "branches.tsv").is_file())

    def test_cli_parser_requires_positive_counts(self):
        args = build_parser().parse_args(["calibrate-exons", "--scenario", "geometry",
            "--rate", "0", "--units", "2", "--replicates", "3", "--seed", "9",
            "--output-dir", "out"])
        self.assertEqual((args.scenario, args.rate, args.units, args.replicates, args.seed),
                         ("geometry", 0., 2, 3, 9))
        for name, value in (("--units", "0"), ("--replicates", "0")):
            with self.subTest(name=name):
                argv = ["calibrate-exons", "--scenario", "geometry", "--rate", ".2",
                        "--units", "1", "--replicates", "1", "--output-dir", "out"]
                argv[argv.index(name) + 1] = value
                with self.assertRaises(SystemExit):
                    build_parser().parse_args(argv)


if __name__ == "__main__":
    unittest.main()
