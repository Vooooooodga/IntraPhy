import json
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

import numpy as np

from intraphy.commands.parser import build_parser
from intraphy.cli import main
from intraphy.inference.genomic_exon_comparison_calibration import (
    _catalogue_spaces,
    _mask_unit,
    _scoreable_native_fit,
    _unit_for_gene,
    _validate_tip_support,
    balanced_tree,
)
from intraphy.inference.configuration_compact import evaluate_compact
from intraphy.inference.configuration_model import RateModel, evaluate_model
from intraphy.inference.exon_resampling import sample_unit_history
from intraphy.inference.genomic_exon_comparison import compare_foreground
from intraphy.structure.edits import EDIT_KINDS
from intraphy.structure.types import Catalogue, ExonConfiguration, ExonSpan, ObservationEvidence
from intraphy.topology import SpeciesTree
from intraphy.verification.exon_comparison_scoring import (
    native_endpoint_scores,
    paired_panel_summary,
    panel_summary,
)


class ExonComparisonCalibrationTests(unittest.TestCase):
    @staticmethod
    def finite_mle_family():
        labels = ("A", "B", "C", "D")
        span = ExonSpan(0, 2)
        observations = tuple(ObservationEvidence(label, (ExonConfiguration(()),))
                             for label in labels)
        catalogue = Catalogue("finite_gene", "single_span", 2, (span,), (),
            observations=observations, discovery="independent_catalogue",
            boundary_candidates=(span,), observation_unit="genomic_exon_spans")
        from intraphy.structure.space import enumerate_space
        space = enumerate_space(catalogue)
        tree_rows = (
            {"node_id": "root", "parent_id": "", "label": "root", "branch_length": 0.},
            {"node_id": "ab", "parent_id": "root", "label": "ab", "branch_length": .5},
            {"node_id": "cd", "parent_id": "root", "label": "cd", "branch_length": .5},
            {"node_id": "A", "parent_id": "ab", "label": "A", "branch_length": .5},
            {"node_id": "B", "parent_id": "ab", "label": "B", "branch_length": .5},
            {"node_id": "C", "parent_id": "cd", "label": "C", "branch_length": .5},
            {"node_id": "D", "parent_id": "cd", "label": "D", "branch_length": .5})
        tree = SpeciesTree(tree_rows)
        state_index = {bool(state.exons): index for index, state in enumerate(space.states)}
        patterns = ((0, 0, 0, 0),) * 6 + ((1, 1, 1, 1),) * 6 + (
                    (1, 0, 0, 0), (1, 0, 0, 0), (1, 0, 0, 0), (1, 0, 0, 0),
                    (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1),
                    (1, 1, 0, 0), (1, 1, 0, 0), (1, 1, 0, 0),
                    (1, 0, 1, 0), (0, 1, 0, 1), (1, 0, 0, 1),
                    (0, 1, 1, 0), (0, 0, 1, 1), (1, 1, 1, 0),
                    (1, 0, 1, 1), (0, 1, 1, 1))
        units = []
        for index, pattern in enumerate(patterns):
            tips = {}
            for label, present in zip(labels, pattern):
                vector = np.zeros(len(space.states), dtype=float)
                vector[state_index[bool(present)]] = 1.
                tips[label] = vector
            units.append({"unit_id": f"unit_{index + 1:02d}", "space": space,
                          "tree": tree, "tips": tips})
        return {"finite_gene": units}

    def test_cli_exposes_fixed_comparison_calibration_arguments(self):
        args = build_parser().parse_args([
            "calibrate-exon-comparison", "--scenario", "geometry",
            "--foreground-multiplier", "1", "--observation-mask", "missing-tip",
            "--gene-rates", ".1", ".2", ".4", ".8", "--taxa", "4",
            "--replicates", "2", "--seed", "17", "--threads", "2",
            "--output-dir", "calibration"])
        self.assertEqual(args.command, "calibrate-exon-comparison")
        self.assertEqual(args.gene_rates, [.1, .2, .4, .8])
        self.assertEqual(args.observation_mask, "missing-tip")

    def test_exact_visible_support_allows_material_zero_and_two_but_not_partial(self):
        catalogues, spaces = _catalogue_spaces("shared-deletion")
        tree_rows, labels = balanced_tree(4)
        tree = SpeciesTree(tree_rows)
        unit = _unit_for_gene(catalogues[0], spaces[0], tree, labels, "gene")
        by_visible = {}
        for index, state in enumerate(unit.space.states):
            visible = (tuple((span.start, span.end) for span in state.exons),
                       tuple(int(value == 1) for value in state.material))
            by_visible.setdefault(visible, {}).setdefault(state.material, []).append(index)
        material_pair = next(value for value in by_visible.values()
                             if (0,) in value and (2,) in value)
        indices = material_pair[(0,)] + material_pair[(2,)]
        tips = {label: np.zeros(len(unit.space.states)) for label in labels}
        for label in labels:
            tips[label][indices] = 1.
        exact = replace(unit, tips=tips)
        _validate_tip_support(exact)
        partial_tips = {label: values.copy() for label, values in tips.items()}
        signatures = {}
        for index, state in enumerate(unit.space.states):
            visible = (tuple((span.start, span.end) for span in state.exons),
                       tuple(int(value == 1) for value in state.material))
            signatures.setdefault(visible, []).append(index)
        visible_indices = list(signatures.values())
        partial_tips[labels[0]][:] = 0.
        partial_tips[labels[0]][visible_indices[0][0]] = 1.
        partial_tips[labels[0]][visible_indices[1][0]] = 1.
        with self.assertRaisesRegex(ValueError, "Partial candidate"):
            _validate_tip_support(replace(unit, tips=partial_tips))

    def test_mask_changes_only_one_whole_tip_to_unknown(self):
        catalogues, spaces = _catalogue_spaces("geometry")
        rows, labels = balanced_tree(4)
        unit = _unit_for_gene(catalogues[0], spaces[0], SpeciesTree(rows), labels, "g")
        masked = _mask_unit(unit, labels[-1])
        np.testing.assert_array_equal(masked.tips[labels[-1]], np.ones(len(unit.space.states)))
        for label in labels[:-1]:
            np.testing.assert_array_equal(masked.tips[label], unit.tips[label])
        _validate_tip_support(masked)

    def test_fit_filter_rejects_flat_nonidentified_and_nonfinite_records(self):
        base = {"status": "estimated_conditional_composite_rate",
                "mu": .3, "log_likelihood": -2.}
        self.assertTrue(_scoreable_native_fit(base))
        self.assertTrue(_scoreable_native_fit({**base, "status": "zero_boundary", "mu": 0.}))
        for bad in (
                {**base, "status": "flat_or_nonidentified", "mu": None},
                {**base, "status": "profile_incomplete"}, {**base, "mu": float("inf")},
                {**base, "log_likelihood": float("nan")}):
            self.assertFalse(_scoreable_native_fit(bad))

    def test_actual_family_fit_profile_supplies_scoreable_rate(self):
        result = compare_foreground(self.finite_mle_family(), {"A"},
                                    profile_multipliers=(.5, 1., 2.))
        alternative = result["alternative"]
        self.assertEqual(alternative["status"], "estimated_conditional_composite_multiplier",
                         msg=f"unexpected alternative fit record: {alternative!r}")
        self.assertTrue(np.isfinite(alternative["foreground_multiplier"]))
        self.assertIsNotNone(result["comparison_statistic"])
        self.assertGreaterEqual(result["comparison_statistic"], -1e-8)
        profile = next(point for point in result["profile"]
                       if point["foreground_multiplier"] == 1.)
        self.assertTrue(_scoreable_native_fit(profile["families"]["finite_gene"]))

    def test_dense_native_endpoint_score_matches_compact_change_probability(self):
        catalogue, spaces = _catalogue_spaces("geometry")
        rows, labels = balanced_tree(4)
        tree = SpeciesTree(rows)
        unit = _unit_for_gene(catalogue[0], spaces[0], tree, labels, "gene")
        foreground = frozenset((labels[0],))
        model = RateModel({kind: .2 for kind in EDIT_KINDS}, foreground=foreground)
        observed, assigned, _origins = sample_unit_history(
            unit, model, np.random.default_rng(83))
        dense_scores = native_endpoint_scores(observed, assigned, .2, foreground, 1.)
        dense = evaluate_model(observed.space, observed.tree, observed.tips, model,
                               posterior=True, counts=False, backend="dense")
        compact = evaluate_compact(observed.space, observed.tree, observed.tips, model,
                                   counts=False)
        compact_rows = {(row["parent"], row["child"]): row
                        for row in compact["branches"]}
        self.assertEqual(set(compact_rows), set(observed.tree.edges()))
        self.assertAlmostEqual(dense["log_likelihood"], compact["log_likelihood"], places=10)
        for score in dense_scores:
            key = (score["parent_node"], score["child_node"])
            probability = compact_rows[key]["observable_endpoint_summary"][
                "probability_exon_structure_change"]
            self.assertAlmostEqual(score["probability_exon_geometry_change"],
                                   probability, places=10)

    def test_panel_denominator_excludes_incomplete_and_keeps_zero_log_status(self):
        complete = {"status": "scored", "branch_count": 2, "mean_endpoint_brier": .2,
                    "mean_true_visible_pair_log_score": None,
                    "zero_probability_true_pairs": 1, "log_score_status": "negative_infinity"}
        genes = [{"joint_history": complete.copy()} for _ in range(4)]
        panel = panel_summary(genes, "joint_history", 2)
        self.assertEqual(panel["status"], "scored")
        self.assertEqual(panel["mean_endpoint_brier"], .2)
        self.assertIsNone(panel["mean_true_visible_pair_log_score"])
        self.assertEqual(panel["log_score_status"], "negative_infinity")
        genes[-1]["joint_history"] = {"status": "failed", "branch_count": 0}
        incomplete = panel_summary(genes, "joint_history", 2)
        self.assertEqual(incomplete["status"], "incomplete")
        self.assertIsNone(incomplete["mean_endpoint_brier"])
        self.assertEqual(incomplete["failed_genes"], 1)
        paired = paired_panel_summary(panel, incomplete)
        self.assertEqual(paired["status"], "incomplete")
        self.assertNotIn("independent_minus_native_brier", paired)

    def test_small_native_calibration_writes_fixed_roster_and_truth_free_fit_input(self):
        with tempfile.TemporaryDirectory() as directory, \
             tempfile.TemporaryDirectory() as parallel_directory:
            def argv_for(path, threads):
                return [
                    "calibrate-exon-comparison", "--scenario", "geometry",
                    "--foreground-multiplier", "1", "--observation-mask", "complete",
                    "--gene-rates", "0", "0", "0", "0", "--taxa", "4",
                    "--replicates", "2", "--seed", "29", "--threads", str(threads),
                    "--output-dir", str(path)]
            self.assertEqual(main(argv_for(directory, 1)), 0)
            self.assertEqual(main(argv_for(parallel_directory, 2)), 0)
            summary = json.loads((Path(directory) / "summary.json").read_text())
            parallel_summary = json.loads((Path(parallel_directory) / "summary.json").read_text())
            output = Path(directory)
            self.assertTrue((output / "calibration_metadata.json").is_file())
            self.assertTrue((output / "replicates.jsonl").is_file())
            self.assertTrue((output / "summary.json").is_file())
            metadata = json.loads((output / "calibration_metadata.json").read_text())
            serial_records = [json.loads(line) for line in
                              (output / "replicates.jsonl").read_text().splitlines()]
            record = serial_records[0]
            self.assertEqual(metadata["worker_count"], 1)
            self.assertEqual(len(record["genes"]), 4)
            self.assertEqual(record["design"]["missing_tip"], None)
            self.assertEqual(record["foreground_comparison"]["p_value"], None)
            fitted_profile = {point["foreground_multiplier"]
                              for point in record["foreground_comparison"].get("profile", [])}
            self.assertTrue({.5, 1., 2.}.issubset(fitted_profile))
            self.assertNotIn("assigned", record)
            self.assertNotIn("origins", record)
            self.assertEqual(summary["requested_replicates"], 2)
            self.assertTrue(all("joint_fit_diagnostic" in gene for gene in record["genes"]))
            self.assertTrue(all("independent_fit" in gene for gene in record["genes"]))
            parallel_metadata = json.loads(
                (Path(parallel_directory) / "calibration_metadata.json").read_text())
            self.assertEqual(parallel_metadata["worker_count"], 2)
            parallel_records = [json.loads(line) for line in
                                (Path(parallel_directory) / "replicates.jsonl").read_text().splitlines()]
            self.assertEqual({row["replicate"] for row in serial_records}, {1, 2})
            self.assertEqual({row["replicate"] for row in parallel_records}, {1, 2})
            def without_timings(records):
                by_id = {}
                for row in records:
                    row.pop("elapsed_seconds", None)
                    row["foreground_comparison"].pop("elapsed_seconds", None)
                    by_id[row["replicate"]] = row
                return by_id
            self.assertEqual(without_timings(serial_records), without_timings(parallel_records))
            self.assertEqual(summary, parallel_summary)


if __name__ == "__main__":
    unittest.main()
