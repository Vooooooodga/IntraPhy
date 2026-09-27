"""Source-only tests for locus CLI routing and persisted fit status."""
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from intraphy.commands.parser import build_parser
from intraphy.commands.preflight import validate_arguments, validate_input_paths
from intraphy.commands.locus import dispatch_locus
from intraphy.inference.locus_io import load_locus_model
from intraphy.inference.locus_run import analyze_locus


def _write_input(directory):
    model = {
        "schema": "intraphy.exon-locus-model/1", "model": "exon-locus-ctmc",
        "branch_length_unit": "declared opportunity units",
        "tree_provenance": "synthetic tree", "provenance": "synthetic model",
        "independence_provenance": "one connected locus",
        "rates": {"provenance": "test values", "groups": {
            "deletion": {"mode": "fit", "initial": 0.2}}},
        "units": [{
            "family": "family", "unit": "unit",
            "catalogue": {
                "provenance": "supported synthetic catalogue",
                "material": [{"id": "m", "start": 10, "end": 20}],
                "copies": [{"id": "copy", "material_ids": ["m"]}],
                "features": [{"id": "exon", "kind": "exon", "required_material": ["m"],
                              "copy_id": "copy", "start": 10, "end": 20}],
                "opportunities": [{"id": "deletion", "outcome_id": "delete",
                                   "kind": "dna_deletion", "rate_group": "deletion",
                                   "material_deletions": ["m"], "interval": [10, 20]}],
            },
            "root": {"provenance": "explicit root", "entries": [
                {"state": {"material": [1], "active_features": ["exon"]}, "weight": 1.0}]},
            "observations": {"provenance": "synthetic tips", "tips": {
                "A": {"material": [1]}, "B": {"material": [0]}}},
        }],
    }
    model_path, tree_path = directory / "model.json", directory / "tree.tsv"
    model_path.write_text(json.dumps(model), encoding="utf-8")
    tree_path.write_text("node_id\tparent_id\tlabel\tbranch_length\n"
                         "root\t\troot\t\nA\troot\tA\t0.2\nB\troot\tB\t0.3\n", encoding="utf-8")
    return model_path, tree_path


def _fit_result(*, success=True, rates=None):
    return SimpleNamespace(
        initial_log_likelihood=-2.0, fitted_log_likelihood=-1.5,
        rates={"deletion": 0.3} if rates is None else rates,
        optimizer_success=success, optimizer_status=0 if success else 2,
        optimizer_message="synthetic optimizer report", iterations=4,
        boundary_rates=(), curvature_eigenvalues=(1.0,), curvature_rank=1,
        curvature_dimension=1, curvature_positive_definite=True,
        curvature_has_negative_eigenvalue=False, locally_flat=False, workers=1)


class LocusCommandTests(unittest.TestCase):
    def test_all_three_commands_accept_locus_inputs_without_prepared_input_dir(self):
        parser = build_parser()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model_path, tree_path = _write_input(root)
            for command in ("analyze", "run", "infer-phylogeny"):
                with self.subTest(command=command):
                    args = parser.parse_args([command, "--locus-model", str(model_path),
                                              "--species-tree", str(tree_path), "--output-dir", str(root / "out"),
                                              "--threads", "3"])
                    self.assertIsNone(getattr(args, "input_dir", None))
                    self.assertEqual(args.model, "exon-locus-ctmc")
                    self.assertEqual(args.parameter_mode, "fit")
                    validate_arguments(args)
                    validate_input_paths(args)
                    self.assertTrue(hasattr(args, "_locus_bundle"))

    def test_required_locus_inputs_and_thread_count_are_preflighted(self):
        parser = build_parser()
        for command in ("analyze", "run", "infer-phylogeny"):
            args = parser.parse_args([command, "--output-dir", "out"])
            with self.subTest(command=command, missing="model"):
                with self.assertRaisesRegex(ValueError, "requires --locus-model"):
                    validate_arguments(args)
            args = parser.parse_args([command, "--locus-model", "model.json", "--output-dir", "out"])
            with self.subTest(command=command, missing="tree"):
                with self.assertRaisesRegex(ValueError, "requires --species-tree"):
                    validate_arguments(args)
        args = parser.parse_args(["analyze", "--locus-model", "model.json", "--species-tree", "tree.tsv",
                                  "--output-dir", "out", "--threads", "0"])
        with self.assertRaisesRegex(ValueError, "--threads must be at least 1"):
            validate_arguments(args)

    def test_expected_edits_is_an_explicit_opt_in(self):
        parser = build_parser()
        for command in ("analyze", "run", "infer-phylogeny"):
            base = [command, "--locus-model", "model.json", "--species-tree", "tree.tsv", "--output-dir", "out"]
            self.assertFalse(parser.parse_args(base).expected_edits)
            self.assertTrue(parser.parse_args(base + ["--expected-edits"]).expected_edits)

    def test_dispatch_passes_expected_count_flag_and_thread_count(self):
        for expected_edits in (False, True):
            args = SimpleNamespace(_locus_bundle=object(), output_dir="out",
                                   parameter_mode="fit", expected_edits=expected_edits, threads=7)
            with self.subTest(expected_edits=expected_edits), patch(
                    "intraphy.inference.locus_run.analyze_locus") as analyze:
                dispatch_locus(args)
            analyze.assert_called_once_with(args._locus_bundle, "out", parameter_mode="fit",
                                            expected_counts=expected_edits, workers=7)

    def test_fixed_and_fitted_parameter_status_are_saved(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model_path, tree_path = _write_input(root)
            bundle = load_locus_model(model_path, tree_path)
            for mode, result, expected in (
                ("fixed", _fit_result(rates={"deletion": 0.2}), "fixed_parameter_evaluation"),
                ("fit", _fit_result(), "converged_local_solution"),
            ):
                output = root / mode
                with self.subTest(mode=mode), patch("intraphy.inference.locus_run.fit_locus_rates", return_value=result):
                    analyze_locus(bundle, output, parameter_mode=mode)
                fit_record = json.loads((output / "locus_fit.json").read_text(encoding="utf-8"))
                history = json.loads((output / "locus_history.json").read_text(encoding="utf-8"))
                self.assertEqual(fit_record["parameter_status"], expected)
                self.assertEqual(history["parameter_status"], expected)
                self.assertEqual(history["parameter_mode"], mode)

    def test_failed_fit_keeps_diagnostics_without_history(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model_path, tree_path = _write_input(root)
            bundle = load_locus_model(model_path, tree_path)
            output = root / "failed"
            with patch("intraphy.inference.locus_run.fit_locus_rates", return_value=_fit_result(success=False)):
                with self.assertRaisesRegex(RuntimeError, "see locus_fit.json"):
                    analyze_locus(bundle, output, parameter_mode="fit")
            record = json.loads((output / "locus_fit.json").read_text(encoding="utf-8"))
            result = json.loads((output / "run_result.json").read_text(encoding="utf-8"))
            self.assertEqual(record["parameter_status"], "optimizer_unsuccessful")
            self.assertFalse(record["optimizer_success"])
            self.assertEqual(result["status"], "failed")
            self.assertFalse((output / "locus_history.json").exists())

    def test_existing_result_is_refused_without_overwrite(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model_path, tree_path = _write_input(root)
            bundle = load_locus_model(model_path, tree_path)
            output = root / "occupied"
            output.mkdir()
            manifest = output / "run_result.json"
            manifest.write_text('{"status":"prior"}\n', encoding="utf-8")
            with patch("intraphy.inference.locus_run.fit_locus_rates") as fit:
                with self.assertRaisesRegex(ValueError, "Refusing existing locus result"):
                    analyze_locus(bundle, output)
            fit.assert_not_called()
            self.assertEqual(manifest.read_text(encoding="utf-8"), '{"status":"prior"}\n')


if __name__ == "__main__":
    unittest.main()
