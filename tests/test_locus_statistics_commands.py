"""Source-only tests for the conditional statistics and evidence CLI surfaces."""
import json
import tempfile
import unittest
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from intraphy.cli import _dispatch
from intraphy.commands.locus_statistics import dispatch_prepare_locus_evidence
from intraphy.commands.parser import build_parser
from intraphy.commands.preflight import validate_arguments, validate_input_paths
from intraphy.inference.locus_statistics_run import run_locus_statistics


@dataclass
class _Fit:
    optimizer_success: bool = True
    unresolved_higher_likelihood: bool = False


@dataclass
class _Comparison:
    null_fit: object
    full_fit: object
    statistic: float | None = 1.2
    status: str = "evaluated_within_numerical_tolerance"
    raw_log_likelihood_difference: float | None = 0.6
    likelihood_tolerance: float = 1e-10
    interpretation_limit: str = "conditional comparison"


def _model():
    return {
        "schema": "intraphy.exon-locus-model/1", "model": "exon-locus-ctmc",
        "branch_length_unit": "declared opportunity units",
        "tree_provenance": "test tree", "provenance": "test model",
        "independence_provenance": "one locus unit",
        "rates": {"provenance": "fixed test value", "groups": {
            "deletion": {"mode": "fit", "initial": 0.2}}},
        "units": [{
            "family": "family", "unit": "relay",
            "catalogue": {
                "provenance": "candidate catalogue",
                "material": [{"id": "copy_A", "start": 10, "end": 20}],
                "copies": [{"id": "A", "material_ids": ["copy_A"]}],
                "features": [],
                "opportunities": [{"id": "delete_A", "outcome_id": "delete",
                                   "kind": "dna_deletion", "rate_group": "deletion",
                                   "material_deletions": ["copy_A"], "interval": [10, 20]}],
            },
            "root": {"provenance": "declared root", "entries": [
                {"state": {"material": [1], "active_features": []}, "weight": 1.0}]},
            "observations": {"provenance": "template calls", "tips": {
                "A": {"material": [None], "features": {}},
                "B": {"material": [None], "features": {}},
            }},
        }],
    }


def _write_inputs(directory):
    model_path = directory / "template.json"
    tree_path = directory / "tree.tsv"
    rows_path = directory / "evidence.json"
    model_path.write_text(json.dumps(_model()), encoding="utf-8")
    tree_path.write_text("node_id\tparent_id\tlabel\tbranch_length\n"
                         "root\t\troot\t\nA\troot\tA\t0.5\nB\troot\tB\t0.5\n",
                         encoding="utf-8")
    rows_path.write_text(json.dumps([
        {"family": "family", "unit": "relay", "tip": "A", "material_id": "copy_A",
         "call": "present", "evidence": "assembly interval", "survey": True,
         "sensitivity": 0.98, "specificity": 0.99},
        {"family": "family", "unit": "relay", "tip": "B", "material_id": "copy_A",
         "call": "unknown", "evidence": "no callable locus", "survey": False},
    ]), encoding="utf-8")
    return model_path, tree_path, rows_path


class LocusStatisticsCommandTests(unittest.TestCase):
    def test_parser_registers_statistics_and_evidence_commands(self):
        parser = build_parser()
        stats = parser.parse_args([
            "locus-statistics", "--locus-model", "m.json", "--species-tree", "t.tsv",
            "--output-dir", "out", "--null-rate-group", "dup", "--null-rate-group", "del",
        ])
        self.assertEqual(stats.model, "exon-locus-ctmc")
        self.assertEqual(stats.start_scales, [0.2, 1.0, 5.0])
        self.assertEqual(stats.null_rate_group, ["dup", "del"])
        self.assertEqual(stats.bootstrap_replicates, 0)
        prepare = parser.parse_args([
            "prepare-locus-evidence", "--locus-model", "m.json", "--species-tree", "t.tsv",
            "--evidence-json", "rows.json", "--output-dir", "out",
        ])
        self.assertEqual(prepare.command, "prepare-locus-evidence")

    def test_statistics_argument_pairs_and_bootstrap_prerequisites(self):
        parser = build_parser()
        args = parser.parse_args([
            "locus-statistics", "--locus-model", "m.json", "--species-tree", "t.tsv",
            "--output-dir", "out", "--null-rate-group", "dup", "--bootstrap-replicates", "4",
        ])
        with self.assertRaisesRegex(ValueError, "--seed is required"):
            validate_arguments(args)
        args.seed = 17
        with self.assertRaisesRegex(ValueError, "sampling-design"):
            validate_arguments(args)
        args.sampling_design = "fixed_catalogue"
        args.profile_rate_group = "dup"
        with self.assertRaisesRegex(ValueError, "supplied together"):
            validate_arguments(args)

    def test_statistics_preflight_loads_model_then_checks_free_group(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model_path, tree_path, _ = _write_inputs(root)
            args = build_parser().parse_args([
                "locus-statistics", "--locus-model", str(model_path), "--species-tree", str(tree_path),
                "--output-dir", str(root / "out"), "--null-rate-group", "deletion",
            ])
            validate_arguments(args)
            validate_input_paths(args)
            self.assertEqual(args._locus_bundle.rates["deletion"]["mode"], "fit")
            args.null_rate_group = ["fixed_group"]
            with self.assertRaisesRegex(ValueError, "declared free"):
                validate_input_paths(args)

    def test_prepare_evidence_preflight_and_output_routing(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model_path, tree_path, rows_path = _write_inputs(root)
            output = root / "prepared"
            args = build_parser().parse_args([
                "prepare-locus-evidence", "--locus-model", str(model_path),
                "--species-tree", str(tree_path), "--evidence-json", str(rows_path),
                "--output-dir", str(output),
            ])
            validate_arguments(args)
            validate_input_paths(args)
            self.assertEqual(args._prepared_model_record["units"][0]["observations"]["tips"]["A"]["material"], [1])
            self.assertEqual(args._prepared_model_record["units"][0]["observations"]["tips"]["B"]["material"], [None])
            output.mkdir()
            dispatch_prepare_locus_evidence(args)
            saved = json.loads((output / "locus_model.json").read_text(encoding="utf-8"))
            result = json.loads((output / "run_result.json").read_text(encoding="utf-8"))
            self.assertEqual(saved["units"][0]["observations"]["tips"]["A"]["material_sensitivity"],
                             {"copy_A": 0.98})
            self.assertEqual(saved["units"][0]["observations"]["tips"]["A"]["features"],
                             {})
            self.assertEqual(result["artifacts"], ["locus_model.json", "species_tree.tsv"])

    def test_cli_routes_both_commands_to_the_new_dispatchers(self):
        stats = SimpleNamespace(command="locus-statistics")
        with patch("intraphy.commands.locus_statistics.dispatch_locus_statistics") as dispatch:
            _dispatch(stats)
        dispatch.assert_called_once_with(stats)
        prepare = SimpleNamespace(command="prepare-locus-evidence")
        with patch("intraphy.commands.locus_statistics.dispatch_prepare_locus_evidence") as dispatch:
            _dispatch(prepare)
        dispatch.assert_called_once_with(prepare)

    def test_statistics_writer_records_json_tree_and_result_ownership(self):
        bundle = SimpleNamespace(
            rates={"dup": {"mode": "fit", "value": 0.2}},
            model_record={"schema": "test"}, tree_rows=({
                "node_id": "root", "parent_id": "", "label": "root", "branch_length": "",
            },), tree_provenance="test tree", branch_length_unit="arbitrary units",
            provenance="test model", units=(), fit_units=lambda: (),
        )
        comparison = _Comparison(_Fit(), _Fit())
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "statistics"
            with patch("intraphy.inference.locus_statistics_run.compare_locus_rates",
                       return_value=comparison):
                run_locus_statistics(bundle, output, null_rate_groups=("dup",))
            record = json.loads((output / "locus_statistics.json").read_text(encoding="utf-8"))
            result = json.loads((output / "run_result.json").read_text(encoding="utf-8"))
            self.assertEqual(record["status"], "completed")
            self.assertEqual(record["comparison"]["statistic"], 1.2)
            self.assertTrue((output / "species_tree.tsv").is_file())
            self.assertEqual(result["artifacts"], ["locus_statistics.json", "species_tree.tsv"])


if __name__ == "__main__":
    unittest.main()
