"""CLI and result-schema contracts for native foreground comparisons."""
import json
from dataclasses import replace
from pathlib import Path
import tempfile
import unittest

from intraphy.cli import main
from intraphy.inference.genomic_exon_calibration import TREE_ROWS, _catalogues
from intraphy.run_result import result_model
from intraphy.storage.tabular import write_tsv
from intraphy.structure.serialization import write_catalogues
from intraphy.structure.types import Catalogue


class GenomicExonForegroundCliTests(unittest.TestCase):
    def test_multiple_catalogues_cli_preserves_requested_excluded_family_and_schema(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            prepared = root / "prepared"
            prepared.mkdir()
            write_tsv(prepared / "species_tree.tsv", TREE_ROWS,
                      ["node_id", "parent_id", "label", "branch_length"])

            eligible_path = root / "eligible.jsonl"
            eligible = _catalogues("shared-deletion")[0]
            observations = tuple(replace(observation, kind="unknown", configurations=(),
                                         material_presence=tuple(None for _ in eligible.material))
                                 if observation.species == "D" else observation
                                 for observation in eligible.observations)
            eligible = replace(eligible, observations=observations)
            write_catalogues(eligible_path, (eligible,))
            excluded = Catalogue("excluded_family", "unannotated", 10, (), (),
                status="unresolved", reasons=("no_annotated_exons",),
                observation_unit="genomic_exon_spans")
            excluded_path = root / "excluded.jsonl"
            write_catalogues(excluded_path, (excluded,))
            foreground_path = root / "foreground.tsv"
            foreground_path.write_text(
                "branch_scope\nab->A\n", encoding="utf-8")
            output = root / "result"

            status = main([
                "compare-exon-foreground", "--input-dir", str(prepared),
                "--exon-configurations", str(eligible_path), str(excluded_path),
                "--foreground-branches", str(foreground_path),
                "--output-dir", str(output), "--threads", "1",
                "--profile-multipliers", "0", "1",
            ])
            self.assertEqual(status, 0)
            result_path = output / "exon_foreground_comparison.json"
            result = json.loads(result_path.read_text(encoding="utf-8"))
            self.assertEqual(result["schema"], "intraphy.exon-foreground-comparison/1")
            self.assertEqual(result["model"], "exon-structure-ctmc")
            self.assertEqual(result["scope"], "single_gene")
            self.assertIsNone(result["p_value"])
            self.assertEqual(result["calibration_status"],
                             "not_calibrated_composite_likelihood")
            self.assertEqual(result["counts"]["requested_families"], 2)
            self.assertEqual(result["counts"]["eligible_families"], 1)
            self.assertEqual(result["counts"]["excluded_families"],
                             ["excluded_family"])
            self.assertEqual(result["foreground"]["canonical_children"], ["A"])
            self.assertTrue(result["input_scope"]["foreground_branch_provenance"])
            expected = {"calibration_family", "excluded_family"}
            for point in [result["null"], result["alternative"], *result["profile"]]:
                self.assertEqual(set(point["families"]), expected)
                self.assertEqual(point["families"]["excluded_family"]["status"],
                                 "no_estimable_units")
            self.assertEqual(result_model(output), "exon-structure-ctmc")
            run_result = json.loads((output / "run_result.json").read_text(encoding="utf-8"))
            self.assertEqual(run_result["status"], "completed_with_unresolved")
            self.assertEqual(run_result["analysis_scope"], "genomic-exon-foreground-comparison")
            self.assertTrue((output / "execution.json").is_file())
            self.assertTrue((output / "intraphy.log").is_file())
            diagnostic_rows = (output / "foreground_unit_diagnostics.tsv").read_text(
                encoding="utf-8").splitlines()
            self.assertIn('"D"', diagnostic_rows[1])
            for filename in (
                "foreground_profile.tsv", "foreground_family_fits.tsv",
                "foreground_unit_diagnostics.tsv", "foreground_branches.tsv",
                "state_space_diagnostics.jsonl", "species_tree.tsv",
                "exon_configurations.jsonl", "run_result.json",
            ):
                self.assertTrue((output / filename).is_file(), filename)

    def test_parser_accepts_nonnegative_profile_points_and_rejects_negative_values(self):
        from intraphy.commands.parser import build_parser
        from intraphy.commands.genomic_exon_comparison import validate_comparison_arguments

        parser = build_parser()
        args = parser.parse_args([
            "compare-exon-foreground", "--input-dir", "prepared",
            "--exon-configurations", "a.jsonl", "b.jsonl",
            "--foreground-branches", "foreground.tsv", "--output-dir", "result",
            "--profile-multipliers", "0", "1.5",
        ])
        self.assertEqual(args.exon_configurations, ["a.jsonl", "b.jsonl"])
        self.assertEqual(args.profile_multipliers, [0., 1.5])
        for value in ("-1", "nan", "inf"):
            args = parser.parse_args([
                "compare-exon-foreground", "--input-dir", "prepared",
                "--exon-configurations", "a.jsonl", "--foreground-branches", "foreground.tsv",
                "--output-dir", "result", "--profile-multipliers", value,
            ])
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, "nonnegative"):
                validate_comparison_arguments(args)

    def test_ambiguous_label_in_branch_scope_is_rejected(self):
        from intraphy.inputs.foreground import read_canonical_foreground
        from intraphy.structure.tree_context import canonical_tree
        from intraphy.topology import SpeciesTree

        tree = SpeciesTree([
            {"node_id": "root", "parent_id": "", "label": "root"},
            {"node_id": "inner", "parent_id": "root", "label": "shared",
             "branch_length": .2},
            {"node_id": "a", "parent_id": "inner", "label": "A", "branch_length": .3},
            {"node_id": "b", "parent_id": "inner", "label": "B", "branch_length": .4},
            {"node_id": "c", "parent_id": "root", "label": "shared",
             "branch_length": .5},
        ])
        with tempfile.TemporaryDirectory() as temporary:
            branch_file = Path(temporary) / "branches.tsv"
            branch_file.write_text("branch_scope\nroot->shared\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "ambiguous node reference"):
                read_canonical_foreground(branch_file, canonical_tree(tree), "supplied")


if __name__ == "__main__":
    unittest.main()
