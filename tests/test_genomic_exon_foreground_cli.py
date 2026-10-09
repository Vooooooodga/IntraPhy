"""CLI and result-schema contracts for native foreground comparisons."""
import json
from dataclasses import replace
from pathlib import Path
import random
import shutil
import tempfile
import unittest
from unittest.mock import patch

from intraphy.cli import main
from intraphy.inference.genomic_exon_calibration import TREE_ROWS, _catalogues
from intraphy.run_result import result_model
from intraphy.structure.prepare import read_preparation_summary, validate_preparation_catalogues
from intraphy.storage.tabular import read_tsv, write_tsv
from intraphy.structure.serialization import read_catalogues, write_catalogues
from intraphy.structure.types import Catalogue
from intraphy.verification.native_cases import build_native_example


ALIGNMENT_TOOLS = bool(shutil.which("mafft") and shutil.which("minimap2"))


class GenomicExonForegroundCliTests(unittest.TestCase):
    def test_empty_catalogue_requires_explicit_allow_empty(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "empty.jsonl"
            path.write_text("", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Configuration collection is empty"):
                read_catalogues(path)
            self.assertEqual(read_catalogues(path, allow_empty=True), ())

    def test_allow_empty_does_not_allow_invalid_catalogue_records(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "invalid.jsonl"
            path.write_text("{not-json}\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Invalid configuration"):
                read_catalogues(path, allow_empty=True)

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
            self.assertEqual(result["input_scope"]["report_scope"],
                             "provided_catalogues_only")
            self.assertIsNone(result["input_scope"]["original_target_families"])
            self.assertEqual(result["input_scope"]["requested_families"],
                             ["calibration_family", "excluded_family"])
            self.assertEqual(result["input_scope"]["supplied_catalogue_families"],
                             ["calibration_family", "excluded_family"])
            self.assertEqual(result["input_scope"]["family_scope_counts"], {
                "original_target_families": None,
                "generated_catalogue_families": None,
                "supplied_catalogue_families": 2,
                "eligible_fit_families": 1,
                "preparation_excluded_families": None,
                "catalogue_not_supplied_families": None,
            })
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

    def test_explicit_summary_reports_all_preparation_failures_without_catalogues(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            prepared = root / "prepared"
            prepared.mkdir()
            write_tsv(prepared / "species_tree.tsv", TREE_ROWS,
                      ["node_id", "parent_id", "label", "branch_length"])
            catalogue_path = root / "empty.jsonl"
            catalogue_path.write_text("", encoding="utf-8")
            preparation_summary = root / "exon_preparation_summary.tsv"
            write_tsv(preparation_summary, [
                {"family_id": "long_family", "status": "unresolved", "units": 0,
                 "qualified_units": 0,
                 "reason": "locus_alignment_budget_exceeded: test",
                 "evidence_directory": "family_00001"},
            ], ["family_id", "status", "units", "qualified_units", "reason", "evidence_directory"])
            foreground_path = root / "foreground.tsv"
            foreground_path.write_text("branch_scope\nab->A\n", encoding="utf-8")
            output = root / "result"

            status = main([
                "compare-exon-foreground", "--input-dir", str(prepared),
                "--exon-configurations", str(catalogue_path),
                "--preparation-summary", str(preparation_summary),
                "--foreground-branches", str(foreground_path),
                "--output-dir", str(output),
            ])

            self.assertEqual(status, 0)
            result = json.loads((output / "exon_foreground_comparison.json").read_text(
                encoding="utf-8"))
            self.assertEqual(result["scope"], "no_estimable_families")
            self.assertEqual(result["input_scope"]["family_scope_counts"], {
                "original_target_families": 1,
                "generated_catalogue_families": 0,
                "supplied_catalogue_families": 0,
                "eligible_fit_families": 0,
                "preparation_excluded_families": 1,
                "catalogue_not_supplied_families": 0,
            })
            self.assertEqual(result["input_scope"]["preparation_excluded_families"],
                             ["long_family"])
            self.assertEqual(result["null"]["families"]["long_family"]["status"],
                             "no_estimable_units")

    def test_prepared_family_can_be_intentionally_omitted_from_supplied_catalogues(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            prepared = root / "prepared"
            prepared.mkdir()
            write_tsv(prepared / "species_tree.tsv", TREE_ROWS,
                      ["node_id", "parent_id", "label", "branch_length"])
            catalogue_path = root / "empty.jsonl"
            catalogue_path.write_text("", encoding="utf-8")
            preparation_summary = root / "exon_preparation_summary.tsv"
            write_tsv(preparation_summary, [
                {"family_id": "prepared_but_omitted", "status": "prepared", "units": 3,
                 "qualified_units": 2, "reason": "", "evidence_directory": "family_00001"},
            ], ["family_id", "status", "units", "qualified_units", "reason", "evidence_directory"])
            foreground_path = root / "foreground.tsv"
            foreground_path.write_text("branch_scope\nab->A\n", encoding="utf-8")
            output = root / "result"

            self.assertEqual(main([
                "compare-exon-foreground", "--input-dir", str(prepared),
                "--exon-configurations", str(catalogue_path),
                "--preparation-summary", str(preparation_summary),
                "--foreground-branches", str(foreground_path),
                "--output-dir", str(output),
            ]), 0)
            result = json.loads((output / "exon_foreground_comparison.json").read_text(
                encoding="utf-8"))
            self.assertEqual(result["input_scope"]["requested_families"],
                             ["prepared_but_omitted"])
            self.assertEqual(result["input_scope"]["supplied_catalogue_families"], [])
            self.assertEqual(result["input_scope"]["family_scope_roster"][0]["fit_status"],
                             "catalogue_not_supplied")
            self.assertEqual(result["input_scope"]["family_scope_counts"][
                "catalogue_not_supplied_families"], 1)

    def test_preparation_summary_rejects_duplicate_targets_and_impossible_catalogue_units(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "summary.tsv"
            fields = ["family_id", "status", "units", "qualified_units", "reason",
                      "evidence_directory"]
            catalogue = _catalogues("shared-deletion")[0]
            row = {"family_id": catalogue.family, "status": "prepared", "units": 2,
                   "qualified_units": 1, "reason": "", "evidence_directory": "family_00001"}
            write_tsv(path, [row, row], fields)
            with self.assertRaisesRegex(ValueError, "duplicate family_id"):
                read_preparation_summary(path)
            write_tsv(path, [row], fields)
            summary = read_preparation_summary(path)
            self.assertEqual(validate_preparation_catalogues(summary, (catalogue,)),
                             (catalogue.family,))
            with self.assertRaisesRegex(ValueError, "absent from --preparation-summary"):
                validate_preparation_catalogues((dict(row, family_id="another_family"),),
                                                 (catalogue,))
            oversized = replace(catalogue, unit=catalogue.unit + "_extra")
            with self.assertRaisesRegex(ValueError, "exceeds preparation summary"):
                validate_preparation_catalogues((dict(row, units=1, qualified_units=1),),
                                                 (catalogue, oversized))
            write_tsv(path, [dict(row, status="skipped")], fields)
            with self.assertRaisesRegex(ValueError, "Invalid preparation status"):
                read_preparation_summary(path)
            write_tsv(path, [dict(row, units=1, qualified_units=2)], fields)
            with self.assertRaisesRegex(ValueError, "Invalid unit counts"):
                read_preparation_summary(path)

    def test_explicit_scope_counts_original_generated_supplied_and_fit_families(self):
        from intraphy.inference.genomic_exon_comparison_run import compare_genomic_exon_foreground

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            prepared = root / "prepared"
            prepared.mkdir()
            write_tsv(prepared / "species_tree.tsv", TREE_ROWS,
                      ["node_id", "parent_id", "label", "branch_length"])
            catalogue = _catalogues("shared-deletion")[0]
            catalogue_path = root / "catalogue.jsonl"
            write_catalogues(catalogue_path, (catalogue,))
            summary_path = root / "exon_preparation_summary.tsv"
            write_tsv(summary_path, [
                {"family_id": catalogue.family, "status": "prepared", "units": 1,
                 "qualified_units": 1, "reason": "", "evidence_directory": "family_00001"},
                {"family_id": "long_family", "status": "unresolved", "units": 0,
                 "qualified_units": 0, "reason": "locus_alignment_budget_exceeded: test",
                 "evidence_directory": "family_00002"},
            ], ["family_id", "status", "units", "qualified_units", "reason", "evidence_directory"])
            foreground = root / "foreground.tsv"
            foreground.write_text("branch_scope\nab->A\n", encoding="utf-8")
            output = root / "result"
            unit_diagnostic = {"family_id": catalogue.family, "unit_id": catalogue.unit,
                "catalogue_status": catalogue.status, "status": "eligible_conditional_unit",
                "probability_status": "not_computed", "state_count": 2,
                "state_space_complete": True, "informative_tips": 4, "unknown_tips": [],
                "state_space_reason": "", "reason": ""}
            unit_input = {"unit_id": catalogue.unit, "tree": object(), "tips": {}}
            family_fit = {catalogue.family: {"mu": 0.1, "log_likelihood": -1.0,
                "status": "conditional_family_fit", "diagnostic": {}}}
            comparison = {
                "null": {"status": "conditional_composite_log_likelihood",
                         "log_likelihood": -1.0, "foreground_multiplier": 1.0,
                         "families": family_fit},
                "alternative": {"status": "estimated_conditional_composite_multiplier",
                                "log_likelihood": -0.9, "foreground_multiplier": 1.2,
                                "best_evaluated_multiplier": 1.2, "families": family_fit},
                "profile": [],
            }
            with patch("intraphy.inference.genomic_exon_comparison_run._diagnose_catalogue",
                       return_value=(unit_diagnostic, None, unit_input)), \
                 patch("intraphy.inference.genomic_exon_comparison_run.compare_foreground",
                       return_value=comparison):
                result = compare_genomic_exon_foreground(
                    prepared, [catalogue_path], foreground, output,
                    preparation_summary=summary_path)

            scope = result["input_scope"]
            self.assertEqual(scope["report_scope"], "explicit_preparation_summary")
            self.assertEqual(scope["requested_families"], [catalogue.family, "long_family"])
            self.assertEqual(scope["supplied_catalogue_families"], [catalogue.family])
            self.assertEqual(scope["family_scope_counts"], {
                "original_target_families": 2,
                "generated_catalogue_families": 1,
                "supplied_catalogue_families": 1,
                "eligible_fit_families": 1,
                "preparation_excluded_families": 1,
                "catalogue_not_supplied_families": 0,
            })

    @unittest.skipUnless(ALIGNMENT_TOOLS, "Native preparation requires MAFFT and minimap2")
    def test_native_preparation_summary_flows_into_comparison_cli(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            raw = root / "raw"
            manifest = build_native_example(raw, "conserved")
            manifest_rows = read_tsv(manifest)
            rng = random.Random(32)
            codons = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"
                      if a + b + c not in {"TAA", "TAG", "TGA"}]
            short_exons = ("ATG" + "".join(rng.choice(codons) for _ in range(59)),
                           "".join(rng.choice(codons) for _ in range(60)),
                           "".join(rng.choice(codons) for _ in range(59)) + "TAA")
            intron = "GT" + "A" * 10 + "AG"
            short_gene = "".join(exon + (intron if index < 2 else "")
                                 for index, exon in enumerate(short_exons))
            for row in manifest_rows:
                species = row["species"]
                fasta = raw / row["genome_fasta"]
                header, sequence = fasta.read_text(encoding="utf-8").splitlines()
                fasta.write_text(f"{header}\n{sequence}{short_gene}\n", encoding="utf-8")
                annotation = raw / row["annotation_file"]
                lines = [f"{species}_chr\tsynthetic\tgene\t{len(sequence) + 1}\t"
                         f"{len(sequence) + len(short_gene)}\t.\t+\t.\tID={species}_short_gene",
                         f"{species}_chr\tsynthetic\tmRNA\t{len(sequence) + 1}\t"
                         f"{len(sequence) + len(short_gene)}\t.\t+\t.\t"
                         f"ID={species}_short_tx;Parent={species}_short_gene"]
                cursor = len(sequence) + 1
                for index, exon in enumerate(short_exons, 1):
                    end = cursor + len(exon) - 1
                    lines.extend((
                        f"{species}_chr\tsynthetic\texon\t{cursor}\t{end}\t.\t+\t.\t"
                        f"ID={species}_short_exon{index};Parent={species}_short_tx",
                        f"{species}_chr\tsynthetic\tCDS\t{cursor}\t{end}\t.\t+\t0\t"
                        f"ID={species}_short_cds{index};Parent={species}_short_tx"))
                    cursor = end + (len(intron) + 1 if index < 3 else 0)
                with annotation.open("a", encoding="utf-8") as handle:
                    handle.write("\n" + "\n".join(lines) + "\n")
            columns = list(manifest_rows[0])
            manifest_rows.extend({**row, "family_id": "short_family",
                                  "gene_id": f"{row['species']}_short_gene"}
                                 for row in manifest_rows.copy())
            write_tsv(manifest, manifest_rows, columns)

            prepared_result = root / "prepared_result"
            self.assertEqual(main([
                "analyze", "--manifest", str(manifest), "--species-tree",
                str(raw / "species_tree.nwk"), "--output-dir", str(prepared_result),
                "--flank", "100", "--max-locus-bases", "800", "--state-space-only",
            ]), 0)
            prepared = prepared_result / "prepared_inputs"
            branches = root / "foreground.tsv"
            tree_rows = read_tsv(prepared / "species_tree.tsv")
            branch = next(row for row in tree_rows if row["parent_id"])
            write_tsv(branches, [{"parent_node": branch["parent_id"],
                                  "child_node": branch["node_id"]}],
                      ["parent_node", "child_node"])
            comparison_result = root / "comparison_result"
            self.assertEqual(main([
                "compare-exon-foreground", "--input-dir", str(prepared),
                "--exon-configurations", str(prepared_result / "exon_configurations.jsonl"),
                "--preparation-summary",
                str(prepared_result / "exon_preparation_summary.tsv"),
                "--foreground-branches", str(branches), "--output-dir",
                str(comparison_result), "--max-states", "1",
            ]), 0)
            result = json.loads((comparison_result / "exon_foreground_comparison.json").read_text(
                encoding="utf-8"))
            scope = result["input_scope"]
            self.assertEqual(scope["requested_families"], ["example_gene", "short_family"])
            self.assertEqual(scope["generated_catalogue_families"], ["short_family"])
            self.assertEqual(scope["preparation_excluded_families"], ["example_gene"])
            self.assertEqual(scope["family_scope_counts"]["supplied_catalogue_families"], 1)

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
        self.assertIsNone(args.preparation_summary)
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
