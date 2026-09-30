import json
from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import Mock, patch

from intraphy.commands.parser import build_parser
from intraphy.commands.genomic_exon_preflight import validate_analysis_options, validate_input_paths
from intraphy.commands.genomic_exons import dispatch_analyze
from intraphy.commands.preflight import required_tools
from intraphy.run_result import RunResult, result_model


class GenomicExonCliTests(unittest.TestCase):
    def test_analyze_defaults_to_genomic_exon_spans_and_keeps_other_models_explicit(self):
        parser = build_parser()
        args = parser.parse_args(["analyze", "--fasta", "genomes", "--gff", "annotations",
                                  "--species-tree", "tree.nwk", "--output-dir", "result"])
        self.assertEqual(args.model, "exon-structure-ctmc")
        self.assertIsNone(args.root_frequency)
        self.assertEqual(parser.parse_args(["analyze", "--model", "dna-presence-ctmc",
                                            "--output-dir", "result"]).root_frequency,
                         "stationary")

    def test_state_space_only_is_limited_to_exon_structure_model(self):
        parser = build_parser()
        args = parser.parse_args(["analyze", "--input-dir", "prepared", "--state-space-only",
                                  "--output-dir", "result"])
        self.assertTrue(args.state_space_only)
        args.model = "dna-presence-ctmc"
        with self.assertRaisesRegex(ValueError, "requires --model exon-structure-ctmc"):
            validate_analysis_options(args)

    def test_state_space_only_rejects_state_limit(self):
        args = build_parser().parse_args(["analyze", "--input-dir", "prepared",
            "--state-space-only", "--max-states", "64", "--output-dir", "result"])
        with self.assertRaisesRegex(ValueError, "omit --max-states"):
            validate_analysis_options(args)

    def test_state_space_only_manifest_is_not_a_renderable_completed_run(self):
        with tempfile.TemporaryDirectory() as directory:
            RunResult("exon-structure-ctmc", "genomic-exon-state-space-diagnostic",
                      "species_tree.tsv", ("state_space_diagnostics.jsonl",),
                      status="state_space_only").write(directory)
            with self.assertRaisesRegex(ValueError, "incomplete run"):
                result_model(directory)

    def test_raw_dispatch_skips_legacy_correspondence_derivation(self):
        selection = Mock()
        selection.rows = [{"family_id": "f"}]
        args = SimpleNamespace(
            output_dir="result", input_dir=None, species_tree="tree.nwk",
            _input_selection=selection,
            exon_configurations=None, exon_rates=None, parameter_mode="fit",
            max_states=None, max_origin_scenarios=None, branch_length_mode="supplied",
            expected_edits=False, alignment_timeout=600, max_locus_bases=100000,
            exon_identity=.7, anchor_bases=12, anchor_identity=.8, threads=1,
            flank=1000, max_extension=10000,
        )
        with patch("intraphy.case.build_case") as build_case, \
                patch("intraphy.inference.genomic_exon_run.infer_genomic_exons", return_value="done") as infer:
            self.assertEqual(dispatch_analyze(args), "done")
        build_case.assert_called_once()
        selection.write.assert_called_once_with(Path("result") / "prepared_inputs")
        self.assertFalse(build_case.call_args.kwargs["derive_correspondence"])
        self.assertTrue(build_case.call_args.kwargs["allow_unannotated"])
        infer.assert_called_once()

    def test_required_tools_follow_catalogue_reuse(self):
        parser = build_parser()
        raw = parser.parse_args(["analyze", "--fasta", "genomes", "--gff", "annotations",
                                 "--species-tree", "tree.nwk", "--output-dir", "result"])
        reused = parser.parse_args(["analyze", "--input-dir", "prepared",
                                    "--exon-configurations", "catalogue.jsonl",
                                    "--output-dir", "result"])
        self.assertEqual(set(required_tools(raw)), {"mafft", "minimap2"})
        self.assertEqual(required_tools(reused), [])

    def test_configuration_reuse_requires_the_genomic_exon_marker_and_prepared_tree(self):
        with tempfile.TemporaryDirectory() as directory:
            catalogue = Path(directory) / "catalogue.jsonl"
            catalogue.write_text(json.dumps({"observation_unit": "transcript_configuration"}) + "\n",
                                 encoding="utf-8")
            args = SimpleNamespace(command="analyze", model="exon-structure-ctmc",
                                   exon_configurations=str(catalogue), input_dir="prepared",
                                   species_tree=None, parameter_mode="fit", exon_rates=None,
                                   root_frequency=None, root_presence=.5,
                                   _root_frequency_explicit=False, _root_presence_explicit=False,
                                   dna_gain_rate=None, dna_loss_rate=None, observation_view="evidence",
                                   max_observation_scenarios=None, survey_min_identity=.7,
                                   survey_min_coverage=.8, survey_max_dp_cells=250000,
                                   intron_anchor_window=15, intron_min_anchor_pairs=8,
                                   origin_root_sensitivity=[], branch_length_mode="supplied",
                                   alignment_timeout=600, max_locus_bases=100000, anchor_bases=12,
                                   max_states=None, max_origin_scenarios=None, exon_identity=.7,
                                   anchor_identity=.8)
            with self.assertRaisesRegex(ValueError, "genomic_exon_spans"):
                validate_analysis_options(args)

    def test_fixed_rate_file_requires_fixed_parameter_mode(self):
        args = SimpleNamespace(command="analyze", model="exon-structure-ctmc",
                               locus_model=None, genomic_evidence_dir=None,
                               input_dir="prepared", fasta=None, manifest=None, gff=None,
                               orthologs=None, _root_frequency_explicit=False,
                               _root_presence_explicit=False, dna_gain_rate=None,
                               dna_loss_rate=None, exon_rates="rates.json",
                               parameter_mode="fit", root_presence=.5, observation_view="evidence",
                               max_observation_scenarios=None, survey_min_identity=.7,
                               survey_min_coverage=.8, survey_max_dp_cells=250000,
                               intron_anchor_window=15, intron_min_anchor_pairs=8,
                               origin_root_sensitivity=[], branch_length_mode="supplied",
                               alignment_timeout=600, max_locus_bases=100000, anchor_bases=12,
                               max_states=None, max_origin_scenarios=None, exon_identity=.7,
                               anchor_identity=.8, species_tree=None)
        with self.assertRaisesRegex(ValueError, "--parameter-mode fixed"):
            validate_analysis_options(args)

    def test_new_model_rejects_explicit_binary_root_controls(self):
        parser = build_parser()
        for flag in (("--root-frequency", "estimated"), ("--root-presence", ".5")):
            with self.subTest(flag=flag):
                args = parser.parse_args(["analyze", "--fasta", "genomes", "--gff", "annotations",
                                          "--species-tree", "tree.nwk", "--output-dir", "result", *flag])
                with self.assertRaisesRegex(ValueError, "root-frequency controls"):
                    validate_analysis_options(args)

    def test_prepared_tree_replacement_uses_exact_tip_panel_and_selected_branch_lengths(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            prepared = root / "prepared"
            prepared.mkdir()
            (prepared / "species_tree.tsv").write_text(
                "node_id\tparent_id\tlabel\tbranch_length\n"
                "root\t\troot\t\nA\troot\tA\t0.1\nB\troot\tB\t0.2\n",
                encoding="utf-8")
            override = root / "override.nwk"
            override.write_text("(A:0.3,B:0.4)R;", encoding="utf-8")
            args = SimpleNamespace(input_dir=str(prepared), exon_configurations="catalogue.jsonl",
                                   species_tree=str(override), branch_length_mode="supplied")
            validate_input_paths(args)
            tree = args._species_tree
            self.assertEqual(tree.length[tree.leaf_by_label["A"]], .3)
            self.assertEqual(tree.length[tree.leaf_by_label["B"]], .4)
            mismatch = root / "mismatch.nwk"
            mismatch.write_text("(A:0.3,C:0.4)R;", encoding="utf-8")
            args.species_tree = str(mismatch)
            with self.assertRaisesRegex(ValueError, "exactly the prepared tree tip panel"):
                validate_input_paths(args)

    def test_unresolved_completed_genomic_exon_result_is_recognized(self):
        with tempfile.TemporaryDirectory() as directory:
            RunResult("exon-structure-ctmc", "single-copy", "species_tree.tsv", (),
                      status="completed_with_unresolved").write(directory)
            self.assertEqual(result_model(directory), "exon-structure-ctmc")


if __name__ == "__main__":
    unittest.main()
