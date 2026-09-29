"""Source-level CLI checks for the automatic genomic presence route."""
import tempfile
import unittest
import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from intraphy.commands.genomic import dispatch_analyze, dispatch_prepare_evidence
from intraphy.commands.parser import build_parser
from intraphy.commands.preflight import validate_arguments, validate_input_paths
from intraphy.cli import main
from intraphy.storage.tabular import read_tsv, write_tsv
from support_intron_fixtures import spliced_fixture


def _occurrence(occurrence_id, species, start, end, *, gene):
    return {"occurrence_id": occurrence_id, "family_id": "fam", "species": species,
            "gene_copy_id": gene, "transcript_id": f"tx_{gene}", "role": "exon",
            "presence_status": "present", "contig": "chr1", "start": start,
            "end": end, "strand": "+", "source_feature_id": f"feature_{occurrence_id}"}


def _edge(left, right):
    return {"query_occurrence_id": left, "subject_occurrence_id": right,
            "position_edge_eligible": "1", "membership_edge_eligible": "1",
            "match_status": "mapped", "candidate_resolution": "resolved",
            "query_genomic_matched_blocks": "chr1:1-4:+",
            "subject_genomic_matched_blocks": "chr1:1-4:+"}


class GenomicCommandTests(unittest.TestCase):
    def _prepare_tree(self, prepared, character_type, tree_text):
        (prepared / "species_tree.tsv").write_text(tree_text, encoding="utf-8")
        args = SimpleNamespace(command="prepare-genomic-evidence", input_dir=str(prepared),
                               character_type=character_type)
        module = ("intraphy.commands.intron_preflight.validate_prepared_input"
                  if character_type == "intron-position"
                  else "intraphy.commands.genomic_preflight.validate_prepared_input")
        with patch(module) as validate_prepared:
            validate_input_paths(args)
        return validate_prepared.call_args.args[1]

    def test_evidence_preparation_accepts_topology_without_branch_lengths(self):
        with tempfile.TemporaryDirectory() as directory:
            prepared = Path(directory)
            topologies = (
                ("blank", "node_id\tparent_id\tlabel\tbranch_length\n"
                 "root\t\troot\t\nA\troot\tA\t\nB\troot\tB\t\n"),
                ("NA", "node_id\tparent_id\tlabel\tbranch_length\n"
                 "root\t\troot\t\nA\troot\tA\tNA\nB\troot\tB\t\n"),
                ("no-column", "node_id\tparent_id\tlabel\n"
                 "root\t\troot\nA\troot\tA\nB\troot\tB\n"),
            )
            for branch_value, tree_text in topologies:
                for character_type in ("dna-presence", "intron-position"):
                    with self.subTest(branch_value=branch_value, character_type=character_type):
                        tree = self._prepare_tree(prepared, character_type, tree_text)
                        self.assertEqual(set(tree.leaf_by_label), {"A", "B"})
                        self.assertIsNone(tree.length["A"])
                        self.assertIsNone(tree.length["B"])

    def test_evidence_preparation_preserves_supplied_lengths_and_rejects_invalid_tree(self):
        with tempfile.TemporaryDirectory() as directory:
            prepared = Path(directory)
            complete = ("node_id\tparent_id\tlabel\tbranch_length\n"
                        "root\t\troot\t\nA\troot\tA\t0.2\nB\troot\tB\t0.3\n")
            for character_type in ("dna-presence", "intron-position"):
                with self.subTest(character_type=character_type):
                    tree = self._prepare_tree(prepared, character_type, complete)
                    self.assertEqual(tree.length["A"], 0.2)
                    self.assertEqual(tree.length["B"], 0.3)

            invalid_trees = (
                ("negative branch length", "node_id\tparent_id\tlabel\tbranch_length\n"
                 "root\t\troot\t\nA\troot\tA\t-0.2\n"),
                ("missing parent", "node_id\tparent_id\tlabel\tbranch_length\n"
                 "root\t\troot\t\nA\tmissing\tA\tNA\n"),
            )
            for message, tree_text in invalid_trees:
                for character_type in ("dna-presence", "intron-position"):
                    module = ("intraphy.commands.intron_preflight.validate_prepared_input"
                              if character_type == "intron-position"
                              else "intraphy.commands.genomic_preflight.validate_prepared_input")
                    with self.subTest(message=message, character_type=character_type), patch(module):
                        (prepared / "species_tree.tsv").write_text(tree_text, encoding="utf-8")
                        args = SimpleNamespace(command="prepare-genomic-evidence",
                                               input_dir=str(prepared), character_type=character_type)
                        with self.assertRaisesRegex(SystemExit, message):
                            validate_input_paths(args)

    def test_cli_evidence_preparation_accepts_topology_and_rejects_roster_mismatch(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            dna_root = root / "dna"
            dna_root.mkdir()
            dna_inputs = self._write_prepared_fixture(dna_root)
            intron_root = root / "intron"
            intron_root.mkdir()
            intron_inputs = spliced_fixture(intron_root)
            topology = ("node_id\tparent_id\tlabel\tbranch_length\n"
                        "root\t\troot\t\nA\troot\tA\tNA\nB\troot\tB\t\n")
            single_tip = ("node_id\tparent_id\tlabel\tbranch_length\n"
                          "root\t\troot\t\nA\troot\tA\tNA\n")
            cases = (("dna-presence", dna_inputs, "dna_observations.tsv"),
                     ("intron-position", intron_inputs, "intron_observations.tsv"))
            for character_type, prepared, observation_name in cases:
                with self.subTest(character_type=character_type):
                    tree_path = prepared / "species_tree.tsv"
                    tree_path.write_text(topology, encoding="utf-8")
                    evidence = root / f"{character_type}_evidence"
                    command = ["prepare-genomic-evidence", "--character-type", character_type,
                               "--input-dir", str(prepared), "--output-dir", str(evidence)]
                    if character_type == "intron-position":
                        with patch("intraphy.commands.preflight.shutil.which", return_value="mafft"), \
                                patch("intraphy.coding_correspondence.alignment.protein_multiple_alignment",
                                      side_effect=lambda records, mode="linsi", threads=1: dict(records)):
                            self.assertEqual(main(command), 0)
                    else:
                        self.assertEqual(main(command), 0)
                    self.assertTrue((evidence / observation_name).is_file())

                    tree_path.write_text(single_tip, encoding="utf-8")
                    rejected_output = root / f"{character_type}_rejected"
                    rejected_command = ["prepare-genomic-evidence", "--character-type", character_type,
                                        "--input-dir", str(prepared), "--output-dir", str(rejected_output)]
                    with patch("intraphy.commands.preflight.shutil.which", return_value="mafft"):
                        self.assertEqual(main(rejected_command), 2)
                    self.assertFalse(rejected_output.exists())

    def test_phylogenetic_inference_still_requires_nonroot_branch_lengths(self):
        with tempfile.TemporaryDirectory() as directory:
            prepared = Path(directory)
            (prepared / "species_tree.tsv").write_text(
                "node_id\tparent_id\tlabel\tbranch_length\n"
                "root\t\troot\t\nA\troot\tA\tNA\nB\troot\tB\t0.3\n",
                encoding="utf-8",
            )
            for model in ("dna-presence-ctmc", "intron-position-ctmc"):
                with self.subTest(model=model):
                    args = SimpleNamespace(command="analyze", model=model,
                                           input_dir=str(prepared), species_tree=None)
                    with self.assertRaisesRegex(SystemExit, "lacks a branch length"):
                        validate_input_paths(args)

    def test_raw_genomic_analyze_defaults_to_exon_structure(self):
        args = build_parser().parse_args([
            "analyze", "--fasta", "genomes", "--gff", "annotations",
            "--species-tree", "tree.nwk", "--output-dir", "result",
        ])
        self.assertEqual(args.model, "exon-structure-ctmc")
        validate_arguments(args)

    def test_locus_input_selects_the_explicit_advanced_model(self):
        args = build_parser().parse_args([
            "analyze", "--locus-model", "model.json", "--species-tree", "tree.nwk",
            "--output-dir", "result",
        ])
        self.assertEqual(args.model, "exon-locus-ctmc")

    def test_fixed_rates_are_paired_and_root_modes_are_explicit(self):
        parser = build_parser()
        args = parser.parse_args([
            "analyze", "--model", "dna-presence-ctmc", "--fasta", "genomes", "--gff", "annotations",
            "--species-tree", "tree.nwk", "--output-dir", "result",
            "--dna-gain-rate", "0.1",
        ])
        with self.assertRaisesRegex(ValueError, "must be supplied together"):
            validate_arguments(args)
        args = parser.parse_args([
            "analyze", "--model", "dna-presence-ctmc", "--fasta", "genomes", "--gff", "annotations",
            "--species-tree", "tree.nwk", "--output-dir", "result",
            "--root-frequency", "fixed", "--root-presence", "0.25",
            "--dna-gain-rate", "0.1", "--dna-loss-rate", "0.2",
        ])
        validate_arguments(args)
        for endpoint in ("0", "1"):
            edge_args = parser.parse_args([
                "analyze", "--model", "dna-presence-ctmc", "--fasta", "genomes", "--gff", "annotations",
                "--species-tree", "tree.nwk", "--output-dir", "result",
                "--root-frequency", "fixed", "--root-presence", endpoint,
            ])
            validate_arguments(edge_args)

    def test_prepared_evidence_command_forwards_survey_parameters(self):
        args = SimpleNamespace(input_dir="prepared", output_dir="evidence", threads=3,
                               survey_min_identity=.75, survey_min_coverage=.85,
                               survey_max_dp_cells=12345)
        with patch("intraphy.commands.genomic.prepare_evidence", return_value=[]) as prepare:
            dispatch_prepare_evidence(args)
        prepare.assert_called_once_with("prepared", "evidence", min_identity=.75,
                                        min_coverage=.85, max_dp_cells=12345, threads=3)

    def test_cli_prepared_evidence_streams_large_segment_match_description(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            prepared = self._write_prepared_fixture(root)
            baseline = root / "baseline_evidence"
            command = ["prepare-genomic-evidence", "--input-dir", str(prepared),
                       "--output-dir", str(baseline)]
            self.assertEqual(main(command), 0)
            expected = read_tsv(baseline / "dna_observations.tsv")

            matches_path = prepared / "segment_matches.tsv"
            matches = read_tsv(matches_path)
            # This explanatory field is part of the prepared match schema and
            # does not affect physical-site eligibility.
            matches[0]["correspondence_basis"] = ('evidence\t"quoted\nline" λ ' * 9000)
            self.assertGreater(len(matches[0]["correspondence_basis"]), 131072)
            write_tsv(matches_path, matches, list(matches[0]))

            actual_dir = root / "large_field_evidence"
            command[command.index(str(baseline))] = str(actual_dir)
            self.assertEqual(main(command), 0)
            self.assertEqual(read_tsv(actual_dir / "dna_observations.tsv"), expected)

    def test_dispatch_uses_fixed_supplied_tree_unit_and_does_not_invent_source_events(self):
        with tempfile.TemporaryDirectory() as temp:
            args = SimpleNamespace(
                output_dir=temp, input_dir="prepared", genomic_evidence_dir="evidence",
                _species_tree=object(), dna_gain_rate=.1, dna_loss_rate=.2,
                root_frequency="fixed", root_presence=.25, expected_edits=True, threads=4,
                survey_min_identity=.7, survey_min_coverage=.8, survey_max_dp_cells=250000,
            )
            with patch("intraphy.commands.genomic.read_tsv", return_value=[{"state": "1"}]), \
                    patch("intraphy.inference.dna_presence.analyze_dna_presence") as analyze:
                dispatch_analyze(args)
            analyze.assert_called_once_with(
                [{"state": "1"}], args._species_tree, Path(temp),
                root_frequency="fixed", root_presence=.25,
                fixed_rates={"gain": .1, "loss": .2},
                branch_length_unit="supplied_tree_units", expected_counts=True, workers=4,
            )

    def _write_prepared_fixture(self, root, *, gene_only_b=False):
        prepared = root / "prepared"
        prepared.mkdir()
        source_sequence = "ACGT" + "TT" + "AAA" + "GGG" + "TGCA"
        target_sequence = "ACGT" + "TT" + "GGG" + "TGCA"
        source_fasta, target_fasta = root / "A.fa", root / "B.fa"
        source_fasta.write_text(f">chr1\n{source_sequence}\n", encoding="utf-8")
        target_fasta.write_text(f">chr1\n{target_sequence}\n", encoding="utf-8")
        occurrences = [
            _occurrence("a_left", "A", 1, 4, gene="geneA"),
            _occurrence("a_focal", "A", 7, 9, gene="geneA"),
            _occurrence("a_right", "A", 13, 16, gene="geneA"),
            _occurrence("b_left", "B", 1, 4, gene="geneB"),
            _occurrence("b_right", "B", 10, 13, gene="geneB"),
        ]
        if gene_only_b:
            occurrences = [_occurrence("a_focal", "A", 7, 9, gene="geneA")]
        write_tsv(prepared / "segment_occurrences.tsv", occurrences, list(occurrences[0]))
        edges = [] if gene_only_b else [_edge("a_left", "b_left"), _edge("a_right", "b_right")]
        write_tsv(prepared / "segment_matches.tsv", edges, list(edges[0]) if edges else list(_edge("q", "s")))
        loci = [
            {"species": "A", "gene_copy_id": "geneA", "contig": "chr1", "strand": "+",
             "search_start": 1, "search_end": len(source_sequence), "genome_fasta": str(source_fasta.resolve())},
            {"species": "B", "gene_copy_id": "geneB", "contig": "chr1", "strand": "+",
             "search_start": 1, "search_end": len(target_sequence), "genome_fasta": str(target_fasta.resolve())},
        ]
        write_tsv(prepared / "gene_loci.tsv", loci, list(loci[0]))
        targets = [{"family_id": "fam", "species": "A", "gene_copy_id": "geneA"},
                   {"family_id": "fam", "species": "B", "gene_copy_id": "geneB"}]
        write_tsv(prepared / "input_targets.tsv", targets, list(targets[0]))
        write_tsv(prepared / "segment_homology.tsv", [], ["homology_id", "occurrence_id"])
        (prepared / "segment_sequences.fasta").write_text(">fixture\nA\n", encoding="utf-8")
        (prepared / "species_tree.tsv").write_text(
            "node_id\tparent_id\tlabel\tbranch_length\n"
            "root\t\troot\t\nA\troot\tA\t0.2\nB\troot\tB\t0.3\n", encoding="utf-8")
        return prepared

    def test_cli_prepared_evidence_fit_and_staged_member_mismatch(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            prepared = self._write_prepared_fixture(root)
            evidence, result = root / "evidence", root / "result"
            self.assertEqual(main(["prepare-genomic-evidence", "--input-dir", str(prepared),
                                   "--output-dir", str(evidence)]), 0)
            self.assertEqual(main(["analyze", "--model", "dna-presence-ctmc", "--input-dir", str(prepared),
                                   "--genomic-evidence-dir", str(evidence),
                                   "--output-dir", str(result), "--root-frequency", "fixed",
                                   "--root-presence", ".5", "--dna-gain-rate", ".1",
                                   "--dna-loss-rate", ".1"]), 0)
            run = json.loads((result / "run_result.json").read_text(encoding="utf-8"))
            history = json.loads((result / "dna_history.json").read_text(encoding="utf-8"))
            self.assertEqual(run["model"], "dna-presence-ctmc")
            self.assertEqual(run["tree_file"], "species_tree.tsv")
            focal = next(row for row in read_tsv(evidence / "genomic_members.tsv")
                         if row["occurrence_ids"] == "a_focal")
            focal_history = next(site for site in history["families"][0]["sites"]
                                 if site["site_id"] == focal["site_id"])
            tip_probabilities = {row["label"]: row["p_present"] for row in focal_history["nodes"]
                                 if row["label"] in {"A", "B"}}
            self.assertEqual(tip_probabilities, {"A": 1.0, "B": 0.0})
            member_rows = read_tsv(evidence / "genomic_members.tsv")
            write_tsv(evidence / "genomic_members.tsv", member_rows[:-1], list(member_rows[0]))
            self.assertEqual(main(["analyze", "--model", "dna-presence-ctmc", "--input-dir", str(prepared),
                                   "--genomic-evidence-dir", str(evidence),
                                   "--output-dir", str(root / "mismatch_result"),
                                   "--dna-gain-rate", ".1", "--dna-loss-rate", ".1"]), 2)

    def test_gene_only_species_is_unknown_and_fit_reports_no_contrast(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            prepared = self._write_prepared_fixture(root, gene_only_b=True)
            evidence, result = root / "gene_only_evidence", root / "gene_only_result"
            self.assertEqual(main(["prepare-genomic-evidence", "--input-dir", str(prepared),
                                   "--output-dir", str(evidence)]), 0)
            observations = read_tsv(evidence / "dna_observations.tsv")
            self.assertEqual({row["species"]: row["state"] for row in observations},
                             {"A": "1", "B": "unknown"})
            self.assertEqual(main(["analyze", "--model", "dna-presence-ctmc", "--input-dir", str(prepared),
                                   "--genomic-evidence-dir", str(evidence),
                                   "--output-dir", str(result)]), 0)
            fit = json.loads((result / "dna_fit.json").read_text(encoding="utf-8"))
            run = json.loads((result / "run_result.json").read_text(encoding="utf-8"))
            self.assertEqual(fit["families"][0]["status"], "no_observed_contrast")
            self.assertEqual(run["status"], "completed_with_unresolved")


if __name__ == "__main__":
    unittest.main()
