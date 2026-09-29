"""Tree import, pruning, and prepared genomic-analysis replacement contracts."""
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from intraphy.commands.preflight import validate_input_paths
from intraphy.cli import main
from intraphy.inference.dna_output import write_binary_presence_outputs
from intraphy.inference.binary_presence import DNA_DOMAIN, INTRON_DOMAIN
from intraphy.inputs.species_tree import (
    prune_species_tree,
    read_species_tree_rows,
    select_species_tree_rows,
)
from intraphy.orthofinder import _write_species_tree
from intraphy.storage.tabular import read_tsv, write_tsv
from intraphy.topology import SpeciesTree
from support_intron_fixtures import spliced_fixture


NEWICK = "((A:1,B:2)X:3,C:4)R:0;"


class SpeciesTreeInputTests(unittest.TestCase):
    def test_newick_subset_preserves_root_unary_path_and_distances(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "tree.nwk"
            source.write_text(NEWICK, encoding="utf-8")
            rows = select_species_tree_rows(source, ["A", "C"], prune=True)
            tree = SpeciesTree(rows)

        self.assertEqual([row["node_id"] for row in rows],
                         ["node_1", "node_2", "node_3", "node_5"])
        self.assertEqual(tree.root, "node_1")
        self.assertEqual(tree.length[tree.root], 0.0)
        self.assertEqual(tree.label[tree.root], "R")
        self.assertEqual(tree.label["node_2"], "X")
        self.assertEqual(tree.children["node_2"], ["node_3"])
        self.assertEqual(tree.length["node_2"], 3.0)
        self.assertEqual(tree.length["node_3"], 1.0)
        self.assertEqual(tree.length["node_5"], 4.0)
        self.assertEqual(tree.branch_length("node_2") + tree.branch_length("node_3"), 4.0)
        self.assertEqual(tree.branch_length("node_5"), 4.0)
        self.assertEqual(tree.branch_length("node_2") + tree.branch_length("node_3")
                         + tree.branch_length("node_5"), 8.0)

    def test_tsv_subset_preserves_root_and_original_branch_lengths(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "tree.tsv"
            source.write_text(
                "node_id\tparent_id\tlabel\tlength\n"
                "r\t\tR\t0\nix\tr\tX\t3\na\tix\tA\t1\n"
                "b\tix\tB\t2\nc\tr\tC\t4\n", encoding="utf-8")
            output = Path(directory) / "selected.tsv"
            prune_species_tree(source, output, ["A", "C"], prune_extra_tips=True)
            rows = read_tsv(output)
            tree = SpeciesTree(rows)

        self.assertEqual([row["node_id"] for row in rows], ["r", "ix", "a", "c"])
        self.assertEqual(tree.length["r"], 0.0)
        self.assertEqual(tree.children["ix"], ["a"])
        self.assertEqual(tree.length["ix"], 3.0)
        self.assertEqual(tree.length["a"], 1.0)
        self.assertEqual(tree.length["c"], 4.0)
        self.assertEqual(tree.branch_length("ix") + tree.branch_length("a")
                         + tree.branch_length("c"), 8.0)

    def test_panel_mismatch_duplicate_and_empty_panel_are_explicit(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "tree.nwk"
            source.write_text(NEWICK, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "explicit pruning"):
                select_species_tree_rows(source, ["A", "C"])
            pruned = select_species_tree_rows(source, ["A", "C"], prune=True)
            self.assertEqual({row["label"] for row in pruned if row["node_id"] in
                              {"node_3", "node_5"}}, {"A", "C"})
            with self.assertRaisesRegex(ValueError, "missing requested tips"):
                select_species_tree_rows(source, ["A", "missing"], prune=True)
            with self.assertRaisesRegex(ValueError, "duplicate"):
                select_species_tree_rows(source, ["A", "A"], prune=True)
            with self.assertRaisesRegex(ValueError, "must not be empty"):
                select_species_tree_rows(source, [], prune=True)

    def test_invalid_nonroot_lengths_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "tree.tsv"
            for value in ("-1", "inf", "-inf", "nan"):
                with self.subTest(value=value):
                    source.write_text(
                        "node_id\tparent_id\tlabel\tbranch_length\n"
                        f"r\t\tR\t0\na\tr\tA\t{value}\n", encoding="utf-8")
                    with self.assertRaises(SystemExit):
                        read_species_tree_rows(source)

    def test_legacy_newick_conversion_keeps_native_node_ids_and_root_length(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "tree.nwk"
            output = Path(directory) / "species_tree.tsv"
            source.write_text("((A:1,B)X:3,C)R;", encoding="utf-8")
            _write_species_tree(source, output)
            rows = read_tsv(output)

        self.assertEqual([row["node_id"] for row in rows],
                         ["node_1", "node_2", "node_3", "node_4", "node_5"])
        self.assertEqual(rows[0]["branch_length"], "0.0")
        self.assertEqual(rows[1]["branch_length"], "3.0")
        self.assertEqual(rows[2]["branch_length"], "1.0")
        self.assertEqual(rows[3]["branch_length"], "NA")
        self.assertEqual(rows[4]["branch_length"], "NA")

    def test_prepared_na_tree_uses_exact_panel_override_after_evidence_validation(self):
        for model in ("dna-presence-ctmc", "intron-position-ctmc"):
            with self.subTest(model=model), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                prepared = _prepared_case(root, "NA")
                override = root / "override.nwk"
                override.write_text("(A:0.5,B:0.75)R:0;", encoding="utf-8")
                evidence = root / "evidence"
                evidence.mkdir()
                observation_file = "dna_observations.tsv" if model.startswith("dna") else "intron_observations.tsv"
                (evidence / observation_file).write_text("species\tstate\nA\t1\n", encoding="utf-8")
                args = _analysis_args(model, prepared, override, evidence)
                validation = ("intraphy.commands.genomic_preflight.validate_reusable_evidence"
                              if model.startswith("dna") else
                              "intraphy.commands.intron_preflight.validate_reusable_evidence")
                prepared_validation = ("intraphy.commands.genomic_preflight.validate_prepared_input"
                                       if model.startswith("dna") else
                                       "intraphy.commands.intron_preflight.validate_prepared_input")
                with patch(prepared_validation) as validate_case, patch(validation) as validate_evidence:
                    validate_input_paths(args)

                checked_case_tree = validate_case.call_args.args[1]
                checked_evidence_tree = validate_evidence.call_args.args[2]
                self.assertEqual(set(checked_case_tree.leaf_by_label), {"A", "B"})
                self.assertIsNone(checked_case_tree.length[checked_case_tree.leaf_by_label["A"]])
                self.assertIs(checked_evidence_tree, checked_case_tree)
                self.assertEqual(args._species_tree.length[args._species_tree.leaf_by_label["A"]], 0.5)
                self.assertEqual(args._species_tree.length[args._species_tree.leaf_by_label["B"]], 0.75)

    def test_na_prepared_tree_without_override_still_fails_inference_preflight(self):
        for model in ("dna-presence-ctmc", "intron-position-ctmc"):
            with self.subTest(model=model), tempfile.TemporaryDirectory() as directory:
                prepared = _prepared_case(Path(directory), "NA")
                args = _analysis_args(model, prepared, None, None)
                prepared_validation = ("intraphy.commands.genomic_preflight.validate_prepared_input"
                                      if model.startswith("dna") else
                                      "intraphy.commands.intron_preflight.validate_prepared_input")
                with patch(prepared_validation):
                    with self.assertRaisesRegex(SystemExit, "lacks a branch length"):
                        validate_input_paths(args)

    def test_override_tip_mismatch_and_invalid_reusable_evidence_are_rejected(self):
        for model in ("dna-presence-ctmc", "intron-position-ctmc"):
            for tips in ("(A:1,C:1)R;", "(A:1,B:1,C:1)R;"):
                with self.subTest(model=model, tips=tips), tempfile.TemporaryDirectory() as directory:
                    root = Path(directory)
                    prepared = _prepared_case(root, "NA")
                    override = root / "override.nwk"
                    override.write_text(tips, encoding="utf-8")
                    args = _analysis_args(model, prepared, override, None)
                    prepared_validation = ("intraphy.commands.genomic_preflight.validate_prepared_input"
                                          if model.startswith("dna") else
                                          "intraphy.commands.intron_preflight.validate_prepared_input")
                    with patch(prepared_validation):
                        with self.assertRaisesRegex(ValueError, "exactly the prepared tree tip panel"):
                            validate_input_paths(args)

            with self.subTest(model=model, evidence="invalid"), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                prepared = _prepared_case(root, "NA")
                override = root / "override.nwk"
                override.write_text("(A:1,B:1)R;", encoding="utf-8")
                evidence = root / "evidence"
                evidence.mkdir()
                observation_file = "dna_observations.tsv" if model.startswith("dna") else "intron_observations.tsv"
                (evidence / observation_file).write_text("species\tstate\nA\t1\n", encoding="utf-8")
                args = _analysis_args(model, prepared, override, evidence)
                prepared_validation = ("intraphy.commands.genomic_preflight.validate_prepared_input"
                                      if model.startswith("dna") else
                                      "intraphy.commands.intron_preflight.validate_prepared_input")
                evidence_validation = ("intraphy.commands.genomic_preflight.validate_reusable_evidence"
                                      if model.startswith("dna") else
                                      "intraphy.commands.intron_preflight.validate_reusable_evidence")
                with patch(prepared_validation), patch(evidence_validation, side_effect=ValueError("bad staged evidence")):
                    with self.assertRaisesRegex(ValueError, "bad staged evidence"):
                        validate_input_paths(args)
                self.assertIsNone(args._species_tree.length[args._species_tree.leaf_by_label["A"]])

    def test_raw_input_route_without_override_has_no_unbound_tree_state(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tree = root / "tree.nwk"
            tree.write_text("(A:1,B:1)R;", encoding="utf-8")
            args = SimpleNamespace(command="analyze", model="dna-presence-ctmc",
                                   input_dir=None, species_tree=str(tree), fasta=["genomes"],
                                   manifest=None, genomic_evidence_dir=None)
            selection = SimpleNamespace(tree_species=("A", "B"))
            with patch("intraphy.inputs.selection.resolve_inputs", return_value=selection):
                validate_input_paths(args)
            self.assertEqual(set(args._species_tree.leaf_by_label), {"A", "B"})

    def test_cli_output_uses_override_without_modifying_prepared_tree(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            prepared = _prepared_case(root, "NA")
            original = (prepared / "species_tree.tsv").read_text(encoding="utf-8")
            override = root / "override.nwk"
            override.write_text("(A:0.5,B:0.75)R:0;", encoding="utf-8")
            output = root / "result"

            def write_using_selected_tree(rows, tree, output_dir, **kwargs):
                write_binary_presence_outputs(output_dir, tree, [], [], DNA_DOMAIN)

            with patch("intraphy.commands.genomic.prepare_evidence", return_value=[]), \
                    patch("intraphy.inference.dna_presence.analyze_dna_presence",
                          side_effect=write_using_selected_tree):
                self.assertEqual(main(["analyze", "--model", "dna-presence-ctmc", "--input-dir", str(prepared),
                                       "--species-tree", str(override), "--output-dir", str(output)]), 0)

            written = SpeciesTree(read_tsv(output / "species_tree.tsv"))
            self.assertEqual(written.length[written.leaf_by_label["A"]], 0.5)
            self.assertEqual(written.length[written.leaf_by_label["B"]], 0.75)
            self.assertEqual((prepared / "species_tree.tsv").read_text(encoding="utf-8"), original)

    def test_cli_reuses_real_dna_and_intron_evidence_with_replacement_tree(self):
        for model in ("dna-presence-ctmc", "intron-position-ctmc"):
            with self.subTest(model=model), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                if model == "dna-presence-ctmc":
                    prepared = _dna_prepared_fixture(root)
                    character_type = "dna-presence"
                    observation_name = "dna_observations.tsv"
                else:
                    prepared = spliced_fixture(root)
                    character_type = "intron-position"
                    observation_name = "intron_observations.tsv"
                original_tree = (
                    "node_id\tparent_id\tlabel\tbranch_length\n"
                    "root\t\troot\t\nA\troot\tA\tNA\nB\troot\tB\tNA\n")
                (prepared / "species_tree.tsv").write_text(original_tree, encoding="utf-8")
                evidence = root / "evidence"
                prep_command = ["prepare-genomic-evidence", "--character-type", character_type,
                                "--input-dir", str(prepared), "--output-dir", str(evidence)]
                if model == "intron-position-ctmc":
                    with patch("intraphy.commands.preflight.shutil.which", return_value="mafft"), \
                            patch("intraphy.coding_correspondence.alignment.protein_multiple_alignment",
                                  side_effect=lambda records, mode="linsi", threads=1: dict(records)):
                        self.assertEqual(main(prep_command), 0)
                else:
                    self.assertEqual(main(prep_command), 0)

                observation_path = evidence / observation_name
                original_observations = observation_path.read_bytes()
                override = root / "override.nwk"
                override.write_text("(A:0.5,B:0.75)R:0;", encoding="utf-8")
                output = root / "accepted"
                args = ["analyze", "--model", model, "--input-dir", str(prepared),
                        "--species-tree", str(override), "--genomic-evidence-dir", str(evidence),
                        "--output-dir", str(output), "--root-frequency", "fixed",
                        "--root-presence", ".5"]
                if model == "dna-presence-ctmc":
                    args += ["--dna-gain-rate", ".1", "--dna-loss-rate", ".2"]
                    analyze_path = "intraphy.inference.dna_presence.analyze_dna_presence"
                    domain = DNA_DOMAIN
                else:
                    args += ["--gain-rate", ".1", "--loss-rate", ".2"]
                    analyze_path = "intraphy.inference.intron_presence.analyze_intron_positions"
                    domain = INTRON_DOMAIN

                def write_selected_tree(rows, tree, output_dir, **kwargs):
                    write_binary_presence_outputs(output_dir, tree, [], [], domain)

                with patch(analyze_path, side_effect=write_selected_tree):
                    self.assertEqual(main(args), 0)
                result_tree = SpeciesTree(read_tsv(output / "species_tree.tsv"))
                self.assertEqual(result_tree.length[result_tree.leaf_by_label["A"]], 0.5)
                self.assertEqual(result_tree.length[result_tree.leaf_by_label["B"]], 0.75)
                self.assertEqual((prepared / "species_tree.tsv").read_text(encoding="utf-8"), original_tree)
                self.assertEqual(observation_path.read_bytes(), original_observations)

                if model == "dna-presence-ctmc":
                    observations = read_tsv(observation_path)
                    observations[0]["min_identity"] = "0.71"
                    write_tsv(observation_path, observations, list(observations[0]))
                else:
                    summaries_path = evidence / "intron_families.tsv"
                    summaries = read_tsv(summaries_path)
                    summaries[0]["anchor_window"] = "14"
                    write_tsv(summaries_path, summaries, list(summaries[0]))
                rejected = root / "rejected"
                args[args.index(str(output))] = str(rejected)
                with patch(analyze_path) as analyze:
                    self.assertEqual(main(args), 2)
                analyze.assert_not_called()


def _prepared_case(root, branch_value):
    prepared = Path(root) / "prepared"
    prepared.mkdir()
    (prepared / "species_tree.tsv").write_text(
        "node_id\tparent_id\tlabel\tbranch_length\n"
        f"root\t\troot\t\nA\troot\tA\t{branch_value}\nB\troot\tB\t{branch_value}\n",
        encoding="utf-8")
    targets = [
        {"family_id": "fam", "species": "A", "gene_copy_id": "copy_A"},
        {"family_id": "fam", "species": "B", "gene_copy_id": "copy_B"},
    ]
    write_tsv(prepared / "input_targets.tsv", targets, list(targets[0]))
    (prepared / "gene_loci.tsv").write_text("species\tgene_copy_id\n", encoding="utf-8")
    occurrences = ["occurrence_id", "family_id", "species", "gene_copy_id", "role",
                   "presence_status", "contig", "start", "end", "strand", "transcript_id",
                   "source_feature_id"]
    (prepared / "segment_occurrences.tsv").write_text("\t".join(occurrences) + "\n", encoding="utf-8")
    (prepared / "segment_homology.tsv").write_text("homology_id\toccurrence_id\n", encoding="utf-8")
    (prepared / "segment_matches.tsv").write_text("query_occurrence_id\tsubject_occurrence_id\n", encoding="utf-8")
    (prepared / "segment_sequences.fasta").write_text(">empty\nA\n", encoding="utf-8")
    return prepared


def _dna_prepared_fixture(root):
    prepared = Path(root) / "prepared"
    prepared.mkdir()
    seq_a = "ACGT" + "TT" + "AAA" + "GGG" + "TGCA"
    seq_b = "ACGT" + "TT" + "GGG" + "TGCA"
    fasta_a, fasta_b = Path(root) / "A.fa", Path(root) / "B.fa"
    fasta_a.write_text(f">chr1\n{seq_a}\n", encoding="utf-8")
    fasta_b.write_text(f">chr1\n{seq_b}\n", encoding="utf-8")
    occurrences = [
        {"occurrence_id": "a_left", "family_id": "fam", "species": "A", "gene_copy_id": "geneA",
         "transcript_id": "txA", "role": "exon", "presence_status": "present", "contig": "chr1",
         "start": 1, "end": 4, "strand": "+", "source_feature_id": "f_a_left"},
        {"occurrence_id": "a_right", "family_id": "fam", "species": "A", "gene_copy_id": "geneA",
         "transcript_id": "txA", "role": "exon", "presence_status": "present", "contig": "chr1",
         "start": 13, "end": 16, "strand": "+", "source_feature_id": "f_a_right"},
        {"occurrence_id": "b_left", "family_id": "fam", "species": "B", "gene_copy_id": "geneB",
         "transcript_id": "txB", "role": "exon", "presence_status": "present", "contig": "chr1",
         "start": 1, "end": 4, "strand": "+", "source_feature_id": "f_b_left"},
        {"occurrence_id": "b_right", "family_id": "fam", "species": "B", "gene_copy_id": "geneB",
         "transcript_id": "txB", "role": "exon", "presence_status": "present", "contig": "chr1",
         "start": 10, "end": 13, "strand": "+", "source_feature_id": "f_b_right"},
    ]
    write_tsv(prepared / "segment_occurrences.tsv", occurrences, list(occurrences[0]))
    edges = [
        {"query_occurrence_id": "a_left", "subject_occurrence_id": "b_left",
         "position_edge_eligible": "1", "membership_edge_eligible": "1", "match_status": "mapped",
         "candidate_resolution": "resolved", "query_genomic_matched_blocks": "chr1:1-4:+",
         "subject_genomic_matched_blocks": "chr1:1-4:+"},
        {"query_occurrence_id": "a_right", "subject_occurrence_id": "b_right",
         "position_edge_eligible": "1", "membership_edge_eligible": "1", "match_status": "mapped",
         "candidate_resolution": "resolved", "query_genomic_matched_blocks": "chr1:13-16:+",
         "subject_genomic_matched_blocks": "chr1:10-13:+"},
    ]
    write_tsv(prepared / "segment_matches.tsv", edges, list(edges[0]))
    loci = [
        {"species": "A", "gene_copy_id": "geneA", "contig": "chr1", "strand": "+",
         "search_start": 1, "search_end": len(seq_a), "genome_fasta": str(fasta_a.resolve())},
        {"species": "B", "gene_copy_id": "geneB", "contig": "chr1", "strand": "+",
         "search_start": 1, "search_end": len(seq_b), "genome_fasta": str(fasta_b.resolve())},
    ]
    write_tsv(prepared / "gene_loci.tsv", loci, list(loci[0]))
    targets = [{"family_id": "fam", "species": "A", "gene_copy_id": "geneA"},
               {"family_id": "fam", "species": "B", "gene_copy_id": "geneB"}]
    write_tsv(prepared / "input_targets.tsv", targets, list(targets[0]))
    write_tsv(prepared / "segment_homology.tsv", [], ["homology_id", "occurrence_id"])
    (prepared / "segment_sequences.fasta").write_text(">fixture\nA\n", encoding="utf-8")
    return prepared


def _analysis_args(model, prepared, override, evidence):
    values = dict(command="analyze", model=model, input_dir=str(prepared),
                  species_tree=str(override) if override else None,
                  genomic_evidence_dir=str(evidence) if evidence else None)
    if model == "dna-presence-ctmc":
        values.update(survey_min_identity=.7, survey_min_coverage=.8, survey_max_dp_cells=250000)
    else:
        values.update(intron_anchor_window=15, intron_min_anchor_pairs=8)
    return SimpleNamespace(**values)


if __name__ == "__main__":
    unittest.main()
