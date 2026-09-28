"""CLI routing tests for the independent intron-position character domain."""
import tempfile
import unittest
import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from intraphy.commands.genomic import dispatch_analyze_introns, dispatch_prepare_evidence
from intraphy.commands.parser import build_parser
from intraphy.commands.preflight import validate_arguments
from intraphy.commands.intron_preflight import validate_reusable_evidence
from intraphy.genomic.introns import prepare_intron_observations
from intraphy.cli import main
from intraphy.topology import SpeciesTree
from intraphy.storage.tabular import read_tsv, write_tsv
from support_intron_fixtures import spliced_fixture


def _shift_contiguous_cds_start_by_one_codon(prepared):
    """Move a continuous CDS path while retaining its translated protein."""
    occurrences = read_tsv(prepared / "segment_occurrences.tsv")
    for row in occurrences:
        if row["species"] == "B":
            row["end"] = "63"
    write_tsv(prepared / "segment_occurrences.tsv", occurrences, list(occurrences[0]))
    paths = read_tsv(prepared / "transcript_paths.tsv")
    for row in paths:
        if row["species"] == "B":
            row["cds_intervals"], row["cds_length"] = "4-63", "60"
    write_tsv(prepared / "transcript_paths.tsv", paths, list(paths[0]))
    features = read_tsv(prepared / "raw_gene_features.tsv")
    for row in features:
        if row["species"] == "B":
            if row["type"].lower() == "cds":
                row["start"], row["end"] = "4", "63"
            else:
                row["start"], row["end"] = "1", "63"
    write_tsv(prepared / "raw_gene_features.tsv", features, list(features[0]))
    loci = read_tsv(prepared / "gene_loci.tsv")
    b_fasta = None
    for row in loci:
        if row["species"] == "B":
            row["annotation_end"] = row["search_end"] = "63"
            b_fasta = Path(row["genome_fasta"])
    write_tsv(prepared / "gene_loci.tsv", loci, list(loci[0]))
    if b_fasta is not None:
        b_fasta.write_text(b_fasta.read_text(encoding="utf-8").replace("ATG" * 20, "ATG" * 21),
                            encoding="utf-8")
    occurrence_fasta = prepared / "segment_sequences.fasta"
    occurrence_text = occurrence_fasta.read_text(encoding="utf-8")
    occurrence_fasta.write_text(occurrence_text.replace(">B1\n" + "ATG" * 20,
                                                         ">B1\n" + "ATG" * 21),
                                encoding="utf-8")


class IntronCliTests(unittest.TestCase):
    def test_model_and_preparation_domain_parse_as_independent_route(self):
        parser = build_parser()
        args = parser.parse_args([
            "analyze", "--model", "intron-position-ctmc", "--fasta", "genomes",
            "--gff", "annotations", "--species-tree", "tree.nwk", "--output-dir", "out",
        ])
        self.assertEqual(args.root_frequency, "stationary")
        validate_arguments(args)
        prep = parser.parse_args([
            "prepare-genomic-evidence", "--character-type", "intron-position",
            "--input-dir", "prepared", "--output-dir", "evidence",
        ])
        self.assertEqual(prep.intron_anchor_window, 15)
        self.assertEqual(prep.intron_min_anchor_pairs, 8)
        validate_arguments(prep)

    def test_domain_specific_controls_are_checked(self):
        parser = build_parser()
        invalid = parser.parse_args([
            "analyze", "--model", "intron-position-ctmc", "--input-dir", "prepared",
            "--output-dir", "out", "--intron-anchor-window", "7",
            "--intron-min-anchor-pairs", "8",
        ])
        with self.assertRaisesRegex(ValueError, "cannot exceed"):
            validate_arguments(invalid)
        invalid = parser.parse_args([
            "analyze", "--model", "dna-presence-ctmc", "--fasta", "genomes",
            "--gff", "annotations", "--species-tree", "tree.nwk", "--output-dir", "out",
            "--intron-anchor-window", "16",
        ])
        with self.assertRaisesRegex(ValueError, "apply only"):
            validate_arguments(invalid)

    def test_prepare_dispatch_forwards_intron_thresholds(self):
        args = SimpleNamespace(input_dir="prepared", output_dir="evidence", threads=3,
                               character_type="intron-position", intron_anchor_window=19,
                               intron_min_anchor_pairs=9)
        with patch("intraphy.genomic.introns.prepare_intron_observations", return_value=[]) as prepare:
            dispatch_prepare_evidence(args)
        prepare.assert_called_once_with("prepared", "evidence", threads=3,
                                        anchor_window=19, min_anchor_pairs=9)

    def test_raw_route_skips_correspondence_derivation(self):
        class Selection:
            rows = [{"family_id": "fam"}]

            def write(self, directory):
                self.destination = directory

        with tempfile.TemporaryDirectory() as directory:
            args = SimpleNamespace(
                output_dir=directory, input_dir=None, genomic_evidence_dir=None,
                _input_selection=Selection(), species_tree="tree.nwk", _species_tree=object(),
                flank=1000, max_extension=10000, threads=2, intron_anchor_window=15,
                intron_min_anchor_pairs=8, dna_gain_rate=None, dna_loss_rate=None,
                root_frequency="stationary", root_presence=.5, expected_edits=False,
            )
            with patch("intraphy.case.build_case") as build_case, \
                    patch("intraphy.genomic.introns.prepare_intron_observations", return_value=[{"state": "unknown"}]), \
                    patch("intraphy.inference.intron_presence.analyze_intron_positions") as analyze:
                dispatch_analyze_introns(args)
            self.assertFalse(build_case.call_args.kwargs["derive_correspondence"])
            self.assertEqual(build_case.call_args.kwargs["allow_unannotated"], True)
            self.assertEqual(analyze.call_args.args[0], [{"state": "unknown"}])
            self.assertEqual(analyze.call_args.kwargs["branch_length_unit"], "supplied_tree_units")

    def test_prepared_staged_route_reads_intron_observations_only(self):
        with tempfile.TemporaryDirectory() as directory:
            prepared = Path(directory) / "prepared"
            evidence = Path(directory) / "evidence"
            output = Path(directory) / "result"
            args = SimpleNamespace(
                output_dir=output, input_dir=prepared, genomic_evidence_dir=evidence,
                _species_tree=object(), dna_gain_rate=.1, dna_loss_rate=.2,
                root_frequency="fixed", root_presence=0.0, expected_edits=True, threads=2,
            )
            rows = [{"observation_type": "intron_position_presence", "state": "unknown"}]
            with patch("intraphy.commands.genomic.read_tsv", return_value=rows) as read, \
                    patch("intraphy.inference.intron_presence.analyze_intron_positions") as analyze:
                dispatch_analyze_introns(args)
            read.assert_called_once_with(evidence / "intron_observations.tsv")
            self.assertEqual(analyze.call_args.args[0], rows)
            self.assertEqual(analyze.call_args.kwargs["fixed_rates"], {"gain": .1, "loss": .2})
            self.assertEqual(analyze.call_args.kwargs["root_presence"], 0.0)

    def test_prepared_evidence_reuse_matches_fresh_fixed_rate_inference(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            prepared = spliced_fixture(root)
            (prepared / "species_tree.tsv").write_text(
                "node_id\tparent_id\tlabel\tbranch_length\n"
                "root\t\troot\t\nA\troot\tA\t0.2\nB\troot\tB\t0.3\n",
                encoding="utf-8",
            )
            evidence, staged_result, fresh_result = (root / name for name in
                                                      ("evidence", "staged_result", "fresh_result"))
            with patch("intraphy.commands.preflight.shutil.which", return_value="mafft"), \
                    patch("intraphy.coding_correspondence.alignment.protein_multiple_alignment",
                          side_effect=lambda records, mode="linsi", threads=1: dict(records)):
                self.assertEqual(main(["prepare-genomic-evidence", "--character-type", "intron-position",
                                       "--input-dir", str(prepared), "--output-dir", str(evidence)]), 0)
                observations = read_tsv(evidence / "intron_observations.tsv")
                self.assertEqual({row["species"]: row["state"] for row in observations}, {"A": "1", "B": "0"})
                common = ["--model", "intron-position-ctmc", "--input-dir", str(prepared),
                          "--root-frequency", "fixed", "--root-presence", ".5",
                          "--gain-rate", ".1", "--loss-rate", ".2"]
                self.assertEqual(main(["analyze", *common, "--genomic-evidence-dir", str(evidence),
                                       "--output-dir", str(staged_result)]), 0)
                self.assertEqual(main(["analyze", *common, "--output-dir", str(fresh_result)]), 0)

            staged_fit = json.loads((staged_result / "intron_fit.json").read_text(encoding="utf-8"))
            fresh_fit = json.loads((fresh_result / "intron_fit.json").read_text(encoding="utf-8"))
            staged_run = json.loads((staged_result / "run_result.json").read_text(encoding="utf-8"))
            fresh_run = json.loads((fresh_result / "run_result.json").read_text(encoding="utf-8"))
            self.assertEqual(staged_fit["model"], "intron-position-ctmc")
            self.assertEqual(staged_fit["families"][0]["log_likelihood"],
                             fresh_fit["families"][0]["log_likelihood"])
            self.assertEqual(staged_run["tree_file"], "species_tree.tsv")
            self.assertNotIn("intron_observations.tsv", staged_run["artifacts"])
            self.assertIn("intron_observations.tsv", fresh_run["artifacts"])
            self.assertFalse((staged_result / "dna_fit.json").exists())
            self.assertFalse((fresh_result / "dna_history.json").exists())

            summary_path = evidence / "intron_families.tsv"
            summary = read_tsv(summary_path)
            summary[0]["species"] = "A"
            write_tsv(summary_path, summary, list(summary[0]))
            bad_roster_result = root / "bad_roster"
            self.assertEqual(main(["analyze", *common, "--genomic-evidence-dir", str(evidence),
                                   "--output-dir", str(bad_roster_result)]), 2)

            summary[0]["species"] = "A;B"
            write_tsv(summary_path, summary, list(summary[0]))
            _shift_contiguous_cds_start_by_one_codon(prepared)
            stale_cut_result = root / "stale_cut"
            with patch("intraphy.commands.preflight.shutil.which", return_value="mafft"):
                self.assertEqual(main(["analyze", *common, "--genomic-evidence-dir", str(evidence),
                                       "--output-dir", str(stale_cut_result)]), 2)


    def test_staged_protein_alignment_rejects_changed_translation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            prepared = spliced_fixture(root)
            tree = SpeciesTree([
                {"node_id": "root", "parent_id": "", "label": "root"},
                {"node_id": "A", "parent_id": "root", "label": "A", "branch_length": "0.2"},
                {"node_id": "B", "parent_id": "root", "label": "B", "branch_length": "0.3"},
            ])
            evidence = root / "evidence"
            with patch("intraphy.coding_correspondence.alignment.protein_multiple_alignment",
                       side_effect=lambda records, mode="linsi", threads=1: dict(records)):
                prepare_intron_observations(prepared, evidence)

            sequence_path = prepared / "segment_sequences.fasta"
            sequence_text = sequence_path.read_text(encoding="utf-8")
            sequence_path.write_text(
                sequence_text.replace("\n" + "ATG" * 20, "\n" + "TTG" + "ATG" * 19, 1),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "family MSA residues disagree with transcript"):
                validate_reusable_evidence(evidence, prepared, tree,
                                           anchor_window=15, min_anchor_pairs=8)


if __name__ == "__main__":
    unittest.main()
