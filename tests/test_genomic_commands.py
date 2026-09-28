"""Source-level CLI checks for the automatic genomic presence route."""
import tempfile
import unittest
import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from intraphy.commands.genomic import dispatch_analyze, dispatch_prepare_evidence
from intraphy.commands.parser import build_parser
from intraphy.commands.preflight import validate_arguments
from intraphy.cli import main
from intraphy.storage.tabular import read_tsv, write_tsv


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
    def test_raw_genomic_analyze_defaults_to_dna_presence(self):
        args = build_parser().parse_args([
            "analyze", "--fasta", "genomes", "--gff", "annotations",
            "--species-tree", "tree.nwk", "--output-dir", "result",
        ])
        self.assertEqual(args.model, "dna-presence-ctmc")
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
            "analyze", "--fasta", "genomes", "--gff", "annotations",
            "--species-tree", "tree.nwk", "--output-dir", "result",
            "--dna-gain-rate", "0.1",
        ])
        with self.assertRaisesRegex(ValueError, "must be supplied together"):
            validate_arguments(args)
        args = parser.parse_args([
            "analyze", "--fasta", "genomes", "--gff", "annotations",
            "--species-tree", "tree.nwk", "--output-dir", "result",
            "--root-frequency", "fixed", "--root-presence", "0.25",
            "--dna-gain-rate", "0.1", "--dna-loss-rate", "0.2",
        ])
        validate_arguments(args)
        for endpoint in ("0", "1"):
            edge_args = parser.parse_args([
                "analyze", "--fasta", "genomes", "--gff", "annotations",
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
            self.assertEqual(main(["analyze", "--input-dir", str(prepared),
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
            self.assertEqual(main(["analyze", "--input-dir", str(prepared),
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
            self.assertEqual(main(["analyze", "--input-dir", str(prepared),
                                   "--genomic-evidence-dir", str(evidence),
                                   "--output-dir", str(result)]), 0)
            fit = json.loads((result / "dna_fit.json").read_text(encoding="utf-8"))
            run = json.loads((result / "run_result.json").read_text(encoding="utf-8"))
            self.assertEqual(fit["families"][0]["status"], "no_observed_contrast")
            self.assertEqual(run["status"], "completed_with_unresolved")


if __name__ == "__main__":
    unittest.main()
