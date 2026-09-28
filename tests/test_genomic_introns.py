from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from intraphy.genomic.introns import prepare_intron_observations
from intraphy.genomic.intron_coordinates import genomic_evidence
from intraphy.genomic.intron_sites import family_intron_rows
from intraphy.storage.tabular import read_tsv, write_tsv
from support_intron_fixtures import competing_transcript_fixture, spliced_fixture


def _base(position, strand="+"):
    return SimpleNamespace(
        genome_interval=SimpleNamespace(start0=position - 1),
        strand=strand,
    )


class IntronEndpointEvidenceTests(unittest.TestCase):
    def _evidence(self, sequence):
        with tempfile.TemporaryDirectory() as directory:
            fasta = Path(directory) / "genome.fa"
            fasta.write_text(f">chr1\n{sequence}\n")
            row = {"start": 3, "end": len(sequence) - 1, "contig": "chr1"}
            locus = {"search_start": "1", "search_end": str(len(sequence)),
                     "genome_fasta": str(fasta)}
            return genomic_evidence(row, _base(2), _base(len(sequence)),
                                    locus, "C", "G")

    def test_interior_n_does_not_invalidate_verified_annotated_boundary(self):
        evidence, reason = self._evidence("ACNNNNNG")
        self.assertEqual(reason, "annotated_intron_interval_and_coding_endpoints_verified")
        self.assertEqual(evidence["evidence_scope"],
                         "coding_endpoints_and_annotated_interval;intron_interior_unassessed")
        self.assertEqual(evidence["intron_length"], 5)

    def test_endpoint_n_remains_unknown(self):
        evidence, reason = self._evidence("ANNNNNNN")
        self.assertIsNone(evidence)
        self.assertEqual(reason, "genomic_boundary_endpoint_non_acgt")

    def test_missing_declared_genome_is_an_operational_error(self):
        row = {"start": 3, "end": 6, "contig": "chr1"}
        locus = {"search_start": "1", "search_end": "8",
                 "genome_fasta": "/path/that/does/not/exist.fa"}
        with self.assertRaisesRegex(FileNotFoundError, "declared genome FASTA"):
            genomic_evidence(row, _base(2), _base(7), locus, "C", "G")

    def test_indexed_out_of_bounds_error_is_unresolved_coordinate_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            fasta = Path(directory) / "genome.fa"
            fasta.write_text(">chr1\nACGTACGT\n")
            row = {"start": 3, "end": 6, "contig": "chr1"}
            locus = {"search_start": "1", "search_end": "8", "genome_fasta": str(fasta)}
            with patch("intraphy.genomic.intron_coordinates.read_fasta_interval",
                       side_effect=SystemExit("invalid FASTA interval for chr1: 2-2")):
                evidence, reason = genomic_evidence(row, _base(2), _base(7), locus, "C", "G")
        self.assertIsNone(evidence)
        self.assertEqual(reason, "genomic_boundary_sequence_unavailable")


class IntronPreparationIntegrationTests(unittest.TestCase):
    def _prepare(self, inputs, output):
        with patch("intraphy.coding_correspondence.alignment.protein_multiple_alignment",
                   side_effect=lambda records, mode="linsi", threads=1: dict(records)):
            return prepare_intron_observations(inputs, output)

    def test_plus_and_minus_introns_retain_all_three_phase_offsets(self):
        for strand in ("+", "-"):
            for offset in (0, 1, 2):
                with self.subTest(strand=strand, right_offset=offset), tempfile.TemporaryDirectory() as directory:
                    root = Path(directory)
                    inputs = spliced_fixture(root, strand=strand, split=30 + offset)
                    rows = self._prepare(inputs, root / "out")
                    states = {row["species"]: row["state"] for row in rows}
                    self.assertEqual(states, {"A": "1", "B": "0"})
                    a_row = next(row for row in rows if row["species"] == "A")
                    self.assertEqual(a_row["coding_right_base_offset"], str(offset))
                    self.assertEqual(a_row["gff_right_phase"], str((3 - ((30 + offset) % 3)) % 3))

    def test_aliases_collapse_and_gene_only_tip_stays_unknown(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inputs = spliced_fixture(root, aliases=2, gene_only=True)
            rows = self._prepare(inputs, root / "out")
            self.assertEqual({row["species"]: row["state"] for row in rows},
                             {"A": "1", "B": "0", "C": "unknown"})
            self.assertEqual(len(read_tsv(root / "out" / "intron_positions.tsv")), 1)
            self.assertEqual(read_tsv(root / "out" / "intron_families.tsv")[0]["n_annotated_introns"], "1")

    def test_unknown_phase_remains_unresolved_and_source_provenance_is_not_copied_to_other_tips(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inputs = spliced_fixture(root, phase_known=False, gene_only=True)
            rows = self._prepare(inputs, root / "out")
            self.assertTrue(rows)
            self.assertTrue(all(row["state"] == "unknown" for row in rows))
            c_row = next(row for row in rows if row["species"] == "C")
            self.assertEqual(c_row["physical_intervals"], "[]")
            self.assertEqual(c_row["annotation_aliases"], "[]")
            self.assertEqual(c_row["intron_length"], "NA")

    def test_missing_required_prepared_table_is_an_input_error(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inputs = root / "inputs"
            inputs.mkdir()
            occurrences = [{"occurrence_id": "a", "family_id": "fam", "species": "A",
                "gene_copy_id": "gA", "contig": "chr1", "start": 1, "end": 3, "strand": "+"}]
            write_tsv(inputs / "segment_occurrences.tsv", occurrences, list(occurrences[0]))
            with self.assertRaisesRegex(FileNotFoundError, "intron_sites.tsv"):
                prepare_intron_observations(inputs, root / "out")

    def test_no_annotated_candidates_create_summary_but_no_phylogenetic_site(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inputs = spliced_fixture(root, annotated_intron=False)
            rows = self._prepare(inputs, root / "out")
            self.assertEqual(rows, [])
            self.assertEqual(read_tsv(root / "out" / "intron_positions.tsv"), [])
            summary = read_tsv(root / "out" / "intron_families.tsv")
            self.assertEqual(summary[0]["status"], "no_annotated_introns")
            self.assertEqual(summary[0]["n_annotated_introns"], "0")

    def test_family_with_no_cds_paths_has_summary_and_real_unresolved_annotation_only(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inputs = spliced_fixture(root, no_paths=True)
            rows = self._prepare(inputs, root / "out")
            self.assertTrue(rows)
            self.assertTrue(all(row["state"] == "unknown" for row in rows))
            summary = read_tsv(root / "out" / "intron_families.tsv")[0]
            self.assertEqual(summary["status"], "no_usable_coding_projection")
            self.assertEqual(summary["n_projected_sites"], "0")

    def test_competing_physical_transcript_boundaries_are_unknown(self):
        key1, key2, key_b = ("fam", "A", "gA", "tx1"), ("fam", "A", "gA", "tx2"), ("fam", "B", "gB", "txB")
        projection = SimpleNamespace(
            aligned_records={"a": "" + "M" * 20, "b": "" + "M" * 20},
            aliases_by_record={"a": (key1, key2), "b": (key_b,)},
            residue_columns={"a": tuple(range(20)), "b": tuple(range(20))},
            record_by_transcript={key1: "a", key2: "a", key_b: "b"},
            mode="linsi",
        )
        projection.aligned_for = lambda key: projection.aligned_records[projection.record_by_transcript[key]]
        index = SimpleNamespace(family_projection=lambda _family: projection)
        p1 = ("fam", "A", "chr1", "+", 30, 34)
        p2 = ("fam", "A", "chr1", "+", 60, 64)
        aliases = {p1: [{"start": 30, "end": 34}], p2: [{"start": 60, "end": 64}]}
        item = lambda physical, tx: {"species": "A", "physical": physical, "kind": "intron",
            "evidence": {"intron_length": 5}, "reason": "verified", "transcript_key": tx,
            "right_offset": 0, "phase": "0", "annotation_aliases": []}
        candidates = {("fam", 29, 30): [item(p1, key1), item(p2, key2)]}
        loci = {
            ("A", "gA"): [{"genome_fasta": "/tmp/a.fa", "contig": "chr1", "strand": "+",
                            "annotation_start": 1, "annotation_end": 100}],
            ("B", "gB"): [{"genome_fasta": "/tmp/b.fa", "contig": "chr1", "strand": "+",
                            "annotation_start": 1, "annotation_end": 100}],
        }
        with patch("intraphy.genomic.intron_sites._collect_candidates", return_value=(candidates, {})), \
             patch("intraphy.genomic.intron_sites.anchor_metrics", return_value=[{"left_pairs": 9, "right_pairs": 9}]):
            _aligned, sites, calls, _errors, _candidates = family_intron_rows(
                "fam", index, [], aliases, {}, loci, {"A", "B"},
                {"fam": {("A", "gA"), ("B", "gB")}}, anchor_window=15, min_anchor_pairs=8)
        self.assertEqual(len(sites), 1)
        self.assertFalse(sites[0]["eligible"])
        self.assertIn("competing_physical_boundaries", sites[0]["reason"])
        self.assertEqual({row["state"] for row in calls}, {"unknown"})

    def test_same_physical_intron_mapped_to_multiple_msa_sites_excludes_both(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inputs = competing_transcript_fixture(root, same_physical=True)
            rows = self._prepare(inputs, root / "out")
            sites = read_tsv(root / "out" / "intron_positions.tsv")
            self.assertEqual(len(sites), 2)
            self.assertTrue(all(not int(row["eligible"]) for row in sites))
            self.assertTrue(all("physical_boundary_maps_to_multiple_msa_positions" in row["reason"]
                                for row in sites))
            self.assertTrue(all(row["state"] == "unknown" for row in rows))

    def test_different_physical_introns_competing_at_same_cut_are_unknown(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inputs = competing_transcript_fixture(root, same_physical=False)
            rows = self._prepare(inputs, root / "out")
            sites = read_tsv(root / "out" / "intron_positions.tsv")
            self.assertEqual(len(sites), 1)
            self.assertFalse(int(sites[0]["eligible"]))
            self.assertIn("competing_physical_boundaries", sites[0]["reason"])
            self.assertTrue(all(row["state"] == "unknown" for row in rows))


if __name__ == "__main__":
    unittest.main()
