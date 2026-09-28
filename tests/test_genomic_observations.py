from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from intraphy.coordinates import CoordinateBlock, Interval0
from intraphy.aligners.types import AlignmentBackendError, AlignmentCandidateSet, AlignmentGap
from intraphy.aligners.short import anchored_short_alignment
from intraphy.genomic.groups import build_physical_sites
from intraphy.genomic.observations import prepare_dna_observations
from intraphy.genomic.survey import (
    _align_bounded,
    _candidate_state,
    _locus_for,
    _placement_consensus,
    _query_interval,
    _span_sequence,
)
from intraphy.storage.tabular import read_tsv, write_tsv


def _occurrence(occurrence_id, species, start, end, *, gene="g1", transcript="tx1", strand="+"):
    return {
        "occurrence_id": occurrence_id, "family_id": "fam", "species": species,
        "gene_copy_id": gene, "transcript_id": transcript, "role": "exon",
        "presence_status": "present", "contig": "chr1", "start": start,
        "end": end, "strand": strand, "source_feature_id": f"f_{transcript}",
    }


def _edge(left, right):
    return {
        "query_occurrence_id": left, "subject_occurrence_id": right,
        "position_edge_eligible": "1", "membership_edge_eligible": "1",
        "match_status": "mapped", "query_genomic_matched_blocks": "chr1:1-4:+",
        "subject_genomic_matched_blocks": "chr1:1-4:+",
    }


class PhysicalSiteGroupingTests(unittest.TestCase):
    def _build(self, occurrences, edges):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        root = Path(directory.name)
        write_tsv(root / "segment_occurrences.tsv", occurrences, list(occurrences[0]))
        write_tsv(root / "segment_matches.tsv", edges, list(edges[0]) if edges else ["match_id"])
        return build_physical_sites(root)

    def test_transcript_aliases_collapse_to_one_physical_member(self):
        rows = [_occurrence("a1", "A", 10, 20, transcript="tx1"),
                _occurrence("a2", "A", 10, 20, transcript="tx2"),
                _occurrence("b1", "B", 30, 40)]
        sites, members, _ = self._build(rows, [_edge("a1", "b1")])
        self.assertEqual(len(sites), 1)
        self.assertTrue(sites[0]["eligible"])
        member_a = next(row for row in members if row["species"] == "A")
        self.assertEqual(len(member_a["occurrence_ids"]), 2)
        self.assertEqual(len(member_a["aliases"]), 2)

    def test_overlapping_unequal_intervals_are_retained_but_ineligible(self):
        rows = [_occurrence("a1", "A", 10, 20), _occurrence("a2", "A", 18, 25),
                _occurrence("b1", "B", 30, 40)]
        sites, members, _ = self._build(rows, [_edge("a1", "b1")])
        self.assertEqual(len(members), 3)
        self.assertTrue(all(not site["eligible"] for site in sites if site["family_id"] == "fam"))
        self.assertTrue(any("overlapping_unequal_physical_intervals" in site["reason"] for site in sites))

    def test_many_copy_component_is_kept_and_excluded(self):
        rows = [_occurrence("a", "A", 1, 10), _occurrence("b1", "B", 20, 30),
                _occurrence("b2", "B", 40, 50)]
        sites, members, _ = self._build(rows, [_edge("a", "b1"), _edge("a", "b2")])
        self.assertEqual(len(members), 3)
        self.assertEqual(len(sites), 1)
        self.assertFalse(sites[0]["eligible"])
        self.assertIn("multiple_physical_intervals_in_species", sites[0]["reason"])

    def test_singleton_focal_candidate_is_retained_but_cannot_anchor(self):
        sites, _, _ = self._build([_occurrence("a", "A", 10, 20)], [])
        self.assertTrue(sites[0]["eligible"])
        self.assertFalse(sites[0]["anchor_eligible"])
        self.assertIn("singleton_physical_candidate", sites[0]["reason"])

    def test_duplicate_ids_and_invalid_coordinates_raise(self):
        row = _occurrence("dup", "A", 10, 20)
        with self.assertRaisesRegex(ValueError, "duplicate occurrence_id"):
            self._build([row, dict(row)], [])
        malformed = _occurrence("bad", "A", 0, 20)
        with self.assertRaisesRegex(ValueError, "invalid occurrence coordinates"):
            self._build([malformed], [])


class BoundedSurveyEvidenceTests(unittest.TestCase):
    def test_identical_alias_loci_collapse_but_distinct_genomes_are_ambiguous(self):
        node = {"species": "A", "contig": "chr1", "strand": "+",
                "gene_copy_ids": ["g1", "g2"]}
        first = {"contig": "chr1", "strand": "+", "search_start": "1",
                 "search_end": "100", "genome_fasta": "/tmp/assembly.fa"}
        same = dict(first)
        different = dict(first, genome_fasta="/tmp/other.fa")
        self.assertEqual(_locus_for(node, {("A", "g1"): first, ("A", "g2"): same}), first)
        self.assertIsNone(_locus_for(node, {("A", "g1"): first, ("A", "g2"): different}))

    def test_minus_strand_interval_projection_uses_native_coordinates(self):
        node = {"start": 103, "end": 104}
        self.assertEqual(_query_interval(node, (100, 109), "-"), (5, 7))
        left = {"contig": "chr1", "strand": "-", "start": 106, "end": 109}
        right = {"contig": "chr1", "strand": "-", "start": 100, "end": 102}
        locus = {"search_start": "100", "search_end": "109", "genome_fasta": "genome.fa"}
        with patch("intraphy.genomic.survey.read_fasta_interval", return_value="AACCGGTTAA"):
            sequence, bounds, status = _span_sequence(left, right, locus)
        self.assertEqual((bounds, status), ((100, 109), "ok"))
        self.assertEqual(len(sequence), 10)

    def test_n_base_in_bounded_sequence_is_unknown(self):
        left = {"contig": "chr1", "strand": "+", "start": 1, "end": 3}
        right = {"contig": "chr1", "strand": "+", "start": 7, "end": 9}
        locus = {"search_start": "1", "search_end": "9", "genome_fasta": "genome.fa"}
        with patch("intraphy.genomic.survey.read_fasta_interval", return_value="AAANNAAAA"):
            sequence, bounds, status = _span_sequence(left, right, locus)
        self.assertIsNone(sequence)
        self.assertIsNone(bounds)
        self.assertEqual(status, "flank_bounded_interval_contains_non_acgt_bases")

    def test_bounded_query_only_gap_with_two_supported_flanks_calls_absent_at_position(self):
        candidate = SimpleNamespace(
            aligned_blocks=(CoordinateBlock(Interval0(0, 2), Interval0(0, 2)),
                            CoordinateBlock(Interval0(4, 6), Interval0(2, 4))),
            gap_blocks=(AlignmentGap(Interval0(2, 4), Interval0(2, 2)),),
        )
        result = _candidate_state(candidate, "AACCGG", "AAGG", (2, 4), (0, 2), (4, 6),
                                  (0, 2), (2, 4), 0.7, 0.8)
        self.assertEqual(result[0], "0")
        self.assertEqual(result[1], "absent_at_homologous_position")

    def test_flanks_aligned_outside_their_homologous_intervals_do_not_support_presence(self):
        candidate = SimpleNamespace(
            aligned_blocks=(CoordinateBlock(Interval0(0, 2), Interval0(4, 6)),
                            CoordinateBlock(Interval0(2, 4), Interval0(0, 2)),
                            CoordinateBlock(Interval0(4, 6), Interval0(2, 4))),
            gap_blocks=(),
        )
        result = _candidate_state(candidate, "AACCGG", "CCAAGG", (2, 4),
                                  (0, 2), (4, 6), (0, 2), (4, 6), 0.7, 0.8)
        self.assertEqual(result[0], "unknown")
        self.assertEqual(result[1], "both_flanking_homologous_intervals_not_supported")

    def test_conflicting_optimal_placements_are_unknown(self):
        deletion = SimpleNamespace(
            aligned_blocks=(CoordinateBlock(Interval0(0, 2), Interval0(0, 2)),
                            CoordinateBlock(Interval0(4, 6), Interval0(2, 4))),
            gap_blocks=(AlignmentGap(Interval0(2, 4), Interval0(2, 2)),),
        )
        unresolved = SimpleNamespace(aligned_blocks=(), gap_blocks=())
        alignments = AlignmentCandidateSet(candidates=[deletion, unresolved], enumeration_complete=True)
        state, reason, _ = _placement_consensus(alignments, "AACCGG", "AAGG", (2, 4),
                                                (0, 2), (4, 6), (0, 2), (2, 4), 0.7, 0.8)
        self.assertEqual((state, reason), ("unknown", "optimal_alignments_disagree"))

    def test_no_hit_and_dp_budget_exhaustion_remain_unknown(self):
        empty = AlignmentCandidateSet(candidates=[], enumeration_complete=True)
        state, reason, _ = _placement_consensus(empty, "AAAA", "TTTT", (1, 3),
                                                (0, 1), (3, 4), (0, 1), (3, 4), 0.7, 0.8)
        self.assertEqual((state, reason), ("unknown", "no_bounded_alignment"))
        with patch("intraphy.genomic.survey.anchored_short_alignment",
                   side_effect=AlignmentBackendError("rejected by the DP resource guard")):
            alignments, status = _align_bounded("AAAA", "TTTT", 1, "q", "t")
        self.assertIsNone(alignments)
        self.assertEqual(status, "resource_limit")


class GenomicObservationIntegrationTests(unittest.TestCase):
    def test_conflicting_duplicate_gene_locus_rows_raise_input_error(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inputs = root / "inputs"
            inputs.mkdir()
            occurrence = _occurrence("a1", "A", 1, 4, gene="geneA")
            write_tsv(inputs / "segment_occurrences.tsv", [occurrence], list(occurrence))
            write_tsv(inputs / "segment_matches.tsv", [], ["query_occurrence_id", "subject_occurrence_id",
                     "position_edge_eligible", "membership_edge_eligible", "match_status",
                     "query_genomic_matched_blocks", "subject_genomic_matched_blocks"])
            loci = [
                {"species": "A", "gene_copy_id": "geneA", "contig": "chr1", "strand": "+",
                 "search_start": 1, "search_end": 100, "genome_fasta": "/tmp/a.fa"},
                {"species": "A", "gene_copy_id": "geneA", "contig": "chr1", "strand": "+",
                 "search_start": 1, "search_end": 100, "genome_fasta": "/tmp/b.fa"},
            ]
            write_tsv(inputs / "gene_loci.tsv", loci, list(loci[0]))
            with self.assertRaisesRegex(ValueError, "conflicting gene_loci resources/spans"):
                prepare_dna_observations(inputs, root / "out")

    def test_gene_only_target_remains_unknown_without_annotation_occurrences(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inputs = root / "inputs"
            inputs.mkdir()
            fasta_a, fasta_b = inputs / "A.fa", inputs / "B.fa"
            fasta_a.write_text(">chr1\nACGTACGT\n")
            fasta_b.write_text(">chr1\nACGTACGT\n")
            occurrence = _occurrence("a1", "A", 2, 5, gene="geneA")
            write_tsv(inputs / "segment_occurrences.tsv", [occurrence], list(occurrence))
            write_tsv(inputs / "segment_matches.tsv", [], ["query_occurrence_id", "subject_occurrence_id",
                     "position_edge_eligible", "membership_edge_eligible", "match_status",
                     "query_genomic_matched_blocks", "subject_genomic_matched_blocks"])
            targets = [
                {"family_id": "fam", "species": "A", "gene_copy_id": "geneA"},
                {"family_id": "fam", "species": "B", "gene_copy_id": "geneB"},
            ]
            write_tsv(inputs / "input_targets.tsv", targets, list(targets[0]))
            loci = [
                {"species": "A", "gene_copy_id": "geneA", "contig": "chr1", "strand": "+",
                 "search_start": 1, "search_end": 8, "genome_fasta": str(fasta_a)},
                {"species": "B", "gene_copy_id": "geneB", "contig": "chr1", "strand": "+",
                 "search_start": 1, "search_end": 8, "genome_fasta": str(fasta_b)},
            ]
            write_tsv(inputs / "gene_loci.tsv", loci, list(loci[0]))
            observations = prepare_dna_observations(inputs, root / "out")
            self.assertEqual(len(observations), 2)
            b_row = next(row for row in observations if row["species"] == "B")
            self.assertEqual(b_row["state"], "unknown")
            self.assertEqual(b_row["eligible"], 1)
            self.assertIn(str(fasta_b.resolve()), b_row["survey_genome_fastas"])

    def _prepare_case(self, root, target_sequence, target_right_start, threads, output_name):
        inputs = root / "inputs"
        inputs.mkdir(exist_ok=True)
        source_sequence = "ACGT" + "TT" + "AAA" + "GGG" + "TGCA"
        source_fasta = inputs / "A.fa"
        target_fasta = inputs / "B.fa"
        source_fasta.write_text(f">chr1\n{source_sequence}\n")
        target_fasta.write_text(f">chr1\n{target_sequence}\n")
        occurrences = [
            _occurrence("a_left", "A", 1, 4, gene="geneA"),
            _occurrence("a_focal", "A", 7, 9, gene="geneA"),
            _occurrence("a_right", "A", 13, 16, gene="geneA"),
            _occurrence("b_left", "B", 1, 4, gene="geneB"),
            _occurrence("b_right", "B", target_right_start, target_right_start + 3, gene="geneB"),
        ]
        write_tsv(inputs / "segment_occurrences.tsv", occurrences, list(occurrences[0]))
        write_tsv(inputs / "segment_matches.tsv", [_edge("a_left", "b_left"),
                 _edge("a_right", "b_right")], list(_edge("a_left", "b_left")))
        loci = [
            {"species": "A", "gene_copy_id": "geneA", "contig": "chr1", "strand": "+",
             "search_start": 1, "search_end": len(source_sequence), "genome_fasta": str(source_fasta)},
            {"species": "B", "gene_copy_id": "geneB", "contig": "chr1", "strand": "+",
             "search_start": 1, "search_end": len(target_sequence), "genome_fasta": str(target_fasta)},
        ]
        write_tsv(inputs / "gene_loci.tsv", loci, list(loci[0]))
        output = root / output_name
        observations = prepare_dna_observations(inputs, output, threads=threads)
        focal_site = next(row["site_id"] for row in read_tsv(output / "genomic_members.tsv")
                          if row["occurrence_ids"] == "a_focal")
        states = {row["species"]: row for row in observations if row["site_id"] == focal_site}
        return states

    def test_singleton_focal_site_uses_cross_species_double_flanks(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            deleted = "ACGT" + "TT" + "GGG" + "TGCA"
            serial = self._prepare_case(root, deleted, 10, 1, "out_serial")
            parallel = self._prepare_case(root, deleted, 10, 2, "out_parallel")
            self.assertEqual(serial["A"]["state"], "1")
            self.assertEqual(serial["B"]["state"], "0")
            self.assertEqual({sp: row["state"] for sp, row in serial.items()},
                             {sp: row["state"] for sp, row in parallel.items()})

    def test_unannotated_homologous_dna_rescues_presence_and_n_stays_unknown(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sequence_with_unannotated_site = "ACGT" + "TT" + "AAA" + "GGG" + "TGCA"
            present = self._prepare_case(root, sequence_with_unannotated_site, 13, 1, "out_present")
            self.assertEqual(present["B"]["state"], "1")
            sequence_with_n = "ACGT" + "TT" + "NNN" + "GGG" + "TGCA"
            ambiguous = self._prepare_case(root, sequence_with_n, 13, 1, "out_n")
            self.assertEqual(ambiguous["B"]["state"], "unknown")

    def test_real_pairwise_aligner_stops_on_repeated_optimal_placements(self):
        alignments = anchored_short_alignment("AAAAA", "AAAA", mode="global", max_alignments=1)
        self.assertFalse(alignments.enumeration_complete)
        state, reason, _ = _placement_consensus(
            alignments, "AAAAA", "AAAA", (1, 2), (0, 1), (4, 5),
            (0, 1), (3, 4), 0.7, 0.8,
        )
        self.assertEqual((state, reason), ("unknown", "optimal_placement_ambiguous"))

    def test_real_pairwise_aligner_supports_unique_bounded_dna_gap(self):
        query, target = "ACGTGTTTTTTGCACT", "ACGTGGCACT"
        alignments = anchored_short_alignment(query, target, mode="global", max_alignments=1)
        self.assertTrue(alignments.enumeration_complete)
        state, reason, _ = _placement_consensus(
            alignments, query, target, (5, 11), (0, 5), (11, 16),
            (0, 5), (5, 10), 0.7, 0.8,
        )
        self.assertEqual((state, reason), ("0", "absent_at_homologous_position"))


if __name__ == "__main__":
    unittest.main()
