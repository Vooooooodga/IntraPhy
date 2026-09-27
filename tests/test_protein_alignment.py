import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from intraphy.alignment import (
    NT_BLASTN_V1,
    AlignmentBackendError,
    anchored_short_alignment,
    global_alignment_stats,
    local_alignment_stats,
    protein_multiple_alignment,
    protein_pair_alignment,
)


class ProteinAlignmentTests(unittest.TestCase):
    def setUp(self):
        availability = patch("intraphy.alignment.shutil.which", return_value="/mock/mafft")
        self.available = availability.start()
        self.addCleanup(availability.stop)
        self.result = SimpleNamespace(returncode=0, stdout=">query\nMK-WX\n>target\nMKQ-X\n", stderr="")
        process = patch("intraphy.alignment.subprocess.run", return_value=self.result)
        self.run = process.start()
        self.addCleanup(process.stop)

    def test_amino_command_preserves_gaps_unknowns_and_query_target_order(self):
        def mafft(command, **kwargs):
            self.assertEqual(Path(command[-1]).read_text(), ">query\nMKWX\n>target\nMKQX\n")
            return self.result

        self.run.side_effect = mafft
        self.result.stdout = ">target description\nMKQ-\nX\n>query\nmk-\nwx\n"
        aligned = protein_pair_alignment("mkwx", "mkqx", threads=4)
        self.assertEqual(aligned, ("MK-WX", "MKQ-X"))
        self.available.assert_called_once_with("mafft")
        self.run.assert_called_once()
        command = self.run.call_args.args[0]
        self.assertEqual(command[:-1], ["/mock/mafft", "--quiet", "--thread", "4", "--amino", "--auto"])

    def test_terminal_alignment_gaps_are_retained(self):
        self.result.stdout = ">query\n--MKWX\n>target\nMKQX--\n"
        self.assertEqual(protein_pair_alignment("MKWX", "MKQX"), ("--MKWX", "MKQX--"))

    def test_nucleotide_mafft_behavior_is_preserved(self):
        self.result.stdout = ">query\nACGT\n>target\nACGT\n"
        stats = global_alignment_stats("ACGT", "ACGT", backend="mafft", threads=2)
        self.assertEqual(stats.backend, "mafft")
        self.assertEqual(stats.identity, 1.0)
        self.assertEqual(stats.coverage, 1.0)
        self.assertEqual(stats.aligned_blocks, [(1, 4, 1, 4)])
        command = self.run.call_args.args[0]
        self.assertEqual(command[:-1], ["/mock/mafft", "--quiet", "--thread", "2", "--nuc", "--auto"])

    def test_missing_mafft_is_explicit(self):
        self.available.return_value = None
        with self.assertRaisesRegex(AlignmentBackendError, "MAFFT was requested but is not available"):
            protein_pair_alignment("MKWX", "MKQX")
        self.run.assert_not_called()

    def test_failed_mafft_is_explicit(self):
        self.result.returncode = 1
        self.result.stderr = "MAFFT could not align the supplied proteins"
        with self.assertRaisesRegex(AlignmentBackendError, "MAFFT could not align"):
            protein_pair_alignment("MKWX", "MKQX")
        self.run.assert_called_once()

    def test_output_must_preserve_both_input_sequences(self):
        outputs = [
            "",
            ">query\nMKWX\n",
            ">query\nMKWX\n>target\nMKQ\n",
            ">query\nMK-W\n>target\nMKQX\n",
            ">query\nMKAX\n>target\nMKQX\n",
            ">query\nMKWX\n>target\nMKQ-\n",
        ]
        for output in outputs:
            with self.subTest(output=output):
                self.result.stdout = output
                with self.assertRaises(AlignmentBackendError):
                    protein_pair_alignment("MKWX", "MKQX")

    def test_caller_owns_stop_and_input_gap_handling(self):
        for query, target in [("", "MKX"), ("MKX", ""), ("MK*", "MKX"), ("M*K", "MKX"), ("MKX", "MK*"), ("MK-X", "MKX")]:
            with self.subTest(query=query, target=target):
                with self.assertRaisesRegex(AlignmentBackendError, "stops normalized by the caller"):
                    protein_pair_alignment(query, target)
        self.run.assert_not_called()

    def test_linsi_family_alignment_uses_stable_ids_and_deterministic_threads(self):
        def mafft(command, **kwargs):
            self.assertEqual(Path(command[-1]).read_text(), ">alpha\nMKWX\n>beta\nMKQX\n")
            return SimpleNamespace(
                returncode=0,
                stdout=">beta description\nMKQ-X\n>alpha\nMK-WX\n",
                stderr="",
            )

        self.run.side_effect = mafft
        aligned = protein_multiple_alignment({"beta": "MKQX", "alpha": "MKWX"}, mode="linsi", threads=8)
        self.assertEqual(list(aligned), ["alpha", "beta"])
        self.assertEqual(aligned, {"alpha": "MK-WX", "beta": "MKQ-X"})
        command = self.run.call_args.args[0]
        self.assertEqual(
            command[:-1],
            [
                "/mock/mafft",
                "--quiet",
                "--thread",
                "8",
                "--threadit",
                "0",
                "--amino",
                "--localpair",
                "--maxiterate",
                "1000",
            ],
        )

    def test_einsi_family_alignment_uses_genafpair(self):
        self.result.stdout = ">alpha\nMK-WX\n>beta\nMKQ-X\n"
        aligned = protein_multiple_alignment(
            [("beta", "MKQX"), ("alpha", "MKWX")],
            mode="einsi",
            threads=3,
        )
        self.assertEqual(list(aligned), ["alpha", "beta"])
        command = self.run.call_args.args[0]
        self.assertEqual(
            command[:-1],
            [
                "/mock/mafft",
                "--quiet",
                "--thread",
                "3",
                "--threadit",
                "0",
                "--amino",
                "--genafpair",
                "--ep",
                "0",
                "--maxiterate",
                "1000",
            ],
        )

    def test_family_alignment_rejects_duplicate_ids_and_residue_changes(self):
        with self.assertRaisesRegex(AlignmentBackendError, "duplicate protein MSA identifier"):
            protein_multiple_alignment([("same", "MKW"), ("same", "MKQ")])
        self.run.assert_not_called()

        self.result.stdout = ">alpha\nMKAX\n>beta\nMKQX\n"
        with self.assertRaisesRegex(AlignmentBackendError, "did not preserve input residues for alpha"):
            protein_multiple_alignment({"alpha": "MKWX", "beta": "MKQX"})


class ShortNucleotideAlignmentTests(unittest.TestCase):
    def test_named_scoring_and_unknown_bases_are_explicit(self):
        mismatch = global_alignment_stats("A", "G", backend="internal")
        self.assertEqual(mismatch.score_scheme, NT_BLASTN_V1)
        self.assertEqual(mismatch.score, -3.0)

        extended_gap = global_alignment_stats("AAAAA", "AA", backend="internal")
        self.assertEqual(extended_gap.score, -7.0)

        unknown = global_alignment_stats("N", "N", backend="internal")
        self.assertEqual(unknown.score_scheme, NT_BLASTN_V1)
        self.assertEqual(unknown.score, 0.0)
        self.assertEqual(unknown.matches, 0)
        self.assertEqual(unknown.unknown_bases, 2)
        self.assertEqual(unknown.known_aligned_pairs, 0)
        self.assertEqual(unknown.unknown_aligned_pairs, 1)
        self.assertEqual(unknown.query_covered_bases, 0)
        self.assertEqual(unknown.target_covered_bases, 0)

    def test_local_api_preserves_distinct_equal_optimal_target_blocks(self):
        result = anchored_short_alignment(
            "AAA",
            "AAATAAA",
            mode="local",
            query_occurrence_id="query_exon",
            left_anchor_id="left_anchor",
            right_anchor_id="right_anchor",
            search_interval={"contig": "chr1", "start": 101, "end": 107, "strand": "+"},
        )
        self.assertTrue(result.enumeration_complete)
        self.assertEqual(result.incomplete_reason, "")
        self.assertEqual(result.score_scheme, NT_BLASTN_V1)
        self.assertEqual(
            {(candidate.target_interval.start0, candidate.target_interval.end0) for candidate in result.candidates},
            {(0, 3), (4, 7)},
        )
        self.assertTrue(all(candidate.query_interval.start0 == 0 for candidate in result.candidates))
        self.assertTrue(all(candidate.query_interval.end0 == 3 for candidate in result.candidates))
        self.assertEqual(result.query_length, 3)
        self.assertEqual(result.target_length, 7)
        self.assertEqual(result.left_anchor_id, "left_anchor")
        self.assertEqual(
            result.search_interval,
            {
                "coordinate_system": "0-based-half-open",
                "contig": "chr1",
                "strand": "+",
                "start0": 100,
                "end0": 107,
            },
        )
        self.assertTrue(result.backend_version)
        candidate_ids = {candidate.candidate_id for candidate in result.candidates}
        self.assertEqual(len(candidate_ids), 2)
        self.assertEqual(set(result.candidate_ids), candidate_ids)
        self.assertEqual(result.hit_count, 2)
        for candidate in result.candidates:
            self.assertEqual(candidate.sequence_kind, "nucleotide")
            self.assertEqual(candidate.nt_identity, 1.0)
            self.assertIsNone(candidate.aa_identity)
            self.assertEqual(candidate.known_aligned_pairs, 3)
            self.assertEqual(candidate.unknown_aligned_pairs, 0)
            self.assertEqual(candidate.query_covered_bases, 3)
            self.assertEqual(candidate.target_covered_bases, 3)
            self.assertEqual(candidate.hit_count, 2)
            self.assertEqual(set(candidate.alternative_candidate_ids), candidate_ids - {candidate.candidate_id})
            self.assertIsNone(candidate.mapping_quality)

        legacy = local_alignment_stats("AAA", "AAATAAA", backend="internal")
        self.assertEqual(legacy.target_start, 1)
        self.assertEqual(legacy.target_end, 3)
        self.assertEqual(legacy.hit_count, 2)
        self.assertEqual(legacy.ambiguous_hit_count, 1)
        self.assertEqual(legacy.alternative_hits[0]["target_start"], 5)
        self.assertEqual(legacy.alternative_hits[0]["target_end"], 7)
        self.assertEqual(legacy.alternative_hits[0]["mapping_quality"], "NA")
        self.assertEqual(legacy.alternative_hits[0]["sequence_kind"], "nucleotide")
        self.assertEqual(legacy.alternative_hits[0]["known_aligned_pairs"], 3)
        self.assertEqual(legacy.alternative_hits[0]["query_length"], 3)
        self.assertEqual(legacy.alternative_hits[0]["target_length"], 7)

    def test_short_query_can_use_longer_anchor_bounded_target(self):
        target = "T" * 20 + "ACGT" + "T" * 20
        result = anchored_short_alignment(
            "ACGT",
            target,
            mode="local",
            left_anchor_id="left",
            right_anchor_id="right",
            search_interval={"contig": "chrB", "start": 50, "end": 93, "strand": "+"},
        )
        self.assertEqual(result.query_length, 4)
        self.assertEqual(result.target_length, 44)
        self.assertGreater(result.target_length, result.query_length)
        self.assertEqual(len(result.candidates), 1)
        candidate = result.primary
        self.assertEqual((candidate.target_interval.start0, candidate.target_interval.end0), (20, 24))
        self.assertEqual(candidate.query_covered_bases, 4)
        self.assertEqual(candidate.target_covered_bases, 4)
        self.assertEqual(candidate.search_interval["start0"], 49)
        self.assertEqual(candidate.search_interval["end0"], 93)

    def test_optimal_alignment_limit_reports_incomplete_enumeration(self):
        result = anchored_short_alignment("AAA", "AAATAAA", mode="local", max_alignments=1)
        self.assertEqual(len(result.candidates), 1)
        self.assertFalse(result.enumeration_complete)
        self.assertEqual(result.incomplete_reason, "optimal_alignment_limit_reached")

    def test_global_api_preserves_alternative_gap_cutpoints(self):
        result = anchored_short_alignment("AAAA", "AAA", mode="global")
        self.assertGreater(len(result.candidates), 1)
        signatures = {
            tuple(
                (block.query.start0, block.query.end0, block.target.start0, block.target.end0)
                for block in candidate.aligned_blocks
            )
            for candidate in result.candidates
        }
        self.assertGreater(len(signatures), 1)
        self.assertEqual(
            {
                (gap.gap_in, gap.query_cut0, gap.target_cut0)
                for candidate in result.candidates
                for gap in candidate.gap_blocks
            },
            {("target", 0, 0), ("target", 1, 1), ("target", 2, 2), ("target", 3, 3)},
        )

    def test_bounded_api_keeps_internal_dp_limit(self):
        with self.assertRaises(AlignmentBackendError) as caught:
            anchored_short_alignment("A" * 501, "A" * 501, mode="global")
        message = str(caught.exception)
        self.assertIn("cells=251001", message)
        self.assertIn("budget=250000", message)
        self.assertIn("complete alignment was rejected without truncating DP states", message)


if __name__ == "__main__":
    unittest.main()
