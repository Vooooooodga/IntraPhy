"""Clustering provider and spawn-process integration contracts."""

import unittest
from unittest.mock import patch

from intraphy.mapping.clustering import cluster_segments
from support_coding_copy_pairs import (
    CodingCopyPairTestSupport,
    coding_rows,
    strong_dna_evidence,
    unchanged_protein_alignment,
    weak_dna_evidence,
)


class CodingCopyPairIntegrationTests(CodingCopyPairTestSupport, unittest.TestCase):
    def test_explicit_provider_preserves_cluster_output_and_dna_evidence(self):
        occurrences, paths, sequences = [], [], {}
        for species, offset in (("A", 0), ("B", 1000)):
            species_paths = []
            for label, amino_acids, start in (("left", 9, 1), ("middle", 2, 101), ("right", 9, 201)):
                name = f"{species}_{label}"
                sequence = "GCT" * amino_acids
                occurrence, path = coding_rows(name, species, start + offset, sequence)
                occurrences.append(occurrence)
                species_paths.append(path)
                sequences[name] = sequence
            paths.extend(species_paths)
            paths.extend({
                **path,
                "transcript_id": "tx_alias",
                "path_id": path["path_id"].replace(":tx:", ":tx_alias:"),
            } for path in species_paths)
        provider = self._provider_for(occurrences, paths, sequences)
        self.addCleanup(provider.close)
        with patch(
            "intraphy.mapping.clustering.match_evidence",
            side_effect=weak_dna_evidence,
        ), patch(
            "intraphy.coding_correspondence.alignment.protein_multiple_alignment",
            side_effect=unchanged_protein_alignment,
        ):
            legacy = cluster_segments(occurrences, sequences, transcript_paths=paths)
            replay = cluster_segments(
                occurrences, sequences, transcript_paths=paths,
                protein_evidence_provider=provider,
            )
        self.assertEqual(replay, legacy)
        self.assertTrue(replay[1])
        focal = next(row for row in replay[1] if {
            row["query_occurrence_id"], row["subject_occurrence_id"],
        } == {"A_middle", "B_middle"})
        self.assertEqual(focal["dna_match_status"], "low_similarity")
        self.assertTrue(focal["protein_hard_observation_eligible"] in {1, "1", True})
        self.assertEqual(focal["correspondence_basis"], "annotated_CDS_protein")

        with patch(
            "intraphy.mapping.clustering.match_evidence",
            side_effect=strong_dna_evidence,
        ), patch(
            "intraphy.coding_correspondence.alignment.protein_multiple_alignment",
            side_effect=unchanged_protein_alignment,
        ):
            legacy_strong = cluster_segments(occurrences, sequences, transcript_paths=paths)
            replay_strong = cluster_segments(
                occurrences, sequences, transcript_paths=paths,
                protein_evidence_provider=provider,
            )
        self.assertEqual(replay_strong, legacy_strong)
        focal_strong = next(row for row in replay_strong[1] if {
            row["query_occurrence_id"], row["subject_occurrence_id"],
        } == {"A_middle", "B_middle"})
        self.assertEqual(focal_strong["dna_match_status"], "mapped")
        self.assertEqual(focal_strong["correspondence_basis"], "DNA_and_annotated_CDS_protein")

    def test_explicit_provider_missing_pair_fails_without_legacy_fallback(self):
        left, left_path = coding_rows("left", "A", 101, "GCT" * 12)
        right, right_path = coding_rows("right", "B", 501, "GCT" * 12)
        occurrences, paths = [left, right], [left_path, right_path]
        sequences = {"left": "GCT" * 12, "right": "GCT" * 12}
        provider = self._provider_for(
            occurrences, paths, sequences, omit=("left", "right"),
        )
        self.addCleanup(provider.close)
        with patch(
            "intraphy.mapping.clustering.match_evidence",
            side_effect=weak_dna_evidence,
        ), patch(
            "intraphy.coding_correspondence.alignment.protein_multiple_alignment",
            side_effect=AssertionError("provider miss must not build a fallback MSA"),
        ):
            with self.assertRaises((KeyError, ValueError)):
                cluster_segments(
                    occurrences, sequences, transcript_paths=paths,
                    protein_evidence_provider=provider,
                )

    def test_spawn_process_pair_scoring_matches_serial_and_two_thread_results(self):
        sequence = "ACGTCAGTGCATGACTGACCTAGGTCATGC"
        altered = sequence[:14] + "A" + sequence[15:]
        left, _ = coding_rows("left", "A", 101, sequence)
        right, _ = coding_rows("right", "B", 501, altered)
        for row in (left, right):
            row["target_interval_bounded"] = True
        occurrences, sequences = [left, right], {"left": sequence, "right": altered}
        serial_rows = []
        threaded_rows = []
        process_rows = []
        serial = cluster_segments(
            occurrences, sequences, aligner="internal", threads=1,
            short_context_max_length=100, pair_scoring_executor="thread",
            match_writer=serial_rows.append,
        )
        threaded = cluster_segments(
            occurrences, sequences, aligner="internal", threads=2,
            short_context_max_length=100, pair_scoring_executor="thread",
            match_writer=threaded_rows.append,
        )
        processed = cluster_segments(
            occurrences, sequences, aligner="internal", threads=2,
            short_context_max_length=100, pair_scoring_executor="process",
            match_writer=process_rows.append,
        )
        self.assertEqual(threaded, serial)
        self.assertEqual(processed, serial)
        self.assertEqual(threaded_rows, serial_rows)
        self.assertEqual(process_rows, serial_rows)
        self.assertEqual(len(serial_rows), 1)
        self.assertEqual(serial_rows[0]["alignment_requested_backend"], "anchored_short_alignment")
        self.assertNotEqual(serial_rows[0]["alignment_cigar"], "NA")
