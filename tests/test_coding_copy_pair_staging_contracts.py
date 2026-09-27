"""Lifecycle, progress, and bounded-return contracts for match staging."""

import sqlite3
import tempfile
import unittest
from pathlib import Path

from intraphy.mapping.clustering import cluster_segments
from intraphy.mapping.match_staging import MatchStagingStore, canonical_copy_pair_key
from support_coding_copy_pair_staging import open_staging_provider, staging_fixture


class CodingCopyPairStagingContractTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)
        self.occurrences, self.paths, self.sequences = staging_fixture()
        self.provider = open_staging_provider(
            self.occurrences, self.paths, self.sequences, self.root,
        )
        self.addCleanup(self.provider.close)

    def _run(self, staging_path, *, collect_matches=True, match_writer=None,
             progress_callback=None):
        return cluster_segments(
            self.occurrences,
            self.sequences,
            transcript_paths=self.paths,
            protein_evidence_provider=self.provider,
            pair_scoring_executor="thread",
            threads=1,
            aligner="internal",
            collect_matches=collect_matches,
            match_writer=match_writer,
            match_staging_path=staging_path,
            progress_callback=progress_callback,
        )

    def test_no_collection_still_writes_all_ordered_rows_and_returns_homology(self):
        collected_rows, streamed_rows = [], []
        collected = self._run(
            self.root / "collect.sqlite",
            collect_matches=True,
            match_writer=collected_rows.append,
        )
        streamed = self._run(
            self.root / "stream.sqlite",
            collect_matches=False,
            match_writer=streamed_rows.append,
        )
        self.assertEqual(streamed[0], collected[0])
        self.assertIsNone(streamed[1])
        self.assertEqual(streamed_rows, collected_rows)
        self.assertEqual(
            [row["match_id"] for row in streamed_rows],
            [f"match_{index:05d}" for index in range(1, len(streamed_rows) + 1)],
        )

    def test_progress_reports_all_staged_phases_and_final_pair_totals(self):
        events = []
        self._run(
            self.root / "progress.sqlite",
            collect_matches=False,
            match_writer=lambda _row: None,
            progress_callback=events.append,
        )
        phases = (
            "occurrence_pair_scoring",
            "unordered_copy_pair_resolution",
            "global_aggregation",
            "global_graph",
            "ordinal_output",
        )
        for phase in phases:
            with self.subTest(phase=phase):
                phase_events = [event for event in events if event["phase"] == phase]
                self.assertTrue(any(event["event"] == "phase_start" for event in phase_events))
                self.assertTrue(any(event["event"] == "phase_complete" for event in phase_events))
                self.assertTrue(all(event["elapsed_seconds"] >= 0 for event in phase_events))
                counts = [event["processed_occurrence_pairs"] for event in phase_events]
                self.assertEqual(counts, sorted(counts))
                self.assertTrue(all(
                    event["total_occurrence_pairs"] == 30 for event in phase_events
                ))
                phase_start = next(event for event in phase_events if event["event"] == "phase_start")
                phase_complete = next(event for event in phase_events if event["event"] == "phase_complete")
                self.assertEqual(
                    phase_complete["total_phase_items"],
                    phase_start["total_phase_items"],
                )
        completed_pairs = [
            event for event in events
            if event["event"] == "phase_complete"
            and event["phase"] == "unordered_copy_pair_resolution"
        ]
        self.assertEqual(completed_pairs[-1]["processed_copy_groups"], 7)
        self.assertEqual(completed_pairs[-1]["total_copy_groups"], 7)

    def test_writer_failure_retains_incomplete_stage_and_existing_path_is_not_reused(self):
        path = self.root / "interrupted.sqlite"

        def fail_writer(_row):
            raise RuntimeError("intentional writer interruption")

        with self.assertRaisesRegex(RuntimeError, "intentional writer interruption"):
            self._run(path, match_writer=fail_writer)
        self.assertTrue(path.is_file())
        with sqlite3.connect(path) as connection:
            metadata = dict(connection.execute("SELECT key, value FROM metadata"))
        self.assertEqual(metadata["schema"], "intraphy.match-stage/1")
        self.assertEqual(metadata["status"], "incomplete")
        with self.assertRaises((FileExistsError, ValueError)):
            self._run(path)

    def test_stage_iteration_can_replace_all_rows_across_multiple_fetch_batches(self):
        store = MatchStagingStore(self.root / "batched-replacement.sqlite")
        self.addCleanup(store.close)
        pair_key = canonical_copy_pair_key(("fam", "A", "gA"), ("fam", "B", "gB"))
        row_count = 513
        for ordinal in range(1, row_count + 1):
            store.add(ordinal, pair_key, {"ordinal": ordinal, "state": "before"})

        observed = []
        for ordinal, row, facts, overlap_pair in store.iter_rows():
            observed.append(ordinal)
            row["state"] = "after"
            store.replace_row(ordinal, row, facts, overlap_pair)

        self.assertEqual(observed, list(range(1, row_count + 1)))
        replayed = [
            (ordinal, row["state"])
            for ordinal, row, _facts, _overlap_pair in store.iter_rows()
        ]
        self.assertEqual(replayed, [(ordinal, "after") for ordinal in range(1, row_count + 1)])
