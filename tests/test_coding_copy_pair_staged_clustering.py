"""Exact output/order contracts for bounded copy-pair staging."""

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from intraphy.mapping.match_staging import MatchStagingStore, canonical_copy_pair_key
from intraphy.mapping.clustering import cluster_segments
from intraphy.mapping import staged_clustering
from support_coding_copy_pairs import unchanged_protein_alignment, weak_dna_evidence
from support_coding_copy_pair_staging import (
    open_staging_provider,
    staging_fixture,
    terminal_partition_alignment,
    terminal_partition_fixture,
)


class CodingCopyPairStagedClusteringTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)
        self.occurrences, self.paths, self.sequences = staging_fixture()
        self.provider = open_staging_provider(
            self.occurrences, self.paths, self.sequences, self.root,
        )
        self.addCleanup(self.provider.close)

    def _cluster(self, *, staging_path=None, executor="thread", threads=1,
                 collect_matches=True, aligner="mafft", use_provider=True):
        written = []
        kwargs = {
            "transcript_paths": self.paths,
            "match_writer": written.append,
            "threads": threads,
            "pair_scoring_executor": executor,
            "collect_matches": collect_matches,
            "aligner": aligner,
            "min_size_ratio": 0.25,
        }
        if use_provider:
            kwargs["protein_evidence_provider"] = self.provider
        if staging_path is not None:
            kwargs["match_staging_path"] = staging_path
        return cluster_segments(self.occurrences, self.sequences, **kwargs), written

    def test_staged_rows_and_components_equal_candidate_memory_path_in_original_order(self):
        with patch(
            "intraphy.coding_correspondence.alignment.protein_multiple_alignment",
            side_effect=unchanged_protein_alignment,
        ), patch(
            "intraphy.mapping.clustering.match_evidence",
            side_effect=weak_dna_evidence,
        ):
            memory_result, memory_written = self._cluster(use_provider=False)
            staged, staged_written = self._cluster(
                staging_path=self.root / "complete.sqlite",
            )

        self.assertEqual(staged, memory_result)
        self.assertEqual(staged_written, memory_written)
        self.assertTrue(staged_written)
        self.assertEqual(
            [row["match_id"] for row in staged_written],
            [f"match_{index:05d}" for index in range(1, len(staged_written) + 1)],
        )
        self.assertEqual(
            [row for row in staged_written if row["match_status"] == "mapped"],
            staged[1],
        )
        self.assertIn("same_copy_alternative_overlap", {
            row["distance_class"] for row in staged_written
        })
        tied_rows = [
            row for row in staged_written
            if {row["query_occurrence_id"], row["subject_occurrence_id"]}
            == {"a1", "b1"}
        ]
        self.assertTrue(tied_rows)
        self.assertEqual(tied_rows[0]["protein_mapping_status"], "ambiguous_mapping")
        self.assertGreaterEqual(int(tied_rows[0]["protein_candidate_mapping_count"]), 2)

    def test_spawn_width_one_and_four_preserve_rows_homology_and_order(self):
        serial, serial_written = self._cluster(
            staging_path=self.root / "process-one.sqlite",
            executor="process", threads=1, aligner="internal",
        )
        parallel, parallel_written = self._cluster(
            staging_path=self.root / "process-four.sqlite",
            executor="process", threads=4, aligner="internal",
        )
        self.assertEqual(parallel, serial)
        self.assertEqual(parallel_written, serial_written)
        self.assertEqual(
            [row["match_id"] for row in parallel_written],
            [f"match_{index:05d}" for index in range(1, len(parallel_written) + 1)],
        )

    def test_complementary_terminal_partition_matches_memory_and_staged_paths(self):
        occurrences, paths, sequences = terminal_partition_fixture()
        provider = open_staging_provider(
            occurrences,
            paths,
            sequences,
            self.root / "terminal-provider",
            protein_alignment=terminal_partition_alignment,
        )
        self.addCleanup(provider.close)
        memory_rows, staged_rows = [], []
        with patch(
            "intraphy.coding_correspondence.alignment.protein_multiple_alignment",
            side_effect=terminal_partition_alignment,
        ), patch(
            "intraphy.mapping.clustering.match_evidence",
            side_effect=weak_dna_evidence,
        ):
            memory_result = cluster_segments(
                occurrences,
                sequences,
                transcript_paths=paths,
                match_writer=memory_rows.append,
                aligner="mafft",
                threads=1,
            )
            staged = cluster_segments(
                occurrences,
                sequences,
                transcript_paths=paths,
                protein_evidence_provider=provider,
                match_writer=staged_rows.append,
                aligner="mafft",
                threads=1,
                match_staging_path=self.root / "terminal-stage.sqlite",
            )
        self.assertEqual(staged, memory_result)
        self.assertEqual(staged_rows, memory_rows)
        resolved = [
            row for row in staged_rows
            if row.get("protein_mapping_status") == "resolved_joint_terminal_partition"
        ]
        self.assertEqual(len(resolved), 2)
        self.assertTrue(all(row["match_status"] == "mapped" for row in resolved))

    def test_one_group_rerun_reapplies_chains_to_false_return_mutated_group(self):
        path = self.root / "global-second-pass.sqlite"
        store = MatchStagingStore(path)
        self.addCleanup(store.close)
        successful = canonical_copy_pair_key(
            ("fam", "A", "gA"), ("fam", "B", "gB"),
        )
        mutated = canonical_copy_pair_key(
            ("fam", "A", "gA2"), ("fam", "C", "gC"),
        )
        store.add(1, successful, {
            "group": "successful-rerun", "match_id": 1,
            "match_status": "candidate_unanchored",
        })
        store.add(2, mutated, {
            "group": "false-return-mutated", "match_id": 2,
            "match_status": "mapped",
        })
        occurrence_by_id = {}
        occurrences = []
        first_pass = []
        second_pass = []

        def resolve(rows, *_args, **_kwargs):
            label = rows[0]["group"]
            first_pass.append(label)
            if label == "successful-rerun":
                rows[0]["match_status"] = "mapped"
                return rows, True
            rows[0]["match_status"] = "candidate_unanchored"
            return rows, False

        def reapply(rows, *_args, **_kwargs):
            second_pass.append([row["group"] for row in rows])
            for row in rows:
                row["second_chain_pass"] = True
                if row["group"] == "false-return-mutated":
                    row["match_status"] = "mapped"

        with patch.object(
            staged_clustering, "resolve_copy_pair_rows", side_effect=resolve,
        ), patch.object(
            staged_clustering, "reapply_ordered_chains", side_effect=reapply,
        ):
            staged_clustering.resolve_staged_pair_groups(
                store,
                occurrence_by_id,
                occurrences,
                transcript_paths=(),
                seqs={},
                gene_loci={"present": True},
                short_alignment_max_dp_cells=1,
            )

        self.assertEqual(first_pass, ["successful-rerun", "false-return-mutated"])
        self.assertEqual(
            second_pass,
            [["successful-rerun"], ["false-return-mutated"]],
        )
        finalized = [row for _ordinal, row, _facts, _overlap in store.iter_rows()]
        self.assertEqual([row["match_id"] for row in finalized], [1, 2])
        self.assertEqual([row["match_status"] for row in finalized], ["mapped", "mapped"])
        self.assertTrue(all(row["second_chain_pass"] for row in finalized))
