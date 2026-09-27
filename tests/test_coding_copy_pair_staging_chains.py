"""Chain, terminal-partition, and global second-pass staging contracts."""

import copy
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from intraphy.coordinates import CoordinateBlock, Interval0
from intraphy.mapping.match_staging import MatchStagingStore, canonical_copy_pair_key
from intraphy.mapping.ordered_chains import _apply_ordered_candidate_chains
from intraphy.mapping.short_alignment import _rerun_anchor_bounded_short_candidates
from intraphy.mapping.staged_clustering import resolve_staged_pair_groups


class CodingCopyPairStagingChainTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)

    def test_feature_flank_chain_rerun_matches_in_memory_all_fields(self):
        occurrences = {
            "short": {"occurrence_id": "short", "family_id": "fam", "species": "A", "gene_copy_id": "gA", "transcript_id": "txA", "contig": "chrA", "start": "101", "end": "112", "strand": "+"},
            "target": {"occurrence_id": "target", "family_id": "fam", "species": "B", "gene_copy_id": "gB", "transcript_id": "txB", "contig": "chrB", "start": "5", "end": "20", "strand": "+"},
            "left_q": {"occurrence_id": "left_q", "family_id": "fam", "species": "A", "gene_copy_id": "gA", "transcript_id": "txA", "contig": "chrA", "start": "1", "end": "4", "strand": "+"},
            "left_t": {"occurrence_id": "left_t", "family_id": "fam", "species": "B", "gene_copy_id": "gB", "transcript_id": "txB", "contig": "chrB", "start": "1", "end": "4", "strand": "+"},
            "right_q": {"occurrence_id": "right_q", "family_id": "fam", "species": "A", "gene_copy_id": "gA", "transcript_id": "txA", "contig": "chrA", "start": "201", "end": "204", "strand": "+"},
            "right_t": {"occurrence_id": "right_t", "family_id": "fam", "species": "B", "gene_copy_id": "gB", "transcript_id": "txB", "contig": "chrB", "start": "21", "end": "24", "strand": "+"},
        }
        left_record = {
            "candidate_id": "left_anchor", "accepted": 1,
            "score": 8.0, "score_scheme": "nucleotide",
            "aligned_blocks": (CoordinateBlock(Interval0(0, 4), Interval0(0, 4)),),
            "query_genomic_blocks0": (Interval0(0, 4),),
            "target_genomic_blocks0": (Interval0(0, 4),),
        }
        right_record = {
            "candidate_id": "right_anchor", "accepted": 1,
            "score": 8.0, "score_scheme": "nucleotide",
            "aligned_blocks": (CoordinateBlock(Interval0(0, 4), Interval0(0, 4)),),
            "query_genomic_blocks0": (Interval0(200, 204),),
            "target_genomic_blocks0": (Interval0(20, 24),),
        }
        focus_record = {
            "candidate_id": "focus_feature_candidate", "accepted": 1,
            "score": 12.0, "score_scheme": "nucleotide",
            "aligned_blocks": (CoordinateBlock(Interval0(0, 12), Interval0(0, 12)),),
            "query_genomic_blocks0": (Interval0(100, 112),),
            "target_genomic_blocks0": (Interval0(4, 16),),
        }
        rows = [
            {"match_id": "left_match", "query_occurrence_id": "left_q", "subject_occurrence_id": "left_t", "match_status": "mapped", "enumeration_complete": 1, "_candidate_records": [left_record]},
            {"match_id": "right_match", "query_occurrence_id": "right_q", "subject_occurrence_id": "right_t", "match_status": "mapped", "enumeration_complete": 1, "_candidate_records": [right_record]},
            {"match_id": "focus_match", "query_occurrence_id": "short", "subject_occurrence_id": "target", "match_status": "mapped", "enumeration_complete": 1, "alignment_input_transposed": 0, "short_context_route": "feature_bounded_candidate", "retained_candidate_ids": "NA", "threshold": "0.7", "_candidate_records": [focus_record]},
        ]
        transcript_paths = [
            {**occurrences[occurrence_id], "path_rank": rank}
            for occurrence_id, rank in (
                ("left_q", 1), ("short", 2), ("right_q", 3),
                ("left_t", 1), ("target", 2), ("right_t", 3),
            )
        ]
        query = "ACGTTGCAAGTC"
        gene_loci = {
            ("B", "gB"): {
                "sequence": "AAAA" + query + "TTTT" + "CCCC",
                "contig": "chrB", "strand": "+", "interval": Interval0(0, 24),
            },
        }
        seqs = {"short": query}

        memory_rows = copy.deepcopy(rows)
        _apply_ordered_candidate_chains(
            memory_rows, occurrences, list(occurrences.values()),
            transcript_paths=transcript_paths,
        )
        self.assertTrue(_rerun_anchor_bounded_short_candidates(
            memory_rows,
            occurrences,
            seqs,
            gene_loci,
            transcript_paths=transcript_paths,
        ))
        rerun_focus = next(
            row for row in memory_rows if row["match_id"] == "focus_match"
        )
        self.assertEqual(rerun_focus["left_anchor_id"], "left_anchor")
        self.assertEqual(rerun_focus["right_anchor_id"], "right_anchor")
        expected_search_interval = {
            "coordinate_system": "1-based-closed",
            "contig": "chrB",
            "start": 5,
            "end": 20,
            "strand": "+",
        }
        self.assertEqual(
            json.loads(rerun_focus["search_interval"]), expected_search_interval,
        )
        self.assertEqual(rerun_focus["search_interval_side"], "target")
        self.assertEqual(rerun_focus["candidate_enumeration_status"], "complete")
        _apply_ordered_candidate_chains(
            memory_rows, occurrences, list(occurrences.values()),
            transcript_paths=transcript_paths,
        )

        store = MatchStagingStore(self.root / "flank-chain.sqlite")
        self.addCleanup(store.close)
        pair_key = canonical_copy_pair_key(("fam", "A", "gA"), ("fam", "B", "gB"))
        for ordinal, row in enumerate(rows, start=1):
            store.add(ordinal, pair_key, copy.deepcopy(row))
        resolve_staged_pair_groups(
            store,
            occurrences,
            list(occurrences.values()),
            transcript_paths=transcript_paths,
            seqs=seqs,
            gene_loci=gene_loci,
            short_alignment_max_dp_cells=55_000_000,
        )
        staged_rows = [row for _ordinal, row, _facts, _overlap in store.iter_rows()]
        self.assertEqual(staged_rows, memory_rows)
        focus = next(row for row in staged_rows if row["match_id"] == "focus_match")
        self.assertEqual(focus["short_context_route"], "anchor_bounded_local")
        self.assertEqual(
            json.loads(focus["search_interval"]), expected_search_interval,
        )
        self.assertEqual(focus["search_interval_side"], "target")
        self.assertEqual(focus["candidate_enumeration_status"], "complete")
        self.assertTrue(focus["candidate_id"].startswith(
            "focus_match.anchor_bounded_candidate_",
        ))

    def test_any_group_rerun_reapplies_chains_to_false_return_mutated_group(self):
        store = MatchStagingStore(self.root / "global-second-pass.sqlite")
        self.addCleanup(store.close)
        successful = canonical_copy_pair_key(("fam", "A", "gA"), ("fam", "B", "gB"))
        mutated = canonical_copy_pair_key(("fam", "A", "gA2"), ("fam", "C", "gC"))
        store.add(1, successful, {"group": "successful-rerun", "match_id": 1, "match_status": "candidate_unanchored"})
        store.add(2, mutated, {"group": "false-return-mutated", "match_id": 2, "match_status": "mapped"})
        first_pass, second_pass = [], []

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

        with patch(
            "intraphy.mapping.staged_clustering.resolve_copy_pair_rows",
            side_effect=resolve,
        ), patch(
            "intraphy.mapping.staged_clustering.reapply_ordered_chains",
            side_effect=reapply,
        ):
            resolve_staged_pair_groups(
                store, {}, [], transcript_paths=(), seqs={}, gene_loci={"present": True},
                short_alignment_max_dp_cells=55_000_000,
            )

        self.assertEqual(first_pass, ["successful-rerun", "false-return-mutated"])
        self.assertEqual(second_pass, [["successful-rerun"], ["false-return-mutated"]])
        finalized = [row for _ordinal, row, _facts, _overlap in store.iter_rows()]
        self.assertEqual([row["match_id"] for row in finalized], [1, 2])
        self.assertEqual([row["match_status"] for row in finalized], ["mapped", "mapped"])
        self.assertTrue(all(row["second_chain_pass"] for row in finalized))
