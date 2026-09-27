"""Candidate IDs are unique within chains and retain backend provenance."""

import copy
import json
import unittest

from intraphy.coordinates import Interval0
from intraphy.mapping.candidate_records import _scope_candidate_record_ids
from intraphy.mapping.match_records import _match_row
from intraphy.mapping.ordered_paths import ordered_candidate_chain
from intraphy.mapping.chain_candidates import _collect_chain_candidates


class CandidateIDScopeTests(unittest.TestCase):
    def _evidence(self, candidate_ids, *, dna_ids=()):
        records = [
            {
                "candidate_id": candidate_id,
                "alternative_candidate_ids": list(candidate_ids[:index] + candidate_ids[index + 1:]),
                "accepted": 1,
                "score": 10.0 - index,
                "score_scheme": "nucleotide",
                "aligned_blocks": (),
                "query_start0": 1 + index * 4,
                "query_end0": 3 + index * 4,
                "target_start0": 1 + index * 4,
                "target_end0": 3 + index * 4,
                "query_genomic_blocks0": (Interval0(1 + index * 4, 3 + index * 4),),
                "target_genomic_blocks0": (Interval0(1 + index * 4, 3 + index * 4),),
                "left_anchor_id": "external_left_anchor",
                "right_anchor_id": "external_right_anchor",
            }
            for index, candidate_id in enumerate(candidate_ids)
        ]
        dna_records = [
            {
                "candidate_id": candidate_id,
                "accepted": 1,
                "score": 10.0,
                "score_scheme": "nucleotide",
                "aligned_blocks": (),
            }
            for candidate_id in dna_ids
        ]
        return {
            "candidate_records": records,
            "dna_candidate_assessments": dna_records,
            "alignment_score": 0.9,
            "coverage_score": 0.9,
            "sequence_score": 0.9,
            "left_context_score": 0.5,
            "right_context_score": 0.5,
            "boundary_score": 0.5,
            "phase_score": 0.5,
            "order_score": 0.5,
            "strand_score": 0.5,
            "splice_score": 0.5,
            "size_ratio": 1.0,
            "total_score": 0.9,
            "alignment_cigar": "2M",
            "alignment_backend": "fixture",
            "alignment_mode": "global",
            "alignment_meaning": "fixture evidence",
            "alignment_requested_backend": "fixture",
        }

    def _row(self, evidence, index):
        return _match_row(
            {"occurrence_id": "query", "start": "1", "end": "10", "species": "A", "gene_copy_id": "a", "family_id": "fam"},
            {"occurrence_id": "target", "start": "1", "end": "10", "species": "B", "gene_copy_id": "b", "family_id": "fam"},
            evidence, 0.9, 0.7, "same_copy", "mapped", index,
        )

    def test_scopes_repeated_backend_ids_and_keeps_provenance_and_references(self):
        evidence_a = self._evidence(
            ["query->target.candidate_001", "query->target.candidate_002"],
            dna_ids=["query->target.candidate_001"],
        )
        evidence_b = self._evidence(
            ["query->target.candidate_001", "query->target.candidate_002"],
            dna_ids=["query->target.candidate_001"],
        )
        before_a, before_b = copy.deepcopy(evidence_a), copy.deepcopy(evidence_b)
        row_a, row_b = self._row(evidence_a, 1), self._row(evidence_b, 2)

        ids_a = [record["candidate_id"] for record in row_a["_candidate_records"]]
        ids_b = [record["candidate_id"] for record in row_b["_candidate_records"]]
        self.assertTrue(set(ids_a).isdisjoint(ids_b))
        self.assertEqual(evidence_a, before_a)
        self.assertEqual(evidence_b, before_b)
        for original, scoped in zip(
            before_a["candidate_records"], row_a["_candidate_records"],
        ):
            for field, value in original.items():
                if field not in {"candidate_id", "alternative_candidate_ids"}:
                    self.assertEqual(scoped[field], value)
        self.assertEqual(
            row_a["_candidate_records"][0]["source_candidate_id"],
            "query->target.candidate_001",
        )
        self.assertEqual(
            row_a["_candidate_records"][0]["alternative_candidate_ids"],
            [ids_a[1]],
        )
        self.assertEqual(row_a["_candidate_records"][0]["left_anchor_id"], "external_left_anchor")
        self.assertEqual(row_a["_candidate_records"][0]["right_anchor_id"], "external_right_anchor")

        public_candidate = json.loads(row_a["candidate_assessments"])[0]
        public_dna = json.loads(row_a["dna_candidate_assessments"])[0]
        self.assertEqual(public_candidate["candidate_id"], ids_a[0])
        self.assertEqual(public_candidate["source_candidate_id"], "query->target.candidate_001")
        self.assertEqual(public_candidate["alternative_candidate_ids"], [ids_a[1]])
        self.assertEqual(public_dna["source_candidate_id"], "query->target.candidate_001")
        self.assertTrue(public_dna["candidate_id"].startswith("match_00001.dna_candidate."))

        occurrences = {
            "query": {"occurrence_id": "query", "start": "1", "end": "10", "strand": "+", "species": "A", "gene_copy_id": "a", "family_id": "fam"},
            "target": {"occurrence_id": "target", "start": "1", "end": "10", "strand": "+", "species": "B", "gene_copy_id": "b", "family_id": "fam"},
        }
        _intervals, _records, _owners, groups, _by_id = _collect_chain_candidates(
            list(occurrences.values()), [row_a, row_b], occurrences, (),
        )
        candidates = next(iter(groups.values()))
        chain = ordered_candidate_chain(candidates, 0.0)
        self.assertEqual(len(chain.retained_ids), 4)
        self.assertEqual(
            sorted(candidate.score for candidate in candidates), [9.0, 9.0, 10.0, 10.0],
        )
        self.assertEqual(
            sorted((candidate.query.start0, candidate.query.end0) for candidate in candidates),
            [(1, 3), (1, 3), (5, 7), (5, 7)],
        )
        self.assertEqual(
            sorted((candidate.target.start0, candidate.target.end0) for candidate in candidates),
            [(1, 3), (1, 3), (5, 7), (5, 7)],
        )

    def test_alternative_reference_container_types_and_external_refs_are_preserved(self):
        records = [
            {"candidate_id": "local_a", "alternative_candidate_ids": ("local_b", "external_ref")},
            {"candidate_id": "local_b", "alternative_candidate_ids": "local_a;external_ref"},
        ]
        _scope_candidate_record_ids(records, "match_00003", "candidate")
        self.assertEqual(
            records[0]["alternative_candidate_ids"],
            ("match_00003.candidate.local_b", "external_ref"),
        )
        self.assertIsInstance(records[0]["alternative_candidate_ids"], tuple)
        self.assertEqual(
            records[1]["alternative_candidate_ids"],
            "match_00003.candidate.local_a;external_ref",
        )
        self.assertIsInstance(records[1]["alternative_candidate_ids"], str)

    def test_missing_ids_keep_rank_format_and_source_na(self):
        evidence = self._evidence([None])
        row = self._row(evidence, 7)
        record = row["_candidate_records"][0]
        self.assertEqual(record["candidate_id"], "match_00007.candidate_001")
        self.assertEqual(record["source_candidate_id"], "NA")

        dna_evidence = self._evidence([], dna_ids=[None])
        dna_row = self._row(dna_evidence, 8)
        dna_record = json.loads(dna_row["dna_candidate_assessments"])[0]
        self.assertEqual(dna_record["candidate_id"], "match_00008.dna_candidate_001")
        self.assertEqual(dna_record["source_candidate_id"], "NA")

    def test_duplicate_source_ids_within_row_raise_precise_error(self):
        evidence = self._evidence(["same_id", "same_id"])
        with self.assertRaisesRegex(ValueError, "ambiguous duplicate source candidate ID.*match_00001.*same_id"):
            self._row(evidence, 1)

        dna_evidence = self._evidence([], dna_ids=["same_dna_id", "same_dna_id"])
        with self.assertRaisesRegex(
            ValueError,
            "ambiguous duplicate source candidate ID within dna_candidate.*match_00001.*same_dna_id",
        ):
            self._row(dna_evidence, 1)


if __name__ == "__main__":
    unittest.main()
