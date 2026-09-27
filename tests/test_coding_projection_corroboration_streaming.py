"""Streaming and semantic contracts for coding-projection corroboration."""

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from intraphy.storage.tabular import iter_tsv, write_tsv
from intraphy.structure.alignment import FamilyAlignment
from intraphy.structure.corroboration import check_coding_projection
from intraphy.structure.native import NativeLocus
from intraphy.structure.types import ExonInstance


class CodingProjectionCorroborationStreamingTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)
        a_exons = (
            ExonInstance("q1", "fam", "A", "gA", "chrA", 0, 5, "+", ("txA",)),
            ExonInstance("q2", "fam", "A", "gA", "chrA", 6, 11, "+", ("txA",)),
        )
        b_exons = (
            ExonInstance("t1", "fam", "B", "gB", "chrB", 0, 5, "+", ("txB",)),
            ExonInstance("t2", "fam", "B", "gB", "chrB", 6, 11, "+", ("txB",)),
        )
        loci = (
            NativeLocus("fam", "A", "gA", "chrA", "+", 1, 12, "ACGTAAACGTAA", a_exons, {"txA": ("q1", "q2")}, False),
            NativeLocus("fam", "B", "gB", "chrB", "+", 1, 12, "ACGTAAATTTAA", b_exons, {"txB": ("t1", "t2")}, False),
        )
        self.alignment = FamilyAlignment(
            loci=loci,
            rows={"A": "ACGTAAACGTAA--", "B": "ACGTAAA--TTTAA"},
            offsets={"A": 0, "B": 0},
            columns={
                "A": tuple(range(12)),
                "B": tuple(range(7)) + tuple(range(9, 14)),
            },
            exons={},
            anomalies={"q1": ("preexisting",)},
            matches=(),
        )

    def _write_matches(self, directory):
        path = Path(directory) / "segment_matches.tsv"
        write_tsv(path, [
            {
                "query_occurrence_id": "q1",
                "subject_occurrence_id": "t1",
                "protein_hard_observation_eligible": "1",
                "protein_projected_blocks": "1-3:1-3",
            },
            {
                "query_occurrence_id": "q2",
                "subject_occurrence_id": "t2",
                "protein_hard_observation_eligible": "1",
                "protein_projected_blocks": "1-3:1-3",
            },
            {
                "query_occurrence_id": "absent",
                "subject_occurrence_id": "t1",
                "protein_hard_observation_eligible": "1",
                "protein_projected_blocks": "1-3:1-3",
            },
            {
                "query_occurrence_id": "q1",
                "subject_occurrence_id": "t1",
                "protein_hard_observation_eligible": "0",
                "protein_projected_blocks": "1-3:1-3",
            },
        ], [
            "query_occurrence_id", "subject_occurrence_id",
            "protein_hard_observation_eligible", "protein_projected_blocks",
        ])
        return path

    def test_one_stream_pass_preserves_records_order_and_anomalies(self):
        path = self._write_matches(self.root)
        calls = []

        def counted_iter(path_arg, **kwargs):
            calls.append((Path(path_arg), kwargs))
            yield from iter_tsv(path_arg, **kwargs)

        with patch(
            "intraphy.structure.corroboration.iter_tsv",
            side_effect=counted_iter,
        ) as reader:
            updated, records = check_coding_projection(self.alignment, self.root)

        self.assertEqual(reader.call_count, 1)
        self.assertEqual(calls, [(path, {"optional": True})])
        self.assertEqual(records, [
            {"query_exon": "q1", "target_exon": "t1", "checked_bases": 3,
             "agreeing_bases": 3, "status": "concordant"},
            {"query_exon": "q2", "target_exon": "t2", "checked_bases": 3,
             "agreeing_bases": 1, "status": "coordinate_conflict"},
        ])
        self.assertEqual(updated.anomalies["q1"], ("preexisting",))
        self.assertEqual(updated.anomalies["q2"], ("protein_genomic_coordinate_conflict",))
        self.assertEqual(updated.anomalies["t2"], ("protein_genomic_coordinate_conflict",))

    def test_missing_optional_match_table_returns_unchanged_alignment(self):
        updated, records = check_coding_projection(self.alignment, self.root)
        self.assertEqual(records, [])
        self.assertEqual(updated, self.alignment)

    def test_malformed_tsv_row_width_still_raises(self):
        path = self.root / "segment_matches.tsv"
        path.write_text(
            "query_occurrence_id\tsubject_occurrence_id\tprotein_hard_observation_eligible\tprotein_projected_blocks\n"
            "q1\tt1\t1\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(ValueError, "wrong number of TSV fields"):
            check_coding_projection(self.alignment, self.root)

    def test_out_of_range_projection_keeps_coordinate_error(self):
        path = self.root / "segment_matches.tsv"
        write_tsv(path, [{
            "query_occurrence_id": "q1",
            "subject_occurrence_id": "t1",
            "protein_hard_observation_eligible": "1",
            "protein_projected_blocks": "30-32:30-32",
        }], [
            "query_occurrence_id", "subject_occurrence_id",
            "protein_hard_observation_eligible", "protein_projected_blocks",
        ])
        with self.assertRaisesRegex(ValueError, "Protein projection exceeds its genomic exon"):
            check_coding_projection(self.alignment, self.root)
