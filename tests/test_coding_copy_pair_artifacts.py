"""Immutable family-MSA and directed copy-pair artifact contracts."""

import sqlite3
import unittest
from unittest.mock import patch

from intraphy.coding.copy_pair import (
    open_copy_pair_evidence,
    write_copy_pair_evidence,
)
from intraphy.coding_correspondence import CodingProjectionIndex
from support_coding_copy_pairs import (
    CodingCopyPairTestSupport,
    INPUT_ID,
    SOURCE_ID,
    copy_key,
    indexed_records,
    unchanged_protein_alignment,
)


class CodingCopyPairArtifactTests(CodingCopyPairTestSupport, unittest.TestCase):
    def _records(self):
        occurrences, paths, sequences = self.split_fixture()
        index = CodingProjectionIndex(occurrences, sequences, paths)
        with patch(
            "intraphy.coding_correspondence.alignment.protein_multiple_alignment",
            side_effect=unchanged_protein_alignment,
        ):
            projection = index._family_alignment("fam")
            records = indexed_records(index, projection, self.root / "artifact-work.sqlite")
        return records

    def test_copy_pair_artifact_replays_typed_evidence_read_only(self):
        records = self._records()
        expected_pairs = tuple(key for key, _ in records)
        path = self.root / "pair.sqlite"
        write_copy_pair_evidence(
            path, family_id="fam", query_copy_key=copy_key("A"),
            target_copy_key=copy_key("B"), source_identity=SOURCE_ID,
            input_identity=INPUT_ID, expected_pairs=expected_pairs, records=iter(records),
        )
        connection = sqlite3.connect(path)
        metadata = dict(connection.execute("SELECT key,value FROM metadata"))
        self.assertEqual(metadata["schema"], "intraphy.copy-pair-evidence/1")
        self.assertEqual(metadata["source_identity"], SOURCE_ID)
        self.assertEqual(metadata["input_identity"], INPUT_ID)
        self.assertEqual(metadata["status"], "complete")
        self.assertEqual(
            {tuple(row) for row in connection.execute(
                "SELECT query_occurrence_id,target_occurrence_id FROM expected",
            )},
            set(expected_pairs),
        )
        connection.close()

        provider = open_copy_pair_evidence(
            [path], expected_source_identity=SOURCE_ID,
            expected_input_identity=INPUT_ID, expected_pairs=expected_pairs,
        )
        self.addCleanup(provider.close)
        expected_by_key = dict(records)
        for key in expected_pairs:
            with self.subTest(pair=key):
                self.assertEqual(provider.get(*key), expected_by_key[key])
        with self.assertRaises(KeyError):
            provider.get("absent-query", "absent-target")
        with self.assertRaises(sqlite3.OperationalError):
            provider._connections[0].execute("INSERT INTO evidence VALUES ('x','y','{}')")
        with self.assertRaises(ValueError):
            open_copy_pair_evidence(
                [path], expected_source_identity="wrong-source",
                expected_input_identity=INPUT_ID, expected_pairs=expected_pairs,
            )
        with self.assertRaises(ValueError):
            open_copy_pair_evidence(
                [path], expected_source_identity=SOURCE_ID,
                expected_input_identity="wrong-input", expected_pairs=expected_pairs,
            )
        with self.assertRaises(ValueError):
            open_copy_pair_evidence(
                [path], expected_source_identity=SOURCE_ID,
                expected_input_identity=INPUT_ID,
                expected_pairs=(*expected_pairs, ("missing", "pair")),
            )
        with self.assertRaises(ValueError):
            open_copy_pair_evidence(
                [path, path], expected_source_identity=SOURCE_ID,
                expected_input_identity={
                    (copy_key("A"), copy_key("B")): INPUT_ID,
                    (copy_key("B"), copy_key("A")): INPUT_ID,
                },
                expected_pairs=expected_pairs,
            )

    def test_copy_pair_writer_rejects_partial_or_duplicate_records(self):
        records = self._records()
        expected_pairs = tuple(key for key, _ in records)
        for name, partial in (("partial.sqlite", records[:-1]),):
            with self.subTest(name=name), self.assertRaises(ValueError):
                write_copy_pair_evidence(
                    self.root / name, family_id="fam",
                    query_copy_key=copy_key("A"), target_copy_key=copy_key("B"),
                    source_identity=SOURCE_ID, input_identity=INPUT_ID,
                    expected_pairs=expected_pairs, records=partial,
                )
        with self.assertRaises(sqlite3.IntegrityError):
            write_copy_pair_evidence(
                self.root / "duplicate.sqlite", family_id="fam",
                query_copy_key=copy_key("A"), target_copy_key=copy_key("B"),
                source_identity=SOURCE_ID, input_identity=INPUT_ID,
                expected_pairs=expected_pairs, records=(*records, records[0]),
            )
