"""Regression tests for the native TSV reader's large-field interface."""
import csv
from contextlib import closing
import tempfile
from pathlib import Path
import unittest

from intraphy.storage.tabular import iter_tsv, read_tsv, write_tsv


class LargeTsvFields(unittest.TestCase):
    def test_large_quoted_fields_roundtrip_plain_and_gzip(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            # Includes a tab, embedded newline, quote, and multibyte Unicode.
            evidence = ('\u8bc1\u636e\t"line one\nline two" \u03bb ' * 9000)
            self.assertGreater(len(evidence), 131072)
            expected = [{"id": "segment-1", "evidence": evidence}]

            for suffix in (".tsv", ".tsv.gz"):
                path = root / f"segments{suffix}"
                write_tsv(path, expected, ["id", "evidence"])
                streamed = iter_tsv(path, required=["id", "evidence"])
                self.assertEqual(next(streamed), expected[0])
                with self.assertRaises(StopIteration):
                    next(streamed)
                self.assertEqual(read_tsv(path), expected)

    def test_interleaved_readers_keep_large_field_limit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            evidence = "\u957f\u5b57\u6bb5\u03bb" * 50000
            paths = [root / "first.tsv", root / "second.tsv"]
            for index, path in enumerate(paths):
                write_tsv(path, [{"id": f"{index}-1", "evidence": evidence},
                                 {"id": f"{index}-2", "evidence": evidence}],
                          ["id", "evidence"])

            previous_limit = csv.field_size_limit()
            csv.field_size_limit(131072)
            self.addCleanup(csv.field_size_limit, previous_limit)
            with closing(iter_tsv(paths[0])) as first, closing(iter_tsv(paths[1])) as second:
                self.assertEqual(next(first), {"id": "0-1", "evidence": evidence})
                self.assertEqual(next(second), {"id": "1-1", "evidence": evidence})
                first.close()
                self.assertEqual(next(second), {"id": "1-2", "evidence": evidence})


if __name__ == "__main__":
    unittest.main()
