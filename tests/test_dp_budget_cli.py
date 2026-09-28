"""CLI plumbing for the internal short-alignment DP budget."""
import csv
import io
import random
import tempfile
import unittest
from contextlib import redirect_stderr
from pathlib import Path
from unittest.mock import patch

from intraphy.aligners.types import MAX_INTERNAL_DP_CELLS
from intraphy.aligners.short import AlignmentBackendError, anchored_short_alignment
from intraphy.cli import _dispatch, main
from intraphy.commands.parser import build_parser
from intraphy.preparation.derive import derive_tables


class ShortAlignmentBudgetCliTests(unittest.TestCase):
    def test_budgeted_and_unlimited_alignment_are_identical_for_255474_cells(self):
        rng = random.Random(20260928)
        query = "".join(rng.choice("ACGT") for _ in range(249))
        target = "N" * 389 + query + "N" * 388
        self.assertEqual(len(query) * len(target), 255_474)

        with self.assertRaisesRegex(AlignmentBackendError, "cells=255474, budget=250000"):
            anchored_short_alignment(query, target, mode="local")

        budgeted = anchored_short_alignment(
            query, target, mode="local", max_dp_cells=300_000
        )
        unlimited = anchored_short_alignment(
            query, target, mode="local", max_dp_cells=None
        )
        self.assertEqual(budgeted.candidates, unlimited.candidates)
        self.assertTrue(budgeted.candidates)
        self.assertTrue(budgeted.enumeration_complete)
        self.assertTrue(unlimited.enumeration_complete)

    def test_parse_default_positive_and_unlimited(self):
        parser = build_parser()
        base = ["derive-tables", "--input-dir", "in"]
        self.assertEqual(
            parser.parse_args(base).short_alignment_max_dp_cells,
            MAX_INTERNAL_DP_CELLS,
        )
        self.assertEqual(
            parser.parse_args(base + ["--short-alignment-max-dp-cells", "55000000"])
            .short_alignment_max_dp_cells,
            55_000_000,
        )
        self.assertIsNone(
            parser.parse_args(base + ["--short-alignment-max-dp-cells", "unlimited"])
            .short_alignment_max_dp_cells
        )

    def test_dispatch_forwards_default_positive_and_unlimited(self):
        parser = build_parser()
        with patch("intraphy.cli.derive_tables") as derive:
            for value, expected in (
                ([], MAX_INTERNAL_DP_CELLS),
                (["--short-alignment-max-dp-cells", "55000000"], 55_000_000),
                (["--short-alignment-max-dp-cells", "unlimited"], None),
            ):
                args = parser.parse_args(["derive-tables", "--input-dir", "in", *value])
                _dispatch(args)
                self.assertEqual(
                    derive.call_args.kwargs["short_alignment_max_dp_cells"], expected
                )

    def test_invalid_values_fail_before_output_creation(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "out"
            for value in ("0", "-1", "1.5", "not-an-integer"):
                with self.subTest(value=value), redirect_stderr(io.StringIO()):
                    with self.assertRaises(SystemExit):
                        main([
                            "derive-tables", "--input-dir", str(Path(tmp) / "in"),
                            "--output-dir", str(output),
                            "--short-alignment-max-dp-cells", value,
                        ])
                    self.assertFalse(output.exists())

    def test_derive_report_records_each_budget(self):
        backend = [{"aligner": "internal", "available": 1, "notes": "internal"}]
        budgets = (
            (MAX_INTERNAL_DP_CELLS, str(MAX_INTERNAL_DP_CELLS)),
            (55_000_000, "55000000"),
            (None, "unlimited"),
        )
        with tempfile.TemporaryDirectory() as tmp:
            with (
                patch("intraphy.preparation.derive.read_tsv", return_value=[]),
                patch("intraphy.preparation.derive.parse_fasta", return_value={}),
                patch("intraphy.preparation.derive._gene_locus_records", return_value=[]),
                patch("intraphy.preparation.derive._species_tree_distances", return_value={}),
                patch("intraphy.preparation.derive.cluster_segments", return_value=([], [])),
                patch("intraphy.preparation.derive.available_alignment_backends", return_value=backend),
            ):
                for index, (budget, expected) in enumerate(budgets):
                    output = Path(tmp) / str(index)
                    derive_tables(
                        tmp,
                        output,
                        short_alignment_max_dp_cells=budget,
                    )
                    with (output / "alignment_backend_report.tsv").open(newline="") as handle:
                        report = next(csv.DictReader(handle, delimiter="\t"))
                    self.assertEqual(report["short_alignment_dp_cell_budget"], expected)
                    self.assertIn("not truncated", report["notes"])


if __name__ == "__main__":
    unittest.main()
