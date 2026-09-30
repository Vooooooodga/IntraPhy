"""Thread forwarding for genomic exon preparation and alignment tools."""
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from intraphy.inference.genomic_exon_run import infer_genomic_exons
from intraphy.structure.alignment import align_family
from intraphy.structure.native import NativeLocus
from intraphy.structure.prepare import prepare_configurations
from intraphy.structure.types import ExonInstance


def _locus(species):
    exon = ExonInstance(f"{species}_exon", "family", species, f"{species}_locus",
                        "contig", 0, 4, "+", (), ())
    return NativeLocus("family", species, f"{species}_locus", "contig", "+",
                       1, 4, "ACGT", (exon,), {}, False)


class GenomicExonThreadTests(unittest.TestCase):
    def test_align_family_passes_threads_to_mafft_and_minimap2(self):
        for requested, expected in ((1, "1"), (4, "4")):
            with self.subTest(threads=requested), tempfile.TemporaryDirectory() as temp:
                paf = Path(temp) / "hits.paf"
                paf.write_text("")
                commands = []

                def record(argv, *args, **kwargs):
                    commands.append(argv)
                    return Path(temp) / "alignment.fa" if argv[0] == "mafft" else paf

                with patch("intraphy.structure.alignment.run_recorded", side_effect=record), \
                     patch("intraphy.structure.alignment.parse_fasta",
                           return_value={"taxon_00000": "ACGT", "taxon_00001": "ACGT"}):
                    align_family((_locus("A"), _locus("B")), temp, threads=requested)

                mafft, minimap = commands
                self.assertEqual(mafft[mafft.index("--thread") + 1], expected)
                self.assertEqual(mafft[mafft.index("--threadit") + 1], "0")
                self.assertNotIn("--threadtb", mafft)
                self.assertEqual(minimap[minimap.index("-t") + 1], expected)

    def test_align_family_defaults_to_one_and_rejects_invalid_threads(self):
        with tempfile.TemporaryDirectory() as temp:
            with patch("intraphy.structure.alignment.run_recorded") as run, \
                 patch("intraphy.structure.alignment.parse_fasta",
                       return_value={"taxon_00000": "ACGT", "taxon_00001": "ACGT"}):
                paf = Path(temp) / "hits.paf"
                paf.write_text("")
                run.side_effect = lambda argv, *args, **kwargs: (
                    Path(temp) / "alignment.fa" if argv[0] == "mafft" else paf)
                align_family((_locus("A"), _locus("B")), temp)
                self.assertEqual(run.call_args_list[0].args[0][
                    run.call_args_list[0].args[0].index("--thread") + 1], "1")
            for invalid in (0, -1, True, 1.5):
                with self.subTest(threads=invalid), self.assertRaisesRegex(ValueError, "positive integer"):
                    align_family((_locus("A"),), temp, threads=invalid)

    def test_prepare_and_infer_forward_threads_and_keep_default(self):
        sentinel = RuntimeError("stop after forwarding")
        fake_locus = SimpleNamespace(family="family", species="A", paths={})
        with tempfile.TemporaryDirectory() as temp:
            with patch("intraphy.structure.prepare.load_native", return_value=(fake_locus,)), \
                 patch("intraphy.structure.prepare.align_family", side_effect=sentinel) as align:
                with self.assertRaisesRegex(RuntimeError, "stop after forwarding"):
                    prepare_configurations(temp, Path(temp) / "prepared", threads=3)
                self.assertEqual(align.call_args.kwargs["threads"], 3)
            with patch("intraphy.structure.prepare.load_native", return_value=(fake_locus,)), \
                 patch("intraphy.structure.prepare.align_family", side_effect=sentinel) as align:
                with self.assertRaisesRegex(RuntimeError, "stop after forwarding"):
                    prepare_configurations(temp, Path(temp) / "prepared_default")
                self.assertEqual(align.call_args.kwargs["threads"], 1)

            with patch("intraphy.inference.genomic_exon_run.prepare_configurations",
                       side_effect=sentinel) as prepare:
                with self.assertRaisesRegex(RuntimeError, "stop after forwarding"):
                    infer_genomic_exons(temp, Path(temp) / "inferred", threads=5)
                self.assertEqual(prepare.call_args.kwargs["threads"], 5)
            with patch("intraphy.inference.genomic_exon_run.prepare_configurations",
                       side_effect=sentinel) as prepare:
                with self.assertRaisesRegex(RuntimeError, "stop after forwarding"):
                    infer_genomic_exons(temp, Path(temp) / "inferred_default")
                self.assertEqual(prepare.call_args.kwargs["threads"], 1)

    def test_prepare_rejects_invalid_threads_before_input_access(self):
        for invalid in (0, -1, True, 2.5):
            with self.subTest(threads=invalid), tempfile.TemporaryDirectory() as temp:
                with patch("intraphy.structure.prepare.load_native") as load:
                    with self.assertRaisesRegex(ValueError, "positive integer"):
                        prepare_configurations(temp, Path(temp) / "prepared", threads=invalid)
                    load.assert_not_called()
                with patch("intraphy.inference.genomic_exon_run.prepare_configurations") as prepare:
                    with self.assertRaisesRegex(ValueError, "positive integer"):
                        infer_genomic_exons(temp, Path(temp) / "inferred", threads=invalid)
                    prepare.assert_not_called()


if __name__ == "__main__":
    unittest.main()
