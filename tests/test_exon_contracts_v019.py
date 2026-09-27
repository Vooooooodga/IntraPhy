"""V19 command, persistence, interval and optional-tool contracts."""
from contextlib import redirect_stdout, redirect_stderr
from dataclasses import replace
from io import StringIO
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import numpy as np

from intraphy.cli import main
from intraphy.commands.parser import build_parser
from intraphy.commands.preflight import required_tools, validate_arguments
from intraphy.structure.serialization import read_catalogues, write_catalogues, encode_catalogue
from intraphy.structure.types import Catalogue, ExonSpan, ExonConfiguration, ObservationEvidence, Material
from intraphy.structure.space import enumerate_space
from intraphy.structure.observations import compatibility
from intraphy.structure.native import NativeLocus
from intraphy.structure.consequences import native_cds
from intraphy.aligners.cesar_adapter import gene_mode_input
from intraphy.aligners.recorded import run_recorded


class ExonContractTests(unittest.TestCase):
    def catalogue(self):
        e = ExonSpan(0, 30)
        return Catalogue("family", "unit", 30, (e,), (), observations=(ObservationEvidence("A", (ExonConfiguration((e,)),)),))

    def test_schema_roundtrip_and_legacy_rejection(self):
        c = self.catalogue()
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/"catalogue.jsonl"
            write_catalogues(p, (c,))
            self.assertEqual(read_catalogues(p), (c,))
            p.write_text('{"schema":"legacy","state":1}\n')
            with self.assertRaises(ValueError): read_catalogues(p)

    def test_duplicate_json_keys_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/"c.jsonl"
            p.write_text('{"schema":"a","schema":"b"}\n')
            with self.assertRaisesRegex(ValueError, "Duplicate JSON"): read_catalogues(p)

    def test_invalid_observation_dimensions_rejected(self):
        c = self.catalogue()
        bad = ObservationEvidence("A", (ExonConfiguration((ExonSpan(0, 40),)),))
        with self.assertRaises(ValueError): replace(c, observations=(bad,))
        with self.assertRaises(ValueError): replace(c, material=(Material("x", 1, 3),))
        with self.assertRaises(ValueError): replace(c, discovery="all_is_known")

    def test_unknown_material_does_not_force_inactivation(self):
        e = ExonSpan(0, 30)
        c = Catalogue("f", "u", 30, (e,), (), material=(Material("x", 0, 30),), boundary_candidates=(e,))
        space = enumerate_space(c)
        o = ObservationEvidence("A", (ExonConfiguration((), (0,)),), "partial",
            material_presence=(None,), unknown_intervals=(e,))
        allowed = compatibility(space, o)
        self.assertTrue(any(s.exons for s, value in zip(space.states, allowed) if value))

    def test_new_default_is_not_legacy_layer_model(self):
        args = build_parser().parse_args(["run", "--input-dir", "in", "--output-dir", "out"])
        self.assertEqual(args.model, "exon-locus-ctmc")
        self.assertEqual(args.observation_view, "evidence")

    def test_legacy_only_flags_are_not_silently_ignored(self):
        for flag, value in [("--annotation-view", "canonical"), ("--root-frequency", "stationary"),
                            ("--ascertainment", "variable-only"), ("--analysis-range", "high-coverage")]:
            args = build_parser().parse_args(["run", "--model", "exon-parsimony", "--input-dir", "in", "--output-dir", "out", flag, value])
            with self.assertRaises(ValueError): validate_arguments(args)

    def test_explicit_catalogue_needs_no_aligner(self):
        args = build_parser().parse_args(["run", "--model", "exon-parsimony", "--input-dir", "in", "--output-dir", "out", "--exon-configurations", "c.jsonl"])
        self.assertEqual(required_tools(args), [])

    def test_rates_are_explicit_and_output_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp, redirect_stderr(StringIO()):
            p = Path(tmp)/"r.json"
            self.assertEqual(main(["exon-rate-template", "--output", str(p)]), 0)
            self.assertEqual(main(["exon-rate-template", "--output", str(p)]), 2)
            self.assertIn("not fitted", json.loads(p.read_text())["provenance"])

    def test_new_guide_uses_computed_configuration_history(self):
        with tempfile.TemporaryDirectory() as tmp, redirect_stderr(StringIO()):
            out = Path(tmp)/"guide"
            self.assertEqual(main(["explain", "--output-dir", str(out)]), 0)
            data = json.loads((out/"illustrative_model.json").read_text())
            self.assertEqual(data["model"], "exon_configuration_v2")
            self.assertIn("V19 exon configuration", (out/"index.html").read_text())
            self.assertEqual(main(["explain", "--output-dir", str(out)]), 2)

    def test_external_error_preserves_stderr_and_returncode(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(RuntimeError):
                run_recorded([sys.executable, "-c", "import sys;print('diagnostic',file=sys.stderr);sys.exit(3)"], tmp, "tool")
            record = json.loads((Path(tmp)/"tool.command.json").read_text())
            self.assertEqual(record["returncode"], 3)
            self.assertIn("diagnostic", (Path(tmp)/"tool.stderr").read_text())

    def test_external_timeout_retains_partial_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(RuntimeError):
                run_recorded([sys.executable, "-c", "import time;print('begin',flush=True);time.sleep(4)"], tmp, "slow", 1)
            record = json.loads((Path(tmp)/"slow.command.json").read_text())
            self.assertEqual(record["status"], "timeout")
            self.assertIn("begin", (Path(tmp)/"slow.stdout").read_text())

    def test_cds_internal_phase_is_not_trimmed_and_exceptions_not_erased(self):
        locus = NativeLocus("f", "s", "g", "chr", "+", 1, 12, "ATGAGGGGATAA", (), {"tx": ()}, False,
                           {"tx": ((0, 4, 0), (8, 12, 2))})
        result = native_cds(locus, "tx")
        self.assertEqual(result["cds"], "ATGAATAA")
        self.assertFalse(result["used_as_evolutionary_filter"])
        self.assertEqual(native_cds(replace(locus, translation_exceptions=("transl_except=x",)), "tx")["status"], "translation_exception_not_modelled")

    def test_cesar_profiles_are_explicit_and_split_codon_bases_preserved(self):
        from intraphy.structure.types import ExonInstance
        a = ExonInstance("e1", "f", "A", "g", "chr", 0, 4, "+", ("tx",))
        b = ExonInstance("e2", "f", "A", "g", "chr", 7, 12, "+", ("tx",))
        ref = NativeLocus("f", "A", "g", "chr", "+", 1, 12, "ATGAGTGAATAA", (a, b), {"tx": ("e1", "e2")}, False,
                          {"tx": ((0, 4, 0), (7, 12, 2))})
        with tempfile.TemporaryDirectory() as tmp:
            profiles = Path(tmp)
            for name in ("acc_profile.txt", "do_profile.txt", "firstCodon_profile.txt", "lastCodon_profile.txt"):
                (profiles/name).write_text("profile\n")
            text = gene_mode_input(ref, ref, "tx", profiles)
            self.assertIn("ATGa\n", text)
            self.assertIn("aaTAA\n", text)
            self.assertIn("####\n>query_genomic_locus", text)
            self.assertNotIn("--clade human", text)


if __name__ == "__main__": unittest.main()
