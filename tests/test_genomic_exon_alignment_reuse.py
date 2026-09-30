import json
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path
from unittest import mock

from intraphy.commands.genomic_exons import dispatch_analyze
from intraphy.commands.genomic_exon_preflight import validate_analysis_options
from intraphy.commands.parser import build_parser
from intraphy.commands.preflight import required_tools
from intraphy.commands.session import reserve_output
from intraphy.structure.alignment import align_family
from intraphy.structure.native import NativeLocus
from intraphy.structure.prepare import prepare_configurations
from intraphy.structure.types import ExonInstance


def _locus(sequence="ACGTACGTACGT"):
    exon = ExonInstance("exonA", "familyA", "beeA", "geneA", "contigA",
                        0, 4, "+", ("txA",), ((0, 4, 0),))
    return NativeLocus("familyA", "beeA", "geneA", "contigA", "+", 1,
        len(sequence), sequence, (exon,), {"txA": ("exonA",)}, False)


def _write_fasta(path, records):
    path.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in records.items()))


def _recorded_evidence(directory, locus):
    directory.mkdir(parents=True)
    old_work = Path("/scratch/previous_mhc_work/alignment_evidence/family_00001")
    raw_name, query_name = "genomic_loci.fa", "exon_queries.fa"
    _write_fasta(directory / raw_name, {"taxon_00000": locus.sequence})
    _write_fasta(directory / query_name, {"exon_000000": locus.sequence[:4]})
    _write_fasta(directory / "genomic_mafft.stdout", {"taxon_00000": locus.sequence})
    _write_fasta(directory / "genomic_alignment.fa", {locus.species: locus.sequence})
    (directory / "exon_minimap2.stdout").write_text("")
    records = (
        ("genomic_mafft", ["/usr/bin/mafft", "--auto", "--inputorder", "--thread", "12",
          "--threadit", "0", str(old_work / raw_name)]),
        ("exon_minimap2", ["/usr/bin/minimap2", "-c", "-x", "asm20", "--secondary=yes",
          "-N", "50", "-t", "16", str(old_work / raw_name), str(old_work / query_name)]),
    )
    for label, command in records:
        (directory / f"{label}.command.json").write_text(json.dumps({
            "command": command, "shell": False, "timeout_seconds": 600,
            "returncode": 0, "status": "completed"}))


class GenomicExonAlignmentReuseTests(unittest.TestCase):
    def test_reuse_skips_external_calls_accepts_relocated_inputs_and_records_provenance(self):
        locus = _locus()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source_root = root / "published" / "alignment_evidence"
            source_family = source_root / "family_00001"
            _recorded_evidence(source_family, locus)
            source_raw = (source_family / "genomic_loci.fa").read_text()
            with mock.patch("intraphy.structure.prepare.load_native", return_value=(locus,)), \
                    mock.patch("intraphy.structure.prepare.native_cds", return_value={}), \
                    mock.patch("intraphy.structure.prepare.check_coding_projection",
                               side_effect=lambda alignment, input_dir: (alignment, [])), \
                    mock.patch("intraphy.structure.prepare.build_catalogues",
                               return_value=([], [], [], [])), \
                    mock.patch("intraphy.structure.alignment.run_recorded",
                               side_effect=AssertionError("external aligner must not run")):
                prepare_configurations(root / "prepared", root / "output",
                    observation_unit="genomic_exon_spans",
                    alignment_evidence_dir=source_root)
            self.assertEqual((source_family / "genomic_loci.fa").read_text(), source_raw)
            policy = json.loads((root / "output" / "exon_evidence_policy.json").read_text())
            reuse = policy["alignment_evidence_reuse"]
            self.assertEqual(reuse["source_directory"], str(source_root.resolve()))
            self.assertFalse(reuse["new_external_commands_executed"])
            commands = reuse["families"][0]["commands"]
            self.assertEqual(Path(commands["mafft"]["command"][-1]).parent,
                             Path("/scratch/previous_mhc_work/alignment_evidence/family_00001"))

    def test_saved_genomic_input_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "family_00001"
            _recorded_evidence(source, _locus())
            with self.assertRaisesRegex(ValueError, "genomic_loci.fa differs"):
                align_family((_locus("ACGTACGTACGA"),), Path(temporary) / "out",
                             evidence_dir=source)

    def test_saved_query_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "family_00001"
            _recorded_evidence(source, _locus())
            _write_fasta(source / "exon_queries.fa", {"exon_000000": "TTTT"})
            with self.assertRaisesRegex(ValueError, "exon_queries.fa differs"):
                align_family((_locus(),), Path(temporary) / "out", evidence_dir=source)

    def test_species_keyed_msa_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "family_00001"
            _recorded_evidence(source, _locus())
            _write_fasta(source / "genomic_alignment.fa", {"another_bee": _locus().sequence})
            with self.assertRaisesRegex(ValueError, "species-keyed MSA differs"):
                align_family((_locus(),), Path(temporary) / "out", evidence_dir=source)

    def test_failed_command_record_is_rejected_without_fallback(self):
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "family_00001"
            _recorded_evidence(source, _locus())
            path = source / "genomic_mafft.command.json"
            record = json.loads(path.read_text())
            record.update(status="failed", returncode=1)
            path.write_text(json.dumps(record))
            with mock.patch("intraphy.structure.alignment.run_recorded",
                            side_effect=AssertionError("must not fall back to alignment")):
                with self.assertRaisesRegex(ValueError, "did not complete successfully"):
                    align_family((_locus(),), Path(temporary) / "out", evidence_dir=source)

    def test_parser_forwarding_preflight_conflict_and_required_tools(self):
        args = build_parser().parse_args(["analyze", "--input-dir", "prepared",
            "--alignment-evidence-dir", "saved/alignment_evidence", "--output-dir", "out"])
        self.assertEqual(args.alignment_evidence_dir, "saved/alignment_evidence")
        self.assertEqual(required_tools(args), [])
        with mock.patch("intraphy.inference.genomic_exon_run.infer_genomic_exons") as infer:
            dispatch_analyze(args)
        self.assertEqual(infer.call_args.kwargs["alignment_evidence_dir"],
                         "saved/alignment_evidence")

        args.exon_configurations = "catalogues.jsonl"
        with self.assertRaisesRegex(ValueError, "mutually exclusive"):
            validate_analysis_options(args)

    def test_alignment_evidence_is_protected_as_input_from_output_ownership(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "result"
            evidence = output / "alignment_evidence"
            evidence.mkdir(parents=True)
            args = Namespace(output_dir=str(output), force=False,
                alignment_evidence_dir=str(evidence))
            with self.assertRaisesRegex(ValueError, "inputs cannot be overwritten"):
                reserve_output(args)

            source_parent = Path(temporary) / "saved_source"
            nested_output = source_parent / "new_result"
            source_parent.mkdir()
            args.output_dir = str(nested_output)
            args.alignment_evidence_dir = str(source_parent)
            with self.assertRaisesRegex(ValueError, "inside --alignment-evidence-dir"):
                reserve_output(args)


if __name__ == "__main__":
    unittest.main()
