"""Raw-input structural counterexamples; truth files never enter preparation."""
from contextlib import redirect_stdout, redirect_stderr
from io import StringIO
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from intraphy.commands.parser import build_parser
from intraphy.commands.exons import dispatch_analyze
from intraphy.inputs.selection import resolve_inputs
from intraphy.inference.configuration_run import infer_configurations
from intraphy.structure.serialization import read_catalogues
from intraphy.verification.exon_cases import write_exon_example


@unittest.skipUnless(shutil.which("mafft") and shutil.which("minimap2"), "MAFFT and minimap2 are both required")
class RawExonTests(unittest.TestCase):
    def run_case(self, scenario, transform=None):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        raw = write_exon_example(root/"raw", scenario)
        (raw/"truth.json").unlink()
        if transform:
            transform(raw)
        args = build_parser().parse_args(["analyze", "--model", "exon-parsimony", "--fasta", str(raw), "--gff", str(raw),
            "--species-tree", str(raw/"species_tree.nwk"), "--output-dir", str(root/"result")])
        args._input_selection = resolve_inputs(args)
        summary = dispatch_analyze(args)
        out = root/"result"
        details = json.loads((out/"exon_history.json").read_text())
        self.assertEqual(details["model"], "exon_configuration_v2")
        self.assertFalse((out/"structural_site_matrix.tsv").exists())
        self.assertFalse((raw/"manifest.tsv").exists())
        return summary, details, out

    def events(self, details, view="evidence"):
        return [e for unit in details["units"] for e in unit.get("views", {}).get(view, {}).get("events", [])]

    def total(self, summary):
        return sum(row.get("minimum_structural_edits") or 0 for row in summary)

    def test_conserved(self):
        rows, details, _ = self.run_case("conserved")
        self.assertEqual(self.total(rows), 0)
        self.assertTrue(all(r["status"] == "conditional_structure_inference" for r in rows))
        self.assertFalse(self.events(details))

    def test_split_is_one_insertion_not_three_exon_events(self):
        rows, details, _ = self.run_case("split_insertion")
        self.assertEqual(self.total(rows), 1)
        required = [e for e in self.events(details) if e["support"] == "required"]
        self.assertEqual(len(required), 1)
        self.assertEqual(required[0]["operation"], "dna_insertion")
        self.assertIn("exon_split", required[0]["consequences"])

    def test_internal_intronization_does_not_require_length_conservation(self):
        rows, details, _ = self.run_case("intronization")
        # This fixture changes GFF only. It cannot verify biological intronization.
        self.assertEqual(self.total(rows), 0)
        self.assertEqual([e["operation"] for e in self.events(details, "annotation")
                          if e["support"] == "required"], ["split"])

    def assert_fusion(self, phase):
        rows, details, out = self.run_case("fusion_phase"+str(phase))
        self.assertEqual(self.total(rows), 1)
        required = [e for e in self.events(details) if e["support"] == "required"]
        self.assertEqual(len(required), 1)
        self.assertEqual(required[0]["operation"], "dna_deletion")
        self.assertIn("exon_fusion", required[0]["consequences"])
        sequences = json.loads((out/"native_cds_consequences.json").read_text())
        self.assertEqual(len({r["cds"] for r in sequences}), 1)
        self.assertEqual(len({r["protein"] for r in sequences}), 1)

    def test_exact_fusion_phase0(self): self.assert_fusion(0)
    def test_exact_fusion_phase1(self): self.assert_fusion(1)
    def test_exact_fusion_phase2(self): self.assert_fusion(2)

    def test_one_exon_deletion(self):
        rows, details, _ = self.run_case("exon_deletion")
        self.assertEqual(self.total(rows), 1)
        required = [e for e in self.events(details) if e["support"] == "required"]
        self.assertEqual([e["operation"] for e in required], ["dna_deletion"])
        self.assertIn("deleted_exons:1", required[0]["consequences"])

    def test_continuous_deletion_affecting_two_exons_is_not_counted_twice(self):
        rows, details, _ = self.run_case("multi_exon_deletion")
        self.assertEqual(self.total(rows), 1)
        required = [e for e in self.events(details) if e["support"] == "required"]
        self.assertEqual(len(required), 1)
        self.assertIn("deleted_exons:2", required[0]["consequences"])

    def test_missing_annotation_is_not_confirmed_exon_death(self):
        rows, details, _ = self.run_case("annotation_dropout")
        self.assertEqual(self.total(rows), 0)
        self.assertFalse(self.events(details))
        self.assertTrue(any(r["annotation_conditional_minima"] != r["evidence_compatible_minima"] for r in rows))

    def test_assembly_gap_is_not_a_deletion(self):
        rows, details, _ = self.run_case("assembly_gap")
        self.assertEqual(self.total(rows), 0)
        self.assertFalse(any(e["operation"] == "dna_deletion" for e in self.events(details)))

    def assert_shift(self, name, kind):
        rows, details, _ = self.run_case(name)
        self.assertEqual(self.total(rows), 0)  # Boundary annotation alternatives remain explicit.
        required = [e for e in self.events(details, "annotation") if e["support"] == "required"]
        self.assertEqual([e["operation"] for e in required], [kind])

    def test_donor_shift_and_annotation_error_are_separate_conditions(self): self.assert_shift("donor_shift", "donor_shift")
    def test_acceptor_shift_and_annotation_error_are_separate_conditions(self): self.assert_shift("acceptor_shift", "acceptor_shift")

    def test_upstream_indel_is_not_a_boundary_shift(self):
        rows, details, _ = self.run_case("neutral_upstream_indel")
        self.assertEqual(self.total(rows), 0)
        self.assertFalse(self.events(details, "annotation"))

    def test_repeat_is_not_complementary_split(self):
        rows, details, _ = self.run_case("duplication")
        self.assertTrue(any(r["status"] == "unresolved" for r in rows))
        self.assertFalse(any(e["operation"] == "split" for e in self.events(details)))

    def test_local_inversion_is_not_a_deletion(self):
        rows, details, _ = self.run_case("inversion")
        self.assertTrue(any(r["status"] == "unresolved" for r in rows))
        self.assertFalse(any(e["operation"] == "dna_deletion" for e in self.events(details)))

    def test_negative_strand(self):
        rows, details, _ = self.run_case("negative_strand")
        self.assertEqual(self.total(rows), 1)
        self.assertTrue(all(e["strand"] == "-" for u in details["units"] for e in u["catalogue"]["exon_instances"]))

    def test_real_coexisting_structures_are_not_a_synthetic_union(self):
        rows, details, _ = self.run_case("coexisting")
        unit = next(u for u in details["units"] if any(o["kind"] == "coexisting" for o in u["catalogue"]["observations"]))
        self.assertEqual(len(unit["views"]["annotation"]["histories"]), 2)
        self.assertTrue(any(r["coexisting_structures"] for r in rows))

    def test_noncoding_exons_are_not_filtered_by_ORF(self):
        rows, _, out = self.run_case("utr")
        self.assertEqual(len(rows), 4)
        self.assertEqual(self.total(rows), 0)
        self.assertTrue(all(r["status"] == "noncoding_or_CDS_unavailable" for r in json.loads((out/"native_cds_consequences.json").read_text())))

    def test_nine_base_microexon_is_not_silently_removed(self):
        rows, _, _ = self.run_case("microexon")
        self.assertEqual(len(rows), 4)
        self.assertEqual(self.total(rows), 0)

    def test_noncanonical_site_is_not_hard_rejected(self):
        rows, _, _ = self.run_case("noncanonical")
        self.assertEqual(len(rows), 4)
        self.assertEqual(self.total(rows), 0)

    def test_gene_without_exon_annotation_is_retained_as_unknown(self):
        def remove_annotation(raw):
            p = raw/"Species_D.gff3"
            p.write_text("\n".join(line for line in p.read_text().splitlines() if line.startswith("#") or "\tgene\t" in line)+"\n")
        rows, details, _ = self.run_case("conserved", remove_annotation)
        self.assertEqual(self.total(rows), 0)
        for unit in details["units"]:
            obs = next(o for o in unit["catalogue"]["observations"] if o["species"] == "Species_D")
            self.assertIn("locus_has_no_exon_annotation", obs["reasons"])

    def test_catalogue_reuse_does_not_call_alignment(self):
        from unittest.mock import patch
        rows, _, out = self.run_case("conserved")
        with patch("intraphy.structure.prepare.align_family", side_effect=AssertionError("Unexpected alignment")):
            again = infer_configurations(out, out.parent/"rerun", configurations=out/"exon_configurations.jsonl")
        self.assertEqual(self.total(rows), self.total(again))


if __name__ == "__main__": unittest.main()
