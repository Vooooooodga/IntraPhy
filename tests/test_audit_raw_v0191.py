"""Raw FASTA/GFF regressions with truth kept outside the analysis directory."""
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from intraphy.verification.exon_audit_cases import write_audit_example, check_audit_result
from intraphy.commands.parser import build_parser
from intraphy.commands.exons import dispatch_analyze
from intraphy.inputs.selection import resolve_inputs
from intraphy.storage.fasta import parse_fasta


@unittest.skipUnless(shutil.which("mafft") and shutil.which("minimap2"), "Real MAFFT and minimap2 are required")
class AuditRawTests(unittest.TestCase):
    def run_case(self,name,extra=()):
        temporary=tempfile.TemporaryDirectory();self.addCleanup(temporary.cleanup);root=Path(temporary.name)
        raw=write_audit_example(root/"raw",name)
        (raw/"truth.json").replace(root/"evaluation_truth.json")
        original={f.name:f.read_bytes() for f in raw.glob("*.gff3")}
        out=root/"result"
        args=build_parser().parse_args(["analyze","--model","exon-parsimony","--fasta",str(raw),"--gff",str(raw),
            "--species-tree",str(raw/"species_tree.nwk"),"--output-dir",str(out),*extra])
        args._input_selection=resolve_inputs(args)
        summary=dispatch_analyze(args);data=json.loads((out/"exon_history.json").read_text())
        self.assertEqual(original,{f.name:f.read_bytes() for f in raw.glob("*.gff3")})
        if not extra: check_audit_result(name,summary,data,out)
        return summary,data,out

    def test_mis_split_is_not_confirmed_by_a_simpler_tree(self): self.run_case("mis_split_unchanged_DNA")
    def test_mis_fusion_has_atomic_whole_exon_alternative(self): self.run_case("mis_fusion_unchanged_DNA")
    def test_single_end_annotation_error(self): self.run_case("one_boundary_unchanged_DNA")
    def test_both_end_annotation_error(self): self.run_case("both_boundaries_unchanged_DNA")
    def test_first_missing_exon_still_in_search(self): self.run_case("first_exon_dropout")
    def test_last_missing_exon_still_in_search(self): self.run_case("last_exon_dropout")
    def test_changed_split_signals_not_blindly_erased(self): self.run_case("split_with_changed_signals")
    def test_changed_boundary_signal_not_blindly_erased(self): self.run_case("shift_with_changed_signal")
    def test_exact_exon_deletion_has_no_none_string_sort_conflict(self): self.run_case("exact_exon_deletion")

    def test_budget_failure_does_not_crop_to_remaining_exons(self):
        summary,_,out=self.run_case("first_exon_dropout",("--max-locus-bases","100"))
        self.assertTrue(all(r["status"]=="unresolved" for r in summary))
        from intraphy.storage.tabular import read_tsv
        windows=read_tsv(out/"alignment_evidence/family_00001/search_windows.tsv")
        self.assertEqual({r["selected_bases"] for r in windows},{"1140"})
        self.assertTrue(all(r["unsearched_available_bases"]=="1140" for r in windows))
        # Selected and actually executed are explicitly different in the budget diagnostic.
        self.assertFalse((out/"alignment_evidence/family_00001/genomic_alignment.fa").exists())
