"""Contracts for endpoint-based genomic exon branch summaries."""
import json
import tempfile
import unittest
from pathlib import Path

from intraphy.inference.genomic_exon_branches import (
    branch_change_rows, classify_exon_change,
)
from intraphy.inference.genomic_exon_output import write_outputs
from intraphy.inference.genomic_exon_run import infer_genomic_exons
from intraphy.storage.tabular import write_tsv
from intraphy.structure.edits import EDIT_KINDS
from intraphy.structure.serialization import write_catalogues, write_json
from intraphy.structure.types import Catalogue, ExonConfiguration, ExonSpan, ObservationEvidence
from intraphy.topology import SpeciesTree


def _tree():
    return SpeciesTree([
        {"node_id": "root", "parent_id": "", "label": "root"},
        {"node_id": "branch", "parent_id": "root", "label": "branch"},
        {"node_id": "A", "parent_id": "branch", "label": "Species A"},
        {"node_id": "B", "parent_id": "branch", "label": "Species B"},
    ])


def _rows(states, endpoint, *, offset=0):
    results = [{"family_id": "f", "units": [{"unit_id": "u", "states": states,
        "ctmc": {"branches": [{"parent": "root", "child": "branch",
                                "endpoint_probabilities": endpoint}]}}]}]
    details = {("f", "u"): {"catalogue": {"alignment_offset": offset}}}
    return branch_change_rows(results, details, _tree())


class GenomicExonBranchTests(unittest.TestCase):
    def test_joint_configuration_mode_is_not_composed_from_marginal_modes(self):
        states = [{"exons": [[0, 2]], "material": []},
                  {"exons": [[4, 6]], "material": []}]
        rows = _rows(states, [[.2, .4], [.4, 0.]])
        self.assertEqual({(r["parent_exons"], r["child_exons"]) for r in rows},
                         {(((0, 2),), ((4, 6),)), (((4, 6),), ((0, 2),))})
        self.assertTrue(all(r["joint_configuration_probability"] == .4 for r in rows))

    def test_unintroduced_and_deleted_material_merge_as_absence(self):
        states = [{"exons": [[0, 2]], "material": [0]},
                  {"exons": [[0, 2]], "material": [2]},
                  {"exons": [[0, 2]], "material": [1]}]
        rows = _rows(states, [[.2, .3, .1], [0., .2, 0.], [0., .1, .1]])
        same_absence = [r for r in rows if r["parent_dna_presence"] == (0,)
                        and r["child_dna_presence"] == (0,)]
        self.assertEqual(len(same_absence), 1)
        self.assertAlmostEqual(same_absence[0]["joint_configuration_probability"], .7)
        self.assertAlmostEqual(same_absence[0]["probability_dna_presence_change"], .2)

    def test_all_exactly_tied_joint_modes_are_retained(self):
        states = [{"exons": [[0, 2]], "material": []},
                  {"exons": [[4, 6]], "material": []}]
        rows = _rows(states, [[.25, .25], [.25, .25]])
        self.assertEqual(len(rows), 4)

    def test_split_fusion_boundary_shift_and_dna_coupling_classifications(self):
        self.assertEqual(classify_exon_change(((0, 10),), ((0, 4), (6, 10))),
                         ("net_split",))
        self.assertEqual(classify_exon_change(((0, 4), (6, 10)), ((0, 10),)),
                         ("net_fusion",))
        self.assertEqual(classify_exon_change(((0, 4),), ((0, 3),)),
                         ("boundary_change",))
        self.assertEqual(classify_exon_change(((0, 10),), ((0, 4), (6, 10)),
                                              dna_changed=True),
                         ("complex_change", "dna_coupled"))
        self.assertEqual(classify_exon_change(((0, 10),), ((0, 4), (6, 10)),
                                              dna_changed=True,
                                              changed_material_intervals=((20, 22),)),
                         ("net_split",))
        self.assertEqual(classify_exon_change(((0, 10),), ((0, 4), (6, 10)),
                                              dna_changed=True,
                                              changed_material_intervals=((4, 6),)),
                         ("complex_change", "dna_coupled"))

    def test_absent_exons_and_missing_ctmc_are_handled(self):
        self.assertEqual(classify_exon_change((), ((0, 3),)), ("span_gain",))
        self.assertEqual(classify_exon_change(((0, 3),), ()), ("span_loss",))
        self.assertEqual(classify_exon_change(((0, 3),), ((0, 3),), dna_changed=True),
                         ("no_exon_structure_change",))
        self.assertEqual(branch_change_rows([{"family_id": "f", "units": [
            {"unit_id": "u", "ctmc": None, "states": []}]}], {}, _tree()), [])

    def test_identical_structure_still_has_a_row_in_local_alignment_coordinates(self):
        rows = _rows([{"exons": [[2, 5]], "material": [0]}], [[1.]], offset=120)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["parent_exons"], ((2, 5),))
        self.assertEqual(rows[0]["coordinate_system"],
                         "local_alignment_interbase_0based_halfopen")
        self.assertEqual(rows[0]["alignment_offset"], 120)
        self.assertEqual(rows[0]["descendant_species"], ["Species A", "Species B"])
        terminal_result = [{"family_id": "f", "units": [{"unit_id": "u",
            "states": [{"exons": [[2, 5]], "material": [0]}],
            "ctmc": {"branches": [{"parent": "branch", "child": "A",
                                    "endpoint_probabilities": [[1.]]}]}}]}]
        terminal = branch_change_rows(terminal_result,
            {("f", "u"): {"catalogue": {"alignment_offset": 120}}}, _tree())
        self.assertEqual(terminal[0]["descendant_species"], ["Species A"])

    def test_default_branch_table_is_registered_without_expected_edit_counts(self):
        with tempfile.TemporaryDirectory() as temp:
            artifacts = write_outputs(Path(temp), tree=[], details=[], summaries=[],
                ancestral=[], branches=[], branch_exon_changes=[], expected_edits=False,
                fit_record={}, diagnostics={})
            self.assertIn("branch_exon_changes.tsv", artifacts)
            table = Path(temp) / "branch_exon_changes.tsv"
            self.assertTrue(table.exists())
            self.assertIn("dna_presence_scope", table.read_text().splitlines()[0].split("\t"))

    def test_small_fixed_catalogue_default_run_writes_branch_changes_without_event_counts(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            span = ExonSpan(2, 5)
            catalogue = Catalogue("f", "u", 10, (span,), (), observations=(
                ObservationEvidence("A", (ExonConfiguration((span,)),)),
                ObservationEvidence("B", (ExonConfiguration(()),))),
                boundary_candidates=(span,), alignment_offset=120,
                observation_unit="genomic_exon_spans")
            write_catalogues(root / "catalogues.jsonl", [catalogue])
            write_tsv(root / "species_tree.tsv", [
                {"node_id": "root", "parent_id": "", "label": "root"},
                {"node_id": "A", "parent_id": "root", "label": "A", "branch_length": .5},
                {"node_id": "B", "parent_id": "root", "label": "B", "branch_length": .5},
                {"node_id": "C", "parent_id": "root", "label": "C", "branch_length": .5}],
                ["node_id", "parent_id", "label", "branch_length"])
            write_json(root / "rates.json", {"schema": "intraphy.exon-rates/1",
                "rates": {kind: .2 for kind in EDIT_KINDS},
                "provenance": "Synthetic test parameters"})
            output = root / "results"
            summaries = infer_genomic_exons(root, output,
                configurations=root / "catalogues.jsonl", rates=root / "rates.json",
                parameter_mode="fixed")
            self.assertEqual(summaries[0]["unknown_tips"], ["C"])
            self.assertTrue((output / "branch_exon_changes.tsv").exists())
            self.assertFalse((output / "branch_exon_events.tsv").exists())
            self.assertIn("branch_exon_changes.tsv",
                          json.loads((output / "run_result.json").read_text())["artifacts"])


if __name__ == "__main__":
    unittest.main()
