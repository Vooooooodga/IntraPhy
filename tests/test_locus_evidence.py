import json
import tempfile
import unittest
from pathlib import Path

from intraphy.inference.locus_evidence import apply_material_evidence
from intraphy.inference.locus_io import load_locus_model


def _model():
    return {
        "schema": "intraphy.exon-locus-model/2", "model": "exon-locus-ctmc",
        "branch_length_unit": "declared opportunity units",
        "tree_provenance": "test tree", "provenance": "test model",
        "independence_provenance": "single locus unit",
        "rates": {"provenance": "fixed test rate", "groups": {
            "delete": {"mode": "fixed", "value": 0.2}}},
        "units": [{
            "family": "Mhc", "unit": "relay",
            "catalogue": {
                "provenance": "synthetic locus",
                "state_model": "irreversible",
                "material": [{"id": "copy_A", "start": 10, "end": 20}],
                "copies": [{"id": "A", "material_ids": ["copy_A"]}],
                "opportunities": [{"id": "delete_A", "outcome_id": "del",
                                   "kind": "dna_deletion", "rate_group": "delete",
                                   "material_deletions": ["copy_A"], "interval": [10, 20]}],
            },
            "root": {"provenance": "fixed present root", "entries": [
                {"state": {"material": [1]}, "weight": 1.0}]},
            "observations": {"provenance": "initial calls", "tips": {
                "A": {"material": [1]},
                "B": {"material": [None]},
            }},
        }],
    }


def _rows():
    return [
        {"family": "Mhc", "unit": "relay", "tip": "A", "material_id": "copy_A",
         "call": "absent", "evidence": "assembly interval search v1", "survey": True,
         "sensitivity": 0.8, "specificity": 0.9},
        {"family": "Mhc", "unit": "relay", "tip": "B", "material_id": "copy_A",
         "call": "unknown", "evidence": "no callable locus evidence", "survey": False},
    ]


class LocusEvidenceBridgeTests(unittest.TestCase):
    def test_complete_table_replaces_dna_calls_and_loads_without_changing_model_states(self):
        source = _model()
        source_tip = source["units"][0]["observations"]["tips"]["A"]
        source_tip.update({"surveyed_material": ["copy_A"],
                           "material_sensitivity": {"copy_A": 0.5},
                           "material_specificity": {"copy_A": 0.5},
                           "material_evidence": {"copy_A": "stale"}})
        result = apply_material_evidence(source, _rows())
        self.assertEqual(source_tip["material"], [1])
        self.assertEqual(source_tip["material_evidence"], {"copy_A": "stale"})
        tip_a = result["units"][0]["observations"]["tips"]["A"]
        tip_b = result["units"][0]["observations"]["tips"]["B"]
        self.assertEqual(tip_a["material"], [0])
        self.assertEqual(tip_a["surveyed_material"], ["copy_A"])
        self.assertEqual(tip_a["material_evidence"], {"copy_A": "assembly interval search v1"})
        self.assertEqual(tip_a["material_sensitivity"], {"copy_A": 0.8})
        self.assertNotIn("features", tip_a)
        self.assertEqual(tip_b["material"], [None])
        self.assertEqual(tip_b["material_evidence"], {"copy_A": "no callable locus evidence"})
        self.assertEqual(tip_b["material_sensitivity"], {})

        with tempfile.TemporaryDirectory() as directory:
            model_path = Path(directory) / "model.json"
            tree_path = Path(directory) / "tree.tsv"
            model_path.write_text(json.dumps(result), encoding="utf-8")
            tree_path.write_text(
                "node_id\tparent_id\tlabel\tbranch_length\n"
                "root\t\troot\t\nA\troot\tA\t1\nB\troot\tB\t1\n",
                encoding="utf-8")
            bundle = load_locus_model(model_path, tree_path)
        process = bundle.units[0].process
        material_state = {state.material[0]: index for index, state in enumerate(process.states)}
        self.assertAlmostEqual(bundle.units[0].tips["A"][material_state[1]], 0.2)
        self.assertAlmostEqual(bundle.units[0].tips["A"][material_state[2]], 0.9)
        self.assertTrue(all(value == 1.0 for value in bundle.units[0].tips["B"]))

    def test_missing_duplicate_or_unknown_rows_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "cover every declared row"):
            apply_material_evidence(_model(), _rows()[:1])
        with self.assertRaisesRegex(ValueError, "Duplicate DNA evidence row"):
            apply_material_evidence(_model(), _rows() + [_rows()[0]])
        extra = _rows() + [{**_rows()[0], "tip": "missing_tip"}]
        with self.assertRaisesRegex(ValueError, "undeclared unit, tip, or material"):
            apply_material_evidence(_model(), extra)

    def test_duplicate_declared_unit_is_rejected(self):
        model = _model()
        duplicate = json.loads(json.dumps(model["units"][0]))
        model["units"].append(duplicate)
        with self.assertRaisesRegex(ValueError, "Duplicate declared family/unit"):
            apply_material_evidence(model, _rows())

    def test_unknown_calls_cannot_carry_probabilities_and_known_calls_need_survey(self):
        rows = _rows()
        rows[1]["sensitivity"] = 0.8
        rows[1]["specificity"] = 0.9
        with self.assertRaisesRegex(ValueError, "Unknown calls"):
            apply_material_evidence(_model(), rows)
        rows = _rows()
        rows[0]["survey"] = False
        with self.assertRaisesRegex(ValueError, "require a surveyed material"):
            apply_material_evidence(_model(), rows)

    def test_probability_pairs_are_validated(self):
        rows = _rows()
        del rows[0]["specificity"]
        with self.assertRaisesRegex(ValueError, "supplied together"):
            apply_material_evidence(_model(), rows)
        rows = _rows()
        rows[0]["sensitivity"] = True
        with self.assertRaisesRegex(ValueError, "finite number"):
            apply_material_evidence(_model(), rows)


if __name__ == "__main__":
    unittest.main()
