"""Source-only contract tests for the strict finite locus-model input format."""
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from intraphy.inference.locus_io import load_locus_model
from intraphy.inference.locus_run import analyze_locus


def _unit(material_id="m", start=10, end=20, unit_id="u"):
    return {
        "family": "family", "unit": unit_id,
        "catalogue": {
            "provenance": "synthetic evidence",
            "state_model": "irreversible",
            "material": [{"id": material_id, "start": start, "end": end}],
            "copies": [{"id": "copy-" + material_id, "material_ids": [material_id]}],
            "opportunities": [{"id": "delete-" + material_id, "outcome_id": "delete",
                               "kind": "dna_deletion", "rate_group": "deletion",
                               "material_deletions": [material_id], "interval": [start, end]}],
        },
        "root": {"provenance": "declared root support", "entries": [
            {"state": {"material": [1]}, "weight": 1.0}]},
        "observations": {"provenance": "synthetic tip calls", "tips": {
            "A": {"material": [None]}, "B": {"material": [0]}}},
    }


def _model(units=None):
    return {
        "schema": "intraphy.exon-locus-model/2", "model": "exon-locus-ctmc",
        "branch_length_unit": "substitutions per declared opportunity",
        "tree_provenance": "synthetic rooted tree", "provenance": "synthetic model",
        "independence_provenance": "one declared locus unit",
        "rates": {"provenance": "explicit initial value", "groups": {
            "deletion": {"mode": "fit", "initial": 0.2}}},
        "units": [_unit()] if units is None else units,
    }


class LocusIOTests(unittest.TestCase):
    def _load(self, model=None, tree_text=None, model_text=None):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        directory = Path(temporary.name)
        model_path, tree_path = directory / "model.json", directory / "tree.tsv"
        if model_text is not None:
            model_path.write_text(model_text, encoding="utf-8")
        else:
            model_path.write_text(json.dumps(_model() if model is None else model), encoding="utf-8")
        tree_path.write_text(tree_text or (
            "node_id\tparent_id\tlabel\tbranch_length\n"
            "root\t\troot\t\nA\troot\tA\t0.2\nB\troot\tB\t0.3\n"), encoding="utf-8")
        return load_locus_model(model_path, tree_path)

    def test_optional_empty_collections_and_root_rate_declarations_load(self):
        bundle = self._load()
        self.assertEqual(bundle.rates["deletion"], {"mode": "fit", "value": 0.2})
        unit = bundle.units[0]
        self.assertEqual(len(unit.root_prior), len(unit.process.states))
        self.assertAlmostEqual(float(unit.root_prior.sum()), 1.0)
        self.assertEqual(tuple(unit.process.catalogue.opportunities[0].preconditions), ())
        self.assertEqual(unit.process.catalogue.opportunities[0].weight, 1.0)
        self.assertEqual(set(unit.tips), {"A", "B"})

    def test_non_array_collections_and_duplicate_json_keys_are_rejected(self):
        record = _model()
        record["units"][0]["catalogue"]["material"] = {}
        with self.assertRaisesRegex(ValueError, "collections must be arrays"):
            self._load(record)
        with self.assertRaisesRegex(ValueError, "duplicate JSON key"):
            self._load(model_text='{"schema":"x","schema":"y"}')

    def test_every_nonroot_tree_edge_requires_a_supplied_branch_length(self):
        tree = ("node_id\tparent_id\tlabel\tbranch_length\n"
                "root\t\troot\t\nA\troot\tA\t\nB\troot\tB\t0.3\n")
        with self.assertRaisesRegex(SystemExit, "lacks a branch length"):
            self._load(tree_text=tree)

    def test_missing_dna_is_uninformative_while_zero_accepts_absent_and_deleted(self):
        unit = self._load().units[0]
        states = unit.process.states
        np.testing.assert_array_equal(unit.tips["A"], np.ones(len(states)))
        supported_zero = [i for i, state in enumerate(states) if state.material[0] in (0, 2)]
        self.assertTrue(supported_zero)
        self.assertTrue(all(unit.tips["B"][i] == 1.0 for i in supported_zero))
        self.assertTrue(all(unit.tips["B"][i] == 0.0 for i, state in enumerate(states)
                            if state.material[0] == 1))

    def test_rate_groups_must_cover_declared_opportunities(self):
        record = _model()
        record["rates"]["groups"] = {"unrelated": {"mode": "fixed", "value": 0.1}}
        with self.assertRaisesRegex(ValueError, "Rate groups missing values"):
            self._load(record)

    def test_zero_event_fixed_evaluation_accepts_an_empty_rate_group_map(self):
        record = _model()
        record["rates"]["groups"] = {}
        record["units"][0]["catalogue"]["opportunities"] = []
        record["units"][0]["observations"]["tips"]["B"]["material"] = [None]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model_path, tree_path = root / "model.json", root / "tree.tsv"
            model_path.write_text(json.dumps(record), encoding="utf-8")
            tree_path.write_text("node_id\tparent_id\tlabel\tbranch_length\n"
                                 "root\t\troot\t\nA\troot\tA\t0.2\nB\troot\tB\t0.3\n", encoding="utf-8")
            bundle = load_locus_model(model_path, tree_path)
            self.assertEqual(bundle.rates, {})
            self.assertEqual(bundle.units[0].process.edges, ())
            analyze_locus(bundle, root / "out", parameter_mode="fixed")
            history = json.loads((root / "out" / "locus_history.json").read_text(encoding="utf-8"))
            self.assertEqual(history["log_likelihood"], 0.0)
            self.assertEqual(history["parameter_status"], "fixed_parameter_evaluation")

    def test_zero_weight_root_entry_and_malformed_precondition_pair_reject(self):
        record = _model()
        record["units"][0]["root"]["entries"].append({
            "state": {"material": [2]}, "weight": 0.0})
        with self.assertRaisesRegex(ValueError, "positive"):
            self._load(record)

        record = _model()
        record["units"][0]["catalogue"]["opportunities"][0]["preconditions"] = [["m"]]
        with self.assertRaisesRegex(ValueError, "two-item arrays"):
            self._load(record)

    def test_overlapping_material_in_same_family_is_not_an_independent_unit(self):
        second = _unit(material_id="m2", unit_id="u2")
        record = _model([_unit(), second])
        with self.assertRaisesRegex(ValueError, "overlap on the shared coordinate axis"):
            self._load(record)

    def test_v1_and_splice_fields_are_rejected(self):
        record = _model()
        record["schema"] = "intraphy.exon-locus-model/1"
        with self.assertRaisesRegex(ValueError, "DNA-only"):
            self._load(record)
        record = _model()
        record["units"][0]["catalogue"]["features"] = []
        with self.assertRaisesRegex(ValueError, "DNA catalogue"):
            self._load(record)
        record = _model()
        record["units"][0]["observations"]["tips"]["A"]["observed_paths"] = []
        with self.assertRaisesRegex(ValueError, "unknown fields"):
            self._load(record)
        record = _model()
        record["units"][0]["root"]["entries"][0]["state"]["active_features"] = []
        with self.assertRaisesRegex(ValueError, "root state"):
            self._load(record)
        record = _model()
        record["units"][0]["catalogue"]["opportunities"][0]["feature_off"] = []
        with self.assertRaisesRegex(ValueError, "event opportunity"):
            self._load(record)

    def test_state_model_is_explicit(self):
        record = _model()
        del record["units"][0]["catalogue"]["state_model"]
        with self.assertRaisesRegex(ValueError, "state_model"):
            self._load(record)

    def test_binary_state_model_has_only_absent_and_present(self):
        record = _model()
        record["units"][0]["catalogue"]["state_model"] = "binary"
        unit = self._load(record).units[0]
        self.assertEqual(unit.process.catalogue.material_states, (0, 1))
        self.assertEqual({state.material[0] for state in unit.process.states}, {0, 1})


if __name__ == "__main__":
    unittest.main()
