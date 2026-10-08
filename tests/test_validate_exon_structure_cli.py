import importlib.util
import csv
import json
from pathlib import Path
import tempfile
import unittest


TOOL_PATH = Path(__file__).resolve().parents[1] / "tools" / "validate_exon_structure.py"
SPEC = importlib.util.spec_from_file_location("validate_exon_structure", TOOL_PATH)
VALIDATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VALIDATOR)


def write_table(path, fields, rows):
    path.write_text("\t".join(fields) + "\n" + "\n".join(
        "\t".join(str(row.get(field, "")) for field in fields) for row in rows) + "\n",
        encoding="utf-8")


def native_result(directory, probability_delta=0.0, branch_delta=None,
                  pair_delta=0.0, pair_total_delta=0.0):
    if branch_delta is None:
        branch_delta = probability_delta
    directory.mkdir()
    artifacts = ["exon_structure_fit.json", "exon_history.json", "exon_structure_summary.tsv",
                 "ancestral_exon_states.tsv", "branch_exon_changes.tsv", "model_diagnostics.json"]
    (directory / "run_result.json").write_text(json.dumps({"status": "completed",
        "model": "exon-structure-ctmc", "artifacts": artifacts}), encoding="utf-8")
    (directory / "exon_structure_fit.json").write_text(json.dumps({"model": "exon-structure-ctmc",
        "parameter_mode": "fit"}), encoding="utf-8")
    (directory / "exon_history.json").write_text(json.dumps({"model": "exon-structure-ctmc",
        "observation_unit": "genomic_exon_spans"}), encoding="utf-8")
    (directory / "model_diagnostics.json").write_text(json.dumps({"model": "exon-structure-ctmc"}), encoding="utf-8")
    write_table(directory / "exon_structure_summary.tsv",
        ["family_id", "unit_id", "status", "log_likelihood", "unknown_tips"],
        [{"family_id": "f", "unit_id": "u", "status": "completed", "log_likelihood": "-1.2", "unknown_tips": "[]"}])
    write_table(directory / "ancestral_exon_states.tsv",
        ["family_id", "unit_id", "node", "state_id", "probability"],
        [{"family_id": "f", "unit_id": "u", "node": "root", "state_id": "C0000", "probability": .25 + probability_delta},
         {"family_id": "f", "unit_id": "u", "node": "root", "state_id": "C0001", "probability": .75 - probability_delta}])
    branch_fields = ["family_id", "unit_id", "parent_node_id", "child_node_id", "parent_exons",
        "child_exons", "parent_dna_presence", "child_dna_presence", "material_ids",
        "changed_material_tracts", "probability_exon_structure_change",
        "probability_dna_presence_change", "joint_configuration_probability",
        "probability_exon_count_increase", "probability_exon_count_decrease",
        "probability_exon_count_unchanged", "exon_count_pair_probabilities"]
    write_table(directory / "branch_exon_changes.tsv", branch_fields,
        [{"family_id": "f", "unit_id": "u", "parent_node_id": "root", "child_node_id": "A",
          "parent_exons": "[[0,10]]", "child_exons": "[[0,10]]", "parent_dna_presence": "[]",
          "child_dna_presence": "[]", "material_ids": "[]", "changed_material_tracts": "[]",
          "probability_exon_structure_change": .2 + branch_delta,
          "probability_dna_presence_change": 0., "joint_configuration_probability": .8 - probability_delta,
          "probability_exon_count_increase": 0., "probability_exon_count_decrease": 0.,
          "probability_exon_count_unchanged": 1.,
          "exon_count_pair_probabilities": json.dumps([
              {"parent_count": 1, "child_count": 1, "probability": .5 - pair_delta},
              {"parent_count": 1, "child_count": 2,
               "probability": .5 + pair_delta - pair_total_delta}])}])
    return VALIDATOR._validate_native_result(directory)


class ValidateExonStructureContractTests(unittest.TestCase):
    def test_probability_validation_rejects_nonfinite_and_out_of_range_values(self):
        self.assertEqual(VALIDATOR._number("0.25", "p", "fixture"), .25)
        for invalid in ("nan", "inf", "-0.01", "1.01"):
            with self.subTest(invalid=invalid), self.assertRaises(AssertionError):
                VALIDATOR._number(invalid, "p", "fixture")

    def test_native_result_requires_nonempty_tables_and_normalized_node_posteriors(self):
        with tempfile.TemporaryDirectory() as temporary:
            result = native_result(Path(temporary) / "result")
            self.assertEqual(len(result["branches"]), 1)
            self.assertEqual(len(result["ancestral"]), 2)

    def test_native_result_rejects_non_normalized_count_pair_probabilities(self):
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(AssertionError, "exon-count-pair probabilities do not normalize"):
                native_result(Path(temporary) / "result", pair_total_delta=.01)

    def test_native_result_rejects_missing_key_branch_probabilities(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary) / "result"
            native_result(directory)
            path = directory / "branch_exon_changes.tsv"
            with path.open(encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                fields, rows = reader.fieldnames, list(reader)
            for field in ("probability_exon_structure_change", "probability_dna_presence_change",
                          "joint_configuration_probability", "exon_count_pair_probabilities"):
                rows[0][field] = ""
            write_table(path, fields, rows)
            with self.assertRaisesRegex(AssertionError, "required native branch probabilities are missing"):
                VALIDATOR._validate_native_result(directory)

    def test_reload_comparison_uses_shared_probability_fields_with_tolerance(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            default = native_result(root / "default")
            reload_result = native_result(root / "reload", probability_delta=1e-10)
            VALIDATOR._compare_native_results(default, reload_result)
            changed = native_result(root / "changed", probability_delta=2e-8)
            with self.assertRaisesRegex(AssertionError, "Reload ancestral probability differs"):
                VALIDATOR._compare_native_results(default, changed)
            changed_branch = native_result(root / "changed_branch", branch_delta=2e-8)
            with self.assertRaisesRegex(AssertionError, "Reload probability_exon_structure_change differs"):
                VALIDATOR._compare_native_results(default, changed_branch)
            changed_pair = native_result(root / "changed_pair", pair_delta=2e-8)
            with self.assertRaisesRegex(AssertionError, "Reload exon-count-pair probability differs"):
                VALIDATOR._compare_native_results(default, changed_pair)

    def test_local_unknown_intervals_are_distinct_from_whole_tip_unknowns(self):
        catalogues = [{"observations": [
            {"species": "D", "kind": "partial", "unknown_intervals": [{"start": 5, "end": 8}]},
            {"species": "A", "kind": "unknown", "unknown_intervals": []},
        ]}]
        local = VALIDATOR._localized_unknown_observations(catalogues, "D")
        self.assertEqual(local[0]["kind"], "partial")
        self.assertTrue(local[0]["unknown_intervals"])
        self.assertEqual(VALIDATOR._localized_unknown_observations(
            [{"observations": [{"species": "D", "kind": "unknown", "unknown_intervals": [{"start": 0, "end": 1}]}]}], "D")[0]["species"], "D")


if __name__ == "__main__":
    unittest.main()
