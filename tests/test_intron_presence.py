import json
import tempfile
import unittest
from pathlib import Path

from intraphy.inference.dna_presence import analyze_dna_presence
from intraphy.inference.intron_presence import analyze_intron_positions
from intraphy.topology import SpeciesTree


def _tree():
    return SpeciesTree([
        {"node_id": "R", "parent_id": "", "label": "root", "branch_length": "NA"},
        {"node_id": "A", "parent_id": "R", "label": "A", "branch_length": 0.7},
        {"node_id": "B", "parent_id": "R", "label": "B", "branch_length": 1.2},
    ])


def _rows(states, observation_type, layer):
    return [{
        "family_id": "f", "site_id": "s", "species": species, "state": state,
        "eligible": True, "reason": "qualified", "evidence": "position evidence",
        "observation_type": observation_type, "layer": layer,
    } for species, state in states.items()]


class IntronPositionPresenceTests(unittest.TestCase):
    def test_relabelled_binary_data_have_same_likelihood_and_posteriors(self):
        rates = {"gain": 0.4, "loss": 0.25}
        dna_rows = _rows({"A": 0, "B": 1}, "homologous_dna_presence", "exon_presence")
        intron_rows = _rows({"A": "absent", "B": "present"}, "intron_position_presence", "intron_position")
        with tempfile.TemporaryDirectory() as directory:
            dna = analyze_dna_presence(
                dna_rows, _tree(), Path(directory) / "dna", root_frequency="fixed",
                root_presence=0.35, fixed_rates=rates,
            )
            intron = analyze_intron_positions(
                intron_rows, _tree(), Path(directory) / "intron", root_frequency="fixed",
                root_presence=0.35, fixed_rates=rates,
            )
        self.assertEqual(dna["families"][0]["log_likelihood"], intron["families"][0]["log_likelihood"])
        dna_nodes = dna["history"][0]["sites"][0]["nodes"]
        intron_nodes = intron["history"][0]["sites"][0]["nodes"]
        self.assertEqual(
            [(row["node_id"], row["p_absent"], row["p_present"]) for row in dna_nodes],
            [(row["node_id"], row["p_no_intron"], row["p_intron_present"])
             for row in intron_nodes],
        )
        self.assertEqual(intron["model"], "intron-position-ctmc")

    def test_output_schema_and_expected_counts_are_domain_specific(self):
        rows = _rows({"A": "absent", "B": "present"}, "intron_position_presence", "intron_position")
        with tempfile.TemporaryDirectory() as directory:
            result = analyze_intron_positions(
                rows, _tree(), directory, fixed_rates={"gain": 0.4, "loss": 0.25},
                expected_counts=True,
            )
            fit = json.loads((Path(directory) / "intron_fit.json").read_text())
            history = json.loads((Path(directory) / "intron_history.json").read_text())
            run_result = json.loads((Path(directory) / "run_result.json").read_text())
        self.assertEqual(result["model"], "intron-position-ctmc")
        self.assertEqual(fit["observation_type"], "intron_position_presence")
        self.assertIn("homologous coding position", fit["assumptions"]["states"]["1"])
        self.assertIn("not modeled", fit["assumptions"]["unobserved_all_zero_sites"])
        self.assertIn("no branch-specific or site-specific", fit["assumptions"]["rate_heterogeneity"])
        edge = history["families"][0]["sites"][0]["edges"][0]
        self.assertIn("expected_intron_position_gain_site_state_change_count", edge)
        self.assertNotIn("expected_gain_site_state_change_count", edge)
        self.assertEqual(run_result["model"], "intron-position-ctmc")
        self.assertIn("intron_history.json", run_result["artifacts"])

    def test_intron_type_is_required_and_other_domains_are_rejected(self):
        rows = _rows({"A": "absent", "B": "present"}, "intron_position_presence", "intron_position")
        with tempfile.TemporaryDirectory() as directory:
            missing_type = [dict(row) for row in rows]
            missing_type[0].pop("observation_type")
            with self.assertRaisesRegex(ValueError, "observation_type"):
                analyze_intron_positions(missing_type, _tree(), directory)
            dna_type = [dict(row, observation_type="homologous_dna_presence") for row in rows]
            with self.assertRaisesRegex(ValueError, "observation_type"):
                analyze_intron_positions(dna_type, _tree(), directory)
            legacy_layer = [dict(row, layer="splice_junction") for row in rows]
            with self.assertRaisesRegex(ValueError, "layer conflicts"):
                analyze_intron_positions(legacy_layer, _tree(), directory)

    def test_boolean_states_are_not_accepted_as_binary_states(self):
        rows = _rows({"A": True, "B": "present"}, "intron_position_presence", "intron_position")
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "invalid intron-position"):
                analyze_intron_positions(rows, _tree(), directory)


if __name__ == "__main__":
    unittest.main()
