import json
import math
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from intraphy.inference.ctmc import _ascertainment_log_probability, _pattern_log_likelihood
from intraphy.inference.dna_fit import _site_log_likelihood, fit_dna_family
from intraphy.inference.dna_presence import analyze_dna_presence
from intraphy.topology import SpeciesTree


def _tree(reverse=False):
    rows = [
        {"node_id": "R", "parent_id": "", "label": "root", "branch_length": "NA"},
        {"node_id": "A", "parent_id": "R", "label": "A", "branch_length": 0.7},
        {"node_id": "B", "parent_id": "R", "label": "B", "branch_length": 1.2},
    ]
    return SpeciesTree(rows[::-1] if reverse else rows)


def _rows(states, *, family="f", site="s", eligible=True):
    return [
        {"family_id": family, "site_id": site, "species": species, "state": state,
         "eligible": eligible, "reason": "", "evidence": ""}
        for species, state in states.items()
    ]


class DnaPresenceTests(unittest.TestCase):
    def test_two_tip_pattern_likelihood_matches_hand_calculation(self):
        tree = _tree()
        gain, loss = 0.3, 0.2
        pi0, pi1 = loss / (gain + loss), gain / (gain + loss)

        def transition(start, end, length):
            changed = 1 - math.exp(-(gain + loss) * length)
            if start == end == 0:
                return pi0 + pi1 * (1 - changed)
            if start == 0 and end == 1:
                return pi1 * changed
            if start == 1 and end == 0:
                return pi0 * changed
            return pi1 + pi0 * (1 - changed)

        expected = pi0 * transition(0, 0, 0.7) * transition(0, 1, 1.2)
        expected += pi1 * transition(1, 0, 0.7) * transition(1, 1, 1.2)
        observed = _pattern_log_likelihood(
            tree, {"A": 0, "B": 1}, gain, loss, 1.0, frozenset(), None
        )
        self.assertAlmostEqual(observed, math.log(expected))

    def test_unknown_tip_is_neutral_and_ascertainment_uses_observed_mask(self):
        tree = _tree()
        gain, loss = 0.3, 0.2
        observations = {"A": 1, "B": "unknown"}
        rho = None
        likelihood = _pattern_log_likelihood(tree, observations, gain, loss, 1.0, frozenset(), rho)
        denominator = _ascertainment_log_probability(
            tree, observations, gain, loss, 1.0, frozenset(), rho, "observed-at-least-one"
        )
        self.assertAlmostEqual(likelihood, math.log(gain / (gain + loss)))
        self.assertAlmostEqual(_site_log_likelihood(tree, observations, gain, loss, rho), 0.0)
        self.assertAlmostEqual(denominator, math.log(gain / (gain + loss)))

    def test_single_known_present_is_unit_likelihood_but_multitip_pattern_is_not(self):
        tree = _tree()
        self.assertAlmostEqual(_site_log_likelihood(
            tree, {"A": 1, "B": "unknown"}, 0.3, 0.2, None
        ), 0.0)
        self.assertLess(_site_log_likelihood(tree, {"A": 1, "B": 1}, 0.3, 0.2, None), 0.0)

    def test_fixed_zero_rates_preserve_diagnostic_output_without_impossible_history(self):
        with tempfile.TemporaryDirectory() as directory:
            result = analyze_dna_presence(
                _rows({"A": 0, "B": 1}), _tree(), directory,
                fixed_rates={"gain": 0.0, "loss": 0.0},
            )
            family = result["families"][0]
            self.assertEqual(family["status"], "optimizer_unresolved")
            self.assertIsNone(family["log_likelihood"])
            history = json.loads((Path(directory) / "dna_history.json").read_text())
            self.assertFalse(history["families"][0]["posterior_available"])
            self.assertTrue((Path(directory) / "run_result.json").exists())

    def test_all_finite_but_unsuccessful_starts_return_diagnostics(self):
        def unsuccessful_minimize(_objective, start, **_kwargs):
            return SimpleNamespace(success=False, status=1, message="iteration limit", nit=1, x=start)

        with patch("intraphy.inference.dna_fit.minimize", side_effect=unsuccessful_minimize):
            result = fit_dna_family(
                _tree(), {"s": {"A": 0, "B": 1}}, root_frequency="stationary",
                root_presence=0.5, fixed_rates=None,
            )
        self.assertFalse(result["optimizer_success"])
        self.assertTrue(math.isfinite(result["log_likelihood"]))
        self.assertEqual(len(result["starts"]), 3)

    def test_result_is_invariant_to_species_and_tree_row_order(self):
        rows = _rows({"A": 0, "B": 1}, site="variable") + _rows(
            {"A": 1, "B": 1}, site="constant"
        )
        with tempfile.TemporaryDirectory() as directory:
            first = analyze_dna_presence(rows, _tree(), Path(directory) / "first")
            second = analyze_dna_presence(
                list(reversed(rows)), _tree(reverse=True), Path(directory) / "second"
            )
        for rate in ("gain_rate", "loss_rate"):
            self.assertAlmostEqual(first["families"][0][rate], second["families"][0][rate])
        self.assertEqual(first["families"][0]["n_constant_present_sites"], 1)

    def test_no_observed_contrast_has_no_estimate_or_history(self):
        with tempfile.TemporaryDirectory() as directory:
            result = analyze_dna_presence(_rows({"A": 1, "B": 1}), _tree(), directory)
            family = result["families"][0]
            self.assertEqual(family["status"], "no_observed_contrast")
            self.assertIsNone(family["gain_rate"])
            self.assertIsNone(family["loss_rate"])
            history = json.loads((Path(directory) / "dna_history.json").read_text())
            self.assertFalse(history["families"][0]["posterior_available"])

    def test_empty_candidate_input_is_explicitly_unresolved(self):
        with tempfile.TemporaryDirectory() as directory:
            result = analyze_dna_presence([], _tree(), directory)
            self.assertEqual(result["families"], [])
            fit = json.loads((Path(directory) / "dna_fit.json").read_text())
            manifest = json.loads((Path(directory) / "run_result.json").read_text())
        self.assertEqual(fit["analysis_status"], "no_candidate_sites")
        self.assertEqual(manifest["status"], "completed_with_unresolved")
        self.assertEqual(manifest["diagnostics"], ["no_candidate_sites"])

    def test_expected_site_state_changes_are_not_molecular_lesion_counts(self):
        rows = _rows({"A": 0, "B": 1}, site="linked_a") + _rows(
            {"A": 0, "B": 1}, site="linked_b"
        )
        with tempfile.TemporaryDirectory() as directory:
            result = analyze_dna_presence(
                rows, _tree(), directory, expected_counts=True,
                fixed_rates={"gain": 0.4, "loss": 0.3},
            )
            self.assertEqual(result["families"][0]["status"], "fixed_rates")
            fit = json.loads((Path(directory) / "dna_fit.json").read_text())
            history = json.loads((Path(directory) / "dna_history.json").read_text())
        self.assertIn("composite likelihood", fit["assumptions"]["independence"])
        self.assertIn("not counts of molecular lesions", fit["assumptions"]["event_interpretation"])
        edges = history["families"][0]["sites"][0]["edges"]
        self.assertIn("expected_gain_site_state_change_count", edges[0])
        self.assertNotIn("molecular_lesion_count", edges[0])

    def test_conflicting_eligibility_and_missing_species_are_rejected(self):
        rows = _rows({"A": 0, "B": 1})
        rows[1]["eligible"] = False
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "conflicting eligible"):
                analyze_dna_presence(rows, _tree(), Path(directory) / "conflict")
            with self.assertRaisesRegex(ValueError, "tip mismatch"):
                analyze_dna_presence(_rows({"A": 0}), _tree(), Path(directory) / "missing")


if __name__ == "__main__":
    unittest.main()
