"""Regression tests for declared coexisting observation scenarios (V19)."""
import json
import itertools
import tempfile
import unittest
from pathlib import Path
import numpy as np

from intraphy.inference.configuration_history import reconstruct
from intraphy.inference.configuration_run import infer_configurations
from intraphy.commands.exons import add_configuration_options
from intraphy.topology import SpeciesTree
from intraphy.storage.tabular import write_tsv
from intraphy.structure.serialization import write_catalogues, write_json, json_safe
from intraphy.structure.types import (Catalogue, ExonConfiguration as C, ExonSpan as E,
                                      ObservationEvidence as O, ConfigurationAlternative)
from intraphy.structure.space import enumerate_space
from intraphy.structure.observations import observation_scenarios, compatibility
from intraphy.structure.edits import EDIT_KINDS


def _catalogue():
    empty, full = C(()), C((E(0, 10),))
    observations = tuple(O(s, (empty, full), "coexisting") for s in "ABCDEFG")
    return Catalogue("g", "u", 10, (E(0, 10),), (), observations=observations,
                     boundary_candidates=(E(0, 10),))


def _tree_rows():
    return [{"node_id": "r", "parent_id": "", "label": "r"}] + [
        {"node_id": s, "parent_id": "r", "label": s, "branch_length": .2} for s in "ABCDEFG"]


class ConfigurationScenarioTests(unittest.TestCase):
    def test_independent_product_has_exact_128_scenarios_and_matches_eager_reference(self):
        catalogue = _catalogue(); space = enumerate_space(catalogue); taxa = tuple("ABCDEFG")
        scenarios = observation_scenarios(space, taxa, "annotation")
        self.assertEqual(scenarios.theoretical_count, 128)
        eager = []
        for choices in itertools.product(*(o.configurations for o in catalogue.observations)):
            selected = dict(zip(taxa, choices))
            tips = {s: compatibility(space, catalogue.observations[i], "annotation", selected[s])
                    for i, s in enumerate(taxa)}
            label = ";".join(f"{s}={selected[s].key}" for s in sorted(selected))
            eager.append((label, tips))
        self.assertEqual(len(eager), 128)
        self.assertEqual([label for label, _ in scenarios], [label for label, _ in eager])
        for (_, lazy), (_, reference) in zip(scenarios, eager):
            for species in taxa:
                np.testing.assert_array_equal(lazy[species], reference[species])

    def test_infer_configurations_reports_all_views_and_selected_counts(self):
        with tempfile.TemporaryDirectory() as directory:
            root, out = Path(directory) / "in", Path(directory) / "out"
            root.mkdir()
            write_tsv(root / "species_tree.tsv", _tree_rows(), ["node_id", "parent_id", "label", "branch_length"])
            path = root / "catalogues.jsonl"; write_catalogues(path, [_catalogue()])
            returned = infer_configurations(root, out, configurations=path, model="exon-parsimony",
                                            observation_view="annotation")
            detail = json.loads((out / "model_diagnostics.json").read_text())
            info = detail["units"][0]
            for view in ("annotation", "evidence"):
                self.assertEqual(info["observation_scenarios"][view]["theoretical_count"], 128)
                self.assertEqual(info["observation_scenarios"][view]["evaluated_count"], 128)
            self.assertEqual(detail["observation_view"], "annotation")
            self.assertEqual(info["observation_scenarios"]["annotation"]["resolved_count"], 128)
            history = json.loads((out / "exon_history.json").read_text())["units"][0]
            self.assertEqual(len(history["views"]["annotation"]["histories"]), 128)
            self.assertEqual(len(history["views"]["evidence"]["histories"]), 128)
            self.assertEqual(returned[0]["observation_scenario_theoretical_count"], 128)
            self.assertEqual(returned[0]["observation_scenario_evaluated_count"], 128)
            reference_space = enumerate_space(_catalogue())
            reference_tree = SpeciesTree(_tree_rows())
            eager = list(observation_scenarios(reference_space, tuple("ABCDEFG"), "annotation"))
            for scenario, (label, reference_tips) in zip(history["views"]["annotation"]["histories"], eager):
                self.assertEqual(scenario["observation_scenario"], label)
                reference = reconstruct(reference_space, reference_tree, reference_tips)
                for field in ("minimum_cost", "events", "nodes", "witness"):
                    self.assertEqual(json_safe(scenario.get(field)), json_safe(reference.get(field)))

    def test_explicit_scenario_cap_fails_without_partial_histories(self):
        catalogue = Catalogue("g", "u", 10, (E(0, 10),), (),
                              observations=tuple(O(s, (C(()), C((E(0, 10),))), "coexisting") for s in "ABCDEFG"),
                              boundary_candidates=(E(0, 10),))
        with self.assertRaisesRegex(ValueError, "observation_scenarios_incomplete"):
            list(observation_scenarios(enumerate_space(catalogue), tuple("ABCDEFG"), "annotation", 64))

    def test_explicit_scenario_cap_integration_reports_zero_evaluated_and_no_histories(self):
        with tempfile.TemporaryDirectory() as directory:
            root, out = Path(directory) / "in", Path(directory) / "out"
            root.mkdir()
            write_tsv(root / "species_tree.tsv", _tree_rows(), ["node_id", "parent_id", "label", "branch_length"])
            path = root / "catalogues.jsonl"; write_catalogues(path, [_catalogue()])
            returned = infer_configurations(root, out, configurations=path, model="exon-parsimony",
                                            observation_view="annotation", max_observation_scenarios=64)
            detail = json.loads((out / "model_diagnostics.json").read_text())["units"][0]
            selected = detail["observation_scenarios"]["annotation"]
            self.assertEqual(selected["theoretical_count"], 128)
            self.assertEqual(selected["evaluated_count"], 0)
            history = json.loads((out / "exon_history.json").read_text())["units"][0]
            self.assertEqual(history["views"], {})
            summary = (out / "exon_structure_summary.tsv").read_text()
            self.assertIn("minimum_structural_edits", summary)
            self.assertIsNone(returned[0]["minimum_structural_edits"])
            self.assertIn("incomplete", returned[0]["status"])

    def test_ctmc_coexisting_is_explicitly_withheld(self):
        with tempfile.TemporaryDirectory() as directory:
            root, out = Path(directory) / "in", Path(directory) / "out"
            root.mkdir()
            write_tsv(root / "species_tree.tsv", _tree_rows(), ["node_id", "parent_id", "label", "branch_length"])
            path = root / "catalogues.jsonl"; write_catalogues(path, [_catalogue()])
            rates = root / "rates.json"
            write_json(rates, {"schema": "intraphy.exon-rates/1", "provenance": "synthetic regression fixture",
                               "rates": {kind: .1 for kind in EDIT_KINDS}})
            infer_configurations(root, out, configurations=path, model="exon-ctmc", rates=rates,
                                 observation_view="annotation", max_observation_scenarios=128)
            diagnostics = json.loads((out / "model_diagnostics.json").read_text())["units"][0]
            self.assertEqual(diagnostics["probability_status"], "coexisting_structures_are_not_a_probability_mixture")
            history = json.loads((out / "exon_history.json").read_text())["units"][0]
            self.assertNotIn("ctmc", history)

    def test_scenario_input_validation_and_alternative_links(self):
        catalogue = _catalogue(); space = enumerate_space(catalogue)
        for limit in (0, -1, True, 1.5):
            with self.assertRaises(ValueError):
                observation_scenarios(space, tuple("ABCDEFG"), "annotation", limit)
        self.assertEqual(len(observation_scenarios(space, tuple("ABCDEFG"), "annotation", None)), 128)
        self.assertEqual(observation_scenarios(space, tuple("ABCDEFG"), "annotation").limit, None)
        incompatible = O("A", (C(()), C((E(0, 10),))), "coexisting")
        with self.assertRaises(ValueError):
            compatibility(space, incompatible, "annotation")
        with self.assertRaises(ValueError):
            compatibility(space, incompatible, "invalid", C(()))
        partial = O("A", (C(()),), "partial", alternative_exons=(E(0, 10),))
        self.assertEqual(int(compatibility(space, partial, "annotation").sum()), 1)
        self.assertEqual(int(compatibility(space, partial, "evidence").sum()), 2)
        bad = Catalogue("g", "u", 10, (E(0, 10),), (),
                        observations=(O("A", (C((E(1, 2),)), C((E(0, 10),))), "coexisting"),),
                        boundary_candidates=(E(0, 10),))
        with self.assertRaises(ValueError):
            observation_scenarios(enumerate_space(bad), ("A",), "annotation")

    def test_native_linked_alternative_is_applied_only_to_its_replaced_configuration(self):
        empty, full = C(()), C((E(0, 10),))
        alternative = ConfigurationAlternative(full, empty.key, "B", "tx", "supported", 1.)
        c = Catalogue("g", "u", 10, (E(0, 10),), (),
                      observations=(O("A", (empty, full), "coexisting", alternatives=(alternative,)),
                                    O("B", (full,))), boundary_candidates=(E(0, 10),))
        space = enumerate_space(c)
        self.assertEqual(int(compatibility(space, c.observations[0], "evidence", empty).sum()), 2)
        self.assertEqual(int(compatibility(space, c.observations[0], "evidence", full).sum()), 1)
        labels = [label for label, _ in observation_scenarios(space, ("A", "B"), "evidence")]
        self.assertEqual(len(labels), 2)

    def test_parser_default_is_uncapped(self):
        import argparse
        parser = argparse.ArgumentParser()
        add_configuration_options(parser)
        self.assertIsNone(parser.parse_args([]).max_observation_scenarios)

    def test_event_aggregation_handles_missing_events_as_zero(self):
        from intraphy.inference.configuration_run import _robust_events
        a = {"parent": "r", "child": "A", "edit_key": "split:0:1:", "minimum_count": 1,
             "maximum_count": 1, "support": "required", "consequences": [], "affected_spans": []}
        result = _robust_events([{"events": [a]}, {"events": []}])
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0]["minimum_count"], 0)
        self.assertEqual(result[0]["maximum_count"], 1)
        self.assertEqual(result[0]["support"], "possible")


if __name__ == "__main__":
    unittest.main()
