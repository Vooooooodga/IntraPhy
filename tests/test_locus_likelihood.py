from types import SimpleNamespace
import unittest

import numpy as np
from scipy.linalg import expm
from scipy.sparse.linalg import expm_multiply

from intraphy.inference.locus_likelihood import evaluate_locus, process_generator
from intraphy.structure.locus_types import (
    CopySlot,
    LocusCatalogue,
    LocusObservation,
    LocusState,
    MaterialTract,
    ProcessEdge,
    SpliceFeature,
)
from intraphy.topology import SpeciesTree


def _process(edges, *, groups=("event",), states=None, catalogue=None):
    states = states or tuple(LocusState((value,), frozenset()) for value in (0, 1))
    if catalogue is None:
        opportunities = tuple(SimpleNamespace(id=group, rate_group=group) for group in groups)
        catalogue = SimpleNamespace(opportunities=opportunities)
    return SimpleNamespace(catalogue=catalogue, states=states, edges=tuple(edges))


def _tree(rows):
    return SpeciesTree(rows)


class LocusLikelihoodTests(unittest.TestCase):
    def test_sparse_exponential_action_matches_dense_asymmetric_generator(self):
        edges = (
            ProcessEdge(0, 1, "a", "a1", "gain", "gain", 1.0),
            ProcessEdge(1, 0, "b", "b1", "loss", "loss", 1.0),
            ProcessEdge(1, 2, "c", "c1", "splice", "splice", 0.7),
        )
        process = _process(edges, groups=("gain", "loss", "splice"),
                           states=tuple(LocusState((i,), frozenset()) for i in range(3)))
        rates = {"gain": 0.2, "loss": 0.8, "splice": 0.3}
        q = process_generator(process, rates)
        vector = np.array([0.3, 0.2, 0.5])
        time = 0.7
        sparse = expm_multiply(q * time, vector)
        dense = expm(q.toarray() * time) @ vector
        np.testing.assert_allclose(sparse, dense, rtol=1e-12, atol=1e-12)

    def test_pruning_and_node_posteriors_match_dense_hidden_state_enumeration(self):
        states = tuple(LocusState((value,), frozenset()) for value in (0, 1, 2))
        process = _process((
            ProcessEdge(0, 1, "gain", "g1", "splice_change", "gain", 1.0),
            ProcessEdge(1, 2, "loss", "l1", "dna_deletion", "loss", 0.7),
            ProcessEdge(2, 0, "return", "r1", "splice_change", "return", 0.4),
            ProcessEdge(1, 0, "reverse", "v1", "splice_change", "reverse", 0.2),
        ), groups=("gain", "loss", "return", "reverse"), states=states)
        rows = [
            {"node_id": "root", "parent_id": "", "label": "root"},
            {"node_id": "inner", "parent_id": "root", "label": "inner", "branch_length": 0.31},
            {"node_id": "outside", "parent_id": "root", "label": "outside", "branch_length": 0.23},
            {"node_id": "left", "parent_id": "inner", "label": "left", "branch_length": 0.47},
            {"node_id": "right", "parent_id": "inner", "label": "right", "branch_length": 0.62},
        ]
        tree = _tree(rows)
        tips = {"left": [0.0, 1.0, 0.2], "right": [0.3, 0.0, 1.0],
                "outside": [1.0, 0.4, 0.0]}
        root_prior = np.array([0.25, 0.5, 0.25])
        rates = {"gain": 0.8, "loss": 0.3, "return": 0.15, "reverse": 0.4}
        q = process_generator(process, rates).toarray()
        transitions = {child: expm(q * tree.branch_length(child))
                       for _, child in tree.edges()}
        joint = np.zeros((3, 3))
        for root_state in range(3):
            for inner_state in range(3):
                for left_state in range(3):
                    for right_state in range(3):
                        for outside_state in range(3):
                            joint[root_state, inner_state] += (
                                root_prior[root_state]
                                * transitions["inner"][root_state, inner_state]
                                * transitions["left"][inner_state, left_state] * tips["left"][left_state]
                                * transitions["right"][inner_state, right_state] * tips["right"][right_state]
                                * transitions["outside"][root_state, outside_state] * tips["outside"][outside_state]
                            )
        total = joint.sum()
        expected_root = joint.sum(axis=1) / total
        expected_inner = joint.sum(axis=0) / total
        result = evaluate_locus(process, tree, tips, root_prior, rates)
        self.assertAlmostEqual(result.log_likelihood, np.log(total), places=11)
        np.testing.assert_allclose(result.node_posteriors["root"], expected_root, rtol=1e-10, atol=1e-12)
        np.testing.assert_allclose(result.node_posteriors["inner"], expected_inner, rtol=1e-10, atol=1e-12)

    def test_irreversible_deletion_likelihood_mle_and_marked_counts(self):
        process = _process((ProcessEdge(1, 0, "deletion", "delete", "dna_deletion", "deletion", 1.0),),
                           groups=("deletion",),
                           states=(LocusState((2,), frozenset()), LocusState((1,), frozenset())))
        tree = _tree([
            {"node_id": "root", "parent_id": "", "label": "root"},
            {"node_id": "present", "parent_id": "root", "label": "present", "branch_length": 1},
            {"node_id": "absent", "parent_id": "root", "label": "absent", "branch_length": 1},
        ])
        tips = {"present": [0, 1], "absent": [1, 0]}
        root_prior = [0, 1]
        evaluated = evaluate_locus(process, tree, tips, root_prior,
                                   {"deletion": 0.4}, counts=True)
        expected_ll = -0.4 + np.log1p(-np.exp(-0.4))
        self.assertAlmostEqual(evaluated.log_likelihood, expected_ll)
        self.assertAlmostEqual(evaluated.branch_event_counts["absent"]["opportunity:deletion"], 1.0)
        self.assertAlmostEqual(evaluated.branch_event_counts["present"]["opportunity:deletion"], 0.0)

        from intraphy.inference.locus_rates import LocusFitUnit, fit_locus_rates

        fitted = fit_locus_rates([LocusFitUnit(process, tree, tips, root_prior)],
                                 {"deletion": 0.4}, maxiter=500)
        self.assertTrue(fitted.optimizer_success)
        self.assertAlmostEqual(fitted.rates["deletion"], np.log(2), delta=2e-3)

    def test_zero_likelihood_is_reported_as_impossible_and_root_prior_is_not_renormalized(self):
        process = _process((), groups=())
        tree = _tree([{"node_id": "tip", "parent_id": "", "label": "tip"}])
        impossible = evaluate_locus(process, tree, {"tip": [1, 0]}, [0, 1], {})
        self.assertEqual(impossible.log_likelihood, -np.inf)
        self.assertIsNone(impossible.node_posteriors)
        with self.assertRaisesRegex(ValueError, "sums to one"):
            evaluate_locus(process, tree, {"tip": [0, 1]}, [0, 2], {})

    def test_detector_emission_uses_surveyed_feature_sensitivity_and_specificity(self):
        catalogue = LocusCatalogue(
            material=(MaterialTract("m", 0, 10),),
            copies=(CopySlot("c", ("m",)),),
            features=(SpliceFeature("f", "exon", ("m",), copy_id="c", start=0, end=10),),
            opportunities=(), provenance="unit test",
        )
        states = (LocusState((1,), frozenset()), LocusState((1,), frozenset({"f"})))
        process = _process((), groups=(), states=states, catalogue=catalogue)
        observation = LocusObservation((1,), (("f", 1),), frozenset({"f"}))
        tree = _tree([{"node_id": "tip", "parent_id": "", "label": "tip"}])
        result = evaluate_locus(process, tree, {"tip": observation}, [0.5, 0.5], {},
                                feature_sensitivity={"f": 0.8},
                                feature_specificity={"f": 0.9})
        self.assertAlmostEqual(result.log_likelihood, np.log(0.45))
        self.assertAlmostEqual(result.node_posteriors["tip"][1], 0.8 / 0.9)

    def test_unknown_tips_on_high_degree_polytomy_have_unit_likelihood(self):
        rows = [{"node_id": "root", "parent_id": "", "label": "root"}]
        tips = {}
        for index in range(12):
            label = f"tip{index}"
            rows.append({"node_id": label, "parent_id": "root", "label": label, "branch_length": 0.5})
            tips[label] = [1, 1]
        tree = _tree(rows)
        process = _process((ProcessEdge(1, 0, "deletion", "delete", "dna_deletion", "event", 1),))
        result = evaluate_locus(process, tree, tips, [0, 1], {"event": 0.8})
        self.assertAlmostEqual(result.log_likelihood, 0.0, delta=1e-12)
