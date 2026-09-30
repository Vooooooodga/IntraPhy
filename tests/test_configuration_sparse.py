"""Exact likelihood-only sparse CTMC backend regression contracts."""
import unittest
from unittest.mock import patch
import numpy as np
from scipy.linalg import expm
from scipy.sparse import csr_matrix

from intraphy.inference.configuration_ctmc import transition_matrix, likelihood
from intraphy.inference.configuration_model import RateModel, evaluate_model, generator
from intraphy.inference.configuration_sparse import likelihood_only, log_action, sparse_generator
from intraphy.inference.kernel_cache import KernelCache
from intraphy.structure.edits import EDIT_KINDS
from intraphy.structure.space import enumerate_space
from intraphy.structure.types import Catalogue, ExonSpan, Material
from intraphy.topology import SpeciesTree


def make_tree(length_a=0.2, length_b=0.3):
    return SpeciesTree([
        {"node_id": "r", "parent_id": "", "label": "r"},
        {"node_id": "a", "parent_id": "r", "label": "A", "branch_length": length_a},
        {"node_id": "b", "parent_id": "r", "label": "B", "branch_length": length_b},
    ])


def make_space(n=1, material=False):
    spans = tuple(ExonSpan(3*i, 3*i+2) for i in range(n))
    catalogue = Catalogue("f", "u", 3*n, spans, (),
        material=(Material("m", 0, 1),) if material else (),
        boundary_candidates=spans)
    return enumerate_space(catalogue)


def make_material_space(n=4):
    materials = tuple(Material(f"m{i}", 2*i, 2*i+1) for i in range(n))
    return enumerate_space(Catalogue("f", "u", 2*n, (), (),
                                     material=materials, boundary_candidates=()))


def model(**kwargs):
    return RateModel({kind: kwargs.get(kind, 0.05) for kind in EDIT_KINDS},
                     foreground_multiplier=kwargs.get("foreground_multiplier", 1.),
                     foreground=frozenset(kwargs.get("foreground", ())))


class SparseLikelihoodTests(unittest.TestCase):
    def test_sparse_generator_matches_dense_generator(self):
        space = make_space(material=True)
        rates = model()
        origins = {"m": "a"}
        dense, _ = generator(space, rates, origins, "a", include_marks=False)
        sparse = sparse_generator(space, rates, origins, "a")
        np.testing.assert_allclose(sparse.toarray(), dense)

    def test_sparse_action_likelihood_matches_dense_on_partial_tips(self):
        space, tree = make_space(), make_tree()
        rates = model()
        tips = {"A": np.array([0.3, 1.0]), "B": np.array([1.0, 0.2])}
        dense_q = {c: generator(space, rates, {}, c, include_marks=False)[0]
                   for _, c in tree.edges()}
        p = {c: transition_matrix(q, tree.branch_length(c))
             for c, q in dense_q.items()}
        root = np.full(len(space.states), 1/len(space.states))
        expected = likelihood(tree, tips, p, root, posterior=False).log_likelihood
        sparse_q = {c: sparse_generator(space, rates, {}, c) for c in dense_q}
        actual = likelihood_only(tree, tips, sparse_q,
                                 {c: tree.branch_length(c) for c in dense_q}, root)
        self.assertAlmostEqual(actual, expected, places=11)

    def test_log_action_handles_matrix_right_sides(self):
        space = make_space()
        q = sparse_generator(space, model(), {}, "a")
        log_values = np.log(np.array([[0.2, 1.0], [1.0, 0.4]]))
        expected = np.log(expm(q.toarray()*0.25) @ np.exp(log_values))
        actual = log_action(q, 0.25, log_values)
        np.testing.assert_allclose(actual, expected, atol=1e-12)

    def test_all_unknown_tips_and_zero_rates_have_unit_likelihood(self):
        space, tree = make_space(), make_tree()
        rates = model(**{kind: 0. for kind in EDIT_KINDS})
        tips = {label: np.ones(len(space.states)) for label in ("A", "B")}
        result = evaluate_model(space, tree, tips, rates, posterior=False,
                                counts=False, backend="sparse")
        self.assertAlmostEqual(result["log_likelihood"], 0., places=12)

    def test_origin_mixture_and_root_weights_match_dense_backend(self):
        space, tree = make_material_space(n=1), make_tree()
        tips = {label: np.ones(len(space.states)) for label in ("A", "B")}
        rates = RateModel({kind: 0.03 for kind in EDIT_KINDS},
                          origin_root_weight=2.)
        dense = evaluate_model(space, tree, tips, rates, posterior=False,
                               counts=False)["log_likelihood"]
        sparse = evaluate_model(space, tree, tips, rates, posterior=False,
                                counts=False, backend="sparse")["log_likelihood"]
        self.assertAlmostEqual(sparse, dense, places=11)

    def test_origin_prior_mass_for_impossible_origin_is_not_renormalized(self):
        space, tree = make_material_space(n=1), make_tree(0.2, 0.3)
        tips = {"A": np.array([0., 1., 0.]), "B": np.array([1., 0., 1.])}
        rates = {kind: 0. for kind in EDIT_KINDS}
        rates["dna_insertion"] = 0.2
        model_with_weight = RateModel(rates, origin_root_weight=2.)
        value = evaluate_model(space, tree, tips, model_with_weight,
                               posterior=False, counts=False,
                               backend="sparse")["log_likelihood"]
        expected = np.log((1. - np.exp(-0.2*0.2)) / 4.)
        self.assertAlmostEqual(value, expected, places=11)

    def test_foreground_polytomy_matches_dense_backend(self):
        space = make_space()
        tree = SpeciesTree([
            {"node_id": "r", "parent_id": "", "label": "r"},
            {"node_id": "a", "parent_id": "r", "label": "A", "branch_length": .2},
            {"node_id": "b", "parent_id": "r", "label": "B", "branch_length": .3},
            {"node_id": "c", "parent_id": "r", "label": "C", "branch_length": .4},
        ])
        tips = {"A": np.array([1., 0.]), "B": np.array([0.2, 1.]),
                "C": np.array([1., 0.4])}
        rates = model(foreground=("a",), foreground_multiplier=3.)
        dense = evaluate_model(space, tree, tips, rates, posterior=False,
                               counts=False)["log_likelihood"]
        sparse = evaluate_model(space, tree, tips, rates, posterior=False,
                                counts=False, backend="sparse")["log_likelihood"]
        self.assertAlmostEqual(sparse, dense, places=11)

    def test_zero_branch_length_is_identity(self):
        space, tree = make_space(), make_tree(0., 0.3)
        tips = {"A": np.array([1., 0.]), "B": np.ones(2)}
        result = evaluate_model(space, tree, tips, model(), posterior=False,
                                counts=False, backend="sparse")
        self.assertTrue(np.isfinite(result["log_likelihood"]))

    def test_zero_generator_preserves_deep_finite_logs(self):
        space = make_space()
        tree = SpeciesTree([
            {"node_id": "r", "parent_id": "", "label": "r"},
            {"node_id": "x", "parent_id": "r", "label": "x", "branch_length": .1},
            {"node_id": "c", "parent_id": "r", "label": "C", "branch_length": .1},
            {"node_id": "a", "parent_id": "x", "label": "A", "branch_length": 0.},
            {"node_id": "b", "parent_id": "x", "label": "B", "branch_length": 0.},
        ])
        tiny = 1e-308
        tips = {"A": np.array([tiny, 1.]), "B": np.array([tiny, 1.]),
                "C": np.array([1., 0.])}
        zero = model(**{kind: 0. for kind in EDIT_KINDS})
        value = evaluate_model(space, tree, tips, zero, posterior=False,
                               counts=False, backend="sparse")["log_likelihood"]
        self.assertAlmostEqual(value, 2*np.log(tiny)-np.log(2.), places=10)

    def test_structurally_impossible_observations_return_negative_infinity(self):
        space, tree = make_space(), make_tree()
        zero = model(**{kind: 0. for kind in EDIT_KINDS})
        tips = {"A": np.array([1., 0.]), "B": np.array([0., 1.])}
        result = evaluate_model(space, tree, tips, zero, posterior=False,
                                counts=False, backend="sparse")
        self.assertEqual(result["log_likelihood"], -np.inf)

    def test_finite_child_support_underflow_is_an_arithmetic_error(self):
        space = make_space()
        tree = SpeciesTree([
            {"node_id": "r", "parent_id": "", "label": "r"},
            {"node_id": "x", "parent_id": "r", "label": "x", "branch_length": .1},
            {"node_id": "c", "parent_id": "r", "label": "C", "branch_length": .1},
            {"node_id": "a", "parent_id": "x", "label": "A", "branch_length": 0.},
            {"node_id": "b", "parent_id": "x", "label": "B", "branch_length": 0.},
        ])
        tiny = 1e-308
        tips = {"A": np.array([tiny, 1.]), "B": np.array([tiny, 1.]),
                "C": np.ones(2)}
        with self.assertRaises(ArithmeticError):
            evaluate_model(space, tree, tips, model(), posterior=False,
                           counts=False, backend="sparse")

    def test_final_negative_infinity_with_graph_feasible_support_is_numeric_error(self):
        tree = make_tree()
        q = csr_matrix([[-1., 1.], [0., 0.]])
        tips = {"A": np.array([0., 1.]), "B": np.ones(2)}
        with patch("intraphy.inference.configuration_sparse.expm_multiply",
                   return_value=np.array([0., 1.])):
            with self.assertRaises(ArithmeticError):
                likelihood_only(tree, tips, {"a": q, "b": q},
                                {"a": .2, "b": .3}, np.array([1., 0.]))

    def test_complete_state_space_larger_than_64_is_supported(self):
        space = make_material_space()
        self.assertGreater(len(space.states), 64)
        tree = make_tree()
        tips = {label: np.ones(len(space.states)) for label in ("A", "B")}
        value = evaluate_model(space, tree, tips, model(), posterior=False,
                               counts=False, backend="sparse")["log_likelihood"]
        self.assertAlmostEqual(value, 0., places=10)

    def test_sparse_backend_rejects_posterior_and_counts(self):
        space, tree = make_space(), make_tree()
        tips = {label: np.ones(len(space.states)) for label in ("A", "B")}
        with self.assertRaises(ValueError):
            evaluate_model(space, tree, tips, model(), posterior=False,
                           backend="sparse")
        with self.assertRaises(ValueError):
            evaluate_model(space, tree, tips, model(), counts=False,
                           backend="sparse")

    def test_sparse_generator_cache_counts_csr_storage_and_evicts(self):
        space = make_material_space(n=1)
        q = sparse_generator(space, model(), {"m0": "a"}, "a")
        self.assertGreater(q.nnz, 0)
        size = KernelCache._size(q)
        self.assertEqual(size, q.data.nbytes + q.indices.nbytes + q.indptr.nbytes)
        cache = KernelCache(maximum_bytes=size)
        first = cache.get_or_compute("first", lambda: q)
        self.assertIs(first, q)
        cache.get_or_compute("second", lambda: q.copy())
        self.assertEqual(len(cache.entries), 1)
        self.assertEqual(cache.bytes, size)

    def test_sparse_likelihood_is_unchanged_when_cache_cannot_retain_generators(self):
        space, tree = make_material_space(n=1), make_tree()
        tips = {label: np.ones(len(space.states)) for label in ("A", "B")}
        args = (space, tree, tips, model())
        cached = evaluate_model(*args, posterior=False, counts=False,
                                backend="sparse")["log_likelihood"]
        with patch("intraphy.inference.kernel_cache.KernelCache",
                   lambda: KernelCache(maximum_bytes=0)):
            evicted = evaluate_model(*args, posterior=False, counts=False,
                                     backend="sparse")["log_likelihood"]
        self.assertAlmostEqual(cached, evicted, places=12)
