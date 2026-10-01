"""Regression contracts for fit-local sparse generator template reuse."""
import unittest
from dataclasses import replace
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from intraphy.inference.configuration_fit_cache import CommonRateKernels
from intraphy.inference.configuration_model import RateModel, evaluate_model
from intraphy.inference.configuration_sparse import sparse_generator
from intraphy.inference.genomic_exon_rates import fit_family_rate
from intraphy.inference.kernel_cache import KernelCache
from intraphy.structure.edits import EDIT_KINDS
from intraphy.structure.space import enumerate_space
from intraphy.structure.types import Catalogue, ExonSpan, Material
from intraphy.topology import SpeciesTree


def make_tree():
    return SpeciesTree([
        {"node_id": "r", "parent_id": "", "label": "r"},
        {"node_id": "a", "parent_id": "r", "label": "A", "branch_length": 0.5},
        {"node_id": "b", "parent_id": "r", "label": "B", "branch_length": 0.5},
    ])


def make_space(with_material=True):
    spans = (ExonSpan(2, 4),)
    materials = (Material("m", 0, 1),) if with_material else ()
    return enumerate_space(Catalogue("f", "u", 4, spans, (), material=materials,
                                     boundary_candidates=spans))


def model(mu, **kwargs):
    return RateModel({kind: mu for kind in EDIT_KINDS}, **kwargs)


def informative_tips(state_space):
    present = np.array([float(bool(state.exons)) for state in state_space.states])
    return {"A": present, "B": np.ones(len(state_space.states))}


class GenomicExonCacheTests(unittest.TestCase):
    def test_templates_match_sparse_generator_across_rates_and_branches(self):
        state_space = make_space()
        origins = {"m": "a"}
        templates = CommonRateKernels(state_space, KernelCache())
        for mu in (0., 0.17):
            fit_model = model(mu, scale=1.7, foreground_multiplier=2.3,
                              foreground=frozenset({"a"}), origin_root_weight=2.)
            for child in ("a", "b"):
                actual = templates.generator(state_space, fit_model, origins, child)
                expected = sparse_generator(state_space, fit_model, origins, child)
                np.testing.assert_allclose(actual.toarray(), expected.toarray(), atol=1e-14)
                self.assertEqual(actual.nnz == 0, mu == 0.)

    def test_sparse_matches_uncached_sparse_and_dense_oracles(self):
        state_space, species_tree = make_space(), make_tree()
        tips = informative_tips(state_space)
        for mode in ("supplied", "unit"):
            for mu in (0., 0.13):
                fit_model = model(mu, origin_root_weight=2.)
                dense = evaluate_model(state_space, species_tree, tips, fit_model,
                    branch_length_mode=mode, posterior=False, counts=False)["log_likelihood"]
                plain = evaluate_model(state_space, species_tree, tips, fit_model,
                    branch_length_mode=mode, posterior=False, counts=False,
                    backend="sparse")["log_likelihood"]
                cached = evaluate_model(state_space, species_tree, tips, fit_model,
                    branch_length_mode=mode, posterior=False, counts=False,
                    backend="sparse", sparse_templates=CommonRateKernels(
                        state_space, KernelCache()))["log_likelihood"]
                self.assertAlmostEqual(plain, dense, places=11)
                self.assertAlmostEqual(cached, plain, places=11)

    def test_zero_rate_conflicting_tips_remain_impossible(self):
        state_space, species_tree = make_space(False), make_tree()
        empty = np.array([float(not state.exons) for state in state_space.states])
        present = np.array([float(bool(state.exons)) for state in state_space.states])
        conflicting = {"A": empty, "B": present}
        unknown = {"A": np.ones(len(state_space.states)),
                   "B": np.ones(len(state_space.states))}
        templates = CommonRateKernels(state_space, KernelCache())
        zero = model(0.)
        impossible = evaluate_model(state_space, species_tree, conflicting, zero,
            posterior=False, counts=False, backend="sparse",
            sparse_templates=templates)["log_likelihood"]
        baseline = evaluate_model(state_space, species_tree, unknown, zero,
            posterior=False, counts=False, backend="sparse",
            sparse_templates=templates)["log_likelihood"]
        self.assertEqual(impossible, -np.inf)
        self.assertAlmostEqual(baseline, 0., places=12)

    def test_builder_reuses_template_across_rates_without_mutating_it(self):
        state_space, origins = make_space(), {"m": "a"}
        cache = KernelCache()
        templates = CommonRateKernels(state_space, cache)
        with patch("intraphy.inference.configuration_fit_cache.sparse_generator",
                   wraps=sparse_generator) as builder:
            first = templates.generator(state_space, model(0.2), origins, "a")
            base = next(iter(cache.entries.values()))[0]
            before = base.copy()
            second = templates.generator(state_space, model(0.4), origins, "a")
        self.assertEqual(builder.call_count, 1)
        np.testing.assert_array_equal(base.toarray(), before.toarray())
        np.testing.assert_allclose(second.toarray(), 2. * first.toarray(), atol=1e-14)

    def test_cache_capacity_eviction_and_unit_namespaces(self):
        state_space, origins = make_space(), {"m": "a"}
        capacity_reference = sparse_generator(state_space,
            RateModel({kind: 1. for kind in EDIT_KINDS}), origins, "a")
        capacity = (capacity_reference.data.nbytes + capacity_reference.indices.nbytes
                    + capacity_reference.indptr.nbytes)
        cache = KernelCache(maximum_bytes=capacity)
        templates = CommonRateKernels(state_space, cache)
        templates.generator(state_space, model(0.1), origins, "a")
        evicted_q = templates.generator(state_space, model(0.1), origins, "b")
        expected_q = sparse_generator(state_space, model(0.1), origins, "b")
        np.testing.assert_allclose(evicted_q.toarray(), expected_q.toarray())
        templates.generator(state_space, model(0.1), origins, "a")
        self.assertEqual(cache.misses, 3)
        self.assertLessEqual(cache.bytes, cache.maximum_bytes)
        self.assertEqual(len(cache.entries), 1)

        shared = KernelCache()
        one = CommonRateKernels(state_space, shared)
        two = CommonRateKernels(state_space, shared)
        q1 = one.generator(state_space, model(0.1), origins, "a")
        q2 = two.generator(state_space, model(0.1), origins, "a")
        self.assertEqual(shared.misses, 2)
        np.testing.assert_allclose(q1.toarray(), q2.toarray())

    def test_distinct_spaces_do_not_contaminate_shared_cache(self):
        first_space, second_space = make_space(), make_space()
        cache = KernelCache()
        first = CommonRateKernels(first_space, cache)
        second = CommonRateKernels(second_space, cache)
        origins = {"m": "a"}
        q1 = first.generator(first_space, model(0.2), origins, "a")
        q2 = second.generator(second_space, model(0.2), origins, "a")
        self.assertEqual(cache.misses, 2)
        np.testing.assert_allclose(q1.toarray(), sparse_generator(
            first_space, model(0.2), origins, "a").toarray())
        np.testing.assert_allclose(q2.toarray(), sparse_generator(
            second_space, model(0.2), origins, "a").toarray())
        with self.assertRaisesRegex(ValueError, "different state space"):
            first.generator(second_space, model(0.2), origins, "a")

    def test_incomplete_space_and_unequal_rates_are_rejected(self):
        state_space = make_space()
        incomplete = replace(state_space, complete=False)
        templates = CommonRateKernels(incomplete, KernelCache())
        with self.assertRaisesRegex(ValueError, "state_space_incomplete"):
            templates.generator(incomplete, model(0.1), {"m": "a"}, "a")
        unequal = RateModel({kind: (0.1 if kind == EDIT_KINDS[0] else 0.2)
                             for kind in EDIT_KINDS})
        with self.assertRaisesRegex(ValueError, "equal edit rates"):
            CommonRateKernels(state_space, KernelCache()).generator(
                state_space, unequal, {"m": "a"}, "a")

    def test_capacity_zero_preserves_likelihood(self):
        state_space, species_tree = make_space(), make_tree()
        tips, fit_model = informative_tips(state_space), model(0.13, origin_root_weight=2.)
        expected = evaluate_model(state_space, species_tree, tips, fit_model,
            posterior=False, counts=False, backend="sparse")["log_likelihood"]
        actual = evaluate_model(state_space, species_tree, tips, fit_model,
            posterior=False, counts=False, backend="sparse",
            sparse_templates=CommonRateKernels(
                state_space, KernelCache(maximum_bytes=0)))["log_likelihood"]
        self.assertAlmostEqual(actual, expected, places=12)

    def test_real_shared_rate_fit_matches_uncached_evaluator(self):
        exon = ExonSpan(0, 2)
        state_space = enumerate_space(Catalogue("f", "fit", 2, (exon,), (),
            boundary_candidates=(exon,)))
        species_tree = make_tree()
        present = np.array([float(bool(state.exons)) for state in state_space.states])
        absent = 1. - present
        same, different = {"A": present, "B": present}, {"A": present, "B": absent}
        units = [{"space": state_space, "tree": species_tree,
                  "tips": same if i < 3 else different} for i in range(4)]
        cached = fit_family_rate(units)
        original = evaluate_model
        def without_templates(*args, **kwargs):
            kwargs.pop("sparse_templates", None)
            return original(*args, **kwargs)
        with patch("intraphy.inference.genomic_exon_rates.evaluate_model",
                   side_effect=without_templates):
            baseline = fit_family_rate(units)
        self.assertEqual(cached["status"], "estimated_conditional_composite_rate")
        self.assertEqual(cached["status"], baseline["status"])
        self.assertTrue(cached["converged"])
        self.assertEqual(cached["upper_tail_diagnostic"]["reason"],
                         baseline["upper_tail_diagnostic"]["reason"])
        self.assertAlmostEqual(cached["mu"], baseline["mu"], places=8)
        self.assertAlmostEqual(cached["log_likelihood"], baseline["log_likelihood"], places=10)

    def test_exact_rate_memo_and_shared_fit_cache_on_analytic_mock(self):
        calls, helpers = [], []
        def analytic(space_arg, tree_arg, tips_arg, fit_model, **kwargs):
            helpers.append(kwargs["sparse_templates"])
            x = fit_model.rates[EDIT_KINDS[0]] * tree_arg.branch_length("a")
            calls.append((id(space_arg), x))
            return {"log_likelihood": -(x - 2.) ** 2 if tips_arg["pattern"] == "finite" else 0.}
        def repeated_optimizer(fun, point, **kwargs):
            fun(point)
            fun(point)
            return SimpleNamespace(x=point, success=True, message="ok", status=0, nit=1)
        unit_tree = SimpleNamespace(edges=lambda: iter((("r", "a"),)),
                                    branch_length=lambda child: 1.)
        units = [{"space": object(), "tree": unit_tree,
                  "tips": {"pattern": "finite" if i == 0 else "constant"}}
                 for i in range(3)]
        with (
            patch("intraphy.inference.genomic_exon_rates.evaluate_model", side_effect=analytic),
            patch("intraphy.inference.genomic_exon_rates.minimize",
                  side_effect=repeated_optimizer),
        ):
            fit_family_rate(units)
        self.assertEqual(len(calls), len(set(calls)))
        self.assertEqual(len({id(helper._cache) for helper in helpers}), 1)

    def test_large_complete_space_matches_uncached_and_origin_limit_still_fails(self):
        materials = tuple(Material(f"m{i}", 2*i, 2*i+1) for i in range(4))
        state_space = enumerate_space(Catalogue("f", "large", 8, (), (),
            material=materials, boundary_candidates=()))
        self.assertTrue(state_space.complete)
        self.assertGreater(len(state_space.states), 64)
        present = np.array([float(all(value == 1 for value in state.material))
                            for state in state_space.states])
        tips = {"A": present, "B": np.ones(len(state_space.states))}
        fit_model = model(0.08, origin_root_weight=2.)
        plain = evaluate_model(state_space, make_tree(), tips, fit_model,
            posterior=False, counts=False, backend="sparse")["log_likelihood"]
        cached = evaluate_model(state_space, make_tree(), tips, fit_model,
            posterior=False, counts=False, backend="sparse",
            sparse_templates=CommonRateKernels(
                state_space, KernelCache()))["log_likelihood"]
        self.assertAlmostEqual(cached, plain, places=11)
        for templates in (None, CommonRateKernels(state_space, KernelCache())):
            with self.assertRaisesRegex(ValueError, "origin_scenarios_incomplete"):
                evaluate_model(state_space, make_tree(), tips, fit_model, max_origins=1,
                    posterior=False, counts=False, backend="sparse",
                    sparse_templates=templates)
