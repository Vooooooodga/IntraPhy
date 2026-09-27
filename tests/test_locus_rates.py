from types import SimpleNamespace
import unittest

import numpy as np

from intraphy.inference.locus_likelihood import evaluate_locus
from intraphy.inference.locus_rates import (
    LocusFitUnit,
    _local_curvature,
    fit_locus_rates,
)
from intraphy.structure.locus_types import LocusState, ProcessEdge
from intraphy.topology import SpeciesTree


def _tree():
    return SpeciesTree([
        {"node_id": "root", "parent_id": "", "label": "root"},
        {"node_id": "a", "parent_id": "root", "label": "a", "branch_length": 1.0},
        {"node_id": "b", "parent_id": "root", "label": "b", "branch_length": 1.0},
    ])


def _deletion_process():
    catalogue = SimpleNamespace(opportunities=(SimpleNamespace(id="del", rate_group="deletion"),))
    edge = ProcessEdge(1, 0, "del", "del_outcome", "dna_deletion", "deletion", 1.0)
    states = (LocusState((2,), frozenset()), LocusState((1,), frozenset()))
    return SimpleNamespace(catalogue=catalogue, states=states, edges=(edge,))


class LocusRateTests(unittest.TestCase):
    def test_analytic_mixed_hessian_stencils_cover_all_boundary_patterns(self):
        def quadratic(point):
            if np.any(point < 0):
                raise AssertionError("curvature stencil evaluated a negative rate")
            first, second = point
            return first * first + 3.0 * first * second + 2.0 * second * second

        expected = np.linalg.eigvalsh(np.array([[2.0, 3.0], [3.0, 4.0]]))
        cases = (
            (np.array([0.0, 2.0]), np.array([True, False])),
            (np.array([2.0, 0.0]), np.array([False, True])),
            (np.array([0.0, 0.0]), np.array([True, True])),
            (np.array([2.0, 3.0]), np.array([False, False])),
            (np.array([1e-6, 2.0]), np.array([False, False])),
        )
        for point, boundary in cases:
            with self.subTest(point=point, boundary=boundary):
                eigenvalues, rank, dimension, positive, negative, flat = _local_curvature(
                    quadratic, point, boundary)
                np.testing.assert_allclose(eigenvalues, expected, rtol=1e-7, atol=1e-7)
                self.assertEqual(rank, dimension)
                self.assertEqual(dimension, 2)
                self.assertFalse(positive)
                self.assertTrue(negative)
                self.assertFalse(flat)

    def test_declared_but_unreachable_free_rate_is_reported_flat(self):
        catalogue = SimpleNamespace(opportunities=(SimpleNamespace(id="unused", rate_group="unused"),))
        states = (LocusState((0,), frozenset()), LocusState((1,), frozenset()))
        process = SimpleNamespace(catalogue=catalogue, states=states, edges=())
        tree = SpeciesTree([{"node_id": "tip", "parent_id": "", "label": "tip"}])
        unit = LocusFitUnit(process, tree, {"tip": [1.0, 1.0]}, [0.25, 0.75])
        result = fit_locus_rates([unit], {"unused": 0.7}, maxiter=100)
        self.assertTrue(result.optimizer_success)
        self.assertAlmostEqual(result.rates["unused"], 0.7)
        self.assertEqual(result.curvature_rank, 0)
        self.assertTrue(result.locally_flat)

    def test_serial_and_independent_unit_worker_fits_agree(self):
        process = _deletion_process()
        tree = _tree()
        tips = {"a": [0, 1], "b": [1, 0]}
        unit = LocusFitUnit(process, tree, tips, [0, 1])
        serial = fit_locus_rates([unit, unit], {"deletion": 0.4}, workers=1, maxiter=500)
        parallel = fit_locus_rates([unit, unit], {"deletion": 0.4}, workers=2, maxiter=500)
        self.assertTrue(serial.optimizer_success and parallel.optimizer_success)
        self.assertEqual(parallel.workers, 2)
        self.assertAlmostEqual(parallel.fitted_log_likelihood, serial.fitted_log_likelihood, delta=1e-9)
        self.assertAlmostEqual(parallel.rates["deletion"], serial.rates["deletion"], delta=1e-7)

    def test_tip_emission_vector_is_a_probability_not_an_arbitrary_scale(self):
        process = _deletion_process()
        tree = SpeciesTree([{"node_id": "tip", "parent_id": "", "label": "tip"}])
        with self.assertRaisesRegex(ValueError, "must not exceed one"):
            evaluate_locus(process, tree, {"tip": [1.2, 0.2]}, [0, 1], {"deletion": 0.2})
