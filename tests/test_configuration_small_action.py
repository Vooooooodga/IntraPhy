"""Contracts for the small-state dense CTMC action backend."""
import math
import unittest
from unittest.mock import patch

import numpy as np
from scipy.linalg import expm
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply

from intraphy.inference.configuration_sparse import _support_after_edge, log_action
from intraphy.inference.configuration_uniformization import log_uniformization_action


def directed_chain(size, rate=1.):
    rows = list(range(size - 1))
    columns = list(range(1, size))
    data = [rate] * (size - 1)
    rows.extend(range(size - 1))
    columns.extend(range(size - 1))
    data.extend([-rate] * (size - 1))
    return csr_matrix((data, (rows, columns)), shape=(size, size))


class SmallActionTests(unittest.TestCase):
    def test_fully_positive_results_skip_support_traversal_for_both_backends(self):
        for size in (2, 17):
            q = directed_chain(size)
            vector_values = np.linspace(-1., 0., size)
            matrix_values = np.column_stack((vector_values, vector_values[::-1]))
            for values in (vector_values, matrix_values):
                with self.subTest(size=size, rhs_ndim=values.ndim):
                    with patch(
                            "intraphy.inference.configuration_sparse._support_after_edge",
                            side_effect=AssertionError("positive result traversed support")):
                        actual = log_action(q, 0.2, values)
                    expected = np.log(expm(q.toarray() * 0.2) @ np.exp(values))
                    np.testing.assert_allclose(actual, expected,
                                               rtol=1e-11, atol=1e-12)

    def test_zero_result_entries_keep_support_traversal_per_column(self):
        for size in (2, 17):
            q = directed_chain(size)
            vector_values = np.full(size, -np.inf)
            vector_values[0] = 0.
            matrix_values = np.column_stack((vector_values, np.zeros(size)))
            for values in (vector_values, matrix_values):
                with self.subTest(size=size, rhs_ndim=values.ndim):
                    with patch(
                            "intraphy.inference.configuration_sparse._support_after_edge",
                            wraps=_support_after_edge) as support:
                        actual = log_action(q, 0.2, values)
                    self.assertEqual(support.call_count, 1)
                    self.assertEqual(actual[-1] if values.ndim == 1
                                     else actual[-1, 0], -np.inf)

    def test_dense_action_matches_matrix_exponential_for_small_sizes_and_orientations(self):
        for size in (2, 3, 16):
            q = directed_chain(size)
            subgenerator = q.copy()
            subgenerator.setdiag(subgenerator.diagonal() - 0.17)
            values = np.linspace(-1.2, -0.1, size)
            matrix_values = np.column_stack((values, values[::-1]))
            for action in (q, q.transpose().tocsr(), subgenerator):
                with self.subTest(size=size, action=action is q):
                    with patch(
                            "intraphy.inference.configuration_sparse.expm_multiply",
                            side_effect=AssertionError("small action used sparse backend")):
                        actual = log_action(action, 0.23, matrix_values)
                    expected = np.log(expm(action.toarray() * 0.23)
                                      @ np.exp(matrix_values))
                    np.testing.assert_allclose(actual, expected, rtol=1e-12,
                                               atol=1e-12)

    def test_dense_action_matches_previous_sparse_action(self):
        for size in (2, 3, 16):
            q = directed_chain(size, rate=0.8)
            subgenerator = q.copy()
            subgenerator.setdiag(subgenerator.diagonal() - 0.13)
            vector_values = np.linspace(-0.8, -0.1, size)
            matrix_values = np.column_stack((vector_values, vector_values[::-1]))
            for action in (q, q.transpose().tocsr(), subgenerator):
                for values in (vector_values, matrix_values):
                    with self.subTest(size=size, rhs_ndim=values.ndim):
                        expected = expm_multiply(action * 0.37, np.exp(values))
                        with np.errstate(divide="ignore"):
                            expected = np.log(expected)
                        actual = log_action(action, 0.37, values)
                        np.testing.assert_allclose(actual, expected,
                                                   rtol=2e-11, atol=2e-12)

    def test_large_action_keeps_sparse_exponential_action(self):
        q = directed_chain(17)
        values = np.linspace(-0.5, 0., 17)
        with patch(
                "intraphy.inference.configuration_sparse.expm_multiply",
                wraps=expm_multiply) as sparse_action:
            actual = log_action(q, 0.2, values)
        expected = np.log(expm(q.toarray() * 0.2) @ np.exp(values))
        np.testing.assert_allclose(actual, expected, rtol=1e-11, atol=1e-12)
        sparse_action.assert_called_once()

    def test_two_state_analytic_transition_and_zero_rate_boundary(self):
        forward, reverse, length = 0.7, 0.2, 0.4
        q = csr_matrix([[-forward, forward], [reverse, -reverse]])
        values = np.array([-np.inf, 0.])
        actual = log_action(q, length, values)
        total = forward + reverse
        p01 = forward / total * (1. - math.exp(-total * length))
        p11 = forward / total + reverse / total * math.exp(-total * length)
        np.testing.assert_allclose(actual, np.log([p01, p11]), atol=1e-13)
        np.testing.assert_array_equal(log_action(q, 0., values), values)
        stationary = log_action(q, 100. / total, values)
        np.testing.assert_allclose(stationary,
                                   np.log([forward / total, forward / total]),
                                   atol=1e-13)

    def test_tiny_reachable_probability_falls_back_to_log_uniformization(self):
        q = directed_chain(3)
        length = 1e-200
        values = np.array([-np.inf, -np.inf, 0.])
        with patch(
                "intraphy.inference.configuration_sparse.log_uniformization_action",
                wraps=log_uniformization_action) as fallback:
            actual = log_action(q, length, values)
        fallback.assert_called_once()
        expected = 2. * math.log(length) - math.log(2.)
        self.assertAlmostEqual(actual[0], expected, delta=1e-10)
        self.assertTrue(np.isfinite(actual).all())

    def test_lost_input_scale_uses_log_domain_fallback(self):
        q = csr_matrix([[-1., 0.], [0., 0.]])
        actual = log_action(q, 1., np.array([-1000., 0.]))
        np.testing.assert_allclose(actual, [-1001., 0.], atol=1e-12)

    def test_reducible_generator_keeps_unreachable_states_at_negative_infinity(self):
        q = csr_matrix([[-1., 1., 0.], [0., 0., 0.], [0., 0., 0.]])
        values = np.array([-np.inf, -np.inf, 0.])
        actual = log_action(q, 0.4, values)
        np.testing.assert_array_equal(actual, values)

    def test_zero_generator_is_identity_for_vector_and_matrix_right_sides(self):
        q = csr_matrix((3, 3))
        values = np.array([[-900., -np.inf], [-np.inf, -2.], [0., -3.]])
        np.testing.assert_array_equal(log_action(q, 1.7, values), values)

    def test_high_norm_small_action_matches_dense_reference(self):
        q = directed_chain(2, rate=1e5)
        values = np.log(np.array([0.3, 0.8]))
        length = 0.02
        actual = log_action(q, length, values)
        expected = np.log(expm(q.toarray() * length) @ np.exp(values))
        np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
