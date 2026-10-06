"""Analytic contracts for log-domain sparse CTMC action fallback."""
import math
import unittest

import numpy as np
from scipy.sparse import csr_matrix

from intraphy.inference.configuration_sparse import log_action


def directed_chain(hops):
    """Unit-rate 0 -> 1 -> ... -> hops chain, with an absorbing endpoint."""
    size = hops + 1
    rows = list(range(hops))
    columns = [index + 1 for index in range(hops)]
    data = [1.] * hops
    rows.extend(range(hops))
    columns.extend(range(hops))
    data.extend([-1.] * hops)
    return csr_matrix((data, (rows, columns)), shape=(size, size))


def basis_log(size, index):
    values = np.full(size, -np.inf)
    values[index] = 0.
    return values


class UniformizationTests(unittest.TestCase):
    def test_two_hop_recovery_for_generator_and_transpose(self):
        q = directed_chain(2)
        t = 1e-18
        expected_two_hop = 2. * math.log(t) - math.log(2.)

        row_action = log_action(q, t, basis_log(3, 2))
        self.assertTrue(np.isfinite(row_action).all())
        self.assertAlmostEqual(row_action[0], expected_two_hop, delta=1e-12)

        transpose = q.transpose().tocsr()
        forward_action = log_action(transpose, t, basis_log(3, 0))
        self.assertTrue(np.isfinite(forward_action).all())
        self.assertAlmostEqual(forward_action[2], expected_two_hop, delta=1e-12)

        unreachable_action = log_action(transpose, t, basis_log(3, 2))
        np.testing.assert_array_equal(unreachable_action, basis_log(3, 2))

    def test_below_float_range_eight_hop_probability_has_finite_log(self):
        hops = 8
        t = 1e-50
        result = log_action(directed_chain(hops), t, basis_log(hops + 1, hops))
        self.assertTrue(np.isfinite(result).all())
        expected = hops * math.log(t) - math.lgamma(hops + 1)
        self.assertAlmostEqual(result[0], expected, delta=1e-10)

    def test_matrix_rhs_keeps_column_support_independent(self):
        q = directed_chain(2)
        t = 1e-18
        values = np.column_stack((basis_log(3, 2), basis_log(3, 0)))
        result = log_action(q, t, values)
        first = log_action(q, t, values[:, 0])
        second = log_action(q, t, values[:, 1])
        np.testing.assert_allclose(result[:, 0], first)
        np.testing.assert_array_equal(result[:, 1], second)
        self.assertTrue(np.isfinite(result[:, 0]).all())
        self.assertTrue(np.isfinite(result[0, 1]))
        self.assertTrue(np.isneginf(result[1:, 1]).all())

    def test_finite_input_far_below_float_range_is_preserved(self):
        # State 0 cannot receive mass from the high-valued isolated state 2.
        q = csr_matrix([[-1., 1., 0.], [0., 0., 0.], [0., 0., 0.]])
        values = np.array([-1000., -np.inf, 0.])
        result = log_action(q, 1., values)
        self.assertAlmostEqual(result[0], -1001., delta=1e-12)
        self.assertEqual(result[1], -np.inf)
        self.assertAlmostEqual(result[2], 0., delta=1e-12)

    def test_zero_length_and_zero_generator_are_identity(self):
        q = directed_chain(2)
        values = np.array([-900., -np.inf, 0.])
        np.testing.assert_array_equal(log_action(q, 0., values), values)
        np.testing.assert_array_equal(log_action(csr_matrix((3, 3)), 1., values), values)


if __name__ == "__main__":
    unittest.main()
