from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

from intraphy.inference.locus_comparison import compare_locus_rates
from intraphy.inference.locus_profiles import profile_locus_rate
from intraphy.inference.locus_rates import LocusFitUnit, fit_locus_rates
from intraphy.structure.locus_types import LocusState, ProcessEdge
from intraphy.topology import SpeciesTree


def _unit(tips=None):
    catalogue = SimpleNamespace(
        opportunities=(SimpleNamespace(id="del", rate_group="deletion"),))
    process = SimpleNamespace(
        catalogue=catalogue,
        states=(LocusState((2,), frozenset()), LocusState((1,), frozenset())),
        edges=(ProcessEdge(1, 0, "del", "delete", "dna_deletion", "deletion", 1.0),),
    )
    tree = SpeciesTree([
        {"node_id": "root", "parent_id": "", "label": "root"},
        {"node_id": "a", "parent_id": "root", "label": "a", "branch_length": 1.0},
        {"node_id": "b", "parent_id": "root", "label": "b", "branch_length": 1.0},
    ])
    return LocusFitUnit(process, tree, tips or {"a": [0.0, 1.0], "b": [0.0, 1.0]}, [0.0, 1.0])


class LocusStatisticsTests(unittest.TestCase):
    def test_multistart_records_attempts_and_marks_skipped_curvature(self):
        result = fit_locus_rates(
            [_unit()], {"deletion": 0.4}, start_scales=(0.2, 1.0, 5.0),
            compute_curvature=False, maxiter=100,
        )
        self.assertTrue(result.optimizer_success)
        self.assertEqual(len(result.start_results), 3)
        self.assertFalse(result.curvature_evaluated)
        self.assertEqual(result.curvature_eigenvalues, ())
        self.assertFalse(result.curvature_positive_definite)
        self.assertFalse(result.unresolved_higher_likelihood)

    def test_profile_retains_requested_boundary_and_conditional_points(self):
        result = profile_locus_rate(
            (_unit() for _ in range(1)), {"deletion": 0.4}, {}, "deletion",
            (value for value in (0.0, 0.5)), start_scales=(0.2, 1.0), maxiter=100,
        )
        self.assertEqual(tuple(point.value for point in result.points), (0.0, 0.5))
        self.assertEqual(len(result.points), 2)
        self.assertIsNotNone(result.points[0].fit)
        self.assertIn(result.points[0].status, {"evaluated", "profile_fit_unsuccessful"})

    def test_deletion_mle_and_profile_match_two_tip_analytic_likelihood(self):
        tips = {"a": [0.0, 1.0], "b": [1.0, 0.0]}
        unit = _unit(tips)
        fitted = fit_locus_rates([unit], {"deletion": 0.4}, maxiter=500)
        self.assertTrue(fitted.optimizer_success)
        self.assertAlmostEqual(fitted.rates["deletion"], np.log(2.0), delta=2e-3)

        profile = profile_locus_rate([unit], {"deletion": 0.4}, {}, "deletion",
                                     (0.0, 0.4, np.log(2.0)), maxiter=500)
        for point in profile.points[1:]:
            expected = -point.value + np.log1p(-np.exp(-point.value))
            self.assertAlmostEqual(point.fit.fitted_log_likelihood, expected, places=8)
        self.assertFalse(profile.points[0].fit.optimizer_success)
        self.assertEqual(profile.points[0].status, "zero_likelihood")
        self.assertEqual(profile.points[0].log_likelihood_difference, -np.inf)
        self.assertEqual(profile.status, "evaluated")

    def test_fixed_zero_nested_comparison_uses_same_inputs_and_has_no_default_pvalue(self):
        result = compare_locus_rates(
            (_unit() for _ in range(1)), {"deletion": 0.4}, {}, (name for name in ("deletion",)),
            start_scales=(0.2, 1.0), maxiter=100, compute_curvature=False,
        )
        self.assertTrue(result.null_fit.optimizer_success)
        self.assertTrue(result.full_fit.optimizer_success)
        self.assertIsNotNone(result.statistic)
        self.assertGreaterEqual(result.statistic, 0.0)
        self.assertLessEqual(result.raw_log_likelihood_difference, result.likelihood_tolerance)
        self.assertFalse(result.full_fit.curvature_evaluated)
        self.assertFalse(hasattr(result, "pvalue"))

    def test_boolean_profile_values_are_rejected(self):
        with self.assertRaises(ValueError):
            profile_locus_rate([_unit()], {"deletion": 0.4}, {}, "deletion", [True])

    def test_higher_likelihood_failed_start_withholds_comparison_statistic(self):
        from intraphy.inference.locus_rates import LocusFitResult

        unit = _unit()
        mocked_optimizers = (
            SimpleNamespace(x=np.array([0.5]), success=True, status=0, message="converged", nit=1),
            SimpleNamespace(x=np.array([0.0]), success=False, status=1, message="iteration limit", nit=1),
        )
        with patch("intraphy.inference.locus_rates.minimize", side_effect=mocked_optimizers):
            unresolved_full = fit_locus_rates([unit], {"deletion": 0.4},
                                              start_scales=(1.0, 5.0), compute_curvature=False)
        self.assertTrue(unresolved_full.optimizer_success)
        self.assertTrue(unresolved_full.unresolved_higher_likelihood)
        self.assertGreater(unresolved_full.start_results[1]["fitted_log_likelihood"],
                           unresolved_full.fitted_log_likelihood)
        null_fit = LocusFitResult(
            0.0, 0.0, {"deletion": 0.0}, True, 0, "all fixed", 0, (), (), 0, 0,
            False, False, False, 1, (), 0, False, False,
        )
        with patch("intraphy.inference.locus_comparison.fit_locus_rates",
                   side_effect=(null_fit, unresolved_full)):
            result = compare_locus_rates([unit], {"deletion": 0.4}, {}, ("deletion",))
        self.assertIsNone(result.statistic)
        self.assertEqual(result.status, "one_or_both_fits_unsuccessful")
        self.assertAlmostEqual(result.raw_log_likelihood_difference, -0.4)

    def test_comparison_clamps_tiny_two_sided_likelihood_difference(self):
        from intraphy.inference.locus_rates import LocusFitResult

        def fit_result(ll):
            return LocusFitResult(
                ll, ll, {"deletion": 0.4}, True, 0, "mock", 1, (), (), 0, 1,
                False, False, False, 1, (), 0, False, False,
            )

        null_fit = fit_result(-1.0)
        full_fit = fit_result(-1.0 + 1e-12)
        with patch("intraphy.inference.locus_comparison.fit_locus_rates",
                   side_effect=(null_fit, full_fit)):
            result = compare_locus_rates([_unit()], {"deletion": 0.4}, {}, ("deletion",))
        self.assertAlmostEqual(result.raw_log_likelihood_difference, 1e-12)
        self.assertEqual(result.statistic, 0.0)


if __name__ == "__main__":
    unittest.main()
