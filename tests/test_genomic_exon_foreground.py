"""Numerical and roster contracts for the native foreground comparison."""
from dataclasses import replace
import math
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from intraphy.inference.configuration_ctmc import likelihood, transition_matrix
from intraphy.inference.configuration_model import RateModel, generator
from intraphy.inference.genomic_exon_comparison import compare_foreground
from intraphy.inference.genomic_exon_rates import FamilyRateLikelihood, fit_family_rate
from intraphy.inference.genomic_exon_calibration import _catalogues, _unit
from intraphy.structure.edits import EDIT_KINDS
from intraphy.structure.origins import origin_scenarios
from intraphy.structure.space import enumerate_space
from intraphy.structure.types import ExonConfiguration, ExonSpan, Material
from intraphy.topology import SpeciesTree


def _family_unit(catalogue):
    space = enumerate_space(catalogue)
    source = _unit(catalogue, space, catalogue.unit)
    return {"unit_id": catalogue.unit, "space": space, "tree": source.tree,
            "tips": source.tips}


def _foreground_origin_unit():
    catalogue = _catalogues("shared-deletion")[0]
    unit = _family_unit(catalogue)
    material_index = 0
    tips = {}
    for species, values in unit["tips"].items():
        allowed = {1} if species == "A" else {0, 2}
        if species == "A":
            supported = np.flatnonzero(values > 0)
            observed_exons = unit["space"].states[int(supported[0])].exons
        else:
            observed_exons = ()
        tips[species] = np.asarray([
            float(state.exons == observed_exons
                  and state.material[material_index] in allowed)
            for state in unit["space"].states
        ])
        assert np.any(tips[species]), f"Synthetic tip {species} has no allowed state"
    unit["tips"] = tips
    return catalogue, unit


def _shift_catalogue(catalogue, offset):
    def shift(span):
        return ExonSpan(span.start + offset, span.end + offset)

    observations = tuple(replace(observation,
        configurations=tuple(ExonConfiguration(
            tuple(shift(span) for span in config.exons), config.material)
            for config in observation.configurations),
        unknown_intervals=tuple(shift(span) for span in observation.unknown_intervals),
        alternative_exons=tuple(shift(span) for span in observation.alternative_exons))
        for observation in catalogue.observations)
    return replace(catalogue, unit=f"{catalogue.unit}_shifted",
        length=catalogue.length + offset,
        spans=tuple(shift(span) for span in catalogue.spans),
        junctions=tuple((donor + offset, acceptor + offset)
                        for donor, acceptor in catalogue.junctions),
        material=tuple(Material(item.id, item.start + offset, item.end + offset)
                       for item in catalogue.material),
        observations=observations,
        boundary_candidates=tuple(shift(span) for span in catalogue.boundary_candidates))


def _dense_origin_log_likelihood(unit, mu, foreground, multiplier, *, root_weight=1.):
    """Literal origin-scenario sum, retaining each scenario's prior mass."""
    model = RateModel({kind: mu for kind in EDIT_KINDS},
                      foreground=frozenset(foreground),
                      foreground_multiplier=multiplier,
                      origin_root_weight=root_weight)
    tree, space, tips = unit["tree"], unit["space"], unit["tips"]
    total = -math.inf
    for origins, root, log_prior in origin_scenarios(
            space, tree, tips=tips, root_weight=root_weight):
        matrices = {}
        for parent, child in tree.edges():
            q, _ = generator(space, model, origins, child, include_marks=False)
            matrices[child] = transition_matrix(q, tree.branch_length(child))
        result = likelihood(tree, tips, matrices, root.astype(float) / root.sum(), False)
        total = float(np.logaddexp(total, result.log_likelihood + log_prior))
    return total


class GenomicExonForegroundTests(unittest.TestCase):
    def test_rho_one_is_nested_null_and_matches_literal_dense_origin_sum(self):
        catalogue = _catalogues("shared-deletion")[0]
        unit = _family_unit(catalogue)
        family_units = {catalogue.family: [unit]}
        result = compare_foreground(family_units, {"A"}, profile_multipliers=(1.,))
        profile_one = next(row for row in result["profile"]
                           if row["foreground_multiplier"] == 1.)
        if result["null"]["log_likelihood"] is None:
            self.assertIsNone(profile_one["log_likelihood"])
        else:
            self.assertAlmostEqual(result["null"]["log_likelihood"],
                                   profile_one["log_likelihood"], delta=1e-10)
        fixed_mu = .2
        workspace = FamilyRateLikelihood([unit])
        fixed_ll = workspace.log_likelihood(
            fixed_mu * workspace.exposure, frozenset({"A"}), 1.)
        self.assertAlmostEqual(fixed_ll,
                               _dense_origin_log_likelihood(unit, fixed_mu, {"A"}, 1.),
                               delta=1e-9)
        self.assertIsNone(result["p_value"])
        self.assertEqual(result["schema"], "intraphy.exon-foreground-comparison/1")

    def test_zero_multiplier_preserves_origin_prior_denominator(self):
        _catalogue, unit = _foreground_origin_unit()
        fixed_mu = .2
        workspace = FamilyRateLikelihood([unit])
        fixed_ll = workspace.log_likelihood(
            fixed_mu * workspace.exposure, frozenset({"A"}), 0.)
        expected = _dense_origin_log_likelihood(unit, fixed_mu, {"A"}, 0.)
        self.assertAlmostEqual(fixed_ll, expected, delta=1e-9)

        # Here an introduction on foreground terminal A is admissible from the
        # tip evidence but impossible when its branch multiplier is zero.
        model = RateModel({kind: fixed_mu for kind in EDIT_KINDS},
                          foreground=frozenset({"A"}), foreground_multiplier=0.)
        terms = []
        for origins, root, log_prior in origin_scenarios(
                unit["space"], unit["tree"], tips=unit["tips"], root_weight=1.):
            matrices = {}
            for parent, child in unit["tree"].edges():
                q, _ = generator(unit["space"], model, origins, child,
                                 include_marks=False)
                matrices[child] = transition_matrix(
                    q, unit["tree"].branch_length(child))
            likelihood_value = likelihood(unit["tree"], unit["tips"], matrices,
                                          root.astype(float) / root.sum(), False).log_likelihood
            terms.append((origins[unit["space"].catalogue.material[0].id] == "A",
                          likelihood_value, log_prior))
        original_mass = float(np.logaddexp.reduce(
            np.asarray([ll + prior for _is_a_origin, ll, prior in terms])))
        finite = [(ll, prior) for _is_a_origin, ll, prior in terms if np.isfinite(ll)]
        surviving_prior_mass = sum(math.exp(prior) for _, prior in finite)
        renormalized = original_mass - math.log(surviving_prior_mass)
        self.assertTrue(any(is_a_origin for is_a_origin, _ll, _prior in terms))
        self.assertTrue(any(is_a_origin and not np.isfinite(ll)
                            for is_a_origin, ll, _prior in terms))
        self.assertAlmostEqual(original_mass, fixed_ll, delta=1e-9)
        self.assertNotAlmostEqual(original_mass, renormalized, delta=1e-6)

    def test_family_likelihood_aggregation_distinguishes_finite_and_nonfinite_values(self):
        unit = _family_unit(_catalogues("geometry")[0])
        for values, expected in (((-.25, -.5), -.75), ((-math.inf, -.5), -math.inf)):
            with self.subTest(values=values), patch(
                    "intraphy.inference.genomic_exon_rates.evaluate_model",
                    side_effect=[{"log_likelihood": value} for value in values]):
                workspace = FamilyRateLikelihood([unit, unit])
                observed = workspace.log_likelihood(1.)
            if expected == -math.inf:
                self.assertEqual(observed, expected)
            else:
                self.assertAlmostEqual(observed, expected, delta=1e-15)

        for values in ((math.nan, -.5), (math.inf, -.5), (1e308, 1e308)):
            with self.subTest(values=values), patch(
                    "intraphy.inference.genomic_exon_rates.evaluate_model",
                    side_effect=[{"log_likelihood": value} for value in values]):
                workspace = FamilyRateLikelihood([unit, unit])
                with self.assertRaises(ArithmeticError):
                    workspace.log_likelihood(1.)

    def test_each_gene_keeps_its_own_background_rate_and_all_profile_rosters(self):
        first = _catalogues("geometry")[0]
        second = replace(first, family="gene_b", unit="unit_b")
        first_unit, second_unit = _family_unit(first), _family_unit(second)
        # Preserve the same tree/tip roster but use a different pattern for gene B.
        states = first_unit["space"].states
        present = np.asarray([float(bool(state.exons)) for state in states])
        second_unit["tips"] = {"A": present, "B": 1. - present,
                               "C": present, "D": 1. - present}
        family_units = {"gene_a": [first_unit], "gene_b": [second_unit],
                        "requested_but_excluded": []}
        result = compare_foreground(family_units, {"A"}, profile_multipliers=(.5, 2.))
        self.assertEqual(result["counts"]["requested_families"], 3)
        self.assertEqual(result["counts"]["eligible_families"], 2)
        self.assertEqual(result["counts"]["excluded_families"],
                         ["requested_but_excluded"])
        expected_names = set(family_units)
        for point in [result["null"], result["alternative"], *result["profile"]]:
            self.assertEqual(set(point["families"]), expected_names)
            self.assertEqual(point["families"]["requested_but_excluded"]["status"],
                             "no_estimable_units")
        for gene, unit in (("gene_a", first_unit), ("gene_b", second_unit)):
            independent = fit_family_rate([unit])
            recorded = result["null"]["families"][gene]
            self.assertEqual(recorded["status"], independent["status"])
            if independent.get("mu") is not None:
                self.assertAlmostEqual(recorded["mu"], independent["mu"], delta=1e-7)
            self.assertNotIn("pooled_mu", result["null"])

        def known_gene_rates(units, *args, family_id, foreground_multiplier=1., **kwargs):
            mu = {"gene_a": .2, "gene_b": .8}[family_id]
            return {"status": "estimated_conditional_composite_rate",
                    "converged": True, "mu": mu,
                    "log_likelihood": -((foreground_multiplier - 1.) ** 2) - mu}

        with patch("intraphy.inference.genomic_exon_comparison.fit_family_rate",
                   side_effect=known_gene_rates):
            independent_rate_case = compare_foreground(
                {"gene_a": [first_unit], "gene_b": [second_unit]}, {"A"},
                profile_multipliers=(.5, 1.))
        self.assertAlmostEqual(independent_rate_case["null"]["families"]["gene_a"]["mu"], .2)
        self.assertAlmostEqual(independent_rate_case["null"]["families"]["gene_b"]["mu"], .8)
        at_half = next(row for row in independent_rate_case["profile"]
                       if row["foreground_multiplier"] == .5)
        self.assertNotEqual(at_half["families"]["gene_a"]["mu"],
                            at_half["families"]["gene_b"]["mu"])

    def test_partial_tip_is_retained_and_serial_parallel_profiles_agree(self):
        catalogue = _catalogues("shared-deletion")[0]
        unit = _family_unit(catalogue)
        shifted_catalogue = _shift_catalogue(catalogue, 100)
        shifted_unit = _family_unit(shifted_catalogue)
        shifted_partial_tips = dict(shifted_unit["tips"])
        shifted_partial_tips["A"] = np.ones_like(shifted_partial_tips["A"])
        shifted_partial = {**shifted_unit, "tips": shifted_partial_tips}
        # Two disjoint units remain a within-family composite likelihood; this
        # checks scheduler parity, not independent-gene replication.
        family_units = {catalogue.family: [unit, shifted_partial]}
        serial = compare_foreground(family_units, {"A"},
            profile_multipliers=(0., .5, 1.), threads=1)
        parallel = compare_foreground(family_units, {"A"},
            profile_multipliers=(0., .5, 1.), threads=2)
        self.assertEqual(serial["counts"]["units"], 2)
        self.assertEqual(parallel["counts"]["units"], 2)
        self.assertEqual(serial["counts"]["requested_families"], 1)
        self.assertEqual(parallel["counts"]["requested_families"], 1)
        self.assertAlmostEqual(serial["null"]["log_likelihood"],
                               parallel["null"]["log_likelihood"], delta=1e-8)
        self.assertEqual(serial["null"]["families"], parallel["null"]["families"])
        self.assertEqual(serial["alternative"]["status"],
                         parallel["alternative"]["status"])
        self.assertEqual(serial["alternative"]["families"],
                         parallel["alternative"]["families"])
        if serial["alternative"]["log_likelihood"] is None:
            self.assertIsNone(parallel["alternative"]["log_likelihood"])
        else:
            self.assertAlmostEqual(serial["alternative"]["log_likelihood"],
                                   parallel["alternative"]["log_likelihood"], delta=1e-8)
        serial_profile = {row["foreground_multiplier"]: row for row in serial["profile"]}
        parallel_profile = {row["foreground_multiplier"]: row for row in parallel["profile"]}
        self.assertEqual(set(serial_profile), set(parallel_profile))
        for multiplier in serial_profile:
            self.assertEqual(serial_profile[multiplier]["status"],
                             parallel_profile[multiplier]["status"])
            self.assertEqual(serial_profile[multiplier]["families"],
                             parallel_profile[multiplier]["families"])
            left, right = (serial_profile[multiplier]["log_likelihood"],
                           parallel_profile[multiplier]["log_likelihood"])
            if left is not None:
                self.assertAlmostEqual(left, right, delta=1e-8)

    def test_canonical_foreground_selection_requires_background_exposure(self):
        catalogue = _catalogues("geometry")[0]
        unit = _family_unit(catalogue)
        children = set(unit["tree"].parent) - {unit["tree"].root}
        for foreground in ({"not-a-child"}, {unit["tree"].root}, children):
            with self.subTest(foreground=foreground), self.assertRaises(ValueError):
                compare_foreground({catalogue.family: [unit]}, foreground)

    def test_null_and_alternative_roster_cannot_diverge_by_tree_or_tip_panel(self):
        catalogue = _catalogues("geometry")[0]
        unit = _family_unit(catalogue)
        changed_tree = SpeciesTree([
            {"node_id": "root", "parent_id": "", "label": "root"},
            {"node_id": "A", "parent_id": "root", "label": "A", "branch_length": .6},
            {"node_id": "B", "parent_id": "root", "label": "B", "branch_length": .8},
        ])
        different_tree = {**unit, "unit_id": "different-tree", "tree": changed_tree}
        with self.assertRaisesRegex(ValueError, "same canonical tree and branch lengths"):
            compare_foreground({"gene_a": [unit], "gene_b": [different_tree]}, {"A"})
        different_tips = {**unit, "unit_id": "different-tips",
                          "tips": {"A": unit["tips"]["A"]}}
        with self.assertRaisesRegex(ValueError, "complete tip roster"):
            compare_foreground({"gene_a": [unit], "gene_b": [different_tips]}, {"A"})

    def test_all_excluded_request_retains_family_across_null_and_profile(self):
        result = compare_foreground({"requested_a": [], "requested_b": []}, {"unused"})
        self.assertEqual(result["counts"]["requested_families"], 2)
        self.assertEqual(result["counts"]["excluded_families"],
                         ["requested_a", "requested_b"])
        for point in [result["null"], result["alternative"], *result["profile"]]:
            self.assertEqual(set(point["families"]), {"requested_a", "requested_b"})
            self.assertEqual(point["status"], "no_estimable_units")
            self.assertEqual(point["unresolved_families"], ["requested_a", "requested_b"])
            self.assertTrue(all(record["status"] == "no_estimable_units"
                                for record in point["families"].values()))

    def test_flat_and_impossible_family_fits_remain_explicit(self):
        catalogue = _catalogues("geometry")[0]
        unit = _family_unit(catalogue)

        def flat_fit(*args, foreground_multiplier=1., **kwargs):
            return {"status": "no_rate_information", "converged": False,
                    "mu": None, "log_likelihood": -2.}

        with patch("intraphy.inference.genomic_exon_comparison.fit_family_rate",
                   side_effect=flat_fit):
            flat = compare_foreground({catalogue.family: [unit]}, {"A"})
        self.assertEqual(flat["null"]["families"][catalogue.family]["status"],
                         "no_rate_information")
        self.assertEqual(flat["alternative"]["status"],
                         "foreground_multiplier_unidentified")
        self.assertIsNone(flat["alternative"]["foreground_multiplier"])

        def impossible_fit(*args, **kwargs):
            return {"status": "impossible_observation_under_model", "converged": False,
                    "mu": None, "log_likelihood": None}

        with patch("intraphy.inference.genomic_exon_comparison.fit_family_rate",
                   side_effect=impossible_fit):
            impossible = compare_foreground({catalogue.family: [unit]}, {"A"})
        self.assertEqual(impossible["null"]["status"],
                         "impossible_observation_under_model")
        self.assertIn(catalogue.family, impossible["null"]["impossible_families"])
        self.assertEqual(impossible["alternative"]["status"], "profile_incomplete")
        self.assertEqual(impossible["alternative"]["resolution_reasons"],
                         ["profile_incomplete"])
        self.assertIsNone(impossible["p_value"])

    def test_unresolved_upper_profile_tail_is_not_reported_as_an_estimate(self):
        catalogue = _catalogues("geometry")[0]
        unit = _family_unit(catalogue)

        def rising_fit(*args, foreground_multiplier=1., **kwargs):
            return {"status": "estimated_conditional_composite_rate",
                    "converged": True, "mu": .2,
                    "log_likelihood": math.log1p(foreground_multiplier)}

        with patch("intraphy.inference.genomic_exon_comparison.fit_family_rate",
                   side_effect=rising_fit), patch(
                "intraphy.inference.genomic_exon_comparison.minimize",
                return_value=SimpleNamespace(x=np.asarray([2.]), success=True,
                    message="synthetic optimizer point", status=0, nit=1)):
            result = compare_foreground({catalogue.family: [unit]}, {"A"})
        self.assertEqual(result["alternative"]["status"], "upper_tail_unresolved")
        self.assertIsNone(result["alternative"]["foreground_multiplier"])
        self.assertEqual(result["alternative"]["upper_tail_diagnostic"]["reason"],
                         "higher_probe")

    def test_zero_rate_impossibility_is_distinct_from_a_nuisance_fit_failure(self):
        catalogue = _catalogues("geometry")[0]
        unit = _family_unit(catalogue)

        def partial_failure(*args, foreground_multiplier=1., **kwargs):
            if foreground_multiplier == 0.:
                return {"status": "impossible_observation_under_model", "converged": False,
                        "mu": None, "log_likelihood": None}
            if foreground_multiplier == .5:
                return {"status": "optimizer_failed", "converged": False,
                        "mu": None, "log_likelihood": None}
            return {"status": "estimated_conditional_composite_rate", "converged": True,
                    "mu": .2, "log_likelihood": -(foreground_multiplier - 1.) ** 2}

        def visit_incomplete_point(objective, initial, **kwargs):
            objective(np.asarray([.5]))
            return SimpleNamespace(x=np.asarray([1.]), success=True,
                message="synthetic optimum", status=0, nit=1)

        with patch("intraphy.inference.genomic_exon_comparison.fit_family_rate",
                   side_effect=partial_failure), patch(
                "intraphy.inference.genomic_exon_comparison.minimize",
                side_effect=visit_incomplete_point):
            result = compare_foreground({catalogue.family: [unit]}, {"A"},
                                        profile_multipliers=(.5,))
        points = {point["foreground_multiplier"]: point for point in result["profile"]}
        self.assertEqual(points[0.]["status"], "impossible_observation_under_model")
        self.assertEqual(points[.5]["status"], "profile_incomplete")
        self.assertEqual(points[.5]["families"][catalogue.family]["status"],
                         "optimizer_failed")
        self.assertIn(.5, result["alternative"]["incomplete_profile_points"])
        self.assertEqual(result["alternative"]["status"], "profile_incomplete")

    def test_mixed_impossible_family_and_failed_nuisance_are_separately_reported(self):
        first, second = _catalogues("geometry")
        units = {"gene_impossible": [_family_unit(first)],
                 "gene_fit_failure": [_family_unit(second)]}

        def mixed_fit(units, *args, family_id, foreground_multiplier=1., **kwargs):
            if family_id == "gene_impossible":
                return {"status": "impossible_observation_under_model", "converged": False,
                        "mu": None, "log_likelihood": None}
            if foreground_multiplier == .5:
                return {"status": "optimizer_failed", "converged": False,
                        "mu": None, "log_likelihood": None}
            return {"status": "estimated_conditional_composite_rate", "converged": True,
                    "mu": .3, "log_likelihood": -(foreground_multiplier - 1.) ** 2}

        def visit_half(objective, initial, **kwargs):
            objective(np.asarray([.5]))
            return SimpleNamespace(x=np.asarray([1.]), success=True,
                message="synthetic optimum", status=0, nit=1)

        with patch("intraphy.inference.genomic_exon_comparison.fit_family_rate",
                   side_effect=mixed_fit), patch(
                "intraphy.inference.genomic_exon_comparison.minimize",
                side_effect=visit_half):
            result = compare_foreground(units, {"A"}, profile_multipliers=(.5,))
        point = next(row for row in result["profile"]
                     if row["foreground_multiplier"] == .5)
        self.assertEqual(point["status"], "profile_incomplete")
        self.assertEqual(point["impossible_families"], ["gene_impossible"])
        self.assertEqual(point["unresolved_families"], ["gene_fit_failure"])

    def test_explicit_zero_boundary_and_unvisited_higher_profile_point_are_reported(self):
        catalogue = _catalogues("geometry")[0]
        unit = _family_unit(catalogue)

        def decreasing_fit(units, *args, foreground_multiplier=1., **kwargs):
            return {"status": "estimated_conditional_composite_rate",
                    "converged": True, "mu": .3,
                    "log_likelihood": -foreground_multiplier}

        with patch("intraphy.inference.genomic_exon_comparison.fit_family_rate",
                   side_effect=decreasing_fit), patch(
                "intraphy.inference.genomic_exon_comparison.minimize",
                return_value=SimpleNamespace(x=np.asarray([0.]), success=True,
                    message="zero optimum", status=0, nit=1)):
            zero = compare_foreground({catalogue.family: [unit]}, {"A"})
        self.assertEqual(zero["alternative"]["status"], "zero_boundary")
        self.assertEqual(zero["alternative"]["foreground_multiplier"], 0.)

        def peaked_fit(units, *args, foreground_multiplier=1., **kwargs):
            return {"status": "estimated_conditional_composite_rate",
                    "converged": True, "mu": .3,
                    "log_likelihood": -((foreground_multiplier - 10.) ** 2)}

        with patch("intraphy.inference.genomic_exon_comparison.fit_family_rate",
                   side_effect=peaked_fit), patch(
                "intraphy.inference.genomic_exon_comparison.minimize",
                return_value=SimpleNamespace(x=np.asarray([1.]), success=True,
                    message="suboptimal candidate", status=0, nit=1)):
            higher = compare_foreground({catalogue.family: [unit]}, {"A"})
        self.assertEqual(higher["alternative"]["status"], "optimization_unresolved")
        self.assertIn("optimization_unresolved",
                      higher["alternative"]["resolution_reasons"])

    def test_finite_mle_produces_resolved_multiplier_and_nested_statistic(self):
        catalogue = _catalogues("geometry")[0]
        unit = _family_unit(catalogue)

        def peaked_fit(units, *args, foreground_multiplier=1., **kwargs):
            return {"status": "estimated_conditional_composite_rate",
                    "converged": True, "mu": .3,
                    "log_likelihood": -((foreground_multiplier - 1.) ** 2)}

        def optimize_at_one(objective, initial, **kwargs):
            objective(np.asarray([1.]))
            return SimpleNamespace(x=np.asarray([1.]), success=True,
                message="finite optimum", status=0, nit=1)

        with patch("intraphy.inference.genomic_exon_comparison.fit_family_rate",
                   side_effect=peaked_fit), patch(
                "intraphy.inference.genomic_exon_comparison.minimize",
                side_effect=optimize_at_one):
            result = compare_foreground({catalogue.family: [unit]}, {"A"})
        self.assertEqual(result["alternative"]["status"],
                         "estimated_conditional_composite_multiplier")
        self.assertAlmostEqual(result["alternative"]["foreground_multiplier"], 1.)
        self.assertAlmostEqual(result["comparison_statistic"], 0., delta=1e-12)


if __name__ == "__main__":
    unittest.main()
