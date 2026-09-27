import unittest
import numpy as np
from intraphy.structure.repertoire import ExonRepertoire
from intraphy.structure.repertoire_process import RepertoireEvent, RepertoireProcess, deletion_process
from intraphy.structure.types import Catalogue, ExonConfiguration as C, ExonSpan as E, Material

class RepertoireProcessTests(unittest.TestCase):
    def setUp(self):
        self.cat = Catalogue("g", "u", 10, (E(0, 10),), (), material=(Material("m", 0, 10),), boundary_candidates=(E(0, 10),))
        self.a = ExonRepertoire((C((E(0, 10),), (1,)),), (1,))
        self.b = ExonRepertoire((C((), (2,)),), (2,))

    def test_generator_aggregates_and_marks(self):
        p = RepertoireProcess((self.a, self.b), (RepertoireEvent(self.a, self.b, "dna_deletion", "m"),), "synthetic", self.cat)
        q, marks = p.generator({"dna_deletion": 2.})
        np.testing.assert_allclose(q, [[-2., 2.], [0., 0.]])
        np.testing.assert_allclose(marks["dna_deletion"], [[0., 2.], [0., 0.]])

    def test_opportunity_weight_and_state_validation(self):
        same = ExonRepertoire((C((), (1,)),), (1,))
        with self.assertRaisesRegex(ValueError, "sum to one"): RepertoireProcess((self.a, same), (RepertoireEvent(self.a, same, "x", "m", .5),), "x", self.cat)
        with self.assertRaises(ValueError): RepertoireProcess((self.a, self.b), (RepertoireEvent(self.a, self.b, "x", "m"),), "", self.cat)

    def test_duplicate_and_conflicting_events(self):
        e = RepertoireEvent(self.a, self.b, "dna_deletion", "m")
        p = RepertoireProcess((self.a, self.b), (e, e), "x", self.cat)
        self.assertEqual(len(p.events), 1)
        with self.assertRaises(ValueError): RepertoireProcess((self.a, self.b), (e, RepertoireEvent(self.a, self.b, "dna_deletion", "m", .5)), "x", self.cat)

    def test_distinct_opportunities_and_weighted_outcomes(self):
        empty = ExonRepertoire((C((), (1,)),), (1,)); split = ExonRepertoire((C((E(0, 10),), (1,)),), (1,))
        cat = self.cat
        events = (RepertoireEvent(empty, split, "split", "o1"), RepertoireEvent(empty, split, "split", "o2"))
        p = RepertoireProcess((empty, split), events, "x", cat)
        np.testing.assert_allclose(p.generator({"split": 2.})[0][0], [-4., 4.])
        left = ExonRepertoire((C((E(0, 10),), (1,)),), (1,)); right = ExonRepertoire((C((), (1,)),), (1,))
        with self.assertRaises(ValueError): RepertoireProcess((left, right), (RepertoireEvent(left, right, "split", "o", .25),), "x", cat)

    def test_one_opportunity_weighted_targets_aggregates(self):
        empty = ExonRepertoire((C((), (1,)),), (1,)); full = ExonRepertoire((C((E(0, 10),), (1,)),), (1,))
        both = ExonRepertoire((C((), (1,)), C((E(0, 10),), (1,))), (1,))
        p = RepertoireProcess((empty, full, both),
                              (RepertoireEvent(empty, full, "x", "o", .25), RepertoireEvent(empty, both, "x", "o", .75)), "x", self.cat)
        np.testing.assert_allclose(p.generator({"x": 2.})[0][0], [-2., .5, 1.5])

    def test_order_and_endpoint_and_lifecycle_validation(self):
        with self.assertRaises(ValueError): RepertoireProcess((self.a, self.a), (), "x", self.cat)
        with self.assertRaises(ValueError): RepertoireProcess((self.a, self.b), (RepertoireEvent(self.b, self.a, "dna_deletion", "m"),), "x", self.cat)
        with self.assertRaises(ValueError): RepertoireProcess((self.a, self.b), (RepertoireEvent(self.a, self.b, "dna_deletion", "m"),), "x", Catalogue("g", "u", 10, (), (), material=(), status="qualified"))

    def test_generator_input_and_overflow_guards(self):
        p = RepertoireProcess((self.a, self.b), (RepertoireEvent(self.a, self.b, "dna_deletion", "m"),), "x", self.cat)
        for rates in ({}, {"dna_deletion": -1.}, {"dna_deletion": float("nan")}, {"dna_deletion": True}):
            with self.assertRaises(ValueError): p.generator(rates)
        with self.assertRaises(ValueError): p.generator({"dna_deletion": 1.}, scale=True)
        with self.assertRaises(TypeError): p.generator({"dna_deletion": 1.}, include_marks=1)
        with self.assertRaises(ArithmeticError): p.generator({"dna_deletion": 1e308}, scale=1e308)

    def test_shared_deletion_closure_seven_materials_and_cap(self):
        materials = tuple(Material(f"m{i}", i * 2, i * 2 + 1) for i in range(7))
        cat = Catalogue("g", "u", 15, tuple(E(i * 2, i * 2 + 1) for i in range(7)), (), material=materials,
                        boundary_candidates=tuple(E(i * 2, i * 2 + 1) for i in range(7)))
        seed = ExonRepertoire((C((), (1,) * 7),), (1,) * 7)
        p = deletion_process(cat, (seed, seed), provenance="x")
        self.assertEqual(len(p.states), 128); self.assertEqual(len(p.events), 448)
        with self.assertRaisesRegex(ValueError, "state_space_incomplete"): deletion_process(cat, (seed,), provenance="x", max_states=64)

    def test_shared_deletion_two_members_is_one_hazard_and_bad_seed_rejected(self):
        split = ExonRepertoire((C((E(0, 10),), (1,)), C((), (1,))), (1,))
        process = deletion_process(self.cat, (split,), provenance="x")
        self.assertEqual(len(process.events), 1)
        bad = ExonRepertoire((C((E(1, 2),), (1,)),), (1,))
        with self.assertRaises(ValueError): deletion_process(self.cat, (bad,), provenance="x")

    def test_deletion_closure_and_cap(self):
        p = deletion_process(self.cat, (self.a,), provenance="declared deletion table")
        self.assertEqual(len(p.states), 2)
        self.assertEqual(len(p.events), 1)
        with self.assertRaisesRegex(ValueError, "state_space_incomplete"):
            deletion_process(self.cat, (self.a,), provenance="x", max_states=1)

if __name__ == "__main__": unittest.main()
