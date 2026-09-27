import unittest
import numpy as np

from intraphy.structure.repertoire import ExonRepertoire, repertoire_compatibility, delete_repertoire_material
from intraphy.structure.space import enumerate_space
from intraphy.structure.types import Catalogue, ExonConfiguration as C, ExonSpan as E, Material, ObservationEvidence as O, ConfigurationAlternative


class RepertoireObservationTests(unittest.TestCase):
    def setUp(self):
        self.empty, self.full = C((), (1,)), C((E(0, 10),), (1,))
        self.catalogue = Catalogue("g", "u", 10, (E(0, 10),), (),
                                   material=(Material("m", 0, 10),),
                                   observations=(), boundary_candidates=(E(0, 10),))
        self.space = enumerate_space(self.catalogue)

    def test_canonical_shared_material_and_duplicate_permutation_invariance(self):
        a = ExonRepertoire((self.full, self.empty, self.full), (1,))
        b = ExonRepertoire((self.empty, self.full), (1,))
        self.assertEqual(a, b)
        with self.assertRaises(ValueError): ExonRepertoire((self.empty, C((), (0,))), (1,))
        with self.assertRaises(ValueError): ExonRepertoire((C((), (True,)),), (1,))
        with self.assertRaises(ValueError): ExonRepertoire((), (1,))

    def test_coexisting_natives_require_distinct_members(self):
        obs = O("A", (self.empty, self.full), "coexisting", material_presence=(None,))
        one = ExonRepertoire((self.full,), (1,))
        two = ExonRepertoire((self.empty, self.full), (1,))
        values = repertoire_compatibility(self.space, (one, two), obs)
        np.testing.assert_array_equal(values, [0, 1])

    def test_duplicate_native_records_are_neutral(self):
        a = O("A", (self.empty, self.full, self.empty), "coexisting", material_presence=(None,))
        b = O("A", (self.empty, self.full), "coexisting", material_presence=(None,))
        candidate = ExonRepertoire((self.empty, self.full), (1,))
        np.testing.assert_array_equal(repertoire_compatibility(self.space, (candidate,), a), repertoire_compatibility(self.space, (candidate,), b))
        np.testing.assert_array_equal(repertoire_compatibility(self.space, (candidate,), a, complete_repertoire=True), repertoire_compatibility(self.space, (candidate,), b, complete_repertoire=True))

    def test_matching_reassigns_and_hall_violation_fails(self):
        alt = ConfigurationAlternative(self.full, self.empty.key, "B", "tx", "supported", 1.)
        obs = O("A", (self.empty, self.full), "coexisting", material_presence=(None,), alternatives=(alt,))
        candidate = ExonRepertoire((self.empty, self.full), (1,))
        np.testing.assert_array_equal(repertoire_compatibility(self.space, (candidate,), obs), [1])
        np.testing.assert_array_equal(repertoire_compatibility(self.space, (ExonRepertoire((self.full,), (1,)),), obs), [0])

    def test_lower_bound_allows_extra_and_complete_forbids_extra(self):
        obs = O("A", (self.full,), "observed", material_presence=(None,))
        repertoire = ExonRepertoire((self.empty, self.full), (1,))
        np.testing.assert_array_equal(repertoire_compatibility(self.space, (repertoire,), obs), [1])
        np.testing.assert_array_equal(repertoire_compatibility(self.space, (repertoire,), obs, complete_repertoire=True), [0])

    def test_observed_and_partial_use_or_compatibility(self):
        candidate = ExonRepertoire((self.full,), (1,))
        for kind in ("observed", "partial"):
            obs = O("A", (self.empty, self.full), kind, material_presence=(None,))
            np.testing.assert_array_equal(repertoire_compatibility(self.space, (candidate,), obs), [1])
        obs = O("A", (self.empty, self.full), "coexisting", material_presence=(None,))
        np.testing.assert_array_equal(repertoire_compatibility(self.space, (candidate,), obs), [0])

    def test_native_linked_alternative_is_preserved(self):
        alt = ConfigurationAlternative(self.full, self.empty.key, "B", "tx", "supported", 1.)
        obs = O("A", (self.empty, self.full), "coexisting", material_presence=(None,), alternatives=(alt,))
        np.testing.assert_array_equal(repertoire_compatibility(self.space, (ExonRepertoire((self.empty, self.full), (1,)),), obs), [1])

    def test_unknown_material_and_all_unknown(self):
        obs = O("A", (), "unknown", material_presence=(1,))
        np.testing.assert_array_equal(repertoire_compatibility(self.space, (ExonRepertoire((self.full,), (1,)),), obs), [1])
        unconstrained = O("A", (), "unknown", material_presence=(None,))
        np.testing.assert_array_equal(repertoire_compatibility(self.space, (ExonRepertoire((self.empty, self.full), (1,)),), unconstrained), [1])
        absent = O("A", (), "unknown", material_presence=(0,))
        np.testing.assert_array_equal(repertoire_compatibility(self.space, (ExonRepertoire((C((), (0,)),), (0,)), ExonRepertoire((C((), (2,)),), (2,))), absent), [1, 1])
        present = O("A", (), "unknown", material_presence=(1,))
        np.testing.assert_array_equal(repertoire_compatibility(self.space, (ExonRepertoire((C((), (1,)),), (1,)), ExonRepertoire((C((), (2,)),), (2,))), present), [1, 0])
        with self.assertRaises(ValueError): ExonRepertoire((self.empty,), (True,))

    def test_rejects_outside_member_and_invalid_complete_mode(self):
        outside = C((E(1, 2),), (1,))
        with self.assertRaisesRegex(ValueError, "outside"):
            repertoire_compatibility(self.space, (ExonRepertoire((outside,), (1,)),), O("A", (), "unknown", material_presence=(None,)))
        with self.assertRaises(TypeError): repertoire_compatibility(self.space, (), O("A", (), "unknown"), complete_repertoire=1)
        with self.assertRaises(ValueError): repertoire_compatibility(self.space, (), O("A", (), "partial", material_presence=(None,)), complete_repertoire=True)
        with self.assertRaises(ValueError): repertoire_compatibility(self.space, (), O("A", (), "excluded", material_presence=(None,)), complete_repertoire=True)

    def test_shared_material_deletion_changes_all_members_once(self):
        repertoire = ExonRepertoire((self.empty, self.full, self.full), (1,))
        deleted = delete_repertoire_material(self.catalogue, repertoire, "m")
        self.assertEqual(deleted.material, (2,))
        self.assertEqual(len(deleted.configurations), 1)
        self.assertEqual(deleted.configurations[0].exons, ())
        for state in (0, 2):
            candidate = C((), (state,))
            with self.assertRaises(ValueError):
                delete_repertoire_material(self.catalogue, ExonRepertoire((candidate,), (state,)), "m")
        with self.assertRaises(ValueError):
            delete_repertoire_material(
                self.catalogue, ExonRepertoire((C((E(1, 2),), (1,)),), (1,)), "m")
        incomplete = Catalogue("g", "u", 30, (E(0, 30),), (),
                               material=(Material("tail", 20, 30),),
                               boundary_candidates=(E(0, 30),))
        with self.assertRaisesRegex(ValueError, "leaves"):
            delete_repertoire_material(
                incomplete, ExonRepertoire((C((E(0, 30),), (1,)),), (1,)), "tail")

    def test_shared_deletion_normalizes_two_members_and_deduplicates(self):
        catalogue = Catalogue("g", "u", 30, (E(0, 10), E(20, 30), E(0, 30)), ((10, 20),),
                              material=(Material("m", 10, 20),), boundary_candidates=(E(0, 10), E(20, 30), E(0, 30)))
        split = C((E(0, 10), E(20, 30)), (1,)); full = C((E(0, 30),), (1,))
        deleted = delete_repertoire_material(catalogue, ExonRepertoire((split, full), (1,)), "m")
        self.assertEqual(deleted.material, (2,))
        self.assertEqual(len(deleted.configurations), 1)
        self.assertEqual(deleted.configurations[0].exons, (E(0, 30),))


if __name__ == "__main__": unittest.main()
