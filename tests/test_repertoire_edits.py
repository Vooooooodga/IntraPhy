import unittest
import numpy as np
from intraphy.structure.repertoire import ExonRepertoire
from intraphy.structure.repertoire_edits import coupled_process, POLICY, _structural_events
from intraphy.structure.repertoire_process import RepertoireEvent
from intraphy.structure.types import Catalogue, ExonConfiguration as C, ExonSpan as E, Material

class CoupledEditTests(unittest.TestCase):
    def setUp(self):
        self.cat = Catalogue("g", "u", 10, (E(0,10), E(0,4), E(6,10)), ((4,6),), boundary_candidates=(E(0,10),E(0,4),E(6,10)))
        self.seed = ExonRepertoire((C((E(0,10),)), C((E(0,4),E(6,10),))), ())

    def test_shared_process_and_option_cardinality(self):
        p = coupled_process(self.cat, (self.seed,), policy=POLICY, option_bank=(), insertion_outcomes=(), provenance="fixture")
        self.assertIn(self.seed, p.states)
        self.assertTrue(any(e.kind == "split" for e in p.events))

    def test_option_gain_rejects_material_mismatch(self):
        cat = Catalogue("m", "u", 10, (E(0,10),), (), material=(Material("m",0,10),), boundary_candidates=(E(0,10),))
        root = ExonRepertoire((C((), (0,)),), (0,))
        option = C((), (1,))
        p = coupled_process(cat, (root,), policy=POLICY, option_bank=(option,), insertion_outcomes=(), provenance="fixture")
        self.assertFalse(any(e.kind == "option_gain" for e in p.events))

    def test_option_gain_and_loss_change_repertoire_cardinality(self):
        option = C((E(0,4),), ())
        p = coupled_process(self.cat, (self.seed,), policy=POLICY, option_bank=(option,), insertion_outcomes=(), provenance="fixture")
        gains = [e for e in p.events if e.kind == "option_gain"]
        losses = [e for e in p.events if e.kind == "option_loss"]
        self.assertTrue(gains)
        self.assertTrue(losses)
        self.assertTrue(all(e.target.configurations for e in losses))

    def test_two_joint_targets_receive_uniform_opportunity_weights(self):
        cat = Catalogue("g", "u", 9, (E(0,5), E(0,7), E(0,9)), (), boundary_candidates=(E(0,5),E(0,7),E(0,9)))
        source = ExonRepertoire((C((E(0,5),)),), ())
        p = coupled_process(cat, (source,), policy=POLICY, option_bank=(), insertion_outcomes=(), provenance="fixture")
        grouped = [e for e in p.events if e.source == source and e.kind == "donor_shift" and e.opportunity == "donor:5"]
        self.assertEqual({next(iter(c.exons)).end for e in grouped for c in e.target.configurations}, {7, 9})
        self.assertEqual([e.weight for e in grouped], [.5, .5])

    def test_explicit_insertion_source_reachable_after_declared_seed_transition(self):
        cat = Catalogue("g", "u", 20, (E(0,5),), (), material=(Material("m",10,20),), boundary_candidates=(E(0,5),))
        seed = ExonRepertoire((C((E(0,5),), (0,)),), (0,))
        intermediate = ExonRepertoire((C((), (0,)),), (0,))
        target = ExonRepertoire((C((), (1,)),), (1,))
        insertion = (RepertoireEvent(intermediate, target, "dna_insertion", "m", 1.),)
        p = coupled_process(cat, (seed,), policy=POLICY, option_bank=(), insertion_outcomes=insertion, provenance="fixture")
        self.assertIn(intermediate, p.states)
        self.assertIn(target, p.states)
        self.assertTrue(any(e.source == intermediate and e.target == target and e.kind == "dna_insertion" for e in p.events))

    def test_shared_donor_event_maps_affected_members_and_preserves_unaffected(self):
        cat = Catalogue("g", "u", 25, (E(0,10),E(2,10),E(0,12),E(2,12),E(20,25)), ((10,12),), boundary_candidates=(E(0,10),E(2,10),E(0,12),E(2,12),E(20,25)))
        source = ExonRepertoire((C((E(0,10),)), C((E(2,10),)), C((E(20,25),))), ())
        events = _structural_events(cat, source)
        targets = [t for (kind, opp), ts in events.items() if kind == "donor_shift" for t in ts]
        self.assertTrue(any(E(0,12) in t.configurations[0].exons and E(2,12) in t.configurations[1].exons and E(20,25) in t.configurations[2].exons for t in targets))

    def test_shared_donor_blocked_when_affected_member_has_no_legal_target(self):
        cat = Catalogue("g", "u", 12, (E(0,10),E(8,10),E(0,5)), ((5,10),), boundary_candidates=(E(0,10),E(8,10),E(0,5)))
        source = ExonRepertoire((C((E(0,10),)), C((E(8,10),))), ())
        self.assertFalse(any(kind == "donor_shift" for kind, _ in _structural_events(cat, source)))

    def test_explicit_insertion_weights_and_invalid_table(self):
        cat = Catalogue("g", "u", 10, (E(0,10),), (), material=(Material("m",0,10),), boundary_candidates=(E(0,10),))
        source = ExonRepertoire((C((), (0,)),), (0,))
        empty = ExonRepertoire((C((), (1,)),), (1,)); full = ExonRepertoire((C((E(0,10),), (1,)),), (1,))
        table = (RepertoireEvent(source, empty, "dna_insertion", "m", .25), RepertoireEvent(source, full, "dna_insertion", "m", .75))
        p = coupled_process(cat, (source,), policy=POLICY, option_bank=(), insertion_outcomes=table, provenance="fixture")
        self.assertEqual(sorted(e.weight for e in p.events if e.kind == "dna_insertion"), [.25, .75])
        for weight in (.5, float("nan"), True):
            bad = (RepertoireEvent(source, empty, "dna_insertion", "m", weight),)
            with self.assertRaises((ValueError, ArithmeticError)): coupled_process(cat, (source,), policy=POLICY, option_bank=(), insertion_outcomes=bad, provenance="fixture")

    def test_shared_deletion_and_large_deletion_only_closure(self):
        cat = Catalogue("g", "u", 10, (E(0,10),), (), material=(Material("m",0,10),), boundary_candidates=(E(0,10),))
        source = ExonRepertoire((C((E(0,10),), (1,)), C((), (1,))), (1,))
        p = coupled_process(cat, (source,), policy=POLICY, option_bank=(), insertion_outcomes=(), provenance="fixture")
        source_deletions = [e for e in p.events if e.source == source and e.kind == "dna_deletion"]
        self.assertEqual(len(source_deletions), 1)
        self.assertEqual(source_deletions[0].weight, 1.)
        self.assertEqual(len(source_deletions[0].target.configurations), 1)
        materials = tuple(Material(f"m{i}", i, i+1) for i in range(7))
        large = Catalogue("g", "u", 7, (), (), material=materials, boundary_candidates=())
        seed = ExonRepertoire((C((), (1,)*7),), (1,)*7)
        q = coupled_process(large, (seed,), policy=POLICY, option_bank=(), insertion_outcomes=(), provenance="fixture")
        self.assertEqual(len(q.states), 128)
        self.assertEqual(sum(e.kind == "dna_deletion" for e in q.events), 448)
        with self.assertRaisesRegex(ValueError, "state_space_incomplete"):
            coupled_process(large, (seed,), policy=POLICY, option_bank=(), insertion_outcomes=(), provenance="fixture", max_states=64)

    def test_invalid_policy_and_cap(self):
        with self.assertRaises(ValueError): coupled_process(self.cat, (self.seed,), policy="", option_bank=(), insertion_outcomes=(), provenance="x")
        with self.assertRaisesRegex(ValueError, "state_space_incomplete"):
            coupled_process(self.cat, (self.seed,), policy=POLICY, option_bank=(), insertion_outcomes=(), provenance="x", max_states=1)

if __name__ == "__main__": unittest.main()
