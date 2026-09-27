import unittest
from dataclasses import replace

from intraphy.structure.locus_events import apply_outcome
from intraphy.structure.locus_process import build_locus_process
from intraphy.structure.locus_types import CopySlot, EventOpportunity, LocusCatalogue, LocusState, MaterialTract


def copy_catalogue(state_model="irreversible"):
    material = (MaterialTract("source_dna", 0, 5), MaterialTract("target_dna", 10, 15))
    copies = (CopySlot("source", ("source_dna",), evidence=("sequence",)),
              CopySlot("target", ("target_dna",), homologous_to=("source",), evidence=("alignment",)))
    duplicate = EventOpportunity(
        "dup", "copy", "copy_duplication", "duplication",
        preconditions=(("source_dna", 1), ("target_dna", 0)),
        material_gains=("target_dna",), source_copy="source", target_copy="target")
    delete = EventOpportunity(
        "del", "target", "dna_deletion", "deletion",
        material_deletions=("target_dna",), interval=(10, 15))
    return LocusCatalogue(material, copies, (duplicate, delete),
                          provenance="test alignment", state_model=state_model)


class LocusEventTests(unittest.TestCase):
    def test_two_models_have_distinct_post_deletion_closure(self):
        for model, expected, edges in (
            ("irreversible", {(1, 0), (1, 1), (1, 2)}, 2),
            ("binary", {(1, 0), (1, 1)}, 2),
        ):
            with self.subTest(model=model):
                catalogue = copy_catalogue(model)
                process = build_locus_process(catalogue, (LocusState((1, 0)),))
                self.assertEqual({state.material for state in process.states}, expected)
                self.assertEqual(len(process.edges), edges)
                self.assertEqual(catalogue.material_states, (0, 1) if model == "binary" else (0, 1, 2))
                self.assertTrue(all(not hasattr(state, "active_features") for state in process.states))
        binary = copy_catalogue("binary")
        self.assertEqual(apply_outcome(binary, LocusState((1, 0)), binary.opportunities[0]),
                         LocusState((1, 1)))
        self.assertEqual(apply_outcome(binary, LocusState((1, 1)), binary.opportunities[1]),
                         LocusState((1, 0)))
        irreversible = copy_catalogue()
        self.assertEqual(apply_outcome(irreversible, LocusState((1, 1)), irreversible.opportunities[1]),
                         LocusState((1, 2)))

    def test_interval_deletes_all_present_intersecting_material_once(self):
        material = (MaterialTract("a", 0, 4), MaterialTract("b", 8, 12),
                    MaterialTract("untouched", 16, 20))
        copies = (CopySlot("copy", ("a", "b", "untouched")),)
        deletion = EventOpportunity("del", "one_interval", "dna_deletion", "deletion",
                                    material_deletions=("a", "b"), interval=(0, 12))
        catalogue = LocusCatalogue(material, copies, (deletion,), provenance="interval evidence")
        process = build_locus_process(catalogue, (LocusState((1, 2, 1)),))
        self.assertEqual(len(process.edges), 1)
        self.assertEqual(process.states[process.edges[0].target], LocusState((2, 2, 1)))
        self.assertEqual(build_locus_process(catalogue, (LocusState((2, 2, 1)),)).edges, ())

    def test_deletion_boundaries_and_skipped_tract_are_rejected(self):
        material = (MaterialTract("a", 0, 3), MaterialTract("b", 5, 8),
                    MaterialTract("c", 10, 13))
        copies = (CopySlot("copy", ("a", "b", "c")),)
        skipped = EventOpportunity("bad", "skip", "dna_deletion", "deletion",
                                   material_deletions=("a", "c"), interval=(0, 13))
        with self.assertRaisesRegex(ValueError, "every supplied material tract"):
            LocusCatalogue(material, copies, (skipped,), provenance="test")
        partial = EventOpportunity("bad", "partial", "dna_deletion", "deletion",
                                   material_deletions=("a",), interval=(1, 3))
        with self.assertRaisesRegex(ValueError, "Split material tracts"):
            LocusCatalogue(material, copies, (partial,), provenance="test")

    def test_full_reachable_closure_has_no_default_cap(self):
        material = [MaterialTract("source", 0, 1)]
        copies = [CopySlot("source", ("source",), evidence=("sequence",))]
        opportunities = []
        for i in range(7):
            mid, cid = f"target_{i}", f"copy_{i}"
            material.append(MaterialTract(mid, 2 * i + 2, 2 * i + 3))
            copies.append(CopySlot(cid, (mid,), homologous_to=("source",), evidence=("alignment",)))
            opportunities.append(EventOpportunity(
                f"dup_{i}", "copy", "copy_duplication", "duplication",
                preconditions=(("source", 1), (mid, 0)),
                material_gains=(mid,), source_copy="source", target_copy=cid))
        catalogue = LocusCatalogue(tuple(material), tuple(copies), tuple(opportunities), provenance="seven alignments")
        process = build_locus_process(catalogue, (LocusState((1,) + (0,) * 7),))
        self.assertEqual(len(process.states), 128)
        self.assertEqual(len(process.edges), 7 * 64)

    def test_partial_alternative_deletions_fail_without_weight_renormalization(self):
        material = (MaterialTract("a", 0, 2), MaterialTract("b", 4, 6))
        copies = (CopySlot("copy", ("a", "b")),)
        outcomes = (
            EventOpportunity("choice", "a", "dna_deletion", "deletion", weight=0.5,
                             material_deletions=("a",), interval=(0, 2)),
            EventOpportunity("choice", "b", "dna_deletion", "deletion", weight=0.5,
                             material_deletions=("b",), interval=(4, 6)),
        )
        catalogue = LocusCatalogue(material, copies, outcomes, provenance="two intervals")
        with self.assertRaisesRegex(ValueError, "partially applicable"):
            build_locus_process(catalogue, (LocusState((1, 0)),))

    def test_duplication_requires_supplied_homology_evidence_and_preconditions(self):
        catalogue = copy_catalogue()
        valid = catalogue.opportunities[0]
        with self.assertRaisesRegex(ValueError, "homology relation"):
            LocusCatalogue(catalogue.material,
                           (catalogue.copies[0], replace(catalogue.copies[1], homologous_to=())),
                           catalogue.opportunities, provenance="test")
        with self.assertRaisesRegex(ValueError, "evidence on both"):
            LocusCatalogue(catalogue.material,
                           (replace(catalogue.copies[0], evidence=()), catalogue.copies[1]),
                           catalogue.opportunities, provenance="test")
        with self.assertRaisesRegex(ValueError, "present source and absent target"):
            LocusCatalogue(catalogue.material, catalogue.copies,
                           (replace(valid, preconditions=()),), provenance="test")
        with self.assertRaisesRegex(ValueError, "selected state model"):
            LocusCatalogue(catalogue.material, catalogue.copies,
                           (replace(valid, preconditions=(("source_dna", 2), ("target_dna", 0))),),
                           provenance="test", state_model="binary")
        with self.assertRaisesRegex(ValueError, "root state is invalid"):
            build_locus_process(copy_catalogue("binary"), (LocusState((1, 2)),))
