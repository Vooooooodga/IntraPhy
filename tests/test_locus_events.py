import unittest
from dataclasses import replace

from intraphy.structure.locus_process import build_locus_process
from intraphy.structure.locus_types import (
    CopySlot, EventOpportunity, LocusCatalogue, LocusState,
    MaterialTract, SpliceFeature,
)


def duplication_catalogue(source_feature_active=True, with_deletion=True):
    material = (MaterialTract("src_dna", 0, 5), MaterialTract("target_dna", 10, 15),
                MaterialTract("flank_dna", 18, 23))
    copies = (CopySlot("source", ("src_dna",), evidence=("collinear block",)),
              CopySlot("target", ("target_dna",), homologous_to=("source",), evidence=("copy alignment",)),
              CopySlot("flank", ("flank_dna",), evidence=("annotation",)))
    features = (SpliceFeature("source_exon", "exon", ("src_dna",), copy_id="source",
                              start=0, end=5, evidence=("source annotation",)),
                SpliceFeature("target_exon", "exon", ("target_dna",), copy_id="target",
                              start=10, end=15, evidence=("target homology",)),
                SpliceFeature("flank_exon", "exon", ("flank_dna",), copy_id="flank",
                              start=18, end=23, evidence=("flank annotation",)),
                SpliceFeature("target_flank_splice", "splice", ("target_dna", "flank_dna"),
                              prerequisites=("target_exon", "flank_exon"), donor=15, acceptor=18))
    opportunities = [EventOpportunity(
        "duplication", "duplicate_target", "copy_duplication", "duplication", 1.0,
        preconditions=(("src_dna", 1), ("target_dna", 0), ("flank_dna", 1)),
        source_copy="source", target_copy="target", material_gains=("target_dna",),
        feature_map=(("source_exon", "target_exon"),), context_features=("flank_exon",),
        graft_features=("target_flank_splice",),
    )]
    if with_deletion:
        opportunities.append(EventOpportunity(
            "target_deletion", "delete_target", "dna_deletion", "deletion", 1.0,
            material_deletions=("target_dna",), interval=(10, 15)))
    catalogue = LocusCatalogue(material, copies, features, tuple(opportunities), provenance="test evidence")
    active = {"flank_exon"}
    if source_feature_active:
        active.add("source_exon")
    return catalogue, LocusState((1, 0, 1), frozenset(active))


class LocusEventTests(unittest.TestCase):
    def test_duplication_maps_features_and_conditionally_grafts_splicing(self):
        catalogue, active_root = duplication_catalogue(True, False)
        _, inactive_root = duplication_catalogue(False, False)
        process = build_locus_process(catalogue, (active_root, inactive_root))
        edges = [edge for edge in process.edges if edge.kind == "copy_duplication"]
        self.assertEqual(len(edges), 2)
        active_target = process.states[edges[0].target]
        inactive_target = process.states[edges[1].target]
        self.assertIn("target_exon", active_target.active_features)
        self.assertIn("target_flank_splice", active_target.active_features)
        self.assertNotIn("target_exon", inactive_target.active_features)
        self.assertNotIn("target_flank_splice", inactive_target.active_features)
        self.assertEqual(active_target.material, (1, 1, 1))
        self.assertEqual(inactive_target.material, (1, 1, 1))

    def test_shared_interval_deletion_is_one_event_and_preserves_unaffected_features(self):
        material = (MaterialTract("left", 0, 4), MaterialTract("right", 8, 12),
                    MaterialTract("untouched", 16, 20))
        copies = (CopySlot("copy", ("left", "right", "untouched")),)
        features = (SpliceFeature("left_exon", "exon", ("left",), copy_id="copy", start=0, end=4),
                    SpliceFeature("right_exon", "exon", ("right",), copy_id="copy", start=8, end=12),
                    SpliceFeature("other_exon", "exon", ("untouched",), copy_id="copy", start=16, end=20))
        deletion = EventOpportunity("del", "one_interval", "dna_deletion", "deletion", 1.0,
                                    material_deletions=("left", "right"), interval=(0, 12))
        catalogue = LocusCatalogue(material, copies, features, (deletion,), provenance="interval evidence")
        root = LocusState((1, 2, 1), frozenset(("left_exon", "other_exon")))
        process = build_locus_process(catalogue, (root,))
        self.assertEqual(len(process.edges), 1)
        target = process.states[process.edges[0].target]
        self.assertEqual(target.material, (2, 2, 1))
        self.assertEqual(target.active_features, frozenset(("other_exon",)))
        self.assertEqual(process.edges[0].kind, "dna_deletion")
        absent_root = LocusState((2, 2, 1), frozenset(("other_exon",)))
        null_process = build_locus_process(catalogue, (absent_root,))
        self.assertEqual(null_process.states, (absent_root,))
        self.assertEqual(null_process.edges, ())

    def test_deletion_rejects_a_skipped_modeled_tract_and_partial_overlap(self):
        material = (MaterialTract("a", 0, 3), MaterialTract("b", 5, 8),
                    MaterialTract("c", 10, 13))
        copy = CopySlot("copy", ("a", "b", "c"))
        features = tuple(SpliceFeature(f"e{i}", "exon", (tract.id,), copy_id="copy",
                                       start=tract.start, end=tract.end)
                        for i, tract in enumerate(material))
        skipped = EventOpportunity("bad", "skip_b", "dna_deletion", "deletion", 1.0,
                                   material_deletions=("a", "c"), interval=(0, 13))
        with self.assertRaisesRegex(ValueError, "every supplied material tract"):
            LocusCatalogue(material, (copy,), features, (skipped,), provenance="test")
        one = (MaterialTract("tract", 0, 10),)
        feature = SpliceFeature("exon", "exon", ("tract",), copy_id="copy", start=0, end=10)
        partial = EventOpportunity("partial", "boundary", "dna_deletion", "deletion", 1.0,
                                   material_deletions=("tract",), interval=(2, 8))
        with self.assertRaisesRegex(ValueError, "Split material tracts"):
            LocusCatalogue(one, (CopySlot("copy", ("tract",)),), (feature,), (partial,), provenance="test")

    def test_closure_has_no_default_cap_for_128_states(self):
        material = (MaterialTract("dna", 0, 20),)
        copy = CopySlot("copy", ("dna",))
        features = tuple(SpliceFeature(f"e{i}", "exon", ("dna",), copy_id="copy",
                                       start=2 * i, end=2 * i + 1) for i in range(7))
        opportunities = []
        for feature in features:
            opportunities.extend((
                EventOpportunity(f"on_{feature.id}", "on", "splice_change", "splice_on", feature_on=(feature.id,)),
                EventOpportunity(f"off_{feature.id}", "off", "splice_change", "splice_off", feature_off=(feature.id,)),
            ))
        catalogue = LocusCatalogue(material, (copy,), features, tuple(opportunities), provenance="toggle fixture")
        process = build_locus_process(catalogue, (LocusState((1,), frozenset()),))
        self.assertEqual(len(process.states), 128)
        self.assertEqual(len({state.active_features for state in process.states}), 128)

    def test_simultaneous_splice_boundary_change_uses_final_feature_set(self):
        material = (MaterialTract("dna", 0, 10),)
        copy = CopySlot("copy", ("dna",))
        features = (SpliceFeature("old_left", "exon", ("dna",), copy_id="copy", start=0, end=3),
                    SpliceFeature("new_left", "exon", ("dna",), copy_id="copy", start=0, end=4),
                    SpliceFeature("right", "exon", ("dna",), copy_id="copy", start=7, end=10),
                    SpliceFeature("old_splice", "splice", ("dna",), prerequisites=("old_left", "right"), donor=3, acceptor=7),
                    SpliceFeature("new_splice", "splice", ("dna",), prerequisites=("new_left", "right"), donor=4, acceptor=7))
        change = EventOpportunity("shift", "boundary", "splice_change", "splice",
                                  feature_on=("new_left", "new_splice"),
                                  feature_off=("old_left", "old_splice"))
        catalogue = LocusCatalogue(material, (copy,), features, (change,), provenance="boundary fixture")
        root = LocusState((1,), frozenset(("old_left", "right", "old_splice")))
        process = build_locus_process(catalogue, (root,))
        self.assertEqual(len(process.edges), 1)
        self.assertEqual(process.states[process.edges[0].target].active_features,
                         frozenset(("new_left", "right", "new_splice")))

    def test_partially_applicable_outcomes_fail_without_renormalization(self):
        material = (MaterialTract("dna", 0, 10), MaterialTract("context_dna", 12, 15),
                    MaterialTract("blocked_dna", 17, 20))
        copy = CopySlot("copy", ("dna", "context_dna", "blocked_dna"))
        features = (SpliceFeature("left", "exon", ("dna",), copy_id="copy", start=0, end=3),
                    SpliceFeature("right", "exon", ("dna",), copy_id="copy", start=7, end=10),
                    SpliceFeature("context", "exon", ("context_dna",), copy_id="copy", start=12, end=15),
                    SpliceFeature("blocked", "splice", ("blocked_dna",), prerequisites=("left", "right"), donor=3, acceptor=7))
        outcomes = (EventOpportunity("choose", "null", "splice_change", "splice", 0.5, feature_on=("left",)),
                    EventOpportunity("choose", "blocked", "splice_change", "splice", 0.5,
                                     feature_on=("blocked",)))
        catalogue = LocusCatalogue(material, (copy,), features, outcomes, provenance="alternative fixture")
        root = LocusState((1, 1, 0), frozenset(("left", "right", "context")))
        with self.assertRaisesRegex(ValueError, "partially applicable"):
            build_locus_process(catalogue, (root,))

    def test_null_alternative_retains_its_weight_without_rate_renormalization(self):
        material = (MaterialTract("dna", 0, 10),)
        copy = CopySlot("copy", ("dna",))
        features = (SpliceFeature("already", "exon", ("dna",), copy_id="copy", start=0, end=3),
                    SpliceFeature("available", "exon", ("dna",), copy_id="copy", start=5, end=8))
        outcomes = (EventOpportunity("choice", "null", "splice_change", "splice", 0.5,
                                     feature_on=("already",)),
                    EventOpportunity("choice", "gain", "splice_change", "splice", 0.5,
                                     feature_on=("available",)))
        catalogue = LocusCatalogue(material, (copy,), features, outcomes, provenance="null outcome fixture")
        root = LocusState((1,), frozenset(("already",)))
        process = build_locus_process(catalogue, (root,))
        edge = next(edge for edge in process.edges if edge.source == process.index[root])
        self.assertEqual(edge.weight, 0.5)
        self.assertEqual(process.states[edge.target].active_features,
                         frozenset(("already", "available")))

    def test_every_alternative_outcome_and_common_gate_are_validated(self):
        catalogue, _ = duplication_catalogue(True, False)
        valid = catalogue.opportunities[0]
        first = replace(valid, weight=0.5)
        malformed = replace(valid, outcome_id="bad_payload", weight=0.5,
                             material_gains=("src_dna",))
        with self.assertRaisesRegex(ValueError, "introduce the disjoint target"):
            LocusCatalogue(catalogue.material, catalogue.copies, catalogue.features,
                           (first, malformed), provenance="bad alternative")
        mismatched_gate = replace(valid, outcome_id="other_context", weight=0.5,
                                  context_features=())
        with self.assertRaisesRegex(ValueError, "share all applicability conditions"):
            LocusCatalogue(catalogue.material, catalogue.copies, catalogue.features,
                           (first, mismatched_gate), provenance="bad gate")

    def test_duplication_cannot_resurrect_deleted_target_material(self):
        catalogue, root = duplication_catalogue(True, True)
        process = build_locus_process(catalogue, (root,))
        for edge in process.edges:
            source = process.states[edge.source]
            target = process.states[edge.target]
            if source.material[1] == 2:
                self.assertNotEqual(target.material[1], 1)
