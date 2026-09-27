import unittest

from intraphy.structure.locus_observations import observation_emission
from intraphy.structure.locus_types import (
    CopySlot, LocusCatalogue, LocusObservation, LocusState,
    MaterialTract, SpliceFeature,
)


def simple_catalogue():
    material = (MaterialTract("dna", 0, 10),)
    copy = CopySlot("copy", ("dna",))
    feature = SpliceFeature("exon", "exon", ("dna",), copy_id="copy", start=0, end=4)
    return LocusCatalogue(material, (copy,), (feature,), (), provenance="observation fixture")


class LocusObservationTests(unittest.TestCase):
    def test_observed_dna_absence_supports_unintroduced_and_deleted_latent_states(self):
        catalogue = simple_catalogue()
        states = (LocusState((0,), frozenset()),
                  LocusState((1,), frozenset(("exon",))),
                  LocusState((2,), frozenset()))
        absent = observation_emission(catalogue, states, LocusObservation((0,)))
        unknown = observation_emission(catalogue, states, LocusObservation((None,)))
        self.assertEqual(absent.tolist(), [1.0, 0.0, 1.0])
        self.assertEqual(unknown.tolist(), [1.0, 1.0, 1.0])

    def test_unsurveyed_annotation_absence_is_uninformative(self):
        catalogue = simple_catalogue()
        states = (LocusState((1,), frozenset()),
                  LocusState((1,), frozenset(("exon",))))
        absent = LocusObservation((1,), features=(("exon", 0),))
        unknown = LocusObservation((1,), features=(("exon", None),))
        surveyed_absent = LocusObservation((1,), features=(("exon", 0),),
                                           surveyed_features=frozenset(("exon",)))
        self.assertEqual(observation_emission(catalogue, states, absent).tolist(), [1.0, 1.0])
        self.assertEqual(observation_emission(catalogue, states, unknown).tolist(), [1.0, 1.0])
        self.assertEqual(observation_emission(catalogue, states, surveyed_absent).tolist(), [1.0, 0.0])

    def test_surveyed_feature_detection_uses_declared_sensitivity_and_specificity(self):
        catalogue = simple_catalogue()
        states = (LocusState((1,), frozenset()),
                  LocusState((1,), frozenset(("exon",))))
        observation = LocusObservation((1,), features=(("exon", 0),),
                                       surveyed_features=frozenset(("exon",)))
        values = observation_emission(catalogue, states, observation,
                                      sensitivity={"exon": 0.8}, specificity={"exon": 0.9})
        self.assertAlmostEqual(float(values[0]), 0.9)
        self.assertAlmostEqual(float(values[1]), 0.2)

    def test_mxe_constraint_applies_to_each_observed_path_without_enumerating_isoforms(self):
        material = (MaterialTract("dna", 0, 12),)
        copy = CopySlot("copy", ("dna",))
        features = (SpliceFeature("A", "exon", ("dna",), copy_id="copy", start=0, end=2, mxe_group="g"),
                    SpliceFeature("C", "exon", ("dna",), copy_id="copy", start=4, end=6),
                    SpliceFeature("B", "exon", ("dna",), copy_id="copy", start=8, end=10, mxe_group="g"),
                    SpliceFeature("AC", "splice", ("dna",), prerequisites=("A", "C"), donor=2, acceptor=4),
                    SpliceFeature("CB", "splice", ("dna",), prerequisites=("C", "B"), donor=6, acceptor=8))
        catalogue = LocusCatalogue(material, (copy,), features, (),
                                   mxe_groups=(("g", ("A", "B")),), provenance="MXE fixture")
        state = LocusState((1,), frozenset(feature.id for feature in features))
        invalid_indirect_path = LocusObservation((1,), observed_paths=(("A", "C", "B"),))
        two_valid_paths = LocusObservation((1,), observed_paths=(("A", "C"), ("C", "B")))
        self.assertEqual(observation_emission(catalogue, (state,), invalid_indirect_path).tolist(), [0.0])
        self.assertEqual(observation_emission(catalogue, (state,), two_valid_paths).tolist(), [1.0])


if __name__ == "__main__":
    unittest.main()
