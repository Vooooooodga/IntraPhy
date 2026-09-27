import unittest

import numpy as np

from intraphy.structure.locus_observations import observation_emission
from intraphy.structure.locus_types import (
    CopySlot, LocusCatalogue, LocusObservation, LocusState, MaterialTract,
)


def _fixture():
    material = (MaterialTract("copy_A", 0, 10), MaterialTract("copy_B", 10, 20))
    copies = (CopySlot("A", ("copy_A",)), CopySlot("B", ("copy_B",)))
    catalogue = LocusCatalogue(material, copies, (), (), provenance="detection fixture")
    states = tuple(LocusState(values, frozenset()) for values in (
        (0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)))
    return catalogue, states


class MaterialDetectionTests(unittest.TestCase):
    def test_hard_calls_unknown_and_soft_detection(self):
        catalogue, states = _fixture()
        hard = observation_emission(catalogue, states, LocusObservation((0, None)))
        self.assertEqual(hard.tolist(), [1.0, 1.0, 0.0, 0.0, 1.0, 1.0])
        unknown = observation_emission(catalogue, states, LocusObservation((None, None)))
        self.assertEqual(unknown.tolist(), [1.0] * len(states))

        soft = observation_emission(
            catalogue, states,
            LocusObservation((0, None), surveyed_material=frozenset({"copy_A"})),
            material_sensitivity={"copy_A": 0.8}, material_specificity={"copy_A": 0.9},
        )
        np.testing.assert_allclose(soft, [0.9, 0.9, 0.2, 0.2, 0.9, 0.9])

    def test_positive_call_uses_false_positive_probability_for_absent_states(self):
        catalogue, states = _fixture()
        result = observation_emission(
            catalogue, states,
            LocusObservation((1, None), surveyed_material=frozenset({"copy_A"})),
            material_sensitivity={"copy_A": 0.8}, material_specificity={"copy_A": 0.9},
        )
        np.testing.assert_allclose(result, [0.1, 0.1, 0.8, 0.8, 0.1, 0.1])

    def test_detection_parameters_require_survey_and_binary_known_call(self):
        catalogue, states = _fixture()
        with self.assertRaisesRegex(ValueError, "surveyed"):
            observation_emission(catalogue, states, LocusObservation((0, None)),
                                 material_sensitivity={"copy_A": 0.8},
                                 material_specificity={"copy_A": 0.9})
        with self.assertRaisesRegex(ValueError, "binary"):
            observation_emission(
                catalogue, states,
                LocusObservation((None, None), surveyed_material=frozenset({"copy_A"})),
                material_sensitivity={"copy_A": 0.8}, material_specificity={"copy_A": 0.9},
            )

    def test_invalid_material_detection_declarations_are_rejected(self):
        catalogue, states = _fixture()
        observation = LocusObservation((0, None), surveyed_material=frozenset({"copy_A"}))
        with self.assertRaisesRegex(ValueError, "same tracts"):
            observation_emission(catalogue, states, observation,
                                 material_sensitivity={"copy_A": 0.8})
        with self.assertRaisesRegex(ValueError, "undeclared"):
            observation_emission(catalogue, states, observation,
                                 material_sensitivity={"unknown": 0.8},
                                 material_specificity={"unknown": 0.9})
        with self.assertRaisesRegex(ValueError, r"\[0,1\]"):
            observation_emission(catalogue, states, observation,
                                 material_sensitivity={"copy_A": True},
                                 material_specificity={"copy_A": 0.9})


if __name__ == "__main__":
    unittest.main()
