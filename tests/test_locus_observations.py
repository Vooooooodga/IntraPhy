import unittest

from intraphy.structure.locus_observations import observation_emission
from intraphy.structure.locus_types import CopySlot, LocusCatalogue, LocusObservation, LocusState, MaterialTract


def simple_catalogue(state_model="irreversible"):
    return LocusCatalogue((MaterialTract("dna", 0, 10),),
                          (CopySlot("copy", ("dna",)),), (),
                          provenance="tip evidence", state_model=state_model)


class LocusObservationTests(unittest.TestCase):
    def test_observed_absence_supports_both_irreversible_latent_absences(self):
        catalogue = simple_catalogue()
        states = (LocusState((0,)), LocusState((1,)), LocusState((2,)))
        self.assertEqual(observation_emission(catalogue, states, LocusObservation((0,))).tolist(),
                         [1.0, 0.0, 1.0])
        self.assertEqual(observation_emission(catalogue, states, LocusObservation((None,))).tolist(),
                         [1.0, 1.0, 1.0])

    def test_binary_observation_uses_same_presence_call(self):
        catalogue = simple_catalogue("binary")
        states = (LocusState((0,)), LocusState((1,)))
        self.assertEqual(observation_emission(catalogue, states, LocusObservation((1,))).tolist(),
                         [0.0, 1.0])
        with self.assertRaisesRegex(ValueError, "does not match locus catalogue"):
            observation_emission(catalogue, states + (LocusState((2,)),), LocusObservation((None,)))

    def test_survey_metadata_does_not_turn_missing_call_into_absence(self):
        catalogue = simple_catalogue()
        states = (LocusState((0,)), LocusState((1,)))
        observation = LocusObservation((None,), surveyed_material=frozenset({"dna"}))
        self.assertEqual(observation_emission(catalogue, states, observation).tolist(), [1.0, 1.0])

    def test_invalid_observation_alphabet_and_unknown_survey_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "0, 1, or None"):
            LocusObservation((2,))
        with self.assertRaisesRegex(ValueError, "undeclared"):
            observation_emission(simple_catalogue(), (LocusState((1,)),),
                                 LocusObservation((1,), surveyed_material=frozenset({"other"})))
