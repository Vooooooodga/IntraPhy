import math
import itertools
import unittest
import numpy as np
from intraphy.topology import SpeciesTree
from intraphy.structure.repertoire import ExonRepertoire
from intraphy.structure.repertoire_process import RepertoireProcess, RepertoireEvent
from intraphy.structure.types import Catalogue, ExonConfiguration, ExonSpan, Material, ObservationEvidence
from intraphy.structure.space import enumerate_space
from intraphy.structure.repertoire import repertoire_compatibility
from intraphy.inference.repertoire_model import evaluate_repertoire_model

RATES = {"dna_deletion": .7}

def tree(left=.4, right=.6):
    return SpeciesTree([{"node_id":"root","parent_id":"","label":"root","branch_length":"0"}, {"node_id":"a","parent_id":"root","label":"A","branch_length":str(left)}, {"node_id":"b","parent_id":"root","label":"B","branch_length":str(right)}])

def setup_states():
    exon = ExonSpan(0, 10); material = Material("m", 0, 10)
    catalogue = Catalogue("g", "u", 10, (exon,), (), material=(material,))
    before = ExonRepertoire((ExonConfiguration((), (1,)), ExonConfiguration((exon,), (1,))), (1,))
    after = ExonRepertoire((ExonConfiguration((), (2,)),), (2,))
    event = RepertoireEvent(before, after, "dna_deletion", "m")
    return catalogue, before, after, event

def evaluate(tips, *, process=None, prior=(1., 0.), tree_value=None, rates=RATES, **kwargs):
    catalogue, before, after, event = setup_states()
    p = process or RepertoireProcess((before, after), (event,), "shared deletion fixture", catalogue)
    return evaluate_repertoire_model(p, tree_value or tree(), tips, rates=rates, root_prior=np.asarray(prior), root_prior_provenance="declared root state", observation_provenance="fixture emissions", rate_provenance="explicit test rates", branch_length_unit="per unit branch length", **kwargs)

class RepertoireModelTests(unittest.TestCase):
    def test_two_way_turnover_matches_poisson_jump_oracle(self):
        # For a two-state symmetric CTMC with rate lambda, the total jump count
        # has Poisson(lambda*t) marginal law; E[N]=lambda*t, while endpoint difference for
        # the reversible pair is (1-exp(-2*lambda*t))/2.
        catalogue = Catalogue("oracle", "u", 10, (ExonSpan(0, 10),), ())
        r0 = ExonRepertoire((ExonConfiguration((), ()),), ())
        r1 = ExonRepertoire((ExonConfiguration((ExonSpan(0, 10),), ()),), ())
        process = RepertoireProcess((r0, r1), (
            RepertoireEvent(r0, r1, "turnover", "forward"),
            RepertoireEvent(r1, r0, "turnover", "backward")), "oracle", catalogue)
        x = .7
        out = evaluate_repertoire_model(process, tree(), {"A": np.ones(2), "B": np.ones(2)},
            rates={"turnover": x}, root_prior=np.array([1., 0.]),
            root_prior_provenance="oracle", observation_provenance="oracle",
            rate_provenance="oracle", branch_length_unit="unit")
        for edge, t in (("a", .4), ("b", .6)):
            item = out["branches"][("root", edge)]
            self.assertAlmostEqual(item["expected_edits"], x*t, places=10)
            self.assertAlmostEqual(item["expected_turnover"], x*t, places=10)
            self.assertAlmostEqual(item["probability_at_least_one_edit"], 1-math.exp(-x*t), places=10)
            self.assertAlmostEqual(item["probability_different_endpoints"], (1-math.exp(-2*x*t))/2, places=10)

    def test_even_poisson_endpoint_conditioning_oracle(self):
        # Conditioning a Poisson process on an even count gives
        # E[N|even]=x*tanh(x), P(any|even)=1-sech(x), and endpoint difference 0.
        catalogue = Catalogue("oracle", "u", 10, (ExonSpan(0, 10),), ())
        r0 = ExonRepertoire((ExonConfiguration((), ()),), ())
        r1 = ExonRepertoire((ExonConfiguration((ExonSpan(0, 10),), ()),), ())
        process = RepertoireProcess((r0, r1), (
            RepertoireEvent(r0, r1, "turnover", "forward"),
            RepertoireEvent(r1, r0, "turnover", "backward")), "oracle", catalogue)
        x, t = .7, .6
        out = evaluate_repertoire_model(process, tree(), {"A": np.ones(2), "B": np.array([1., 0.])},
            rates={"turnover": x}, root_prior=np.array([1., 0.]),
            root_prior_provenance="oracle", observation_provenance="oracle",
            rate_provenance="oracle", branch_length_unit="unit")
        item = out["branches"][("root", "b")]
        z = x*t
        self.assertAlmostEqual(item["expected_edits"], z*math.tanh(z), places=10)
        self.assertAlmostEqual(item["expected_turnover"], z*math.tanh(z), places=10)
        self.assertAlmostEqual(item["probability_at_least_one_edit"], 1-1/math.cosh(z), places=10)
        self.assertAlmostEqual(item["probability_different_endpoints"], 0., places=10)

    def test_shared_deletion_likelihood_and_marked_count(self):
        out = evaluate({"A":np.array([1.,0.]),"B":np.array([0.,1.])}, posterior=True, counts=True)
        expected = math.exp(-.7*.4) * (1-math.exp(-.7*.6))
        self.assertAlmostEqual(out["log_likelihood"], math.log(expected), places=12)
        self.assertAlmostEqual(out["branches"][("root","b")]["expected_dna_deletion"], 1., places=12)

    def test_all_unknown_and_all_zero(self):
        self.assertAlmostEqual(evaluate({"A":np.ones(2),"B":np.ones(2)})["log_likelihood"], 0., places=12)
        out = evaluate({"A":np.zeros(2),"B":np.zeros(2)})
        self.assertEqual(out["log_likelihood"], -math.inf); self.assertEqual(out["nodes"], {})

    def test_zero_length_conflict_is_impossible(self):
        out = evaluate({"A":np.array([1.,0.]),"B":np.array([0.,1.])}, tree_value=tree(0., 0.), posterior=False)
        self.assertEqual(out["log_likelihood"], -math.inf)

    def test_singleton_state_process(self):
        catalogue, before, _, _ = setup_states()
        p = RepertoireProcess((before,), (), "singleton fixture", catalogue)
        self.assertEqual(evaluate({"A":np.ones(1),"B":np.ones(1)}, process=p, prior=(1.,), rates={})["log_likelihood"], 0.)

    def test_invalid_prior_emission_labels_scale_and_provenance(self):
        with self.assertRaises(ValueError): evaluate({"A":np.ones(2),"B":np.ones(2)}, prior=(1.,))
        with self.assertRaises(ValueError): evaluate({"A":np.array([1.1,0.]),"B":np.ones(2)})
        with self.assertRaises(ValueError): evaluate({"A":np.ones(2),"C":np.ones(2)})
        with self.assertRaises(ValueError): evaluate({"A":np.ones(2),"B":np.ones(2)}, scale=0.)
        catalogue, before, after, _ = setup_states()
        p = RepertoireProcess((before, after), (), "fixture", catalogue)
        with self.assertRaises(ValueError): evaluate_repertoire_model(p, tree(), {"A":np.ones(2),"B":np.ones(2)}, rates=RATES, root_prior=np.array([1.,0.]), root_prior_provenance="", observation_provenance="x", rate_provenance="x", branch_length_unit="x")

    def test_permutation_preserves_likelihood(self):
        a = evaluate({"A":np.array([1.,0.]),"B":np.array([0.,1.])})
        catalogue, before, after, event = setup_states()
        p = RepertoireProcess((after, before), (RepertoireEvent(before, after, "dna_deletion", "m"),), "shared deletion fixture", catalogue)
        b = evaluate({"A":np.array([0.,1.]),"B":np.array([1.,0.])}, process=p, prior=(0.,1.))
        self.assertAlmostEqual(a["log_likelihood"], b["log_likelihood"], places=12)
        self.assertTrue(np.allclose(a["nodes"]["root"], b["nodes"]["root"][::-1]))

    def test_invalid_emissions_lengths_and_missing_branch(self):
        for value in (np.array([-1., 0.]), np.array([np.nan, 0.])):
            with self.subTest(value=value):
                with self.assertRaises(ValueError): evaluate({"A":value,"B":np.ones(2)})
        bad = tree(); bad.length["b"] = None
        with self.assertRaises(SystemExit): evaluate({"A":np.ones(2),"B":np.ones(2)}, tree_value=bad)
        for length in (-1., np.nan):
            bad = tree(); bad.length["b"] = length
            with self.subTest(length=length):
                with self.assertRaises(ValueError): evaluate({"A":np.ones(2),"B":np.ones(2)}, tree_value=bad)

    def test_invalid_root_vectors(self):
        for prior in ((-1., 2.), (np.nan, 1.), (0.5, 0.5000001)):
            with self.subTest(prior=prior):
                with self.assertRaises(ValueError): evaluate({"A":np.ones(2),"B":np.ones(2)}, prior=prior)

    def test_independent_two_state_enumeration_oracle(self):
        t = SpeciesTree([{"node_id":"r","parent_id":"","label":"r","branch_length":"0"}, {"node_id":"x","parent_id":"r","label":"x","branch_length":"0.2"}, {"node_id":"a","parent_id":"x","label":"A","branch_length":"0.4"}, {"node_id":"b","parent_id":"x","label":"B","branch_length":"0.6"}, {"node_id":"c","parent_id":"r","label":"C","branch_length":"0.5"}])
        emissions = {"A":np.array([1.,0.]),"B":np.array([0.,1.]),"C":np.ones(2)}
        out = evaluate(emissions, tree_value=t)
        nodes = t.preorder(); p = lambda length: np.array([[math.exp(-.7*length),1-math.exp(-.7*length)],[0.,1.]])
        total = 0.; node_sum = {n:np.zeros(2) for n in nodes}; edge_sum = {e:np.zeros((2,2)) for e in t.edges()}
        for assignment in itertools.product(range(2), repeat=len(nodes)):
            prob = .0 + (1. if assignment[0] == 0 else 0.)
            for parent, child in t.edges(): prob *= p(t.branch_length(child))[assignment[nodes.index(parent)], assignment[nodes.index(child)]]
            for leaf in t.leaves: prob *= emissions[t.label[leaf]][assignment[nodes.index(leaf)]]
            total += prob
            for n in nodes: node_sum[n][assignment[nodes.index(n)]] += prob
            for e in t.edges(): edge_sum[e][assignment[nodes.index(e[0])],assignment[nodes.index(e[1])]] += prob
        self.assertAlmostEqual(out["log_likelihood"], math.log(total), places=12)
        for n in nodes: self.assertTrue(np.allclose(out["nodes"][n], node_sum[n]/total, atol=1e-10))
        for e in t.edges(): self.assertTrue(np.allclose(out["branches"][e]["endpoint_matrix"], edge_sum[e]/total, atol=1e-10))

    def test_native_observation_constraints_feed_joint_emissions(self):
        catalogue, before, after, event = setup_states()
        space = enumerate_space(catalogue)
        obs_a = ObservationEvidence("A", before.configurations, "coexisting", material_presence=(1,))
        obs_b = ObservationEvidence("B", (), "unknown", material_presence=(0,))
        vectors_a = repertoire_compatibility(space, (before, after), obs_a)
        vectors_b = repertoire_compatibility(space, (before, after), obs_b)
        self.assertTrue(np.array_equal(vectors_a, [1., 0.]))
        self.assertTrue(np.array_equal(vectors_b, [0., 1.]))
        p = RepertoireProcess((before, after), (event,), "shared deletion fixture", catalogue)
        out = evaluate_repertoire_model(p, tree(), {"A":vectors_a,"B":vectors_b}, rates=RATES, root_prior=np.array([1.,0.]), root_prior_provenance="declared root state", observation_provenance="native compatibility vectors", rate_provenance="explicit test rates", branch_length_unit="per unit branch length")
        self.assertAlmostEqual(out["log_likelihood"], math.log(math.exp(-.7*.4) * (1-math.exp(-.7*.6))), places=12)
        self.assertAlmostEqual(out["branches"][("root","b")]["expected_dna_deletion"], 1., places=12)
