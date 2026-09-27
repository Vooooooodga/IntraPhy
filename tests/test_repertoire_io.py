"""Strict repertoire-model I/O fixtures; intentionally not run in this phase."""
import json
import tempfile
import unittest
from pathlib import Path
import numpy as np

from intraphy.inference.repertoire_inputs import read_repertoire_model, encode_repertoire, SCHEMA
from intraphy.structure.repertoire import ExonRepertoire
from intraphy.structure.repertoire_process import deletion_process, RepertoireProcess
from dataclasses import replace
from intraphy.structure.types import Catalogue, ExonConfiguration as C, ExonSpan as E, Material, ObservationEvidence as O
from intraphy.topology import SpeciesTree


def _fixture():
    c = Catalogue("g", "u", 10, (E(0, 10),), (), material=(Material("m", 0, 10),), boundary_candidates=(E(0, 10),))
    seed = ExonRepertoire((C((E(0, 10),), (1,)),), (1,))
    tree = SpeciesTree([{"node_id":"r","parent_id":"","label":"r"}, {"node_id":"A","parent_id":"r","label":"A","branch_length":1.}])
    return c, seed, tree

def _record(c, process, root_states, rates, observation):
    states=[encode_repertoire(s) for s in process.states]
    events=[]
    for e in process.events:
        events.append({"source":process.states.index(e.source),"target":process.states.index(e.target),"kind":e.kind,"opportunity":e.opportunity,"weight":e.weight})
    return {"schema":SCHEMA,"model":"exon-repertoire-ctmc","branch_length_unit":"unit","scale":1.,"provenance":"fixture","units":[{"family":c.family,"unit":c.unit,"process":{"mode":"explicit","states":states,"events":events,"provenance":process.provenance},"rates":rates,"rate_provenance":"r","root":{"entries":[{"state":states[i],"weight":1/len(root_states)} for i in root_states],"provenance":"root"},"observation":observation}]}

def _read(record, c, tree):
    with tempfile.TemporaryDirectory() as d:
        p=Path(d)/"model.json"; p.write_text(json.dumps(record)); return read_repertoire_model(p,(c,),tree)


class RepertoireIOTests(unittest.TestCase):
    def test_strict_deletion_roundtrip_and_root_order(self):
        c, seed, tree = _fixture(); deleted = deletion_process(c, (seed,), provenance="synthetic")
        states = [encode_repertoire(s) for s in deleted.states]
        record = {"schema":SCHEMA,"model":"exon-repertoire-ctmc","branch_length_unit":"unit",
                  "scale":1.,"provenance":"fixture","units":[{"family":"g","unit":"u",
                  "process":{"mode":"explicit","states":states,"events":[{"source":0,"target":1,"kind":"dna_deletion","opportunity":"m","weight":1.}],"provenance":"p"},
                  "rates":{"dna_deletion":1.},"rate_provenance":"r","root":{"entries":[{"state":states[1],"weight":0.25},{"state":states[0],"weight":0.75}],"provenance":"root"},
                  "observation":{"mode":"explicit","provenance":"userfixedemissions","tips":{"A":[1.,0.]}}}]}
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/"model.json"; p.write_text(json.dumps(record)); bundle=read_repertoire_model(p,(c,),tree)
        self.assertEqual(bundle.units[0].root_prior.tolist(), [0.75,0.25])

    def test_unknown_json_keys_duplicate_keys_and_nonfinite_rejected(self):
        c, seed, tree = _fixture()
        for payload in ["{\"schema\":1,\"schema\":2}", "{\"schema\":\"x\",\"extra\":1}"]:
            with tempfile.TemporaryDirectory() as d:
                p=Path(d)/"bad.json"; p.write_text(payload)
                with self.assertRaises(ValueError): read_repertoire_model(p,(c,),tree)

    def test_constraints_missing_leaf_emits_ones(self):
        c, seed, tree = _fixture(); process=deletion_process(c,(seed,),provenance="p")
        obs={"mode":"constraints","provenance":"declared","view":"annotation","missing_policy":"uninformative","complete_species":[]}
        bundle=_read(_record(c,process,range(len(process.states)),{"dna_deletion":1.},obs),c,tree)
        np.testing.assert_array_equal(bundle.units[0].tips["A"],np.ones(len(process.states)))

    def test_deletion_mode_decodes_two_states_one_event(self):
        c, seed, tree = _fixture(); process=deletion_process(c,(seed,),provenance="p")
        rec=_record(c,process,range(len(process.states)),{"dna_deletion":1.},{"mode":"constraints","provenance":"declared","view":"annotation","missing_policy":"uninformative","complete_species":[]})
        rec["units"][0]["process"]={"mode":"deletion","seeds":[encode_repertoire(seed)],"provenance":"p"}; rec["units"][0]["rates"]={"dna_deletion":1.}
        out=_read(rec,c,tree); self.assertEqual(len(out.units[0].process.states),2); self.assertEqual(len(out.units[0].process.events),1)

    def test_malformed_containers_complete_species_and_huge_numbers_reject(self):
        c, seed, tree = _fixture(); process=deletion_process(c,(seed,),provenance="p")
        base=_record(c,process,range(len(process.states)),{"dna_deletion":1.},{"mode":"constraints","provenance":"declared","view":"annotation","missing_policy":"uninformative","complete_species":[]})
        cases=[]
        bad=json.loads(json.dumps(base)); bad["units"][0]["process"]["states"]=None; cases.append(bad)
        bad=json.loads(json.dumps(base)); bad["units"][0]["observation"]["complete_species"]=[[]]; cases.append(bad)
        bad=json.loads(json.dumps(base)); bad["scale"]=10**400; cases.append(bad)
        for field, values in (("configurations", (None, {})), ("material", (None, {}))):
            for value in values:
                bad=json.loads(json.dumps(base)); bad["units"][0]["process"]["states"][0][field]=value; cases.append(bad)
        for field, values in (("exons", (None, {})), ("material", (None, {}))):
            for value in values:
                bad=json.loads(json.dumps(base)); bad["units"][0]["process"]["states"][0]["configurations"][0][field]=value; cases.append(bad)
        for value in ([], {}):
            bad=json.loads(json.dumps(base)); bad["units"][0]["process"]["mode"]=value; cases.append(bad)
        for field, value in (("mode", []), ("view", {}), ("missing_policy", [])):
            bad=json.loads(json.dumps(base)); bad["units"][0]["observation"][field]=value; cases.append(bad)
        for record in cases:
            with self.assertRaises(ValueError): _read(record,c,tree)

    def test_complete_species_forbids_unannotated_extra_member(self):
        c, seed, tree = _fixture(); full=seed.configurations[0]; empty=C((),(1,))
        c=replace(c, observations=(O("A", (full,), material_presence=(None,)),))
        states=(ExonRepertoire((full,), (1,)), ExonRepertoire((full,empty), (1,)))
        process=RepertoireProcess(states, (), "p", c)
        base={
            "schema":SCHEMA,
            "model":"exon-repertoire-ctmc",
            "branch_length_unit":"unit",
            "scale":1.,
            "provenance":"fixture",
            "units":[{
                "family":"g",
                "unit":"u",
                "process":{
                    "mode":"explicit",
                    "states":[encode_repertoire(s) for s in states],
                    "events":[],
                    "provenance":"p",
                },
                "rates":{},
                "rate_provenance":"r",
                "root":{
                    "entries":[{"state":encode_repertoire(states[0]),"weight":1.}],
                    "provenance":"root",
                },
            }],
        }
        for complete in ([], ["A"]):
            record=json.loads(json.dumps(base)); record["units"][0]["observation"]={"mode":"constraints","provenance":"declared","view":"annotation","missing_policy":"uninformative","complete_species":complete}
            out=_read(record,c,tree)
            np.testing.assert_array_equal(out.units[0].tips["A"], [1., 1.] if not complete else [1., 0.])
            self.assertEqual(out.observation_metadata[("g","u")]["complete_species"], complete)

    def test_coupled_insertion_weights_are_preserved(self):
        c, _seed, tree = _fixture(); c=replace(c, material=(Material("m",0,10),))
        source=ExonRepertoire((C((),(0,)),),(0,)); empty=ExonRepertoire((C((),(1,)),),(1,)); full=ExonRepertoire((C((E(0,10),),(1,)),),(1,))
        states=[encode_repertoire(source),encode_repertoire(empty),encode_repertoire(full)]
        event=lambda target,weight:{"source":encode_repertoire(source),"target":encode_repertoire(target),"kind":"dna_insertion","opportunity":"m","weight":weight}
        record={"schema":SCHEMA,"model":"exon-repertoire-ctmc","branch_length_unit":"unit","scale":1.,"provenance":"fixture","units":[{"family":"g","unit":"u","process":{"mode":"coupled","policy":"shared_coordinate_opportunities_v1","seeds":[states[0]],"option_bank":[],"insertion_outcomes":[event(empty,.25),event(full,.75)],"provenance":"p"},"rates":{"dna_insertion":1.,"dna_deletion":1.,"exon_inactivation":1.,"exonization":1.},"rate_provenance":"r","root":{"entries":[{"state":states[0],"weight":1.}],"provenance":"root"},"observation":{"mode":"constraints","provenance":"declared","view":"annotation","missing_policy":"uninformative","complete_species":[]}}]}
        out=_read(record,c,tree); events=[e for e in out.units[0].process.events if e.kind=="dna_insertion"]
        self.assertEqual(sorted(e.weight for e in events), [.25,.75]); self.assertEqual(out.units[0].root_prior.tolist()[0],1.); np.testing.assert_array_equal(out.units[0].tips["A"],np.ones(len(out.units[0].process.states)))

    def test_json_nonfinite_literals_rejected(self):
        c, seed, tree = _fixture(); process=deletion_process(c,(seed,),provenance="p")
        record=_record(c,process,range(len(process.states)),{"dna_deletion":1.},{"mode":"explicit","provenance":"fixed","tips":{"A":[1.,0.]}})
        for literal in ("NaN","Infinity","1e999"):
            text=json.dumps(record).replace("1.0",literal,1)
            with tempfile.TemporaryDirectory() as d:
                p=Path(d)/"bad.json"; p.write_text(text)
                with self.assertRaises(ValueError): read_repertoire_model(p,(c,),tree)


if __name__ == "__main__": unittest.main()
