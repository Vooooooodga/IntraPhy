import unittest
import math
import json
import numpy as np
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
from intraphy.commands.parser import build_parser
from intraphy.inference.repertoire_rates import RepertoireInferenceUnit
from intraphy.inference.repertoire_inputs import RepertoireBundle
from intraphy.structure.repertoire_process import RepertoireProcess, RepertoireEvent
from intraphy.structure.repertoire import ExonRepertoire
from intraphy.structure.types import Catalogue, ExonConfiguration, ExonSpan, Material
from intraphy.topology import SpeciesTree

class RepertoireCommandParserTests(unittest.TestCase):
    def _bundle(self, families=("g",), rate_values=None):
        units=[]; rates={}; rprov={}; metadata={}
        for family in families:
            exon=ExonSpan(0,10); cat=Catalogue(family,"u",10,(exon,),(),material=(Material("m",0,10),))
            before=ExonRepertoire((ExonConfiguration((),(1,)),), (1,)); after=ExonRepertoire((ExonConfiguration((),(2,)),), (2,))
            process=RepertoireProcess((before,after),(RepertoireEvent(before,after,"dna_deletion","m"),),"process",cat)
            tree=SpeciesTree([{"node_id":"r","parent_id":"","label":"r"},{"node_id":"a","parent_id":"r","label":"A","branch_length":"1"}])
            units.append(RepertoireInferenceUnit(process,tree,{"A":np.array([1.,0.])},np.array([1.,0.]),"root","observation"))
            key=(family,"u"); rates[key]={"dna_deletion":.7 if rate_values is None else rate_values[len(units)-1]}; rprov[key]="rate"; metadata[key]={"A":"native"}
        return RepertoireBundle(tuple(units),rates,rprov,"unit",1.,"fixture",metadata,{"schema":"model/1"})
    def test_infer_requires_explicit_model_inputs(self):
        parser = build_parser()
        args = parser.parse_args(["infer-exon-repertoires", "--exon-configurations", "c.jsonl", "--species-tree", "tree.tsv", "--repertoire-model", "model.json", "--output-dir", "out"])
        self.assertEqual(args.command, "infer-exon-repertoires")
        self.assertFalse(args.expected_edits)

    def test_fit_requires_bounds_and_collection_provenance(self):
        parser = build_parser()
        args = parser.parse_args(["fit-exon-repertoire-rates", "--exon-configurations", "c.jsonl", "--species-tree", "tree.tsv", "--repertoire-model", "model.json", "--output-dir", "out", "--log-scale-bounds", "-2", "2", "--collection-provenance", "independent genes"])
        self.assertEqual(args.log_scale_bounds, [-2., 2.])
        self.assertEqual(args.collection_provenance, "independent genes")

    def test_required_model_inputs_are_enforced(self):
        parser = build_parser()
        with self.assertRaises(SystemExit):
            parser.parse_args(["infer-exon-repertoires", "--output-dir", "out"])

    def test_history_runner_serializes_payload_and_refuses_existing_target(self):
        from intraphy.inference import repertoire_run
        exon = ExonSpan(0, 10); cat = Catalogue("g", "u", 10, (exon,), (), material=(Material("m", 0, 10),))
        before = ExonRepertoire((ExonConfiguration((), (1,)),), (1,)); after = ExonRepertoire((ExonConfiguration((), (2,)),), (2,))
        process = RepertoireProcess((before, after), (RepertoireEvent(before, after, "dna_deletion", "m"),), "process", cat)
        tree = SpeciesTree([{"node_id":"r","parent_id":"","label":"r"},{"node_id":"a","parent_id":"r","label":"A","branch_length":"1"}])
        unit = RepertoireInferenceUnit(process, tree, {"A": [1., 0.]}, np.array([1., 0.]), "root", "observation")
        bundle = self._bundle()
        fake = {"log_likelihood": 0., "scale": 1., "root_prior": np.array([1.,0.]), "root_prior_provenance":"root", "observation_provenance":"observation", "process_provenance":"process", "nodes":{"r":np.array([1.,0.])}, "branches":{("r","a"):{"endpoint_matrix":np.eye(2),"probability_at_least_one_edit":0.,"probability_different_endpoints":0.}}}
        with TemporaryDirectory() as tmp:
            with patch.object(repertoire_run, "write_json") as writer, patch.object(repertoire_run, "evaluate_repertoire_model", return_value=fake):
                repertoire_run.infer_repertoires(bundle, tmp)
                payload = writer.call_args.args[1]
                json.dumps(payload, allow_nan=False)
                unit_record = payload["units"][0]
                self.assertTrue(unit_record["catalogue"] and unit_record["states"] and unit_record["events"])
                self.assertEqual(unit_record["events"][0]["source"], 0)
                self.assertTrue(unit_record["tree"] and unit_record["tips"] and unit_record["observation_metadata"])
                self.assertEqual(payload["model_specification"], {"schema":"model/1"})
                self.assertIsInstance(unit_record["branches"], list)
            Path(tmp, "joint_repertoire_history.json").touch()
            with self.assertRaises(ValueError): repertoire_run.infer_repertoires(bundle, tmp)

    def test_impossible_likelihood_is_explicitly_statused(self):
        from intraphy.inference import repertoire_run
        with TemporaryDirectory() as tmp:
            exon = ExonSpan(0, 10); cat = Catalogue("g", "u", 10, (exon,), (), material=(Material("m",0,10),))
            before = ExonRepertoire((ExonConfiguration((),(1,)),), (1,)); after = ExonRepertoire((ExonConfiguration((),(2,)),), (2,)); process = RepertoireProcess((before,after),(RepertoireEvent(before,after,"dna_deletion","m"),),"p",cat)
            tree = SpeciesTree([{"node_id":"r","parent_id":"","label":"r"},{"node_id":"a","parent_id":"r","label":"A","branch_length":"1"}]); unit = RepertoireInferenceUnit(process,tree,{"A":[0.,1.]},np.array([1.,0.]),"r","o")
            bundle = self._bundle()
            fake = {"log_likelihood":-math.inf,"scale":1.,"root_prior":np.array([1.,0.]),"root_prior_provenance":"r","observation_provenance":"o","process_provenance":"p","nodes":{},"branches":{}}
            with patch.object(repertoire_run, "write_json") as writer, patch.object(repertoire_run, "evaluate_repertoire_model", return_value=fake):
                repertoire_run.infer_repertoires(bundle, tmp)
                record = writer.call_args.args[1]["units"][0]
                self.assertIsNone(record["log_likelihood"]); self.assertEqual(record["log_likelihood_status"], "impossible_observation")

    def test_conflicting_relative_rates_fail_before_fit(self):
        from intraphy.inference import repertoire_run
        bundle = self._bundle(("g1", "g2"), rate_values=(.7, .8))
        with TemporaryDirectory() as tmp, patch.object(repertoire_run, "fit_repertoire_scale") as fit:
            with self.assertRaises(ValueError): repertoire_run.fit_repertoire_rates(bundle, tmp, log_scale_bounds=(-1.,1.), collection_provenance="independent")
            fit.assert_not_called()

    def test_fewer_than_two_genes_status_is_preserved(self):
        from intraphy.inference import repertoire_run
        bundle = self._bundle()
        result = {"status":"fewer_than_two_genes", "converged":False, "scale":None}
        with TemporaryDirectory() as tmp, patch.object(repertoire_run, "fit_repertoire_scale", return_value=result):
            with patch.object(repertoire_run, "write_json") as writer:
                repertoire_run.fit_repertoire_rates(bundle, tmp, log_scale_bounds=(-1.,1.), collection_provenance="one gene")
                self.assertEqual(writer.call_args.args[1]["fit"]["status"], "fewer_than_two_genes")

    def test_invalid_model_fails_before_output_reservation(self):
        from intraphy.commands import preflight
        from intraphy.commands import session
        args = type("Args", (), {"command":"infer-exon-repertoires", "exon_configurations":"c", "species_tree":"t", "repertoire_model":"m", "expected_edits":False, "output_dir":"out"})()
        with patch.object(preflight.Path, "is_file", return_value=True), patch("intraphy.structure.serialization.read_catalogues", return_value=()), patch("intraphy.inference.repertoire_run.load_repertoire_bundle", side_effect=ValueError("invalid model")), patch.object(session, "reserve_output") as reserve:
            with self.assertRaises(ValueError):
                with session.command_session(args): pass
            reserve.assert_not_called()
