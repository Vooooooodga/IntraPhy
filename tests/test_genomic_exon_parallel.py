"""Contracts for ordered local-unit threading in genomic exon inference."""
from concurrent.futures import ThreadPoolExecutor
from threading import Barrier, Lock, get_ident
import unittest
from unittest.mock import patch

import numpy as np

from intraphy.inference.configuration_fit_cache import CommonRateKernels
from intraphy.inference.configuration_model import RateModel
from intraphy.inference.genomic_exon_family import run_family_tasks
from intraphy.inference.genomic_exon_rates import fit_family_rate, fixed_family_fit
from intraphy.inference.kernel_cache import KernelCache
from intraphy.structure.edits import EDIT_KINDS
from intraphy.structure.space import enumerate_space
from intraphy.structure.types import Catalogue, ExonSpan, Material
from intraphy.topology import SpeciesTree


def _fixture():
    exon = ExonSpan(0, 2)
    space = enumerate_space(Catalogue("f", "u", 2, (exon,), (),
        boundary_candidates=(exon,), observation_unit="genomic_exon_spans"))
    tree = SpeciesTree([{"node_id": "r", "parent_id": "", "label": "r"},
        {"node_id": "a", "parent_id": "r", "label": "A", "branch_length": .5},
        {"node_id": "b", "parent_id": "r", "label": "B", "branch_length": .5}])
    present = np.asarray([float(bool(state.exons)) for state in space.states])
    absent = 1. - present
    same, different = {"A": present, "B": present}, {"A": present, "B": absent}
    units = [{"unit_id": f"u{i}", "space": space, "tree": tree,
              "tips": same if i < 3 else different} for i in range(4)]
    payload = ("f", units, {"parameter_mode": "fit", "rates": None,
        "branch_length_mode": "supplied", "max_origins": None,
        "expected_edits": False})
    return space, tree, units, payload


def _assert_nested_close(test, left, right, tolerance=1e-10):
    if isinstance(left, dict):
        test.assertEqual(set(left), set(right))
        for key in left:
            _assert_nested_close(test, left[key], right[key], tolerance)
    elif isinstance(left, list):
        test.assertEqual(len(left), len(right))
        for a, b in zip(left, right):
            _assert_nested_close(test, a, b, tolerance)
    elif isinstance(left, (float, np.floating)):
        test.assertAlmostEqual(float(left), float(right), delta=tolerance)
    else:
        test.assertEqual(left, right)


class GenomicExonParallelTests(unittest.TestCase):
    def test_threaded_family_matches_serial_fit_posterior_and_unit_order(self):
        _, _, _, payload = _fixture()
        with self.assertLogs("intraphy", level="INFO") as captured:
            serial, serial_family_workers, serial_unit_workers = run_family_tasks([payload], 1)
            threaded, family_workers, unit_workers = run_family_tasks([payload], 2)
        self.assertEqual((serial_family_workers, serial_unit_workers), (1, 1))
        self.assertEqual((family_workers, unit_workers), (1, 2))
        self.assertEqual([unit["unit_id"] for unit in threaded[0]["units"]],
                         ["u0", "u1", "u2", "u3"])
        self.assertEqual(threaded[0]["fit"]["status"], serial[0]["fit"]["status"])
        self.assertAlmostEqual(threaded[0]["fit"]["mu"], serial[0]["fit"]["mu"], delta=1e-8)
        self.assertAlmostEqual(threaded[0]["fit"]["log_likelihood"],
                               serial[0]["fit"]["log_likelihood"], delta=1e-10)
        for threaded_unit, serial_unit in zip(threaded[0]["units"], serial[0]["units"]):
            self.assertAlmostEqual(threaded_unit["log_likelihood"],
                                   serial_unit["log_likelihood"], delta=1e-10)
            _assert_nested_close(self, threaded_unit["ctmc"], serial_unit["ctmc"])
        messages = "\n".join(captured.output)
        self.assertIn("Family fit started", messages)
        self.assertIn("Family fit completed", messages)
        self.assertIn("Family unit likelihood started", messages)
        self.assertIn("Family unit likelihood completed", messages)
        self.assertIn("Family unit posterior started", messages)
        self.assertIn("Family unit posterior completed", messages)
        self.assertIn("Family analysis completed", messages)

    def test_scheduler_uses_one_parallel_layer(self):
        class FakePool:
            def __init__(self, max_workers):
                self.max_workers = max_workers

            def __enter__(self):
                return self

            def __exit__(self, *args):
                return False

            def map(self, function, items):
                self.items = list(items)
                return []

        payload = ("f", [{"unit_id": f"u{i}"} for i in range(4)], {})
        with patch("intraphy.inference.genomic_exon_family.ThreadPoolExecutor",
                   side_effect=FakePool) as threads, \
             patch("intraphy.inference.genomic_exon_family.ProcessPoolExecutor",
                   side_effect=AssertionError("nested process pool")), \
             patch("intraphy.inference.genomic_exon_family.analyze_family",
                   return_value={"family_id": "f"}) as analyze:
            result, family_workers, unit_workers = run_family_tasks([payload], 4)
        self.assertEqual(result, [{"family_id": "f"}])
        self.assertEqual((family_workers, unit_workers), (1, 4))
        self.assertEqual(threads.call_args.kwargs["max_workers"], 4)
        self.assertIsNotNone(analyze.call_args.kwargs["unit_executor"])

        with patch("intraphy.inference.genomic_exon_family.ProcessPoolExecutor",
                   side_effect=FakePool) as processes, \
             patch("intraphy.inference.genomic_exon_family.ThreadPoolExecutor",
                   side_effect=AssertionError("nested thread pool")), \
             patch("intraphy.inference.genomic_exon_family._fit_family_safely",
                   return_value={"family_id": "f"}):
            result, family_workers, unit_workers = run_family_tasks([payload, payload], 4)
        self.assertEqual(result, [])
        self.assertEqual((family_workers, unit_workers), (2, 1))
        self.assertEqual(processes.call_args.kwargs["max_workers"], 2)

    def test_fixed_parameters_and_unknown_tip_keep_compact_posterior(self):
        _, _, units, payload = _fixture()
        unknown_units = []
        for unit in units[:2]:
            tips = dict(unit["tips"])
            tips["B"] = np.ones_like(tips["B"])
            unknown_units.append({**unit, "tips": tips})
        rates = RateModel({kind: .08 for kind in EDIT_KINDS})
        fixed_payload = ("fixed", unknown_units, {"parameter_mode": "fixed", "rates": rates,
            "branch_length_mode": "supplied", "max_origins": None,
            "expected_edits": False})
        serial, _, _ = run_family_tasks([fixed_payload], 1)
        threaded, _, workers = run_family_tasks([fixed_payload], 2)
        self.assertEqual(workers, 2)
        self.assertEqual([row["unit_id"] for row in threaded[0]["units"]], ["u0", "u1"])
        for parallel, single in zip(threaded[0]["units"], serial[0]["units"]):
            self.assertAlmostEqual(parallel["log_likelihood"], single["log_likelihood"], delta=1e-10)
            _assert_nested_close(self, parallel["ctmc"], single["ctmc"])

    def test_fixed_fit_unit_evaluations_overlap_on_two_threads(self):
        _, _, units, _ = _fixture()
        barrier = Barrier(2)
        observed_threads = set()
        observed_lock = Lock()

        def evaluator(*args, **kwargs):
            barrier.wait(timeout=5)
            with observed_lock:
                observed_threads.add(get_ident())
            return {"log_likelihood": -1.}

        with patch("intraphy.inference.genomic_exon_rates.evaluate_model",
                   side_effect=evaluator), ThreadPoolExecutor(max_workers=2) as pool:
            fit = fixed_family_fit(units[:2], RateModel(
                {kind: .1 for kind in EDIT_KINDS}), unit_map=pool.map)
        self.assertEqual(fit["unit_log_likelihoods"], [-1., -1.])
        self.assertEqual(len(observed_threads), 2)

    def test_two_family_process_workers_preserve_order_and_propagate_failure(self):
        _, _, _, payload = _fixture()
        second = ("g", payload[1], payload[2])
        results, family_workers, unit_workers = run_family_tasks([payload, second], 2)
        self.assertEqual((family_workers, unit_workers), (2, 1))
        self.assertEqual([result["family_id"] for result in results], ["f", "g"])
        invalid = ("bad", [], {"parameter_mode": "unknown"})
        with self.assertRaises(KeyError):
            run_family_tasks([payload, invalid], 2)

    def test_origin_limit_remains_explicit_in_parallel_fit(self):
        exon = ExonSpan(0, 2)
        material = Material("m", 2, 3)
        space = enumerate_space(Catalogue("f", "m", 4, (exon,), (),
            material=(material,), boundary_candidates=(exon,),
            observation_unit="genomic_exon_spans"))
        tree = SpeciesTree([{"node_id": "r", "parent_id": "", "label": "r"},
            {"node_id": "a", "parent_id": "r", "label": "A", "branch_length": .5},
            {"node_id": "b", "parent_id": "r", "label": "B", "branch_length": .5}])
        present = np.asarray([float(state.material == (1,)) for state in space.states])
        unit = {"space": space, "tree": tree,
                "tips": {"A": present, "B": 1. - present}}
        unrestricted = fit_family_rate([unit])
        limited = fit_family_rate([unit], max_origins=2)
        self.assertEqual(unrestricted["status"], limited["status"])
        self.assertAlmostEqual(unrestricted["log_likelihood"], limited["log_likelihood"])
        with ThreadPoolExecutor(max_workers=2) as pool:
            with self.assertRaisesRegex(ValueError, "origin_scenarios_incomplete"):
                fit_family_rate([unit], max_origins=1, unit_map=pool.map)

    def test_shared_template_cache_is_guarded_during_concurrent_lookup(self):
        space, tree, units, _ = _fixture()
        lock = Lock()
        cache = KernelCache()

        class LockCheckingCache:
            def get_or_compute(self, key, compute):
                if not lock.locked():
                    raise AssertionError("cache lookup escaped shared lock")
                return cache.get_or_compute(key, compute)

        guarded_cache = LockCheckingCache()
        template = CommonRateKernels(space, guarded_cache, lock)
        templates = [template for _ in units]
        model_rates = {kind: .1 for kind in EDIT_KINDS}
        origins = {}

        def build(item):
            index, template = item
            return template.generator(space, RateModel(model_rates), origins, "a")

        with ThreadPoolExecutor(max_workers=2) as pool:
            matrices = list(pool.map(build, enumerate(templates)))
        for matrix in matrices[1:]:
            np.testing.assert_array_equal(matrix.toarray(), matrices[0].toarray())
        self.assertGreater(cache.hits, 0)


if __name__ == "__main__":
    unittest.main()
