"""Family-level genomic exon fitting and posterior work scheduling."""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import logging
from time import perf_counter

import numpy as np

from ..structure.edits import EDIT_KINDS
from ..structure.serialization import json_safe
from .configuration_compact import evaluate_compact
from .configuration_model import RateModel, evaluate_model
from .genomic_exon_output import state_rows
from .genomic_exon_rates import fit_family_rate, fixed_family_fit


logger = logging.getLogger("intraphy")


def analyze_family(payload, *, unit_executor=None, unit_workers=1):
    family, units, parameters = payload
    ordered_map = map if unit_executor is None else unit_executor.map
    fit_started = perf_counter()
    logger.info("Family fit started: family=%s mode=%s unit_workers=%d",
                family, parameters["parameter_mode"], unit_workers)
    try:
        if parameters["parameter_mode"] == "fit":
            fit = fit_family_rate(units, parameters["branch_length_mode"],
                parameters["max_origins"], unit_map=ordered_map, family_id=family)
            rate_model = None if not fit.get("converged") or fit.get("mu") is None else RateModel(
                {kind: fit["mu"] for kind in EDIT_KINDS}, scale=1.)
        else:
            rate_model = parameters["rates"]
            fit = fixed_family_fit(units, rate_model, parameters["branch_length_mode"],
                parameters["max_origins"], unit_map=ordered_map, family_id=family)
    except BaseException:
        logger.exception("Family fit failed: family=%s elapsed_seconds=%.6f",
                         family, perf_counter() - fit_started)
        raise
    logger.info("Family fit completed: family=%s elapsed_seconds=%.6f",
                family, perf_counter() - fit_started)

    def evaluate_posterior(item):
        index, unit = item
        unit_id = unit.get("unit_id", f"unit-{index}")
        record = {"family_id": family, "unit_id": unit_id,
                  "status": fit["status"] if rate_model is None else "conditional_on_declared_catalogue",
                  "log_likelihood": None,
                  "ctmc": None, "states": state_rows(unit["space"])}
        if rate_model is None:
            return record
        started = perf_counter()
        logger.info("Family unit posterior started: family=%s unit=%s", family, unit_id)
        try:
            if parameters["expected_edits"]:
                result = evaluate_model(unit["space"], unit["tree"], unit["tips"], rate_model,
                    max_origins=parameters["max_origins"], posterior=True, counts=True,
                    branch_length_mode=parameters["branch_length_mode"])
            else:
                result = evaluate_compact(unit["space"], unit["tree"], unit["tips"], rate_model,
                    max_origins=parameters["max_origins"], counts=False,
                    branch_length_mode=parameters["branch_length_mode"])
        except BaseException:
            logger.exception("Family unit posterior failed: family=%s unit=%s elapsed_seconds=%.6f",
                             family, unit_id, perf_counter() - started)
            raise
        logger.info("Family unit posterior completed: family=%s unit=%s elapsed_seconds=%.6f",
                    family, unit_id, perf_counter() - started)
        record["log_likelihood"] = result["log_likelihood"]
        record["ctmc"] = json_safe(result)
        if not np.isfinite(result["log_likelihood"]):
            record["status"] = "zero_probability_under_supplied_or_fitted_parameters"
        else:
            record["status"] = fit["status"]
        return record

    family_started = fit_started
    try:
        unit_results = list(ordered_map(evaluate_posterior, enumerate(units)))
    except BaseException:
        logger.exception("Family analysis failed: family=%s elapsed_seconds=%.6f",
                         family, perf_counter() - family_started)
        raise
    logger.info("Family analysis completed: family=%s elapsed_seconds=%.6f",
                family, perf_counter() - family_started)
    return {"family_id": family, "fit": fit, "units": unit_results}


def _fit_family_safely(payload, *, unit_executor=None, unit_workers=1):
    try:
        return analyze_family(payload, unit_executor=unit_executor,
                              unit_workers=unit_workers)
    except ValueError as exc:
        if "origin_scenarios_incomplete" not in str(exc):
            raise
        family, units, _ = payload
        return {"family_id": family, "fit": {"status": "origin_scenarios_incomplete",
                    "converged": False, "mu": None, "reason": str(exc)},
                "units": [{"family_id": family, "unit_id": unit["unit_id"],
                           "status": "origin_scenarios_incomplete", "log_likelihood": None,
                           "ctmc": None, "states": state_rows(unit["space"])}
                          for unit in units]}


def run_family_tasks(tasks, threads):
    """Use process workers across families and threads within one family only."""
    effective_family_workers = min(threads, len(tasks)) if tasks else 1
    if len(tasks) == 1:
        effective_unit_workers = min(threads, len(tasks[0][1])) if tasks[0][1] else 1
        if effective_unit_workers > 1:
            with ThreadPoolExecutor(max_workers=effective_unit_workers) as pool:
                results = [_fit_family_safely(tasks[0], unit_executor=pool,
                    unit_workers=effective_unit_workers)]
        else:
            results = [_fit_family_safely(tasks[0])]
    elif effective_family_workers > 1:
        effective_unit_workers = 1
        with ProcessPoolExecutor(max_workers=effective_family_workers) as pool:
            results = list(pool.map(_fit_family_safely, tasks))
    else:
        effective_unit_workers = 1
        results = [_fit_family_safely(task) for task in tasks]
    return results, effective_family_workers, effective_unit_workers
