"""Small fixed-catalogue calibration for the genomic exon-span CTMC."""
from __future__ import annotations

from dataclasses import asdict, replace
import json
import math
from pathlib import Path

import numpy as np

from ..structure.edits import EDIT_KINDS
from ..structure.space import enumerate_space
from ..structure.types import Catalogue, ExonSpan, Material, ObservationEvidence, ExonConfiguration
from ..topology import SpeciesTree
from ..storage.tabular import write_tsv
from ..structure.serialization import json_safe, write_json
from .configuration_model import RateModel
from .exon_rates import InferenceUnit
from .exon_resampling import sample_unit_history
from .genomic_exon_rates import fit_family_rate
from ..verification.genomic_exon_calibration import score_replicate, summarize_arm

TIP_LABELS = ("A", "B", "C", "D")
TREE_ROWS = (
    {"node_id": "root", "parent_id": "", "label": "root", "branch_length": 0},
    {"node_id": "ab", "parent_id": "root", "label": "ab", "branch_length": .7},
    {"node_id": "cd", "parent_id": "root", "label": "cd", "branch_length": .9},
    {"node_id": "A", "parent_id": "ab", "label": "A", "branch_length": .5},
    {"node_id": "B", "parent_id": "ab", "label": "B", "branch_length": .8},
    {"node_id": "C", "parent_id": "cd", "label": "C", "branch_length": .6},
    {"node_id": "D", "parent_id": "cd", "label": "D", "branch_length": 1.1},
)


def _catalogues(scenario):
    if scenario == "geometry":
        specs = (("geom_a", 12, ((0, 12), (0, 5), (7, 12)), ((5, 7),), ()),
                 ("geom_b", 14, ((0, 14), (0, 6), (8, 14)), ((6, 8),), ()))
    elif scenario == "shared-deletion":
        specs = (("delete_a", 12, ((0, 4), (8, 12)), ((4, 8),), ((0, 12),)),
                 ("delete_b", 14, ((0, 5), (9, 14)), ((5, 9),), ((0, 14),)))
    else:
        raise ValueError("scenario must be 'geometry' or 'shared-deletion'")
    catalogues = []
    for unit, length, raw_spans, junctions, raw_material in specs:
        spans = tuple(ExonSpan(*span) for span in raw_spans)
        materials = tuple(Material(f"tract_{i+1}", *span) for i, span in enumerate(raw_material))
        observations = []
        for index, species in enumerate(TIP_LABELS):
            structure = spans[index % len(spans)]
            material_state = (1,) * len(materials)
            observations.append(ObservationEvidence(species,
                (ExonConfiguration((structure,), material_state),),
                material_presence=material_state))
        catalogues.append(Catalogue("calibration_family", unit, length, spans,
            junctions, materials, tuple(observations), discovery="independent_catalogue",
            boundary_candidates=spans, observation_unit="genomic_exon_spans"))
    return tuple(catalogues)


def _visible(state):
    return (tuple((e.start, e.end) for e in state.exons),
            tuple(int(value == 1) for value in state.material))


def _exact_tip(space, configuration):
    visible = _visible(configuration)
    return np.asarray([float(_visible(state) == visible) for state in space.states])


def _unit(catalogue, space, unit_id):
    tip_states = {}
    for index, observation in enumerate(catalogue.observations):
        tip_states[observation.species] = _exact_tip(space, observation.configurations[0])
    return InferenceUnit(catalogue.family, unit_id, space, SpeciesTree(TREE_ROWS), tip_states)


def calibrate_genomic_exons(output_dir, scenario, rate, units, replicates, seed):
    """Run fixed-catalogue simulations and compare fixed and estimated-rate inference."""
    if not math.isfinite(rate) or rate < 0 or type(units) is not int or units < 1:
        raise ValueError("rate must be finite and nonnegative; units must be positive")
    if type(replicates) is not int or replicates < 1:
        raise ValueError("replicates must be a positive integer")
    catalogues = _catalogues(scenario)
    spaces = [enumerate_space(c) for c in catalogues]
    if any(not space.complete for space in spaces):
        raise ValueError("Built-in calibration catalogue did not enumerate completely")
    base_units = tuple(_unit(c, s, c.unit) for c, s in zip(catalogues, spaces))
    truth_model = RateModel({kind: rate for kind in EDIT_KINDS})
    rng = np.random.default_rng(seed)
    replicate_rows, branches = [], []
    for replicate_index in range(1, replicates + 1):
        sampled = []
        latent = []
        for index in range(units):
            source = base_units[index % len(base_units)]
            template = replace(source, unit=f"{source.unit}_draw{index + 1:03d}")
            observed, assigned, origins = sample_unit_history(template, truth_model, rng)
            sampled.append({"family": observed.family, "unit": observed.unit,
                "space": observed.space, "tree": observed.tree, "tips": observed.tips})
            latent.append((observed, assigned, origins))
        fit_error = None
        try:
            fit = fit_family_rate(sampled)
        except (ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError) as exc:
            fit = {"status": "numerical_fit_error", "converged": False, "mu": None}
            fit_error = exc
        mu = fit.get("mu")
        valid_mu = (isinstance(mu, (int, float, np.number)) and not isinstance(mu, bool)
                    and np.isfinite(mu) and mu >= 0)
        fit_status = ("scored" if fit.get("converged") and valid_mu else
                      "nonidentified" if fit.get("status") in {"flat_or_nonidentified", "no_rate_information"}
                      else "failed")
        if fit_status != "scored":
            mu = None
        models = (("fixed", truth_model), ("estimated", RateModel(
            {kind: mu for kind in EDIT_KINDS}) if mu is not None else None))
        replicate = {"replicate": replicate_index, "fit_detail": fit.get("status"),
            "estimated_mu": mu, "fit_converged": fit.get("converged", False),
            "fit_valid_for_resampling": fit.get("valid_for_resampling", mu is not None),
            "fit_record": json_safe(fit),
            "fit_error_type": type(fit_error).__name__ if fit_error else None,
            "fit_error_message": str(fit_error) if fit_error else None}
        replicate_rows.append(replicate)
        for arm, model in models:
            if model is None:
                replicate[f"{arm}_status"] = fit_status if arm == "estimated" else "failed"
                continue
            try:
                scored_rows = score_replicate(arm, replicate["replicate"], sampled, latent, model)
            except (ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError) as exc:
                scored_rows = []
                replicate[f"{arm}_error_type"] = type(exc).__name__
                replicate[f"{arm}_error_message"] = str(exc)
            if scored_rows:
                branches.extend(scored_rows)
                replicate[f"{arm}_status"] = "scored"
            else:
                replicate[f"{arm}_status"] = "failed"
        # Do not keep latent assignments in the fit input or output record.
    summary_by_arm = {arm: summarize_arm(arm, replicate_rows,
        [row for row in branches if row["arm"] == arm], replicates, units,
        sum(bool(row["parent_id"]) for row in TREE_ROWS))
        for arm in ("fixed", "estimated")}
    out = Path(output_dir)
    metadata = {"scenario": scenario, "rate": rate, "units_per_replicate": units,
        "replicates": replicates, "seed": seed, "tree": TREE_ROWS,
        "catalogues": [asdict(c) for c in catalogues],
        "catalogue_discovery": "independent_catalogue",
        "observation_unit": "genomic_exon_spans",
        "generating_model": asdict(truth_model),
        "material_origin_sampling": "one origin per material tract; root weight is origin_root_weight and each other canonical node has weight 1",
        "root_state_sampling": "uniform among states whose material mask is 1 for root-origin tracts and 0 otherwise",
        "calibration_scope": "fixed_catalogue_conditional_calibration",
        "discovery_pipeline_calibrated": False,
        "rate_fit_inputs": "simulated tip observations, state spaces, and fixed tree only",
        "fit_initialization": "fit_family_rate defaults; generated rate is not supplied"}
    metadata["probability_clamp"] = "values within 1e-12 of [0,1] are clipped for roundoff; larger excursions fail the arm"
    out.mkdir(parents=True, exist_ok=True)
    write_json(out / "calibration_metadata.json", metadata)
    write_json(out / "replicates.json", replicate_rows)
    write_json(out / "branches.json", branches)
    write_json(out / "summary.json", summary_by_arm)
    for arm in ("fixed", "estimated"):
        write_tsv(out / f"reliability_{arm}.tsv", summary_by_arm[arm]["reliability_bins"],
                  ["bin", "lower", "upper", "count", "mean_predicted", "actual_frequency"])
    write_tsv(out / "replicates.tsv", replicate_rows,
        ["replicate", "fit_detail", "estimated_mu", "fit_converged", "fit_valid_for_resampling",
         "fit_error_type", "fit_error_message", "fixed_status", "fixed_error_type",
         "fixed_error_message", "estimated_status", "estimated_error_type", "estimated_error_message"])
    branch_table = [{**row,
        "truth_parent_exons": json.dumps(row["truth_parent_exons"]),
        "truth_child_exons": json.dumps(row["truth_child_exons"]),
        "truth_parent_dna_presence": json.dumps(row["truth_parent_dna_presence"]),
        "truth_child_dna_presence": json.dumps(row["truth_child_dna_presence"]),
        "modal_pairs": json.dumps(json_safe(row["modal_pairs"]))} for row in branches]
    write_tsv(out / "branches.tsv", branch_table,
        ["replicate", "arm", "unit_id", "parent_node_id", "child_node_id", "truth_changed",
         "probability_change", "modal_tie_count", "modal_hit", "modal_credit", "brier", "accuracy",
         "truth_parent_exons", "truth_child_exons", "truth_parent_dna_presence",
         "truth_child_dna_presence", "modal_pairs"])
    return summary_by_arm
