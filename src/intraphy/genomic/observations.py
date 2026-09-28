"""Prepare physical homologous-DNA site observations from mapping outputs."""
from __future__ import annotations

from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path

from .groups import build_physical_sites, family_species_roster
from .survey import _locus_for, _locus_resource_key, survey_absent_site
from intraphy.storage.tabular import iter_tsv, write_tsv


OBSERVATION_TYPE = "homologous_dna_presence"
SITE_FIELDS = ("site_id", "family_id", "species", "member_node_ids", "member_count",
               "eligible", "anchor_eligible", "reason", "observation_type")
MEMBER_FIELDS = ("site_id", "node_id", "family_id", "species", "contig", "start", "end",
                 "strand", "genome_fasta", "aliases", "occurrence_ids", "overlap_ambiguous")
OBSERVATION_FIELDS = ("family_id", "site_id", "species", "state", "eligible", "reason",
                      "evidence", "genome_fasta", "survey_genome_fastas",
                      "min_identity", "min_coverage",
                      "max_dp_cells", "discovery_rule", "survey_scope", "observation_type")


_SURVEY_CONTEXT = None


def _initialize_survey_worker(context):
    global _SURVEY_CONTEXT
    _SURVEY_CONTEXT = context


def _survey_site_worker(site):
    return site["site_id"], survey_absent_site(site, *_SURVEY_CONTEXT)


def prepare_dna_observations(input_dir, output_dir, *, min_identity=0.7,
                             min_coverage=0.8, max_dp_cells=250000, threads=1) -> list[dict]:
    """Write DNA material calls at cross-species homologous physical positions.

    Calls describe homologous DNA material only. They do not infer exon use,
    splice state, evolutionary direction, or molecular mechanism.
    """
    if not 0 <= min_identity <= 1 or not 0 <= min_coverage <= 1:
        raise ValueError("identity and coverage thresholds must be within [0, 1]")
    if type(max_dp_cells) is not int or max_dp_cells < 1:
        raise ValueError("max_dp_cells must be a positive integer")
    if type(threads) is not int or threads < 1:
        raise ValueError("threads must be a positive integer")
    input_dir, output_dir = Path(input_dir), Path(output_dir)
    sites, members, occurrences = build_physical_sites(input_dir)
    locus_rows = list(iter_tsv(input_dir / "gene_loci.tsv", required=(
        "species", "gene_copy_id", "contig", "strand", "search_start", "search_end", "genome_fasta",
    )))
    loci = {}
    for row in locus_rows:
        key = (row["species"], row["gene_copy_id"])
        prior = loci.get(key)
        if prior is not None and _locus_resource_key(prior) != _locus_resource_key(row):
            raise ValueError(f"conflicting gene_loci resources/spans for {key[0]}:{key[1]}")
        loci[key] = row
    family_species, family_copies = family_species_roster(input_dir, occurrences)
    fasta_by_family_species = defaultdict(set)
    for family, copies in family_copies.items():
        for species, copy_id in copies:
            locus = loci.get((species, copy_id))
            if locus and locus.get("genome_fasta") not in {None, "", "NA"}:
                fasta_by_family_species[(family, species)].add(
                    str(Path(locus["genome_fasta"]).resolve())
                )
    loci_by_copy = {"loci": loci, "species_by_family": family_species}

    member_by_site = defaultdict(dict)
    repeated_species = set()
    for member in members:
        # Components with repeated species remain in the catalogue but cannot
        # provide a unique member for inference.
        if member["species"] in member_by_site[member["site_id"]]:
            repeated_species.add((member["site_id"], member["species"]))
        else:
            member_by_site[member["site_id"]][member["species"]] = member
    for site_id, species in repeated_species:
        member_by_site[site_id][species] = None

    site_eligible = {site["site_id"]: bool(site["eligible"]) for site in sites}
    calls_by_site = {site["site_id"]: {} for site in sites if not site_eligible[site["site_id"]]}
    survey_sites = [site for site in sites if site_eligible[site["site_id"]]]
    context = (members, sites, member_by_site, loci_by_copy,
               min_identity, min_coverage, max_dp_cells)
    if threads > 1 and len(survey_sites) > 1:
        with ProcessPoolExecutor(max_workers=min(threads, len(survey_sites)),
                                 initializer=_initialize_survey_worker,
                                 initargs=(context,)) as pool:
            calls_by_site.update(dict(pool.map(_survey_site_worker, survey_sites)))
    else:
        calls_by_site.update((site["site_id"], survey_absent_site(site, *context))
                             for site in survey_sites)

    observations = []
    for site in sites:
        for species in sorted(family_species[site["family_id"]]):
            eligible = site_eligible[site["site_id"]]
            member = member_by_site[site["site_id"]].get(species)
            if not eligible:
                state, reason, evidence = "unknown", site["reason"], []
            elif member:
                state, reason, evidence = calls_by_site[site["site_id"]].get(
                    species, ("unknown", "genome_resource_unavailable", []))
            else:
                state, reason, evidence = calls_by_site[site["site_id"]].get(
                    species, ("unknown", "no_concordant_bounded_reference_survey", []))
            observations.append({
                "family_id": site["family_id"], "site_id": site["site_id"],
                "species": species, "state": state, "eligible": int(eligible),
                "reason": reason,
                "evidence": json.dumps(evidence, sort_keys=True, separators=(",", ":")),
                "genome_fasta": ";".join(sorted(fasta_by_family_species[(site["family_id"], species)])) or "NA",
                "survey_genome_fastas": json.dumps(sorted({path for (fam, _sp), paths in fasta_by_family_species.items()
                    if fam == site["family_id"] for path in paths}), separators=(",", ":")),
                "min_identity": min_identity, "min_coverage": min_coverage,
                "max_dp_cells": max_dp_cells,
                "discovery_rule": "unique_optimal_flank_consistent_placement;physical_interval_deduplication",
                "survey_scope": "unique_optimal_placement_within_flank_bounded_sequence_pair",
                "observation_type": OBSERVATION_TYPE,
            })

    site_rows = []
    for site in sites:
        site_rows.append({"site_id": site["site_id"], "family_id": site["family_id"],
                          "species": ";".join(site["species"]),
                          "member_node_ids": ";".join(site["member_node_ids"]),
                          "member_count": len(site["member_node_ids"]),
                          "eligible": int(site["eligible"]),
                          "anchor_eligible": int(site["anchor_eligible"]),
                          "reason": site["reason"], "observation_type": OBSERVATION_TYPE})
    member_rows = [{"site_id": row["site_id"], "node_id": row["node_id"],
                    "family_id": row["family_id"], "species": row["species"],
                    "contig": row["contig"], "start": row["start"], "end": row["end"],
                    "strand": row["strand"],
                    "genome_fasta": (str(Path(resource["genome_fasta"]).resolve())
                        if (resource := _locus_for(row, loci)) else "NA"),
                    "aliases": json.dumps(row["aliases"]),
                    "occurrence_ids": ";".join(row["occurrence_ids"]),
                    "overlap_ambiguous": int(row["overlap_ambiguous"])} for row in members]
    write_tsv(output_dir / "genomic_sites.tsv", site_rows, SITE_FIELDS)
    write_tsv(output_dir / "genomic_members.tsv", member_rows, MEMBER_FIELDS)
    write_tsv(output_dir / "dna_observations.tsv", observations, OBSERVATION_FIELDS)
    return observations
