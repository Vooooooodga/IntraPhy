"""Validation specific to prepared genomic DNA-presence evidence."""
import json
from pathlib import Path

from ..storage.tabular import read_tsv


def validate_prepared_input(directory, tree):
    """Check prepared files and exact family/tree panel agreement."""
    directory = Path(directory)
    for name in ("gene_loci.tsv", "segment_occurrences.tsv", "segment_homology.tsv", "segment_matches.tsv", "segment_sequences.fasta"):
        if not (directory / name).is_file():
            raise FileNotFoundError(f"Prepared genomic input lacks {name}")
    occurrences = read_tsv(directory / "segment_occurrences.tsv", ["family_id", "species", "gene_copy_id"])
    from ..genomic import family_species_roster
    family_species, _ = family_species_roster(directory, {str(i): row for i, row in enumerate(occurrences)})
    tips = set(tree.leaf_by_label)
    for family, species in family_species.items():
        if species != tips:
            raise ValueError(f"Prepared family {family!r} does not match the supplied tree tips.")


def validate_reusable_evidence(evidence_dir, input_dir, tree, *, min_identity, min_coverage, max_dp_cells):
    """Require staged evidence to match physical membership, tree and provenance."""
    from ..genomic.groups import build_physical_sites
    from ..genomic.survey import _locus_resource_key

    evidence_dir, input_dir = Path(evidence_dir), Path(input_dir)
    site_rows, member_rows, occurrences = build_physical_sites(input_dir)
    from ..genomic import family_species_roster
    family_species, family_copies = family_species_roster(input_dir, occurrences)
    expected_sites = {row["site_id"]: row for row in site_rows}
    expected_members = {row["node_id"]: row for row in member_rows}
    site_for_node = {row["node_id"]: row["site_id"] for row in member_rows}
    sites = read_tsv(evidence_dir / "genomic_sites.tsv", ["site_id", "family_id", "species", "member_node_ids", "member_count", "eligible", "anchor_eligible", "reason", "observation_type"])
    members = read_tsv(evidence_dir / "genomic_members.tsv", ["site_id", "node_id", "family_id", "species", "contig", "start", "end", "strand", "aliases", "occurrence_ids", "overlap_ambiguous", "genome_fasta"])
    observations = read_tsv(evidence_dir / "dna_observations.tsv", [
        "family_id", "site_id", "species", "state", "eligible", "reason", "evidence",
        "discovery_rule", "survey_scope", "observation_type", "survey_genome_fastas", "min_identity",
        "min_coverage", "max_dp_cells",
    ])
    if len(sites) != len(expected_sites) or len({r["site_id"] for r in sites}) != len(sites) or set(r["site_id"] for r in sites) != set(expected_sites):
        raise ValueError("Staged DNA evidence does not match the prepared physical-site catalogue.")
    for row in sites:
        expected = expected_sites[row["site_id"]]
        if (row["family_id"] != expected["family_id"] or
                set(filter(None, row["species"].split(";"))) != set(expected["species"]) or
                set(filter(None, row["member_node_ids"].split(";"))) != set(expected["member_node_ids"]) or
                int(row["member_count"]) != len(expected["member_node_ids"]) or
                (row["eligible"] in {"1", "true", "True"}) != bool(expected["eligible"]) or
                (row["anchor_eligible"] in {"1", "true", "True"}) != bool(expected["anchor_eligible"]) or
                row["reason"] != expected["reason"] or row["observation_type"] != "homologous_dna_presence"):
            raise ValueError("Staged DNA evidence does not match the prepared physical-site catalogue.")
    if len(members) != len(expected_members) or len({r["node_id"] for r in members}) != len(members) or set(r["node_id"] for r in members) != set(expected_members):
        raise ValueError("Staged DNA evidence does not contain the exact prepared physical-member set.")
    loci = read_tsv(input_dir / "gene_loci.tsv", ["species", "gene_copy_id", "genome_fasta"])
    locus_by_copy = {(row["species"], row["gene_copy_id"]): row for row in loci}
    for row in members:
        expected = expected_members[row["node_id"]]
        source_resources = {
            _locus_resource_key(locus_by_copy[(expected["species"], copy)])
            for copy in expected["gene_copy_ids"]
            if (expected["species"], copy) in locus_by_copy
            and locus_by_copy[(expected["species"], copy)].get("contig") == expected["contig"]
            and locus_by_copy[(expected["species"], copy)].get("strand") == expected["strand"]
        }
        source_resources.discard(None)
        expected_fasta = next(iter(source_resources))[0] if len(source_resources) == 1 else "NA"
        try:
            same = (row["site_id"] == site_for_node[row["node_id"]]
                    and row["family_id"] == expected["family_id"] and row["species"] == expected["species"]
                    and row["contig"] == expected["contig"] and int(row["start"]) == expected["start"]
                    and int(row["end"]) == expected["end"] and row["strand"] == expected["strand"]
                    and json.loads(row["aliases"]) == expected["aliases"]
                    and set(filter(None, row["occurrence_ids"].split(";"))) == set(expected["occurrence_ids"])
                    and (row["overlap_ambiguous"] in {"1", "True", "true"}) == bool(expected["overlap_ambiguous"])
                    and (str(Path(row["genome_fasta"]).resolve()) == expected_fasta
                         if row["genome_fasta"] not in {"", "NA"} else expected_fasta == "NA"))
        except (KeyError, TypeError, ValueError):
            same = False
        if not same:
            raise ValueError("Staged DNA evidence physical members or genome paths differ from the prepared case.")
    tips = set(tree.leaf_by_label)
    if any(species != tips for species in family_species.values()):
        raise ValueError("Prepared families do not match the supplied tree tips.")
    expected_keys = {(site["family_id"], site["site_id"], species)
                     for site in site_rows for species in family_species[site["family_id"]]}
    family_fastas = {}
    species_fastas = {}
    for family in family_species:
        per_species = {}
        for species, copy_id in family_copies[family]:
            locus = locus_by_copy.get((species, copy_id))
            if locus and locus.get("genome_fasta") not in {None, "", "NA"}:
                per_species.setdefault(species, set()).add(str(Path(locus["genome_fasta"]).resolve()))
        fastas = set().union(*per_species.values()) if per_species else set()
        family_fastas[family] = json.dumps(sorted(fastas), separators=(",", ":"))
        species_fastas[family] = {species: ";".join(sorted(paths)) or "NA"
                                  for species, paths in per_species.items()}
    actual_keys = set()
    for row in observations:
        key = (row["family_id"], row["site_id"], row["species"])
        actual_keys.add(key)
        site = expected_sites.get(row["site_id"])
        if (site is None or row["family_id"] != site["family_id"]
                or row["observation_type"] != "homologous_dna_presence"
                or (row["eligible"] in {"1", "true", "True"}) != bool(site["eligible"])
                or row["discovery_rule"] != "unique_optimal_flank_consistent_placement;physical_interval_deduplication"
                or row["survey_scope"] != "unique_optimal_placement_within_flank_bounded_sequence_pair"
                or row["survey_genome_fastas"] != family_fastas.get(row["family_id"])
                or row["genome_fasta"] != species_fastas.get(row["family_id"], {}).get(row["species"], "NA")
                or float(row["min_identity"]) != min_identity or float(row["min_coverage"]) != min_coverage
                or int(row["max_dp_cells"]) != max_dp_cells
                or row["eligible"] not in {"0", "1", "true", "false", "True", "False"}
                or row["state"] not in {"0", "1", "unknown"}):
            raise ValueError("Staged DNA observations have incompatible input provenance or survey parameters.")
    if actual_keys != expected_keys or len(observations) != len(expected_keys):
        raise ValueError("Staged DNA observations do not match the prepared families, sites, and tree tips.")
