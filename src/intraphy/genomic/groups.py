"""Group physical exon intervals using explicit resolved positional evidence."""
from __future__ import annotations

from collections import defaultdict

from intraphy.elements import EXON_LIKE_ROLES
from intraphy.mapping.member_intervals import _hard_position_eligible
from intraphy.storage.tabular import iter_tsv


def _truth(value):
    return value in {True, 1, "1"}


def family_species_roster(input_dir, occurrences):
    """Return family species and (species, gene-copy) rosters for prepared inputs.

    Explicit targets preserve gene-only species that have no segment occurrences.
    Older prepared cases without input_targets.tsv use the occurrence roster.
    """
    from pathlib import Path

    targets_path = Path(input_dir) / "input_targets.tsv"
    family_species, family_copies = defaultdict(set), defaultdict(set)
    if targets_path.is_file():
        for row in iter_tsv(targets_path, required=("family_id", "species", "gene_copy_id")):
            key = (row["family_id"], row["species"], row["gene_copy_id"])
            if not all(key):
                raise ValueError("input_targets.tsv contains an empty family/species/gene_copy_id")
            family_species[key[0]].add(key[1])
            family_copies[key[0]].add((key[1], key[2]))
        if not family_species:
            raise ValueError("input_targets.tsv contains no selected targets")
        allowed = {(family, species, copy_id) for family, pairs in family_copies.items()
                   for species, copy_id in pairs}
        for row in occurrences.values():
            key = (row["family_id"], row["species"], row["gene_copy_id"])
            if key not in allowed:
                raise ValueError(f"segment occurrence is inconsistent with input_targets.tsv: {key}")
    else:
        for row in occurrences.values():
            family_species[row["family_id"]].add(row["species"])
            family_copies[row["family_id"]].add((row["species"], row["gene_copy_id"]))
    return dict(family_species), dict(family_copies)


def build_physical_sites(input_dir):
    """Return (sites, members, occurrences); coordinates are native GFF coordinates."""
    from pathlib import Path

    occurrences = {}
    physical = {}
    for row in iter_tsv(Path(input_dir) / "segment_occurrences.tsv", required=(
        "occurrence_id", "family_id", "species", "gene_copy_id", "role", "presence_status",
        "contig", "start", "end", "strand", "transcript_id", "source_feature_id",
    )):
        if not row["occurrence_id"]:
            raise ValueError("empty occurrence_id in segment_occurrences.tsv")
        if row["occurrence_id"] in occurrences:
            raise ValueError(f"duplicate occurrence_id: {row['occurrence_id']}")
        occurrences[row["occurrence_id"]] = row
        try:
            start, end = int(row["start"]), int(row["end"])
        except (TypeError, ValueError) as exc:
            raise ValueError(f"invalid occurrence coordinates for {row['occurrence_id']}") from exc
        if start < 1 or end < start:
            raise ValueError(f"invalid occurrence coordinates for {row['occurrence_id']}")
        if row.get("presence_status") != "present" or row.get("role") not in EXON_LIKE_ROLES:
            continue
        try:
            key = (row["family_id"], row["species"], row["contig"], row["strand"],
                   start, end)
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"invalid physical exon interval for occurrence {row.get('occurrence_id')}") from exc
        if key[3] not in {"+", "-"}:
            continue
        if key[4] < 1 or key[5] < key[4]:
            raise ValueError(f"invalid physical exon interval for occurrence {row['occurrence_id']}")
        item = physical.setdefault(key, {"aliases": set(), "occurrence_ids": set(), "gene_copy_ids": set()})
        item["aliases"].add("|".join((row.get("gene_copy_id", "NA"), row.get("transcript_id", "NA"), row.get("source_feature_id", "NA"))))
        item["occurrence_ids"].add(row["occurrence_id"])
        item["gene_copy_ids"].add(row.get("gene_copy_id", "NA"))

    nodes = {}
    node_for_occurrence = {}
    for index, (key, info) in enumerate(sorted(physical.items()), start=1):
        node_id = f"physical_{index:06d}"
        family, species, contig, strand, start, end = key
        nodes[node_id] = {"node_id": node_id, "family_id": family, "species": species,
                          "contig": contig, "strand": strand, "start": start, "end": end,
                          "aliases": sorted(info["aliases"]),
                          "occurrence_ids": sorted(info["occurrence_ids"]),
                          "gene_copy_ids": sorted(info["gene_copy_ids"]), "overlap_ambiguous": False}
        for occurrence_id in info["occurrence_ids"]:
            node_for_occurrence[occurrence_id] = node_id

    by_locus = defaultdict(list)
    for node_id, node in nodes.items():
        by_locus[(node["family_id"], node["species"], node["contig"], node["strand"])].append(node_id)
    for node_ids in by_locus.values():
        ordered = sorted(node_ids, key=lambda n: (nodes[n]["start"], nodes[n]["end"]))
        for i, left_id in enumerate(ordered):
            left = nodes[left_id]
            for right_id in ordered[i + 1:]:
                right = nodes[right_id]
                if right["start"] > left["end"]:
                    break
                if (left["start"], left["end"]) != (right["start"], right["end"]):
                    left["overlap_ambiguous"] = right["overlap_ambiguous"] = True

    adjacency = {node_id: set() for node_id in nodes}
    matches_path = Path(input_dir) / "segment_matches.tsv"
    for match in iter_tsv(matches_path, optional=True):
        # Require both explicit current-schema qualification and the shared hard predicate.
        if not _truth(match.get("position_edge_eligible")) or not _hard_position_eligible(match):
            continue
        if match.get("query_genomic_matched_blocks") in {None, "", "NA"} or match.get("subject_genomic_matched_blocks") in {None, "", "NA"}:
            continue
        left_id = node_for_occurrence.get(match.get("query_occurrence_id"))
        right_id = node_for_occurrence.get(match.get("subject_occurrence_id"))
        if not left_id or not right_id or left_id == right_id:
            continue
        left, right = nodes[left_id], nodes[right_id]
        if left["family_id"] != right["family_id"] or left["species"] == right["species"]:
            continue
        adjacency[left_id].add(right_id)
        adjacency[right_id].add(left_id)

    components = []
    unseen = set(nodes)
    while unseen:
        seed = min(unseen)
        stack, component = [seed], set()
        while stack:
            current = stack.pop()
            if current in component:
                continue
            component.add(current)
            unseen.discard(current)
            stack.extend(adjacency[current] - component)
        components.append(sorted(component))

    sites, members = [], []
    for index, component in enumerate(components, start=1):
        site_id = f"site_{index:06d}"
        component_nodes = [nodes[node_id] for node_id in component]
        species_counts = defaultdict(int)
        for node in component_nodes:
            species_counts[node["species"]] += 1
        reasons = []
        if any(species_counts[s] > 1 for s in species_counts):
            reasons.append("multiple_physical_intervals_in_species")
        if any(node["overlap_ambiguous"] for node in component_nodes):
            reasons.append("overlapping_unequal_physical_intervals")
        eligible = not reasons
        anchor_eligible = eligible and len(species_counts) > 1
        family = component_nodes[0]["family_id"]
        if len(species_counts) == 1 and not reasons:
            reason = "singleton_physical_candidate_requires_two_cross_species_anchor_groups"
        else:
            reason = ";".join(reasons) if reasons else "unique_cross_species_physical_group"
        sites.append({"site_id": site_id, "family_id": family,
                      "member_node_ids": component,
                      "species": sorted(species_counts), "eligible": eligible,
                      "anchor_eligible": anchor_eligible,
                      "reason": reason,
                      "observation_type": "homologous_dna_presence"})
        for node in component_nodes:
            members.append({"site_id": site_id, **node})
    return sites, members, occurrences
