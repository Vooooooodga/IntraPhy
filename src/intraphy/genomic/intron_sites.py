"""Family MSA boundary collection and physical-placement consensus."""
from __future__ import annotations

from collections import defaultdict
import json

from intraphy.genomic.intron_coordinates import anchor_metrics, coding_boundaries
from intraphy.genomic.intron_coordinates import fasta_locus, genomic_evidence, physical_key
from intraphy.genomic.intron_coordinates import multiple_physical_copy_species
from intraphy.genomic.intron_coordinates import prepared_base_matches


OBSERVATION_TYPE = "intron_position_presence"
ALIGNMENT_SCOPE = "family_protein_msa_conditional_position_not_orthology_certification"
POSITION_FIELDS = (
    "family_id", "site_id", "site_kind", "eligible", "reason", "left_msa_nt0",
    "right_msa_nt0", "species", "physical_intervals", "annotation_aliases",
    "anchor_window", "min_anchor_pairs", "alignment_mode", "alignment_scope",
    "intron_length", "observation_type",
)
OBSERVATION_FIELDS = (
    "family_id", "site_id", "species", "state", "eligible", "reason", "evidence",
    "physical_intervals", "left_msa_nt0", "right_msa_nt0", "annotation_aliases",
    "coding_right_base_offset", "gff_right_phase", "anchor_window", "min_anchor_pairs",
    "alignment_mode", "alignment_scope", "intron_length", "observation_type",
)
ALIGNMENT_FIELDS = (
    "family_id", "record_id", "aligned_protein", "transcript_keys",
    "residue_columns0", "msa_mode", "alignment_scope",
)


def _alias_records(aliases, physical):
    return [{"intron_id": row.get("intron_id", "NA"),
             "transcript_id": row.get("transcript_id", "NA"),
             "gene_copy_id": row.get("gene_copy_id", "NA"),
             "left_feature_id": row.get("left_feature_id", "NA"),
             "right_feature_id": row.get("right_feature_id", "NA"),
             "right_phase": row.get("right_phase", "NA")}
            for row in aliases[physical]]


def _collect_candidates(family, index, intron_rows, aliases, sequences, copies):
    candidates, transcript_errors = defaultdict(list), {}
    if index is None:
        return candidates, transcript_errors
    projection = index.family_projection(family)
    for transcript_key, transcript in index.transcripts.items():
        if transcript.unavailable_reason:
            transcript_errors[transcript_key] = transcript.unavailable_reason
            continue
        transcript_introns = [row for row in intron_rows
            if tuple(row.get(k) for k in ("family_id", "species", "gene_copy_id", "transcript_id")) == transcript_key]
        for left_nt, right_nt, kind, intron, left, right, right_offset, left_nt_base, right_nt_base in coding_boundaries(
            transcript, projection.columns_for(transcript_key), transcript_introns, sequences,
        ):
            if kind == "intron":
                physical = physical_key(intron)
                if prepared_base_matches(left, left_nt_base, sequences) and prepared_base_matches(right, right_nt_base, sequences):
                    dna, reason = genomic_evidence(intron, left, right, fasta_locus(intron, copies),
                                                   left_nt_base, right_nt_base)
                else:
                    dna, reason = None, "prepared_occurrence_dna_disagrees_with_coding_projection"
                if dna:
                    dna.update({"left_occurrence_id": left.occurrence_id,
                        "left_occurrence_offset0": left.occurrence_interval.start0,
                        "right_occurrence_id": right.occurrence_id,
                        "right_occurrence_offset0": right.occurrence_interval.start0})
                item = {"species": intron["species"], "physical": physical,
                    "kind": kind, "evidence": dna, "reason": reason,
                    "transcript_key": transcript_key, "right_offset": right_offset,
                    "phase": intron.get("right_phase", "NA"),
                    "annotation_aliases": _alias_records(aliases, physical)}
            else:
                physical = (family, transcript_key[1], left.contig, left.strand,
                    min(left.genome_interval.start0, right.genome_interval.start0) + 1,
                    max(left.genome_interval.start0, right.genome_interval.start0) + 1)
                item = {"species": transcript_key[1], "physical": physical,
                    "kind": kind, "evidence": None, "reason": "pending_genomic_check",
                    "transcript_key": transcript_key, "right_offset": right_offset,
                    "phase": "NA", "bases": (left, right, left_nt_base, right_nt_base)}
            candidates[(family, left_nt, right_nt)].append(item)
    return candidates, transcript_errors


def _complete_contiguous_evidence(item, family, copies, sequences):
    left, right, left_nt, right_nt = item["bases"]
    pseudo = {"family_id": family, "species": item["species"],
        "gene_copy_id": item["transcript_key"][2], "contig": left.contig,
        "strand": left.strand, "start": min(left.genome_interval.start0, right.genome_interval.start0) + 1,
        "end": max(left.genome_interval.end0, right.genome_interval.end0)}
    if not (prepared_base_matches(left, left_nt, sequences) and prepared_base_matches(right, right_nt, sequences)):
        item["reason"] = "prepared_occurrence_dna_disagrees_with_coding_projection"
        return
    item["evidence"], item["reason"] = genomic_evidence(
        pseudo, left, right, fasta_locus(pseudo, copies), left_nt, right_nt)
    if item["evidence"]:
        item["evidence"].update({"left_occurrence_id": left.occurrence_id,
            "left_occurrence_offset0": left.occurrence_interval.start0,
            "right_occurrence_id": right.occurrence_id,
            "right_occurrence_offset0": right.occurrence_interval.start0})


def family_intron_rows(family, index, intron_rows, aliases, sequences, copies,
                       family_species, family_copies, *, anchor_window, min_anchor_pairs,
                       site_ordinal_start=0):
    """Build alignment, site catalogue and species calls for one family."""
    alignment_rows, position_rows, observations = [], [], []
    projection = index.family_projection(family) if index else None
    if projection:
        for record_id, sequence in sorted(projection.aligned_records.items()):
            alignment_rows.append({"family_id": family, "record_id": record_id,
                "aligned_protein": sequence,
                "transcript_keys": json.dumps(projection.aliases_by_record[record_id], separators=(",", ":")),
                "residue_columns0": json.dumps(projection.residue_columns[record_id], separators=(",", ":")),
                "msa_mode": projection.mode,
                "alignment_scope": ALIGNMENT_SCOPE})
    candidates, transcript_errors = _collect_candidates(family, index, intron_rows, aliases, sequences, copies)
    intron_keys = {key for key, items in candidates.items() if any(item["kind"] == "intron" for item in items)}
    physical_keys = defaultdict(set)
    for key in intron_keys:
        for item in candidates[key]:
            if item["kind"] == "intron":
                physical_keys[item["physical"]].add(key)
    multiple_copy_species = multiple_physical_copy_species(family_copies, family, copies)
    for ordinal, site_key in enumerate(sorted(intron_keys), site_ordinal_start + 1):
        _fam, left_nt, right_nt = site_key
        site_id = f"IP_{family}_{ordinal:06d}"
        by_species = defaultdict(list)
        for item in candidates[site_key]:
            if item["kind"] == "continuous":
                _complete_contiguous_evidence(item, family, copies, sequences)
            key = item["transcript_key"]
            item["anchors"] = anchor_metrics(projection, key, left_nt, right_nt,
                                               anchor_window, min_anchor_pairs) if projection else []
            if not item["anchors"]:
                item["reason"] = "local_amino_acid_flanks_unsupported"
            by_species[item["species"]].append(item)

        conflicts, physical_members = [], set()
        for species, items in by_species.items():
            physical = {item["physical"] for item in items}
            kinds = {item["kind"] for item in items}
            if len(physical) > 1:
                conflicts.append(f"{species}:competing_physical_boundaries")
            if len(kinds) > 1:
                conflicts.append(f"{species}:conflicting_intron_vs_contiguous_paths")
            if species in multiple_copy_species:
                conflicts.append(f"{species}:multiple_gene_copies_in_family")
            physical_members.update(physical)
        if any(len(physical_keys[item["physical"]]) > 1 for item in candidates[site_key]
               if item["kind"] == "intron"):
            conflicts.append("physical_boundary_maps_to_multiple_msa_positions")
        has_anchor = any(item["anchors"] for item in candidates[site_key])
        if not has_anchor:
            conflicts.append("local_amino_acid_flanks_unsupported")
        eligible = not conflicts and bool(by_species)
        aliases_here = [alias for item in candidates[site_key] if item["kind"] == "intron"
                        for alias in item.get("annotation_aliases", [])]
        lengths = sorted({str(int(alias["end"]) - int(alias["start"]) + 1)
                          for item in candidates[site_key] if item["kind"] == "intron"
                          for alias in aliases[item["physical"]]})
        position_rows.append({"family_id": family, "site_id": site_id,
            "site_kind": "annotated_intron_boundary", "eligible": int(eligible),
            "reason": ";".join(conflicts) if conflicts else "conditional_family_msa_boundary",
            "left_msa_nt0": left_nt, "right_msa_nt0": right_nt,
            "species": ";".join(sorted(by_species)),
            "physical_intervals": json.dumps(sorted(physical_members), separators=(",", ":")),
            "annotation_aliases": json.dumps(aliases_here, sort_keys=True, separators=(",", ":")),
            "anchor_window": anchor_window, "min_anchor_pairs": min_anchor_pairs,
            "alignment_mode": "linsi", "alignment_scope": ALIGNMENT_SCOPE,
            "intron_length": ";".join(lengths) or "NA",
            "observation_type": OBSERVATION_TYPE})
        for species in sorted(family_species):
            items = by_species.get(species, [])
            evidence_items = [item for item in items if item.get("evidence") and item["anchors"]]
            evidence = [{key: value for key, value in item.items() if key != "bases"}
                        for item in evidence_items]
            states = {"1" if item["kind"] == "intron" else "0" for item in evidence_items}
            if not eligible:
                state, reason = "unknown", ";".join(conflicts)
            elif len(states) == 1 and len({item["physical"] for item in evidence_items}) == 1:
                state = next(iter(states))
                reason = ("annotated_intron_interval_and_coding_endpoints_verified" if state == "1"
                          else "consecutive_genomic_coding_bases_verified")
            elif not items:
                state, reason = "unknown", "no_coding_boundary_observation"
            elif not evidence:
                state, reason = "unknown", ";".join(sorted({item["reason"] for item in items}))
            else:
                state, reason = "unknown", "boundary_observations_disagree_or_compete"
            intron_items = [item for item in items if item["kind"] == "intron"]
            observations.append({"family_id": family, "site_id": site_id, "species": species,
                "state": state, "eligible": int(eligible), "reason": reason,
                "evidence": json.dumps(evidence, sort_keys=True, separators=(",", ":")),
                "physical_intervals": json.dumps(sorted({item["physical"] for item in items}), separators=(",", ":")),
                "annotation_aliases": json.dumps([alias for item in items for alias in item.get("annotation_aliases", [])],
                                                   sort_keys=True, separators=(",", ":")),
                "left_msa_nt0": left_nt, "right_msa_nt0": right_nt,
                "coding_right_base_offset": ";".join(sorted({str(item["right_offset"]) for item in items})) or "NA",
                "gff_right_phase": ";".join(sorted({str(item["phase"]) for item in items})) or "NA",
                "anchor_window": anchor_window, "min_anchor_pairs": min_anchor_pairs,
                "alignment_mode": "linsi", "alignment_scope": ALIGNMENT_SCOPE,
                "intron_length": ";".join(sorted({str(int(alias["end"]) - int(alias["start"]) + 1)
                    for item in intron_items for alias in aliases[item["physical"]]})) or "NA",
                "observation_type": OBSERVATION_TYPE})
    return alignment_rows, position_rows, observations, transcript_errors, candidates


def unresolved_intron_rows(family, introns, projected_physical, aliases, species,
                           transcript_errors, *, anchor_window, min_anchor_pairs,
                           site_ordinal_start=0):
    positions, observations = [], []
    for ordinal, intron in enumerate(introns, site_ordinal_start + 1):
        physical = physical_key(intron)
        if physical in projected_physical:
            continue
        site_id = f"IP_{family}_UNRESOLVED_{ordinal:06d}"
        reason = transcript_errors.get((family, intron["species"], intron["gene_copy_id"],
                                        intron.get("transcript_id")),
            "annotated_intron_not_mapped_between_adjacent_known_coding_bases")
        alias_rows = _alias_records(aliases, physical)
        length = int(intron["end"]) - int(intron["start"]) + 1
        positions.append({"family_id": family, "site_id": site_id,
            "site_kind": "unresolved_annotated_intron", "eligible": 0, "reason": reason,
            "left_msa_nt0": "NA", "right_msa_nt0": "NA", "species": intron["species"],
            "physical_intervals": json.dumps([physical], separators=(",", ":")),
            "annotation_aliases": json.dumps(alias_rows, sort_keys=True, separators=(",", ":")),
            "anchor_window": anchor_window, "min_anchor_pairs": min_anchor_pairs,
            "alignment_mode": "linsi", "alignment_scope": ALIGNMENT_SCOPE,
            "intron_length": length, "observation_type": OBSERVATION_TYPE})
        for tip in sorted(species):
            own_tip = tip == intron["species"]
            observations.append({"family_id": family, "site_id": site_id, "species": tip,
                "state": "unknown", "eligible": 0, "reason": reason, "evidence": "[]",
                "physical_intervals": json.dumps([physical] if own_tip else [], separators=(",", ":")),
                "annotation_aliases": json.dumps(alias_rows if own_tip else [], sort_keys=True, separators=(",", ":")),
                "left_msa_nt0": "NA", "right_msa_nt0": "NA",
                "coding_right_base_offset": "NA", "gff_right_phase": "NA",
                "anchor_window": anchor_window, "min_anchor_pairs": min_anchor_pairs,
                "alignment_mode": "linsi", "alignment_scope": ALIGNMENT_SCOPE,
                "intron_length": length if own_tip else "NA", "observation_type": OBSERVATION_TYPE})
    return positions, observations
