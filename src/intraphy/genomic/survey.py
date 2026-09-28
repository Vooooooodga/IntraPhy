"""Bounded, two-flank DNA presence observations for physical interval groups."""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from intraphy.aligners.columns import revcomp
from intraphy.aligners.short import anchored_short_alignment
from intraphy.aligners.types import AlignmentBackendError
from intraphy.storage.fasta import read_fasta_interval


def _as_int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _locus_for(node, loci):
    compatible = {}
    for copy_id in node.get("gene_copy_ids", []):
        locus = loci.get((node["species"], copy_id))
        if locus and locus.get("contig") == node["contig"] and locus.get("strand") == node["strand"]:
            lo, hi = _as_int(locus.get("search_start")), _as_int(locus.get("search_end"))
            fasta = locus.get("genome_fasta")
            if lo is None or hi is None or not fasta or fasta == "NA":
                continue
            resource = (str(Path(fasta).resolve()), locus["contig"], locus["strand"], lo, hi)
            compatible[resource] = locus
    return next(iter(compatible.values())) if len(compatible) == 1 else None


def _locus_resource_key(locus):
    lo, hi = _as_int(locus.get("search_start")), _as_int(locus.get("search_end"))
    fasta = locus.get("genome_fasta")
    if lo is None or hi is None or not fasta or fasta == "NA":
        return None
    return str(Path(fasta).resolve()), locus.get("contig"), locus.get("strand"), lo, hi


def _locus_failure_reason(node, loci):
    keys = {_locus_resource_key(loci[(node["species"], copy_id)])
            for copy_id in node.get("gene_copy_ids", [])
            if (node["species"], copy_id) in loci
            and loci[(node["species"], copy_id)].get("contig") == node["contig"]
            and loci[(node["species"], copy_id)].get("strand") == node["strand"]}
    keys.discard(None)
    return "genome_resource_ambiguous" if len(keys) > 1 else "genome_resource_unavailable"


def _span_sequence(left, right, locus):
    if left["contig"] != right["contig"] or left["strand"] != right["strand"]:
        return None, None, "flanks_not_on_same_physical_path"
    if left["strand"] == "+":
        start, end = left["start"], right["end"]
    else:
        start, end = right["start"], left["end"]
    if start > end:
        return None, None, "flanking_sites_not_ordered"
    lo, hi = _as_int(locus.get("search_start")), _as_int(locus.get("search_end"))
    if lo is None or hi is None or start < lo or end > hi:
        return None, None, "flank_bounded_interval_outside_locus_sequence"
    try:
        sequence = read_fasta_interval(locus["genome_fasta"], left["contig"], start, end)
    except SystemExit:
        return None, None, "flank_bounded_interval_sequence_unavailable"
    if any(base not in "ACGT" for base in sequence):
        return None, None, "flank_bounded_interval_contains_non_acgt_bases"
    if left["strand"] == "-":
        sequence = revcomp(sequence)
    return sequence, (start, end), "ok"


def _query_interval(node, bounds, strand):
    start, end = bounds
    if strand == "+":
        return node["start"] - start, node["end"] - start + 1
    return end - node["end"], end - node["start"] + 1


def _interval_stats(query, target, blocks, interval, target_interval=None):
    start, end = interval
    paired = matches = 0
    for block in blocks:
        left = max(start, block.query.start0)
        right = min(end, block.query.end0)
        if left >= right:
            continue
        offset = left - block.query.start0
        target_start = block.target.start0 + offset
        for qpos in range(left, right):
            tpos = target_start + qpos - left
            if target_interval and not (target_interval[0] <= tpos < target_interval[1]):
                continue
            paired += 1
            matches += query[qpos] == target[tpos]
    return paired / max(1, end - start), matches / max(1, paired), paired


def _candidate_state(candidate, query, target, expected, left_q, right_q,
                     left_t, right_t,
                     min_identity, min_coverage):
    qstart, qend = expected
    left_cov, left_id, _ = _interval_stats(query, target, candidate.aligned_blocks, left_q, left_t)
    right_cov, right_id, _ = _interval_stats(query, target, candidate.aligned_blocks, right_q, right_t)
    if min(left_cov, right_cov) < min_coverage or min(left_id, right_id) < min_identity:
        return "unknown", "both_flanking_homologous_intervals_not_supported"

    coverage, identity, _ = _interval_stats(query, target, candidate.aligned_blocks, expected)
    if coverage >= min_coverage and identity >= min_identity:
        return "1", "expected_interval_aligned_with_identity_and_coverage_between_supported_flanks"

    for gap in candidate.gap_blocks:
        if gap.query.start0 > qstart or gap.query.end0 < qend or gap.target.length != 0:
            continue
        before = next((b for b in candidate.aligned_blocks if b.query.end0 == gap.query.start0), None)
        after = next((b for b in candidate.aligned_blocks if b.query.start0 == gap.query.end0), None)
        if before and after and before.target.end0 == after.target.start0:
            return "0", "absent_at_homologous_position"
    return "unknown", "bounded_alignment_does_not_resolve_expected_interval"


def _placement_consensus(alignments, query, target, expected, left_q, right_q,
                        left_t, right_t,
                        min_identity, min_coverage):
    if not alignments.enumeration_complete:
        reason = ("optimal_placement_ambiguous"
                  if alignments.incomplete_reason == "optimal_alignment_limit_reached"
                  else "optimal_alignment_enumeration_incomplete")
        return "unknown", reason, []
    if not alignments.candidates:
        return "unknown", "no_bounded_alignment", []
    states = [_candidate_state(c, query, target, expected, left_q, right_q, left_t, right_t,
                               min_identity, min_coverage)
              for c in alignments.candidates]
    calls = {state for state, _ in states}
    if len(calls) != 1:
        return "unknown", "optimal_alignments_disagree", states
    state = states[0][0]
    return state, states[0][1], states


def _align_bounded(query, target, max_dp_cells, query_id, target_id):
    try:
        return anchored_short_alignment(
            query, target, mode="global", max_alignments=1,
            max_dp_cells=max_dp_cells, query_occurrence_id=query_id,
            target_occurrence_id=target_id,
        ), "ok"
    except AlignmentBackendError as exc:
        if "resource guard" in str(exc):
            return None, "resource_limit"
        raise


def _site_order(sites, member_by_site, species):
    rows = [member_by_site[site_id].get(species) for site_id in sites]
    rows = [row for row in rows if row]
    return sorted(rows, key=lambda n: (
        n["contig"], n["strand"], n["start"] if n["strand"] == "+" else -n["end"],
    ))


def survey_absent_site(site, members, all_sites, member_by_site, loci_by_copy,
                       min_identity, min_coverage, max_dp_cells):
    """Survey one site in species lacking a physical member, from all references."""
    family, site_id = site["family_id"], site["site_id"]
    own = member_by_site[site_id]
    known = {}
    eligible_ids = {s["site_id"] for s in all_sites if s["family_id"] == family and s["anchor_eligible"]}
    ordered_site_ids = [s["site_id"] for s in all_sites if s["site_id"] in eligible_ids]
    for target_species in sorted(loci_by_copy.get("species_by_family", {}).get(family, set())):
        attempts = []
        if target_species in own:
            node = own[target_species]
            locus = _locus_for(node, loci_by_copy["loci"])
            if locus is None:
                known[target_species] = ("unknown", _locus_failure_reason(node, loci_by_copy["loci"]), [])
            else:
                try:
                    seq = read_fasta_interval(locus["genome_fasta"], node["contig"], node["start"], node["end"])
                except SystemExit:
                    seq = ""
                state = "1" if seq and all(b in "ACGT" for b in seq) else "unknown"
                reason = "annotated_physical_interval_sequence" if state == "1" else "annotated_interval_contains_non_acgt_or_unavailable_sequence"
                known[target_species] = (state, reason, [{"interval": [node["start"], node["end"]], "evidence": reason}])
            continue
        ref_calls = []
        for ref_species, ref_node in own.items():
            ref_locus = _locus_for(ref_node, loci_by_copy["loci"])
            if not ref_locus:
                attempts.append({"reference_species": ref_species,
                                 "status": _locus_failure_reason(ref_node, loci_by_copy["loci"])})
                continue
            ordered = _site_order(ordered_site_ids, member_by_site, ref_species)
            ordered = [row for row in ordered
                       if row["contig"] == ref_node["contig"] and row["strand"] == ref_node["strand"]]
            if all(row["node_id"] != ref_node["node_id"] for row in ordered):
                ordered.append(ref_node)
                ordered.sort(key=lambda n: n["start"] if n["strand"] == "+" else -n["end"])
            pos = next((i for i, row in enumerate(ordered) if row["node_id"] == ref_node["node_id"]), None)
            if pos is None or pos == 0 or pos + 1 >= len(ordered):
                continue
            left_ref, right_ref = ordered[pos - 1], ordered[pos + 1]
            target_left, target_right = member_by_site.get(left_ref["site_id"], {}).get(target_species), member_by_site.get(right_ref["site_id"], {}).get(target_species)
            if not target_left or not target_right:
                continue
            target_locus = _locus_for(target_left, loci_by_copy["loci"])
            right_locus = _locus_for(target_right, loci_by_copy["loci"])
            if not target_locus or not right_locus:
                ambiguous = any(_locus_failure_reason(node, loci_by_copy["loci"]) == "genome_resource_ambiguous"
                                for node in (target_left, target_right))
                attempts.append({"reference_species": ref_species,
                                 "status": "genome_resource_ambiguous" if ambiguous else "genome_resource_unavailable"})
                continue
            if _locus_resource_key(right_locus) != _locus_resource_key(target_locus):
                attempts.append({"reference_species": ref_species, "status": "flank_locus_resources_disagree"})
                continue
            if (target_left["contig"] != target_right["contig"] or target_left["strand"] != target_right["strand"]):
                continue
            if target_left["strand"] == "+" and target_left["end"] >= target_right["start"]:
                continue
            if target_left["strand"] == "-" and target_left["start"] <= target_right["end"]:
                continue

            qseq, qbounds, qstatus = _span_sequence(left_ref, right_ref, ref_locus)
            tseq, tbounds, tstatus = _span_sequence(target_left, target_right, target_locus)
            if qstatus != "ok" or tstatus != "ok":
                attempts.append({"reference_species": ref_species, "status": qstatus if qstatus != "ok" else tstatus})
                continue
            exp = _query_interval(ref_node, qbounds, ref_node["strand"])
            left_q = _query_interval(left_ref, qbounds, ref_node["strand"])
            right_q = _query_interval(right_ref, qbounds, ref_node["strand"])
            left_t = _query_interval(target_left, tbounds, target_left["strand"])
            right_t = _query_interval(target_right, tbounds, target_right["strand"])
            alignments, alignment_status = _align_bounded(
                qseq, tseq, max_dp_cells, ref_node["node_id"], f"{target_species}:{site_id}",
            )
            if alignment_status != "ok":
                attempts.append({"reference_species": ref_species, "status": alignment_status})
                continue
            call, reason, states = _placement_consensus(
                alignments, qseq, tseq, exp, left_q, right_q, left_t, right_t,
                min_identity, min_coverage,
            )
            placements = []
            for candidate, (placement_state, placement_reason) in zip(alignments.candidates, states):
                coverage, identity, paired = _interval_stats(qseq, tseq, candidate.aligned_blocks, exp)
                placements.append({"state": placement_state, "reason": placement_reason,
                                   "candidate_id": candidate.candidate_id or "NA",
                                   "cigar": candidate.cigar, "expected_identity": identity,
                                   "expected_coverage": coverage, "expected_paired_bases": paired})
            ref_calls.append((call, reason, {
                "reference_species": ref_species, "target_species": target_species,
                "survey_state": call, "survey_reason": reason,
                "alignment_count_observed": len(alignments.candidates),
                "alignment_incomplete_reason": alignments.incomplete_reason or "NA",
                "placements": placements,
                "alignment_scope": "unique_optimal_placement_within_flank_bounded_sequence_pair",
                "enumeration_complete": bool(alignments.enumeration_complete),
                "reference_genome_fasta": str(Path(ref_locus["genome_fasta"]).resolve()),
                "target_genome_fasta": str(Path(target_locus["genome_fasta"]).resolve()),
                "reference_span": list(qbounds), "target_span": list(tbounds),
                "expected_interval": [ref_node["start"], ref_node["end"]],
                "left_flank_site_id": left_ref["site_id"],
                "right_flank_site_id": right_ref["site_id"],
                "target_left_flank": [target_left["start"], target_left["end"]],
                "target_right_flank": [target_right["start"], target_right["end"]],
            }))
        calls = {call for call, _, _ in ref_calls if call != "unknown"}
        if any(attempt.get("status") == "resource_limit" for attempt in attempts):
            state, reason = "unknown", "resource_limit"
        elif len(calls) == 1 and all(call != "unknown" for call, _, _ in ref_calls):
            state = next(iter(calls))
            reason = "absent_at_homologous_position" if state == "0" else "present_at_homologous_position"
        else:
            state = "unknown"
            if not ref_calls and attempts and all(a["status"] == "resource_limit" for a in attempts):
                reason = "resource_limit"
            else:
                reason = "no_concordant_bounded_reference_survey" if not ref_calls else "reference_surveys_disagree_or_unresolved"
        known[target_species] = (state, reason, [detail for _, _, detail in ref_calls] + attempts)
    return known
