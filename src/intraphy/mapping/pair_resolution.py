"""Exact within-copy-pair terminal and ordered-candidate resolution."""
from __future__ import annotations

from intraphy.coordinates import parse_legacy_blocks
from intraphy.mapping.match_records import _projection_record
from intraphy.mapping.ordered_chains import _apply_ordered_candidate_chains
from intraphy.mapping.policies import occurrence_copy_key
from intraphy.mapping.short_alignment import _rerun_anchor_bounded_short_candidates


def resolve_joint_terminal_partitions(rows, occurrence_by_id):
    """Apply the legacy joint-terminal partition rule to a complete row set."""
    terminal_groups = {}
    for row in rows:
        if (
            row.get("protein_mapping_status") != "supported_unanchored"
            or row.get("protein_candidate_coordinate_consensus") not in {1, "1", True}
            or row.get("protein_terminal_side") not in {"left", "right"}
            or row.get("protein_projected_blocks") in {None, "", "NA"}
        ):
            continue
        query = occurrence_by_id.get(row.get("query_occurrence_id"), {})
        subject = occurrence_by_id.get(row.get("subject_occurrence_id"), {})
        key = (
            tuple(sorted((occurrence_copy_key(query), occurrence_copy_key(subject)))),
            row.get("protein_best_query_transcript", "NA"),
            row.get("protein_best_target_transcript", "NA"),
        )
        terminal_groups.setdefault(key, []).append(row)

    for group in terminal_groups.values():
        if {row.get("protein_terminal_side") for row in group} != {"left", "right"}:
            continue
        shared_occurrences = set.intersection(*(
            {row.get("query_occurrence_id"), row.get("subject_occurrence_id")}
            for row in group
        ))
        if len(shared_occurrences) != 1:
            continue
        reference_id = next(iter(shared_occurrences))
        reference_intervals = []
        parsed_by_row = {}
        for row in group:
            try:
                blocks = tuple(parse_legacy_blocks(row["protein_projected_blocks"]))
            except (TypeError, ValueError):
                blocks = tuple()
            if not blocks:
                break
            parsed_by_row[id(row)] = blocks
            reference_intervals.extend(
                block.query if row.get("query_occurrence_id") == reference_id else block.target
                for block in blocks
            )
        else:
            ordered = sorted(reference_intervals)
            if any(left.overlaps(right) for left, right in zip(ordered, ordered[1:])):
                continue
            for row in group:
                row["protein_mapping_status"] = "resolved_joint_terminal_partition"
                row["protein_membership_eligible"] = 1
                row["protein_position_eligible"] = 1
                row["protein_hard_observation_eligible"] = 1
                row["match_status"] = "mapped"
                row["correspondence_basis"] = "annotated_CDS_protein"
                row["correspondence_score"] = (
                    0.70 * float(row["protein_aa_identity"])
                    + 0.30 * min(float(row["protein_query_cds_coverage"]),
                                 float(row["protein_target_cds_coverage"]))
                )
                blocks = parsed_by_row[id(row)]
                evidence = dict(row, protein_projected_blocks=blocks)
                left_id = row["query_occurrence_id"]
                right_id = row["subject_occurrence_id"]
                facts = dict(row.get("_staged_projection_facts", ()))
                facts[(left_id, right_id)] = _projection_record(evidence, "target")
                facts[(right_id, left_id)] = _projection_record(evidence, "query")
                row["_staged_projection_facts"] = tuple(facts.items())


def resolve_copy_pair_rows(
    rows, occurrence_by_id, occurrences, *, transcript_paths, seqs, gene_loci,
    short_alignment_max_dp_cells,
):
    """Resolve one complete unordered copy pair while retaining every row."""
    resolve_joint_terminal_partitions(rows, occurrence_by_id)

    _apply_ordered_candidate_chains(
        rows, occurrence_by_id, occurrences, transcript_paths=transcript_paths,
    )
    changed = bool(gene_loci and _rerun_anchor_bounded_short_candidates(
        rows,
        occurrence_by_id,
        seqs,
        gene_loci,
        transcript_paths=transcript_paths,
        short_alignment_max_dp_cells=short_alignment_max_dp_cells,
    ))
    return rows, changed


def reapply_ordered_chains(rows, occurrence_by_id, occurrences, transcript_paths=None):
    """Apply the legacy global second chain pass to one complete pair group."""
    return _apply_ordered_candidate_chains(
        rows, occurrence_by_id, occurrences, transcript_paths=transcript_paths,
    )
