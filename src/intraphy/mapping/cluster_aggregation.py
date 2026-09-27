"""Global, ordinal-stable aggregation of finalized staged match rows."""
from __future__ import annotations

from collections import defaultdict

from intraphy.mapping.match_records import _projection_compatibility
from intraphy.mapping.policies import (
    graph_components,
    inferred_source_label,
    known_source_labels,
    occurrence_copy_key,
)
from intraphy.storage.values import to_float


def _position_pair(row):
    return frozenset((row["query_occurrence_id"], row["subject_occurrence_id"]))


def _protein_hard(row):
    return (
        "annotated_CDS_protein" in str(row.get("correspondence_basis", ""))
        and row.get("protein_hard_observation_eligible") in {1, "1", True}
        and row.get("protein_mapping_status") in {
            "resolved_local", "resolved_joint_terminal_partition",
        }
    )


def fold_global_match_facts(records, occurrence_by_id, position_pairs, *, on_record=None):
    """Fold finalized rows in input order into the legacy graph facts."""
    accepted_edges = []
    score_by_occ = defaultdict(list)
    source_support = defaultdict(lambda: defaultdict(float))
    projection_by_occ_ref = {}
    genomic_overlap_compatible = set()
    for index, (row, projection_facts, overlap_pair) in enumerate(records, start=1):
        if overlap_pair is not None:
            genomic_overlap_compatible.add(frozenset(overlap_pair))
        if _protein_hard(row):
            row["match_status"] = "mapped"
            row["membership_edge_eligible"] = 1
            row["membership_edge_reason"] = "resolved_annotated_CDS_protein_membership"
            row["_membership_edge_eligible"] = True
        if row.get("_membership_edge_eligible"):
            left_id = row["query_occurrence_id"]
            right_id = row["subject_occurrence_id"]
            edge_score = to_float(row.get("correspondence_score"), 0.0)
            accepted_edges.append((left_id, right_id, edge_score))
            score_by_occ[left_id].append(edge_score)
            score_by_occ[right_id].append(edge_score)
            left_sources = known_source_labels(occurrence_by_id.get(left_id, {}))
            right_sources = known_source_labels(occurrence_by_id.get(right_id, {}))
            if left_sources and not right_sources:
                for source in left_sources:
                    source_support[right_id][source] += edge_score
            if right_sources and not left_sources:
                for source in right_sources:
                    source_support[left_id][source] += edge_score
        if _position_pair(row) in position_pairs:
            for key, value in projection_facts:
                projection_by_occ_ref[key] = value
        if on_record is not None:
            on_record(index)
    return (
        accepted_edges, score_by_occ, source_support,
        projection_by_occ_ref, genomic_overlap_compatible,
    )


def aggregate_staged_matches(
    store, occurrences, occurrence_by_id, species_distances, *,
    match_writer=None, collect_matches=True, progress=None,
):
    """Fold cross-pair evidence in legacy ordinal order and build components."""
    position_pairs = set()
    for _ordinal, row, _facts, _overlap in store.iter_rows():
        if row.get("_position_edge_eligible") or (
            "annotated_CDS_protein" in str(row.get("correspondence_basis", ""))
            and row.get("protein_hard_observation_eligible") in {1, "1", True}
            and row.get("protein_mapping_status") in {
                "resolved_local", "resolved_joint_terminal_partition",
            }
        ):
            position_pairs.add(_position_pair(row))

    row_total = store.row_count()
    if progress is not None:
        progress.start("global_aggregation", total_phase_items=row_total,
                       processed_occurrence_pairs=progress.total_occurrence_pairs)
    def stage_row_updates():
        for ordinal, row, projection_facts, overlap_pair in store.iter_rows():
            yield row, projection_facts, overlap_pair
            store.replace_row(ordinal, row, projection_facts, overlap_pair)

    def report_record(processed):
        if progress is not None:
            progress.update(processed, row_total,
                            processed_occurrence_pairs=progress.total_occurrence_pairs)

    (
        accepted_edges, score_by_occ, source_support,
        projection_by_occ_ref, genomic_overlap_compatible,
    ) = fold_global_match_facts(
        stage_row_updates(), occurrence_by_id, position_pairs,
        on_record=report_record,
    )
    if progress is not None:
        progress.complete(processed_occurrence_pairs=progress.total_occurrence_pairs)

    if progress is not None:
        progress.start("global_graph", total_phase_items=len(accepted_edges),
                       processed_occurrence_pairs=progress.total_occurrence_pairs)
    same_copy_compatible, cross_copy_compatible = _projection_compatibility(
        projection_by_occ_ref, occurrence_by_id,
    )
    same_copy_compatible |= genomic_overlap_compatible
    nodes = [row["occurrence_id"] for row in occurrences]
    components = graph_components(
        nodes,
        accepted_edges,
        occurrence_by_id,
        species_distances,
        same_copy_compatible,
        cross_copy_compatible,
    )
    if progress is not None:
        progress.complete(processed_occurrence_pairs=progress.total_occurrence_pairs)

    homology = []
    for idx, occ_ids in enumerate(sorted(components, key=lambda vals: vals[0]), start=1):
        component_id = f"HC_{idx:04d}"
        for occ_id in occ_ids:
            scores = score_by_occ.get(occ_id, [])
            confidence = sum(scores) / len(scores) if scores else 0.5
            homology.append({
                "homology_id": component_id,
                "occurrence_id": occ_id,
                "support_type": "ordered_sequence_correspondence_graph",
                "confidence": f"{confidence:.6g}",
                "source_label": inferred_source_label(
                    occurrence_by_id.get(occ_id, {}), source_support.get(occ_id, {}),
                ),
            })

    if progress is not None:
        progress.start("ordinal_output", total_phase_items=row_total,
                       processed_occurrence_pairs=progress.total_occurrence_pairs)
    matches = [] if collect_matches else None
    processed = 0
    for _ordinal, row, _facts, _overlap in store.iter_rows():
        public_row = {key: value for key, value in row.items() if not key.startswith("_")}
        if match_writer is not None:
            match_writer(public_row)
        if collect_matches and public_row.get("match_status") == "mapped":
            matches.append(public_row)
        processed += 1
        if progress is not None:
            progress.update(processed, row_total,
                            processed_occurrence_pairs=progress.total_occurrence_pairs)
    if progress is not None:
        progress.complete(processed_occurrence_pairs=progress.total_occurrence_pairs)
    return homology, matches
