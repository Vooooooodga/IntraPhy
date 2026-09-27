"""Disk-backed, copy-pair staged candidate resolution."""
from __future__ import annotations

from intraphy.mapping.cluster_aggregation import aggregate_staged_matches
from intraphy.mapping.pair_resolution import (
    reapply_ordered_chains,
    resolve_copy_pair_rows,
)


def _attach_facts(records):
    rows = []
    for _ordinal, row, facts, overlap_pair in records:
        row["_staged_projection_facts"] = tuple(facts)
        if overlap_pair is not None:
            row["_staged_overlap_pair"] = overlap_pair
        rows.append(row)
    return rows


def _detach_facts(records, rows):
    packed = []
    for (ordinal, _old_row, old_facts, old_overlap), row in zip(records, rows):
        facts = tuple(row.pop("_staged_projection_facts", old_facts))
        overlap_pair = row.pop("_staged_overlap_pair", old_overlap)
        packed.append((ordinal, row, facts, overlap_pair))
    return packed


def resolve_staged_pair_groups(
    store, occurrence_by_id, occurrences, *, transcript_paths, seqs, gene_loci,
    short_alignment_max_dp_cells, progress=None,
):
    """Resolve all copy-pair groups, preserving the legacy global second pass."""
    total_groups = store.pair_count()
    if progress is not None:
        progress.start(
            "unordered_copy_pair_resolution",
            total_phase_items=total_groups,
            processed_occurrence_pairs=progress.total_occurrence_pairs,
            total_copy_groups=total_groups,
        )
    changed_anywhere = False
    processed_groups = 0
    for pair_key in store.pair_keys():
        records = store.pair_rows(pair_key)
        rows = _attach_facts(records)
        rows, changed = resolve_copy_pair_rows(
            rows,
            occurrence_by_id,
            occurrences,
            transcript_paths=transcript_paths,
            seqs=seqs,
            gene_loci=gene_loci,
            short_alignment_max_dp_cells=short_alignment_max_dp_cells,
        )
        changed_anywhere = changed_anywhere or changed
        store.replace_pair_rows(pair_key, _detach_facts(records, rows))
        processed_groups += 1
        if progress is not None:
            progress.update(
                processed_groups,
                total_groups,
                processed_occurrence_pairs=progress.total_occurrence_pairs,
                processed_copy_groups=processed_groups,
            )
    if progress is not None:
        progress.complete(
            processed_occurrence_pairs=progress.total_occurrence_pairs,
            processed_copy_groups=processed_groups,
        )

    if changed_anywhere:
        if progress is not None:
            progress.start(
                "global_second_chain_pass",
                total_phase_items=total_groups,
                processed_occurrence_pairs=progress.total_occurrence_pairs,
                total_copy_groups=total_groups,
            )
        processed_groups = 0
        for pair_key in store.pair_keys():
            records = store.pair_rows(pair_key)
            rows = [row for _ordinal, row, _facts, _overlap in records]
            reapply_ordered_chains(
                rows, occurrence_by_id, occurrences, transcript_paths=transcript_paths,
            )
            store.replace_pair_rows(pair_key, records)
            processed_groups += 1
            if progress is not None:
                progress.update(
                    processed_groups,
                    total_groups,
                    processed_occurrence_pairs=progress.total_occurrence_pairs,
                    processed_copy_groups=processed_groups,
                )
        if progress is not None:
            progress.complete(
                processed_occurrence_pairs=progress.total_occurrence_pairs,
                processed_copy_groups=processed_groups,
            )


def finish_staged_clustering(
    store, occurrence_by_id, occurrences, species_distances, *, transcript_paths,
    seqs, gene_loci, short_alignment_max_dp_cells, match_writer,
    collect_matches, progress,
):
    store.set_status("resolving")
    resolve_staged_pair_groups(
        store,
        occurrence_by_id,
        occurrences,
        transcript_paths=transcript_paths,
        seqs=seqs,
        gene_loci=gene_loci,
        short_alignment_max_dp_cells=short_alignment_max_dp_cells,
        progress=progress,
    )
    store.set_status("aggregating")
    return aggregate_staged_matches(
        store,
        occurrences,
        occurrence_by_id,
        species_distances,
        match_writer=match_writer,
        collect_matches=collect_matches,
        progress=progress,
    )
