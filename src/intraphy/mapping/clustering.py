"""mapping / clustering: explicit implementation ownership."""
from __future__ import annotations

from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from intraphy.coding_correspondence import CodingProjectionIndex
from intraphy.coding.copy_pair import expected_all_copy_pair_keys
from intraphy.aligners.short import _validate_max_dp_cells
from intraphy.aligners.types import MAX_INTERNAL_DP_CELLS
from intraphy.mapping.alternative import alternative_overlap_evidence
from intraphy.mapping.match_context import cheap_match_evidence
from intraphy.mapping.match_context import copy_order_context
from intraphy.mapping.match_context import load_distance_table
from intraphy.mapping.match_context import should_align_pair
from intraphy.mapping.match_records import _match_row
from intraphy.mapping.match_records import _projection_compatibility
from intraphy.mapping.match_records import _projection_record
from intraphy.mapping.ordered_chains import _apply_ordered_candidate_chains
from intraphy.mapping.pairwise_matches import match_evidence
from intraphy.mapping.policies import graph_components
from intraphy.mapping.policies import inferred_source_label
from intraphy.mapping.policies import known_source_labels
from intraphy.mapping.policies import occurrence_copy_key
from intraphy.mapping.short_alignment import _rerun_anchor_bounded_short_candidates
from intraphy.mapping.scoring import process_scored_pairs
from intraphy.mapping.match_staging import MatchStagingStore
from intraphy.mapping.match_staging import canonical_copy_pair_key
from intraphy.mapping.progress import ProgressReporter
from intraphy.mapping.progress import print_progress_event
from intraphy.mapping.staged_clustering import finish_staged_clustering
from intraphy.mapping.pair_resolution import resolve_joint_terminal_partitions
from intraphy.mapping.cluster_aggregation import fold_global_match_facts
from intraphy.mapping.score_assembly import assemble_scored_match
from intraphy.preparation.annotation_index import parse_attributes
from itertools import islice


def _cluster_segments_impl(occurrences, seqs, identity_threshold=0.7, distance_table=None, aligner="mafft", threads=1, min_size_ratio=0.25, species_distances=None, match_writer=None, context_aligner="minimap2", transcript_paths=None, raw_features=None, coding_msa_mode="linsi", short_context_max_length=300, gene_loci=None, *, protein_evidence_provider=None, pair_scoring_executor="thread", short_alignment_max_dp_cells=MAX_INTERNAL_DP_CELLS, _staging_store=None, progress_callback=None, collect_matches=True):
    short_alignment_max_dp_cells = _validate_max_dp_cells(short_alignment_max_dp_cells)
    if pair_scoring_executor not in {"thread", "process"}:
        raise ValueError("pair_scoring_executor must be 'thread' or 'process'")
    if not isinstance(collect_matches, bool):
        raise TypeError("collect_matches must be Boolean")
    if protein_evidence_provider is not None:
        expected_protein_keys = frozenset(expected_all_copy_pair_keys(occurrences))
        provider_keys = getattr(protein_evidence_provider, "keys", None)
        if provider_keys is None or frozenset(provider_keys) != expected_protein_keys:
            raise ValueError("protein evidence provider does not cover the exact expected occurrence-pair set")
    context = copy_order_context(occurrences, transcript_paths)
    distance_lookup = load_distance_table(distance_table)
    occurrence_by_id = {row["occurrence_id"]: row for row in occurrences}
    occurrence_row_index = {row["occurrence_id"]: index for index, row in enumerate(occurrences)}
    matches = []
    pending_match_rows = []
    accepted_edges = []
    projection_by_occ_ref = {}
    score_by_occ = defaultdict(list)
    source_support = defaultdict(lambda: defaultdict(float))
    genomic_overlap_compatible = set()
    protein_index = None
    if protein_evidence_provider is None and aligner in {"mafft", "auto"} and transcript_paths:
        parsed_features = [
            {**row, "attrs": parse_attributes(row["attrs"]) if isinstance(row.get("attrs"), str) else row.get("attrs", {})}
            for row in raw_features or []
        ]
        protein_index = CodingProjectionIndex(
            occurrences,
            seqs,
            transcript_paths,
            parsed_features,
            threads=threads,
            msa_mode=coding_msa_mode,
        )

    def emit_match(row, *, projection_facts=(), overlap_pair=None):
        if _staging_store is None:
            if not projection_facts:
                left_id = row["query_occurrence_id"]
                right_id = row["subject_occurrence_id"]
                projection_facts = (
                    ((left_id, right_id), projection_by_occ_ref.get((left_id, right_id))),
                    ((right_id, left_id), projection_by_occ_ref.get((right_id, left_id))),
                )
            row["_staged_projection_facts"] = tuple(projection_facts)
            pending_match_rows.append(row)
            return
        left = occurrence_by_id[row["query_occurrence_id"]]
        right = occurrence_by_id[row["subject_occurrence_id"]]
        pair_key = canonical_copy_pair_key(occurrence_copy_key(left), occurrence_copy_key(right))
        ordinal = int(str(row["match_id"]).split("_")[-1])
        _staging_store.add(ordinal, pair_key, row, projection_facts, overlap_pair)

    def iter_pairs():
        index = 0
        by_family = defaultdict(list)
        for occurrence in occurrences:
            by_family[occurrence["family_id"]].append(occurrence)
        for family_rows in by_family.values():
            for left_index, left in enumerate(family_rows):
                for right in family_rows[left_index + 1 :]:
                    if occurrence_copy_key(left) == occurrence_copy_key(right):
                        continue
                    if left.get("role") == "intron" and right.get("role") == "intron":
                        continue
                    yield index, left, right
                    index += 1

    progress = None
    scored_total = None
    if _staging_store is not None or progress_callback is not None:
        scored_total = sum(1 for _ in iter_pairs())
        callback = progress_callback
        if callback is None and _staging_store is not None:
            callback = print_progress_event
        progress = ProgressReporter(callback, scored_total)

    match_count = 0
    by_copy = defaultdict(list)
    for occurrence in occurrences:
        by_copy[occurrence_copy_key(occurrence)].append(occurrence)
    for copy_rows in by_copy.values():
        for left_index, left in enumerate(copy_rows):
            for right in copy_rows[left_index + 1 :]:
                evidence = alternative_overlap_evidence(left, right, context)
                if evidence is None:
                    continue
                score = evidence["sequence_score"]
                pair = frozenset((left["occurrence_id"], right["occurrence_id"]))
                if _staging_store is None:
                    accepted_edges.append((left["occurrence_id"], right["occurrence_id"], score))
                    genomic_overlap_compatible.add(pair)
                    projection_by_occ_ref[(left["occurrence_id"], right["occurrence_id"])] = _projection_record(evidence, "target")
                    projection_by_occ_ref[(right["occurrence_id"], left["occurrence_id"])] = _projection_record(evidence, "query")
                    score_by_occ[left["occurrence_id"]].append(score)
                    score_by_occ[right["occurrence_id"]].append(score)
                match_count += 1
                facts = (
                    ((left["occurrence_id"], right["occurrence_id"]), _projection_record(evidence, "target")),
                    ((right["occurrence_id"], left["occurrence_id"]), _projection_record(evidence, "query")),
                )
                row = _match_row(left, right, evidence, score, 1.0, "same_copy_alternative_overlap", "mapped", match_count)
                emit_match(
                    row,
                    projection_facts=facts,
                    overlap_pair=(left["occurrence_id"], right["occurrence_id"]),
                )

    def score_pair(item):
        idx, left, right = item
        should_align, prefilter_status = should_align_pair(left, right, seqs, min_size_ratio)
        if should_align:
            evidence = match_evidence(
                left,
                right,
                seqs,
                context,
                aligner=aligner,
                threads=1,
                context_aligner=context_aligner,
                short_context_max_length=short_context_max_length,
                short_alignment_max_dp_cells=short_alignment_max_dp_cells,
            )
        else:
            evidence = cheap_match_evidence(left, right, context, alignment_backend=prefilter_status)
        return idx, left, right, evidence, prefilter_status

    worker_count = max(1, int(threads or 1))
    if pair_scoring_executor == "process":
        pool = None
        scored_pairs = process_scored_pairs(
            iter_pairs(), occurrences, occurrence_row_index, seqs, context,
            workers=worker_count, aligner=aligner, context_aligner=context_aligner,
            min_size_ratio=min_size_ratio,
            short_context_max_length=short_context_max_length,
            short_alignment_max_dp_cells=short_alignment_max_dp_cells,
        )
    elif worker_count > 1:
        pool = ThreadPoolExecutor(max_workers=worker_count)
        pair_iterator = iter(iter_pairs())

        def bounded_scores():
            batch_size = max(8, worker_count * 4)
            while True:
                batch = list(islice(pair_iterator, batch_size))
                if not batch:
                    break
                yield from pool.map(score_pair, batch)

        scored_pairs = bounded_scores()
    else:
        pool = None
        scored_pairs = map(score_pair, iter_pairs())

    if progress is not None:
        progress.start("occurrence_pair_scoring", total_phase_items=scored_total)
    scored_processed = 0

    try:
        for _idx, left, right, evidence, prefilter_status in scored_pairs:
            score, status, mapped, row, facts = assemble_scored_match(
                left,
                right,
                evidence,
                prefilter_status,
                identity_threshold,
                distance_lookup,
                protein_index=protein_index,
                protein_evidence_provider=protein_evidence_provider,
                aligner=aligner,
                match_index=match_count + 1,
            )
            if mapped:
                if _staging_store is None:
                    accepted_edges.append((left["occurrence_id"], right["occurrence_id"], score))
                    projection_by_occ_ref[(left["occurrence_id"], right["occurrence_id"])] = _projection_record(evidence, "target")
                    projection_by_occ_ref[(right["occurrence_id"], left["occurrence_id"])] = _projection_record(evidence, "query")
                    score_by_occ[left["occurrence_id"]].append(score)
                    score_by_occ[right["occurrence_id"]].append(score)
                    left_sources = known_source_labels(left)
                    right_sources = known_source_labels(right)
                    if left_sources and not right_sources:
                        for source in left_sources:
                            source_support[right["occurrence_id"]][source] += score
                    if right_sources and not left_sources:
                        for source in right_sources:
                            source_support[left["occurrence_id"]][source] += score
            match_count += 1
            emit_match(row, projection_facts=facts)
            scored_processed += 1
            if progress is not None:
                progress.update(scored_processed, scored_total, processed_occurrence_pairs=scored_processed)
    finally:
        close_scored = getattr(scored_pairs, "close", None)
        try:
            if close_scored is not None:
                close_scored()
        finally:
            if pool is not None:
                pool.shutdown(wait=True)

    if progress is not None:
        progress.complete(processed_occurrence_pairs=scored_processed)

    if _staging_store is not None:
        homology, staged_matches = finish_staged_clustering(
            _staging_store,
            occurrence_by_id,
            occurrences,
            species_distances,
            transcript_paths=transcript_paths,
            seqs=seqs,
            gene_loci=gene_loci,
            short_alignment_max_dp_cells=short_alignment_max_dp_cells,
            match_writer=match_writer,
            collect_matches=collect_matches,
            progress=progress,
        )
        return homology, staged_matches

    resolve_joint_terminal_partitions(pending_match_rows, occurrence_by_id)
    for row in pending_match_rows:
        projection_by_occ_ref.update(dict(row.get("_staged_projection_facts", ())))
    _apply_ordered_candidate_chains(
        pending_match_rows,
        occurrence_by_id,
        occurrences,
        transcript_paths=transcript_paths,
    )
    if gene_loci and _rerun_anchor_bounded_short_candidates(
        pending_match_rows,
        occurrence_by_id,
        seqs,
        gene_loci,
        transcript_paths=transcript_paths,
        short_alignment_max_dp_cells=short_alignment_max_dp_cells,
    ):
        _apply_ordered_candidate_chains(
            pending_match_rows,
            occurrence_by_id,
            occurrences,
            transcript_paths=transcript_paths,
        )
    position_pairs = {
        frozenset((row["query_occurrence_id"], row["subject_occurrence_id"]))
        for row in pending_match_rows
        if row.get("_position_edge_eligible") or (
            "annotated_CDS_protein" in str(row.get("correspondence_basis", ""))
            and row.get("protein_hard_observation_eligible") in {1, "1", True}
            and row.get("protein_mapping_status") in {
                "resolved_local", "resolved_joint_terminal_partition",
            }
        )
    }
    (
        accepted_edges, score_by_occ, source_support,
        projection_by_occ_ref, _staged_overlap_compatible,
    ) = fold_global_match_facts(
        ((row, row.get("_staged_projection_facts", ()), None) for row in pending_match_rows),
        occurrence_by_id,
        position_pairs,
    )
    for row in pending_match_rows:
        public_row = {key: value for key, value in row.items() if not key.startswith("_")}
        if match_writer is not None:
            match_writer(public_row)
        if collect_matches and public_row.get("match_status") == "mapped":
            matches.append(public_row)
    nodes = [row["occurrence_id"] for row in occurrences]
    same_copy_compatible, cross_copy_compatible = _projection_compatibility(projection_by_occ_ref, occurrence_by_id)
    same_copy_compatible |= genomic_overlap_compatible
    components = graph_components(
        nodes,
        accepted_edges,
        occurrence_by_id,
        species_distances,
        same_copy_compatible,
        cross_copy_compatible,
    )
    homology = []
    for idx, occ_ids in enumerate(sorted(components, key=lambda vals: vals[0]), start=1):
        component_id = f"HC_{idx:04d}"
        for occ_id in occ_ids:
            scores = score_by_occ.get(occ_id, [])
            confidence = sum(scores) / len(scores) if scores else 0.5
            homology.append(
                {
                    "homology_id": component_id,
                    "occurrence_id": occ_id,
                    "support_type": "ordered_sequence_correspondence_graph",
                    "confidence": f"{confidence:.6g}",
                    "source_label": inferred_source_label(occurrence_by_id.get(occ_id, {}), source_support.get(occ_id, {})),
                }
            )
    return homology, matches if collect_matches else None


def cluster_segments(
    occurrences, seqs, identity_threshold=0.7, distance_table=None, aligner="mafft",
    threads=1, min_size_ratio=0.25, species_distances=None, match_writer=None,
    context_aligner="minimap2", transcript_paths=None, raw_features=None,
    coding_msa_mode="linsi", short_context_max_length=300, gene_loci=None, *,
    protein_evidence_provider=None, pair_scoring_executor="thread",
    short_alignment_max_dp_cells=MAX_INTERNAL_DP_CELLS, match_staging_path=None,
    progress_callback=None, collect_matches=True,
):
    """Cluster occurrences, optionally staging exact match rows in SQLite.

    ``match_staging_path`` is a fresh SQLite filename owned by this task. Its
    database is retained with a complete/incomplete status for auditability.
    ``collect_matches=False`` suppresses the returned mapped-row list; writers
    continue receiving every row in original order.
    """
    if match_staging_path is None:
        return _cluster_segments_impl(
            occurrences, seqs, identity_threshold, distance_table, aligner,
            threads, min_size_ratio, species_distances, match_writer,
            context_aligner, transcript_paths, raw_features, coding_msa_mode,
            short_context_max_length, gene_loci,
            protein_evidence_provider=protein_evidence_provider,
            pair_scoring_executor=pair_scoring_executor,
            short_alignment_max_dp_cells=short_alignment_max_dp_cells,
            progress_callback=progress_callback,
            collect_matches=collect_matches,
        )
    store = MatchStagingStore(match_staging_path)
    try:
        result = _cluster_segments_impl(
            occurrences, seqs, identity_threshold, distance_table, aligner,
            threads, min_size_ratio, species_distances, match_writer,
            context_aligner, transcript_paths, raw_features, coding_msa_mode,
            short_context_max_length, gene_loci,
            protein_evidence_provider=protein_evidence_provider,
            pair_scoring_executor=pair_scoring_executor,
            short_alignment_max_dp_cells=short_alignment_max_dp_cells,
            _staging_store=store,
            progress_callback=progress_callback,
            collect_matches=collect_matches,
        )
        store.set_status("complete")
    except BaseException as error:
        try:
            store.set_status("incomplete", f"{type(error).__name__}: {error}")
        except BaseException:
            pass
        try:
            store.close()
        except BaseException:
            pass
        raise
    else:
        store.close()
    return result
