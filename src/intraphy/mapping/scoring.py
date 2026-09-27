"""Spawn-safe worker state for independent DNA occurrence-pair scoring."""
from __future__ import annotations


_WORKER_CONFIG = None


def initialize_pair_scoring_worker(occurrences, seqs, context, config):
    """Install read-only scoring inputs once in a spawned worker."""
    global _WORKER_CONFIG
    import os

    for name in ("http_proxy", "https_proxy", "all_proxy", "HTTP_PROXY", "HTTPS_PROXY", "ALL_PROXY"):
        os.environ.pop(name, None)
    _WORKER_CONFIG = (tuple(occurrences), seqs, context, config)


def score_pair_worker(item):
    """Score one DNA pair using worker-local immutable sequence/context data."""
    if _WORKER_CONFIG is None:
        raise RuntimeError("pair-scoring worker was not initialized")
    from intraphy.mapping.match_context import cheap_match_evidence, should_align_pair
    from intraphy.mapping.pairwise_matches import match_evidence

    occurrences, seqs, context, config = _WORKER_CONFIG
    idx, left_index, right_index = item
    left, right = occurrences[left_index], occurrences[right_index]
    should_align, prefilter_status = should_align_pair(
        left, right, seqs, config["min_size_ratio"],
    )
    if should_align:
        evidence = match_evidence(
            left, right, seqs, context,
            aligner=config["aligner"], threads=1,
            context_aligner=config["context_aligner"],
            short_context_max_length=config["short_context_max_length"],
            short_alignment_max_dp_cells=config["short_alignment_max_dp_cells"],
        )
    else:
        evidence = cheap_match_evidence(left, right, context, alignment_backend=prefilter_status)
    return idx, left_index, right_index, evidence, prefilter_status


def process_scored_pairs(pairs, occurrences, row_index, seqs, context, *, workers,
                         aligner, context_aligner, min_size_ratio,
                         short_context_max_length, short_alignment_max_dp_cells):
    """Yield bounded, ordered process results as original parent row objects."""
    from concurrent.futures import ProcessPoolExecutor
    from itertools import islice
    from multiprocessing import get_context

    config = {
        "aligner": aligner,
        "context_aligner": context_aligner,
        "min_size_ratio": min_size_ratio,
        "short_context_max_length": short_context_max_length,
        "short_alignment_max_dp_cells": short_alignment_max_dp_cells,
    }
    pool = ProcessPoolExecutor(
        max_workers=workers,
        mp_context=get_context("spawn"),
        initializer=initialize_pair_scoring_worker,
        initargs=(occurrences, seqs, context, config),
    )
    iterator = iter(pairs)
    try:
        batch_size = max(8, workers * 4)
        while True:
            batch = list(islice(iterator, batch_size))
            if not batch:
                break
            indexed = [(idx, row_index[left["occurrence_id"]],
                        row_index[right["occurrence_id"]])
                       for idx, left, right in batch]
            for _idx, left_index, right_index, evidence, prefilter_status in pool.map(
                    score_pair_worker, indexed):
                yield (_idx, occurrences[left_index], occurrences[right_index],
                       evidence, prefilter_status)
    finally:
        pool.shutdown(wait=True)
