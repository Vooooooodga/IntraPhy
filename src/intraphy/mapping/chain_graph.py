"""mapping / chain graph: explicit implementation ownership."""
from __future__ import annotations


def _parent_order_compatible(left_order, right_order, left_parent, right_parent):
    # Equal ranks are only allowed for complementary, disjoint subintervals of
    # the same explicitly identified parent. Geometry is checked by the caller.
    return left_order < right_order or (
        left_order == right_order and bool(left_parent) and left_parent == right_parent
    )


def _ordered_before(left, right, left_path=None, right_path=None):
    path_compatible = True
    if left_path is not None or right_path is not None:
        path_compatible = (
            left_path is not None
            and right_path is not None
            and left_path.context == right_path.context
            and _parent_order_compatible(left_path.query_order, right_path.query_order,
                                         left_path.query_parent_id, right_path.query_parent_id)
            and _parent_order_compatible(left_path.target_order, right_path.target_order,
                                         left_path.target_parent_id, right_path.target_parent_id)
        )
    return (
        path_compatible
        and left.relative_strand == right.relative_strand == "+"
        and left.query.end0 <= right.query.start0
        and left.target.end0 <= right.target.start0
    )


def _path_membership_for_context(candidate, context):
    if context is None:
        return None
    return next(
        (membership for membership in candidate.path_memberships if membership.context == context),
        None,
    )


def _top_complete_paths(candidates, candidate_indices, predecessors, terminal_indices, start_ids):
    """Return the two highest-scoring distinct complete paths in one DAG."""

    def preference(item):
        score, predecessor = item
        predecessor_key = predecessor if predecessor is not None else (-1, -1)
        return -score, predecessor_key

    def retain_two(items):
        retained = []
        seen = set()
        for item in items:
            if item[1] in seen:
                continue
            seen.add(item[1])
            retained.append(item)
            retained.sort(key=preference)
            if len(retained) > 2:
                retained.pop()
        return retained

    paths_by_index = {index: [] for index in candidate_indices}
    path_predecessor = {}
    for index in candidate_indices:
        candidate = candidates[index]
        options = []
        if not start_ids or candidate.candidate_id in start_ids:
            options.append((float(candidate.score), None))
        for predecessor in predecessors[index]:
            for score, path_id in paths_by_index[predecessor]:
                options.append((score + float(candidate.score), path_id))
        selected = retain_two(options)
        for rank, (score, predecessor_id) in enumerate(selected):
            path_id = (index, rank)
            path_predecessor[path_id] = predecessor_id
            paths_by_index[index].append((score, path_id))

    complete = []
    for index in terminal_indices:
        complete.extend(paths_by_index[index])
    selected_complete = retain_two(complete)

    def materialize(path_id):
        path = []
        while path_id is not None:
            candidate_index, _rank = path_id
            path.append(candidates[candidate_index].candidate_id)
            path_id = path_predecessor[path_id]
        return tuple(reversed(path))

    return tuple(
        (score, materialize(path_id))
        for score, path_id in selected_complete
    )


def _solve_context(candidates, context, start_ids, end_ids, score_delta, *, prepared=None):
    negative_infinity = float("-inf")
    if prepared is None:
        path_by_index = [
            _path_membership_for_context(candidate, context)
            for candidate in candidates
        ]
        candidate_indices = [
            index
            for index in range(len(candidates))
            if context is None or path_by_index[index] is not None
        ]
        if start_ids and not any(candidates[index].candidate_id in start_ids for index in candidate_indices):
            return None
        if end_ids and not any(candidates[index].candidate_id in end_ids for index in candidate_indices):
            return None

        predecessors = {index: [] for index in candidate_indices}
        successors = {index: [] for index in candidate_indices}
        for offset, left_index in enumerate(candidate_indices):
            left = candidates[left_index]
            left_path = path_by_index[left_index]
            for right_index in candidate_indices[offset + 1 :]:
                right = candidates[right_index]
                right_path = path_by_index[right_index]
                if _ordered_before(left, right, left_path, right_path):
                    predecessors[right_index].append(left_index)
                    successors[left_index].append(right_index)

        negative_infinity = float("-inf")
        forward = {index: negative_infinity for index in candidate_indices}
        for index in candidate_indices:
            candidate = candidates[index]
            starts_here = not start_ids or candidate.candidate_id in start_ids
            previous = max((forward[item] for item in predecessors[index]), default=negative_infinity)
            if starts_here:
                previous = max(0.0, previous)
            if previous != negative_infinity:
                forward[index] = float(candidate.score) + previous

        terminal_indices = [
            index for index in candidate_indices
            if not end_ids or candidates[index].candidate_id in end_ids
        ]
        reachable_terminal = [index for index in terminal_indices if forward[index] != negative_infinity]
        if not reachable_terminal:
            return None
        best_score = max(forward[index] for index in reachable_terminal)

        backward = {index: negative_infinity for index in candidate_indices}
        for index in reversed(candidate_indices):
            candidate = candidates[index]
            ends_here = not end_ids or candidate.candidate_id in end_ids
            following = max((backward[item] for item in successors[index]), default=negative_infinity)
            if ends_here:
                following = max(0.0, following)
            if following != negative_infinity:
                backward[index] = float(candidate.score) + following

    else:
        candidate_indices, predecessors, successors, forward, backward, reachable_terminal, best_score = prepared

    tolerance = 1e-12 * max(1.0, abs(best_score))
    threshold = best_score - float(score_delta) - tolerance
    through = {
        index: forward[index] + backward[index] - float(candidates[index].score)
        if forward[index] != negative_infinity and backward[index] != negative_infinity
        else negative_infinity
        for index in candidate_indices
    }
    retained = {
        candidates[index].candidate_id for index, value in through.items() if value >= threshold
    }
    best = {
        candidates[index].candidate_id
        for index, value in through.items()
        if abs(value - best_score) <= tolerance
    }
    retained_edges = {
        (candidates[left].candidate_id, candidates[right].candidate_id)
        for left in candidate_indices
        for right in successors[left]
        if forward[left] != negative_infinity
        and backward[right] != negative_infinity
        and forward[left] + backward[right] >= threshold
    }
    complete_paths = _top_complete_paths(
        candidates,
        candidate_indices,
        predecessors,
        reachable_terminal,
        start_ids,
    )
    best_path_count = sum(
        1 for score, _path in complete_paths if abs(score - best_score) <= tolerance
    )
    near_path_count = sum(1 for score, _path in complete_paths if score >= threshold)
    # A node is mandatory iff removing it eliminates *all* near-optimal paths.
    # Work on the existing DAG: do not introduce new edges across a removed node.
    shared_near_ids = set()
    for excluded in candidate_indices:
        candidate_id = candidates[excluded].candidate_id
        if candidate_id not in retained:
            continue
        without = {}
        for index in candidate_indices:
            if index == excluded:
                without[index] = negative_infinity
                continue
            previous = max((without.get(item, negative_infinity) for item in predecessors[index]),
                           default=negative_infinity)
            if not start_ids or candidates[index].candidate_id in start_ids:
                previous = max(0.0, previous)
            without[index] = previous + float(candidates[index].score)
        alternative = max((without.get(index, negative_infinity) for index in reachable_terminal),
                          default=negative_infinity)
        if alternative < threshold:
            shared_near_ids.add(candidate_id)
    return {
        "_prepared": (candidate_indices, predecessors, successors, forward, backward, reachable_terminal, best_score),
        "best_score": best_score,
        "retained": retained,
        "best": best,
        "retained_edges": retained_edges,
        "complete_paths": complete_paths,
        "best_path_count_capped": min(2, best_path_count),
        "near_optimal_path_count_capped": min(2, near_path_count),
        "ambiguous_ids": retained - shared_near_ids,
    }
