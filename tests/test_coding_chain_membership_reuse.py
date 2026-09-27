"""Regression contracts for per-candidate transcript-context lookup."""

import unittest
from unittest.mock import patch

from intraphy.coordinates import Interval0
from intraphy.mapping import chain_graph
from intraphy.mapping.chain_graph import _solve_context
from intraphy.mapping.chain_types import ChainCandidate, ChainPathMembership
from intraphy.mapping.ordered_paths import ordered_candidate_chain


def membership(context, query_order, target_order):
    query_path, target_path, query_contig, target_contig, query_strand, target_strand = context
    return ChainPathMembership(
        query_path, target_path, query_order, target_order,
        query_contig, target_contig, query_strand, target_strand,
    )


def candidate(candidate_id, start, score, memberships=()):
    interval = Interval0(start, start + 2)
    return ChainCandidate(
        candidate_id, interval, interval, score, "nt_test",
        path_memberships=tuple(memberships),
    )


class ChainMembershipReuseTests(unittest.TestCase):
    def setUp(self):
        self.context = ("q-path", "t-path", "chrQ", "chrT", "+", "+")
        self.other_context = ("q-other", "t-other", "chrQ", "chrT", "+", "+")
        self.start = candidate("start", 0, 2, [membership(self.context, 1, 1)])
        # The first duplicate-context membership must retain legacy next(...)
        # behavior; the second membership would break the anchored chain.
        self.middle = candidate("middle", 2, 3, [
            membership(self.context, 2, 2), membership(self.context, 0, 0),
        ])
        self.end = candidate("end", 4, 4, [membership(self.context, 3, 3)])
        self.outside = candidate(
            "outside", 6, 50, [membership(self.other_context, 4, 4)],
        )

    def test_membership_is_computed_once_and_reused_with_anchors_and_edges(self):
        candidates = (self.start, self.middle, self.end, self.outside)
        with patch.object(
            chain_graph, "_path_membership_for_context",
            wraps=chain_graph._path_membership_for_context,
        ) as lookup:
            solved = _solve_context(
                candidates, self.context, {"start"}, {"end"}, 0.0,
            )
        self.assertEqual(lookup.call_count, len(candidates))
        self.assertEqual(solved["_prepared"][0], [0, 1, 2])
        self.assertEqual(solved["retained"], {"start", "middle", "end"})
        self.assertEqual(solved["retained_edges"], {
            ("start", "middle"), ("middle", "end"),
        })
        # Complete-path reporting keeps the top two paths; retained_edges
        # separately reflects the near-optimal threshold and excludes start→end.
        self.assertEqual(solved["complete_paths"], (
            (9.0, ("start", "middle", "end")),
            (6.0, ("start", "end")),
        ))

        with patch.object(
            chain_graph, "_path_membership_for_context",
            wraps=chain_graph._path_membership_for_context,
        ) as lookup:
            prepared = _solve_context(
                candidates, self.context, {"start"}, {"end"}, 0.0,
                prepared=solved["_prepared"],
            )
        lookup.assert_not_called()
        self.assertEqual(prepared, solved)

        public = ordered_candidate_chain(
            candidates, 0.0,
            context_anchor_ids={self.context: ({"start"}, {"end"})},
        )
        self.assertEqual(public.start_anchor_ids, frozenset({"start"}))
        self.assertEqual(public.end_anchor_ids, frozenset({"end"}))
        self.assertEqual(public.retained_edges, frozenset({
            ("start", "middle"), ("middle", "end"),
        }))
        self.assertEqual(public.best_path_member_ids, frozenset({
            "start", "middle", "end",
        }))

    def test_context_none_keeps_all_candidates_and_passes_none_paths(self):
        candidates = (self.start, self.outside)
        with patch.object(
            chain_graph, "_path_membership_for_context",
            wraps=chain_graph._path_membership_for_context,
        ) as lookup, patch.object(
            chain_graph, "_ordered_before",
            wraps=chain_graph._ordered_before,
        ) as ordered:
            solved = _solve_context(
                candidates, None, set(), set(), 0.0,
            )
        self.assertEqual(lookup.call_count, len(candidates))
        self.assertEqual(solved["_prepared"][0], [0, 1])
        self.assertEqual(solved["retained_edges"], {("start", "outside")})
        self.assertTrue(all(
            call.args[2:] == (None, None) for call in ordered.call_args_list
        ))


if __name__ == "__main__":
    unittest.main()
