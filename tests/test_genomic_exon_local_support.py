"""Local exon support preserves only strict shared-span evidence."""
from types import SimpleNamespace
import unittest

from intraphy.structure.genomic_observations import genomic_span_observation
from intraphy.structure.observations import compatibility
from intraphy.structure.space import enumerate_space
from intraphy.structure.types import Catalogue, ExonInstance, ExonSpan, Material


def _evaluate(rows, features, *, extra_spans=(), material=(), presence=(),
              anchor_bases=4, minimum_identity=.7, locus_options=None):
    length = len(next(iter(rows.values())))
    spans, instances, loci = {}, {}, []
    locus_options = locus_options or {}
    for species, entries in features.items():
        for identifier, span in entries:
            spans[identifier] = span
            # Empty CDS records model sequence-only support for UTR/noncoding exons.
            instances[identifier] = ExonInstance(identifier, "family", species,
                "locus", "contig", span.start, span.end, "+", ("tx",), cds=())
    for species, row in rows.items():
        options = dict(locus_options.get(species, {}))
        default_paths = {"tx": tuple(identifier for identifier, _ in features.get(species, ())) }
        values = dict(species=species, partial=False, partial_location_known=False,
            partial_intervals=(), strand="+", search_start=1,
            search_end=sum(base != "-" for base in row),
            paths=default_paths if default_paths["tx"] else {})
        values.update(options)
        loci.append(SimpleNamespace(**values))
    candidates = tuple(sorted(set(spans.values()) | set(extra_spans)))
    proto = Catalogue("family", "unit", length, candidates, (), material=material,
        boundary_candidates=candidates, observation_unit="genomic_exon_spans")
    alignment = SimpleNamespace(rows=rows, loci=tuple(loci),
        columns={species: tuple(i for i, base in enumerate(row) if base != "-")
                 for species, row in rows.items()})
    target = next(iter(rows))
    locus = next(value for value in loci if value.species == target)
    observation = genomic_span_observation(species=target, proto=proto, locus=locus,
        ids=tuple(spans), spans=spans, instances=instances, alignment=alignment,
        material_presence=presence, all_material_presence={}, start=0,
        minimum_identity=minimum_identity, anchor_bases=anchor_bases,
        anchor_identity=.8)
    space = enumerate_space(proto)
    return observation, space, compatibility(space, observation)


class GenomicExonLocalSupportTests(unittest.TestCase):
    def test_noncoding_shared_exon_localizes_uncertainty_and_state_compatibility(self):
        core, neighbor = ExonSpan(10, 30), ExonSpan(60, 70)
        rows = {"A": "A" * 100, "B": "C" * 10 + "A" * 20 + "C" * 70}
        observation, space, allowed = _evaluate(rows, {
            "A": (("a_core_utr", core),),
            "B": (("b_core_utr", core), ("b_neighbor", neighbor)),
        }, extra_spans=(ExonSpan(10, 14), ExonSpan(16, 30)))

        self.assertEqual(observation.kind, "partial")
        self.assertEqual(observation.configurations[0].exons, (core,))
        self.assertIn("unresolved_sequence_correspondence", observation.reasons)
        self.assertIn("local_shared_exon_correspondence", observation.reasons)
        self.assertNotIn(ExonSpan(0, 100), observation.unknown_intervals)
        self.assertIn(ExonSpan(60, 70), observation.unknown_intervals)

        accepted = {state.exons for state, keep in zip(space.states, allowed) if keep}
        split = (ExonSpan(10, 14), ExonSpan(16, 30))
        self.assertIn(split, {state.exons for state in space.states})
        self.assertIn((core,), accepted)
        self.assertIn((core, neighbor), accepted)
        # A fusion extending from supported exon sequence into masked sequence remains unresolved.
        self.assertIn((ExonSpan(10, 70),), accepted)
        # An internal split changes geometry inside the independently supported exon.
        self.assertNotIn(split, accepted)

    def test_unlocalized_or_weak_same_boundary_matches_do_not_create_local_support(self):
        span = ExonSpan(10, 30)
        cases = {
            "low_identity": "C" * 100,
            "ambiguous_base": "C" * 10 + "A" * 5 + "N" + "A" * 14 + "C" * 70,
            "alignment_gap": "C" * 10 + "A" * 5 + "-" + "A" * 14 + "C" * 70,
        }
        for label, source in cases.items():
            with self.subTest(label=label):
                observation, _, allowed = _evaluate({"A": "A" * 100, "B": source}, {
                    "A": (("a", span),), "B": (("b", span),)})
                self.assertIn(ExonSpan(0, 100), observation.unknown_intervals)
                self.assertNotIn("local_shared_exon_correspondence", observation.reasons)
                self.assertTrue(all(allowed))

        short = ExonSpan(10, 13)
        observation, _, allowed = _evaluate({"A": "A" * 100,
            "B": "C" * 10 + "A" * 3 + "C" * 87}, {
                "A": (("a", short),), "B": (("b", short),)}, anchor_bases=4)
        self.assertIn(ExonSpan(0, 100), observation.unknown_intervals)
        self.assertNotIn("local_shared_exon_correspondence", observation.reasons)
        self.assertTrue(all(allowed))

    def test_partial_pathless_and_conflicting_counterparts_cannot_supply_fallback(self):
        core = ExonSpan(10, 30)
        base_rows = {"A": "A" * 100, "B": "C" * 10 + "A" * 20 + "C" * 70}
        cases = (
            ({"B": {"partial": True}}, (("b", core),)),
            ({"B": {"paths": {}}}, (("b", core),)),
            ({"B": {"paths": {"tx": ()}}}, (("b", core),)),
            ({}, (("b", core), ("b_adjacent", ExonSpan(30, 40)))),
            ({}, (("b", core), ("b_overlap", ExonSpan(29, 40)))),
        )
        for options, source_features in cases:
            with self.subTest(options=options, source_features=source_features):
                observation, _, allowed = _evaluate(base_rows, {
                    "A": (("a", core),), "B": source_features,
                }, locus_options=options)
                self.assertIn(ExonSpan(0, 100), observation.unknown_intervals)
                self.assertNotIn("local_shared_exon_correspondence", observation.reasons)
                self.assertTrue(all(allowed))

    def test_existing_partial_mask_is_preserved_over_local_support(self):
        core = ExonSpan(10, 30)
        rows = {"A": "A" * 100, "B": "C" * 10 + "A" * 20 + "C" * 70}
        observation, space, allowed = _evaluate(rows, {
            "A": (("a", core),),
            "B": (("b", core), ("b_neighbor", ExonSpan(60, 70))),
        }, extra_spans=(ExonSpan(10, 14), ExonSpan(16, 30)), locus_options={
            "A": {"partial": True, "partial_location_known": True,
                  "partial_intervals": ((14, 16),)}})
        self.assertIn(ExonSpan(14, 16), observation.unknown_intervals)
        self.assertIn("localized_partial_annotation_boundary", observation.reasons)
        accepted = {state.exons for state, keep in zip(space.states, allowed) if keep}
        self.assertIn((ExonSpan(10, 14), ExonSpan(16, 30)), accepted)

        unlocated, _, _ = _evaluate(rows, {
            "A": (("a", core),),
            "B": (("b", core), ("b_neighbor", ExonSpan(60, 70))),
        }, locus_options={"A": {"partial": True, "partial_location_known": False}})
        self.assertIn(ExonSpan(0, 100), unlocated.unknown_intervals)

    def test_global_pass_is_unchanged_and_no_shared_interval_stays_whole_unknown(self):
        core = ExonSpan(10, 30)
        features = {"A": (("a", core),), "B": (("b", core),)}
        observed, observed_space, observed_allowed = _evaluate(
            {"A": "A" * 100, "B": "A" * 100}, features)
        self.assertEqual(observed.kind, "observed")
        self.assertEqual(observed.unknown_intervals, ())
        self.assertNotIn("unresolved_sequence_correspondence", observed.reasons)
        self.assertEqual({state.exons for state, keep in zip(
            observed_space.states, observed_allowed) if keep}, {(core,)})

        unsupported, _, unsupported_allowed = _evaluate({"A": "A" * 100, "B": "C" * 100}, {
            "A": (("a", core),), "B": (("b", ExonSpan(60, 70)),)})
        self.assertIn(ExonSpan(0, 100), unsupported.unknown_intervals)
        self.assertNotIn("local_shared_exon_correspondence", unsupported.reasons)
        self.assertTrue(all(unsupported_allowed))

    def test_material_absence_does_not_spread_to_an_unqualified_neighbor(self):
        deletion, neighbor = ExonSpan(40, 50), ExonSpan(70, 80)
        target = "A" * 40 + "-" * 10 + "A" * 50
        observation, space, allowed = _evaluate({"A": target, "B": "C" * 100}, {
            "A": (), "B": (("b_deleted", deletion), ("b_neighbor", neighbor)),
        }, material=(Material("tract", 40, 50),), presence=(0,))

        self.assertEqual(observation.kind, "unknown")
        self.assertNotIn("all_candidate_spans_qualified_absent", observation.reasons)
        self.assertIn("genomic_absence_unqualified", observation.reasons)
        accepted = [state for state, keep in zip(space.states, allowed) if keep]
        self.assertTrue(any(neighbor in state.exons and state.material != (1,)
                            for state in accepted))
        self.assertFalse(any(state.material == (1,) for state in accepted))


if __name__ == "__main__":
    unittest.main()
