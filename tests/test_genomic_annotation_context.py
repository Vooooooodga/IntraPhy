"""Contracts for annotation-conditioned non-exonic observations."""
import unittest
from types import SimpleNamespace

from intraphy.structure.annotation_context import annotation_supported_nonexonic_interval
from intraphy.structure.genomic_observations import genomic_span_observation
from intraphy.structure.native import NativeLocus
from intraphy.structure.types import Catalogue, ExonInstance, ExonSpan, Material
from intraphy.structure.alignment import FamilyAlignment
from intraphy.structure.build import build_catalogues


def _instance(identifier, species, start, end, strand="+"):
    return ExonInstance(identifier, "family", species, f"{species}_locus", "ctg",
                        start, end, strand, ("tx",))


def _setup(*, strand="+", partial=False, partial_location_known=True,
           partial_intervals=(), paths=None, extra_exons=(), anomalies=None,
           extra_span_map=None, target_row=None, source_row=None):
    if strand == "+":
        left_coords, right_coords, positive_coords = (10, 20), (40, 50), (52, 56)
    else:
        left_coords, right_coords, positive_coords = (40, 50), (10, 20), (4, 8)
    left = _instance("left", "target", *left_coords, strand=strand)
    right = _instance("right", "target", *right_coords, strand=strand)
    positive = _instance("positive", "target", *positive_coords, strand=strand)
    source_exon = _instance("candidate", "source", 25, 32)
    target_exons = (left, right, *extra_exons)
    target_paths = paths if paths is not None else {"tx": ("left", "right")}
    target = NativeLocus("family", "target", "target_locus", "ctg", strand,
        1, 60, "A" * 60, target_exons, target_paths, partial,
        partial_intervals=partial_intervals,
        partial_location_known=partial_location_known)
    source = NativeLocus("family", "source", "source_locus", "ctg", "+",
        1, 60, "A" * 60, (source_exon,), {"tx": ("candidate",)}, False,
        partial_location_known=True)
    spans = {"candidate": ExonSpan(25, 32), "positive": ExonSpan(52, 56)}
    exon_spans = {"left": ExonSpan(10, 20), "right": ExonSpan(40, 50),
                  "candidate": ExonSpan(25, 32), "positive": ExonSpan(52, 56)}
    exon_spans.update(extra_span_map or {})
    for exon in extra_exons:
        # Test helper uses IDs and alignment intervals supplied by the caller.
        if exon.id == "cover":
            exon_spans[exon.id] = ExonSpan(27, 29)
    rows = {"target": target_row or "A" * 60,
            "source": source_row or "A" * 60}
    alignment = SimpleNamespace(loci=(target, source), rows=rows,
        exons=exon_spans, anomalies=anomalies or {},
        columns={name: tuple(i for i, base in enumerate(row) if base != "-")
                 for name, row in rows.items()})
    return target, source, spans, exon_spans, alignment


def _qualified(**options):
    target, source, spans, _, alignment = _setup(**options)
    return annotation_supported_nonexonic_interval(species="target", source_species="source",
        candidate=spans["candidate"], source_id="candidate", locus=target,
        alignment=alignment, start=0, anchor_bases=2,
        minimum_identity=.7, anchor_identity=.8)


class GenomicAnnotationContextTests(unittest.TestCase):
    def test_retained_candidate_between_adjacent_annotated_exons_is_supported(self):
        self.assertTrue(_qualified())

    def test_strand_oriented_native_paths_support_both_directions(self):
        self.assertTrue(_qualified(strand="+"))
        self.assertTrue(_qualified(strand="-"))

    def test_annotation_absence_and_unlocated_partial_do_not_prove_intronic_status(self):
        self.assertFalse(_qualified(paths={}))
        self.assertFalse(_qualified(partial=True, partial_location_known=False))

    def test_legacy_records_without_whole_locus_context_remain_unresolved(self):
        target = SimpleNamespace(species="target", partial=False)
        alignment = SimpleNamespace(rows={"target": "A" * 60})
        self.assertFalse(annotation_supported_nonexonic_interval(
            species="target", source_species="source", candidate=ExonSpan(25, 32),
            source_id="candidate", locus=target, alignment=alignment, start=0,
            anchor_bases=2, minimum_identity=.7, anchor_identity=.8))

    def test_candidate_or_boundary_partial_blocks_but_distant_partial_does_not(self):
        self.assertFalse(_qualified(partial=True, partial_intervals=((18, 20),)))
        self.assertFalse(_qualified(partial=True, partial_intervals=((40, 42),)))
        self.assertFalse(_qualified(partial=True, partial_intervals=((26, 27),)))
        self.assertFalse(_qualified(partial=True, partial_intervals=((24, 25),)))
        self.assertTrue(_qualified(partial=True, partial_intervals=((5, 7),)))
        # Genomic intervals reverse into those same path boundaries on minus strand.
        self.assertFalse(_qualified(strand="-", partial=True,
                                    partial_intervals=((40, 42),)))
        self.assertFalse(_qualified(strand="-", partial=True,
                                    partial_intervals=((18, 20),)))
        self.assertTrue(_qualified(strand="-", partial=True,
                                   partial_intervals=((53, 55),)))

    def test_another_clean_path_pair_can_support_when_one_pair_boundary_is_partial(self):
        left2 = _instance("left2", "target", 8, 17)
        right2 = _instance("right2", "target", 42, 48)
        self.assertTrue(_qualified(partial=True, partial_intervals=((18, 20),),
            extra_exons=(left2, right2),
            extra_span_map={"left2": ExonSpan(8, 17), "right2": ExonSpan(42, 48)},
            paths={"tx1": ("left", "right"), "tx2": ("left2", "right2")}))

    def test_another_provided_transcript_exon_overlapping_candidate_vetoes_call(self):
        covering = _instance("cover", "target", 27, 29)
        self.assertFalse(_qualified(extra_exons=(covering,),
                                    paths={"tx": ("left", "right"),
                                           "alt": ("cover",)}))

    def test_intact_dna_and_source_candidate_and_both_anchors_are_required(self):
        self.assertFalse(_qualified(target_row="A" * 25 + "N" + "A" * 34))
        self.assertFalse(_qualified(target_row="A" * 25 + "-" * 7 + "A" * 28))
        weak_candidate = "A" * 25 + "C" * 7 + "A" * 28
        self.assertFalse(_qualified(source_row=weak_candidate))
        weak_anchor = "A" * 23 + "C" * 2 + "A" * 35
        self.assertFalse(_qualified(source_row=weak_anchor))

    def test_relevant_alignment_anomalies_block_annotation_conditioning(self):
        self.assertFalse(_qualified(anomalies={"candidate": ("competing_copy_candidates",)}))
        self.assertFalse(_qualified(anomalies={"left": ("reverse_correspondence_candidate",)}))

    def test_no_local_exon_uses_provided_annotation_context_and_auditable_reason(self):
        target, _, spans, _, alignment = _setup()
        proto = Catalogue("family", "unit", 60, (spans["candidate"],), (),
            boundary_candidates=(spans["candidate"],), observation_unit="genomic_exon_spans")
        observation = genomic_span_observation(species="target", proto=proto, locus=target,
            ids=("candidate",), spans={"candidate": spans["candidate"]},
            instances={"candidate": _instance("candidate", "source", 25, 32)},
            alignment=alignment, material_presence=(), all_material_presence={}, start=0,
            minimum_identity=.7, anchor_bases=2, anchor_identity=.8)
        self.assertEqual(observation.kind, "observed")
        self.assertEqual(observation.configurations[0].exons, ())
        self.assertIn("annotation_supported_nonexonic_interval", observation.reasons)
        self.assertNotIn("all_candidate_spans_qualified_absent", observation.reasons)

    def test_build_catalogues_preserves_whole_locus_context_outside_local_unit(self):
        target, source, _, exon_spans, alignment = _setup()
        native = FamilyAlignment(alignment.loci, alignment.rows,
            {species: 0 for species in alignment.rows}, alignment.columns,
            exon_spans, alignment.anomalies, ())
        catalogues, _, _, _ = build_catalogues(native,
            observation_unit="genomic_exon_spans", anchor_bases=2, anchor_identity=.8)
        candidate_catalogue = next(catalogue for catalogue in catalogues
            if any(exon.id == "candidate" for exon in catalogue.exon_instances))
        observation = next(value for value in candidate_catalogue.observations
                           if value.species == "target")
        self.assertEqual(candidate_catalogue.alignment_offset, 25)
        self.assertEqual(candidate_catalogue.spans, (ExonSpan(0, 7),))
        self.assertEqual(observation.kind, "observed")
        self.assertEqual(observation.configurations[0].exons, ())
        self.assertIn("annotation_supported_nonexonic_interval", observation.reasons)

    def test_positive_local_exons_can_coexist_with_independent_nonexonic_candidate(self):
        target, _, spans, _, alignment = _setup(
            paths={"tx": ("left", "right", "positive")},
            extra_exons=(_instance("positive", "target", 52, 56),))
        proto = Catalogue("family", "unit", 60,
            (spans["candidate"], spans["positive"]), (),
            boundary_candidates=(spans["candidate"], spans["positive"]),
            observation_unit="genomic_exon_spans")
        instances = {"candidate": _instance("candidate", "source", 25, 32),
                     "positive": _instance("positive", "target", 52, 56)}
        observation = genomic_span_observation(species="target", proto=proto, locus=target,
            ids=tuple(instances), spans=spans, instances=instances, alignment=alignment,
            material_presence=(), all_material_presence={}, start=0,
            minimum_identity=.7, anchor_bases=2, anchor_identity=.8)
        self.assertEqual(observation.configurations[0].exons, (spans["positive"],))
        self.assertIn("annotation_supported_nonexonic_interval", observation.reasons)

    def test_remote_partial_context_remains_usable_while_unknown_is_preserved(self):
        self.assertTrue(_qualified(partial=True, partial_intervals=((5, 7),)))
        target, _, spans, _, alignment = _setup(
            partial=True, partial_intervals=((5, 7),))
        proto = Catalogue("family", "unit", 60, (spans["candidate"],), (),
            boundary_candidates=(spans["candidate"],), observation_unit="genomic_exon_spans")
        observation = genomic_span_observation(species="target", proto=proto, locus=target,
            ids=("candidate",), spans={"candidate": spans["candidate"]},
            instances={"candidate": _instance("candidate", "source", 25, 32)},
            alignment=alignment, material_presence=(), all_material_presence={}, start=0,
            minimum_identity=.7, anchor_bases=2, anchor_identity=.8)
        self.assertEqual(observation.kind, "unknown")

    def test_dna_deletion_remains_a_distinct_qualified_absence(self):
        target, _, spans, _, alignment = _setup(target_row="A" * 25 + "-" * 7 + "A" * 28)
        material = Material("tract", 25, 32)
        proto = Catalogue("family", "unit", 60, (spans["candidate"],), (),
            material=(material,), boundary_candidates=(spans["candidate"],),
            observation_unit="genomic_exon_spans")
        observation = genomic_span_observation(species="target", proto=proto, locus=target,
            ids=("candidate",), spans={"candidate": spans["candidate"]},
            instances={"candidate": _instance("candidate", "source", 25, 32)},
            alignment=alignment, material_presence=(0,), all_material_presence={}, start=0,
            minimum_identity=.7, anchor_bases=2, anchor_identity=.8)
        self.assertEqual(observation.kind, "observed")
        self.assertIn("all_candidate_spans_qualified_absent", observation.reasons)
        self.assertNotIn("annotation_supported_nonexonic_interval", observation.reasons)
        self.assertEqual(observation.configurations[0].material, (0,))


if __name__ == "__main__":
    unittest.main()
