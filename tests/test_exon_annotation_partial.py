"""NCBI partial-boundary handling for physical genomic exon observations."""
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

from intraphy.partial_boundaries import range_value
from intraphy.preparation.transcripts import _path_role_record
from intraphy.structure.genomic_observations import genomic_span_observation
from intraphy.structure.native import NativeLocus, _partial_intervals, load_native
from intraphy.structure.observations import compatibility
from intraphy.structure.space import enumerate_space
from intraphy.structure.types import (Catalogue, ExonConfiguration, ExonInstance,
                                      ExonSpan, ObservationEvidence)


class PartialAnnotationTests(unittest.TestCase):
    def test_genomic_start_range_maps_to_transcript_boundary_by_strand(self):
        for strand, expected in (("+", ("1", "0")), ("-", ("0", "1"))):
            feature = {"id": "e", "type": "exon", "start": 10, "end": 20,
                       "strand": strand, "attrs": {"partial": "true",
                                                       "start_range": ".,12"}}
            record = _path_role_record("tx", 1, [feature], 0, feature)
            self.assertEqual((record["partial_start"], record["partial_end"]), expected)

    def test_propagated_generic_partial_uses_sibling_range_not_every_exon(self):
        attrs = {"partial": "true"}
        first = {"start": "100", "end": "150", "original_attributes":
                 "ID=e1;partial=true;start_range=.,105"}
        second = {"start": "300", "end": "350", "original_attributes":
                  "ID=e2;partial=true"}
        intervals, localized, generic, unresolved = _partial_intervals(
            {"tx": [first, second]}, 1, 500)
        self.assertEqual(intervals, ((0, 104),))
        self.assertTrue(localized)
        self.assertTrue(generic)
        self.assertFalse(unresolved)
        self.assertEqual(range_value(attrs, "start"), "")

    def test_internal_and_terminal_ranges_are_clipped_to_neighbor_exons(self):
        rows = [
            {"start": "100", "end": "150", "original_attributes": "ID=e1"},
            {"start": "300", "end": "350", "original_attributes":
             "ID=e2;start_range=285,300;end_range=345,."},
            {"start": "400", "end": "450", "original_attributes": "ID=e3"},
        ]
        intervals, localized, _, unresolved = _partial_intervals({"tx": rows}, 1, 500)
        self.assertEqual(intervals, ((284, 299), (345, 399)))
        self.assertTrue(localized)
        self.assertFalse(unresolved)

    def test_noncanonical_or_unbounded_ranges_remain_unlocalized(self):
        for attributes in ("ID=e;start_range=.,.", "ID=e;start_range=100,.",
                           "ID=e;end_range=.,200"):
            row = {"start": "100", "end": "150", "original_attributes": attributes}
            intervals, localized, partial, unresolved = _partial_intervals(
                {"tx": [row]}, 1, 500)
            self.assertEqual(intervals, ())
            self.assertFalse(localized)
            self.assertTrue(partial)
            self.assertTrue(unresolved)

    def test_old_path_attributes_recover_range_and_generic_only_stays_conservative(self):
        located = {"start": "100", "end": "150",
                   "original_attributes": "ID=e;partial=true;start_range=.,105"}
        intervals, localized, _, unresolved = _partial_intervals({"tx": [located]}, 1, 500)
        self.assertEqual(intervals, ((0, 104),))
        self.assertTrue(localized)
        self.assertFalse(unresolved)

    def test_load_native_recovers_selected_raw_transcript_range_only(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "gene_loci.fasta").write_text(">A|g|locus\n" + "A" * 500 + "\n")
            (root / "gene_loci.tsv").write_text(
                "species\tgene_copy_id\tcontig\tstrand\tsearch_start\tsearch_end\n"
                "A\tg\tchr1\t+\t1\t500\n")
            (root / "transcript_paths.tsv").write_text(
                "family_id\tspecies\tgene_copy_id\toccurrence_id\ttranscript_id\trole\tstart\tend\t"
                "annotation_source\tcds_intervals\tcds_phase\toriginal_attributes\tpartial_start\tpartial_end\n"
                "f\tA\tg\te\ttx\texon\t100\t150\tfixture\t\t.\tID=e\t0\t0\n")
            (root / "raw_gene_features.tsv").write_text(
                "family_id\tspecies\tgene_copy_id\tgene_id\tseqid\tsource\ttype\tstart\tend\tstrand\t"
                "phase\tid\tname\tparent\tparents\townership\tattrs\n"
                "f\tA\tg\tg1\tchr1\ttest\texon\t100\t150\t+\t.\te\tNA\ttx\ttx\t"
                "target_gene_descendant\tID=e\n"
                "f\tA\tg\tg1\tchr1\ttest\tmRNA\t100\t200\t+\t.\ttx\tNA\tg1\tg1\t"
                "target_gene_descendant\tID=tx;partial=true;start_range=.,105\n"
                "f\tA\tg\tg1\tchr1\ttest\tmRNA\t250\t300\t+\t.\tother\tNA\tg1\tg1\t"
                "overlapping_context\tID=other;partial=true;start_range=.,260\n")
            locus = load_native(root)[0]
            self.assertTrue(locus.partial)
            self.assertTrue(locus.partial_location_known)
            self.assertEqual(locus.partial_intervals, ((0, 104),))
            raw_path = root / "raw_gene_features.tsv"
            raw_path.write_text(raw_path.read_text().replace(
                "ID=tx;partial=true;start_range=.,105", "ID=tx;partial=true"))
            generic_parent = load_native(root)[0]
            self.assertTrue(generic_parent.partial)
            self.assertFalse(generic_parent.partial_location_known)
            self.assertEqual(generic_parent.partial_intervals, ())

    def test_direct_legacy_partial_locus_defaults_to_unlocated(self):
        locus = NativeLocus("f", "A", "g", "chr1", "+", 1, 4, "AAAA", (), {"tx": ()}, True)
        self.assertFalse(locus.partial_location_known)
        span = ExonSpan(0, 2)
        instance = ExonInstance("e", "f", "A", "g", "chr1", 1, 3, "+", ("tx",))
        proto = Catalogue("f", "u", 4, (span,), (), observation_unit="genomic_exon_spans")
        observation = genomic_span_observation(species="A", proto=proto, locus=locus,
            ids=("e",), spans={"e": span}, instances={"e": instance},
            alignment=SimpleNamespace(rows={"A": "AAAA"}, loci=(locus,)),
            material_presence=(), all_material_presence={}, start=0,
            minimum_identity=.7, anchor_bases=2, anchor_identity=.8)
        self.assertIn(ExonSpan(0, 4), observation.unknown_intervals)

    def test_partial_uncertainty_maps_on_reverse_alignment_axis_and_stays_unknown(self):
        span = ExonSpan(5, 15)
        instance = ExonInstance("e", "f", "A", "g", "chr", 10, 20, "-", ("tx",))
        proto = Catalogue("f", "u", 20, (span,), (), observation_unit="genomic_exon_spans")
        locus = SimpleNamespace(species="A", strand="-", search_start=100,
            search_end=119, partial=True, partial_location_known=True,
            partial_intervals=((109, 113),), paths={"tx": ("e",)})
        alignment = SimpleNamespace(rows={"A": "A" * 20}, loci=(locus,),
            columns={"A": tuple(range(20))})
        observation = genomic_span_observation(species="A", proto=proto, locus=locus,
            ids=("e",), spans={"e": span}, instances={"e": instance},
            alignment=alignment, material_presence=(), all_material_presence={}, start=0,
            minimum_identity=.7, anchor_bases=2, anchor_identity=.8)
        self.assertEqual(observation.configurations[0].exons, (span,))
        self.assertIn(ExonSpan(6, 10), observation.unknown_intervals)
        self.assertIn("localized_partial_annotation_boundary", observation.reasons)

        other_proto = Catalogue("f", "other", 5, (ExonSpan(0, 5),), (),
                                 observation_unit="genomic_exon_spans")
        other = genomic_span_observation(species="A", proto=other_proto, locus=locus,
            ids=("e",), spans={"e": ExonSpan(0, 5)}, instances={"e": instance},
            alignment=SimpleNamespace(rows={"A": "A" * 32}, loci=(locus,),
                                      columns={"A": tuple(range(32))}),
            material_presence=(), all_material_presence={}, start=12,
            minimum_identity=.7, anchor_bases=2, anchor_identity=.8)
        self.assertEqual(other.unknown_intervals, ())
        self.assertNotIn("localized_partial_annotation_boundary", other.reasons)

    def test_unlocated_legacy_partial_does_not_become_qualified_absence(self):
        span = ExonSpan(5, 15)
        instance = ExonInstance("e", "f", "B", "g", "chr", 10, 20, "+", ("tx",))
        proto = Catalogue("f", "u", 20, (span,), (), observation_unit="genomic_exon_spans")
        locus = SimpleNamespace(species="A", partial=True, partial_location_known=False,
                                paths={})
        alignment = SimpleNamespace(rows={"A": "A" * 5 + "-" * 10 + "A" * 5,
                                          "B": "A" * 20},
                                   loci=(locus, SimpleNamespace(species="B")))
        observation = genomic_span_observation(species="A", proto=proto, locus=locus,
            ids=("e",), spans={"e": span}, instances={"e": instance},
            alignment=alignment, material_presence=(), all_material_presence={}, start=0,
            minimum_identity=.7, anchor_bases=2, anchor_identity=.8)
        self.assertEqual(observation.kind, "unknown")
        self.assertNotIn("all_candidate_spans_qualified_absent", observation.reasons)

        remote = SimpleNamespace(species="A", strand="+", search_start=1,
            search_end=24, partial=True, partial_location_known=True,
            partial_intervals=((20, 24),), paths={"tx": ()})
        remote_proto = Catalogue("f", "remote", 20, (ExonSpan(2, 18),), (),
                                 observation_unit="genomic_exon_spans")
        remote_alignment = SimpleNamespace(
            rows={"A": "A" * 12 + "-" * 16 + "A" * 12, "B": "A" * 40},
            loci=(remote, SimpleNamespace(species="B")),
            columns={"A": (*range(12), *range(28, 40))})
        localized = genomic_span_observation(species="A", proto=remote_proto, locus=remote,
            ids=("e",), spans={"e": ExonSpan(2, 18)}, instances={"e": instance},
            alignment=remote_alignment, material_presence=(), all_material_presence={}, start=10,
            minimum_identity=.7, anchor_bases=2, anchor_identity=.8)
        self.assertEqual(localized.kind, "observed")
        self.assertIn("all_candidate_spans_qualified_absent", localized.reasons)

    def test_uncertain_endpoint_does_not_mask_known_exon_interior(self):
        catalogue = Catalogue("f", "u", 12,
            (ExonSpan(0, 10), ExonSpan(0, 12)), (), observation_unit="genomic_exon_spans")
        space = enumerate_space(catalogue)
        observation = ObservationEvidence("A", (ExonConfiguration((ExonSpan(0, 10),)),),
            "partial", unknown_intervals=(ExonSpan(0, 3),))
        allowed = compatibility(space, observation)
        accepted = {state.exons for state, value in zip(space.states, allowed) if value}
        self.assertIn((ExonSpan(0, 10),), accepted)
        self.assertNotIn((ExonSpan(0, 12),), accepted)


if __name__ == "__main__":
    unittest.main()
