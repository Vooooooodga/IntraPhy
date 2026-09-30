"""Bounded conditioning of within-exon sequence indels in genomic catalogues."""
from types import SimpleNamespace
import unittest

from intraphy.structure.alignment import FamilyAlignment
from intraphy.structure.build import build_catalogues
from intraphy.structure.space import enumerate_space
from intraphy.structure.types import ExonInstance, ExonSpan


def _build(rows, exon_spans, *, observation_unit="genomic_exon_spans", anomalies=None,
           partial=None):
    instances, spans, loci = {}, {}, []
    for species, physical in exon_spans.items():
        exon_ids = []
        for number, span in enumerate(physical):
            identifier = f"{species}_e{number}"
            exon_ids.append(identifier)
            instances[identifier] = ExonInstance(identifier, "family", species, "locus",
                "contig", span.start, span.end, "+", ("tx1",))
            spans[identifier] = span
        attributes = (partial or {}).get(species, {})
        loci.append(SimpleNamespace(family="family", species=species,
            exons=tuple(instances[eid] for eid in exon_ids),
            paths={"tx1": tuple(exon_ids)} if exon_ids else {},
            partial=attributes.get("partial", False),
            partial_intervals=attributes.get("partial_intervals", ()),
            partial_location_known=attributes.get("partial_location_known", False),
            strand=attributes.get("strand", "+"),
            search_start=attributes.get("search_start", 1),
            search_end=attributes.get("search_end", len(rows[species].replace("-", "")))))
    columns = {species: tuple(i for i, base in enumerate(row) if base != "-")
               for species, row in rows.items()}
    alignment = FamilyAlignment(tuple(loci), rows, {species: 0 for species in rows},
        columns, spans, anomalies or {}, ())
    return build_catalogues(alignment, observation_unit=observation_unit,
                            anchor_bases=2, anchor_identity=.8)[0]


class StructuralMaterialConditioningTests(unittest.TestCase):
    def test_internal_indel_is_projected_out_and_unknown_evidence_is_preserved(self):
        rows = {"A": "A"*40, "B": "A"*10 + "-"*4 + "A"*26,
                "C": "A"*10 + "N"*4 + "A"*26}
        catalogues = _build(rows, {species: (ExonSpan(0, 40),) for species in rows})
        self.assertEqual(len(catalogues), 1)
        catalogue = catalogues[0]
        self.assertEqual(catalogue.material, ())
        self.assertIn("conditioned_internal_sequence_indels", catalogue.reasons)
        for observation in catalogue.observations:
            self.assertTrue(all(not config.material for config in observation.configurations))
        ambiguous = next(o for o in catalogue.observations if o.species == "C")
        self.assertEqual(ambiguous.kind, "partial")
        self.assertIn(ExonSpan(10, 14), ambiguous.unknown_intervals)
        self.assertEqual(ambiguous.material_presence, ())
        space = enumerate_space(catalogue)
        self.assertTrue(space.complete)
        self.assertEqual(space.diagnostics["material_combinations"], 1)

    def test_internal_indel_with_known_dna_is_projected_out_too(self):
        rows = {"A": "A"*40, "B": "A"*10 + "-"*4 + "A"*26}
        catalogues = _build(rows, {species: (ExonSpan(0, 40),) for species in rows})
        self.assertEqual(catalogues[0].material, ())
        self.assertEqual(tuple(o.material_presence for o in catalogues[0].observations), ((), ()))

    def test_low_support_full_geometry_remains_unknown_after_projection(self):
        rows = {"A": "A"*40, "B": "A"*10 + "-"*4 + "A"*26, "C": "C"*40}
        catalogues = _build(rows, {species: (ExonSpan(0, 40),) for species in rows})
        self.assertEqual(catalogues[0].material, ())
        unsupported = next(o for o in catalogues[0].observations if o.species == "C")
        self.assertEqual(unsupported.kind, "partial")
        self.assertIn(ExonSpan(0, 40), unsupported.unknown_intervals)
        self.assertIn("unresolved_sequence_correspondence", unsupported.reasons)

    def test_boundary_touching_and_split_spacer_tracts_are_retained(self):
        rows = {"A": "A"*40, "B": "A"*10 + "-"*4 + "A"*26}
        catalogues = _build(rows, {"A": (ExonSpan(0, 40),),
            "B": (ExonSpan(0, 10), ExonSpan(14, 40))})
        self.assertEqual(len(catalogues[0].material), 1)

    def test_whole_exon_and_shared_deletion_tracts_are_retained(self):
        rows = {"A": "A"*40, "B": "A"*10 + "-"*4 + "A"*26}
        whole_exon = _build(rows, {"A": (ExonSpan(10, 14),), "B": ()})
        self.assertEqual(len(whole_exon[0].material), 1)

        deletion_rows = {"A": "A"*40, "B": "A"*4 + "-"*32 + "A"*4}
        shared = _build(deletion_rows, {"A": (ExonSpan(5, 15), ExonSpan(25, 35)), "B": ()})
        self.assertEqual(len(shared), 1)
        self.assertEqual(len(shared[0].material), 1)
        self.assertEqual(enumerate_space(shared[0]).diagnostics["material_combinations"], 3)

    def test_missing_or_conflicted_native_geometry_and_alignment_anomaly_retain_tract(self):
        rows = {"A": "A"*40, "B": "A"*10 + "-"*4 + "A"*26}
        missing = _build(rows, {"A": (ExonSpan(0, 40),), "B": ()})
        overlap = _build(rows, {"A": (ExonSpan(0, 40),),
            "B": (ExonSpan(0, 20), ExonSpan(18, 35))})
        anomalous = _build(rows, {species: (ExonSpan(0, 40),) for species in rows},
            anomalies={"A_e0": ("alignment_anomaly",)})
        for catalogue in (missing[0], overlap[0], anomalous[0]):
            self.assertEqual(len(catalogue.material), 1)

    def test_duplicate_physical_aliases_do_not_prevent_conditioning(self):
        rows = {"A": "A"*40, "B": "A"*10 + "-"*4 + "A"*26}
        catalogues = _build(rows, {"A": (ExonSpan(0, 40), ExonSpan(0, 40)),
                                   "B": (ExonSpan(0, 40),)})
        self.assertEqual(catalogues[0].material, ())
        self.assertEqual(next(o for o in catalogues[0].observations if o.species == "A")
                         .configurations[0].exons, (ExonSpan(0, 40),))

    def test_unlocated_or_overlapping_partial_window_retains_tract(self):
        rows = {"A": "A"*40, "B": "A"*10 + "-"*4 + "A"*26}
        exons = {species: (ExonSpan(0, 40),) for species in rows}
        generic = _build(rows, exons, partial={"A": {"partial": True}})
        overlapping = _build(rows, exons, partial={"A": {"partial": True,
            "partial_location_known": True, "partial_intervals": ((9, 12),)}})
        distant_terminal = _build(rows, exons, partial={"A": {"partial": True,
            "partial_location_known": True, "partial_intervals": ((0, 5),)}})
        self.assertEqual(len(generic[0].material), 1)
        self.assertEqual(len(overlapping[0].material), 1)
        self.assertEqual(distant_terminal[0].material, ())

    def test_legacy_transcript_catalogue_is_unchanged(self):
        rows = {"A": "A"*40, "B": "A"*10 + "-"*4 + "A"*26}
        catalogues = _build(rows, {species: (ExonSpan(0, 40),) for species in rows},
                            observation_unit="transcript_configuration")
        self.assertEqual(len(catalogues[0].material), 1)
        self.assertNotIn("conditioned_internal_sequence_indels", catalogues[0].reasons)


if __name__ == "__main__":
    unittest.main()
