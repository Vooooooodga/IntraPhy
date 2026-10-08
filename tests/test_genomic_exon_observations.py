"""Regression contracts for physical-span observations (written, not executed here)."""
from types import SimpleNamespace
from pathlib import Path
import json
import shutil
import tempfile
import unittest

from intraphy.structure.genomic_observations import genomic_span_observation
from intraphy.structure.types import Catalogue, ExonInstance, ExonSpan, Material
from intraphy.structure.build import build_catalogues
from intraphy.structure.alignment import FamilyAlignment
from intraphy.structure.serialization import decode_catalogue, encode_catalogue


def _instance(identifier, species, start, end, transcripts=("tx1",)):
    return ExonInstance(identifier, "family", species, "locus", "contig", start, end,
                        "+", tuple(transcripts))


def _call(species, spans, instances, rows, *, materials=(), presence=(), partial=False, paths=None):
    proto = Catalogue("family", "unit", len(next(iter(rows.values()))), tuple(spans.values()), (),
        tuple(materials), boundary_candidates=tuple(spans.values()), observation_unit="genomic_exon_spans")
    loci = tuple(SimpleNamespace(species=key) for key in rows)
    locus = SimpleNamespace(species=species, partial=partial,
                            paths=paths if paths is not None else {"tx1": ()})
    alignment = SimpleNamespace(rows=rows, loci=loci)
    return genomic_span_observation(species=species, proto=proto, locus=locus,
        ids=tuple(spans), spans=spans, instances=instances, alignment=alignment,
        material_presence=presence, all_material_presence={}, start=0,
        minimum_identity=.7, anchor_bases=2, anchor_identity=.8)


class GenomicExonObservationTests(unittest.TestCase):
    def _build_alignment_catalogues(self, rows, exon_spans):
        instances = {key: _instance(key, species, span.start, span.end)
                     for species, values in exon_spans.items()
                     for key, span in values.items()}
        spans = {key: span for values in exon_spans.values() for key, span in values.items()}
        loci = tuple(SimpleNamespace(family="family", species=species,
            exons=tuple(instances[key] for key in values),
            paths={"tx1": tuple(values)} if values else {}, partial=False)
            for species, values in exon_spans.items())
        alignment = FamilyAlignment(loci, rows, {species: 0 for species in rows},
            {species: tuple(i for i, base in enumerate(row) if base != "-")
             for species, row in rows.items()}, spans, {}, ())
        return build_catalogues(alignment, observation_unit="genomic_exon_spans",
                                anchor_bases=2, anchor_identity=.8)[0]

    def test_complete_exon_spans_remain_separate_from_cds_subintervals(self):
        spans = {"coding": ExonSpan(0, 12), "noncoding": ExonSpan(20, 25)}
        instances = {
            "coding": ExonInstance("coding", "family", "A", "locus", "contig",
                100, 112, "+", ("tx1",), cds=((104, 109, 0),)),
            "noncoding": _instance("noncoding", "A", 120, 125),
        }
        observation = _call("A", spans, instances, {"A": "A"*60, "B": "A"*60})
        self.assertEqual(observation.configurations[0].exons,
                         (ExonSpan(0, 12), ExonSpan(20, 25)))
        self.assertEqual(instances["coding"].cds, ((104, 109, 0),))

    def test_physical_span_dedup_ignores_transcript_path_count(self):
        spans = {"e1": ExonSpan(10, 20), "e1_copy": ExonSpan(10, 20),
                 "e2": ExonSpan(30, 40)}
        instances = {"e1": _instance("e1", "A", 10, 20, ("tx1", "tx2")),
                     "e1_copy": _instance("e1_copy", "A", 10, 20, ("tx3",)),
                     "e2": _instance("e2", "A", 30, 40, ("tx1",))}
        observation = _call("A", spans, instances, {"A": "A"*60, "B": "A"*60})
        self.assertEqual(observation.configurations[0].exons,
                         (ExonSpan(10, 20), ExonSpan(30, 40)))

    def test_overlapping_physical_spans_are_unknown_and_unmerged(self):
        spans = {"e1": ExonSpan(5, 25), "e2": ExonSpan(20, 35)}
        instances = {key: _instance(key, "A", span.start, span.end)
                     for key, span in spans.items()}
        observation = _call("A", spans, instances, {"A": "A"*60, "B": "A"*60})
        self.assertEqual(observation.kind, "unknown")
        self.assertEqual(observation.configurations, ())
        self.assertEqual(len(instances), 2)

    def test_split_and_fusion_geometry_is_retained_without_spacer_masking(self):
        spans = {"fusion": ExonSpan(0, 40), "left": ExonSpan(0, 20),
                 "right": ExonSpan(25, 40)}
        instances = {"fusion": _instance("fusion", "A", 0, 40),
                     "left": _instance("left", "B", 0, 20),
                     "right": _instance("right", "B", 25, 40)}
        rows = {"A": "A"*50, "B": "A"*50}
        fusion = _call("A", spans, instances, rows)
        split = _call("B", spans, instances, rows)
        self.assertEqual(fusion.configurations[0].exons, (ExonSpan(0, 40),))
        self.assertEqual(split.configurations[0].exons, (ExonSpan(0, 20), ExonSpan(25, 40)))
        self.assertNotIn(ExonSpan(20, 25), split.unknown_intervals)

    def test_ambiguous_window_keeps_local_positive_exon(self):
        spans = {"local": ExonSpan(0, 10), "projected": ExonSpan(20, 30)}
        instances = {"local": _instance("local", "A", 0, 10),
                     "projected": _instance("projected", "B", 20, 30)}
        rows = {"A": "A"*20 + "N"*10 + "A"*30, "B": "A"*60}
        observation = _call("A", spans, instances, rows)
        self.assertEqual(observation.kind, "partial")
        self.assertEqual(observation.configurations[0].exons, (ExonSpan(0, 10),))
        self.assertIn(ExonSpan(20, 30), observation.unknown_intervals)

    def test_no_annotation_can_be_absent_when_material_covers_supported_gap(self):
        spans = {"candidate": ExonSpan(12, 20)}
        instances = {"candidate": _instance("candidate", "B", 12, 20)}
        rows = {"A": "A"*12 + "-"*8 + "A"*20, "B": "A"*40}
        observation = _call("A", spans, instances, rows,
            materials=(Material("gap", 12, 20),), presence=(0,), paths={})
        self.assertEqual(observation.kind, "observed")
        self.assertEqual(observation.configurations[0].exons, ())
        self.assertIn("all_candidate_spans_qualified_absent", observation.reasons)

    def test_no_annotation_with_unqualified_sequence_remains_unknown(self):
        spans = {"candidate": ExonSpan(12, 20)}
        instances = {"candidate": _instance("candidate", "B", 12, 20)}
        observation = _call("A", spans, instances, {"A": "A"*40, "B": "A"*40}, paths={})
        self.assertEqual(observation.kind, "unknown")
        self.assertIn("annotation_missing_or_no_physical_exon_spans", observation.reasons)
        self.assertIn("genomic_absence_unqualified", observation.reasons)

    def test_catalogue_roundtrip_preserves_new_marker_and_old_default(self):
        old = Catalogue("f", "old", 10, (), (), observations=())
        new = Catalogue("f", "new", 10, (), (), observation_unit="genomic_exon_spans")
        self.assertEqual(decode_catalogue(encode_catalogue(old)).observation_unit,
                         "transcript_configuration")
        self.assertEqual(decode_catalogue(encode_catalogue(new)).observation_unit,
                         "genomic_exon_spans")

    def test_build_uses_physical_union_independent_of_transcript_partition(self):
        spans = {"a1": ExonSpan(5, 10), "a2": ExonSpan(20, 25),
                 "b1": ExonSpan(5, 10), "b2": ExonSpan(20, 25),
                 "c": ExonSpan(5, 25)}
        instances = {"a1": _instance("a1", "A", 5, 10), "a2": _instance("a2", "A", 20, 25),
                     "b1": _instance("b1", "B", 5, 10), "b2": _instance("b2", "B", 20, 25),
                     "c": _instance("c", "C", 5, 25)}
        rows = {"A": "A"*40, "B": "A"*40}
        rows["C"] = "A"*40
        exons = {key: value for key, value in spans.items()}

        def build(partitioned):
            a_paths = ({"tx1": ("a1", "a2")} if not partitioned else
                       {"tx1": ("a1",), "tx2": ("a2",)})
            loci = (
                SimpleNamespace(family="f", species="A", exons=(instances["a1"], instances["a2"]), paths=a_paths, partial=False),
                SimpleNamespace(family="f", species="B", exons=(instances["b1"], instances["b2"]), paths={"tx1": ("b1", "b2")}, partial=False),
                SimpleNamespace(family="f", species="C", exons=(instances["c"],), paths={"tx1": ("c",)}, partial=False),
            )
            return build_catalogues(FamilyAlignment(loci, rows,
                {species: 0 for species in rows}, {species: tuple(range(40)) for species in rows},
                exons, {}, ()), observation_unit="genomic_exon_spans", anchor_bases=2)

        first, first_rows, _, first_coordinates = build(False)
        second, second_rows, _, second_coordinates = build(True)
        self.assertEqual(len(first), 1)
        self.assertEqual(first[0].observation_unit, "genomic_exon_spans")
        self.assertEqual(first[0].junctions, ((5, 15),))
        self.assertEqual(first[0].junctions, second[0].junctions)
        self.assertEqual(first[0].observations, second[0].observations)
        observation_a = next(o for o in first[0].observations if o.species == "A")
        self.assertEqual(observation_a.configurations[0].exons,
                         (ExonSpan(0, 5), ExonSpan(15, 20)))
        self.assertEqual(len(first_rows), 5)
        self.assertEqual(first_rows, second_rows)
        self.assertEqual(len(first_coordinates), 5)
        self.assertEqual(first_coordinates, second_coordinates)

    def test_unsupported_gap_does_not_connect_units_and_remains_unknown(self):
        rows = {"A": "A"*40, "B": "C"*10 + "-"*22 + "C"*8}
        catalogues = self._build_alignment_catalogues(rows, {
            "A": {"left": ExonSpan(5, 15), "right": ExonSpan(30, 35)},
            "B": {},
        })
        self.assertEqual(len(catalogues), 2)
        for catalogue in catalogues:
            observation = next(o for o in catalogue.observations if o.species == "B")
            self.assertEqual(observation.kind, "unknown")
            self.assertIn("unsupported_alignment_gap", observation.reasons)
            self.assertTrue(observation.unknown_intervals)
            self.assertEqual(catalogue.material, ())
            self.assertNotIn("all_candidate_spans_qualified_absent", observation.reasons)
        self.assertEqual(next(o for o in catalogues[0].observations if o.species == "B").unknown_intervals,
                         (ExonSpan(5, 10),))
        self.assertEqual(next(o for o in catalogues[1].observations if o.species == "B").unknown_intervals,
                         (ExonSpan(0, 2),))

    def test_unsupported_gap_marks_observed_tip_partial_and_retains_exon(self):
        rows = {"A": "A"*40, "B": "A"*8 + "--" + "C"*30,
                "C": "A"*8 + "--" + "C"*30}
        catalogues = self._build_alignment_catalogues(rows, {
            "A": {"a": ExonSpan(5, 15)},
            "B": {"b": ExonSpan(5, 15)},
            "C": {"c": ExonSpan(5, 15)},
        })
        self.assertEqual(len(catalogues), 1)
        observation = next(o for o in catalogues[0].observations if o.species == "B")
        self.assertEqual(observation.kind, "partial")
        self.assertIn(ExonSpan(3, 5), observation.unknown_intervals)
        self.assertEqual(observation.configurations[0].exons, (ExonSpan(0, 10),))

    def test_qualified_shared_deletion_keeps_spanning_exons_in_one_unit(self):
        rows = {"A": "A"*40, "B": "A"*4 + "-"*32 + "A"*4}
        catalogues = self._build_alignment_catalogues(rows, {
            "A": {"left": ExonSpan(5, 15), "right": ExonSpan(25, 35)},
            "B": {},
        })
        self.assertEqual(len(catalogues), 1)
        catalogue = catalogues[0]
        self.assertEqual(catalogue.status, "qualified")
        self.assertEqual(len(catalogue.material), 1)
        self.assertEqual((catalogue.material[0].start, catalogue.material[0].end), (0, 32))
        presence = {o.species: o.material_presence for o in catalogue.observations}
        self.assertEqual(presence, {"A": (1,), "B": (0,)})
        absent = next(o for o in catalogue.observations if o.species == "B")
        self.assertEqual(absent.kind, "observed")
        self.assertEqual(absent.configurations[0].exons, ())

    def test_overlapping_qualified_gaps_remain_unresolved(self):
        rows = {"A": "A"*40, "B": "A"*8 + "-"*17 + "A"*15,
                "C": "A"*15 + "-"*15 + "A"*10}
        catalogues = self._build_alignment_catalogues(rows, {
            "A": {"left": ExonSpan(5, 10), "right": ExonSpan(20, 35)},
            "B": {}, "C": {},
        })
        self.assertEqual(len(catalogues), 1)
        self.assertEqual(catalogues[0].status, "unresolved")
        self.assertIn("overlapping_indel_tracts_unresolved", catalogues[0].reasons)
        self.assertEqual(catalogues[0].material, ())


@unittest.skipUnless(shutil.which("mafft") and shutil.which("minimap2"),
                     "MAFFT and minimap2 are required for raw-input end-to-end coverage")
class GenomicExonDefaultEndToEnd(unittest.TestCase):
    def test_default_raw_conserved_and_annotation_dropout_write_genomic_results(self):
        from intraphy.cli import _dispatch
        from intraphy.commands.parser import build_parser
        from intraphy.inputs.selection import resolve_inputs
        from intraphy.structure.serialization import read_catalogues
        from intraphy.verification.exon_cases import write_exon_example

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for scenario in ("conserved", "annotation_dropout"):
                raw = write_exon_example(root/scenario/"raw", scenario)
                (raw/"truth.json").unlink()
                output = root/scenario/"result"
                args = build_parser().parse_args(["analyze", "--fasta", str(raw), "--gff", str(raw),
                    "--species-tree", str(raw/"species_tree.nwk"), "--output-dir", str(output)])
                args._input_selection = resolve_inputs(args)
                self.assertEqual(args.model, "exon-structure-ctmc")
                _dispatch(args)
                result = json.loads((output/"run_result.json").read_text())
                diagnostics = json.loads((output/"model_diagnostics.json").read_text())
                self.assertEqual(result["model"], "exon-structure-ctmc")
                self.assertTrue((output/"exon_structure_fit.json").exists())
                self.assertTrue((output/"exon_history.json").exists())
                self.assertFalse((output/"dna_presence_history.json").exists())
                catalogues = read_catalogues(output/"exon_configurations.jsonl")
                self.assertTrue(catalogues)
                self.assertTrue(all(c.observation_unit == "genomic_exon_spans" for c in catalogues))
                self.assertEqual(diagnostics["method_scope"],
                                 "constrained_elementary_edit_graph_conditional_composite_likelihood")
                self.assertEqual(diagnostics["ascertainment_correction"], "not_applied")
                if scenario == "annotation_dropout":
                    conditioned = [(c, o) for c in catalogues for o in c.observations
                                   if o.species == "Species_D"
                                   and "annotation_supported_nonexonic_interval" in o.reasons]
                    self.assertTrue(conditioned)
                    for catalogue, observation in conditioned:
                        self.assertEqual(observation.kind, "observed")
                        self.assertTrue(catalogue.spans)
                        self.assertTrue(all(span not in observation.configurations[0].exons
                                            for span in catalogue.spans))
                    # Inputs do not distinguish genuine non-exonic structure from an unmarked
                    # omitted exon annotation; truth is deliberately not supplied to inference.


if __name__ == "__main__":
    unittest.main()
