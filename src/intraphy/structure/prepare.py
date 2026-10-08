"""Prepare V19 exon configuration evidence using existing FASTA/GFF extraction."""
from __future__ import annotations
from collections import defaultdict
from pathlib import Path
from ..storage.tabular import write_tsv
from .native import load_native
from .consequences import native_cds
from .alignment import align_family
from .corroboration import check_coding_projection
from .build import build_catalogues
from .serialization import write_catalogues, write_json


def _table(path, rows, fields):
    write_tsv(path, rows, list(dict.fromkeys(fields + [k for row in rows for k in row])))


def prepare_configurations(input_dir: str | Path, output_dir: str | Path, *,
                           timeout=600, max_locus_bases=100000,
                           minimum_identity=.7, anchor_bases=12, anchor_identity=.8,
                           observation_unit="transcript_configuration", threads=1,
                           alignment_evidence_dir=None):
    if type(threads) is not int or threads < 1:
        raise ValueError("threads must be a positive integer")
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    grouped = defaultdict(list)
    for locus in load_native(input_dir):
        grouped[locus.family].append(locus)
    catalogues, correspondences, candidates, coords, reports, coding = [], [], [], [], [], []
    coding_projection_records = 0
    reused_families = []
    for number, (family, loci) in enumerate(sorted(grouped.items()), 1):
        # Filesystem names do not use untrusted gene/annotation strings.
        directory = out/"alignment_evidence"/f"family_{number:05d}"
        for locus in loci:
            for tx in sorted(locus.paths):
                coding.append({"family_id": family, "species": locus.species, "transcript": tx,
                               **native_cds(locus, tx)})
        try:
            alignment = align_family(tuple(loci), directory, timeout=timeout,
                                     max_bases=max_locus_bases, threads=threads,
                                     evidence_dir=(Path(alignment_evidence_dir) / f"family_{number:05d}"
                                                   if alignment_evidence_dir is not None else None))
            if alignment.reuse_provenance is not None:
                reused_families.append({"family_id": family,
                    "source_family_directory": f"family_{number:05d}",
                    "commands": alignment.reuse_provenance})
            alignment, protein = check_coding_projection(alignment, input_dir)
            coding_projection_records += len(protein)
            cs, matches, predictions, coordinates = build_catalogues(alignment,
                minimum_identity=minimum_identity, anchor_bases=anchor_bases,
                anchor_identity=anchor_identity, observation_unit=observation_unit)
            _table(directory/"nucleotide_search.tsv", alignment.matches, ["query_exon", "target_species", "coverage", "identity"])
            _table(directory/"coding_projection_check.tsv", protein, ["query_exon", "target_exon", "checked_bases", "agreeing_bases", "status"])
            catalogues.extend(cs)
            correspondences.extend(matches)
            candidates.extend(predictions)
            coords.extend(coordinates)
            # Save only actual paired genomic blocks for read-only ribbon plots.
            mapping = []
            for locus in alignment.loci:
                cols = alignment.columns[locus.species]
                mapping.append({"species": locus.species, "source_contig": locus.contig,
                    "strand": locus.strand, "source_positions0": [locus.source_position(i+alignment.offsets[locus.species]) for i in range(len(cols))],
                    "alignment_columns0": cols})
            write_json(directory/"coordinate_map.json", mapping)
            reports.append({"family_id": family, "status": "prepared", "units": len(cs),
                            "qualified_units": sum(c.status == "qualified" for c in cs), "reason": "", "evidence_directory": str(directory.relative_to(out))})
        except ValueError as exc:
            # Scientific/resource qualification is not a successful negative call.
            # External executable failure propagates rather than disappearing.
            if not any(str(exc).startswith(prefix) for prefix in ("locus_alignment_budget_exceeded", "no_annotated_exons")):
                raise
            reports.append({"family_id": family, "status": "unresolved", "units": 0,
                            "qualified_units": 0, "reason": str(exc), "evidence_directory": str(directory.relative_to(out))})
    if not grouped:
        raise ValueError("Prepared inputs contain no identifiable gene loci")
    write_json(out/"native_cds_consequences.json", coding)
    write_catalogues(out/"exon_configurations.jsonl", catalogues)
    _table(out/"exon_correspondence.tsv", correspondences, ["family_id", "unit_id", "species", "exon_id", "status"])
    _table(out/"annotation_structure_candidates.tsv", candidates, ["family_id", "unit_id", "species", "source_species", "source_transcript", "status"])
    _table(out/"exon_coordinates.tsv", coords, ["family_id", "unit_id", "species", "exon_id", "coordinate_system"])
    _table(out/"exon_preparation_summary.tsv", reports, ["family_id", "status", "units", "qualified_units", "reason", "evidence_directory"])
    policy = {"observation_unit": observation_unit,
        "minimum_nucleotide_identity": minimum_identity,
        "gap_anchor_bases": anchor_bases, "gap_anchor_identity": anchor_identity,
        "protein_projection_agreement": .95, "thresholds_biologically_calibrated": False,
        "genomic_alignment": "MAFFT --auto --thread N --threadit 0; threadtb left at MAFFT default",
        "nucleotide_search": "minimap2 -t N",
        "search_window": "entire_supplied_or_extracted_locus; no annotation-based cropping",
        "annotation_alternatives": ("not_used_for_genomic_exon_spans" if observation_unit == "genomic_exon_spans"
                                     else "atomic_source_configuration; paired_flanks; exact_4base_changed_cut_context"),
        "protein_projection_check": {"status": "records_used" if coding_projection_records else "not_available",
                                     "records": coding_projection_records,
                                     "applies_only_to_retained_segment_matches": True},
        "native_cds_consequences": "annotation_descriptive_only_not_a_state_or_rate_observation",
        "annotation_conditioned_nonexonic_evidence": ({
            "scope": "provided_annotation_paths_only",
            "requirements": "adjacent_exons_on_one_provided_path; no_provided_exon_overlaps_candidate; retained_homologous_DNA_with_paired_candidate_and_boundary_anchors; localized_partial_windows_avoid_candidate_and_boundaries",
            "transcript_path_role": "annotation_context_only; transcript_ids_do_not_define_states_or_RNA_usage",
            "thresholds_biologically_calibrated": False,
        } if observation_unit == "genomic_exon_spans" else None),
        "boundary_context_is_splice_function_evidence": False,
        "copy_and_orientation_audit": "minimap2 plus exact short-exon search",
        "discovery": "annotation_discovered", "reference_is_ancestor": False,
        "transcript_ids_are_provenance_only": observation_unit == "genomic_exon_spans",
        "required_inputs": "genomic sequence, exon annotations and species tree; no RNA measurements"}
    policy["alignment_evidence_reuse"] = (None if alignment_evidence_dir is None else {
        "source_directory": str(Path(alignment_evidence_dir).resolve()),
        "reuse_only": True,
        "new_external_commands_executed": False,
        "families": reused_families})
    write_json(out/"exon_evidence_policy.json", policy)
    return tuple(catalogues), reports
