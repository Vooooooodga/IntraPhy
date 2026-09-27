"""preparation / derive: explicit implementation ownership."""
from __future__ import annotations

from intraphy.aligners.short import _validate_max_dp_cells
from intraphy.aligners.runner import available_alignment_backends
from intraphy.aligners.types import MAX_INTERNAL_DP_CELLS
from intraphy.mapping.clustering import cluster_segments
from intraphy.mapping.fields import MATCH_FIELDS
from intraphy.mapping.fields import SEGMENT_FIELDS
from intraphy.mapping.flank_context import _gene_locus_records
from intraphy.mapping.policies import _species_tree_distances
from intraphy.preparation.copy_context import make_adjacencies
from intraphy.preparation.copy_context import make_copy_context
from intraphy.preparation.copy_context import make_copy_relationships
from intraphy.storage.fasta import parse_fasta
from intraphy.storage.tabular import read_tsv
from intraphy.storage.tabular import write_tsv
from pathlib import Path
import csv


def derive_tables(input_dir, output_dir=None, identity_threshold=0.7, distance_table=None, aligner="mafft", threads=1, min_size_ratio=0.25, context_aligner="minimap2", coding_msa_mode="linsi", short_context_max_length=300, *, protein_evidence_provider=None, pair_scoring_executor="thread", short_alignment_max_dp_cells=MAX_INTERNAL_DP_CELLS, match_staging_path=None, progress_callback=None, collect_matches=True):
    short_alignment_max_dp_cells = _validate_max_dp_cells(short_alignment_max_dp_cells)
    input_dir = Path(input_dir)
    output_dir = Path(output_dir or input_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    occurrences = read_tsv(input_dir / "segment_occurrences.tsv", SEGMENT_FIELDS)
    transcript_paths = read_tsv(input_dir / "transcript_paths.tsv", optional=True)
    raw_features = read_tsv(input_dir / "raw_gene_features.tsv", optional=True)
    seqs = parse_fasta(input_dir / "segment_sequences.fasta")
    gene_loci = _gene_locus_records(input_dir)
    species_distances = _species_tree_distances(input_dir / "species_tree.tsv")
    match_path = output_dir / "segment_matches.tsv"
    with match_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, MATCH_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        homology, matches = cluster_segments(
            occurrences, seqs, identity_threshold, distance_table, aligner=aligner,
            threads=threads, min_size_ratio=min_size_ratio, species_distances=species_distances,
            match_writer=lambda row: writer.writerow({field: row.get(field, "NA") for field in MATCH_FIELDS}),
            context_aligner=context_aligner,
            transcript_paths=transcript_paths,
            raw_features=raw_features,
            coding_msa_mode=coding_msa_mode,
            short_context_max_length=short_context_max_length,
            gene_loci=gene_loci,
            protein_evidence_provider=protein_evidence_provider,
            pair_scoring_executor=pair_scoring_executor,
            short_alignment_max_dp_cells=short_alignment_max_dp_cells,
            match_staging_path=match_staging_path,
            progress_callback=progress_callback,
            collect_matches=collect_matches,
        )
    write_tsv(output_dir / "segment_homology.tsv", homology, ["homology_id", "occurrence_id", "support_type", "confidence", "source_label"])
    backend_rows = []
    for row in available_alignment_backends():
        selected_exon = row["aligner"] == aligner
        selected_context = row["aligner"] == context_aligner
        backend_rows.append(
            {
                **row,
                "selected": int(selected_exon or selected_context),
                "selected_exon": int(selected_exon),
                "selected_context": int(selected_context),
                "alignment_mode": (
                    "overlap_projection" if selected_exon and aligner in {"mafft", "auto"}
                    else "local" if selected_exon or selected_context else "NA"
                ),
                "threads": threads,
                "min_size_ratio": f"{min_size_ratio:.6g}",
                "short_context_max_length": int(short_context_max_length),
                "short_alignment_dp_cell_budget": (
                    "unlimited" if short_alignment_max_dp_cells is None
                    else short_alignment_max_dp_cells
                ),
                "coding_msa_mode": coding_msa_mode,
                "notes": (
                    f"{row['notes']}; min_size_ratio is restricted to non-exon-like prefiltering"
                    if selected_exon or selected_context
                    else row["notes"]
                ) + "; short-alignment DP budget rejects the complete alignment when exceeded; DP states are not truncated",
            }
        )
    write_tsv(
        output_dir / "alignment_backend_report.tsv",
        backend_rows,
        [
            "aligner",
            "available",
            "selected",
            "threads",
            "min_size_ratio",
            "short_context_max_length",
            "short_alignment_dp_cell_budget",
            "coding_msa_mode",
            "notes",
            "selected_exon",
            "selected_context",
            "alignment_mode",
        ],
    )
    write_tsv(output_dir / "physical_adjacencies.tsv", make_adjacencies(occurrences), ["adjacency_id", "family_id", "species", "gene_copy_id", "left_occurrence_id", "right_occurrence_id", "adjacency_status"])
    write_tsv(output_dir / "copy_context.tsv", make_copy_context(occurrences), ["family_id", "species", "gene_copy_id", "copy_class", "copy_subclass", "copy_span"])
    write_tsv(output_dir / "copy_relationships.tsv", make_copy_relationships(occurrences), ["family_id", "species", "query_copy_id", "subject_copy_id", "relationship_class", "synteny_score", "distance_bp", "evidence"])
    if not (output_dir / "sequence_synteny_evidence.tsv").exists():
        evidence = []
        for row in occurrences:
            evidence.append(
                {
                    "evidence_id": f"ev_{row['occurrence_id']}",
                    "family_id": row["family_id"],
                    "species": row["species"],
                    "gene_copy_id": row["gene_copy_id"],
                    "homology_id": "NA",
                    "annotation_status": "annotated",
                    "evidence_status": "supports_annotation",
                    "inferred_role": row["role"],
                    "contig": row.get("contig", "NA"),
                    "start": row.get("start", "NA"),
                    "end": row.get("end", "NA"),
                    "strand": row.get("strand", "NA"),
                    "sequence_score": "1.0",
                    "left_synteny_score": "1.0",
                    "right_synteny_score": "1.0",
                    "splice_motif_score": row.get("splice_motif_score", "0.5"),
                    "phase_compatibility": "compatible" if row.get("phase") not in {".", "NA", ""} else "unknown",
                    "inferred_event": "annotated_segment",
                    "frame_status": row.get("frame_status", "unknown"),
                }
            )
        write_tsv(
            output_dir / "sequence_synteny_evidence.tsv",
            evidence,
            [
                "evidence_id",
                "family_id",
                "species",
                "gene_copy_id",
                "homology_id",
                "annotation_status",
                "evidence_status",
                "inferred_role",
                "contig",
                "start",
                "end",
                "strand",
                "sequence_score",
                "left_synteny_score",
                "right_synteny_score",
                "splice_motif_score",
                "phase_compatibility",
                "inferred_event",
                "frame_status",
            ],
        )
    return homology, matches
