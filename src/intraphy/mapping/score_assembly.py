"""Build exact public match rows from completed nucleotide/protein evidence."""
from __future__ import annotations

from intraphy.mapping.candidate_coordinates import _coordinate_block0
from intraphy.mapping.candidate_coordinates import _format_alignment_blocks
from intraphy.mapping.candidate_coordinates import _genomic_matched_blocks
from intraphy.mapping.fields import EXON_LIKE_ROLES
from intraphy.mapping.match_context import pair_threshold
from intraphy.mapping.match_records import _match_row
from intraphy.mapping.match_records import _projection_record
from intraphy.mapping.policies import _candidate_sequence_accepted
from intraphy.mapping.policies import roles_compatible
from intraphy.mapping.policies import sequence_supported_mapping
from intraphy.coordinates import ClosedInterval1
from intraphy.coordinates import Interval0


def assemble_scored_match(
    left, right, evidence, prefilter_status, identity_threshold, distance_lookup,
    *, protein_index, protein_evidence_provider, aligner, match_index,
):
    """Qualify scored evidence and assemble its legacy-compatible match row."""
    threshold, distance_class = pair_threshold(
        left, right, identity_threshold, distance_lookup,
    )
    score = evidence["sequence_score"]
    compatible = roles_compatible(left, right)
    candidate_records = evidence.get("candidate_records", ())
    short_context = evidence.get("short_context_route") in {
        "feature_bounded_candidate", "anchor_bounded_local",
    }
    for record in candidate_records:
        if record.get("source") not in {None, "", "NA"}:
            record.setdefault("alignment_source", record["source"])
        record["source"] = "nucleotide_alignment"
        record.setdefault("score_scheme", evidence.get("score_scheme", "unspecified"))
        record["accepted"] = int(
            _candidate_sequence_accepted(record, threshold, short_context)
        )
        record["acceptance_threshold"] = f"{threshold:.6g}"
    sequence_ok = any(record.get("accepted") == 1 for record in candidate_records)
    if not candidate_records:
        sequence_ok = sequence_supported_mapping(evidence, threshold)
    mapped = compatible and sequence_ok
    evidence["dna_match_status"] = (
        "mapped" if mapped else
        prefilter_status if prefilter_status != "aligned_candidate" else "low_similarity"
    )
    exon_pair = left.get("role") in EXON_LIKE_ROLES and right.get("role") in EXON_LIKE_ROLES
    evidence["dna_candidate_assessments"] = [
        dict(record) for record in candidate_records
    ]
    if exon_pair and (protein_index is not None or protein_evidence_provider is not None):
        protein_evidence = (
            protein_evidence_provider.get(left["occurrence_id"], right["occurrence_id"])
            if protein_evidence_provider is not None
            else protein_index.evidence(left["occurrence_id"], right["occurrence_id"])
        )
        evidence.update(protein_evidence)
        protein_hard = bool(evidence.get("protein_hard_observation_eligible"))
        protein_position = bool(evidence.get("protein_position_eligible"))
        if mapped and protein_hard:
            evidence["correspondence_basis"] = "DNA_and_annotated_CDS_protein"
        elif compatible and protein_hard:
            mapped = True
            score = (
                0.70 * float(evidence["protein_aa_identity"])
                + 0.30 * min(
                    float(evidence["protein_query_cds_coverage"]),
                    float(evidence["protein_target_cds_coverage"]),
                )
            )
            evidence["correspondence_basis"] = "annotated_CDS_protein"
        if mapped and protein_position:
            protein_blocks = tuple(
                _coordinate_block0(block)
                for block in evidence.get("protein_projected_blocks", ())
            )
            evidence["projected_reference_blocks"] = _format_alignment_blocks(protein_blocks)
            evidence["matched_blocks"] = evidence["projected_reference_blocks"]
            if protein_blocks:
                query_interval = Interval0(
                    min(block.query.start0 for block in protein_blocks),
                    max(block.query.end0 for block in protein_blocks),
                )
                target_interval = Interval0(
                    min(block.target.start0 for block in protein_blocks),
                    max(block.target.end0 for block in protein_blocks),
                )
                public_query = ClosedInterval1.from_interval0(query_interval)
                public_target = ClosedInterval1.from_interval0(target_interval)
                evidence["query_alignment_start"] = public_query.start
                evidence["query_alignment_end"] = public_query.end
                evidence["target_alignment_start"] = public_target.start
                evidence["target_alignment_end"] = public_target.end
                evidence["query_genomic_matched_blocks"] = _genomic_matched_blocks(
                    left, protein_blocks, "query",
                )
                evidence["subject_genomic_matched_blocks"] = _genomic_matched_blocks(
                    right, protein_blocks, "subject",
                )
            evidence["score_scheme"] = "blosum62_cds_projection"
            evidence["raw_alignment_score"] = evidence.get("protein_blosum62_score", "NA")
            evidence["alignment_backend"] = "family_protein_msa_projection"
            evidence["alignment_mode"] = "coding_projection"
            evidence["alignment_meaning"] = (
                "family protein MSA projected through transcript codons to CDS bases"
            )
            evidence["alignment_strand"] = "+"
            evidence["candidate_records"] = []
            evidence["candidate_accepted"] = True
            evidence["enumeration_complete"] = (
                evidence.get("protein_candidate_details") not in {None, "", "NA"}
            )
            evidence["candidate_enumeration_status"] = (
                "complete" if evidence["enumeration_complete"] else "incomplete"
            )
            evidence["incomplete_reason"] = (
                "NA" if evidence["enumeration_complete"]
                else "protein_candidate_set_unavailable"
            )
        elif evidence.get("protein_candidate_evidence_available"):
            evidence["incomplete_reason"] = "protein_candidate_without_hard_coordinates"
    elif not mapped and exon_pair and protein_index is None:
        evidence["protein_status"] = (
            "unavailable" if aligner in {"mafft", "auto"} else "disabled_backend"
        )
        evidence["protein_unavailable_reason"] = (
            "no_CDS_transcript_path" if aligner in {"mafft", "auto"}
            else "MAFFT_exon_backend_required"
        )

    status = "mapped" if mapped else (
        prefilter_status if prefilter_status != "aligned_candidate" else "low_similarity"
    )
    projection_facts = ()
    if mapped:
        projection_facts = (
            ((left["occurrence_id"], right["occurrence_id"]), _projection_record(evidence, "target")),
            ((right["occurrence_id"], left["occurrence_id"]), _projection_record(evidence, "query")),
        )
    row = _match_row(
        left, right, evidence, score, threshold, distance_class, status, match_index,
    )
    return score, status, mapped, row, projection_facts
