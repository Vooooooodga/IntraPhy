"""coding_correspondence: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from collections import OrderedDict, defaultdict
from intraphy.preparation.features import FeatureHierarchy
from intraphy import alignment
from intraphy.coding.projection import _aligned_occurrence_pairs
from intraphy.coding.transcripts import _copy_key
from intraphy.coding.transcripts import build_coding_transcript
from intraphy.coding.types import CodingResidueProjection
from intraphy.coding.types import FamilyCodingProjection


class CodingProjectionIndex:
    """Family MSA index with transcript-specific CDS and genome projections."""

    MAX_CACHED_PAIRS = 32
    MAX_CACHED_BASE_PAIRS = 1_000_000

    def __init__(
        self, occurrences, sequences, transcript_paths, raw_features=(), threads=1,
        msa_mode="linsi",
    ):
        self.threads = threads
        self.msa_mode = msa_mode
        self.cache = OrderedDict()
        self.cached_bases = 0
        self.transcripts = {}
        self.by_occurrence = defaultdict(list)
        self.family_alignments = {}
        self.occurrences = {row["occurrence_id"]: row for row in occurrences}
        self.occurrences_by_copy = defaultdict(list)
        for row in occurrences:
            self.occurrences_by_copy[_copy_key(row)].append(row["occurrence_id"])
        by_copy = defaultdict(list)
        for row in raw_features:
            if row.get("ownership") == "target_gene_descendant":
                by_copy[_copy_key(row)].append(row)
        by_transcript = defaultdict(list)
        for row in transcript_paths:
            by_transcript[(*_copy_key(row), row["transcript_id"])].append(row)
        by_id = self.occurrences
        feature_indexes = {key: FeatureHierarchy(rows) for key, rows in by_copy.items()}
        for key, paths in sorted(by_transcript.items()):
            transcript = build_coding_transcript(key, paths, by_id, sequences, by_copy[key[:-1]], feature_index=feature_indexes.get(key[:-1]))
            self.transcripts[key] = transcript
            for occurrence_id in {row["occurrence_id"] for row in paths}:
                self.by_occurrence[occurrence_id].append(key)

    def _family_alignment(self, family_id):
        if family_id in self.family_alignments:
            return self.family_alignments[family_id]
        family_transcripts = [
            transcript for key, transcript in sorted(self.transcripts.items())
            if key[0] == family_id and not transcript.unavailable_reason and transcript.protein
        ]
        aliases_by_protein = defaultdict(list)
        for transcript in family_transcripts:
            aliases_by_protein[transcript.protein].append(transcript.key)
        records = []
        aliases_by_record = {}
        record_by_transcript = {}
        for index, protein in enumerate(sorted(aliases_by_protein), 1):
            record_id = f"coding_{index:06d}"
            aliases = tuple(sorted(aliases_by_protein[protein]))
            records.append((record_id, protein))
            aliases_by_record[record_id] = aliases
            for transcript_key in aliases:
                record_by_transcript[transcript_key] = record_id
        aligned_records = dict(alignment.protein_multiple_alignment(
            tuple(records), mode=self.msa_mode, threads=self.threads,
        )) if records else {}
        expected_ids = {record_id for record_id, _protein in records}
        if set(aligned_records) != expected_ids:
            raise alignment.AlignmentBackendError("family protein MSA returned a different record set")
        lengths = {len(sequence) for sequence in aligned_records.values()}
        if len(lengths) > 1:
            raise alignment.AlignmentBackendError("family protein MSA returned unequal alignment lengths")
        residue_columns = {}
        for record_id, protein in records:
            aligned = aligned_records[record_id]
            if aligned.replace("-", "") != protein:
                raise alignment.AlignmentBackendError("family protein MSA did not preserve input residues")
            residue_columns[record_id] = tuple(
                column0 for column0, amino_acid in enumerate(aligned) if amino_acid != "-"
            )
        projection = FamilyCodingProjection(
            family_id=family_id, mode=self.msa_mode,
            aligned_records=aligned_records,
            aliases_by_record=aliases_by_record,
            record_by_transcript=record_by_transcript,
            residue_columns=residue_columns,
        )
        self.family_alignments[family_id] = projection
        return projection

    def use_family_projection(self, projection):
        """Install a strictly validated full-family MSA for this two-copy index."""
        if not isinstance(projection, FamilyCodingProjection):
            raise TypeError("FamilyCodingProjection required")
        if projection.mode != self.msa_mode:
            raise ValueError("family MSA mode does not match coding index")
        if any(key[0] != projection.family_id for key in projection.record_by_transcript):
            raise ValueError("family MSA contains an alias from another family")
        for key, transcript in self.transcripts.items():
            if key[0] != projection.family_id or transcript.unavailable_reason or not transcript.protein:
                continue
            record_id = projection.record_by_transcript.get(key)
            if record_id is None:
                raise ValueError(f"family MSA lacks usable transcript alias {key}")
            aligned = projection.aligned_for(key)
            if aligned.replace("-", "") != transcript.protein:
                raise ValueError(f"family MSA residues disagree with transcript {key}")
            columns = tuple(i for i, aa in enumerate(aligned) if aa != "-")
            if projection.columns_for(key) != columns:
                raise ValueError(f"family MSA residue columns disagree with transcript {key}")
        previous = self.family_alignments.get(projection.family_id)
        if previous is not projection:
            for key in tuple(self.cache):
                if key[0][0] == projection.family_id:
                    _pairs, _columns, _left, _right, size = self.cache.pop(key)
                    self.cached_bases -= size
        self.family_alignments[projection.family_id] = projection
        return projection

    def _canonical_pair_projection(self, query_key, target_key, family_projection=None, *, cache=True):
        """Return one canonical projection, optionally participating in the legacy LRU."""
        key = tuple(sorted((query_key, target_key)))
        if not cache and key in self.cache:
            _pairs, _columns, _left, _right, old_size = self.cache.pop(key)
            self.cached_bases -= old_size
        if cache and key in self.cache:
            self.cache.move_to_end(key)
            pairs, known_columns, aligned_left, aligned_right, _size = self.cache[key]
        else:
            left, right = (self.transcripts[item] for item in key)
            if left.key[0] != right.key[0]:
                return key, {}, (), "", ""
            family = family_projection or self._family_alignment(left.key[0])
            if family_projection is not None and self.family_alignments.get(left.key[0]) is not family_projection:
                self.use_family_projection(family_projection)
            aligned_left = family.aligned_for(left.key)
            aligned_right = family.aligned_for(right.key)
            pairs, known_columns = _aligned_occurrence_pairs(
                left, right, aligned_left, aligned_right,
            )
            size = sum(len(record["positions0"]) for record in pairs.values())
            if cache and size <= self.MAX_CACHED_BASE_PAIRS:
                while self.cache and (len(self.cache) >= self.MAX_CACHED_PAIRS or
                                      self.cached_bases + size > self.MAX_CACHED_BASE_PAIRS):
                    _old_key, (_old_pairs, _old_columns, _left, _right, old_size) = self.cache.popitem(last=False)
                    self.cached_bases -= old_size
                self.cache[key] = (pairs, known_columns, aligned_left, aligned_right, size)
                self.cached_bases += size
        return key, pairs, known_columns, aligned_left, aligned_right

    def _pair(self, query_key, target_key):
        key, pairs, known_columns, aligned_left, aligned_right = self._canonical_pair_projection(query_key, target_key)
        if not aligned_left and not aligned_right:
            return {}, (), "", "", False
        inverted = query_key != key[0]
        if not inverted:
            return pairs, known_columns, aligned_left, aligned_right, False
        reversed_pairs = {}
        for (left_occurrence, right_occurrence), record in pairs.items():
            reversed_pairs[(right_occurrence, left_occurrence)] = {
                **record,
                "positions0": [(right0, left0) for left0, right0 in record["positions0"]],
            }
        reversed_columns = tuple(
            (column0, right_aa, left_aa, right_codon, left_codon)
            for column0, left_aa, right_aa, left_codon, right_codon in known_columns
        )
        return reversed_pairs, reversed_columns, aligned_right, aligned_left, True

    def project(self, occurrence_id):
        """Return residue-to-MSA-to-CDS/genome projections for one occurrence."""
        projections = []
        for transcript_key in sorted(self.by_occurrence.get(occurrence_id, ())):
            transcript = self.transcripts[transcript_key]
            if transcript.unavailable_reason:
                continue
            family = self._family_alignment(transcript_key[0])
            columns = family.columns_for(transcript_key)
            for residue_index0, codon in enumerate(transcript.codon_sources):
                if not any(base.occurrence_id == occurrence_id for base in codon):
                    continue
                projections.append(CodingResidueProjection(
                    transcript_key=transcript_key,
                    residue_index0=residue_index0,
                    msa_column0=columns[residue_index0],
                    amino_acid=transcript.protein[residue_index0],
                    cds_bases=codon,
                ))
        return tuple(projections)

    def family_projection(self, family_id):
        """Expose the run-local family MSA and its complete alias mapping."""
        return self._family_alignment(family_id)

    def candidates(self, query_occurrence, target_copy_id):
        """Return local coding evidence for every occurrence in a target copy."""
        query = self.occurrences.get(query_occurrence)
        if query is None:
            return tuple()
        if isinstance(target_copy_id, tuple):
            copy_key = target_copy_id
        else:
            matches = [
                key for key in self.occurrences_by_copy
                if key[0] == query.get("family_id") and key[2] == target_copy_id
            ]
            if len(matches) != 1:
                return tuple()
            copy_key = matches[0]
        return tuple(
            {"target_occurrence_id": target, **self.evidence(query_occurrence, target)}
            for target in sorted(self.occurrences_by_copy.get(copy_key, ()))
        )

    def _competing_occurrences(
        self, query_key, target_key, query_occurrence, target_occurrence,
        query_positions0, target_positions0,
    ):
        """Find alternative occurrence placements supported by compatible paths."""
        competitors = set()
        query_copy = query_key[:3]
        target_copy = target_key[:3]
        target_transcripts = [
            key for key, transcript in self.transcripts.items()
            if key[:3] == target_copy and not transcript.unavailable_reason
        ]
        query_transcripts = [
            key for key, transcript in self.transcripts.items()
            if key[:3] == query_copy and not transcript.unavailable_reason
        ]
        for alternative_target_key in target_transcripts:
            pairs, _columns, _aligned_query, _aligned_target, _inverted = self._pair(
                query_key, alternative_target_key,
            )
            for (other_query, other_target), record in pairs.items():
                if other_query != query_occurrence or other_target == target_occurrence:
                    continue
                other_query_positions0 = {
                    query0 for query0, _target0 in record["positions0"]
                }
                if query_positions0 & other_query_positions0:
                    competitors.add(f"target:{other_target}")
        for alternative_query_key in query_transcripts:
            pairs, _columns, _aligned_query, _aligned_target, _inverted = self._pair(
                alternative_query_key, target_key,
            )
            for (other_query, other_target), record in pairs.items():
                if other_target != target_occurrence or other_query == query_occurrence:
                    continue
                other_target_positions0 = {
                    target0 for _query0, target0 in record["positions0"]
                }
                if target_positions0 & other_target_positions0:
                    competitors.add(f"query:{other_query}")
        return tuple(sorted(competitors))

    def evidence(self, query_occurrence, target_occurrence):
        from intraphy.coding.projection_evidence import evidence
        return evidence(self, query_occurrence, target_occurrence)

    def evidence_for_copy_pair(self, query_copy_key, target_copy_key, family_projection=None,
                               *, work_db_path):
        from intraphy.coding.copy_pair_projection import evidence_for_copy_pair
        return evidence_for_copy_pair(
            self, query_copy_key, target_copy_key, family_projection,
            work_db_path=work_db_path,
        )



# Public entry points; implementations have a single owner.
from intraphy.coding.projection import (
    _blosum62_score,
    _aligned_occurrence_pairs,
    _coordinate_blocks0,
    _anchor_metrics,
    _candidate_json,
)
from intraphy.coding.transcripts import (
    _copy_key,
    _transcript_features,
    _translation_reason,
    _synthetic_transcript_features,
    build_coding_transcript,
)
from intraphy.coding.types import (
    _STANDARD_CODE,
    _BLOSUM62,
    _EMPTY,
    _ANCHOR_WINDOW,
    _MIN_ANCHOR_PAIRS,
    CodingBaseProjection,
    CodingResidueProjection,
    CodingTranscript,
    FamilyCodingProjection,
    CodingProjectionCandidate,
    CodingProjectionCandidateSet,
)
import json
from intraphy import alignment
