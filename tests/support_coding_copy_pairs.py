"""Shared bounded fixtures for copy-pair cache regression tests."""

import tempfile
from pathlib import Path
from unittest.mock import patch

from intraphy.coding.copy_pair import open_copy_pair_evidence, write_copy_pair_evidence
from intraphy.coding_correspondence import CodingProjectionIndex
from support_coding_correspondence_coding_correspondence import (
    CodingCorrespondenceTestsSupport,
    coding_rows,
    unchanged_protein_alignment,
    weak_dna_evidence,
)


SOURCE_ID = "synthetic-source-v1"
INPUT_ID = "synthetic-input-v1"


def copy_key(species):
    return ("fam", species, f"g{species}")


def indexed_records(index, projection, work_path, query_species="A", target_species="B"):
    return tuple(index.evidence_for_copy_pair(
        copy_key(query_species), copy_key(target_species), projection,
        work_db_path=work_path,
    ))


def strong_dna_evidence(left, right, sequences, context, **kwargs):
    evidence = weak_dna_evidence(left, right, sequences, context, **kwargs)
    evidence.update(
        alignment_score=0.95, coverage_score=1.0, sequence_score=0.965,
        total_score=evidence["total_score"] + 0.34 * 0.55,
    )
    return evidence


def disjoint_protein_alignment(records, mode="linsi", threads=1):
    aligned = {}
    for record_id, protein in dict(records).items():
        if protein and set(protein) == {"A"}:
            aligned[record_id] = protein + "-" * 12
        else:
            aligned[record_id] = "-" * 12 + protein
    return aligned


class CodingCopyPairTestSupport(CodingCorrespondenceTestsSupport):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)

    def _provider_for(self, occurrences, paths, sequences, *, omit=None):
        index = CodingProjectionIndex(occurrences, sequences, paths)
        with patch(
            "intraphy.coding_correspondence.alignment.protein_multiple_alignment",
            side_effect=unchanged_protein_alignment,
        ):
            projection = index._family_alignment("fam")
            records = indexed_records(index, projection, self.root / "provider-work.sqlite")
        if omit is not None:
            records = tuple(row for row in records if row[0] != omit)
        expected_pairs = tuple(key for key, _ in records)
        path = self.root / ("provider-" + str(len(list(self.root.iterdir()))) + ".sqlite")
        write_copy_pair_evidence(
            path, family_id="fam", query_copy_key=copy_key("A"),
            target_copy_key=copy_key("B"), source_identity=SOURCE_ID,
            input_identity=INPUT_ID, expected_pairs=expected_pairs,
            records=iter(records),
        )
        return open_copy_pair_evidence(
            [path], expected_source_identity=SOURCE_ID,
            expected_input_identity=INPUT_ID, expected_pairs=expected_pairs,
        )
