"""Small deterministic fixtures for disk-staged clustering regressions."""

from itertools import combinations
from pathlib import Path
from unittest.mock import patch

from intraphy.coding.copy_pair import (
    expected_all_copy_pair_keys,
    open_copy_pair_evidence,
    write_copy_pair_evidence,
)
from intraphy.coding_correspondence import CodingProjectionIndex
from support_coding_copy_pairs import (
    INPUT_ID,
    SOURCE_ID,
    coding_rows,
    unchanged_protein_alignment,
)


def staging_fixture():
    """Return equal-length complete proteins and a same-copy exon alternative."""
    occurrences, paths, sequences = [], [], {}
    copies = (
        ("A", "gA1", (("a1", 101, "txA1", 1), ("a2", 301, "txA1", 2),
                       ("a2_alt", 301, "txA_competitor", 2))),
        ("A", "gA2", (("a3", 701, "txA2", 1), ("a4", 901, "txA2", 2))),
        ("B", "gB1", (("b1", 1101, "txB1", 1), ("b2", 1301, "txB1", 2))),
        ("B", "gB2", (("b3", 1701, "txB2", 1), ("b4", 1901, "txB2", 2))),
    )
    for species, gene_copy_id, members in copies:
        for name, start, transcript_id, path_rank in members:
            occurrence, path = coding_rows(name, species, start, "GCT" * 15)
            occurrence["gene_copy_id"] = gene_copy_id
            occurrence["transcript_id"] = transcript_id
            path["gene_copy_id"] = gene_copy_id
            path["transcript_id"] = transcript_id
            path["path_rank"] = str(path_rank)
            path["path_id"] = f"{species}:{gene_copy_id}:{transcript_id}:{name}"
            path["cds_intervals"] = f"{start}-{start + 35}"
            path["cds_length"] = "36"
            occurrences.append(occurrence)
            paths.append(path)
            sequences[name] = "GCT" * 15
            if name == "a1":
                # A second complete transcript reuses a1 with an in-frame
                # offset CDS, then substitutes the overlapping a2_alt exon.
                alternate = {
                    **path,
                    "transcript_id": "txA_competitor",
                    "path_id": f"{species}:{gene_copy_id}:txA_competitor:{name}",
                    "cds_intervals": f"{start + 9}-{start + 44}",
                }
                paths.append(alternate)
    return occurrences, paths, sequences


def terminal_partition_fixture():
    left, left_path = coding_rows("split_left", "A", 101, "GCT" * 9)
    right, right_path = coding_rows("split_right", "A", 301, "GAT" * 9)
    whole, whole_path = coding_rows("whole", "B", 701, "GCT" * 9 + "GAT" * 9)
    occurrences = [left, right, whole]
    paths = [left_path, right_path, whole_path]
    sequences = {
        "split_left": "GCT" * 9,
        "split_right": "GAT" * 9,
        "whole": "GCT" * 9 + "GAT" * 9,
    }
    return occurrences, paths, sequences


def terminal_partition_alignment(records, mode="linsi", threads=1):
    """Align complementary short proteins to opposite whole-protein ends."""
    aligned = {}
    for record_id, protein in dict(records).items():
        if protein == "A" * 9:
            aligned[record_id] = "A" * 9 + "-" * 9
        elif protein == "D" * 9:
            aligned[record_id] = "-" * 9 + "D" * 9
        elif protein == "A" * 9 + "D" * 9:
            aligned[record_id] = protein
        else:
            raise AssertionError(f"unexpected terminal-partition protein: {protein}")
    return aligned


def open_staging_provider(occurrences, paths, sequences, work_root, *,
                          protein_alignment=unchanged_protein_alignment):
    """Build immutable pair artifacts covering the complete copy-pair set."""
    work_root = Path(work_root)
    work_root.mkdir(parents=True, exist_ok=True)
    index = CodingProjectionIndex(occurrences, sequences, paths)
    with patch(
        "intraphy.coding_correspondence.alignment.protein_multiple_alignment",
        side_effect=protein_alignment,
    ):
        projection = index.family_projection("fam")
        copy_keys = sorted({
            (row["family_id"], row["species"], row["gene_copy_id"])
            for row in occurrences
        })
        paths_by_pair = []
        input_identity_by_pair = {}
        for pair_index, (left, right) in enumerate(combinations(copy_keys, 2)):
            records = tuple(index.evidence_for_copy_pair(
                left,
                right,
                projection,
                work_db_path=work_root / f"projection-{pair_index:03d}.sqlite",
            ))
            artifact = work_root / f"pair-{pair_index:03d}.sqlite"
            write_copy_pair_evidence(
                artifact,
                family_id="fam",
                query_copy_key=left,
                target_copy_key=right,
                source_identity=SOURCE_ID,
                input_identity=INPUT_ID,
                expected_pairs=tuple(key for key, _ in records),
                records=iter(records),
            )
            paths_by_pair.append(artifact)
            input_identity_by_pair[(left, right)] = INPUT_ID
    provider = open_copy_pair_evidence(
        paths_by_pair,
        expected_source_identity=SOURCE_ID,
        expected_input_identity=input_identity_by_pair,
        expected_pairs=expected_all_copy_pair_keys(occurrences),
    )
    return provider
