"""Copy-pair projection equivalence, coordinate, and cache contracts."""

import json
import unittest
from dataclasses import replace
from unittest.mock import patch

from intraphy.coding.copy_pair import read_family_projection, write_family_projection
from intraphy.coding_correspondence import CodingProjectionIndex
from support_coding_copy_pairs import (
    CodingCopyPairTestSupport,
    copy_key,
    disjoint_protein_alignment,
    indexed_records,
    unchanged_protein_alignment,
    coding_rows,
)


class CodingCopyPairProjectionTests(CodingCopyPairTestSupport, unittest.TestCase):
    def test_family_projection_round_trip_and_identity_contract(self):
        occurrences, paths, sequences = self.split_fixture()
        index = CodingProjectionIndex(occurrences, sequences, paths)
        with patch(
            "intraphy.coding_correspondence.alignment.protein_multiple_alignment",
            side_effect=unchanged_protein_alignment,
        ):
            projection = index._family_alignment("fam")
        path = self.root / "family.json"
        write_family_projection(path, projection, source_identity="source", input_identity="input")
        replay = read_family_projection(
            path, expected_family_id="fam", expected_mode="linsi",
            source_identity="source", input_identity="input",
        )
        self.assertEqual(replay, projection)
        self.assertEqual(replay.aligned_records, projection.aligned_records)
        self.assertEqual(replay.aliases_by_record, projection.aliases_by_record)
        self.assertEqual(replay.record_by_transcript, projection.record_by_transcript)
        self.assertEqual(replay.residue_columns, projection.residue_columns)
        for kwargs in (
            {"expected_family_id": "other", "expected_mode": "linsi", "source_identity": "source", "input_identity": "input"},
            {"expected_family_id": "fam", "expected_mode": "einsi", "source_identity": "source", "input_identity": "input"},
            {"expected_family_id": "fam", "expected_mode": "linsi", "source_identity": "other-source", "input_identity": "input"},
            {"expected_family_id": "fam", "expected_mode": "linsi", "source_identity": "source", "input_identity": "other-input"},
        ):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                read_family_projection(path, **kwargs)

    def test_family_projection_rejects_malformed_aliases_and_columns(self):
        occurrences, paths, sequences = self.split_fixture()
        index = CodingProjectionIndex(occurrences, sequences, paths)
        with patch(
            "intraphy.coding_correspondence.alignment.protein_multiple_alignment",
            side_effect=unchanged_protein_alignment,
        ):
            projection = index._family_alignment("fam")
        path = self.root / "malformed-family.json"
        write_family_projection(path, projection, source_identity="source", input_identity="input")
        payload = json.loads(path.read_text(encoding="utf-8"))
        record_id = next(iter(payload["records"]))
        original_record = dict(payload["records"][record_id])
        original_record["aliases"] = list(original_record["aliases"])
        original_record["residue_columns"] = list(original_record["residue_columns"])
        payload["records"][record_id]["aliases"] = [
            original_record["aliases"][0], original_record["aliases"][0],
        ]
        path.write_text(json.dumps(payload), encoding="utf-8")
        with self.assertRaises(ValueError):
            read_family_projection(
                path, expected_family_id="fam", expected_mode="linsi",
                source_identity="source", input_identity="input",
            )
        payload["records"][record_id] = {
            **original_record,
            "aliases": list(original_record["aliases"]),
            "residue_columns": [999],
        }
        path.write_text(json.dumps(payload), encoding="utf-8")
        with self.assertRaises(ValueError):
            read_family_projection(
                path, expected_family_id="fam", expected_mode="linsi",
                source_identity="source", input_identity="input",
            )

    def test_copy_pair_index_matches_legacy_for_every_directed_occurrence_pair(self):
        occurrences, paths, sequences = self.split_fixture()
        indexed = CodingProjectionIndex(occurrences, sequences, paths)
        legacy = CodingProjectionIndex(occurrences, sequences, paths)
        with patch(
            "intraphy.coding_correspondence.alignment.protein_multiple_alignment",
            side_effect=unchanged_protein_alignment,
        ):
            projection = indexed._family_alignment("fam")
            expected = {
                (query["occurrence_id"], target["occurrence_id"]):
                    legacy.evidence(query["occurrence_id"], target["occurrence_id"])
                for query in occurrences for target in occurrences
                if query["species"] != target["species"]
            }
            project = indexed._canonical_pair_projection
            with patch.object(indexed, "_canonical_pair_projection", wraps=project) as counted:
                observed = dict(indexed_records(
                    indexed, projection, self.root / "projection-work.sqlite",
                ))
        self.assertEqual(set(observed), set(expected))
        expected_count = sum(
            1 for query_key in indexed.transcripts
            if query_key[:3] == copy_key("A")
            and not indexed.transcripts[query_key].unavailable_reason
            and indexed.transcripts[query_key].protein
            for target_key in indexed.transcripts
            if target_key[:3] == copy_key("B")
            and not indexed.transcripts[target_key].unavailable_reason
            and indexed.transcripts[target_key].protein
        )
        self.assertEqual(counted.call_count, expected_count)
        self.assertEqual(indexed.cache, {})
        self.assertEqual(indexed.cached_bases, 0)
        for key in sorted(expected):
            with self.subTest(pair=key):
                self.assertEqual(observed[key], expected[key])

    def test_installing_changed_family_projection_evicts_stale_pair_cache(self):
        occurrences, paths, sequences = self.split_fixture()
        index = CodingProjectionIndex(occurrences, sequences, paths)
        with patch(
            "intraphy.coding_correspondence.alignment.protein_multiple_alignment",
            side_effect=unchanged_protein_alignment,
        ):
            original = index._family_alignment("fam")
            before = index.evidence("left", "ref")
        self.assertTrue(index.cache)
        self.assertGreater(index.cached_bases, 0)
        shifted = replace(
            original,
            aligned_records={key: "-" + value for key, value in original.aligned_records.items()},
            residue_columns={key: tuple(column + 1 for column in columns)
                             for key, columns in original.residue_columns.items()},
        )
        index.use_family_projection(shifted)
        self.assertEqual(index.cache, {})
        self.assertEqual(index.cached_bases, 0)
        after = index.evidence("left", "ref")
        self.assertEqual(
            after["protein_candidate_set"].candidates[0].coordinate_blocks,
            before["protein_candidate_set"].candidates[0].coordinate_blocks,
        )
        self.assertEqual(after["protein_msa_column_start0"], before["protein_msa_column_start0"] + 1)

    def test_two_copy_worker_accepts_full_family_projection_aliases(self):
        occurrences, paths, sequences = self.split_fixture()
        third, third_path = coding_rows("third", "C", 901, "GCT" * 3)
        occurrences.append(third)
        paths.append(third_path)
        sequences["third"] = "GCT" * 3
        full_family = CodingProjectionIndex(occurrences, sequences, paths)
        worker = CodingProjectionIndex(
            occurrences[:3], {key: sequences[key] for key in ("left", "right", "ref")}, paths[:3],
        )
        legacy = CodingProjectionIndex(
            occurrences[:3], {key: sequences[key] for key in ("left", "right", "ref")}, paths[:3],
        )
        with patch(
            "intraphy.coding_correspondence.alignment.protein_multiple_alignment",
            side_effect=unchanged_protein_alignment,
        ):
            projection = full_family._family_alignment("fam")
            self.assertTrue(any(key[1] == "C" for key in projection.record_by_transcript))
            observed = dict(worker.evidence_for_copy_pair(
                copy_key("A"), copy_key("B"), projection,
                work_db_path=self.root / "full-family-aliases.sqlite",
            ))
            expected = {
                (left["occurrence_id"], right["occurrence_id"]):
                    legacy.evidence(left["occurrence_id"], right["occurrence_id"])
                for left in occurrences[:3] for right in occurrences[:3]
                if left["species"] != right["species"]
            }
        self.assertEqual(observed, expected)

    def test_transcript_aliases_with_distinct_codon_sources_stay_distinct(self):
        query, path_a = coding_rows("q", "A", 101, "GCTGCTGCTGCT", transcript="tx_a")
        path_b = {
            **path_a, "transcript_id": "tx_b", "path_id": "A:tx_b:q",
            "cds_intervals": "104-112", "cds_length": "9",
        }
        path_a = {**path_a, "cds_intervals": "101-109", "cds_length": "9"}
        target, target_path = coding_rows("t", "B", 501, "GCTGCTGCT")
        occurrences, paths = [query, target], [path_a, path_b, target_path]
        sequences = {"q": "GCTGCTGCTGCT", "t": "GCTGCTGCT"}
        indexed = CodingProjectionIndex(occurrences, sequences, paths)
        legacy = CodingProjectionIndex(occurrences, sequences, paths)
        with patch(
            "intraphy.coding_correspondence.alignment.protein_multiple_alignment",
            side_effect=unchanged_protein_alignment,
        ):
            projection = indexed._family_alignment("fam")
            observed = dict(indexed_records(
                indexed, projection, self.root / "alias-work.sqlite",
            ))[("q", "t")]
            expected = legacy.evidence("q", "t")
        self.assertEqual(observed, expected)
        self.assertEqual(observed["protein_mapping_status"], "ambiguous_mapping")
        candidates = observed["protein_candidate_set"].candidates
        self.assertEqual(
            [(item.query_transcript_id, item.target_transcript_id) for item in candidates],
            [("tx_a", "tx"), ("tx_b", "tx")],
        )
        self.assertEqual(observed["protein_best_query_transcript"], "tx_a")
        self.assertEqual(observed["protein_best_target_transcript"], "tx")
        self.assertEqual(observed["protein_supporting_transcripts"], "tx_a>tx;tx_b>tx")
        self.assertFalse(observed["protein_hard_observation_eligible"])

    def test_reverse_asymmetric_projection_matches_legacy_and_swaps_axes(self):
        occurrences, paths, sequences = self.split_fixture(strand="-")
        indexed, legacy = CodingProjectionIndex(occurrences, sequences, paths), CodingProjectionIndex(occurrences, sequences, paths)
        with patch(
            "intraphy.coding_correspondence.alignment.protein_multiple_alignment",
            side_effect=unchanged_protein_alignment,
        ):
            projection = indexed._family_alignment("fam")
            forward = dict(indexed.evidence_for_copy_pair(
                copy_key("A"), copy_key("B"), projection,
                work_db_path=self.root / "reverse-forward.sqlite",
            ))
            reverse = dict(indexed.evidence_for_copy_pair(
                copy_key("B"), copy_key("A"), projection,
                work_db_path=self.root / "reverse-backward.sqlite",
            ))
            for query in occurrences:
                for target in occurrences:
                    if query["species"] == target["species"]:
                        continue
                    key = (query["occurrence_id"], target["occurrence_id"])
                    self.assertEqual(forward[key], legacy.evidence(*key))
                    self.assertEqual(reverse[key], legacy.evidence(*key))
        blocks = forward[("left", "ref")]["protein_projected_coordinate_blocks"]
        reverse_blocks = forward[("ref", "left")]["protein_projected_coordinate_blocks"]
        self.assertTrue(blocks)
        reverse_coordinates = tuple(sorted(
            (block.query.start0, block.query.end0, block.target.start0, block.target.end0)
            for block in reverse_blocks
        ))
        expected_reversed = tuple(sorted(
            (block.target.start0, block.target.end0, block.query.start0, block.query.end0)
            for block in blocks
        ))
        self.assertEqual(reverse_coordinates, expected_reversed)

    def test_competing_projection_checks_all_alternative_transcripts(self):
        occurrences, paths, sequences = [], [], {}
        for index in range(4):
            occurrence, path = coding_rows(
                f"q{index}", "A", 101 + index * 100, "GCT" * 12,
                transcript=f"tx{index}",
            )
            occurrences.append(occurrence)
            paths.append(path)
            sequences[occurrence["occurrence_id"]] = "GCT" * 12
        target, target_path = coding_rows("t", "B", 901, "GCT" * 12)
        occurrences.append(target)
        paths.append(target_path)
        sequences["t"] = "GCT" * 12
        indexed, legacy = CodingProjectionIndex(occurrences, sequences, paths), CodingProjectionIndex(occurrences, sequences, paths)
        with patch(
            "intraphy.coding_correspondence.alignment.protein_multiple_alignment",
            side_effect=unchanged_protein_alignment,
        ):
            projection = indexed._family_alignment("fam")
            observed = dict(indexed_records(
                indexed, projection, self.root / "competitor-work.sqlite",
            ))[("q0", "t")]
            expected = legacy.evidence("q0", "t")
        self.assertEqual(observed, expected)
        self.assertEqual(observed["protein_competing_occurrences"], "query:q1;query:q2;query:q3")
        self.assertEqual(observed["protein_mapping_status"], "ambiguous_mapping")
        self.assertFalse(observed["protein_position_eligible"])

    def test_unavailable_and_unaligned_cds_states_match_legacy(self):
        query, query_path = coding_rows("q", "A", 101, "GCT" * 4)
        target, target_path = coding_rows("t", "B", 501, "GCT" * 4)
        no_cds_path = {**query_path, "cds_intervals": ".", "cds_length": "0", "cds_phase": "."}
        no_cds = CodingProjectionIndex(
            [query, target], {"q": "GCT" * 4, "t": "GCT" * 4}, [no_cds_path, target_path],
        )
        legacy = CodingProjectionIndex(
            [query, target], {"q": "GCT" * 4, "t": "GCT" * 4}, [no_cds_path, target_path],
        )
        with patch(
            "intraphy.coding_correspondence.alignment.protein_multiple_alignment",
            side_effect=unchanged_protein_alignment,
        ):
            projection = no_cds._family_alignment("fam")
            observed = dict(no_cds.evidence_for_copy_pair(
                copy_key("A"), copy_key("B"), projection,
                work_db_path=self.root / "no-cds-work.sqlite",
            ))[("q", "t")]
            expected = legacy.evidence("q", "t")
        self.assertEqual(observed, expected)
        self.assertIn("no_annotated_CDS", observed["protein_unavailable_reason"])

        query_l, query_path_l = coding_rows("ql", "A", 101, "GCT" * 4)
        target_l, target_path_l = coding_rows("tl", "B", 501, "TTA" * 4)
        aligned = CodingProjectionIndex(
            [query_l, target_l], {"ql": "GCT" * 4, "tl": "TTA" * 4}, [query_path_l, target_path_l],
        )
        legacy_aligned = CodingProjectionIndex(
            [query_l, target_l], {"ql": "GCT" * 4, "tl": "TTA" * 4}, [query_path_l, target_path_l],
        )
        with patch(
            "intraphy.coding_correspondence.alignment.protein_multiple_alignment",
            side_effect=disjoint_protein_alignment,
        ):
            projection = aligned._family_alignment("fam")
            observed = dict(aligned.evidence_for_copy_pair(
                copy_key("A"), copy_key("B"), projection,
                work_db_path=self.root / "unaligned-work.sqlite",
            ))[("ql", "tl")]
            expected = legacy_aligned.evidence("ql", "tl")
        self.assertEqual(observed, expected)
        self.assertEqual(observed["protein_status"], "no_aligned_CDS")
        self.assertFalse(observed["protein_candidate_evidence_available"])
