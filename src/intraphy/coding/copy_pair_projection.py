"""Public copy-pair projection entry point over a prepared coding index."""
from __future__ import annotations

from collections import OrderedDict
from dataclasses import replace
import json
from pathlib import Path
import sqlite3

from intraphy.coding.copy_pair_codec import decode_candidate, encode_candidate
from intraphy.coding.evidence import finalize_evidence
from intraphy.coding.projection import _anchor_metrics, _coordinate_blocks0
from intraphy.coding.types import CodingProjectionCandidate
from intraphy.mapping.fields import EXON_LIKE_ROLES


def evidence_for_copy_pair(index, query_copy_key, target_copy_key, family_projection=None,
                           *, work_db_path):
    """Build exact directed evidence for one pair using index-owned MSA state."""
    query_copy_key, target_copy_key = tuple(query_copy_key), tuple(target_copy_key)
    if len(query_copy_key) != 3 or len(target_copy_key) != 3 or query_copy_key == target_copy_key:
        raise ValueError("two distinct three-part copy keys are required")
    if query_copy_key[0] != target_copy_key[0]:
        raise ValueError("copy-pair evidence requires copies from one family")
    if family_projection is not None:
        index.use_family_projection(family_projection)
    query_keys = tuple(sorted(key for key in index.transcripts if key[:3] == query_copy_key))
    target_keys = tuple(sorted(key for key in index.transcripts if key[:3] == target_copy_key))
    query_occurrences = tuple(sorted(occ for occ in index.occurrences_by_copy.get(query_copy_key, ())
        if index.occurrences[occ].get("role") in EXON_LIKE_ROLES))
    target_occurrences = tuple(sorted(occ for occ in index.occurrences_by_copy.get(target_copy_key, ())
        if index.occurrences[occ].get("role") in EXON_LIKE_ROLES))
    if not query_occurrences or not target_occurrences:
        raise ValueError("copy pair has no cross-copy exon-like occurrence pairs")
    work_db_path = Path(work_db_path)
    work_db_path.parent.mkdir(parents=True, exist_ok=True)
    if work_db_path.exists():
        raise FileExistsError(work_db_path)
    db = None
    try:
        db = sqlite3.connect(work_db_path)
        db.execute("PRAGMA journal_mode=DELETE")
        db.execute("PRAGMA synchronous=OFF")
        db.execute("CREATE TABLE query_bits (transcript_key TEXT NOT NULL, focal TEXT NOT NULL, partner TEXT NOT NULL, bits BLOB NOT NULL, PRIMARY KEY(transcript_key,focal,partner)) WITHOUT ROWID")
        db.execute("CREATE TABLE target_bits (transcript_key TEXT NOT NULL, focal TEXT NOT NULL, partner TEXT NOT NULL, bits BLOB NOT NULL, PRIMARY KEY(transcript_key,focal,partner)) WITHOUT ROWID")
        db.execute("CREATE TABLE candidate (query_occ TEXT NOT NULL, target_occ TEXT NOT NULL, query_key TEXT NOT NULL, target_key TEXT NOT NULL, payload TEXT NOT NULL, query_bits BLOB NOT NULL, target_bits BLOB NOT NULL, PRIMARY KEY(query_occ,target_occ,query_key,target_key)) WITHOUT ROWID")
    except Exception:
        if db is not None:
            db.close()
        work_db_path.unlink(missing_ok=True)
        raise

    def add_mask(table, transcript_key, focal, partner, mask):
        key = json.dumps(transcript_key, separators=(",", ":"))
        row = db.execute(f"SELECT bits FROM {table} WHERE transcript_key=? AND focal=? AND partner=?",
                         (key, focal, partner)).fetchone()
        if row:
            mask |= int.from_bytes(row[0], "little")
        blob = mask.to_bytes(max(1, (mask.bit_length() + 7) // 8), "little")
        db.execute(f"INSERT OR REPLACE INTO {table} VALUES (?,?,?,?)", (key, focal, partner, blob))

    def candidate_for(qkey, tkey, qocc, tocc, record, columns, aq, at):
        positions = set(record["positions0"])
        anchors = _anchor_metrics(record, columns, aq, at, qocc, tocc)
        ordered = sorted(positions)
        qpos = {q0 for q0, _t0 in positions}
        tpos = {t0 for _q0, t0 in positions}
        candidate = CodingProjectionCandidate(
            query_transcript_key=qkey, target_transcript_key=tkey,
            coordinate_blocks=_coordinate_blocks0(positions),
            aa_identity=record["aa_matches"] / record["aa_pairs"],
            query_cds_coverage=len(positions) / index.transcripts[qkey].coding_lengths[qocc],
            target_cds_coverage=len(positions) / index.transcripts[tkey].coding_lengths[tocc],
            known_aa_pairs=record["aa_pairs"], blosum62_score=record["blosum62_score"],
            gap_fraction=anchors["gap_fraction"], left_anchor_pairs=anchors["left_pairs"],
            right_anchor_pairs=anchors["right_pairs"], left_anchor_score=anchors["left_score"],
            right_anchor_score=anchors["right_score"], left_anchor_supported=anchors["left_supported"],
            right_anchor_supported=anchors["right_supported"], terminal_side=anchors["terminal_side"],
            msa_column_interval=anchors["column_interval"], anchor_resolved=anchors["resolved"],
            position_monotonic=all(right[0] > left[0] and right[1] > left[1]
                                    for left, right in zip(ordered, ordered[1:])),
            competing_occurrences=(),
        )
        return candidate, sum(1 << pos for pos in qpos), sum(1 << pos for pos in tpos)

    try:
        qocc_set, tocc_set = set(query_occurrences), set(target_occurrences)
        for qkey in query_keys:
            qtx = index.transcripts[qkey]
            if qtx.unavailable_reason or not qtx.protein:
                continue
            for tkey in target_keys:
                ttx = index.transcripts[tkey]
                if ttx.unavailable_reason or not ttx.protein:
                    continue
                canonical_key, canonical_pairs, known_columns, aligned_left, aligned_right = index._canonical_pair_projection(
                    qkey, tkey, cache=False)
                if qkey == canonical_key[0]:
                    pairs, columns, aligned_q, aligned_t = canonical_pairs, known_columns, aligned_left, aligned_right
                else:
                    pairs = {(right, left): {**value,
                              "positions0": [(b, a) for a, b in value["positions0"]]}
                             for (left, right), value in canonical_pairs.items()}
                    columns = tuple((c, ta, qa, tc, qc) for c, qa, ta, qc, tc in known_columns)
                    aligned_q, aligned_t = aligned_right, aligned_left
                reversed_columns = tuple((column, target_aa, query_aa,
                                          target_codon, query_codon)
                                         for column, query_aa, target_aa,
                                             query_codon, target_codon in columns)
                for (qocc, tocc), record in pairs.items():
                    if qocc not in qocc_set or tocc not in tocc_set:
                        continue
                    if not qtx.coding_lengths.get(qocc) or not ttx.coding_lengths.get(tocc):
                        continue
                    candidate, qmask, tmask = candidate_for(
                        qkey, tkey, qocc, tocc, record, columns, aligned_q, aligned_t)
                    add_mask("query_bits", qkey, qocc, tocc, qmask)
                    add_mask("target_bits", tkey, tocc, qocc, tmask)
                    db.execute("INSERT INTO candidate VALUES (?,?,?,?,?,?,?)",
                        (qocc, tocc, json.dumps(qkey, separators=(",", ":")),
                         json.dumps(tkey, separators=(",", ":")),
                         json.dumps(encode_candidate(candidate), sort_keys=True, separators=(",", ":")),
                         qmask.to_bytes(max(1, (qmask.bit_length() + 7) // 8), "little"),
                         tmask.to_bytes(max(1, (tmask.bit_length() + 7) // 8), "little")))
                    reversed_record = {**record, "positions0": [(b, a) for a, b in record["positions0"]]}
                    reverse_candidate, reverse_qmask, reverse_tmask = candidate_for(
                        tkey, qkey, tocc, qocc, reversed_record, reversed_columns, aligned_t, aligned_q)
                    add_mask("query_bits", tkey, tocc, qocc, reverse_qmask)
                    add_mask("target_bits", qkey, qocc, tocc, reverse_tmask)
                    db.execute("INSERT INTO candidate VALUES (?,?,?,?,?,?,?)",
                        (tocc, qocc, json.dumps(tkey, separators=(",", ":")),
                         json.dumps(qkey, separators=(",", ":")),
                         json.dumps(encode_candidate(reverse_candidate), sort_keys=True, separators=(",", ":")),
                         reverse_qmask.to_bytes(max(1, (reverse_qmask.bit_length() + 7) // 8), "little"),
                         reverse_tmask.to_bytes(max(1, (reverse_tmask.bit_length() + 7) // 8), "little")))
                db.commit()
                del canonical_pairs, known_columns, aligned_left, aligned_right
                del pairs, columns, aligned_q, aligned_t, reversed_columns
    except Exception:
        db.close()
        work_db_path.unlink(missing_ok=True)
        raise

    def records():
        query_cache, target_cache = OrderedDict(), OrderedDict()

        def competing(table, transcript_key, focal, partner, mask, label, cache):
            group = (json.dumps(transcript_key, separators=(",", ":")), focal)
            if group not in cache:
                cache[group] = {other: int.from_bytes(bits, "little")
                    for other, bits in db.execute(
                        f"SELECT partner,bits FROM {table} WHERE transcript_key=? AND focal=?", group)}
                if len(cache) > 8:
                    cache.popitem(last=False)
            cache.move_to_end(group)
            return [f"{label}:{other}" for other, bits in cache[group].items()
                    if other != partner and bits & mask]

        try:
            for qocc in query_occurrences:
                for tocc in target_occurrences:
                    for q, t in ((qocc, tocc), (tocc, qocc)):
                        qkeys = sorted(index.by_occurrence.get(q, ()))
                        tkeys = sorted(index.by_occurrence.get(t, ()))
                        unavailable, protein_status, candidates = set(), "unavailable", []
                        for qkey in qkeys:
                            query = index.transcripts[qkey]
                            for tkey in tkeys:
                                target = index.transcripts[tkey]
                                for side, transcript in (("query", query), ("target", target)):
                                    if transcript.unavailable_reason:
                                        unavailable.add(f"{side}:{transcript.key[-1]}:{transcript.unavailable_reason}")
                                if query.unavailable_reason or target.unavailable_reason:
                                    continue
                                if not query.coding_lengths.get(q) or not target.coding_lengths.get(t):
                                    unavailable.add("no_CDS_in_requested_occurrence")
                                    continue
                                protein_status = "no_aligned_CDS"
                        for qkey, tkey, payload, qblob, tblob in db.execute(
                                "SELECT query_key,target_key,payload,query_bits,target_bits FROM candidate WHERE query_occ=? AND target_occ=? ORDER BY query_key,target_key", (q, t)):
                            qkey, tkey = tuple(json.loads(qkey)), tuple(json.loads(tkey))
                            candidate = decode_candidate(json.loads(payload))
                            qmask, tmask = int.from_bytes(qblob, "little"), int.from_bytes(tblob, "little")
                            competitors = set(competing("query_bits", qkey, q, t, qmask, "target", query_cache))
                            competitors.update(competing("target_bits", tkey, t, q, tmask, "query", target_cache))
                            candidates.append(replace(candidate, competing_occurrences=tuple(sorted(competitors))))
                        base = {"protein_status": protein_status, "protein_mapping_status": "uncovered",
                            "protein_unavailable_reason": "no_CDS_transcript_path",
                            "protein_membership_eligible": False, "protein_position_eligible": False,
                            "protein_hard_observation_eligible": False,
                            "protein_candidate_evidence_available": False, "_unavailable": unavailable}
                        yield (q, t), finalize_evidence(base, candidates, index.transcripts, index.msa_mode)
        finally:
            db.close()
            work_db_path.unlink(missing_ok=True)

    yield from records()
