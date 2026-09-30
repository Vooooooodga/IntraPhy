"""preparation / transcripts: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

from intraphy.aligners.scoring import splice_motif_score
from intraphy.partial_boundaries import generic_partial, parse_attributes, range_declared, range_value
from intraphy.preparation.annotation_index import overlaps_gene
from intraphy.storage.tabular import read_tsv
from pathlib import Path


def sequence_slice(seqs, contig, start, end, strand):
    seq = seqs.get(contig, "")
    if not seq:
        return ""
    bounds = (seqs.get("__bounds__") or {}).get(contig)
    if bounds:
        lower, upper = bounds
        if start < lower or end > upper:
            raise SystemExit(f"requested interval {contig}:{start}-{end} outside loaded FASTA window {lower}-{upper}")
        local_start = start - lower + 1
        local_end = end - lower + 1
    else:
        local_start = start
        local_end = end
    sub = seq[local_start - 1 : local_end]
    if strand == "-":
        table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
        sub = sub.translate(table)[::-1]
    return sub.upper()


def translate_cds(sequence):
    sequence = (sequence or "").upper().replace("U", "T")
    try:
        from Bio.Seq import Seq
    except ImportError as exc:
        raise SystemExit("Biopython is required to translate reconstructed CDS sequences") from exc
    usable = sequence[: len(sequence) - (len(sequence) % 3)]
    return str(Seq(usable).translate()) if usable else ""


def transcript_features(features, gene, gene_ids):
    transcript_types = {"mrna", "transcript", "lnc_rna", "ncrna", "rrna", "trna"}
    transcripts = []
    for row in features:
        if row["type"].lower() not in transcript_types:
            continue
        if not overlaps_gene(row, gene):
            continue
        if set(row.get("parents", [])) & gene_ids or row.get("id") in gene_ids:
            transcripts.append(row)
    if not transcripts:
        transcripts = [
            {
                "seqid": gene["seqid"],
                "source": gene.get("source", "IntraPhy"),
                "type": "synthetic_transcript",
                "start": gene["start"],
                "end": gene["end"],
                "score": ".",
                "strand": gene["strand"],
                "phase": ".",
                "attrs": {"ID": f"{gene.get('id', 'gene')}.synthetic_tx"},
                "id": f"{gene.get('id', 'gene')}.synthetic_tx",
                "parent": gene.get("id", ""),
                "parents": [gene.get("id", "")],
                "name": f"{gene.get('id', 'gene')}.synthetic_tx",
            }
        ]
    return sorted(transcripts, key=lambda row: (row["start"], row["end"], row.get("id", "")))


def feature_role(feature):
    ftype = feature["type"].lower()
    if ftype == "cds":
        return "CDS"
    if ftype in {"utr", "five_prime_utr", "three_prime_utr"}:
        return "UTR"
    if ftype == "intron":
        return "intron"
    if ftype == "gene_body":
        return "exon"
    return "exon"


def _merge_exonic_intervals(features, transcript):
    """Reconstruct exon intervals when an annotation omits explicit exon rows."""
    intervals = sorted(features, key=lambda row: (row["start"], row["end"]))
    merged = []
    for feature in intervals:
        if not merged or feature["start"] > merged[-1]["end"] + 1:
            merged.append(
                {
                    **feature,
                    "type": "exon",
                    "id": f"{transcript.get('id', 'tx')}.reconstructed_exon_{len(merged) + 1}",
                    "attrs": {"ID": f"{transcript.get('id', 'tx')}.reconstructed_exon_{len(merged) + 1}"},
                    "components": [feature],
                    "reconstructed": True,
                }
            )
        else:
            merged[-1]["end"] = max(merged[-1]["end"], feature["end"])
            merged[-1]["components"].append(feature)
    return merged


def _annotate_exon(exon, child_rows):
    item = dict(exon)
    overlapping = [
        row for row in child_rows
        if row["start"] <= item["end"] and row["end"] >= item["start"]
    ]
    cds = sorted(
        [row for row in overlapping if row["type"].lower() == "cds"],
        key=lambda row: transcript_sort_key(row, item.get("strand", "+")),
    )
    utr = [row for row in overlapping if row["type"].lower() in {"utr", "five_prime_utr", "three_prime_utr"}]
    item["cds_intervals"] = [(row["start"], row["end"]) for row in cds]
    item["utr_intervals"] = [(row["start"], row["end"]) for row in utr]
    item["utr_types"] = [row["type"] for row in utr]
    item["cds_length"] = sum(row["end"] - row["start"] + 1 for row in cds)
    item["cds_start"] = min((row["start"] for row in cds), default=None)
    item["cds_end"] = max((row["end"] for row in cds), default=None)
    item["cds_phase"] = cds[0].get("phase", ".") if cds else "."
    item["coding_status"] = "coding" if cds else "noncoding"
    item["utr_status"] = "contains_utr" if utr else "no_annotated_utr"
    item["phase"] = item["cds_phase"]
    return item


def _format_intervals(intervals):
    if not intervals:
        return "NA"
    return ";".join(f"{int(start)}-{int(end)}" for start, end in sorted(intervals))


def _feature_parents(feature):
    parents = set(feature.get("parents", []))
    parent = feature.get("parent")
    if parent:
        parents.add(parent)
    return sorted(parent for parent in parents if parent)


def _format_attrs(attrs):
    return ";".join(f"{key}={attrs[key]}" for key in sorted(attrs)) if attrs else "NA"


def _partial_boundary(feature, boundary, *, localized_range=False):
    attrs = parse_attributes(feature.get("attrs", {}) or {})
    keys = (
        ("partial_start", "start_partial", "partial5")
        if boundary == "start"
        else ("partial_end", "end_partial", "partial3")
    )
    for key in keys:
        if key in attrs:
            return attrs[key]
    strand = feature.get("strand", "+")
    genomic_boundary = ("start" if (boundary == "start") == (strand != "-") else "end")
    if range_value(attrs, genomic_boundary):
        return "1"
    if range_declared(attrs, genomic_boundary):
        return "1"
    if generic_partial(attrs) and not localized_range:
        return "1"
    return "0"


def _coding_role(feature):
    if feature_role(feature) == "intron":
        return "unknown"
    cds = feature.get("cds_intervals", ())
    utr_types = {str(value).lower() for value in feature.get("utr_types", ())}
    if cds and utr_types:
        return "mixed"
    if cds:
        return "CDS"
    if utr_types == {"five_prime_utr"}:
        return "five_prime_UTR"
    if utr_types == {"three_prime_utr"}:
        return "three_prime_UTR"
    if utr_types:
        return "mixed" if len(utr_types) > 1 else "unknown"
    return "noncoding_exon"


def _position_role(path_features, index):
    feature = path_features[index]
    if feature_role(feature) == "intron":
        return "internal"
    exonic_indices = [
        item_index
        for item_index, item in enumerate(path_features)
        if feature_role(item) != "intron"
    ]
    if len(exonic_indices) == 1:
        return "single"
    if index == exonic_indices[0]:
        return "first"
    if index == exonic_indices[-1]:
        return "last"
    return "internal"


def _path_role_record(transcript_id, rank, path_features, index, feature, transcript=None):
    position_role = _position_role(path_features, index)
    path_has_range = any(
        range_value(item.get("attrs", {}), boundary)
        for item in (*path_features, *((transcript,) if transcript else ()))
        for boundary in ("start", "end")
    )
    localized_range = path_has_range and generic_partial(feature.get("attrs", {}))
    partial_start = _partial_boundary(feature, "start", localized_range=localized_range)
    partial_end = _partial_boundary(feature, "end", localized_range=localized_range)
    if transcript is not None and position_role in {"first", "single"} and partial_start == "0":
        parent_localized = path_has_range and generic_partial(transcript.get("attrs", {}))
        partial_start = _partial_boundary(transcript, "start", localized_range=parent_localized)
    if transcript is not None and position_role in {"last", "single"} and partial_end == "0":
        parent_localized = path_has_range and generic_partial(transcript.get("attrs", {}))
        partial_end = _partial_boundary(transcript, "end", localized_range=parent_localized)
    # A complete alternative transcript may start/end within the gene envelope.
    # Preserve explicit partial flags; an unmarked end means unassessed completeness,
    # not biological truncation inferred from another transcript's coordinates.
    return {
        "transcript_id": transcript_id,
        "path_rank": rank,
        "source_feature_id": feature.get("id", "NA") or "NA",
        "source_feature_type": feature.get("type", "NA") or "NA",
        "path_role": "intronic" if feature_role(feature) == "intron" else "exonic",
        "coding_role": _coding_role(feature),
        "position_role": position_role,
        "annotation_source": feature.get("source", (transcript or {}).get("source", "NA")) or "NA",
        "original_attributes": _format_attrs(feature.get("attrs", {})),
        "partial_start": partial_start,
        "partial_end": partial_end,
    }


def _interval_union_length(intervals):
    cleaned = sorted((int(start), int(end)) for start, end in intervals if int(end) >= int(start))
    if not cleaned:
        return 0
    merged = []
    for start, end in cleaned:
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return sum(end - start + 1 for start, end in merged)


def child_features_for_transcript(features, gene, gene_ids, transcript):
    tx_id = transcript.get("id", "")
    selected = []
    allowed = {"exon", "cds", "utr", "five_prime_utr", "three_prime_utr"}
    synthetic = transcript["type"] == "synthetic_transcript"
    for row in features:
        if row["type"].lower() not in allowed or not overlaps_gene(row, gene):
            continue
        parents = set(row.get("parents", []))
        if (tx_id and tx_id in parents) or (synthetic and (parents & gene_ids or row.get("id") in gene_ids)):
            selected.append(row)
    if not selected:
        return []
    explicit_exons = [row for row in selected if row["type"].lower() == "exon"]
    components = [row for row in selected if row["type"].lower() != "exon"]
    exons = explicit_exons or _merge_exonic_intervals(components, transcript)
    return sorted(
        [_annotate_exon(exon, components) for exon in exons],
        key=lambda row: transcript_sort_key(row, gene["strand"]),
    )


def transcript_sort_key(row, strand):
    key = (int(row["start"]), int(row["end"]))
    if strand == "-":
        key = (-int(row["end"]), -int(row["start"]))
    return key


def transcript_cds_length(tx_features):
    intervals = []
    for row in tx_features:
        intervals.extend(row.get("cds_intervals") or [])
        if row.get("type", "").lower() == "cds":
            intervals.append((row["start"], row["end"]))
    return _interval_union_length(intervals)


def select_transcripts(transcripts, features_by_tx, transcript_policy="canonical", canonical_rule="longest_cds"):
    if transcript_policy == "all":
        return transcripts
    if not transcripts:
        return []
    if canonical_rule == "longest_span":
        return [max(transcripts, key=lambda tx: (tx["end"] - tx["start"] + 1, tx.get("id", "")))]
    return [max(transcripts, key=lambda tx: (transcript_cds_length(features_by_tx[tx["id"]]), tx["end"] - tx["start"] + 1, tx.get("id", "")))]


def introns_from_path(path_features, gene, transcript_id, seqs=None):
    exonic = [row for row in path_features if feature_role(row) in {"CDS", "UTR", "exon"}]
    merged = []
    for exon in sorted(exonic, key=lambda row: (row["start"], row["end"])):
        if not merged or exon["start"] > merged[-1]["end"] + 1:
            merged.append({"start": exon["start"], "end": exon["end"], "left": exon, "right": exon})
        else:
            merged[-1]["end"] = max(merged[-1]["end"], exon["end"])
            merged[-1]["right"] = exon
    introns = []
    for idx, left in enumerate(merged[:-1], start=1):
        right = merged[idx]
        start = left["end"] + 1
        end = right["start"] - 1
        if start <= end:
            intron_seq = sequence_slice(seqs or {}, gene["seqid"], start, end, gene["strand"])
            motif_score, donor, acceptor = splice_motif_score(intron_seq)
            upstream, downstream = left["right"], right["left"]
            if gene["strand"] == "-":
                upstream, downstream = downstream, upstream
            introns.append(
                {
                    "seqid": gene["seqid"],
                    "type": "intron",
                    "start": start,
                    "end": end,
                    "strand": gene["strand"],
                    "phase": ".",
                    "id": f"{transcript_id}.intron_{idx}",
                    "parent": transcript_id,
                    "parents": [transcript_id],
                    "attrs": {"ID": f"{transcript_id}.intron_{idx}"},
                    "left_feature_id": upstream.get("id", "NA") or "NA",
                    "right_feature_id": downstream.get("id", "NA") or "NA",
                    "left_phase": upstream.get("phase", "."),
                    "right_phase": downstream.get("phase", "."),
                    "left_cds_length": upstream.get("cds_length", 0),
                    "splice_motif_score": motif_score,
                    "splice_donor": donor,
                    "splice_acceptor": acceptor,
                }
            )
    return introns


def exon_features_for_gene(features, gene, gene_ids):
    transcripts = transcript_features(features, gene, gene_ids)
    features_by_tx = {tx["id"]: child_features_for_transcript(features, gene, gene_ids, tx) for tx in transcripts}
    selected = select_transcripts(transcripts, features_by_tx)
    rows = []
    for tx in selected:
        rows.extend(features_by_tx[tx["id"]])
    return sorted(rows, key=lambda row: (row["start"], row["end"], feature_role(row)))


def load_existing_segments(output_dir):
    path = Path(output_dir) / "segment_occurrences.tsv"
    if not path.exists():
        return []
    return read_tsv(path, optional=True)


def _feature_key(feature):
    return (
        feature["seqid"],
        int(feature["start"]),
        int(feature["end"]),
        feature.get("strand", "."),
        feature_role(feature),
        feature.get("phase", "."),
    )


def _occurrence_id(species, gene_copy_id, index, role):
    safe_role = role.replace("/", "_")
    return f"{species}_{gene_copy_id}_{index:03d}_{safe_role}"
