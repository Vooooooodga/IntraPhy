"""Prepared genomic sources to complete exon instances, without tree inference."""
from __future__ import annotations
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from ..storage.fasta import parse_fasta
from ..storage.tabular import read_tsv
from ..partial_boundaries import generic_partial, parse_attributes, range_declared, range_value
from .types import ExonInstance


@dataclass(frozen=True)
class NativeLocus:
    family: str
    species: str
    locus: str
    contig: str
    strand: str
    search_start: int
    search_end: int
    sequence: str
    exons: tuple[ExonInstance, ...]
    paths: dict[str, tuple[str, ...]]
    partial: bool
    coding_by_transcript: dict[str, tuple[tuple[int, int, int | None], ...]] = field(default_factory=dict)
    translation_exceptions: tuple[str, ...] = ()
    partial_intervals: tuple[tuple[int, int], ...] = ()
    partial_location_known: bool = False

    def oriented_interval(self, exon: ExonInstance) -> tuple[int, int]:
        if self.strand == "+":
            return exon.start-(self.search_start-1), exon.end-(self.search_start-1)
        return self.search_end-exon.end, self.search_end-exon.start

    def source_position(self, oriented_base: int) -> int:
        """Return 0-based source-genomic base position."""
        return self.search_start-1+oriented_base if self.strand == "+" else self.search_end-1-oriented_base


def _cds_intervals(row):
    value = row.get("cds_intervals", "")
    phase = row.get("cds_phase", ".")
    phase = int(phase) if str(phase) in {"0", "1", "2"} else None
    result = []
    for part in value.split(";"):
        if "-" in part:
            left, right = part.split("-")
            result.append((int(left)-1, int(right), phase))
    return tuple(result)


def _attrs(row):
    return parse_attributes(row.get("original_attributes", row.get("attrs", "")))


def _has_partial_metadata(rows):
    return (any(generic_partial(_attrs(row)) for row in rows) or
            any(row.get("partial_start") == "1" or row.get("partial_end") == "1"
                for row in rows) or
            any(range_declared(_attrs(row), boundary) for row in rows
                for boundary in ("start", "end")))


def _partial_intervals(transcript_rows, search_start, search_end):
    """Return 0-based genomic uncertainty windows for selected transcript paths."""
    intervals = set()
    localized = False
    unlocated_generic = False
    selected_rows = [row for group in transcript_rows.values() for row in group]
    for tx_rows in transcript_rows.values():
        ordered = sorted(tx_rows, key=lambda row: (int(row["start"]), int(row["end"])))
        has_range = any(range_value(_attrs(row), boundary)
                        for row in ordered for boundary in ("start", "end"))
        has_unlocated_range = any(range_declared(_attrs(row), boundary) and
                                  not range_value(_attrs(row), boundary)
                                  for row in ordered for boundary in ("start", "end"))
        has_partial = _has_partial_metadata(ordered)
        if has_partial and (not has_range or has_unlocated_range):
            unlocated_generic = True
        for index, row in enumerate(ordered):
            attrs = _attrs(row)
            start_range = range_value(attrs, "start")
            end_range = range_value(attrs, "end")
            for name, value in (("start", start_range), ("end", end_range)):
                if not value:
                    continue
                localized = True
                bounds = [part.strip() for part in value.split(",")]
                feature_start, feature_end = int(row["start"]), int(row["end"])
                if name == "start":
                    left = (max((int(other["end"]) for other in ordered[:index]
                                 if int(other["end"]) < feature_start), default=search_start - 1)
                            if bounds[0] == "." else int(bounds[0]) - 1)
                    right = feature_start - 1 if bounds[1] == "." else int(bounds[1]) - 1
                else:
                    left = feature_end if bounds[0] == "." else int(bounds[0])
                    right = (min((int(other["start"]) - 1 for other in ordered[index + 1:]
                                  if int(other["start"]) > feature_end), default=search_end)
                             if bounds[1] == "." else int(bounds[1]))
                left, right = max(search_start - 1, left), min(search_end, right)
                if left < right:
                    intervals.add((left, right))
    has_partial = _has_partial_metadata(selected_rows)
    return tuple(sorted(intervals)), localized, has_partial, unlocated_generic


def load_native(input_dir: str | Path) -> tuple[NativeLocus, ...]:
    root = Path(input_dir)
    loci = read_tsv(root / "gene_loci.tsv")
    paths = read_tsv(root / "transcript_paths.tsv")
    sequences = parse_fasta(root / "gene_loci.fasta")
    raw = read_tsv(root / "raw_gene_features.tsv", optional=True)
    metadata = defaultdict(list)
    for r in raw:
        metadata[(r["species"], r["gene_copy_id"])].append(r)
    exon_rows = [r for r in paths if r.get("role") == "exon"]
    by_locus = defaultdict(list)
    for row in exon_rows:
        by_locus[(row["species"], row["gene_copy_id"])].append(row)
    result = []
    for locus in loci:
        rows = by_locus.get((locus["species"], locus["gene_copy_id"]), [])
        raw_locus = metadata[(locus["species"], locus["gene_copy_id"])]
        by_tx_rows = defaultdict(list)
        for row in rows:
            by_tx_rows[row["transcript_id"]].append(row)
        raw_tx = {row.get("id", ""): row for row in raw_locus
                  if row.get("id") and row.get("type", "").lower() in
                  {"mrna", "transcript", "lnc_rna", "ncrna", "rrna", "trna"} and
                  row.get("ownership") == "target_gene_descendant"}
        for tx, tx_row in raw_tx.items():
            if tx in by_tx_rows:
                by_tx_rows[tx].append(tx_row)
        for raw_exon in raw_locus:
            if (raw_exon.get("type", "").lower() != "exon" or
                    raw_exon.get("ownership") != "target_gene_descendant"):
                continue
            parents = set(raw_exon.get("parents", "").split(";"))
            parents.update(raw_exon.get("parent", "").replace(",", ";").split(";"))
            for tx in parents & by_tx_rows.keys():
                if any(int(path_row["start"]) == int(raw_exon["start"]) and
                       int(path_row["end"]) == int(raw_exon["end"])
                       for path_row in by_tx_rows[tx]):
                    by_tx_rows[tx].append(raw_exon)
        families = {r["family_id"] for r in (rows or raw_locus)}
        if len(families) != 1:
            raise ValueError("Prepared locus belongs to multiple families")
        prefix = f"{locus['species']}|{locus['gene_copy_id']}|"
        keys = [key for key in sequences if key.startswith(prefix)]
        if len(keys) != 1:
            raise ValueError(f"Expected one prepared genomic sequence for {prefix}")
        by_span = defaultdict(list)
        for row in rows:
            by_span[(int(row["start"])-1, int(row["end"]))].append(row)
        instances = []
        by_tx = defaultdict(list)
        for (start, end), evidence in sorted(by_span.items()):
            occurrence = sorted({r["occurrence_id"] for r in evidence})[0]
            transcripts = tuple(sorted({r["transcript_id"] for r in evidence}))
            cds = tuple(sorted({p for r in evidence for p in _cds_intervals(r)}, key=lambda p: (p[0], p[1], str(p[2]))))
            instance = ExonInstance(occurrence, next(iter(families)), locus["species"], locus["gene_copy_id"],
                locus["contig"], start, end, locus["strand"], transcripts, cds,
                ";".join(sorted({r.get("annotation_source", "provided") for r in evidence})))
            instances.append(instance)
            for tx in transcripts:
                by_tx[tx].append(instance)
        reverse = locus["strand"] == "-"
        ordered_paths = {tx: tuple(e.id for e in sorted(es, key=lambda e: e.start, reverse=reverse)) for tx, es in by_tx.items()}
        coding = {}
        for tx in ordered_paths:
            coding[tx] = tuple(sorted({p for r in rows if r["transcript_id"] == tx
                                      for p in _cds_intervals(r)}, key=lambda p: (p[0], p[1], str(p[2])), reverse=reverse))
        exceptions = set()
        for r in raw_locus:
            attrs = r.get("attrs", "")
            if any(token in attrs for token in ("transl_except=", "translation_exception=", "exception=")):
                exceptions.add(attrs)
            for token in attrs.split(";"):
                if token.startswith(("transl_table=", "translation_table=")) and token.split("=", 1)[1] != "1":
                    exceptions.add(token)
        sequence = sequences[keys[0]]
        expected = int(locus["search_end"])-int(locus["search_start"])+1
        if len(sequence) != expected:
            raise ValueError("Prepared genomic sequence length does not match its coordinate record")
        partial_intervals, localized_partial, generic_partial_seen, unlocated_generic = _partial_intervals(
            by_tx_rows, int(locus["search_start"]), int(locus["search_end"]))
        partial = (generic_partial_seen or bool(partial_intervals) or
                   any(r.get("partial_start") == "1" or r.get("partial_end") == "1" for r in rows))
        partial_location_known = not partial or (localized_partial and not unlocated_generic)
        result.append(NativeLocus(next(iter(families)), locus["species"], locus["gene_copy_id"],
            locus["contig"], locus["strand"], int(locus["search_start"]), int(locus["search_end"]),
            sequence, tuple(instances), ordered_paths, partial, coding, tuple(sorted(exceptions)),
            partial_intervals, partial_location_known))
    keys = [(r.family, r.species) for r in result]
    if len(keys) != len(set(keys)):
        raise ValueError("V19 exon configurations require one orthologous locus per family and species")
    return tuple(result)
