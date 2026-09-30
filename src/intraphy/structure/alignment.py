"""Audited common coordinates for local exon configurations.

The genomic MSA is a correspondence hypothesis. Exon matches and genomic gap
calls are separately qualified. Weak intronic alignment is not called homology.
Existing protein-projection evidence is retained in the preparation directory.
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from collections import defaultdict
import re

from ..aligners.recorded import run_recorded
from ..aligners.columns import revcomp
from ..storage.fasta import parse_fasta
from ..storage.tabular import write_tsv
from .native import NativeLocus
from .types import ExonSpan


@dataclass(frozen=True)
class FamilyAlignment:
    loci: tuple[NativeLocus, ...]
    rows: dict[str, str]
    offsets: dict[str, int]
    columns: dict[str, tuple[int, ...]]
    exons: dict[str, ExonSpan]
    anomalies: dict[str, tuple[str, ...]]
    matches: tuple[dict, ...]
    reuse_provenance: dict | None = None


def _write_fasta(path: Path, records: dict[str, str]):
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in records.items()))


def align_family(loci: tuple[NativeLocus, ...], output_dir: str | Path, *,
                 timeout: int = 600, max_bases: int = 100000,
                 threads: int = 1, evidence_dir: str | Path | None = None) -> FamilyAlignment:
    if type(threads) is not int or threads < 1:
        raise ValueError("threads must be a positive integer")
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    ordered = tuple(sorted(loci, key=lambda l: l.species))
    ids, records, offsets = {}, {}, {}
    anomalies = defaultdict(set)
    for i, locus in enumerate(ordered):
        identifier = f"taxon_{i:05d}"
        ids[identifier] = locus.species
        # Annotation incompleteness must not remove already supplied terminal DNA.
        a, b = 0, len(locus.sequence)
        records[identifier] = locus.sequence[a:b]
        offsets[locus.species] = a
    queries, query_exon = {}, {}
    for locus in ordered:
        for exon in locus.exons:
            identifier = f"exon_{len(queries):06d}"
            a, b = locus.oriented_interval(exon)
            queries[identifier] = locus.sequence[a:b]
            query_exon[identifier] = exon
    window_rows = [{"species": l.species, "source_contig": l.contig,
        "oriented_start0": 0, "oriented_end0": len(l.sequence),
        "available_bases": len(l.sequence), "selected_bases": len(l.sequence),
        "source_first_base0": l.source_position(0),
        "source_last_base0": l.source_position(len(l.sequence)-1),
        "strand": l.strand, "policy": "entire_supplied_or_extracted_locus",
        "annotation_based_cropping": False, "searched_bases": 0,
        "unsearched_available_bases": len(l.sequence), "alignment_status": "not_completed"}
        for l in ordered]
    if window_rows:
        write_tsv(directory / "search_windows.tsv", window_rows, list(window_rows[0]))
    if not any(locus.exons for locus in ordered):
        raise ValueError("no_annotated_exons: no target structure can be defined")
    if max(map(len, records.values())) > max_bases or sum(map(len, records.values())) > max_bases*30:
        for record in window_rows:
            record["alignment_status"] = "budget_exceeded"
        write_tsv(directory / "search_windows.tsv", window_rows, list(window_rows[0]))
        raise ValueError("locus_alignment_budget_exceeded: no genomic sequence was silently trimmed")
    raw = directory / "genomic_loci.fa"
    _write_fasta(raw, records)
    query_file = directory / "exon_queries.fa"
    _write_fasta(query_file, queries)
    reuse_provenance = None
    paf_text = None
    if evidence_dir is not None:
        from .alignment_evidence import load_alignment_evidence
        aligned, paf_text, reuse_provenance = load_alignment_evidence(
            evidence_dir, records, queries, ids)
    elif len(records) == 1:
        aligned = records
    else:
        output = run_recorded(["mafft", "--auto", "--inputorder", "--thread", str(threads),
                               "--threadit", "0", str(raw.resolve())],
                              directory, "genomic_mafft", timeout)
        aligned = parse_fasta(output)
    if set(aligned) != set(records) or len({len(s) for s in aligned.values()}) != 1:
        raise ValueError("MAFFT output changed record identifiers or alignment lengths")
    for identifier, sequence in aligned.items():
        if sequence.replace("-", "").upper() != records[identifier].upper():
            raise ValueError("MAFFT output changed genomic residues")
    for record in window_rows:
        record.update(alignment_status="completed", searched_bases=record["selected_bases"], unsearched_available_bases=0)
    write_tsv(directory / "search_windows.tsv", window_rows, list(window_rows[0]))
    rows = {ids[k]: v.upper() for k, v in aligned.items()}
    columns = {s: tuple(i for i, b in enumerate(seq) if b != "-") for s, seq in rows.items()}
    spans = {}
    for locus in ordered:
        for exon in locus.exons:
            a, b = locus.oriented_interval(exon)
            ia, ib = a-offsets[locus.species], b-offsets[locus.species]
            spans[exon.id] = ExonSpan(columns[locus.species][ia], columns[locus.species][ib-1]+1)
    # An independent nucleotide search retains inversions/competing local copies.
    if evidence_dir is None:
        paf = run_recorded(["minimap2", "-c", "-x", "asm20", "--secondary=yes", "-N", "50",
                            "-t", str(threads),
                            str(raw.resolve()), str(query_file.resolve())], directory, "exon_minimap2", timeout)
        paf_text = paf.read_text()
    matches = []
    hit_groups = defaultdict(list)
    for line in paf_text.splitlines():
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        q, target = fields[0], fields[5]
        if q not in query_exon or target not in ids:
            raise ValueError("Unexpected minimap2 record identifier")
        coverage = (int(fields[3])-int(fields[2]))/max(1, int(fields[1]))
        identity = int(fields[9])/max(1, int(fields[10]))
        item = {"query_exon": query_exon[q].id, "target_species": ids[target],
                "query_start0": int(fields[2]), "query_end0": int(fields[3]), "strand": fields[4],
                "target_start0": int(fields[7]), "target_end0": int(fields[8]),
                "coverage": coverage, "identity": identity, "mapq": int(fields[11]),
                "cigar": next((v[5:] for v in fields[12:] if v.startswith("cg:Z:")), "unreported")}
        matches.append(item)
        if coverage >= .8 and identity >= .7:
            hit_groups[(q, target)].append(item)
    for (q, target), hits in hit_groups.items():
        if any(h["strand"] == "-" for h in hits):
            anomalies[query_exon[q].id].add("reverse_correspondence_candidate")
        locations = {(h["target_start0"], h["target_end0"], h["strand"]) for h in hits}
        if len(locations) > 1:
            anomalies[query_exon[q].id].add("competing_copy_candidates")
    # Exact short exon repetitions may lie below minimap2's seeding limit.
    for q, seq in queries.items():
        if len(seq) < 9 or set(seq) - set("ACGT"):
            continue
        for target_seq in records.values():
            if len(list(re.finditer(f"(?={seq})", target_seq))) > 1:
                anomalies[query_exon[q].id].add("repeated_exact_exon_sequence")
            reverse = revcomp(seq)
            if reverse != seq and reverse in target_seq and seq not in target_seq:
                anomalies[query_exon[q].id].add("reverse_correspondence_candidate")
    (directory / "genomic_alignment.fa").write_text("".join(f">{s}\n{row}\n" for s, row in rows.items()))
    return FamilyAlignment(ordered, rows, offsets, columns, spans,
                           {k: tuple(sorted(v)) for k, v in anomalies.items()}, tuple(matches),
                           reuse_provenance)


def paired_support(row_a: str, row_b: str, start: int, end: int):
    pairs = [(a, b) for a, b in zip(row_a[start:end], row_b[start:end]) if a in "ACGT" and b in "ACGT"]
    return len(pairs), sum(a == b for a, b in pairs)/max(1, len(pairs))


def gap_supported(rows: dict[str, str], species: str, start: int, end: int,
                  flank: int = 12, identity: float = .8) -> bool:
    target = rows[species]
    if start < flank or end+flank > len(target) or set(target[start:end]) != {"-"}:
        return False
    if any(c not in "ACGT" for c in target[start-flank:start]+target[end:end+flank]):
        return False
    for other, row in rows.items():
        if other == species or not any(c in "ACGT" for c in row[start:end]):
            continue
        left = paired_support(target, row, start-flank, start)
        right = paired_support(target, row, end, end+flank)
        if left[0] == flank and right[0] == flank and min(left[1], right[1]) >= identity:
            return True
    return False
