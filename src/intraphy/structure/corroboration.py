"""Check nucleotide coordinates against retained family protein projections."""
from __future__ import annotations
from collections import defaultdict
from dataclasses import replace
from pathlib import Path
import re
from ..storage.tabular import iter_tsv
from .alignment import FamilyAlignment


def check_coding_projection(alignment: FamilyAlignment, input_dir: str | Path):
    occurrences = {e.id: (l, e) for l in alignment.loci for e in l.exons}
    support = defaultdict(list)
    anomalies = {k: set(v) for k, v in alignment.anomalies.items()}
    records = []
    for row in iter_tsv(Path(input_dir)/"segment_matches.tsv", optional=True):
        qid, tid = row.get("query_occurrence_id"), row.get("subject_occurrence_id")
        if qid not in occurrences or tid not in occurrences:
            continue
        if row.get("protein_hard_observation_eligible") != "1":
            continue
        ql, qe = occurrences[qid]
        tl, te = occurrences[tid]
        qstart, _ = ql.oriented_interval(qe)
        tstart, _ = tl.oriented_interval(te)
        blocks = row.get("protein_projected_blocks", "")
        checked, agreements = 0, 0
        for block in blocks.split(";"):
            match = re.fullmatch(r"(\d+)-(\d+):(\d+)-(\d+)", block)
            if not match:
                continue
            a, b, c, d = map(int, match.groups())
            if b-a != d-c:
                continue
            for j in range(b-a+1):
                qi = qstart+a-1+j-alignment.offsets[ql.species]
                ti = tstart+c-1+j-alignment.offsets[tl.species]
                if not 0 <= qi < len(alignment.columns[ql.species]) or not 0 <= ti < len(alignment.columns[tl.species]):
                    raise ValueError("Protein projection exceeds its genomic exon")
                checked += 1
                agreements += alignment.columns[ql.species][qi] == alignment.columns[tl.species][ti]
        if checked:
            fraction = agreements/checked
            status = "concordant" if fraction >= .95 else "coordinate_conflict"
            records.append({"query_exon": qid, "target_exon": tid, "checked_bases": checked,
                            "agreeing_bases": agreements, "status": status})
            for eid in (qid, tid):
                support[eid].append(status)
                if status == "coordinate_conflict":
                    anomalies.setdefault(eid, set()).add("protein_genomic_coordinate_conflict")
    return replace(alignment, anomalies={k: tuple(sorted(v)) for k, v in anomalies.items()}), records
