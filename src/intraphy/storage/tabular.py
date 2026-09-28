"""storage / tabular: extracted responsibilities; see docs/architecture.md."""
from __future__ import annotations

import ctypes
from pathlib import Path
import csv
import gzip


def open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def iter_tsv(path, required=None, optional=False):
    """Yield TSV rows while the handle is open; reject malformed row widths."""
    # csv.field_size_limit is process-global. Keep one platform-safe maximum
    # instead of restoring it when this generator closes: TSV readers may be
    # interleaved, and restoring a smaller limit could break another reader.
    csv.field_size_limit((1 << (ctypes.sizeof(ctypes.c_long) * 8 - 1)) - 1)
    path = Path(path)
    if optional and not path.exists():
        return
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [field for field in (required or []) if field not in (reader.fieldnames or [])]
        if missing:
            raise SystemExit(f"{path} missing required fields: {','.join(missing)}")
        if reader.fieldnames and len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError(f"{path}: duplicate TSV column names")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"{path}:{reader.line_num}: wrong number of TSV fields")
            yield row


def read_tsv(path, required=None, optional=False):
    """Materialise rows for callers that rely on repeated traversal."""
    return list(iter_tsv(path, required=required, optional=optional))


def write_tsv(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "NA") for field in fields})
