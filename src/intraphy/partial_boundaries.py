"""Small shared parser for NCBI GFF3 partial-boundary attributes."""
from __future__ import annotations


def parse_attributes(value):
    if isinstance(value, dict):
        return {str(key): str(item) for key, item in value.items()}
    attributes = {}
    for field in str(value or "").split(";"):
        if "=" in field:
            key, item = field.split("=", 1)
            attributes[key.strip()] = item.strip().strip('"')
    return attributes


def range_value(attributes, genomic_boundary):
    """Return a declared genomic-column range, excluding empty placeholders."""
    value = parse_attributes(attributes).get(f"{genomic_boundary}_range", "")
    return value if range_is_partial(value, genomic_boundary) else ""


def range_declared(attributes, genomic_boundary):
    """Whether a range attribute exists, even if its location is malformed."""
    return bool(parse_attributes(attributes).get(f"{genomic_boundary}_range", ""))


def range_is_partial(value, genomic_boundary=None):
    """Recognize canonical NCBI ranges with a locatable uncertain boundary."""
    parts = str(value or "").split(",")
    if len(parts) != 2 or not all(part.strip() == "." or part.strip().isdigit()
                                  for part in parts):
        return False
    left, right = (part.strip() for part in parts)
    if left == right or (left == right == "."):
        return False
    if left == ".":
        return genomic_boundary != "end" and right.isdigit()
    if right == ".":
        return genomic_boundary != "start" and left.isdigit()
    return int(left) < int(right)


def generic_partial(attributes):
    return parse_attributes(attributes).get("partial", "").lower() in {"1", "true", "yes"}
