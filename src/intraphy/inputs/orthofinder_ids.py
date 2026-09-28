"""Resolve optional, explicitly declared OrthoFinder member-ID prefixes."""
from __future__ import annotations


def member_id_candidates(member_id, declared_prefix=None):
    """Return the source ID and, only on an exact prefix match, its suffix."""
    member_id = str(member_id or "")
    prefix = str(declared_prefix or "")
    if prefix in {"NA", "."} or not prefix:
        return (member_id,)
    if member_id.startswith(prefix) and len(member_id) > len(prefix):
        return member_id, member_id[len(prefix):]
    return (member_id,)


def member_alias_headers(member_ids, aliases):
    """Collect aliases, preserving the established full-ID lookup priority."""
    headers = set()
    for member_id in member_ids:
        direct = set(aliases.get(member_id, set()))
        if direct:
            headers.update(direct)
            continue
        primary = str(member_id or "").strip().lstrip(">").split()
        if primary:
            headers.update(aliases.get(primary[0], set()))
    return headers


def member_lookup_headers(member_id, declared_prefix, aliases):
    """Resolve raw and alias headers without losing exact-prefix conflicts."""
    member_ids = member_id_candidates(member_id, declared_prefix)
    source_headers = member_alias_headers(member_ids, aliases)
    alias_ids = set()
    alias_prefix_matched = False
    for header in source_headers:
        candidates = member_id_candidates(header, declared_prefix)
        alias_ids.update(candidates)
        alias_prefix_matched = alias_prefix_matched or len(candidates) > 1

    member_prefix_matched = len(member_ids) > 1
    if member_prefix_matched:
        return source_headers | alias_ids | set(member_ids), bool(source_headers)
    if alias_prefix_matched:
        return alias_ids, bool(source_headers)
    return source_headers or {str(member_id)}, bool(source_headers)
