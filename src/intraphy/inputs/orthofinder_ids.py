"""Parse OrthoFinder member IDs and map them through GFF ancestry."""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from ..preparation.annotation_index import load_annotation_index


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


def single_gene(value):
    return [part.strip() for part in str(value or "").split(",") if part.strip()]


def name_key(value):
    return str(value or "").strip().replace(" ", "_")


def split_attr_values(value):
    tokens = []
    for part in str(value or "").replace(",", ";").split(";"):
        part = part.strip().strip('"')
        if part:
            tokens.append(part)
    return tokens


def token_variants(token):
    token = str(token or "").strip().strip(">").strip('"')
    if not token:
        return set()
    variants = {token}
    for prefix in (
        "gene-", "rna-", "cds-", "transcript-", "protein-",
        "gene:", "transcript:", "CDS:", "cds:", "protein:",
    ):
        if token.startswith(prefix) and len(token) > len(prefix):
            variants.add(token[len(prefix):])
    return {variant for variant in variants if variant}


def member_token_groups(member_id):
    parts = str(member_id or "").strip().lstrip(">").split()
    if not parts:
        return set(), set()
    metadata = set()
    for part in parts[1:]:
        for key in ("gene", "gene_id", "transcript", "transcript_id", "protein", "protein_id"):
            for separator in ("=", ":"):
                prefix = key + separator
                if part.startswith(prefix):
                    metadata.update(split_attr_values(part[len(prefix):]))
    return {parts[0]}, metadata


def sequence_ids_file(orthofinder_dir):
    root = Path(orthofinder_dir)
    candidates = [
        root / "WorkingDirectory" / "SequenceIDs.txt",
        root / "SequenceIDs.txt",
    ]
    for path in candidates:
        if path.exists():
            return path
    return None


def read_sequence_id_aliases(orthofinder_dir):
    path = sequence_ids_file(orthofinder_dir)
    if not path:
        return {}
    aliases = defaultdict(set)
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line or ":" not in line:
                continue
            left, right = line.split(":", 1)
            left = left.strip()
            right = right.strip()
            if not left or not right:
                continue
            source_id = right.split()[0]
            for token in (left, right, source_id, source_id.replace(":", "_")):
                aliases[token].add(right)
    return aliases


def feature_index_tokens(feature):
    attrs = feature.get("attrs", {})
    values = [feature.get("id", "")]
    keys = ["ID", "Alias", "transcript_id", "protein_id"]
    if feature.get("type", "").lower() == "gene":
        keys.append("gene_id")
    for key in keys:
        values.extend(split_attr_values(attrs.get(key, "")))
    exact = {value for value in values if value}
    version = str(attrs.get("version", "")).strip()
    if version:
        for key in keys:
            if key == "Alias":
                continue
            for value in split_attr_values(attrs.get(key, "")):
                if "." not in value or not value.rsplit(".", 1)[-1].isdigit():
                    exact.add(f"{value}.{version}")
    derived = set()
    for value in exact:
        derived.update(token_variants(value))
    for key in ("Name", "gene_name"):
        derived.update(split_attr_values(attrs.get(key, "")))
    if feature.get("name"):
        derived.add(feature["name"])
    return exact, derived - exact


def build_gff_locus_index(annotation_file):
    index = load_annotation_index(annotation_file)
    features = index.rows
    features_by_id = index.by_id

    gene_ids_to_loci = defaultdict(set)
    gene_alias_to_loci = defaultdict(set)
    gene_names_to_loci = defaultdict(set)
    for feature in features:
        if feature.get("type", "").lower() != "gene":
            continue
        locus = feature.get("id") or feature.get("name")
        if not locus:
            continue
        gene_ids_to_loci[locus].add(locus)
        exact, derived = feature_index_tokens(feature)
        for token in exact:
            gene_alias_to_loci[token].add(locus)
        for token in derived:
            gene_names_to_loci[token].add(locus)

    locus_cache = {}

    def feature_loci(feature, visiting=None):
        visiting = set(visiting or set())
        key = (
            feature.get("seqid", ""), feature.get("type", ""),
            feature.get("start", ""), feature.get("end", ""),
            feature.get("id", ""), feature.get("parent", ""),
        )
        if key in locus_cache:
            return set(locus_cache[key])
        if key in visiting:
            return set()
        visiting.add(key)
        if feature.get("type", "").lower() == "gene":
            locus = feature.get("id") or feature.get("name")
            loci = {locus} if locus else set()
            locus_cache[key] = tuple(sorted(loci))
            return loci
        loci = set()
        for parent in feature.get("parents", []):
            for parent_feature in features_by_id.get(parent, []):
                loci.update(feature_loci(parent_feature, visiting))
        # Explicit parent links define the locus; dangling links stay unresolved.
        attrs = feature.get("attrs", {})
        if not loci and not attrs.get("Parent"):
            for gene_key in ("gene_id", "gene"):
                for value in split_attr_values(attrs.get(gene_key, "")):
                    loci.update(
                        gene_ids_to_loci.get(value)
                        or gene_alias_to_loci.get(value)
                        or gene_names_to_loci.get(value, set())
                    )
        locus_cache[key] = tuple(sorted(loci))
        return loci

    exact_index = defaultdict(set)
    normalized_index = defaultdict(set)
    derived_index = defaultdict(set)
    for feature in features:
        loci = feature_loci(feature)
        if not loci:
            continue
        exact, derived = feature_index_tokens(feature)
        for token in exact:
            exact_index[token].update(loci)
            if ":" in token:
                normalized_index[token.replace(":", "_")].update(loci)
        for token in derived:
            derived_index[token].update(loci)
    for token, loci in gene_ids_to_loci.items():
        exact_index[token] = set(loci)
    return exact_index, normalized_index, derived_index


def resolve_members_to_locus(
    members,
    annotation_file,
    sequence_aliases=None,
    *,
    locus_index=None,
    member_id_prefix=None,
):
    sequence_aliases = sequence_aliases or {}
    exact_index, normalized_index, derived_index = (
        locus_index if locus_index is not None else build_gff_locus_index(annotation_file)
    )
    member_rows = []
    species_loci = set()
    for member in members:
        token_hits = []
        loci = set()
        headers, has_source_aliases = member_lookup_headers(
            member, member_id_prefix, sequence_aliases
        )
        for header in sorted(headers):
            exact_tokens, metadata_tokens = member_token_groups(header)
            header_loci = set()
            for tokens in (exact_tokens, metadata_tokens):
                for token in sorted(tokens):
                    hits = set(exact_index.get(token, set()))
                    if not has_source_aliases:
                        hits.update(normalized_index.get(token, set()))
                    if hits:
                        token_hits.append(f"{token}:{','.join(sorted(hits))}")
                        header_loci.update(hits)
                if header_loci:
                    break
            if not header_loci:
                fallback = set().union(*(token_variants(token) for token in exact_tokens | metadata_tokens))
                for token in sorted(fallback):
                    hits = exact_index.get(token) or derived_index.get(token, set())
                    if hits:
                        token_hits.append(f"{token}:{','.join(sorted(hits))}")
                        header_loci.update(hits)
            loci.update(header_loci)
        if not loci:
            member_rows.append({
                "source_member_id": member, "gene_id": "NA",
                "matched_tokens": "NA", "mapping_status": "unresolved_member",
            })
            continue
        if len(loci) > 1:
            member_rows.append({
                "source_member_id": member, "gene_id": ";".join(sorted(loci)),
                "matched_tokens": ";".join(token_hits), "mapping_status": "ambiguous_member",
            })
            continue
        locus = next(iter(loci))
        species_loci.add(locus)
        member_rows.append({
            "source_member_id": member, "gene_id": locus,
            "matched_tokens": ";".join(token_hits), "mapping_status": "mapped",
        })
    statuses = {row["mapping_status"] for row in member_rows}
    if "ambiguous_member" in statuses:
        return None, "ambiguous_member_id", member_rows
    if "unresolved_member" in statuses:
        return None, "unresolved_member_id", member_rows
    if len(species_loci) != 1:
        return None, "orthofinder_members_map_to_multiple_loci", member_rows
    return next(iter(species_loci)), "only_gene_unique", member_rows


def members_for_species_from_txt(all_members, resource_rows, sequence_aliases):
    species_members = {resource["species"]: [] for resource in resource_rows}
    member_hits = defaultdict(list)
    for resource in resource_rows:
        locus_index = build_gff_locus_index(resource["annotation_file"])
        for member in all_members:
            locus_id, mapping_status, _ = resolve_members_to_locus(
                [member], resource["annotation_file"], sequence_aliases,
                locus_index=locus_index, member_id_prefix=resource.get("member_id_prefix"),
            )
            if locus_id is not None:
                species_members[resource["species"]].append(member)
                member_hits[member].append(resource["species"])
            elif mapping_status == "ambiguous_member_id":
                species_members[resource["species"]].append(member)
                member_hits[member].append(resource["species"])
    ambiguous = {member: species for member, species in member_hits.items() if len(species) > 1}
    if ambiguous:
        details = "; ".join(f"{member}=>{','.join(sorted(species))}" for member, species in sorted(ambiguous.items()))
        raise SystemExit(f"Orthogroups.txt member-to-species mapping is ambiguous: {details}")
    return species_members
