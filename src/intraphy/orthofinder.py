"""Import one single-copy orthogroup without re-inferring gene homology."""

from __future__ import annotations

from pathlib import Path
from collections import defaultdict
from intraphy.preparation.annotation_index import load_annotation_index

from .io import read_tsv, write_tsv
from .preprocess import read_annotation
from .inputs.orthofinder_ids import member_lookup_headers


MANIFEST_FIELDS = [
    "case_id",
    "species",
    "family_id",
    "gene_id",
    "gene_copy_id",
    "genome_fasta",
    "annotation_file",
    "assembly",
    "annotation",
    "source_url",
    "release",
    "orthofinder_member_ids",
    "orthofinder_member_count",
    "orthofinder_locus_count",
    "orthofinder_mapping_status",
    "notes",
]


def _orthogroups_table(orthofinder_dir):
    root = Path(orthofinder_dir)
    candidates = [
        root / "Orthogroups" / "Orthogroups.tsv",
        root / "Orthogroups.tsv",
        root / "Orthogroups" / "Orthogroups.txt",
        root / "Orthogroups.txt",
    ]
    for path in candidates:
        if path.exists():
            return path
    raise SystemExit(f"Orthogroups.tsv or Orthogroups.txt not found below {root}")


def _read_orthogroups(path):
    path = Path(path)
    if path.suffix.lower() == ".tsv":
        return read_tsv(path)
    rows = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line or ":" not in line:
                continue
            group_id, members = line.split(":", 1)
            rows.append({"Orthogroup": group_id.strip(), "__members__": members.strip()})
    return rows


def _single_gene(value):
    genes = [part.strip() for part in str(value or "").split(",") if part.strip()]
    return genes


def _name_key(value):
    return str(value or "").strip().replace(" ", "_")


def _split_attr_values(value):
    tokens = []
    for part in str(value or "").replace(",", ";").split(";"):
        part = part.strip().strip('"')
        if part:
            tokens.append(part)
    return tokens


def _token_variants(token):
    token = str(token or "").strip().strip(">").strip('"')
    if not token:
        return set()
    variants = {token}
    for prefix in (
        "gene-", "rna-", "cds-", "transcript-", "protein-",
        "gene:", "transcript:", "CDS:", "cds:", "protein:",
    ):
        if token.startswith(prefix) and len(token) > len(prefix):
            variants.add(token[len(prefix) :])
    return {variant for variant in variants if variant}


def _member_token_groups(member_id):
    parts = str(member_id or "").strip().lstrip(">").split()
    if not parts:
        return set(), set()
    metadata = set()
    for part in parts[1:]:
        for key in ("gene", "gene_id", "transcript", "transcript_id", "protein", "protein_id"):
            for separator in ("=", ":"):
                prefix = key + separator
                if part.startswith(prefix):
                    metadata.update(_split_attr_values(part[len(prefix):]))
    return {parts[0]}, metadata


def _sequence_ids_file(orthofinder_dir):
    root = Path(orthofinder_dir)
    candidates = [
        root / "WorkingDirectory" / "SequenceIDs.txt",
        root / "SequenceIDs.txt",
    ]
    for path in candidates:
        if path.exists():
            return path
    return None


def _read_sequence_id_aliases(orthofinder_dir):
    path = _sequence_ids_file(orthofinder_dir)
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


def _feature_index_tokens(feature):
    attrs = feature.get("attrs", {})
    values = [feature.get("id", "")]
    keys = ["ID", "Alias", "transcript_id", "protein_id"]
    if feature.get("type", "").lower() == "gene":
        keys.append("gene_id")
    for key in keys:
        values.extend(_split_attr_values(attrs.get(key, "")))
    exact = {value for value in values if value}
    version = str(attrs.get("version", "")).strip()
    if version:
        for key in keys:
            if key == "Alias":
                continue
            for value in _split_attr_values(attrs.get(key, "")):
                if "." not in value or not value.rsplit(".", 1)[-1].isdigit():
                    exact.add(f"{value}.{version}")
    derived = set()
    for value in exact:
        derived.update(_token_variants(value))
    for key in ("Name", "gene_name"):
        derived.update(_split_attr_values(attrs.get(key, "")))
    if feature.get("name"):
        derived.add(feature["name"])
    return exact, derived - exact


def _build_gff_locus_index(annotation_file):
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
        exact, derived = _feature_index_tokens(feature)
        for token in exact:
            gene_alias_to_loci[token].add(locus)
        for token in derived:
            gene_names_to_loci[token].add(locus)

    locus_cache = {}

    def feature_loci(feature, visiting=None):
        visiting = set(visiting or set())
        key = (
            feature.get("seqid", ""),
            feature.get("type", ""),
            feature.get("start", ""),
            feature.get("end", ""),
            feature.get("id", ""),
            feature.get("parent", ""),
        )
        if key in locus_cache:
            return set(locus_cache[key])
        if key in visiting:
            return set()
        visiting.add(key)
        ftype = feature.get("type", "").lower()
        if ftype == "gene":
            locus = feature.get("id") or feature.get("name")
            loci = {locus} if locus else set()
            locus_cache[key] = tuple(sorted(loci))
            return loci
        loci = set()
        parents = feature.get("parents", [])
        for parent in parents:
            for parent_feature in features_by_id.get(parent, []):
                loci.update(feature_loci(parent_feature, visiting))
        # Explicit parent links define the locus; dangling links stay unresolved.
        attrs = feature.get("attrs", {})
        if not loci and not attrs.get("Parent"):
            for gene_key in ("gene_id", "gene"):
                for value in _split_attr_values(attrs.get(gene_key, "")):
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
        exact, derived = _feature_index_tokens(feature)
        for token in exact:
            exact_index[token].update(loci)
            if ":" in token:
                normalized_index[token.replace(":", "_")].update(loci)
        for token in derived:
            derived_index[token].update(loci)
    for token, loci in gene_ids_to_loci.items():
        exact_index[token] = set(loci)
    return exact_index, normalized_index, derived_index


def _resolve_members_to_locus(
    members,
    annotation_file,
    sequence_aliases=None,
    *,
    locus_index=None,
    member_id_prefix=None,
):
    sequence_aliases = sequence_aliases or {}
    exact_index, normalized_index, derived_index = (
        locus_index if locus_index is not None else _build_gff_locus_index(annotation_file)
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
            exact_tokens, metadata_tokens = _member_token_groups(header)
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
                fallback = set().union(*(_token_variants(token) for token in exact_tokens | metadata_tokens))
                for token in sorted(fallback):
                    hits = exact_index.get(token) or derived_index.get(token, set())
                    if hits:
                        token_hits.append(f"{token}:{','.join(sorted(hits))}")
                        header_loci.update(hits)
            loci.update(header_loci)
        if not loci:
            member_rows.append(
                {
                    "source_member_id": member,
                    "gene_id": "NA",
                    "matched_tokens": "NA",
                    "mapping_status": "unresolved_member",
                }
            )
            continue
        if len(loci) > 1:
            member_rows.append(
                {
                    "source_member_id": member,
                    "gene_id": ";".join(sorted(loci)),
                    "matched_tokens": ";".join(token_hits),
                    "mapping_status": "ambiguous_member",
                }
            )
            continue
        locus = next(iter(loci))
        species_loci.add(locus)
        member_rows.append(
            {
                "source_member_id": member,
                "gene_id": locus,
                "matched_tokens": ";".join(token_hits),
                "mapping_status": "mapped",
            }
        )
    statuses = {row["mapping_status"] for row in member_rows}
    if "ambiguous_member" in statuses:
        return None, "ambiguous_member_id", member_rows
    if "unresolved_member" in statuses:
        return None, "unresolved_member_id", member_rows
    if len(species_loci) != 1:
        return None, "orthofinder_members_map_to_multiple_loci", member_rows
    return next(iter(species_loci)), "only_gene_unique", member_rows


def _members_for_species_from_txt(all_members, resource_rows, sequence_aliases):
    species_members = {resource["species"]: [] for resource in resource_rows}
    member_hits = defaultdict(list)
    for resource in resource_rows:
        locus_index = _build_gff_locus_index(resource["annotation_file"])
        for member in all_members:
            locus_id, mapping_status, _ = _resolve_members_to_locus(
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


def _write_species_tree(tree_path, output_path):
    path = Path(tree_path)
    if path.suffix.lower() == ".tsv":
        rows = read_tsv(path, ["node_id", "parent_id", "label"])
        fields = ["node_id", "parent_id", "label", "branch_length"]
        write_tsv(output_path, rows, fields)
        return
    from .inputs.species_tree import read_species_tree_rows

    rows = read_species_tree_rows(path)
    for row in rows:
        if not row["parent_id"]:
            row["branch_length"] = 0.0
        elif row["branch_length"] is None or row["branch_length"] == "":
            row["branch_length"] = "NA"
    write_tsv(output_path, rows, ["node_id", "parent_id", "label", "branch_length"])


def import_orthofinder(
    orthofinder_dir,
    orthogroup,
    genome_manifest,
    output_dir,
    species_tree=None,
    *,
    on_unresolved="error",
    prune_species_tree=False,
):
    """Require one actual gene locus per selected species, retaining all source IDs.

    Multiple protein or transcript members may represent the same gene locus.
    Every member must resolve uniquely; members from distinct loci are excluded.
    """
    if on_unresolved not in {"error", "exclude"}:
        raise SystemExit("--on-unresolved must be 'error' or 'exclude'")
    if prune_species_tree and not species_tree:
        raise SystemExit("--prune-species-tree requires --species-tree")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    resource_rows = read_tsv(genome_manifest, ["species", "genome_fasta", "annotation_file"])
    table = _orthogroups_table(orthofinder_dir)
    orthogroups = _read_orthogroups(table)
    sequence_aliases = _read_sequence_id_aliases(orthofinder_dir)
    selected = [row for row in orthogroups if row.get("Orthogroup") == orthogroup]
    if len(selected) != 1:
        raise SystemExit(f"orthogroup {orthogroup} was not found exactly once in {table}")
    row = selected[0]
    manifest_rows = []
    excluded = []
    mapping_rows = []
    txt_members_by_species = None
    if "__members__" in row:
        txt_members_by_species = _members_for_species_from_txt(_single_gene(row["__members__"].replace(" ", ",")), resource_rows, sequence_aliases)
    for resource in resource_rows:
        if txt_members_by_species is None:
            species_column = None
            for column in row:
                if column != "Orthogroup" and _name_key(column) == _name_key(resource["species"]):
                    species_column = column
                    break
            if species_column is None:
                raise SystemExit(f"selected species {resource['species']} from genome manifest is missing from Orthofinder table")
            value = row.get(species_column, "")
            members = _single_gene(value)
        else:
            members = txt_members_by_species.get(resource["species"], [])
        if not members:
            excluded.append(
                {
                    "family_id": orthogroup,
                    "species": resource["species"],
                    "member_count": 0,
                    "locus_count": 0,
                    "source_member_ids": "NA",
                    "gene_ids": "NA",
                    "reason": "orthofinder_species_has_no_member",
                }
            )
            continue
        locus_id, mapping_status, member_rows = _resolve_members_to_locus(
            members,
            resource["annotation_file"],
            sequence_aliases,
            member_id_prefix=resource.get("member_id_prefix"),
        )
        locus_ids = sorted({
            locus for item in member_rows if item["gene_id"] != "NA"
            for locus in item["gene_id"].split(";")
        })
        for item in member_rows:
            mapping_rows.append(
                {
                    "family_id": orthogroup,
                    "species": resource["species"],
                    "source_member_id": item["source_member_id"],
                    "gene_id": item["gene_id"],
                    "matched_tokens": item["matched_tokens"],
                    "mapping_status": item["mapping_status"],
                }
            )
        if locus_id is None:
            excluded.append(
                {
                    "family_id": orthogroup,
                    "species": resource["species"],
                    "member_count": len(members),
                    "locus_count": len(locus_ids),
                    "source_member_ids": ";".join(members),
                    "gene_ids": ";".join(locus_ids) or "NA",
                    "reason": mapping_status,
                }
            )
            continue
        manifest_rows.append(
            {
                "case_id": orthogroup,
                "species": resource["species"],
                "family_id": orthogroup,
                "gene_id": locus_id,
                "gene_copy_id": locus_id,
                "genome_fasta": resource["genome_fasta"],
                "annotation_file": resource["annotation_file"],
                "assembly": resource.get("assembly", "NA"),
                "annotation": resource.get("annotation", "NA"),
                "source_url": resource.get("source_url", "NA"),
                "release": resource.get("release", "NA"),
                "orthofinder_member_ids": ";".join(members),
                "orthofinder_member_count": len(members),
                "orthofinder_locus_count": 1,
                "orthofinder_mapping_status": mapping_status,
                "notes": "imported_from_orthofinder_only_gene_unique_locus",
            }
        )
    if mapping_rows:
        write_tsv(
            output_dir / "orthofinder_member_mapping.tsv",
            mapping_rows,
            ["family_id", "species", "source_member_id", "gene_id", "matched_tokens", "mapping_status"],
        )
    if excluded:
        write_tsv(
            output_dir / "excluded_families.tsv",
            excluded,
            ["family_id", "species", "member_count", "locus_count", "source_member_ids", "gene_ids", "reason"],
        )
        if on_unresolved == "error":
            raise SystemExit(f"orthogroup {orthogroup} does not map to exactly one annotated gene locus in every species")
    if on_unresolved == "exclude" and len(manifest_rows) < 2:
        raise SystemExit("--on-unresolved exclude requires at least two species with uniquely mapped loci")
    if not manifest_rows:
        raise SystemExit(f"orthogroup {orthogroup} contains no usable species")
    if species_tree:
        if prune_species_tree:
            from .inputs.species_tree import prune_species_tree as prune_tree
            try:
                prune_tree(
                    species_tree,
                    output_dir / "species_tree.tsv",
                    [row["species"] for row in manifest_rows],
                    prune_extra_tips=True,
                )
            except ValueError as exc:
                raise SystemExit(str(exc)) from exc
        else:
            _write_species_tree(species_tree, output_dir / "species_tree.tsv")
    write_tsv(output_dir / "manifest.tsv", manifest_rows, MANIFEST_FIELDS)
    return manifest_rows
