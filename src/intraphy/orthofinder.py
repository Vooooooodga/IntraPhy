"""Import one single-copy orthogroup without re-inferring gene homology."""

from __future__ import annotations

from pathlib import Path

from .io import read_tsv, write_tsv
from .inputs.orthofinder_ids import (
    members_for_species_from_txt,
    name_key,
    read_sequence_id_aliases,
    resolve_members_to_locus as _resolve_members_to_locus,
    single_gene,
)


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
    sequence_aliases = read_sequence_id_aliases(orthofinder_dir)
    selected = [row for row in orthogroups if row.get("Orthogroup") == orthogroup]
    if len(selected) != 1:
        raise SystemExit(f"orthogroup {orthogroup} was not found exactly once in {table}")
    row = selected[0]
    manifest_rows = []
    excluded = []
    mapping_rows = []
    txt_members_by_species = None
    if "__members__" in row:
        txt_members_by_species = members_for_species_from_txt(single_gene(row["__members__"].replace(" ", ",")), resource_rows, sequence_aliases)
    for resource in resource_rows:
        if txt_members_by_species is None:
            species_column = None
            for column in row:
                if column != "Orthogroup" and name_key(column) == name_key(resource["species"]):
                    species_column = column
                    break
            if species_column is None:
                raise SystemExit(f"selected species {resource['species']} from genome manifest is missing from Orthofinder table")
            value = row.get(species_column, "")
            members = single_gene(value)
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
