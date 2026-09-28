"""Focused integration coverage for OrthoFinder panel import controls."""
import csv
import io
import tempfile
import unittest
from contextlib import redirect_stderr
from pathlib import Path
from unittest.mock import patch

from intraphy.cli import main as cli_main
from intraphy.io import read_tsv
from intraphy.topology import SpeciesTree


def _write_gff(path, loci):
    lines = []
    start = 1
    for gene_id, transcript_ids in loci:
        end = start + 99
        lines.append(
            f"chr1\ttest\tgene\t{start}\t{end}\t.\t+\t.\tID={gene_id}"
        )
        for transcript_id in transcript_ids:
            lines.append(
                f"chr1\ttest\tmRNA\t{start}\t{end}\t.\t+\t.\t"
                f"ID={transcript_id};Parent={gene_id}"
            )
        start = end + 1
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_import_inputs(
    root,
    species_members,
    loci_by_species,
    *,
    member_prefixes=None,
    sequence_ids=None,
):
    member_prefixes = member_prefixes or {}
    orthofinder_dir = root / "orthofinder"
    orthogroups_dir = orthofinder_dir / "Orthogroups"
    orthogroups_dir.mkdir(parents=True)
    species = list(species_members)
    with (orthogroups_dir / "Orthogroups.tsv").open("w", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Orthogroup", *species])
        writer.writerow([
            "OG0001",
            *[", ".join(species_members[name]) for name in species],
        ])
    if sequence_ids:
        with (orthofinder_dir / "SequenceIDs.txt").open("w", encoding="utf-8") as handle:
            for sequence_id, header in sequence_ids:
                handle.write(f"{sequence_id}: {header}\n")

    manifest = root / "genomes.tsv"
    with manifest.open("w", encoding="utf-8") as handle:
        fields = ["species", "genome_fasta", "annotation_file", "member_id_prefix"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for name in species:
            annotation = root / f"{name}.gff3"
            _write_gff(annotation, loci_by_species.get(name, []))
            writer.writerow({
                "species": name,
                "genome_fasta": root / f"{name}.fa",
                "annotation_file": annotation,
                "member_id_prefix": member_prefixes.get(name, ""),
            })
    return orthofinder_dir, manifest


def _run_import(orthofinder_dir, manifest, output, *options):
    stderr = io.StringIO()
    argv = [
        "import-orthofinder",
        "--orthofinder-dir", str(orthofinder_dir),
        "--orthogroup", "OG0001",
        "--genome-manifest", str(manifest),
        "--output-dir", str(output),
        "--quiet",
        *map(str, options),
    ]
    with patch("intraphy.commands.session.environment_report", return_value={}), \
            redirect_stderr(stderr):
        status = cli_main(argv)
    return status, stderr.getvalue()


class OrthoFinderPanelImportTests(unittest.TestCase):
    def test_exact_prefix_and_sequence_alias_map_to_locus_and_keep_source_id(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            orthofinder, manifest = _write_import_inputs(
                root,
                {"A": ["Apis_GCF_rna-A"], "B": ["tx-B"]},
                {"A": [("gene-A", ["rna-A"])], "B": [("gene-B", ["tx-B"]) ]},
                member_prefixes={"A": "Apis_GCF_"},
            )
            output = root / "direct-prefix"
            status, _ = _run_import(orthofinder, manifest, output)
            self.assertEqual(status, 0)
            manifest_rows = read_tsv(output / "manifest.tsv")
            row_a = next(row for row in manifest_rows if row["species"] == "A")
            self.assertEqual(row_a["gene_id"], "gene-A")
            self.assertEqual(row_a["orthofinder_member_ids"], "Apis_GCF_rna-A")
            mapping = read_tsv(output / "orthofinder_member_mapping.tsv")
            member_a = next(row for row in mapping if row["species"] == "A")
            self.assertEqual(member_a["source_member_id"], "Apis_GCF_rna-A")
            self.assertIn("rna-A:gene-A", member_a["matched_tokens"])

            alias_root = root / "aliased"
            alias_root.mkdir()
            orthofinder, manifest = _write_import_inputs(
                alias_root,
                {"A": ["seq-001"], "B": ["tx-B"]},
                {"A": [("gene-A", ["rna-A"])], "B": [("gene-B", ["tx-B"]) ]},
                member_prefixes={"A": "Apis_GCF_"},
                sequence_ids=[("seq-001", "Apis_GCF_rna-A gene=gene-A")],
            )
            output = alias_root / "sequence-alias"
            status, _ = _run_import(orthofinder, manifest, output)
            self.assertEqual(status, 0)
            row_a = next(row for row in read_tsv(output / "manifest.tsv") if row["species"] == "A")
            self.assertEqual(row_a["gene_id"], "gene-A")
            self.assertEqual(row_a["orthofinder_member_ids"], "seq-001")

    def test_prefix_conflict_between_original_and_stripped_alias_is_ambiguous(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            orthofinder, manifest = _write_import_inputs(
                root,
                {"A": ["0"], "B": ["tx-B"]},
                {
                    "A": [
                        ("gene-original", ["Apis_GCF_rna-A"]),
                        ("gene-stripped", ["rna-A"]),
                    ],
                    "B": [("gene-B", ["tx-B"])],
                },
                member_prefixes={"A": "Apis_GCF_"},
                sequence_ids=[("0", "Apis_GCF_rna-A")],
            )
            output = root / "conflict"
            status, _ = _run_import(orthofinder, manifest, output)
            self.assertEqual(status, 2)
            excluded = read_tsv(output / "excluded_families.tsv")
            row_a = next(row for row in excluded if row["species"] == "A")
            self.assertEqual(row_a["reason"], "ambiguous_member_id")
            mapping = read_tsv(output / "orthofinder_member_mapping.tsv")
            member_a = next(row for row in mapping if row["species"] == "A")
            self.assertEqual(member_a["source_member_id"], "0")
            self.assertEqual(
                set(member_a["gene_id"].split(";")),
                {"gene-original", "gene-stripped"},
            )
            self.assertFalse((output / "manifest.tsv").exists())

    def test_without_matching_prefix_existing_alias_priority_is_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            orthofinder, manifest = _write_import_inputs(
                root,
                {"A": ["shared"], "B": ["tx-B"]},
                {
                    "A": [
                        ("gene-raw", ["shared"]),
                        ("gene-alias", ["tx-alias"]),
                    ],
                    "B": [("gene-B", ["tx-B"])],
                },
                member_prefixes={"A": "No_match_"},
                sequence_ids=[("shared", "tx-alias")],
            )
            output = root / "alias-priority"
            status, _ = _run_import(orthofinder, manifest, output)
            self.assertEqual(status, 0)
            row_a = next(row for row in read_tsv(output / "manifest.tsv") if row["species"] == "A")
            self.assertEqual(row_a["gene_id"], "gene-alias")

    def test_unresolved_default_errors_and_exclude_records_missing_and_multilocus(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            members = {"A": ["tx-A"], "B": [], "C": ["tx-C1", "tx-C2"], "D": ["tx-D"]}
            loci = {
                "A": [("gene-A", ["tx-A"])],
                "B": [("gene-B", ["tx-B"])],
                "C": [("gene-C1", ["tx-C1"]), ("gene-C2", ["tx-C2"])],
                "D": [("gene-D", ["tx-D"])],
            }
            orthofinder, manifest = _write_import_inputs(root, members, loci)

            error_output = root / "default-error"
            status, _ = _run_import(orthofinder, manifest, error_output)
            self.assertEqual(status, 2)
            self.assertFalse((error_output / "manifest.tsv").exists())
            error_rows = read_tsv(error_output / "excluded_families.tsv")
            reasons = {row["species"]: row["reason"] for row in error_rows}
            self.assertEqual(reasons["B"], "orthofinder_species_has_no_member")
            self.assertEqual(reasons["C"], "orthofinder_members_map_to_multiple_loci")

            exclude_output = root / "exclude"
            status, _ = _run_import(
                orthofinder, manifest, exclude_output, "--on-unresolved", "exclude",
            )
            self.assertEqual(status, 0)
            imported = read_tsv(exclude_output / "manifest.tsv")
            self.assertEqual({row["species"] for row in imported}, {"A", "D"})
            excluded = read_tsv(exclude_output / "excluded_families.tsv")
            self.assertEqual(
                {row["species"]: row["reason"] for row in excluded},
                {
                    "B": "orthofinder_species_has_no_member",
                    "C": "orthofinder_members_map_to_multiple_loci",
                },
            )

    def test_exclude_requires_two_uniquely_mapped_species(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            orthofinder, manifest = _write_import_inputs(
                root,
                {"A": ["tx-A"], "B": []},
                {"A": [("gene-A", ["tx-A"])], "B": [("gene-B", ["tx-B"]) ]},
            )
            output = root / "one-retained"
            status, stderr = _run_import(
                orthofinder, manifest, output, "--on-unresolved", "exclude",
            )
            self.assertEqual(status, 2)
            self.assertIn("at least two species", stderr)
            self.assertFalse((output / "manifest.tsv").exists())

    def test_pruning_matches_manifest_and_default_keeps_full_tree(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            orthofinder, manifest = _write_import_inputs(
                root,
                {"A": ["tx-A"], "B": ["tx-B"]},
                {"A": [("gene-A", ["tx-A"])], "B": [("gene-B", ["tx-B"]) ]},
            )
            tree = root / "tree.nwk"
            tree.write_text("(A:0.1,(B:0.2,C:0.3):0.4);\n", encoding="utf-8")

            retained_output = root / "retained-full-tree"
            status, _ = _run_import(
                orthofinder, manifest, retained_output, "--species-tree", tree,
            )
            self.assertEqual(status, 0)
            full_tree = SpeciesTree(read_tsv(retained_output / "species_tree.tsv"))
            self.assertEqual(set(full_tree.leaf_by_label), {"A", "B", "C"})

            pruned_output = root / "pruned-tree"
            status, _ = _run_import(
                orthofinder,
                manifest,
                pruned_output,
                "--species-tree", tree,
                "--prune-species-tree",
            )
            self.assertEqual(status, 0)
            pruned_tree = SpeciesTree(read_tsv(pruned_output / "species_tree.tsv"))
            manifest_species = {row["species"] for row in read_tsv(pruned_output / "manifest.tsv")}
            self.assertEqual(set(pruned_tree.leaf_by_label), manifest_species)

    def test_pruning_rejects_a_tree_missing_a_successful_manifest_species(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            orthofinder, manifest = _write_import_inputs(
                root,
                {"A": ["tx-A"], "B": ["tx-B"]},
                {"A": [("gene-A", ["tx-A"])], "B": [("gene-B", ["tx-B"]) ]},
            )
            tree = root / "tree.nwk"
            tree.write_text("(A:0.1,C:0.2);\n", encoding="utf-8")
            output = root / "missing-tip"
            status, stderr = _run_import(
                orthofinder,
                manifest,
                output,
                "--species-tree", tree,
                "--prune-species-tree",
            )
            self.assertEqual(status, 2)
            self.assertIn("missing requested tips: B", stderr)
            self.assertFalse((output / "manifest.tsv").exists())


if __name__ == "__main__":
    unittest.main()
