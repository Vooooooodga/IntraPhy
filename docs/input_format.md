> **V18 legacy documentation / retained reference.** The default V19 model and
> commands are described in [the README](../README.md) and
> [exon_structure_model.md](exon_structure_model.md). Do not interpret the old
> independent-layer analyses as exon configuration inference.

# Input reference

For the standard FASTA/GFF/tree interface, start with [file inputs](inputs.md).
The table interface below is optional for advanced pipelines; users do not need
to create a manifest for ordinary analysis.


## Raw gene manifest

Required TSV columns are `species`, `family_id`, `gene_id`, `genome_fasta` and
`annotation_file`. Resource paths may be absolute or relative to the manifest.
Optional `gene_copy_id` defaults to `gene_id`; `case_id` defaults to `family_id`.
Source and copy-role fields are optional metadata and are never inferred from
gene names. Formal input contains one distinct gene locus per species and family.

Genome FASTA names must match GFF3/GTF sequence IDs. Annotation must preserve
feature parentage and strand. All selected supplied transcripts are retained;
`--transcript-policy canonical` is an explicit preparation restriction.
Search flanks and extension bounds do not redefine the annotated gene boundary.

## Rooted tree

A Newick tree can be supplied to `build-case` or `import-orthofinder`. Species
labels must match species file stems (or species in the optional manifest). The prepared `species_tree.tsv` uses
`node_id`, `parent_id`, `label` and optional `branch_length`. The root has no parent.
Topology is validated independently of the inference method. CTMC requires finite
nonnegative non-root branch lengths, unless `--branch-length-mode unit` is explicit.
Parsimony does not require branch lengths for its cost.

## OrthoFinder

The resource manifest requires `species`, `genome_fasta` and
`annotation_file`. An optional `member_id_prefix` declares the exact literal
prefix to remove from that species' OrthoFinder member IDs before annotation
lookup. No prefix is inferred from underscores or species names. The original
member ID is retained in the mapping report, and conflicting mappings from the
original and stripped ID remain ambiguous.

Supply an upstream orthogroup and the species tree:

```bash
intraphy import-orthofinder \
  --orthofinder-dir OrthoFinder/Results_run \
  --orthogroup OG0001 --genome-manifest genomes.tsv \
  --species-tree species.nwk --prune-species-tree \
  --output-dir prepared/OG0001
```

Exact annotation ID resolution determines distinct gene loci. Multiple isoform
IDs at one locus do not become extra gene copies. By default, any missing,
ambiguous, or multi-locus species mapping is recorded and stops the import.
`--on-unresolved exclude` retains only species with a unique annotated locus;
the excluded species and reasons are written to `excluded_families.tsv`, and at
least two species must remain. It does not assign a state or select among
paralogs. `--prune-species-tree` requires `--species-tree` and trims extra tips
to the retained manifest species; missing required tips remain an error. The
root, branch lengths, and retained paths are preserved. By default, the full
supplied tree is retained, including tips outside the resource manifest. Use
the pruning flag when the prepared tree must match the retained species panel.
This import does not infer orthology.

## Frozen structural matrix

Use a generated `structural_site_matrix.tsv` rather than fabricating negative
states. Schema-v3 rows identify `family_id`, `layer`, `site_id`, `species`,
`state`, `state_0`, `state_1`, annotation view, discovery rule, observation mask,
applicability and evidence. Each character requires exactly one row per tree tip.
`linked_group_id` records character dependence, not a compound-event conclusion.

An independent complete-universe catalogue can contain observed all-zero sites
and explicit unknown tips with missing masks and reasons. Unknown tips must not
be replaced with zero to satisfy a validator. A discovered positive-only catalogue
cannot be relabeled independent after filtering.

Internal tables for occurrences, transcript paths, matches, element membership
and position projection retain source coordinates. Their definitions and readers
are in `observations/schema.py`, `mapping/fields.py` and `storage/`. Stable
character identity is conditional on frozen input and algorithm settings; IDs
are not permanent biological accession numbers.
