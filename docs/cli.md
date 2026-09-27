# Command-line guide

## Genomic evidence

`check`, `build-case`, and `derive-tables` are the staged interface from
genomic FASTA and existing GFF3/GFF/GTF annotations to correspondence evidence.

```bash
intraphy check --fasta genomes --gff annotations \
  --orthologs families --species-tree tree.nwk
intraphy build-case --fasta genomes --gff annotations \
  --orthologs families --species-tree tree.nwk \
  --output-dir prepared_case --threads 8
intraphy derive-tables --input-dir prepared_case \
  --output-dir correspondence_tables --threads 8
```

The `--orthologs` selection is needed when annotations contain multiple
genes; it does not itself establish orthology. Omit it for inputs containing
exactly one selected gene locus per species. `check` resolves identifiers,
tree tips, and coordinate bounds. `build-case` prepares loci and annotation
evidence. `derive-tables` calculates correspondence information. See
[genomic input contracts](inputs.md) for file naming, compressed input,
coordinate handling, and optional `extract-loci`.

Annotation exons and introns are evidence at native coordinates. Repeated
transcript records covering the same DNA must be deduplicated by physical
interval when defining copy units. Lack of an annotation alone does not
establish DNA absence. Ambiguous alignment or unsurveyed sequence is unknown.

## Qualified DNA-copy model

```bash
intraphy analyze --locus-model model.json \
  --species-tree tree.nwk --output-dir locus_result
```

`model.json` uses `intraphy.exon-locus-model/2`. It supplies the fixed rooted
tree's branch-length unit and provenance, a root distribution, tip DNA
evidence, rate groups, and a finite event catalogue. Every catalogue names
`state_model: binary` or `state_model: irreversible`. The current model
contains DNA material, copy positions, source-qualified `copy_duplication`,
and continuous `dna_deletion` opportunities. It has no feature, transcript
path, or mutually exclusive exon field. Correspondence tables are not
automatically converted into this catalogue.

`--parameter-mode fixed` evaluates supplied rates; the default fits groups
declared `fit`. `--threads` parallelizes independent units, while one
connected unit is serial. `--expected-edits` adds conditional branch counts
and can be more expensive. Results include `locus_fit.json`,
`locus_history.json`, and a run ownership record. State and event summaries
are conditional on the declared tree, root, mapping, catalogue, and
observation design.

`prepare-locus-evidence` applies a complete caller-supplied tip-by-material
DNA evidence array to a model template. Known presence or absence calls need
survey and evidence information; unknown is retained as unknown. Detection
sensitivity/specificity, if supplied, are fixed input probabilities.
`locus-statistics` provides same-data nested rate comparisons, profiles, and
an optional fixed-catalogue DNA bootstrap. See
[conditional statistics](locus_statistics.md).

Other command routes retained in the package serve separate historical or
experimental analyses. They do not supply the qualified DNA-copy catalogue or
extend the current locus model to intron-boundary evolution.
