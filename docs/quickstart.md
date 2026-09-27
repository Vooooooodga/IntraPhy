# Quick start

The biological inputs are genomic FASTA sequences, matching gene annotations,
and a supplied rooted species tree. For a whole-genome annotation, provide an
ortholog-family FASTA to select one gene locus per species. Existing gene
annotations give exon and intron coordinates, boundaries, CDS phase, and
strand. Multiple transcript records can describe the same DNA; they do not
multiply the number of physical copies.

```bash
intraphy check --fasta genomes --gff annotations \
  --orthologs family.fa --species-tree species_tree.nwk
intraphy build-case --fasta genomes --gff annotations \
  --orthologs family.fa --species-tree species_tree.nwk \
  --output-dir prepared_case
intraphy derive-tables --input-dir prepared_case \
  --output-dir correspondence_tables
```

`check` reports resolved locus selection and coordinate problems.
`build-case` prepares selected genomic loci and annotation evidence;
`derive-tables` calculates correspondence tables. See
[input requirements](inputs.md) and the
[CLI guide](cli.md).

The DNA phylogenetic model needs a separate, reviewed JSON catalogue of
homologous copy positions and qualified event opportunities. Mapping output
does not automatically establish duplication or deletion. An unresolved match
or missing annotation remains unknown.

For the bundled synthetic teaching example:

```bash
intraphy analyze --locus-model examples/exon_locus/dscam_like_model.json \
  --species-tree examples/exon_locus/species_tree.tsv \
  --output-dir locus_result
```

The model schema is `intraphy.exon-locus-model/2`. Its catalogue selects a
`binary` or `irreversible` DNA state model and supplies ordered material
tracts, copy positions, source-qualified duplications, continuous deletion
intervals, root distribution, tip calls, rate groups, and provenance. Rates
marked `fit` are estimated by maximum likelihood. Use `--parameter-mode fixed`
to evaluate supplied values. Optional `--expected-edits` adds conditional
branch event counts. `locus_fit.json` and `locus_history.json` report the fit
and DNA-copy history. Binary regain requires a present qualified source and
does not imply a reversible full process. The irreversible model permits
independent introduction on separate branches, so it is not a strict
stochastic Dollo process. See the
[model reference](exon_structure_model.md).
