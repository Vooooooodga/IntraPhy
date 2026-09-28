# Quick start

The biological inputs are genomic FASTA sequences, matching gene annotations,
and a supplied rooted species tree. For a whole-genome annotation, provide an
ortholog-family FASTA to select one gene locus per species. Existing gene
annotations give exon and intron coordinates, boundaries, CDS phase, and
strand. Multiple transcript records can describe the same DNA; they do not
multiply the number of physical copies.

```bash
intraphy analyze --fasta genomes --gff annotations \
  --orthologs family.fa --species-tree species_tree.nwk \
  --output-dir dna_result
```

`analyze` prepares selected genomic loci, derives correspondence evidence,
surveys qualified DNA positions, and fits their presence histories. It
describes DNA-position presence gains and losses, not source-directed
duplication mechanisms. Ambiguous placements, failed surveys, and unqualified
candidates remain unknown. DNA presence alone cannot infer intron-loss exon
fusion, boundary shifts, exonization with retained DNA, transcript usage, or
complete gene-structure evolution. See
[input requirements](inputs.md) and the
[CLI guide](cli.md).

For staged use, build a prepared case and run `prepare-genomic-evidence`; the
survey can be reused only with matching prepared inputs, tree panel, and survey
settings.

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
