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
  --output-dir exon_result
```

`analyze` maps physical exon spans from the selected genomic annotations and
fits the `exon-structure-ctmc` model. Repeated transcript aliases for the same
native interval are deduplicated. Conflicting overlapping annotations,
missing annotation, and ambiguous mapping can leave local observations
unknown. The model includes terminal and UTR exons when annotated and does not
infer RNA transcript use or splicing. See [input requirements](inputs.md), the
[genomic exon-span model](genomic_exon_model.md), and the [CLI guide](cli.md).

The explicit `--model intron-position-ctmc` route remains available as a
separate analysis; see the [CLI guide](cli.md) and
[intron-position method](intron_position_model.md).

`--model dna-presence-ctmc` selects the separate binary DNA-position model.
Its staged evidence route is described in the [CLI guide](cli.md).

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
