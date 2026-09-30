# Hdac3 Genomic Case

The manifest supplies HDAC3 in D. ananassae, D. melanogaster, D. simulans and
D. yakuba. It combines the newly downloaded D. ananassae assembly with the
three existing Drosophila genome/GFF resources. All four numeric GeneIDs and
their gene/mRNA/CDS records were read from these exact local GFF files.

| Species | Manifest GeneID | GFF gene ID | Annotated transcript |
|---|---|---|---|
| D. ananassae | 6500490 | gene-LOC6500490 | XM_001953300.4 |
| D. melanogaster | 44446 | gene-Dmel_CG2128 | NM_143721.4 |
| D. simulans | 6726814 | gene-LOC6726814 | XM_002102185.4 |
| D. yakuba | 6535413 | gene-LOC6535413 | XM_002096026.4 |

Use `manifest.tsv` and `species_tree.tsv` as the two inputs to `build-case`.
The `gene_id` values resolve through the GFF `Dbxref=GeneID:` attributes.
The three comparator protein-file cells are empty: their genome/GFF resources
supply the required sequence and annotation inputs. Original annotation
evidence fields are preserved; no RNA sequencing data enter this case.
These are curated gene mappings supplied to the method; an OrthoFinder analysis
has not been performed for this case in this preparation.

The GFFs name the three non-melanogaster genes `LOC6500490`, `LOC6726814` and
`LOC6535413`; their mRNA/CDS product is histone deacetylase 3. D. melanogaster
also carries the original gene-level `exception=dicistronic gene`. Retain that
metadata and use the HDAC3 transcript/CDS hierarchy specified above.

The topology `(ananassae,(yakuba,(melanogaster,simulans)))` is the restriction
of the [Drosophila 12 Genomes Consortium 2007 tree](https://doi.org/10.1038/nature06341)
to these four species. In `species_tree.tsv`, the supplied non-root branch
lengths are all one and provide an arbitrary common scale. No calibrated divergence-time or
substitution lengths are supplied in this tree. Native CTMC rates and posterior
histories are conditional on this scale, the supplied topology, the qualified
evidence, and this four-species sample. Directional support is limited by that
sampling and should be read from the conditional branch posterior.

## Fourfold-site branch-length comparison

`species_tree_fourfold.tsv` is an alternate four-taxon input for a separate
conditional-ML comparison. It restricts the published 12-species fourfold-site
tree to the four species in this manifest, roots the restriction at their MRCA,
and sums branch lengths across collapsed unary paths. The lengths are expected
substitutions per site estimated from fourfold-degenerate sites in the
historical comparative-genomics dataset. They are not divergence times, rates
for structural evolution, or estimates re-fit on the assemblies used here.
The native analysis keeps the same model and re-fits its parameters; improved
fit or stronger branch support is not guaranteed.

Source: [Stark et al. 2007 supplementary tree](https://compbio.mit.edu/flies/stree/flies.fourfold.tree), linked from the [author-hosted supplement](https://compbio.mit.edu/flies/). In this four-taxon restriction, the D. yakuba length is the unary-path sum 0.038942 + 0.087380 = 0.126322, and the D. simulans length is 0.020119 + 0.022856 = 0.042975; all other retained branch lengths are unchanged. These path sums are documented so the input can be reviewed without a conversion step.

Published expectations and the distinction between structural correspondence
and branch direction are in [truth_events.tsv](truth_events.tsv), for external
scoring only. Keep that file outside the generated analysis input directory.
See the [genomic exon-span model](../../../docs/genomic_exon_model.md) for
model scope and coordinate interpretation.
