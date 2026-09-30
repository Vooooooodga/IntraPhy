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
to these four species. The supplied non-root branch lengths are all one and
provide an arbitrary common scale. No calibrated divergence-time or
substitution lengths are supplied in this tree. Native CTMC rates and posterior
histories are conditional on this scale, the supplied topology, the qualified
evidence, and this four-species sample. Directional support is limited by that
sampling and should be read from the conditional branch posterior.

Published expectations and the distinction between structural correspondence
and branch direction are in [truth_events.tsv](truth_events.tsv), for external
scoring only. Keep that file outside the generated analysis input directory.
See the [genomic exon-span model](../../../docs/genomic_exon_model.md) for
model scope and coordinate interpretation.
