# SPBC29A10.02 Genomic Case

The manifest supplies the spo5 gene group in four Schizosaccharomyces species.
All four numeric GeneIDs and their gene/mRNA/CDS records were read directly from
the assembly-matched local GFF3 files.

| Species | Manifest GeneID | GFF gene ID | Annotated transcript |
|---|---|---|---|
| S. pombe | 2540637 | gene-SPOM_SPBC29A10.02 | NM_001021957.3 |
| S. japonicus | 7052533 | gene-SJAG_04271 | XM_002175352.2 |
| S. octosporus | 25033683 | gene-SOCG_04721 | XM_013165313.1 |
| S. cryophilus | 25034398 | gene-SPOG_00066 | XM_013167573.1 |

Use `manifest.tsv` and `species_tree.tsv` as the two inputs to `build-case`.
The `gene_id` values resolve through the GFF `Dbxref=GeneID:` attributes. Raw
FASTA/GFF files stay at the flat database paths recorded in the manifest.
Transcript accessions identify provided annotations; no transcriptome data
are required. The manifest provides curated gene mappings for this case.

The rooted topology is `(japonicus,(pombe,(octosporus,cryophilus)))`, following
[Zhu and Niu 2013, Figure 1](https://doi.org/10.1371/journal.pone.0061683.g001).
Every non-root branch has unit length for qualitative parsimony. These values
carry no divergence-time or substitution-rate estimates.

Published expectations are exclusively in [truth_events.tsv](truth_events.tsv)
for external scoring. Keep that file outside the generated analysis input
directory and pass no event label, expected branch or ancestral state to
`build-case` or `run`. Assembly accessions and FASTA/GFF source paths are listed
in [manifest.tsv](manifest.tsv).
