# ERG7 Candidate Genomic Case

This four-species case prepares the annotated Erg7 loci for an exploratory
genomic exon-structure analysis. [Zhu and Niu](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0061683) report an imprecise intron loss
shared by the *S. cryophilus* gene SPOG_00241 and *S. octosporus* gene
SOCG_04299, leaving 3 nt of intron sequence on the upstream exon, and propose
that it occurred in their common ancestor (Fig. 4D). This published observation
motivates the case; the four-way mapping and correspondence to current assembly
coordinates remain to be assessed.

The manifest selects RefSeq-annotated records whose product is lanosterol
synthase Erg7. The paper explicitly identifies the *S. cryophilus* and
*S. octosporus* pair as orthologous. The *S. pombe* and *S. japonicus* records
are selected from the same annotated product; their orthology to that pair has
not been independently established here by reciprocal comparison or a gene
tree. These are curated locus inputs, not orthology inferred by `build-case`.
The native workflow retains annotated transcript records and their partial
flags; no transcript usage is inferred. No truth file is supplied, and the
published event is not passed to preparation or analysis. Exact correspondence
between the published sequence tract and the current assembly annotations has
not been confirmed.

The rooted topology follows Zhu and Niu 2013, Figure 1:
`(japonicus,(pombe,(octosporus,cryophilus)))`. Every non-root branch in
`species_tree.tsv` has length 1 as an arbitrary scale. Fitted parameters and
posterior summaries are conditional on this tree and scale; branch lengths are
not measured divergence times or substitution lengths.

```bash
intraphy build-case --manifest manifest.tsv --species-tree species_tree.tsv \
  --output-dir prepared_case --threads 1 \
  --short-alignment-max-dp-cells unlimited

intraphy analyze --input-dir prepared_case --output-dir analysis \
  --model exon-structure-ctmc --parameter-mode fit \
  --branch-length-mode supplied --threads 1
```

The case is an exploratory comparison motivated by a published hypothesis. It
does not establish the event from current assemblies, validate orthology for
all four records, or support a general accuracy or novelty claim.

The unlimited setting removes the short-alignment DP-cell budget guard; it
does not truncate alignment states. Memory and execution time still bound the
preparation.

## Native analysis

The native fit qualified one unit with 277 complete states and converged.
Scry and Soct supplied informative partial observations; Spom and Sjap were
unknown across the full unit because sequence correspondence remained
unresolved. The focal branch `pombe_cry_octo` → `cry_octo` has exon-count
endpoint probabilities 0.2726608222 decrease, 0.3669871303 increase, and
0.3603520475 unchanged (`branch_exon_changes.tsv`). The analysis did not
recover a reliable ancestral branch fusion, and the published 3-nt remnant
remains unconfirmed in these assembly-matched inputs. This summarizes this
case's inputs and fit only; it is not an accuracy claim.
