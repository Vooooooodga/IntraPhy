> **V18 legacy documentation / retained reference.** The default V19 model and
> commands are described in [the README](../README.md) and
> [exon_structure_model.md](exon_structure_model.md). Do not interpret the old
> independent-layer analyses as exon configuration inference.

# Biological scope and analysis coverage

The formal analysis concerns local structure in supplied single-copy orthologs.
It can report evidence insufficient to decide. Supporting an input feature type
is not a guarantee that its homology, role or history can be recovered.

## Coverage is not conservation

Known state 0 is a supported negative observation. Known state 1 is a supported
positive. Unknown is neither. Seven present and three explicitly absent tips have
100% callable coverage. Seven present and three unknown tips have 70%. The
coverage fraction uses the full supplied tree panel, not the fraction with state 1.

Default `all` does not filter by this fraction. Optional `high-coverage` is a
sensitivity subset; 0.70 is an engineering setting without a calibrated biological
interpretation. Identical fractions can have different phylogenetic information
when called tips cluster in one clade. The output therefore records called species,
MRCA and representation of root subtrees.

Character selection leaves coordinates, full annotation and transcript paths
unchanged. A–unknown B–C does not become an A–C splice junction. A drawing split
never creates additional character replicates. Ascertainment and state-dependent
observation failure remain separate statistical issues.

## Twenty-eight scope cases

These cases are overlapping biological and observation situations, not 28 states
or a claim of 28 experimentally validated capabilities.

| No. | Situation | Supported interpretation and boundary |
|---|---|---|
| 01 | Conserved structure | Compare explicit local states; unknown regions are not conserved by default. |
| 02 | Unalignable intronic DNA with corresponding splice positions | Position characters may be informative without complete intron sequence homology. |
| 03 | Missing exon annotation | Retain DNA/projection candidates; do not automatically assert exon use. |
| 04 | Missing transcript annotation | Known-locus evidence is possible; no complete de novo gene prediction. |
| 05 | Supported DNA deletion | Code absence only with resolved location and adequate sequence evidence. |
| 06 | Assembly gaps or contig ends | Retain unknown and the observation reason. |
| 07 | Shifted splice boundary | Record exact positions; no generic one-step boundary-shift history model. |
| 08 | One-to-many complementary exon correspondence | Compare mapped blocks and elementary cutpoints; no compound-event count. |
| 09 | Many-to-one complementary exon correspondence | Same constraints in the reverse direction; reference does not imply ancestor. |
| 10 | Exonization | Compare DNA and annotation-conditional role separately; function is not inferred. |
| 11 | Intronization | Retained DNA and loss of exonic role remain distinct. |
| 12 | Exon skipping | Retain supplied paths; do not infer absent paths or usage frequency. |
| 13 | Intron retention | Describe annotated paths; no inclusion rate or NMD inference. |
| 14 | Alternative first or last exons | Retain path-specific endpoints and explicit completeness limitations. |
| 15 | CDS versus UTR identity | Preserve labels; no separate coding-identity evolutionary model. |
| 16 | Microexons | Require appropriate local evidence; short-sequence thresholds are not calibrated probabilities. |
| 17 | Internal repeats | Distinguish repeated reference coverage from complementary fragments; abstain on ambiguous instances. |
| 18 | Overlapping exons | Retain path ownership and actual intervals; drawing fragments are not replicates. |
| 19 | Antisense or nested genes | Preserve strand and ownership; position alone does not assign function. |
| 20 | Intronic noncoding RNA | Retain annotation context; no RNA-class-specific history model. |
| 21 | Transposable elements and exonization | Retain repeat labels and role observations; no TE origin or activity mechanism. |
| 22 | Unannotated noncoding sequence | Local DNA correspondence may be retained; no exhaustive element discovery. |
| 23 | Inversion or order rearrangement | Report incompatible order/ambiguous correspondence; no complete rearrangement history. |
| 24 | Frameshift or premature stop | Retain coding anomalies; no automatic pseudogene or adaptation claim. |
| 25 | Annotation micro-intron or exception | Preserve original exception metadata; platform-specific cases are not exhausted. |
| 26 | uORF, SECIS or RNA secondary structure | Retain supplied annotations; no translation-control or covariation inference. |
| 27 | Circular RNA or trans-splicing | Outside the formal linear single-locus model. |
| 28 | Multicopy families or gene fusion | Outside formal single-copy inference; historical experimental code remains separate. |

No molecular mechanism, selection coefficient, phenotypic cause or independent
mutation count is inferred. [Software tests](validation.md) do not establish
biological accuracy or statistical calibration.
