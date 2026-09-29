# Research foundations and scope

This page situates IntraPhy's structural observations and phylogenetic analyses
among relevant comparative studies. It summarizes the main texts of the cited
papers; supplementary materials were not comprehensively audited. The papers
motivate distinct questions and are not interchangeable benchmarks for one
another.

## What the studies establish

Glick et al. analyzed genome-wide traits across 590 eukaryotic genomes, selecting
a longest-CDS representative for gene-level summaries and using phylogenetic
comparative methods including RRphylo and PGLS. This is a broad comparative
framework for genome traits. Such summaries do not identify evolutionary
histories at aligned, homologous structural positions.
[Glick et al. 2024, *Molecular Biology and Evolution* 41:msae248](https://doi.org/10.1093/molbev/msae248).

Csűrös et al. modeled intron positions in orthologous genes as binary characters,
with presence, absence and ambiguous observations. Their likelihood treatment
accounts for ascertainment of observed intron sites, and the analyses consider
variation in gain/loss processes across lineages and sites. These are
position-specific intron histories; they do not model exon-sequence homology or
the full annotation process.
[Csűrös et al. 2011, *PLoS Computational Biology* 7:e1002150](https://doi.org/10.1371/journal.pcbi.1002150).

Ma et al. studied intron gain and loss in nematodes using orthology and flanking
sequence evidence to qualify candidate positions, with broad and dense sampling
for different phylogenetic depths. Their MALIN analyses address rate
heterogeneity, and reported candidate gains were evaluated with additional
sequence and transcript evidence. The study illustrates how orthologue
qualification and sampling affect inferred histories; it does not remove
annotation and homology-detection uncertainty from all datasets.
[Ma et al. 2022, *Biology Direct*](https://doi.org/10.1186/s13062-022-00328-8).

Jakt et al. treated intron length as a separate evolutionary quantity from
intron presence. They compared orthologous teleost introns and reconstructed
length change using discretized size states and Sankoff parsimony. Presence at
a homologous position can persist while its length evolves, so length
trajectories should not be read as gain/loss histories.
[Jakt et al. 2022, *BMC Genomics* 23:628](https://doi.org/10.1186/s12864-022-08760-w).

Lee et al. examined arthropod *Dscam* exon clusters through exon-family
relationships, cluster order and ancestral structural reconstruction. Their
interpretation of exon expansion draws on those complementary signals and
proposes staggered homologous recombination as a plausible process. A cluster
history is a copy/exon-family problem, distinct from binary presence at a
single-copy homologous position on a species tree.
[Lee et al. 2010, *RNA* 16:91–105](https://doi.org/10.1261/rna.1812710).

Zhang et al.'s ExonEvo compares exon organization across angiosperms using
annotated transcripts, hierarchical orthogroups, exon similarities and
phylogenetic placement on species relationships. The main text describes
selecting a representative transcript per gene and using conserved exon and
flanking-exon evidence to define comparisons and map inferred changes to
ancestral branches. Representative-transcript and annotation choices bound
which alternative exons can be observed; a mapped change remains conditional
on those inputs and homology rules.
[Zhang et al. 2026, *Nature Communications* 17:106](https://doi.org/10.1038/s41467-025-66816-3).

Dong et al. combined comparative *Dscam1* cluster analysis with CRISPR-mediated
deletions and measured fitness and immune phenotypes in *Drosophila*. The
functional experiments address consequences of cluster variation in that
system; they do not test whether Hymenoptera sociality caused Dscam evolution.
[Dong et al. 2025, *PLOS Biology* 23:e3003383](https://doi.org/10.1371/journal.pbio.3003383).

Yu et al. combined RNA-seq evidence from ten plants with conserved microexon
clusters and searched genomic sequence across a much broader panel using
coding-context and splice-boundary features. This design highlights a key
distinction: a genomic sequence match can support homology without confirming
exon use in a transcript. Training-set composition, intron-length bounds and
RNA evidence constrain the scope of such annotation methods.
[Yu et al. 2022, *Nature Communications* 13:820](https://doi.org/10.1038/s41467-022-28449-8).

## Implications for IntraPhy

The default `exon-structure-ctmc` observation unit is a mapped physical genomic
exon span. Homology is supported by genomic sequence and local geometry using
the selected gene annotations. Repeated transcript aliases for one interval
are deduplicated, while distinct overlapping annotations remain distinct and
can have unresolved local structure. A missing annotation does not establish
absence. The model follows complete exon spans, including terminal and UTR
exons. Exon-start and exon-end displacement describe boundaries on a common
transcriptional axis; terminal ends do not automatically carry splice donor or
acceptor function. Genomic sequence, coordinates, strand, and local order
support correspondence. CDS phase is retained as annotation metadata; it does
not score or qualify correspondence.

ExOrthist uses genomic coordinates together with exon sequence and flanking
exon context to support exon-orthology inference. This provides methodological
context for IntraPhy's genome and annotation based correspondence strategy; it
does not validate the exon-structure CTMC, its priors, or its rate assumptions.
[Márquez et al. 2021, *Genome Biology* 22:239](https://link.springer.com/article/10.1186/s13059-021-02441-9).

The exon-structure CTMC conditions on the supplied rooted tree, its branch
length scale, the mapped finite catalogue, and local-unit dependence
assumptions. It uses a uniform distribution over valid root exon geometries
conditional on material origin. Separately, material-origin opportunities at
the root and on each branch have unit weight. Its
default rate fit is one nonnegative scalar shared across elementary edit
opportunities within a family, using the composite likelihood over linked local
units. These conventions define a baseline for the declared edit graph; they
do not establish equal biological event rates or provide a significance test.
Annotation discovery is not corrected for ascertainment. Conditional expected
transition counts are model changes, not counts of physical lesions.

IntraPhy also retains explicit `dna-presence-ctmc` and
`intron-position-ctmc` analyses. They represent separate estimands: binary DNA
presence at qualified homologous intervals, and annotated intron presence at
aligned coding positions. They are selected explicitly and are not combined
with exon-span histories in a joint likelihood. The intron model's observation
rules and limitations are described in the [intron-position method](intron_position_model.md).

The species-tree analysis of single-copy families and a copy/exon-tree analysis
of duplicated exon families concern different evolutionary units. IntraPhy's
supplied species tree is fixed for the former; the latter requires an explicit
copy-family history and should not be conflated with a species-tree state
reconstruction. Likewise, orthology pipelines, exon-family trees, and
annotation-aware homology qualification in cited studies provide methodological
context, not automatically equivalent inputs or results.

## Implemented and prospective scope

Implemented analyses report annotation-supported exon-span observations and
conditional histories under the selected tree and finite edit graph. Linked
units contribute a composite likelihood. DNA-presence and intron-position
analyses retain their own observation rules when explicitly selected.
Unresolved and insufficient-information observations remain unknown; they do
not support a fabricated finite estimate or ancestral history.

The advanced `intraphy.exon-locus-model/2` route evaluates a supplied catalogue
of source-qualified copy duplications and continuous shared deletion
opportunities. Fully automatic genealogy inference, unified copy and exon
boundary evolution, integrated models of annotation error, and full
reproductions of ExonEvo or MALIN remain prospective.
The cited papers motivate careful homology definitions, explicit observation
rules, sampling and rate diagnostics. They do not by themselves establish
empirical novelty or comparative performance for IntraPhy.
