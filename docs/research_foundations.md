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

Homology defines the unit to which a structural state refers. A DNA-presence
observation concerns DNA at a proposed homologous interval; an intron-position
observation concerns a corresponding splice boundary; and intron-length
analysis concerns the size of sequence between boundaries. These estimands
must remain separate. A match of flanking exons alone does not establish
homology of all intervening bases.

IntraPhy's genomic observation layer records evidence and unresolved states for
candidate DNA intervals and structural features. Its binary DNA-presence CTMC
uses state 0 for supported absence and state 1 for supported presence at a
homologous position. A 0-to-1 change is a structural material gain at that
position; the model does not identify duplication, transfer, or another source
mechanism. Intron boundaries and lengths are genomic annotation descriptors
here; the DNA-presence history does not infer their evolution.

The DNA-presence rate analysis conditions on the supplied species tree and
branch lengths. It fits family-shared gain/loss rates over eligible sites, with
an observed-presence ascertainment correction for discovered sites and
unknown tip states treated as compatible with either state. This conditions on
the observation mask; it does not estimate missingness or correct all
nonrandom annotation, transcript-selection, mapping, or discovery biases.
Eligibility and evidence review therefore remain central. A site-state change
count, including an expected CTMC transition count, is not a count of physical
molecular lesions. Linked sites may also make the product likelihood a
composite rather than a fully independent-site likelihood.

The species-tree analysis of single-copy families and a copy/exon-tree analysis
of duplicated exon families concern different evolutionary units. IntraPhy's
supplied species tree is fixed for the former; the latter requires an explicit
copy-family history and should not be conflated with a species-tree state
reconstruction. Likewise, orthology pipelines, exon-family trees, and
annotation-aware homology qualification in cited studies provide methodological
context, not automatically equivalent inputs or results.

## Implemented and prospective scope

Implemented analyses report evidence-qualified structural observations and
conditional histories under the selected tree, state model and ascertainment
rule. DNA-presence CTMC rate estimates are conditional on included eligible
sites, fixed branch lengths and the stated root assumption; rate diagnostics
and data limitations must accompany interpretation. Unresolved and
insufficient-information families do not support a fabricated finite estimate
or ancestral history. Site-level changes are not summed into inferred mutation
events.

The advanced `intraphy.exon-locus-model/2` route already evaluates a supplied
catalogue of source-qualified copy duplications and continuous shared deletion
opportunities. Fully automatic genealogy inference, or a unified model of
boundary, copy and length evolution, remains prospective. So do full
reproductions of ExonEvo or MALIN and integrated models of annotation error.
The cited papers motivate careful homology definitions, explicit observation
rules, sampling and rate diagnostics. They do not by themselves establish
empirical novelty or comparative performance for IntraPhy.
