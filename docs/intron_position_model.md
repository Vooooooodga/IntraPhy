# Genomic intron-position histories

`intron-position-ctmc` estimates histories for annotated intron presence at
homologous coding positions. It is an independent character domain from
`dna-presence-ctmc`; each command fits one domain and emits its own fit and
history files.

## Observation unit and states

For each selected family, available coding paths are projected through a
protein multiple-sequence alignment. Adjacent coding nucleotides define
candidate positions, including positions within a codon. A coding nucleotide's
projected coordinate is `3 × MSA_column0 + codon_base_offset`. A position is retained when
an annotated physical intron in at least one species projects there. Repeated
transcript records are represented as aliases of physical intervals. Competing
physical boundaries or incompatible transcript paths make a candidate
ineligible; unresolved annotated introns remain catalogued with unknown states.

The tip states are:

- `1`: the annotated genomic intron interval exactly fills the gap between
  adjacent aligned coding bases, and both coding endpoint bases agree with the
  strand-oriented genomic sequence. The intron interior is not assessed for
  sequence conservation; internal `N` bases alone do not invalidate the call.
- `0`: the corresponding coding bases are consecutive in genomic sequence.
- `unknown`: the coding boundary cannot be projected or verified, or the
  physical evidence conflicts. An endpoint `N` or endpoint mismatch remains
  unknown; an absent declared FASTA resource or a backend failure raises an
  input/runtime error.

The protein MSA supplies conditional position homology for coding paths in the
selected family. It does not certify orthology or transcript use. The anchor
settings count amino-acid alignment columns on both flanks of the focal
boundary; unsupported local flanks remain unknown. Raw `analyze` prepares the
selected loci and coding paths, then runs this alignment route without deriving
the separate all-vs-all DNA correspondence tables. Staged evidence can be
prepared with `prepare-genomic-evidence --character-type intron-position`.

CDS phase in GFF is retained as an annotation field. `coding_right_base_offset`
identifies the right coding base's offset within its codon at the projected
boundary; it is separate from the GFF phase. At an internal junction in a
consistent CDS path, the downstream GFF CDS phase is
`(3 - coding_right_base_offset) % 3`. GFF3 defines phase relative to the 5′ end
of each CDS feature, with strand orientation taken into account; see the
[Sequence Ontology GFF3 specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md).

## Binary model

The latent states are intron absence and intron presence at the aligned coding
position; annotation and genome sequence provide the observations. The
continuous-time generator is

```text
              no intron     intron
no intron     -gain         gain
intron         loss        -loss
```

Gain and loss rates are shared among eligible sites within each family. Rates
are estimated by maximum likelihood, or both can be supplied with
`--gain-rate` and `--loss-rate`. The default root distribution is stationary;
`--root-frequency fixed --root-presence P` supplies a fixed root state
probability. All estimates use the branch lengths and units in the supplied
rooted species tree. The likelihood conditions on the observed-at-least-one
discovery rule and fixed observation mask. Unknown tips are integrated as
unobserved states.

Outputs are `intron_fit.json`, `intron_history.json`, and `run_result.json`.
The preparation artifacts are `intron_families.tsv`, `intron_positions.tsv`,
`intron_observations.tsv`, and `intron_alignment.tsv`. The family summary keeps
zero-candidate roster/count/status information without adding a phylogenetic
site. Catalogued introns
without an analyzable position are retained as unresolved evidence; a family
with no annotated introns has no phylogenetic site and does not receive a
fabricated history.

The model does not estimate intron-length evolution, ancestral intron length,
transcript usage, rate heterogeneity, or genome-wide intron density. Inferred
state changes describe changes in intron absence/presence at an aligned coding
position, with annotation and genomic sequence supplying the tip evidence; they
are not counts of molecular lesions or splice outcomes. Linked sites may make the
product likelihood composite. Annotation, orthology selection, MSA placement,
genomic resource agreement, and sampling remain conditions on interpretation.

## Related comparative methods

Csűrös et al. modeled binary intron positions with ambiguous observations and
ascertainment correction; their framework provides context for position-level
histories, while its assumptions do not establish this mapping process.
[Csűrös et al. 2011](https://doi.org/10.1371/journal.pcbi.1002150).
Glick et al. compared species-level summaries of exon/intron structure across
eukaryotes using phylogenetic comparative analyses; those summaries estimate a
different unit than aligned intron sites.
[Glick et al. 2024](https://doi.org/10.1093/molbev/msae248). ExonEvo compares
homologous exon/domain organization in angiosperms and is a related structural
baseline, not validation of this model.
[Zhang et al. 2026](https://doi.org/10.1038/s41467-025-66816-3).
