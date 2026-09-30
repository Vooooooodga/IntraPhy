# Genomic exon-span model

## Observation unit

The `exon-structure-ctmc` model follows annotated physical exon spans on a
supplied rooted species tree. Whole-locus genomic alignment and nucleotide
searches support correspondence using sequence, coordinates, strand, and local
order. Existing staged protein-projection evidence, when present, can provide
corroboration. CDS phase is retained as annotation metadata; it does not score
or qualify correspondence. Full exon spans are retained, including terminal
and UTR exons.
An exon start or end is a geometric boundary on the transcriptional axis; a
terminal boundary does not automatically have splice donor or acceptor
function.

Multiple transcripts can identify the same physical interval. Those aliases
are deduplicated for the observation unit. Distinct overlapping annotated spans
remain separate evidence records. If their local structures conflict, the
affected observation is unresolved. A missing annotation, ambiguous mapping,
or insufficient sequence support is unknown. Annotation records are evidence
and can contain errors. The method does not infer RNA abundance, isoform use,
splicing, or an RNA phenotype.

Each catalogue has `observation_unit="genomic_exon_spans"`. The CLI requires
this marker when `--exon-configurations` reuses a prepared catalogue. A
transcript-configuration catalogue has a different observation unit and cannot
be relabelled for this model.
`--threads` defaults to 1 and sets per-family MAFFT and minimap2 thread counts
during genomic alignment preparation, then controls parallel family fits.
MAFFT receives `--thread N --threadit 0`; `threadtb` remains at its default.

`analyze --state-space-only` runs genomic exon catalogue preparation, complete
state enumeration, and annotation-view observation eligibility assessment. It
does not fit family rates or evaluate CTMC probabilities. Qualified catalogues
are enumerated without a state-count limit, and unqualified catalogues are
recorded as not enumerated. Its result manifest uses status `state_space_only`;
the flushed `state_space_diagnostics.jsonl` records each unit before
enumeration and after eligibility assessment. An interrupted run retains the
latest completed phase. This output is a diagnostic record rather than a
completed inference result.

## State and transition process

A local unit represents an ordered exon-span geometry together with the
presence state of declared local DNA material tracts. Its finite edit graph
contains eligible exon split, adjacent exon fusion, exon-boundary displacement,
exonization, exon inactivation, DNA insertion, and DNA deletion transitions.
A DNA deletion can remove material shared by multiple exons in one transition;
its resulting geometry is recorded once. Transition opportunities are those
generated from the declared sequence and annotation catalogue. Omitted
opportunities have zero rate.

For the genomic exon-span observation unit, alignment gaps connect exon spans
into one local event unit and enter the declared material catalogue only when
they pass the stated flanking-anchor support test. Unsupported gaps do not
imply material absence: their intervals are retained as species-specific
unknown regions, and affected observed tips are marked partial. Distinct
overlapping qualified material tracts remain unresolved.

Each identified variable material tract has states 0 (not introduced), 1
(present), and 2 (deleted). Its origin is assigned once, either at the root or
on one branch. State 2 cannot return to state 1 in that positional catalogue.
Tip absence is compatible with states 0 and 2, while tip presence requires
state 1. DNA presence alone does not require an exon annotation at that span.
This single-origin material process is distinct from the advanced
`exon-locus-ctmc` model, whose declared copy opportunities can reintroduce a
position independently on multiple branches.

The root geometry is uniform over valid geometries conditional on its
material-origin state. Material-origin scenarios have a separate prior: the
root opportunity and each possible first-introduction branch have weight one.
This prior does not estimate root-state frequencies. On branches, elementary
edit opportunities have unit weight subject to the edit graph's normalization
for multiple destinations.

By default, one nonnegative scalar rate is fitted per family and shared across
all elementary edit kinds and linked local units. Equal rates form a baseline
for the declared edit graph; the data and design do not establish equal
biological rates. A fixed rate file can specify rates explicitly and is used
with `--parameter-mode fixed`. Linked units contribute a composite likelihood.
The likelihood is conditional on the supplied tree, its branch-length unit,
the finite catalogue, and the stated independence structure. Annotation
discovery is not corrected for ascertainment. The model does not provide a
significance test for exon evolution.
Before accepting a fitted rate, the fitter probes once at twice the larger of
the selected dimensionless rate and its largest start; a non-finite or
numerically indistinguishable/higher probe withholds the estimate as
`upper_tail_unresolved`. This finite probe is a safeguard, not an asymptotic or
global-optimum test; accepted estimates retain the existing multi-start fit
scope.

`calibrate-exons` provides a small fixed-catalogue conditional check for
endpoint-change probabilities. It simulates histories on two built-in
independent catalogues and compares inference at the generating rate with
inference after fitting the shared rate from tips, state spaces, and tree.
Scores report scored, failed, and nonidentified replicate denominators for each
inference arm. This calibration does not cover FASTA/GFF catalogue discovery.

The finite geometry catalogue is derived from observed boundaries and declared
material-tract cuts; the model does not enumerate every possible exon boundary
in sequence space. There is no default state-count cap. `--max-states` and
`--max-origin-scenarios` are optional resource-abort limits. Reaching a supplied
limit stops inference without truncating or renormalizing the catalogue or
state space.
`--expected-edits` reports conditional model transition counts, which are not
counts of molecular lesions.

Every run writes `branch_exon_changes.tsv` from the existing CTMC joint
endpoint probabilities. It summarizes the probability of exon-span structure
change and declared material-tract DNA-presence change on each branch, and
lists the joint modal observable parent/child configurations. Hidden material
states 0 (unintroduced) and 2 (deleted) are both reported as absent. Every
exactly tied joint mode is retained; no probability cutoff is applied. The
reported configurations are endpoint net differences and do not estimate the
number of transitions or a molecular mechanism. DNA-presence probabilities
cover only material tracts declared in the local catalogue, not all DNA or
coding sequence.

`change_classification` describes endpoint geometry. `net_split` and
`net_fusion` require matching outer span boundaries for a one-to-many or
many-to-one overlap group. `boundary_change`, `span_gain`, and `span_loss`
describe the corresponding span differences. When a declared DNA tract
changes with a split/fusion group, the classification is `complex_change` and
`dna_coupled`. These labels summarize endpoint differences and do not infer
the path or number of evolutionary changes. Tied modal rows are alternative
joint endpoint explanations and must not be added as event counts.

Exon spans use local alignment coordinates, 0-based and half-open, with
`alignment_offset` locating the local unit on its shared alignment axis. For
observed native exon coordinates, join the family/unit and local span through
`exon_coordinates.tsv` and `exon_correspondence.tsv`.

## Interpretation

Inferred boundary changes describe shifts in exon-span geometry on a common
transcriptional axis. The registry retains the labels `acceptor_shift` and
`donor_shift` for the two directional boundary edit classes; for terminal
exons these remain exon-start or exon-end shifts, without inferred splice-site
function. Split, fusion, insertion, and deletion histories depend on the
declared candidates and material tracts. A model transition does not establish
a molecular mechanism, transcript usage, or selection.

Results remain conditional on gene-family selection, annotation quality,
sequence correspondence, the finite catalogue, the root and branch priors,
and the supplied tree. Local unknown states are retained as unknown; they are
not converted to absence.
