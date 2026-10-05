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

If whole-unit correspondence fails, an exact physical exon interval can retain
local support when another complete, path-associated, nonconflicting annotated
locus has the same interval in genomic-alignment coordinates, and every
alignment column contains A/C/G/T in both sequences with sequence identity at
or above the configured threshold. The span must also meet the anchor-length
requirement.
The supported exon interior constrains compatible states; its complement and
any pre-existing unknown windows remain unknown, so boundary extensions and
fusion histories through those regions are not resolved. This local nucleotide
match does not establish gene-family orthology and does not require protein
projection evidence.

For NCBI GFF3 annotations, `start_range` and `end_range` locate uncertainty at
their respective genomic-coordinate boundaries, regardless of strand. These
ranges are mapped to local alignment windows and clipped to the current unit;
known sequence inside an exon remains observed. A propagated `partial=true`
flag is localized using explicit sibling or transcript ranges when available.
Generic partial flags without retained boundary locations remain conservatively
unknown over the affected unit, including older prepared inputs that lack the
original attributes.

Each catalogue has `observation_unit="genomic_exon_spans"`. The CLI requires
this marker when `--exon-configurations` reuses a prepared catalogue. A
transcript-configuration catalogue has a different observation unit and cannot
be relabelled for this model.
`--threads` defaults to 1 and sets per-family MAFFT/minimap2 preparation
threads (MAFFT uses `--thread N --threadit 0`; `threadtb` stays at its
default). During inference, multiple families run in parallel; for one family,
threads schedule local-unit likelihood and posterior work, capped by eligible
unit count. The pools are not nested, and the longest unit remains serial.
Keep OpenMP/BLAS threads at 1. `intraphy.log` records family-fit,
unit-likelihood, and unit-posterior stages; the markers aid diagnosis and do
not provide restart or recovery.

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

In the default automatic genomic catalogue, a qualified sequence-gap tract is
not declared as a DNA material history when every present-day locus places it
strictly inside one annotated exon and it has no exon-boundary or overlap
conflict. Such a tract creates no DNA states 0/1/2 or new catalogue cut; its
MSA columns, source coordinates, and existing unknown windows remain intact.
Cross-boundary tracts, tracts spanning a whole exon, and shared-deletion
opportunities remain eligible. This conditional structural model therefore
excludes ancestor-only boundary histories at unobserved internal-gap endpoints
and does not model the history of those internal indels. Legacy and
explicit-catalogue inputs retain their existing behavior.

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
marginalizes the complete endpoint distribution into physical exon-count pairs
and reports the probabilities of count increase, decrease, and no net
count change. The full finite count-pair domain is emitted without a probability
cutoff. A count can remain unchanged while exon geometry changes, including
boundary shifts or compensating changes. A count change does not identify a
split, fusion, or event count. These are
unit-specific endpoint summaries, not whole-gene extrapolations. The count
marginals repeat on each modal-pair row and must not be summed across rows. The
table also lists the joint modal observable parent/child configurations. Hidden material states
0 (unintroduced) and 2 (deleted) are both reported as absent. Joint
modes within the recorded absolute and relative tolerances of the maximum are
retained; no probability mass cutoff is applied. The
reported configurations are endpoint net differences and do not estimate the
number of transitions or a molecular mechanism. DNA-presence probabilities
cover only material tracts declared in the local catalogue, not all DNA or
coding sequence.

Likelihood fitting uses the sparse edit generator and exponential-action
pruning. By default, posterior inference also uses exponential actions and
contracts states that share the same observable exon geometry and declared
DNA-presence pattern. It retains state-order node marginals and reports
observable joint endpoint modes and exon-count pair marginals without storing
a full state-pair endpoint matrix. Count-pair probabilities marginalize all
states and material-origin scenarios under the same fitted CTMC; they describe
parent/child endpoint counts rather than transition histories. Mode rows use
absolute and relative tie tolerances of `1e-12`; these values apply to both
compact and dense branch summaries and are recorded in
`branch_exon_changes.tsv`. Compact rows also record the number of parent
groups examined by the modal search; dense rows record the full observable
group count. Complete state catalogues are evaluated under the declared root,
origin, and transition models.

`--expected-edits` requests transition-count expectations and uses the dense
posterior evaluator for that calculation. Its state-pair endpoint and
transition matrices require memory proportional to the square of the state
count. The default path reports the probability of at least one edit and the
observable endpoint summaries without those dense matrices.

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
