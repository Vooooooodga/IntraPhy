# Planned DNA copy model ablations

## Question and status

This protocol asks which mapping evidence and DNA-copy event assumptions
affect reconstruction under the advanced `exon-locus-ctmc` process, using
supplied `intraphy.exon-locus-model/2` catalogues. It is a planned study. No
simulation, empirical benchmark, calibration result, performance threshold,
or resource estimate is reported here. The contrasts assess DNA-copy histories
for these declared catalogue inputs; they do not test the default
`exon-structure-ctmc` or its genomic-exon mapping and preparation pipeline. The
default model follows physical genomic exon spans and models boundary
displacement, exon split and fusion, exonization, inactivation, and DNA
insertion and deletion. Its scope and assumptions are described in the
[genomic exon-span model](genomic_exon_model.md).

Contrast 1 varies how sequence and genomic context are used to map homologous
DNA-copy positions and qualify copy opportunities. Its context/order arm
examines local intron/exon context, CDS phase, collinearity, and copy order as
planned mapping evidence; these are experiment-specific assumptions. CDS
phase remains annotation metadata and does not score or qualify correspondence
in the default exon-span model. Multiple transcript records describing one
physical DNA interval are deduplicated. Missing annotation and unresolved
sequence correspondence are unknown; deletion requires surveyed sequence
evidence. These are correctness requirements in every arm, not ablated
options. The likelihood in this experiment infers DNA-copy presence and
events. Annotation boundary and intron-position differences are mapping
outputs, with no separate phylogenetic event process in this experiment.

## Prespecified contrasts

Exactly three main contrasts are planned:

| Contrast | Common target | Changed assumption |
| --- | --- | --- |
| 1. Sequence-only versus sequence plus context/order mapping | The same true homologous genomic intervals and copy positions | Candidate correspondence scores and qualification use sequence alone or additionally use local intron/exon context, CDS phase, collinearity and order. |
| 2. Independent DNA changes versus shared-interval deletion | The same observed DNA tracts, tree, and root conditions | Deletion opportunities act on tracts separately or as one continuous event across all intersecting tracts. |
| 3. Binary versus irreversible state assumption | The same labelled positions, observations, event opportunities, tree, and root observables | A deleted position can later regain presence through a qualified source in `binary`; state 2 cannot regain in `irreversible`. |

Contrasts 1 and 2 require prespecified mapping/catalogue construction and
simulation configurations; the CLI does not automatically generate an ablation
experiment. Contrast 3 is an explicit `state_model` selection in the current
DNA locus schema. The independent-change arm of contrast 2 must use eligible
single-tract deletion intervals in the same declared DNA opportunity
framework. When coupling changes, calibrate its deletion rate so the
predefined per-tract marginal deletion intensity is comparable over the
relevant root and occupancy conditions. Record any residual difference
caused by state-dependent eligibility.

## Designs and data partitions

Prespecify the homologous DNA units, rooted tree and branch-length unit,
observable root conditions, locus discovery rule, surveyed regions, unknown
calls, and detection process before generating data. Keep these identical
within each paired contrast. Represent root observables consistently across
binary and irreversible states; if both latent 0 and 2 are possible at the
root in an irreversible design, declare their mixture rather than silently
changing the root distribution. Use separate simulation seeds, loci, or
families for parameter tuning and final evaluation. Fix the primary metrics
and analysis code before examining the final evaluation set.

Simulated truth must include both independent and compound deletion processes,
and both recurrent and irreversible position histories. Vary copy-number
ambiguity, tract spacing, branch-length scale, detection coverage, and
sequence similarity within a prespecified, tractable design. Include complete
surveys, fixed unknown masks, and unknown masks generated independently of
latent presence as distinct designs. Add a separate, explicitly labelled
missingness mechanism if real annotation incompleteness may depend on the
state. Use an independently specified generator or independent
implementation to avoid evaluating only internal consistency. Keep catalogue
discovery rules fixed and document which simulated units are unobservable or
excluded by that rule.

Evaluate mapping against true correspondence first. Evaluate phylogenetic
inference initially with true correspondence supplied to all methods, so
mapping error does not confound event-model accuracy. Then propagate each
mapping arm into the same downstream DNA model under a separately reported
end-to-end analysis. The two stages answer different questions and should not
be pooled.

An independent real-data set of annotated genomes will assess agreement with
held-out sequence context, conserved position/order, and manually reviewable
genomic alignments. Keep these families separate from method development.
Species-tree source, assembly/annotation versions, locus selection, callable
regions, and coordinate transforms must be recorded. Agreement with annotation
is evidence about correspondence, not a known evolutionary truth. Mhc or
Dscam observations chosen for biological interest cannot alone calibrate
false-positive rates, probability calibration, or P values.

## Endpoints and interpretation

Primary endpoints are correspondence precision/recall with an unresolved
fraction for contrast 1, and event-history accuracy and false-positive rate
for contrasts 2 and 3. Report event localization error, event-count error,
ancestral-state probability calibration, and uncertainty intervals across
independent loci as secondary endpoints. Stratify by callability, copy
ambiguity, tract overlap, branch, and generative truth class. Report
identifiability diagnostics, failed/unresolved fits, CPU time, peak RAM, and
state-space size for every arm.

Compare inference methods on matched observed DNA units with the same tree,
survey design, and observable root conditions. When coding or discovery
changes the data space, do not use the larger raw likelihood as a winner.
Use truth-based error metrics and, where meaningful, proper predictive scores
on a common held-out observation space. Separate rate estimation from
posterior event reconstruction. A low error under one simulation family
supports only that family and its declared discovery/detection assumptions.
No full factorial design or unreported significance threshold is implied.
