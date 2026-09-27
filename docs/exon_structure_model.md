# Exon-copy and splice-feature likelihood model

## Biological question and scope

The model estimates event rates and reconstructs ancestral states for supplied
exon-copy material and splice-feature availability along a rooted species tree. It is conditional on
the topology, branch lengths, root distribution, homology catalogue, tip evidence
and finite event-opportunity catalogue. It does not estimate transcript usage,
RNA expression, protein abundance, selection, unrestricted copy birth or copy
genealogy.

## Model objects

Each independent unit has a `LocusCatalogue` containing:

- `MaterialTract`: one nonoverlapping segment on a shared ordered coordinate axis.
- `CopySlot`: one labelled copy position, its material IDs, explicit homology and
  collinearity links, orientation and evidence.
- `SpliceFeature`: exon availability or a directed splice link. Exons require
  material and a copy; splice links require two exon features and their endpoint
  order. Prerequisites form a directed acyclic graph.
- `EventOpportunity`: one event opportunity and its weighted alternative
  outcomes, including material/feature preconditions, source and target copy,
  and source-to-target feature mappings where required.
- `mxe_groups`: feature sets whose alternatives constrain compatible paths.

The latent state is the material lifecycle vector and the set of active
availability features. Material states are 0 (not introduced), 1 (present), and
2 (deleted). Feature availability is separate from DNA material. A sequence
absence at a tip only establishes binary DNA absence; it does not identify
whether the tract was never introduced or was deleted.

Write a state as $x=(m,a)$, where $m$ is the vector of material lifecycle
values and $a$ is the set of available features. Along one lineage, deleted
material (state 2) cannot be reintroduced: supplied duplication opportunities
can change 0 to 1, and deletion opportunities can change 1 to 2.

The process is compiled from the explicitly supplied root-state support. All
states reachable under the declared opportunities are retained. There is no
default state-count cap. The compiled space can grow rapidly with the number of
copies, material segments, features and opportunities; full closure may require
substantial memory and time.

Material coordinates use one shared axis, increasing in transcriptional
direction across the modeled units of a family. Inputs in different species'
native chromosome coordinates are not automatically converted to that axis.
Orientation and collinearity are recorded evidence metadata; inversion events
are outside this process.

## Event process

The current typed events are:

1. **Copy duplication:** a supplied source copy can introduce a distinct target
   slot only when their homology and evidence are declared, source material is
   present, target material is unintroduced, and a feature map is supplied.
2. **DNA deletion:** one continuous declared interval changes its present
   material to deleted in one event. Dependent exon and splice features follow
   the declared prerequisite cleanup.
3. **Splice change:** context-specific feature availability changes while DNA
   copy material remains present. Availability does not imply use.

The finite opportunity catalogue defines the event process. Omitted events are
unavailable in the fitted model. Alternative outcomes for one opportunity have
weights summing to one. An opportunity applies only when its shared preconditions
hold; partial applicability fails model compilation rather than redistributing
weights. A no-change outcome keeps its probability mass and is not counted as a
jump.

The model does not invent homology from transcript similarity, infer root states
from the union of tip observations, or infer unrestricted duplication history.
The Dscam-like example uses explicitly labelled synthetic alternatives. Primary
studies report arthropod Dscam alternative splicing and copy diversification
[Lee et al. 2010](https://doi.org/10.1261/rna.1812710); those results motivate a
use case and do not establish any synthetic history.

## Likelihood and rates

For each unit, a continuous-time Markov chain is evaluated on the supplied tree
with sparse exponential actions and pruning. The root probability vector is
explicit. Global rate-group labels tie opportunity intensities across units.
Values marked `fit` provide positive starting values for nonnegative maximum
likelihood; values marked `fixed` remain fixed. `--parameter-mode fixed` treats
all declared values as fixed. No arbitrary upper rate bound is imposed.

For distinct states $i\ne j$, the generator entry is
$q_{ij}=\sum_e r_{g(e)}w_e$ over eligible declared events from $i$ to $j;
diagonal entries are the negative outgoing rate. Here $r_{g(e)}$ is the
rate-group intensity and $w_e$ the supplied outcome weight. Along a branch of
length $t$, $P(t)=\exp(Qt)$. For unit $u$ on tree nodes $V$ and edges $B$, the
likelihood is

$$L_u=\sum_{x_V}\pi(x_{\mathrm{root}})
\prod_{(p,c)\in B}P_{x_p,x_c}(t_{pc})
\prod_{v\in\mathrm{tips}}E_v(x_v),$$

where $\pi$ is the supplied root distribution and $E_v$ is the tip emission
likelihood. Fitting maximizes $\sum_u\log L_u$ over declared free rate
groups. Branch lengths and their units are fixed inputs.

The branch-length unit is required in the model, as is provenance for the tree
and branch-length basis. Rates have units of events per declared opportunity
unit per tree branch-length unit. Tree scale and rate scale are confounded unless
branch lengths and units are fixed as supplied. Fit diagnostics report optimizer
status, boundary estimates and a local curvature/rank summary. These describe
the fitted objective near the returned solution; they do not certify a global
maximum or broad identifiability. Single-gene fits are allowed, with limited
information exposed in those diagnostics.

With `--expected-edits`, branch event counts are conditional on the fitted or
fixed rates and observations. They count jumps in the supplied model, not
molecular mutations. Marked count calculations can take substantially longer;
the implementation does not derive p-values from these counts.

## Tip observation model

Each tree tip has a binary material presence vector with values 0, 1 or unknown. A
present call constrains latent material to state 1. An absent call permits latent
state 0 or 2. Unknown material contributes no discrimination. Feature values
record evidence about availability. Positive observations support availability;
negative observations constrain it only where the feature was surveyed. A
missing feature record is unknown. Optional sensitivity and specificity apply
only to features explicitly surveyed at that tip; without those parameters the
observation is a hard compatibility constraint. Surveyed-feature detection
factors are conditionally independent given the latent state.

Optional observed transcript paths record ordered exon-feature paths and their
connected splice links. Multiple paths may coexist at a tip. The model does not
enumerate every path allowed by a state or estimate usage/frequency weights.
Alternative MXE features can remain simultaneously available where no observed
path requires both alternatives.

## Input and output

The strict JSON schema is `intraphy.exon-locus-model/1`. It records catalogue
objects, root distribution, observation provenance, global rate groups,
conditional independence scope, model provenance, branch-length unit and tree
provenance. The CLI accepts rooted species trees as TSV or Newick. Every nonroot
branch needs a finite nonnegative length; the exact input tree remains fixed.

`locus_fit.json` records log likelihoods, fitted rates, optimizer status,
iterations, boundary groups, curvature diagnostics, workers and rate/tree/unit
provenance. `locus_history.json` records reachable states, node marginals for
material lifecycle, complete copy-material retention (`copy_intact`), and
feature availability, plus the opportunity
records. When `--expected-edits` is selected, it also records expected branch
event counts by opportunity and rate group. Marginals are conditional on all
declared inputs.

The declared units are assumed conditionally independent given shared global
rates. This assumption and its provenance are stored with the analysis. IntraPhy
does not infer independence from the number of genes; overlapping material or
event support should be represented in one connected unit.

The likelihood is conditional on the supplied catalogue and observed loci. It
does not generally correct for catalogue discovery/ascertainment or integrate
unobserved extinct copy positions. Labelled copy slots encode supplied position
and homology hypotheses, not a reconstructed copy genealogy; a branch-homogeneous
process permits the same target slot to arise independently on separate
lineages. Availability alone does not imply transcript usage or full-length
isoform frequency.

## Numerical implementation references

The likelihood uses SciPy's sparse matrix-exponential action
[`expm_multiply`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.expm_multiply.html),
which computes the action of a matrix exponential on a vector without requiring
an explicitly materialized transition matrix. The phylogenetic maximum
likelihood framework follows established software practice described by
[HyPhy 2.5](https://doi.org/10.1093/molbev/msz197); this citation is methodological
context, not a claim of using HyPhy code.
