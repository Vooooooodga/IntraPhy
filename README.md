# IntraPhy

IntraPhy addresses how exon-copy DNA and splice-feature availability changed
along a species tree. It estimates event rates and reconstructs ancestral
material, copy-position and feature states from a supplied root distribution and
tip evidence. The model is a joint continuous-time Markov process: DNA retention,
copy duplication, splice availability and shared-interval deletion contribute to
the same latent history.

## Quick start

The included Dscam-like example is synthetic. It illustrates duplicated exon
copies, splice-link availability and mutually exclusive alternatives without
claiming an observed biological history.

```bash
intraphy analyze --locus-model examples/exon_locus/dscam_like_model.json \
  --species-tree examples/exon_locus/species_tree.tsv \
  --output-dir results/dscam_like
```

The default `--model exon-locus-ctmc` fits groups marked `fit` by maximum
likelihood and holds groups marked `fixed` constant. Use
`--parameter-mode fixed` to evaluate all supplied rate values without
optimization. The tree and branch lengths stay fixed. `--threads` parallelizes
independent declared units; a single connected locus is evaluated serially.
Add `--expected-edits` for expected event counts by branch and opportunity;
marked calculations can take substantially longer.

## Model objects and events

Each unit declares ordered genomic material tracts, labelled copy positions,
availability features and typed event opportunities. A material tract has latent
lifecycle state 0 (not introduced), 1 (present) or 2 (deleted). A copy slot is
supported by explicit homology/evidence records. Exon and splice-link features
are governed by material requirements and a directed prerequisite graph.
Mutually exclusive exon groups constrain compatible paths while allowing
multiple alternative features to remain available in the latent state.

The supported opportunities are evidence-qualified source-to-target copy
duplication, deletion of one continuous shared DNA interval, and
context-specific splice-feature changes. Deletion removes a continuous supplied
interval once; dependent features are removed through the declared dependency
rules. The catalogue is finite and conditional: omitted opportunities are
unavailable. Each labelled copy is a supported position and homology hypothesis;
the model does not reconstruct copy genealogy. Under a branch-homogeneous
process, the same labelled target can arise independently on different
branches.

## Likelihood and observations

The process starts from the supplied root distribution and compiles its full
reachable closure under the listed opportunities. For state i and j, the
generator entry is qᵢⱼ = Σ(rate[group] × outcome weight) over eligible declared
events from i to j, with diagonal entries set to minus the outgoing rate. The
branch transition matrix is P(t) = exp(Qt); pruning sums over unobserved states
to obtain the tip likelihood. Sparse matrix-exponential actions evaluate this
likelihood. Across conditionally independent units, log-likelihood contributions
are summed. Free global rates are optimized by nonnegative maximum likelihood;
the optimizer reports convergence, boundary rates and local curvature
diagnostics. These are local fit diagnostics, not a proof of global optimality
or parameter identifiability. A single gene can be fitted with the same
diagnostics.

Tip DNA observations report binary presence/absence or unknown. Absence does not
distinguish a never-introduced tract from a deleted tract. Feature observations
record availability evidence, not expression, usage, PSI, protein abundance or
selection. An omitted feature annotation is unknown. Optional
sensitivity/specificity values apply only to explicitly surveyed features.
Observed transcript paths constrain connected exon/splice links. Alternative
features can be jointly available when they do not conflict within an observed
path. IntraPhy does not estimate transcript usage or enumerate every possible
full transcript.

## Inputs and outputs

`--locus-model` requires schema `intraphy.exon-locus-model/1`; identifiers and
declared record consistency are checked, while biological evidence and
assumptions remain user-supplied. `--species-tree`
accepts the existing rooted TSV or Newick conversion route. Every non-root
branch requires a finite nonnegative length. `branch_length_unit` and
`tree_provenance` are required in the JSON model. Omitted biological evidence is
not filled from a tip union or inferred annotation.

The output directory contains:

- `locus_fit.json`: fit status, log likelihood, rate estimates, boundary status,
  local curvature/flatness diagnostics, and conditional unit-independence scope.
- `locus_history.json`: model assumptions and provenance, complete reachable
  state records and ancestral copy/material/feature marginal probabilities.
  Add `--expected-edits` for branch expected event counts by opportunity and
  rate group; marked counts can take substantially longer.
- `run_result.json`, `execution.json`, and the command log: model ownership and
  command-session status.

Expected counts describe jumps in the supplied opportunity catalogue, not
molecular mutations; expected counts alone are not significance tests. The model
does not supply a catalogue-discovery correction, integrate extinct unobserved
copies, or test sociality. Conditional rate tests are described below.

## Conditional rate statistics and tip evidence

`locus-statistics` compares a full rate fit with the same-data nested fit that
fixes one or more declared free rate groups at zero. Repeat
`--null-rate-group GROUP` for each tested group. Optional
`--profile-rate-group GROUP --profile-values VALUE...` profiles one free rate
while refitting nuisance rates. These likelihood comparisons are conditional
on the supplied locus units, tree, root distributions, observations, and event
catalogue; they do not report asymptotic P values. A parametric bootstrap is
available for a fixed catalogue with DNA-only tip observations; latent splice
features may be declared, while tip feature/path observations are excluded. Use
`--bootstrap-replicates B --seed N --sampling-design fixed_catalogue`. The
bootstrap keeps the supplied tree, root, catalogue, and missingness/detection
design fixed. Its calibration is plug-in and conditional, and an incomplete
replicate set has no reported bootstrap P value.

`prepare-locus-evidence` applies a JSON array covering every declared
family/unit/tip/material combination to an unchanged model template. Calls are
`present`, `absent`, or `unknown`; unknown calls remain null. Known calls require
an explicit survey and evidence source. Optional per-material `sensitivity` and
`specificity` values are caller-supplied fixed detection probabilities and are
accepted only for surveyed known calls. No detection probability is estimated
by this command. Without these fields, binary DNA observations remain hard
constraints. The command writes a validated model and parsed species tree,
preserving the template's catalogue, root, opportunities, and feature records.

Both commands write auditable JSON summaries and a `run_result.json` ownership
record. Their results remain conditional on declared inputs and assumptions;
they do not resolve copy homology or infer a locus-discovery correction.

## Scope and interpretation

Results are conditional on the topology, branch-length scale, root distribution,
homology hypotheses, observation process, conditional independence of units and
finite opportunity catalogue. There is no general correction for selecting
observed exons or loci. Surveyed-feature detection factors assume conditional
independence given the latent state. The method does not infer RNA usage or
validate a supplied molecular mechanism. Dscam studies motivate the synthetic
example; they do not validate its states or imply a sociality conclusion.

Raw FASTA/GFF inputs can be checked, selected and prepared through `check`,
`build-case` and `derive-tables`. This preparation does not automatically create
the qualified copy/homology catalogue required by `--locus-model`; users supply
that model explicitly. Prepared-input `run` and `infer-phylogeny` retain
`--model exon-parsimony` and `--model exon-ctmc` as explicit compatibility paths.
The explicit
`infer-exon-repertoires` command is a separate conditional engine. These routes
do not silently substitute for the exon-locus model.

## References

- Pond SLK et al. 2020. HyPhy 2.5—A Customizable Platform for Evolutionary
  Hypothesis Testing Using Phylogenies. *Molecular Biology and Evolution*
  37:295–299. [doi:10.1093/molbev/msz197](https://doi.org/10.1093/molbev/msz197).
- Lee C, Kim N, Roy M, Graveley BR. 2010. Massive expansions of Dscam splicing
  diversity via staggered homologous recombination during arthropod evolution.
  *RNA* 16:91–105. [doi:10.1261/rna.1812710](https://doi.org/10.1261/rna.1812710).
- SciPy `expm_multiply` reference:
  [SciPy sparse linear algebra documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.expm_multiply.html).
