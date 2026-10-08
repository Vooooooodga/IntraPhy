# 0.20.0 — exon structure and DNA-copy evolution models

- Map physical exon spans across species from genomic FASTA and gene annotations,
  then fit the default exon-structure CTMC on the supplied rooted species tree.
- Fit one nonnegative rate per family, shared across elementary edits and linked
  local units; report ancestral exon configurations and branch endpoint changes.
- Add source-directed DNA-copy evolution through the advanced
  `--locus-model` route, using a supplied `intraphy.exon-locus-model/2` catalogue.
- Preserve raw FASTA/GFF preparation and named compatibility model paths.

# 0.19.1 — exon structure model corrections (2026-09-21)

- Qualify annotation alternatives using complete source paths and local sequence evidence.
- Search all supplied locus DNA for missing-exon candidates and preserve search-range
  information in prepared results.
- Reject overlapping physical regions across input catalogues and rate collections.
- Preserve species-tree branch information in parsimony, likelihood and simulation results.
- Expose root material-origin assumptions and corresponding simulation settings.
- Count geometry before enumeration, cache state indices, and bound numerical cache reuse.
- Represent multi-exon insertions as explicit, evidence-qualified opportunities.
- Support repeated-gene bootstrap samples and report unresolved model comparisons.
- Render ancestral configuration probabilities and branch summaries separately from
  minimum-change histories.
- Configuration model/schema becomes v2. Old catalogues/results are rejected with a
  regenerate-from-original-inputs message, never silently relabeled.

---

# 0.19.0 — exon configurations and elementary edits (2026-09-21)

- New native exon configuration schema; biological units are not P/R/J columns.
- One edit registry for split/fusion, donor/acceptor shifts, exon appearance/
  inactivation and source-aware genomic interval insertion/deletion.
- Complete bounded candidate-state enumeration including ancestral intermediates;
  incomplete spaces stop probability inference. Irreversible deletion and explicit
  introduction-opportunity scenarios preserve material-source identity.
- All-optimal generalized Sankoff histories, representative compatible paths and
  non-additive possible placements. The legacy binary wrapper uses this kernel.
- Arbitrary-state CTMC, log-space likelihood and ancestor/endpoint probabilities;
  optional marked expected edit counts and probability of at least one edit.
- Explicit rate files; pooled scalar/foreground fitting across genes; whole-gene
  bootstrap and conditional parametric testing with discovery/validity gates.
- Manifest-free `analyze`; inherited file inputs/flanks/coordinate export remain.
- Genomic MSA, coding-projection corroboration, exon correspondence, copy/orientation
  rejection, unknown gene-only loci and separate annotation-alternative scenarios.
- Raw structural counterexamples include phase 0/1/2 fusion, insertion, deletion,
  intronization, boundary changes, repeats, negative strand and observation damage.
- Native CDS diagnostics preserve split codons, exceptions and noncoding exons.
- Optional explicit CESAR2 prediction export; predictions never replace observations.
- Read-only tree/exon reports and synthetic guide generated through the real engine.
- V18 binary models remain explicit legacy baselines. No RNA usage analysis added.

No independent biological benchmark or genome-wide discovery-aware statistical
calibration is claimed by this release.

---

# 0.18.0 — structural characters and elementary changes (2026-09-20)

- Remove the previous import namespace and executable alias; require Python 3.10+.
- Add character catalogues, coordinate evidence, dependence metadata and applicability checks.
- Remove formal compound-event summaries; report per-layer minimum character changes and a compatible optimum.
- Block all independent-character likelihood inference when known linked included sites are present.
- Permit explicit unknown tips in independently specified complete catalogues.
- Correct complementary-fragment repeat classification using actual reference coverage.
- Restore conservative CI assertions and add raw FASTA/GFF3 integration regressions.
- Add CLI preflight, explicit output ownership, preserved backups, failure logs and version provenance.
- Split large implementation modules and provide English documentation and a rendered vector methods figure.

No biological-accuracy benchmark or finite-sample statistical calibration is claimed.
Earlier release notes and unchanged historical validation remain in repository history.


### V18 file-input and explanatory-figure completion

- Accept genomic FASTA/GFF files or directories and a supplied rooted tree without
  a hand-written table; exact ortholog-FASTA identifiers select target loci.
- Support combined genomic locus FASTA with distinct species sequence IDs.
- Export paired locus FASTA/GFF with automatic flanks, strand/phase preservation
  and reversible source-coordinate mapping.
- Add optional, explicit AGAT normalization with command/failure provenance.
- Add original tree/structure/ribbon teaching plates and `intraphy explain`;
  keep synthetic guide and data-derived target figures clearly separate.
- Test file contracts, coordinate round-trips, conservative rejection and the
  illustrated parsimony examples. No new biological/statistical calibration.

### Biology-facing method figures and explicit model explanation

- Replace the previous five teaching plates with three focused figures: overview,
  homology inference, and gene structures to a phylogenetic probability model.
- Show the species tree and gene structures together, including sequence loss,
  intron loss, missing annotation and uncertainty rather than only exon splitting.
- Explain compatible-match chain selection and distinguish complementary coverage
  from competing copies. Use named biological characters instead of internal keys.
- Calculate illustrative ancestral-state probabilities and test them against the
  production pruning and posterior algorithms; fixed teaching parameters are explicit.
- Document likelihood construction, missing-state treatment, rate-sharing scope,
  character dependence and conditional inference in `docs/model_bridge.md`.
