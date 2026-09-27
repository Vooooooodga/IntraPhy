# Quick start

IntraPhy estimates event rates and reconstructs ancestral exon-copy material and
splice-feature availability on a fixed rooted species tree. A run requires a
qualified, user-supplied finite homology/opportunity catalogue in
`intraphy.exon-locus-model/1` JSON and one rooted tree in TSV or Newick format.

Run the small synthetic Dscam-like example from the repository root:

```bash
intraphy analyze --locus-model examples/exon_locus/dscam_like_model.json \
  --species-tree examples/exon_locus/species_tree.tsv \
  --expected-edits \
  --output-dir results/dscam_like
```

Groups declared `fit` use their positive supplied values as starting values for
nonnegative maximum likelihood. Groups declared `fixed` remain constant. Use
`--parameter-mode fixed` to evaluate every supplied value as fixed. The tree,
branch lengths and their units are fixed inputs. `--threads` applies to
conditionally independent locus units; one connected unit is serial.

The model file explicitly supplies the root distribution, homology/copy slots,
DNA material, exon and splice-link features, typed event opportunities,
observation evidence, global rate groups, and provenance. IntraPhy compiles the
full state closure from root support. It records fit diagnostics in
`locus_fit.json` and ancestral material/copy/feature marginals in
`locus_history.json`. The command requests branch event counts with
`--expected-edits`; this calculation is optional and can take substantially
longer.

Raw genomic inputs remain available for upstream preparation. See [input
notes](inputs.md) for FASTA/GFF selection, `check`, `build-case` and
`derive-tables`. Those operations do not automatically infer the qualified copy
homology and opportunity hypotheses required by the locus JSON model. Prepared
inputs can use explicit compatibility choices `--model exon-parsimony` or
`--model exon-ctmc`.

The [model reference](exon_structure_model.md) defines the state space, event
process, likelihood, observation model and inferential scope. The
[architecture](architecture.md) maps those concepts to the code.
