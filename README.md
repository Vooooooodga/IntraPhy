# IntraPhy

IntraPhy maps homologous intragenic DNA structure across species and estimates
how qualified copy positions changed on a supplied species tree. Its evidence
starts with genomic sequence and existing gene annotations. Exon and intron
coordinates, boundaries, CDS phase, and copy order help establish
correspondence. Multiple transcript records can provide annotation evidence
for one locus; overlapping records of the same physical DNA do not create
extra DNA copies. A missing annotation or unresolved alignment remains unknown
until there is evidence for absence.

The default genomic-input analysis fits binary presence histories for qualified
cross-species DNA position-homology groups. These histories describe DNA-site
presence gains and losses; they do not identify a molecular duplication source
or mechanism. Genomic exon/intron boundaries remain annotation descriptors.
Presence histories cannot detect intron-loss exon fusion, boundary shifts, or
exonization when homologous DNA remains. IntraPhy does not infer transcript use
or the full evolution of intragenic architecture.

The advanced `--locus-model` route retains explicit, source-directed DNA-copy
duplication and continuous DNA-deletion opportunities in an
`intraphy.exon-locus-model/2` model.

## From genomic inputs to correspondence evidence

Provide per-species genomic FASTA and matching GFF3/GFF/GTF files, with a rooted
species tree. File stems identify species and must match tree-tip names. For a
whole-genome annotation, `--orthologs` selects one gene locus per species for
each family; orthology is supplied by the user. For already selected
annotations containing one gene locus per species, omit `--orthologs`.

```bash
intraphy analyze --fasta genomes --gff annotations \
  --orthologs families --species-tree species_tree.nwk \
  --output-dir dna_result --threads 8
```

`analyze` resolves inputs, prepares selected genomic loci, derives
correspondence evidence, surveys qualified DNA positions, and infers presence
histories. Use `check`, `build-case`, and `prepare-genomic-evidence` when staged
inspection or evidence reuse is useful.
See [genomic input requirements](docs/inputs.md) and the
[command guide](docs/cli.md) for options and output interpretation.

For a separate homologous intron-position history, select
`--model intron-position-ctmc`. This analyzes annotated genomic introns at
coding positions aligned through a family protein MSA; it does not combine
intron states with DNA-material states. GFF phase is retained as annotation
metadata and is distinct from the coding-base offset at the aligned boundary.
The intron model does not infer RNA transcript use, ancestral intron length,
rate heterogeneity, or genome-wide intron density. Zero-candidate families are
reported in `intron_families.tsv` without a fabricated site. See the
[intron-position method](docs/intron_position_model.md).

The automatic mapping route only admits position-qualified homology groups and
keeps unresolved surveys unknown. Its CTMC histories are conditional on these
groups and the supplied rooted tree. A candidate sequence match, unannotated
exon, or unmapped region does not automatically establish an evolutionary
event.

## Advanced source-directed DNA-copy inference

The locus model uses schema `intraphy.exon-locus-model/2`. Each connected unit
declares ordered material tracts, copy positions, an explicit `state_model`,
qualified event opportunities, a root distribution, and tip DNA observations.
The `binary` model has states 0 (absent) and 1 (present); a deletion returns a
tract to 0 and a later source-qualified duplication may reintroduce it. The
`irreversible` model has states 0 (not introduced), 1 (present), and 2
(deleted); state 2 cannot regain presence. The choice is made per catalogue.
Binary regain requires an eligible source, so the full process need not be
reversible.
Both permit the same position to be introduced independently on separate
branches, so the irreversible option is not a stochastic Dollo model.

A copy duplication requires a present, qualified source and an eligible target.
One DNA deletion removes all present tracts intersected by its continuous
declared interval in a single event. The catalogue is finite; omitted
opportunities have zero modeled rate. Full reachable-state closure has no
arbitrary count cap and may become costly for large catalogues.

```bash
intraphy analyze --locus-model examples/exon_locus/dscam_like_model.json \
  --species-tree examples/exon_locus/species_tree.tsv \
  --output-dir locus_result
```

The bundled Dscam-like locus and tree are synthetic teaching inputs. They do
not represent an inferred history for Dscam. Rate groups marked `fit` are
estimated by nonnegative maximum likelihood; groups marked `fixed` retain
their supplied values. `--parameter-mode fixed` evaluates all supplied rates
without optimization. The rooted tree, branch lengths, branch-length unit, and
root distribution are supplied conditions. `--threads` distributes
independent declared units; one connected locus is evaluated serially.
Optional `--expected-edits` calculates branch event counts and can take
longer.

The likelihood sums over latent ancestral DNA states on the fixed tree.
`locus_fit.json` reports fit status, likelihood, rate estimates, and local
diagnostics. `locus_history.json` reports ancestral DNA and copy-presence
marginals, with event counts when requested. Estimates remain conditional on
the supplied correspondence, catalogue, root, tree, observations, and
independence assumptions. A DNA absence call does not establish its historical
cause; unresolved evidence is recorded as unknown. See the
[model reference](docs/exon_structure_model.md) and
[conditional statistics](docs/locus_statistics.md).

## Study status

See [research foundations](docs/research_foundations.md) for literature
context and scope of the structural characters and phylogenetic analyses.

The [ablation protocol](docs/ablation_study.md) specifies planned comparisons
for mapping context, shared deletions, and state-model assumptions. It reports
no completed validation or empirical performance estimate. Existing annotated
genomes can provide independent evidence for correspondence; synthetic
examples and selected observations alone cannot calibrate error rates or P
values.
