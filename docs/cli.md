# Command-line guide

## Genomic evidence

`check`, `build-case`, and `derive-tables` are the staged interface from
genomic FASTA and existing GFF3/GFF/GTF annotations to correspondence evidence.

```bash
intraphy check --fasta genomes --gff annotations \
  --orthologs families --species-tree tree.nwk
intraphy build-case --fasta genomes --gff annotations \
  --orthologs families --species-tree tree.nwk \
  --output-dir prepared_case --threads 8
```

The `--orthologs` selection is needed when annotations contain multiple
genes; it does not itself establish orthology. Omit it for inputs containing
exactly one selected gene locus per species. `check` resolves identifiers,
tree tips, and coordinate bounds. `build-case` prepares loci and correspondence
evidence. Run `derive-tables` separately only when recomputing tables into
another output directory. See
[genomic input contracts](inputs.md) for file naming, compressed input,
coordinate handling, and optional `extract-loci`.

`derive-tables --short-alignment-max-dp-cells N` sets the maximum dynamic-
programming cell budget for each internal short alignment; the default is
250000. An alignment exceeding the budget fails in full, and DP states are
never truncated. Use `--short-alignment-max-dp-cells unlimited` to remove
IntraPhy's DP cell guard; available memory and other system resources still
limit the alignment.

Annotation exons and introns are evidence at native coordinates. Repeated
transcript records covering the same DNA must be deduplicated by physical
interval when defining copy units. Lack of an annotation alone does not
establish DNA absence. Ambiguous alignment or unsurveyed sequence is unknown.

## Default exon-structure histories

```bash
intraphy analyze --fasta genomes --gff annotations \
  --orthologs families --species-tree tree.nwk --output-dir exon_result
```

The default model is `exon-structure-ctmc`. It maps annotated physical exon
spans using genomic sequence similarity and local genomic geometry, then fits
the finite edit graph over the supplied tree. Identical physical intervals
listed by multiple transcripts are deduplicated. Distinct overlapping
annotations remain separately represented; conflicting evidence can leave a
local observation unknown. Missing annotations are unknown. This includes
annotated terminal and UTR exons; an exon end alone does not imply splice
function.

The default fit estimates one nonnegative scalar rate per family, shared across
all elementary edit kinds and linked local units. Linked units contribute a
composite likelihood. The root geometry prior is uniform over valid geometries
conditional on material origin. Separately, each root or branch
material-origin opportunity has weight one. Rates are conditional on the finite
mapped catalogue and supplied tree. Equal rates provide a constrained baseline
for this edit graph; they do not establish equal biological rates. The model
does not estimate root frequencies or perform a significance test. Annotation
discovery is not corrected for ascertainment, and inferred counts are model
transitions rather than molecular lesions. See the
[genomic exon-span model](genomic_exon_model.md).

Every default genomic exon run writes `branch_exon_changes.tsv`. It reports
branch probabilities for exon-span structure change and declared material-tract
DNA-presence change, along with joint modal endpoint configurations and every
exact tie. It uses no probability threshold and does not count branch events.
Reported spans are local alignment coordinates (0-based, half-open); use
`alignment_offset` with `exon_coordinates.tsv` and
`exon_correspondence.tsv` to connect them to native exon coordinates. DNA
presence is summarized only for declared material tracts. `--expected-edits`
additionally writes `branch_exon_events.tsv` with conditional model transition
counts.

No state-count cap is imposed by default. `--max-states` and
`--max-origin-scenarios` provide explicit stopping limits. To reuse a genomic
exon catalogue, supply `--exon-configurations` with `--input-dir`; every record
must declare `observation_unit="genomic_exon_spans"`. Legacy transcript-path
catalogues are rejected for this route. A fixed rate file requires
`--exon-rates FILE --parameter-mode fixed`.
`--threads` defaults to 1. During genomic catalogue preparation it sets the
per-family MAFFT and minimap2 thread counts; MAFFT receives `--thread N
--threadit 0`, with `threadtb` left at its default. The same value controls
parallel independent family fits.

```bash
intraphy analyze --input-dir prepared_case \
  --exon-configurations exon_configurations.jsonl \
  --output-dir exon_result
```

An optional `--species-tree` replaces the prepared tree only when it has the
same tip labels. Its topology and branch lengths condition inference. Without
catalogue reuse, MAFFT and minimap2 are required for family alignment and
nucleotide correspondence.

## Exon endpoint calibration

```bash
intraphy calibrate-exons --scenario geometry --rate 0.2 --units 2 \
  --replicates 20 --seed 101 --output-dir calibration
```

The `shared-deletion` scenario uses a material tract spanning two exons. Each
replicate contains independent local units sharing a generated rate. Outputs
compare fixed-rate reference inference with inference after estimating the rate
from tips, state spaces, and tree. Failed or nonidentified fits remain counted
per arm. This calibration is conditional on its fixed catalogues; FASTA/GFF
catalogue discovery is not calibrated.

## Separate DNA-presence histories

```bash
intraphy analyze --model dna-presence-ctmc --fasta genomes --gff annotations \
  --orthologs families --species-tree tree.nwk --output-dir dna_result
```

This explicit model infers binary histories for qualified cross-species DNA
position-homology groups, with family-shared gain and loss rates and an
observed-at-least-one ascertainment condition. Unqualified candidates and
failed or over-budget surveys remain unknown; a second optimal placement is
ambiguous and also remains unknown. The bounded survey covers selected
sequence pairs, not genome-wide search space. The default root distribution is
stationary; `--root-frequency fixed --root-presence P` supplies a fixed root
presence probability. `--dna-gain-rate` and `--dna-loss-rate` must be supplied
together to evaluate fixed rates. Rates use the branch-length units of the
supplied rooted tree, whose non-root branches must have explicit lengths.
`--expected-edits` optionally reports conditional transition counts.

These histories describe DNA-position presence changes and do not identify a
molecular duplication source or mechanism. They do not estimate exon boundary
displacement, transcript use, or splicing.

For staged reuse:

```bash
intraphy prepare-genomic-evidence --input-dir prepared_case \
  --output-dir dna_evidence --threads 8
intraphy analyze --model dna-presence-ctmc --input-dir prepared_case --genomic-evidence-dir dna_evidence \
  --output-dir dna_result
```

Reuse checks the prepared physical-site/member catalogue, genome resource paths,
tree tip panel, survey thresholds, and evidence schema. It does not repeat the
alignment survey; preserve and review staged evidence and its source files.
Evidence preparation validates the rooted tree topology and tip panel but does
not require branch lengths; phylogenetic inference still requires explicit
non-root branch lengths.

For `analyze --input-dir`, `--species-tree` may explicitly replace the prepared
tree when it has exactly the same tip labels. The replacement rooted topology
and branch lengths condition the inference; prepared loci, species roster, and
staged observations remain unchanged. Case and reusable-evidence checks are
performed against the prepared tree before inference uses the replacement.
The chosen tree is written to the result's `species_tree.tsv`, and the supplied
tree path is recorded in `execution.json` arguments.

[Glick et al. 2024](https://doi.org/10.1093/molbev/msae248) analyze species-level structural summaries,
whereas this workflow models qualified per-site observations. Binary intron-position reconstruction such as
[Csűrös et al. 2011](https://doi.org/10.1371/journal.pcbi.1002150) and homologous exon/domain structure comparison such as
[ExonEvo](https://doi.org/10.1038/s41467-025-66816-3) provide related
methodological context, not validation of this observation process or event
interpretation.

## Intron-position histories

```bash
intraphy analyze --model intron-position-ctmc --fasta genomes --gff annotations \
  --orthologs families --species-tree tree.nwk --output-dir intron_result
```

This is a separate binary model of annotated genomic intron presence at coding
positions projected through a protein-family MSA. Its states refer to an
annotated intron at the aligned position or consecutive genomic coding bases;
unsupported projections and conflicting physical boundaries remain unknown.
The MSA is conditional on the selected family and available coding paths, and
does not certify orthology. Rates are family-shared, and the stationary root
distribution is the default; `--root-frequency fixed --root-presence P` sets a
fixed root probability. `--gain-rate` and `--loss-rate` evaluate paired fixed
rates in the supplied tree's branch-length units. The retained `--dna-gain-rate`
and `--dna-loss-rate` spellings are aliases for the selected binary model.

Prepare reusable observations with
`intraphy prepare-genomic-evidence --character-type intron-position`; the
result has separate `intron_families.tsv`, `intron_positions.tsv`,
`intron_observations.tsv`, and `intron_alignment.tsv` files. The family summary
records zero-site families without a synthetic character. The fit and history
are `intron_fit.json` and `intron_history.json`. DNA-presence and intron-position observations are fitted
separately, with no joint likelihood. The intron model does not estimate RNA
transcript usage, ancestral intron length, rate heterogeneity, or genome-wide
intron density. See [intron-position method](intron_position_model.md) for
state, phase, and evidence details.

For staged `analyze --input-dir`, an explicit `--species-tree` replacement is
allowed only with the exact prepared tip panel. Prepared case and reusable
intron-evidence coherence are checked against the original tree first; the
replacement tree then supplies the rooted topology and branch lengths for
inference. Prepared loci and intron observations are unchanged. The chosen
tree is written to the result's `species_tree.tsv`, and its supplied path is
recorded in `execution.json` arguments.

## Advanced source-directed DNA-copy model

```bash
intraphy analyze --locus-model model.json \
  --species-tree tree.nwk --output-dir locus_result
```

`model.json` uses `intraphy.exon-locus-model/2`. It supplies the fixed rooted
tree's branch-length unit and provenance, a root distribution, tip DNA
evidence, rate groups, and a finite event catalogue. Every catalogue names
`state_model: binary` or `state_model: irreversible`. The current model
contains DNA material, copy positions, source-qualified `copy_duplication`,
and continuous `dna_deletion` opportunities. It has no feature, transcript
path, or mutually exclusive exon field. Correspondence tables are not
automatically converted into this catalogue.

`--parameter-mode fixed` evaluates supplied rates; the default fits groups
declared `fit`. `--threads` parallelizes independent units, while one
connected unit is serial. `--expected-edits` adds conditional branch counts
and can be more expensive. Results include `locus_fit.json`,
`locus_history.json`, and a run ownership record. State and event summaries
are conditional on the declared tree, root, mapping, catalogue, and
observation design.

`prepare-locus-evidence` applies a complete caller-supplied tip-by-material
DNA evidence array to a model template. Known presence or absence calls need
survey and evidence information; unknown is retained as unknown. Detection
sensitivity/specificity, if supplied, are fixed input probabilities.
`locus-statistics` provides same-data nested rate comparisons, profiles, and
an optional fixed-catalogue DNA bootstrap. See
[conditional statistics](locus_statistics.md).

Additional command routes implement separate experimental analyses. They do
not supply the qualified DNA-copy catalogue or
extend the current locus model to intron-boundary evolution.
