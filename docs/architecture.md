# IntraPhy architecture

Genomic preparation starts with FASTA, existing GFF3/GFF/GTF annotations, and
a supplied species tree. It resolves selected loci, extracts annotation and
sequence evidence, and prepares correspondence tables. The DNA inference path
then consumes a separately qualified `intraphy.exon-locus-model/2` JSON
catalogue and fixed rooted tree. Preparation does not automatically decide
which copy positions, duplications, or deletion intervals are evolutionarily
qualified. Domain objects compile the declared catalogue into its complete
root-seeded reachable closure; sparse CTMC likelihood, fitting, and posterior
summaries remain separate from parsing and serialization.

| Responsibility | Owner |
|---|---|
| FASTA/GFF locus selection, preparation and correspondence evidence | `inputs/`, `case.py`, `preprocess.py`, `mapping/`, `commands/parser.py` |
| CLI model routing and preflight | `commands/parser.py`, `commands/locus.py`, `commands/preflight.py`, `cli.py` |
| Conditional rate statistics and evidence-table preparation | `commands/locus_statistics.py`, `inference/locus_statistics_run.py`, `inference/locus_evidence.py` |
| Strict JSON records | `inference/locus_codec.py` |
| Model assembly and fixed-tree input | `inference/locus_io.py` |
| Labelled copy positions, DNA material, state model and opportunities | `structure/locus_types.py` |
| Reachable finite event process | `structure/locus_process.py`, `locus_events.py` |
| Tip observation compatibility and surveyed detection | `structure/locus_observations.py` |
| Cross-record biological consistency checks | `structure/locus_validation.py` |
| Sparse likelihood, posteriors and marked counts | `inference/locus_likelihood.py` |
| Shared nonnegative rate fitting and diagnostics | `inference/locus_rates.py` |
| Model-specific history and fit outputs | `inference/locus_run.py` |
| Additional historical or compatibility analyses | `commands/exons.py`, other existing inference modules |

The model preflight loads and validates the user-qualified DNA catalogue before
it reserves an output directory. The `check`, `build-case`, and
`derive-tables` route handles genomic inputs separately. Annotation exon
boundaries, intron positions, and CDS phase remain mapping evidence; the
current locus likelihood reconstructs DNA-copy states and events. Other
command routes retained in the package have distinct assumptions.

The main result artifacts are `locus_fit.json` and `locus_history.json`; their
meaning and conditional assumptions are defined in
[the model documentation](exon_structure_model.md). A locus-specific graphic
renderer is not yet provided; `visualize` reports this clearly for these results.
