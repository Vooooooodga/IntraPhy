# IntraPhy architecture

The default inference path consumes a strict JSON model and a fixed rooted tree.
Domain objects compile the declared catalogue into its complete root-seeded
reachable closure; sparse CTMC likelihood, fitting and posterior summaries remain
separate from input parsing and output serialization.

| Responsibility | Owner |
|---|---|
| CLI model routing and preflight | `commands/parser.py`, `commands/locus.py`, `commands/preflight.py`, `cli.py` |
| Conditional rate statistics and evidence-table preparation | `commands/locus_statistics.py`, `inference/locus_statistics_run.py`, `inference/locus_evidence.py` |
| Strict JSON records | `inference/locus_codec.py` |
| Model assembly and fixed-tree input | `inference/locus_io.py` |
| Labelled copies, material, features and opportunities | `structure/locus_types.py` |
| Reachable finite event process | `structure/locus_process.py`, `locus_events.py` |
| Tip observation compatibility and surveyed detection | `structure/locus_observations.py` |
| Cross-record biological consistency checks | `structure/locus_validation.py` |
| Sparse likelihood, posteriors and marked counts | `inference/locus_likelihood.py` |
| Shared nonnegative rate fitting and diagnostics | `inference/locus_rates.py` |
| Model-specific history and fit outputs | `inference/locus_run.py` |
| Raw genomic preparation and compatibility models | `inputs/`, `preparation/`, `commands/exons.py`, existing inference modules |

The CLI preflight loads and validates a qualified locus model before it reserves
an output directory. This route does not invoke sequence alignment tools or
prepare FASTA/GFF inputs. Raw genomic input preparation remains a separate,
explicit path and the old model names remain available for those prepared inputs.

The main result artifacts are `locus_fit.json` and `locus_history.json`; their
meaning and conditional assumptions are defined in
[the model documentation](exon_structure_model.md). A locus-specific graphic
renderer is not yet provided; `visualize` reports this clearly for these results.
