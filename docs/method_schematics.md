# Method schematics

The default model is the genomic exon-span CTMC: one local physical exon-span
geometry, including annotated UTR and terminal exons, combined with declared
local DNA-material states. It does not model RNA use, transcript repertoires or
splicing as hidden states. One rate is fitted per family and shared across its
edit kinds and linked units. There is no default state-count cap. For the
implemented default, see the [genomic exon-span model](genomic_exon_model.md).
The publication figures below are scoped illustrations, not a complete diagram
of this default model or of outputs from an empirical run.

## Mapping evidence

![Schematic genomic and coding-projection tracks followed by ordered candidate correspondence and separate symbols for unknown and supported absence.](figures/publication/exon_correspondence.svg)

[SVG](figures/publication/exon_correspondence.svg) ·
[PDF](figures/publication/exon_correspondence.pdf) ·
[PNG](figures/publication/exon_correspondence.png)

Sequence, native coordinates, strand and local order support correspondence;
existing staged protein-projection evidence can provide corroboration. CDS
phase is retained as annotation metadata and does not score or qualify a match.
Unknown or conflicting matches remain unresolved. A survey-qualified DNA
absence call may constrain the separate DNA-copy model. Multiple transcript
annotations of the same physical exon interval are deduplicated as one physical
exon observation.

This figure illustrates correspondence evidence relevant to the default
genomic exon-span model; its symbols are not model output fields.

## Broader conceptual figures

- [Method scheme](figures/publication/method_scheme.svg): a broad conceptual
  sequence from loci to structural summaries. Parsimony and other panels are
  not a depiction of the default fitted analysis.
- [Structural states](figures/publication/structural_states.svg): its local
  repertoire of coexisting configurations and whole-annotation alternatives
  belongs to a separate repertoire scope, not the default single-geometry
  genomic exon-span state.
- [Evolutionary events](figures/publication/evolutionary_events.svg): its
  split, fusion, boundary and material-change motifs illustrate possible
  transitions. They do not represent an inferred path, event count, or
  molecular mechanism. Default branch output is endpoint-based net change.
- [Phylogenetic inference](figures/publication/phylogenetic_inference.svg):
  general pruning illustration. Its parsimony panel is a separate summary;
  pooled one-scale fitting across two or more independent genes is not the
  default per-family fit; expected transition counts are optional.

The advanced source-directed DNA-copy CTMC is a separate model with a distinct
catalogue and interpretation; see the
[DNA-copy model reference](exon_structure_model.md). For default branch
endpoint and optional expected-edit semantics, use the
[genomic exon-span model](genomic_exon_model.md).
