# Method schematics

The existing publication figures illustrate correspondence and several broader
structural hypotheses. Their panels showing transcript repertoires,
splice-feature availability, split/fusion/boundary edits, or transcript
evolution are conceptual context and do not depict the current DNA-copy
likelihood. The current inference has source-qualified copy duplication and
continuous DNA deletion only. Figure labels or state symbols should not be
read as output fields of `intraphy.exon-locus-model/2`.

## Mapping evidence

![Schematic genomic and coding-projection tracks followed by ordered candidate correspondence and separate symbols for unknown and supported absence.](figures/publication/exon_correspondence.svg)

[SVG](figures/publication/exon_correspondence.svg) ·
[PDF](figures/publication/exon_correspondence.pdf) ·
[PNG](figures/publication/exon_correspondence.png)

Genomic sequence blocks, exon/intron coordinates, CDS phase, strand and copy
order provide correspondence evidence. An unknown or conflicting match stays
unresolved. A survey-qualified DNA absence call can constrain the locus
model. Multiple annotation transcript records covering the same physical
interval remain one DNA unit.

## Broader conceptual figures

- [Method scheme](figures/publication/method_scheme.svg): includes local exon
  configurations and edit transitions beyond the current likelihood.
- [Structural states](figures/publication/structural_states.svg): depicts
  transcript repertoires; these are outside the current DNA state space.
- [Evolutionary events](figures/publication/evolutionary_events.svg): only its
  continuous shared DNA deletion motif directly corresponds to a current
  event type. Other pictured structural edits are outside this engine.
- [Phylogenetic inference](figures/publication/phylogenetic_inference.svg):
  shows general tree pruning, with parsimony and transcript-configuration
  panels outside the current DNA-copy likelihood.

For the implemented state definitions, events, and result interpretation, use
the [DNA-copy model reference](exon_structure_model.md).
