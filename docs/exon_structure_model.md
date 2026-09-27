# DNA-copy locus model

## Biological object

The model represents homologous DNA-copy positions within a gene locus on a
fixed rooted species tree. Genomic sequence and existing annotations identify
candidate correspondence: exon and intron coordinates, boundaries, CDS phase,
strand, and local order are evidence. They do not themselves define a gain or
loss. One physical DNA interval is represented once even if several transcript
records annotate it. Missing annotation, failed alignment, and ambiguous copy
assignment are unknown observations; sequence-supported absence requires an
explicitly adequate survey.

This engine estimates the DNA-copy component of intragenic evolution. Annotation
boundaries and intron positions remain important mapping outputs, but are not
separate latent event types in this process. It makes no claim about transcript
use, splice availability, or alternative splicing.

## Catalogue and states

The strict input schema is `intraphy.exon-locus-model/2`. Each conditionally
independent connected unit has ordered, nonoverlapping material tracts on a
shared local coordinate axis and labelled copy positions supported by supplied
homology evidence. Native chromosome coordinates require explicit mapping to
that axis. Orientation and order are correspondence evidence; the model does
not infer inversion events or copy genealogy.

Each catalogue declares one `state_model`:

| State model | Values | Consequence after deletion |
| --- | --- | --- |
| `binary` | 0 absent, 1 present | A tract returns to 0; a later source-qualified duplication can restore presence. |
| `irreversible` | 0 not introduced, 1 present, 2 deleted | A tract enters 2 and cannot regain presence on that lineage. |

The irreversible state definition concerns a labelled position along a
lineage. Different branches can independently introduce the same position, so
the model does not impose a single origin over the tree. Binary regain still
requires a present, qualified source; the full binary process need not be
reversible. Tip presence means
state 1. Tip absence is state 0 in the binary model and is compatible with 0 or
2 in the irreversible model. Unknown tip calls impose no presence constraint.
Known calls may use caller-supplied fixed sensitivity and specificity when the
survey design supports them; these probabilities are not estimated from the
same call.

The root distribution is supplied explicitly. The full set of states reachable
from its support under the declared opportunities is compiled without an
arbitrary state-count cap. Memory and computation can grow steeply with tract
number and event coupling.

## Events and likelihood

A `copy_duplication` opportunity names a qualified source and target position.
Its source material must be present and target material eligible for gain under
the selected state model. A `dna_deletion` opportunity covers one continuous
interval and changes every intersecting present tract together. The event is
counted once even when several tracts are lost. These typed events define a
finite, conditional catalogue; opportunities not declared cannot occur in the
fitted process. Catalogue qualification rests on sequence correspondence and
independent annotation/context evidence, not on a favorable inferred history.

For distinct reachable states \(i,j\), the generator is
\(q_{ij}=\sum_{e:i\to j} r_{g(e)}w_e\), summing rates of eligible declared
event outcomes. The diagonal is minus the outgoing sum. Branch transitions are
\(P(t)=\exp(Qt)\). Pruning sums over unobserved internal states and applies tip
DNA emission probabilities to obtain the likelihood. Log likelihoods of
declared independent units are summed. Overlapping tracts or coupled events
belong in one unit.

Global rate groups can be fixed or estimated by nonnegative maximum likelihood.
Branch lengths and their unit, tree topology, and root distribution are fixed
inputs. Rate values depend on that branch-length scale and on the number of
declared opportunities. Fit status, boundary estimates, and local curvature
diagnostics describe the returned fit; they do not prove global optimality or
identifiability. Conditional branch counts, requested by `--expected-edits`,
count model transitions rather than individual molecular lesions.

## Interpretation and limits

`locus_fit.json` contains likelihood and rate diagnostics. `locus_history.json`
contains ancestral tract/copy DNA marginals and optional branch event counts.
There are no feature, path, or splice-state marginals. Inference remains
conditional on the supplied tree, root, homology and position hypotheses,
survey design, opportunity catalogue, and independence of units. The method
does not automatically correct for discovering the catalogue from the same
observed loci or integrate unobserved extinct copy positions. A likelihood
comparison across different codings or sampled loci requires a common
observation space and explicit treatment of discovery; a larger raw likelihood
alone is not evidence that one biological model is better.

The likelihood is an application of established discrete-state phylogenetic
methods. [Malin](https://doi.org/10.1093/bioinformatics/btn226) applies
maximum likelihood to corresponding intron positions; the
[ExOrthist](https://doi.org/10.1186/s13059-021-02441-9) framework motivates
using exon sequence and genomic context when establishing correspondence.
