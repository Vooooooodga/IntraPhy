# RpL32 Genomic Exon Case

This native genomic-exon case uses the five species listed in `manifest.tsv`
and `species_tree.tsv`: *D. melanogaster*, *D. simulans*, *D. erecta*,
*D. yakuba*, and *D. teissieri*. It prepares evidence from the supplied local
genome and annotation files, then fits the default genomic exon-span CTMC.

```bash
intraphy build-case \
  --manifest examples/real_cases/rpl32_control/manifest.tsv \
  --species-tree examples/real_cases/rpl32_control/species_tree.tsv \
  --output-dir work/rpl32_control_case \
  --threads 1

intraphy analyze \
  --input-dir work/rpl32_control_case \
  --output-dir results/rpl32_control \
  --model exon-structure-ctmc \
  --parameter-mode fit \
  --branch-length-mode supplied \
  --threads 1
```

The case retains all annotated transcript records. Identical physical exon
intervals are deduplicated; conflicting overlapping annotations can remain
unknown. The workflow does not infer transcript usage. The gene identifiers
and locus mappings are curated inputs in the manifest; this case command does
not establish orthology.

All non-root lengths in the supplied tree equal one, so fitted parameters and
posterior histories are conditional on that arbitrary scale, the supplied
topology, the annotations, and this five-species sample. Protein conservation
does not establish invariance of physical exon structure. This is a
comparative locus case, with no prespecified near-zero rate or true event
expectation. `truth_events.tsv` contains no event records and supplies no
scoring denominator.
