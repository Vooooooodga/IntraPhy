The default one-command route is `intraphy analyze --fasta ... --gff ...
--species-tree ... --output-dir ...`; add `--orthologs` for multi-gene whole-
genome annotations. Gene-only selected loci remain unknown, not exons. For
staged `build-case`, use `--allow-unannotated-loci` explicitly when retaining
such loci.

# Genomic FASTA, annotation and a rooted tree

## Standard file inputs

Provide per-species genomic FASTA and GFF3/GFF/GTF files, plus a rooted Newick or
TSV species tree. Plain and gzip files are accepted. File stems match tree tips:
`Species_A.fa`, `Species_A.gff3`, and the tree tip `Species_A`. Direct files and
flat directories both work. Directory searches are not recursive. Extra species,
duplicated resources and missing counterparts are rejected.

```bash
intraphy check --fasta genomes --gff annotations \
  --orthologs OG0001.fa --species-tree species_tree.nwk
```

`check` validates inputs without requiring an alignment executable. An example
FASTA header is `>rna123 gene=gene123`; an exact `rna123` GFF transcript or its
explicit `gene123` metadata must resolve to one gene. Protein IDs explicitly
present in GFF attributes are also supported. When the same ID occurs in several
species, use `>Species_A|rna123`. No gene-symbol guesses, prefix-only matching or
sequence search silently determine the selected orthologs.

An ortholog FASTA is an upstream **selection**, not a genomic substitute. Its
sequence can be CDS, transcript, protein or genomic sequence. IntraPhy does not
infer or certify gene orthology from this sequence. One FASTA file defines one
family, named after the file; `--family-id` overrides a single family's name.
A directory of FASTA files selects several families. Each family must resolve to
one distinct locus in every supplied species. Multiple isoforms of that same
locus are collapsed for locus counting, preserving their IDs. A locus assigned
to two input families is rejected.

Explicit, consistent gene records and descendant relationships are required.
GTF/GFF lacking gene-level records must first be normalized (for example with
AGAT) and inspected; IntraPhy does not silently infer a replacement gene model.
Without `--orthologs`, every GFF must contain exactly one target gene record.
This route is suitable for already selected locus annotations. A whole-genome
GFF containing many genes requires the ortholog selection. An automatic
`input_targets.tsv` records the resolved genes and absolute resource paths.
Users do not need to construct this table.

## One combined genomic FASTA

A single multi-record genomic locus FASTA may be paired with separate species
GFFs. Species come from GFF file stems; each selected species must use different
FASTA sequence IDs. GFF column 1 matches the FASTA record ID, and GFF coordinates
refer to that record. Duplicate FASTA record IDs are rejected.

```bash
intraphy build-case --fasta orthologous_genomic_loci.fa \
  --gff locus_annotations --species-tree species_tree.nwk \
  --output-dir work/one_gene
```

Do not concatenate genome FASTAs with duplicate contig names such as `chr1` and
expect the program to determine their species. Use the per-species file route
in that situation. Identical input sequence content alone does not establish or
refute orthology; identifiers and original genomic coordinates still matter.

## Genomic sequence and coordinate consistency

A mature transcript/CDS FASTA lacks introns. A protein FASTA lacks genomic bases.
Neither permits recovery of missing introns, their exact lengths, or upstream and
downstream sequence. Use such files with `--orthologs`, together with genomic
DNA under `--fasta`. Source sequence IDs, bounds and a target DNA-alphabet probe are
checked. These checks cannot establish that an incorrect-but-plausible assembly
and annotation actually belong together; that remains an input-quality issue.

Genomic input can be a whole assembly or a contiguous locus already extracted
from it. A cropped sequence requires matching rebased annotation. Changing only
a FASTA header does not correct chromosome-based GFF coordinates. Minus-strand
genes must retain a consistent sequence orientation and GFF strand/phase.

## Intron-position evidence

The `intron-position-ctmc` route uses selected genomic loci, annotated CDS paths,
and intron intervals. A candidate intron position must project between adjacent
coding bases in the family protein MSA and have supported local amino-acid
flanks. State 1 requires the intron interval to exactly span adjacent coding
endpoints whose bases agree with the strand-oriented genome; intron interior
sequence is not assessed. State 0 requires adjacent coding bases in the genome.
An endpoint `N` or mismatch remains unknown; a missing declared FASTA resource
or backend failure stops preparation. Missing CDS,
unsupported projection, competing physical boundaries, or unavailable sequence
yields `unknown`. The route does not infer an exon for a gene-only annotation.

GFF CDS phase is annotation metadata used while parsing coding segments. The
reported `coding_right_base_offset` is the position of the right coding base
within its codon at the projected MSA boundary; it is distinct from the GFF
phase field. At a consistent internal CDS junction, downstream phase is
`(3 - coding_right_base_offset) % 3`, following the
[Sequence Ontology GFF3 specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md).
The model reports intron presence at a homologous coding position,
not intron-length change or transcript usage. See the
[intron-position method](intron_position_model.md).

The existing interval reader uses a supplied `.fai` when available and otherwise
streams the FASTA. It does not require a user-created index or read every complete
genome into a multi-species in-memory sequence dictionary. Repeated unindexed
access to large genomes can nevertheless be I/O-intensive.

## Automatic flanks and optional locus export

`build-case` extracts the gene interval and available flanking sequence. It uses
all selected transcript paths when determining the target bounds. `--flank`
defaults to 1000 bases on each side; `--max-extension` bounds subsequent search
extension. These are explicit search limits with no claimed universal biological
optimum. The program cannot extend beyond the supplied contig or locus crop.

To export small portable inputs:

```bash
intraphy extract-loci --fasta genomes --gff annotations --orthologs OG0001.fa \
  --species-tree species_tree.nwk --flank 1000 --output-dir loci

intraphy check --fasta loci/OG0001 --gff loci/OG0001 \
  --species-tree loci/OG0001/species_tree.nwk
```

The export retains the entire selected locus including introns, always in the
forward genomic orientation. It rebases the GFF to a new local sequence ID,
retains negative-strand and phase annotations, and writes
`locus_coordinate_map.tsv`. The conversion is:

```text
source_position = exported_position + source_start - 1
```

Requested and actual upstream/downstream bases are reported in transcriptional
orientation. A truncated flank is reported explicitly. Only the selected gene's
annotation hierarchy is exported; overlapping neighboring genes are not included
in this portable target-only GFF. Use full genomic input when that background is
important. Both original and exported coordinates remain auditable.

## Optional AGAT preparation

AGAT is useful for normalizing GFF/GTF and extracting annotated sequences. It is
not bundled, silently installed, or used as an independent biological validator.
The opt-in adapter calls the real executable:

```bash
intraphy normalize-annotation --gff annotations --output-dir normalized_annotations
```

This invokes `agat_convert_sp_gxf2gxf.pl --gff INPUT --output OUTPUT` per file.
An optional `--config` passes an explicit AGAT configuration. Commands, source and
output paths, stdout/stderr and return status are retained. No fake executable or
internal substitute is used when AGAT is missing. Normalization can change IDs,
parent relationships, boundaries and inferred missing features. Rerun `check`
against the resulting annotation and inspect changes before inference.

An independently installed AGAT can also prepare the upstream sequence files:

```bash
# Example: extract annotated CDS for an upstream orthology workflow.
agat_sp_extract_sequences.pl --gff Species_A.gff3 --fasta Species_A.fa \
  --type cds --output Species_A.cds.fa
```

AGAT `gene=` FASTA metadata is recognized by the exact-ID resolver. A FASTA
containing all CDS is not automatically an ortholog family: the user still
supplies the upstream homolog set. AGAT ordinarily reverse-complements features
on the minus strand; pairing that output with original chromosome-coordinate
GFF is unsafe. Use IntraPhy's coordinate-safe locus export or an explicitly
consistent external conversion. `--type exon --merge` creates spliced sequence
and must not be supplied as genomic DNA for structure analysis.

The AGAT adapter has interface tests. A real AGAT executable was not present in
the local development environment; no real-AGAT integration or biological
validation is claimed by those mocked-process tests.

Official documentation:
[normalization](https://agat.readthedocs.io/en/latest/tools/agat_convert_sp_gxf2gxf.html),
[sequence extraction](https://agat.readthedocs.io/en/latest/tools/agat_sp_extract_sequences.html).

## Advanced table and OrthoFinder import

The existing `--manifest` and OrthoFinder import interfaces remain available to
pipeline developers. They are optional, mutually exclusive alternatives to direct
file input; a hand-written table is not the standard user interface. See the
[table/schema reference](input_format.md). No runtime alias for the old project
name is retained.
