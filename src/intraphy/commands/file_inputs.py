"""Shared command-line definitions for manifest-free genomic inputs."""


def add_file_inputs(command, *, required=True):
    command.add_argument("--allow-unannotated-loci", action="store_true",
                         help="Keep gene-only loci as unknown; never treat a gene span as an exon.")
    source = command.add_mutually_exclusive_group(required=required)
    source.add_argument("--fasta", nargs="+", metavar="PATH",
                        help="Genomic FASTA files or a flat directory; Species.fa matches Species.gff3.")
    source.add_argument("--manifest", help="Optional advanced TSV input; ordinary runs do not need a table.")
    command.add_argument("--gff", nargs="+", metavar="PATH",
                         help="GFF3/GFF/GTF files or directory with species stems matching --fasta.")
    command.add_argument("--orthologs", nargs="+", metavar="FASTA",
                         help="Upstream homolog FASTA files (one family/file). IDs select GFF loci; "
                              "the sequences may be CDS/protein but are not genomic absence evidence.")
    command.add_argument("--family-id", help="Family label for a single selected family; "
                         "defaults to the ortholog filename or target_gene for one-gene GFF inputs.")
