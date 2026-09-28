"""Prepared-input fixtures for genomic intron observation tests."""
from intraphy.aligners.columns import revcomp
from intraphy.storage.tabular import write_tsv


def spliced_fixture(root, *, strand="+", split=30, phase_known=True, aliases=1,
                    gene_only=False, internal_n=True, annotated_intron=True, no_paths=False):
    inputs = root / "inputs"
    inputs.mkdir()
    coding = "ATG" * 20
    intron = "N" * 5 if internal_n else "C" * 5
    species_data = {"A": coding[:split] + intron + coding[split:], "B": coding}
    occurrences, paths, features, loci, introns, fasta_lines, targets = [], [], [], [], [], [], []
    for species, genomic_tx in species_data.items():
        locus_sequence = revcomp(genomic_tx) if species == "A" and strand == "-" else genomic_tx
        fasta = root / f"{species}.fa"
        fasta.write_text(f">chr1\n{locus_sequence}\n")
        gene = f"g{species}"
        targets.append({"family_id": "fam", "species": species, "gene_copy_id": gene})
        loci.append({"species": species, "gene_copy_id": gene, "contig": "chr1",
            "strand": strand if species == "A" else "+", "annotation_start": 1,
            "annotation_end": len(locus_sequence), "search_start": 1,
            "search_end": len(locus_sequence), "genome_fasta": str(fasta)})
        if species == "B":
            occurrence = {"occurrence_id": "B1", "family_id": "fam", "species": "B",
                "gene_copy_id": gene, "transcript_id": "txB", "role": "exon",
                "presence_status": "present", "contig": "chr1", "start": 1,
                "end": len(coding), "strand": "+"}
            occurrences.append(occurrence)
            paths.append({**occurrence, "path_id": "B:txB:B1", "cds_intervals": f"1-{len(coding)}",
                "cds_length": len(coding), "cds_phase": "0"})
            features.append({"family_id": "fam", "species": "B", "gene_copy_id": gene,
                "id": "txB", "parent": gene, "ownership": "target_gene_descendant",
                "type": "mRNA", "seqid": "chr1", "start": 1, "end": len(coding),
                "strand": "+", "phase": ".", "attrs": "ID=txB;Parent=gB"})
            features.append({"family_id": "fam", "species": "B", "gene_copy_id": gene,
                "id": "cdsB", "parent": "txB", "ownership": "target_gene_descendant",
                "type": "CDS", "seqid": "chr1", "start": 1, "end": len(coding),
                "strand": "+", "phase": "0", "attrs": "ID=cdsB;Parent=txB"})
            fasta_lines.append(("B1", coding))
            continue
        left_seq, right_seq = coding[:split], coding[split:]
        gap_left, gap_right = split + 1, split + 5
        if strand == "-":
            left_start, left_end = len(locus_sequence) - split + 1, len(locus_sequence)
            right_start, right_end = 1, len(locus_sequence) - split - 5
        else:
            left_start, left_end = 1, split
            right_start, right_end = split + 6, len(locus_sequence)
        phase = str((3 - (split % 3)) % 3) if phase_known else "."
        for alias in range(aliases):
            tx = f"txA{alias}"
            features.append({"family_id": "fam", "species": "A", "gene_copy_id": gene,
                "id": tx, "parent": gene, "ownership": "target_gene_descendant",
                "type": "mRNA", "seqid": "chr1", "start": 1, "end": len(locus_sequence),
                "strand": strand, "phase": ".", "attrs": f"ID={tx};Parent={gene}"})
            for suffix, start, end, seq, part_phase in (("L", left_start, left_end, left_seq, "0"),
                                                          ("R", right_start, right_end, right_seq, phase)):
                oid = f"A{alias}{suffix}"
                occurrence = {"occurrence_id": oid, "family_id": "fam", "species": "A",
                    "gene_copy_id": gene, "transcript_id": tx, "role": "exon",
                    "presence_status": "present", "contig": "chr1", "start": start,
                    "end": end, "strand": strand}
                occurrences.append(occurrence)
                paths.append({**occurrence, "path_id": f"A:{tx}:{oid}",
                    "cds_intervals": f"{start}-{end}", "cds_length": len(seq),
                    "cds_phase": part_phase})
                features.append({"family_id": "fam", "species": "A", "gene_copy_id": gene,
                    "id": f"cds{oid}", "parent": tx, "ownership": "target_gene_descendant",
                    "type": "CDS", "seqid": "chr1", "start": start, "end": end,
                    "strand": strand, "phase": part_phase, "attrs": f"ID=cds{oid};Parent={tx}"})
                fasta_lines.append((oid, seq))
            intron_start, intron_end = (right_end + 1, left_start - 1) if strand == "-" else (gap_left, gap_right)
            if annotated_intron:
                introns.append({"family_id": "fam", "species": "A", "gene_copy_id": gene,
                    "contig": "chr1", "strand": strand, "start": intron_start,
                    "end": intron_end, "intron_id": f"intron{alias}", "transcript_id": tx,
                    "left_feature_id": f"cds{tx}L", "right_feature_id": f"cds{tx}R",
                    "right_phase": phase if phase_known else "."})
    if gene_only:
        targets.append({"family_id": "fam", "species": "C", "gene_copy_id": "gC"})
    write_tsv(inputs / "segment_occurrences.tsv", occurrences, list(occurrences[0]))
    path_fields = list(paths[0]) if paths else ["occurrence_id", "family_id", "species", "gene_copy_id",
        "transcript_id", "path_id", "cds_intervals", "cds_length", "cds_phase"]
    write_tsv(inputs / "transcript_paths.tsv", [] if no_paths else paths, path_fields)
    write_tsv(inputs / "raw_gene_features.tsv", features, list(features[0]))
    write_tsv(inputs / "gene_loci.tsv", loci, list(loci[0]))
    intron_fields = ["family_id", "species", "gene_copy_id", "contig", "strand", "start", "end",
        "intron_id", "transcript_id", "left_feature_id", "right_feature_id", "right_phase"]
    write_tsv(inputs / "intron_sites.tsv", introns, intron_fields)
    write_tsv(inputs / "input_targets.tsv", targets, list(targets[0]))
    (inputs / "segment_sequences.fasta").write_text("".join(f">{key}\n{sequence}\n" for key, sequence in fasta_lines))
    return inputs


def competing_transcript_fixture(root, *, same_physical=False):
    inputs = root / "inputs"
    inputs.mkdir()
    dna_a, dna_b = "GCT" * 62, "GCT" * 60
    occurrences, paths, features, loci, introns, targets, fasta_records = [], [], [], [], [], [], []
    for species, gene, dna in (("A", "gA", dna_a), ("B", "gB", dna_b)):
        fasta = root / f"{species}.fa"
        fasta.write_text(f">chr1\n{dna}\n")
        targets.append({"family_id": "fam", "species": species, "gene_copy_id": gene})
        loci.append({"species": species, "gene_copy_id": gene, "contig": "chr1", "strand": "+",
            "annotation_start": 1, "annotation_end": len(dna), "search_start": 1,
            "search_end": len(dna), "genome_fasta": str(fasta)})
        oid, default_tx = f"{species}1", f"tx{species}"
        occurrence = {"occurrence_id": oid, "family_id": "fam", "species": species,
            "gene_copy_id": gene, "transcript_id": default_tx, "role": "exon",
            "presence_status": "present", "contig": "chr1", "start": 1,
            "end": len(dna), "strand": "+"}
        occurrences.append(occurrence)
        fasta_records.append((oid, dna))
        if species == "B":
            paths.append({**occurrence, "path_id": "B:txB", "cds_intervals": "1-180",
                "cds_length": 180, "cds_phase": "0"})
            features.append({"family_id": "fam", "species": "B", "gene_copy_id": gene,
                "id": "txB", "parent": gene, "ownership": "target_gene_descendant",
                "type": "mRNA", "seqid": "chr1", "start": 1, "end": 180,
                "strand": "+", "phase": ".", "attrs": "ID=txB;Parent=gB"})
            features.append({"family_id": "fam", "species": "B", "gene_copy_id": gene,
                "id": "cdsB", "parent": "txB", "ownership": "target_gene_descendant",
                "type": "CDS", "seqid": "chr1", "start": 1, "end": 180,
                "strand": "+", "phase": "0", "attrs": "ID=cdsB;Parent=txB"})
            continue
        for ordinal, intervals in enumerate(("1-60;64-183",
                "4-60;64-186" if same_physical else "1-60;67-186")):
            tx = f"txA{ordinal}"
            paths.append({**occurrence, "transcript_id": tx, "path_id": f"A:{tx}",
                "cds_intervals": intervals, "cds_length": 180, "cds_phase": "0"})
            parts = [tuple(map(int, part.split("-"))) for part in intervals.split(";")]
            features.append({"family_id": "fam", "species": "A", "gene_copy_id": gene,
                "id": tx, "parent": gene, "ownership": "target_gene_descendant",
                "type": "mRNA", "seqid": "chr1", "start": 1, "end": len(dna),
                "strand": "+", "phase": ".", "attrs": f"ID={tx};Parent={gene}"})
            for part_i, (start, end) in enumerate(parts):
                features.append({"family_id": "fam", "species": "A", "gene_copy_id": gene,
                    "id": f"cds{ordinal}_{part_i}", "parent": tx,
                    "ownership": "target_gene_descendant", "type": "CDS", "seqid": "chr1",
                    "start": start, "end": end, "strand": "+", "phase": "0",
                    "attrs": f"ID=cds{ordinal}_{part_i};Parent={tx}"})
            first_end, second_start = parts[0][1], parts[1][0]
            introns.append({"family_id": "fam", "species": "A", "gene_copy_id": gene,
                "contig": "chr1", "strand": "+", "start": first_end + 1,
                "end": second_start - 1, "intron_id": f"intron{ordinal}",
                "transcript_id": tx, "left_feature_id": "NA", "right_feature_id": "NA",
                "right_phase": "0"})
    write_tsv(inputs / "segment_occurrences.tsv", occurrences, list(occurrences[0]))
    write_tsv(inputs / "transcript_paths.tsv", paths, list(paths[0]))
    write_tsv(inputs / "raw_gene_features.tsv", features, list(features[0]))
    write_tsv(inputs / "gene_loci.tsv", loci, list(loci[0]))
    write_tsv(inputs / "input_targets.tsv", targets, list(targets[0]))
    write_tsv(inputs / "intron_sites.tsv", introns, ["family_id", "species", "gene_copy_id", "contig",
        "strand", "start", "end", "intron_id", "transcript_id", "left_feature_id",
        "right_feature_id", "right_phase"])
    (inputs / "segment_sequences.fasta").write_text(
        "".join(f">{record_id}\n{sequence}\n" for record_id, sequence in fasta_records))
    return inputs
