"""Explicit read-only reuse of recorded genomic alignment evidence."""
from __future__ import annotations

import json
from pathlib import Path

from ..storage.fasta import parse_fasta


def _read_command(directory, label, executable, expected, inputs):
    path = directory / f"{label}.command.json"
    stdout = directory / f"{label}.stdout"
    if not path.is_file() or not stdout.is_file():
        raise ValueError(f"--alignment-evidence-dir lacks {label} command record or stdout.")
    try:
        record = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"--alignment-evidence-dir has an invalid {label} command record.") from exc
    command = record.get("command") if isinstance(record, dict) else None
    if (not isinstance(record, dict) or record.get("status") != "completed"
            or record.get("returncode") != 0
            or record.get("shell") is not False or not isinstance(command, list)
            or not command or any(not isinstance(item, str) for item in command)
            or Path(command[0]).name != executable):
        raise ValueError(f"--alignment-evidence-dir {label} command did not complete successfully.")
    args = command[1:]
    recorded_inputs = args[-len(inputs):]
    compatible = (len(args) == len(expected) + len(inputs)
        and all(wanted is None or actual == wanted
                for actual, wanted in zip(args[:len(expected)], expected))
        and all(Path(actual).is_absolute() and Path(actual).name == item.name
                for actual, item in zip(recorded_inputs, inputs)))
    if compatible:
        thread_values = [args[index] for index, item in enumerate(expected) if item is None]
        compatible = all(value.isdigit() and int(value) > 0 for value in thread_values)
    if not compatible:
        raise ValueError(f"--alignment-evidence-dir {label} command used incompatible inputs or algorithm options.")
    return record, stdout, recorded_inputs


def load_alignment_evidence(directory, records, queries, identifier_species):
    """Load prior MAFFT/minimap2 output only after direct input/provenance checks."""
    directory = Path(directory)
    if not directory.is_dir():
        raise ValueError(f"--alignment-evidence-dir family evidence is missing: {directory}")
    raw_path, query_path = directory / "genomic_loci.fa", directory / "exon_queries.fa"
    species_alignment = directory / "genomic_alignment.fa"
    for path in (raw_path, query_path, species_alignment):
        if not path.is_file():
            raise ValueError(f"--alignment-evidence-dir is missing required evidence file: {path.name}")
    try:
        saved_records = parse_fasta(raw_path)
        saved_queries = parse_fasta(query_path)
    except (OSError, ValueError) as exc:
        raise ValueError("--alignment-evidence-dir contains invalid saved FASTA inputs.") from exc
    if saved_records != records:
        raise ValueError("--alignment-evidence-dir genomic_loci.fa differs from current prepared loci.")
    if saved_queries != queries:
        raise ValueError("--alignment-evidence-dir exon_queries.fa differs from current exon queries.")

    mafft_record, mafft_stdout, mafft_inputs = _read_command(directory, "genomic_mafft", "mafft",
        ["--auto", "--inputorder", "--thread", None, "--threadit", "0"],
        [raw_path])
    minimap_record, minimap_stdout, minimap_inputs = _read_command(directory, "exon_minimap2", "minimap2",
        ["-c", "-x", "asm20", "--secondary=yes", "-N", "50", "-t", None],
        [raw_path, query_path])
    if (mafft_inputs[0] != minimap_inputs[0]
            or Path(minimap_inputs[0]).parent != Path(minimap_inputs[1]).parent):
        raise ValueError("--alignment-evidence-dir recorded commands do not share the same historical family inputs.")

    aligned = parse_fasta(mafft_stdout)
    if set(aligned) != set(records) or len({len(sequence) for sequence in aligned.values()}) != 1:
        raise ValueError("--alignment-evidence-dir saved MAFFT output changed record identifiers or alignment lengths.")
    for identifier, sequence in aligned.items():
        if sequence.replace("-", "").upper() != records[identifier].upper():
            raise ValueError("--alignment-evidence-dir saved MAFFT output changed genomic residues.")
    by_species = parse_fasta(species_alignment)
    expected_species = {identifier_species[key]: sequence.upper()
                        for key, sequence in aligned.items()}
    if by_species != expected_species:
        raise ValueError("--alignment-evidence-dir species-keyed MSA differs from saved MAFFT records.")
    paf = minimap_stdout.read_text(encoding="utf-8")
    query_ids, target_ids = set(queries), set(records)
    for line in paf.splitlines():
        fields = line.split("\t")
        if len(fields) >= 12 and (fields[0] not in query_ids or fields[5] not in target_ids):
            raise ValueError("--alignment-evidence-dir minimap2 output has unexpected query/target identifiers.")
    return aligned, paf, {"mafft": mafft_record, "minimap2": minimap_record}
