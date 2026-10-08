"""Resolve upstream FASTA member identifiers through explicit GFF ancestry."""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from Bio import SeqIO

from .orthofinder_ids import build_gff_locus_index, member_token_groups
from ..storage.tabular import open_text
from .resources import SpeciesFiles


@dataclass(frozen=True)
class Target:
    family: str
    species: str
    gene: str
    members: tuple[str, ...]
    source: str


def read_headers(path: Path) -> tuple[str, ...]:
    """Sequences identify an input set; they are not genomic absence evidence."""
    headers: list[str] = []
    identifiers: set[str] = set()
    with open_text(path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if not len(record.seq):
                raise ValueError(f"Empty ortholog FASTA record: {record.id} in {path}")
            if record.id in identifiers:
                raise ValueError(f"Duplicate FASTA identifier {record.id} in {path}; use Species|ID")
            identifiers.add(record.id)
            headers.append(record.description)
    if not headers:
        raise ValueError(f"Ortholog FASTA contains no records: {path}")
    return tuple(headers)


def _candidate_headers(header: str, species: set[str]) -> tuple[str | None, str]:
    first, *rest = header.split(maxsplit=1)
    prefix, separator, member = first.partition("|")
    # Many accession formats contain '|'. Only an exact known species prefix
    # has special meaning; arbitrary prefixes are never guessed.
    if separator and prefix in species:
        if not member:
            raise ValueError(f"Missing member identifier in FASTA header: {header}")
        return prefix, member + (" " + rest[0] if rest else "")
    return None, header


def targets_from_fastas(resources: tuple[SpeciesFiles, ...],
                        families: dict[str, Path]) -> tuple[Target, ...]:
    """Use exact ID/alias/metadata matches; a hit must resolve to one locus.

    Indices are processed one species at a time, so whole-genome annotations
    from all species are not retained simultaneously. No sequence search or
    gene-symbol fallback is used to assign upstream ortholog membership.
    """
    species = {item.species for item in resources}
    records = [(family, path, header, *_candidate_headers(header, species))
               for family, path in families.items() for header in read_headers(path)]
    hits: dict[int, set[tuple[str, str]]] = defaultdict(set)
    for resource in resources:
        exact, _, _ = build_gff_locus_index(resource.gff)
        for index, (_, _, _, declared_species, header) in enumerate(records):
            if declared_species and declared_species != resource.species:
                continue
            primary, metadata = member_token_groups(header)
            for tokens in (primary, metadata):
                matched = set().union(*(exact.get(token, set()) for token in tokens))
                if matched:
                    hits[index].update((resource.species, locus) for locus in matched)
                    break
    grouped: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    for index, (family, path, header, _, _) in enumerate(records):
        candidates = hits[index]
        if len(candidates) != 1:
            detail = "unresolved" if not candidates else f"ambiguous: {sorted(candidates)}"
            raise ValueError(
                f"Ortholog member {header!r} in {path} is {detail}. "
                "Use an exact GFF gene/transcript/protein ID, or Species|ID; "
                "AGAT gene= metadata is accepted. Orthology is not inferred."
            )
        taxon, gene = next(iter(candidates))
        grouped[family, taxon, gene].add(header)
    result = [Target(family, taxon, gene, tuple(sorted(members)), str(families[family]))
              for (family, taxon, gene), members in sorted(grouped.items())]
    counts: dict[tuple[str, str], int] = defaultdict(int)
    for target in result:
        counts[target.family, target.species] += 1
    for family in families:
        invalid = {taxon: counts[family, taxon] for taxon in sorted(species)
                   if counts[family, taxon] != 1}
        if invalid:
            raise ValueError(f"Family {family} requires one gene locus per species; counts={invalid}. "
                             "Multiple isoforms of one locus are allowed; missing taxa and paralogs are not.")
    return tuple(result)
