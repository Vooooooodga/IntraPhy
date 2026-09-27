"""Explicit joint-repertoire command-line surface."""
from __future__ import annotations

def add_repertoire_commands(sub):
    for name, help_text in (("infer-exon-repertoires", "Evaluate an explicit joint exon-repertoire CTMC."),
                            ("fit-exon-repertoire-rates", "Fit one explicit shared repertoire-rate scale.")):
        command = sub.add_parser(name, help=help_text)
        command.add_argument("--exon-configurations", required=True)
        command.add_argument("--species-tree", required=True)
        command.add_argument("--repertoire-model", required=True)
        command.add_argument("--output-dir", required=True)
        if name == "infer-exon-repertoires":
            command.add_argument("--expected-edits", action="store_true")
        else:
            command.add_argument("--log-scale-bounds", nargs=2, type=float, required=True)
            command.add_argument("--collection-provenance", required=True)
