"""Public CTMC API for homologous intron-position presence histories."""
from __future__ import annotations

from intraphy.inference.binary_presence import INTRON_DOMAIN
from intraphy.inference.binary_presence import analyze_binary_presence


def analyze_intron_positions(
    rows, tree, output_dir, *, root_frequency="stationary", root_presence=0.5,
    fixed_rates=None, branch_length_unit="supplied_tree_units", expected_counts=False,
    workers=1,
):
    """Fit intron-position presence histories from explicitly typed observations.

    Every row must declare ``observation_type='intron_position_presence'``.
    ``absent`` and ``present`` encode coding continuity or intron interruption
    at homologous coding positions; annotation supplies the evidence for those
    states. Unknown or inapplicable observations remain unobserved.
    """
    encoded = []
    for index, original in enumerate(rows):
        row = dict(original)
        if row.get("observation_type") != INTRON_DOMAIN["observation_type"]:
            raise ValueError(
                f"row {index + 1} observation_type must be {INTRON_DOMAIN['observation_type']!r}"
            )
        if row.get("layer") not in (None, "", "intron_position"):
            raise ValueError(f"row {index + 1} layer conflicts with intron-position observation_type")
        state = row.get("state")
        if isinstance(state, bool):
            raise ValueError(f"invalid intron-position presence state in row {index + 1}: {state!r}")
        if state in (0, "0", "absent"):
            row["state"] = 0
        elif state in (1, "1", "present"):
            row["state"] = 1
        elif state == "unknown":
            row["state"] = "unknown"
        else:
            raise ValueError(f"invalid intron-position presence state in row {index + 1}: {state!r}")
        row["layer"] = "intron_position"
        encoded.append(row)
    return analyze_binary_presence(
        encoded, tree, output_dir, domain=INTRON_DOMAIN,
        root_frequency=root_frequency, root_presence=root_presence,
        fixed_rates=fixed_rates, branch_length_unit=branch_length_unit,
        expected_counts=expected_counts, workers=workers,
    )
