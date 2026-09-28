"""Public DNA-material presence CTMC API."""
from __future__ import annotations

from intraphy.inference.binary_presence import DNA_DOMAIN
from intraphy.inference.binary_presence import analyze_binary_presence


def analyze_dna_presence(
    rows, tree, output_dir, *, root_frequency="stationary", root_presence=0.5,
    fixed_rates=None, branch_length_unit="supplied_tree_units", expected_counts=False,
    workers=1,
):
    """Fit homologous DNA-material presence histories with family-shared rates."""
    return analyze_binary_presence(
        rows, tree, output_dir, domain=DNA_DOMAIN,
        root_frequency=root_frequency, root_presence=root_presence,
        fixed_rates=fixed_rates, branch_length_unit=branch_length_unit,
        expected_counts=expected_counts, workers=workers,
    )
