"""Reference values for the slab-waveguide mode-matching comparison.

The numerical values below are the analytical references supplied in the project
brief and are intentionally left as fixed validation data. They should not be
modified to force agreement with FDTD output.
"""

from __future__ import annotations

from config import ANALYTICAL_EFFECTIVE_INDICES, MODE_MATCHING_REFERENCE, WIDTH_INPUT, WIDTH_OUTPUT


def analytical_mode_table(width: float) -> dict:
    """Return the analytical effective indices for a given guide width."""
    return dict(ANALYTICAL_EFFECTIVE_INDICES[width])


if __name__ == "__main__":
    print("Input guide analytical TE effective indices:")
    for mode, neff in analytical_mode_table(WIDTH_INPUT).items():
        print(f"  {mode}: {neff:.6f}")
    print("\nOutput guide analytical TE effective indices:")
    for mode, neff in analytical_mode_table(WIDTH_OUTPUT).items():
        print(f"  {mode}: {neff:.6f}")
    print("\nGuided-mode truncation reference for 0.90 -> 0.60 um step:")
    for key, value in MODE_MATCHING_REFERENCE.items():
        print(f"  {key}: {value:.6f}")
