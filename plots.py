"""Create comparison plots for the slab-waveguide results."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt

from config import MODE_MATCHING_REFERENCE, RESULTS_DIR


def load_csv(path: Path) -> list[dict]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def make_step_comparison_plot() -> None:
    path = RESULTS_DIR / "width_step_modal_coefficients.csv"
    data = load_csv(path)
    if not data:
        raise FileNotFoundError(f"No width-step modal results found at {path}")

    # Keep the comparison honest: the current Meep basis is band-ordered, not
    # TE_m-labeled, so we compare only the dominant guided band against the
    # analytical TE0 reference and the total guided power against the mode-matching
    # guided totals.
    ref_te0 = MODE_MATCHING_REFERENCE["transmitted_TE0"]
    ref_total = MODE_MATCHING_REFERENCE["guided_transmitted_total"]

    band_1 = next((float(row["transmitted_power"]) for row in data if row["mode_label"] == "band_1"), 0.0)
    band_1_ref = next((float(row["reflected_power"]) for row in data if row["mode_label"] == "band_1"), 0.0)

    labels = ["TE0 transmitted", "guided transmitted total"]
    refp = [ref_te0, ref_total]
    fdtp = [band_1, sum(float(row["transmitted_power"]) for row in data)]

    fig, ax = plt.subplots(figsize=(8, 5))
    x = range(len(labels))
    ax.bar(x, refp, width=0.35, label="mode-matching reference")
    ax.bar([i + 0.35 for i in x], fdtp, width=0.35, label="FDTD")
    ax.set_xticks([i + 0.175 for i in x])
    ax.set_xticklabels(labels, rotation=18, ha="right")
    ax.set_ylabel("Power")
    ax.set_title("Width-step guided power comparison")
    ax.legend()
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / "width_step_modal_power_comparison.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    make_step_comparison_plot()
    print(f"Saved comparison plot to {RESULTS_DIR / 'width_step_modal_power_comparison.png'}")
