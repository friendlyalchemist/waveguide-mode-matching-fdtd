"""Convergence sweep for the width-step simulation."""

from __future__ import annotations

import csv

from config import RESULTS_DIR
from fdtd_step import run_step_simulation


def convergence_sweep(resolutions=(20, 30, 40)) -> list[dict]:
    rows: list[dict] = []
    for resolution in resolutions:
        result = run_step_simulation(resolution=resolution, save_prefix=f"step_res_{resolution}")
        rows.append({
            "resolution": resolution,
            "guided_reflected_total": result["guided_reflected_total"],
            "guided_transmitted_total": result["guided_transmitted_total"],
            "unaccounted_power": result["unaccounted_power"],
            "band_1_reflected": result["per_mode"].get("band_1", {}).get("reflected", 0.0),
            "band_1_transmitted": result["per_mode"].get("band_1", {}).get("transmitted", 0.0),
            "band_2_reflected": result["per_mode"].get("band_2", {}).get("reflected", 0.0),
            "band_2_transmitted": result["per_mode"].get("band_2", {}).get("transmitted", 0.0),
        })
    return rows


def save_convergence(rows: list[dict]) -> None:
    path = RESULTS_DIR / "convergence.csv"
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "resolution",
                "guided_reflected_total",
                "guided_transmitted_total",
                "unaccounted_power",
                "band_1_reflected",
                "band_1_transmitted",
                "band_2_reflected",
                "band_2_transmitted",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    rows = convergence_sweep()
    save_convergence(rows)
    print(f"Convergence results saved to {RESULTS_DIR / 'convergence.csv'}")
    for row in rows:
        print(row)
