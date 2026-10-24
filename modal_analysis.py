
"""Modal analysis for the Orange slab-waveguide FDTD project.

The Meep eigenmode solver returns the dominant propagation constant in the
wavevector units used by MPB. For this TE slab, the effective index is

    n_eff = kdom.y * wavelength

because kdom.y is in units of 2*pi / um and the wavelength is in um.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import meep as mp
import numpy as np

from config import (
    ANALYTICAL_EFFECTIVE_INDICES,
    CELL_WIDTH,
    FREQUENCY,
    N_CLAD,
    N_CORE,
    PML_THICKNESS,
    PROJECT_DIR,
    RESULTS_DIR,
    RESOLUTION,
    WAVELENGTH,
    WIDTH_INPUT,
    WIDTH_OUTPUT,
)


def slab_modes_for_width(width: float, bands: list[int] | tuple[int, ...] = (1, 2, 3, 4)) -> list[dict]:
    """Return the solved guided modes for a uniform slab of a given width."""
    sim = mp.Simulation(
        cell_size=mp.Vector3(CELL_WIDTH, 12.0, 0),
        geometry=[
            mp.Block(
                size=mp.Vector3(width, mp.inf, mp.inf),
                center=mp.Vector3(0, 0, 0),
                material=mp.Medium(index=N_CORE),
            )
        ],
        boundary_layers=[mp.PML(PML_THICKNESS)],
        default_material=mp.Medium(index=N_CLAD),
        resolution=RESOLUTION,
        dimensions=2,
    )
    sim.init_sim()

    rows: list[dict] = []
    mode_region = mp.Volume(center=mp.Vector3(0, 0, 0), size=mp.Vector3(CELL_WIDTH, 0, 0))
    for band in bands:
        try:
            mode = sim.get_eigenmode(
                frequency=FREQUENCY,
                where=mode_region,
                band_num=band,
                direction=mp.Y,
                kpoint=mp.Vector3(0, 1, 0),
            )
        except Exception:
            continue

        kdom_y = float(mode.kdom.y)
        neff = kdom_y * WAVELENGTH
        is_guided = neff > N_CLAD + 1e-8
        rows.append(
            {
                "width": width,
                "band": band,
                "freq": float(mode.freq),
                "kdom_y": kdom_y,
                "n_eff": neff,
                "group_velocity": float(mode.group_velocity),
                "is_guided": is_guided,
                "reference": ANALYTICAL_EFFECTIVE_INDICES[width].get(f"TE{band-1}", None),
            }
        )

    sim.reset_meep()
    return rows


def save_mode_csv(rows: list[dict], filename: str) -> None:
    """Write a CSV of the solved modal results."""
    path = RESULTS_DIR / filename
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["width", "band", "freq", "kdom_y", "n_eff", "group_velocity", "is_guided", "reference"],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def mode_parity_and_confinement(profile: np.ndarray, xs: np.ndarray, width: float) -> tuple[str, float, float]:
    """Estimate parity, node count, and core confinement from the transverse field profile."""
    core_mask = np.abs(xs) <= width / 2
    core_power = float(np.sum(np.abs(profile[core_mask]) ** 2))
    total_power = float(np.sum(np.abs(profile) ** 2))
    if total_power <= 0:
        parity = "undetermined"
        node_count = 0
        confinement = 0.0
    else:
        node_count = int(np.sum(np.diff(np.sign(profile[core_mask])) != 0))
        center_val = float(profile[np.argmin(np.abs(xs))])
        parity = "even" if abs(center_val) > 0.1 * np.max(np.abs(profile)) else "odd"
        confinement = core_power / total_power
    return parity, float(node_count), float(confinement)


def plot_mode_profiles(width: float, rows: list[dict]) -> None:
    """Plot the transverse field profile for the guided modes of a slab."""
    sim = mp.Simulation(
        cell_size=mp.Vector3(CELL_WIDTH, 12.0, 0),
        geometry=[
            mp.Block(
                size=mp.Vector3(width, mp.inf, mp.inf),
                center=mp.Vector3(0, 0, 0),
                material=mp.Medium(index=N_CORE),
            )
        ],
        boundary_layers=[mp.PML(PML_THICKNESS)],
        default_material=mp.Medium(index=N_CLAD),
        resolution=RESOLUTION,
        dimensions=2,
    )
    sim.init_sim()

    xs = np.linspace(-CELL_WIDTH / 2, CELL_WIDTH / 2, 1001)
    fig, ax = plt.subplots(figsize=(8, 6))
    for row in rows:
        if not row["is_guided"]:
            continue
        band = row["band"]
        mode = sim.get_eigenmode(
            frequency=FREQUENCY,
            where=mp.Volume(center=mp.Vector3(0, 0, 0), size=mp.Vector3(CELL_WIDTH, 0, 0)),
            band_num=band,
            direction=mp.Y,
            kpoint=mp.Vector3(0, 1, 0),
        )
        profile = np.array([
            np.real(mode.amplitude(mp.Vector3(x, 0, 0), mp.Ez)) for x in xs
        ])
        parity, node_count, confinement = mode_parity_and_confinement(profile, xs, width)
        ax.plot(xs, profile, label=f"band {band}: n_eff={row['n_eff']:.4f}, parity={parity}, nodes={int(node_count)}, core={confinement:.3f}")

    ax.axvline(-width / 2, color="k", linestyle="--", linewidth=1, alpha=0.7)
    ax.axvline(width / 2, color="k", linestyle="--", linewidth=1, alpha=0.7)
    ax.text(-width / 2 - 0.08, 0.9 * ax.get_ylim()[1], "core edge", rotation=90, va="top", ha="center", fontsize=8)
    ax.text(width / 2 + 0.08, 0.9 * ax.get_ylim()[1], "core edge", rotation=90, va="top", ha="center", fontsize=8)
    ax.set_title(f"Uniform slab mode profiles (width={width:.2f} um)")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("Re{Ez}")
    ax.grid(alpha=0.2)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(RESULTS_DIR / f"modal_profile_width_{width:.2f}um.png", dpi=200)
    plt.close(fig)
    sim.reset_meep()


def save_mode_validation_csv(rows: list[dict], filename: str) -> None:
    """Save a validation table with parity, node count, and core confinement."""
    path = RESULTS_DIR / filename
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["width", "band", "n_eff", "kdom_y", "group_velocity", "parity", "node_count", "core_confinement", "is_guided"])
        for row in rows:
            if not row["is_guided"]:
                continue
            # Re-solve the mode to extract the profile-based parity/confinement data.
            sim = mp.Simulation(
                cell_size=mp.Vector3(CELL_WIDTH, 12.0, 0),
                geometry=[
                    mp.Block(
                        size=mp.Vector3(row["width"], mp.inf, mp.inf),
                        center=mp.Vector3(0, 0, 0),
                        material=mp.Medium(index=N_CORE),
                    )
                ],
                boundary_layers=[mp.PML(PML_THICKNESS)],
                default_material=mp.Medium(index=N_CLAD),
                resolution=RESOLUTION,
                dimensions=2,
            )
            sim.init_sim()
            mode = sim.get_eigenmode(
                frequency=FREQUENCY,
                where=mp.Volume(center=mp.Vector3(0, 0, 0), size=mp.Vector3(CELL_WIDTH, 0, 0)),
                band_num=row["band"],
                direction=mp.Y,
                kpoint=mp.Vector3(0, 1, 0),
            )
            xs = np.linspace(-CELL_WIDTH / 2, CELL_WIDTH / 2, 1001)
            profile = np.array([np.real(mode.amplitude(mp.Vector3(x, 0, 0), mp.Ez)) for x in xs])
            parity, node_count, core_confinement = mode_parity_and_confinement(profile, xs, float(row["width"]))
            sim.reset_meep()
            writer.writerow([
                row["width"],
                row["band"],
                row["n_eff"],
                row["kdom_y"],
                row["group_velocity"],
                parity,
                node_count,
                core_confinement,
                row["is_guided"],
            ])


def main() -> None:
    widths = [WIDTH_INPUT, WIDTH_OUTPUT]
    all_rows: list[dict] = []
    for width in widths:
        rows = slab_modes_for_width(width)
        all_rows.extend(rows)
        plot_mode_profiles(width, rows)
        save_mode_validation_csv(rows, f"mode_validation_{width:.2f}um.csv")
    save_mode_csv(all_rows, "modal_results.csv")

    print("\nModal results (validated effective-index extraction):")
    for row in all_rows:
        print(
            f"width={row['width']:.2f} um, band={row['band']}, "
            f"n_eff={row['n_eff']:.6f}, group_velocity={row['group_velocity']:.6f}, "
            f"guided={row['is_guided']}"
        )
    print(f"\nCSV written to: {RESULTS_DIR / 'modal_results.csv'}")


if __name__ == "__main__":
    main()
