"""Width-step FDTD simulation for a 0.90 to 0.60 um dielectric slab.

The input guide fills the region y < 0 and the output guide fills y > 0,
with the discontinuity at z = 0. The launched source is the fundamental TE0
mode of the input guide. Fields are decomposed into guided modes of the input
and output slab at monitors placed away from the interface and away from the
PML.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import meep as mp
import numpy as np

from config import (
    CELL_WIDTH,
    CELL_LENGTH,
    FREQUENCY,
    N_CLAD,
    N_CORE,
    PML_THICKNESS,
    RESULTS_DIR,
    RESOLUTION,
    SOURCE_Y,
    WIDTH_INPUT,
    WIDTH_OUTPUT,
    WAVELENGTH,
)


def solve_input_mode(width: float = WIDTH_INPUT, band_num: int = 1):
    """Solve the guided TE mode in a uniform slab of the given width."""
    sim = mp.Simulation(
        cell_size=mp.Vector3(CELL_WIDTH, CELL_LENGTH, 0),
        geometry=[
            mp.Block(
                size=mp.Vector3(width, mp.inf, mp.inf),
                center=mp.Vector3(0, 0, 0),
                material=mp.Medium(index=N_CORE),
            )
        ],
        default_material=mp.Medium(index=N_CLAD),
        boundary_layers=[mp.PML(PML_THICKNESS)],
        resolution=RESOLUTION,
        dimensions=2,
    )
    sim.init_sim()
    mode = sim.get_eigenmode(
        frequency=FREQUENCY,
        where=mp.Volume(center=mp.Vector3(0, 0, 0), size=mp.Vector3(CELL_WIDTH, 0, 0)),
        band_num=band_num,
        direction=mp.Y,
        kpoint=mp.Vector3(0, 1, 0),
    )
    sim.reset_meep()
    return mode


def build_step_geometry():
    input_block = mp.Block(
        size=mp.Vector3(WIDTH_INPUT, 6.0, 0),
        center=mp.Vector3(0, -3.0, 0),
        material=mp.Medium(index=N_CORE),
    )
    output_block = mp.Block(
        size=mp.Vector3(WIDTH_OUTPUT, 6.0, 0),
        center=mp.Vector3(0, 3.0, 0),
        material=mp.Medium(index=N_CORE),
    )
    return [input_block, output_block]


def mode_power_dictionary(alpha: np.ndarray, bands: list[int]) -> dict[str, float]:
    """Map Meep band numbers to forward/backward coefficients.

    We deliberately avoid assuming that Meep band number equals a TE_m label.
    The physical mode ordering must be confirmed from the actual transverse
    field profile, not from the band index alone.
    """
    out: dict[str, float] = {}
    for idx, band in enumerate(bands):
        out[f"band_{band}"] = {
            "forward": float(np.abs(alpha[idx, 0, 0]) ** 2),
            "backward": float(np.abs(alpha[idx, 0, 1]) ** 2),
        }
    return out


def run_step_simulation(resolution: int = RESOLUTION, save_prefix: str = "width_step") -> dict:
    """Run the width-step simulation and return guided-mode powers."""
    sim = mp.Simulation(
        cell_size=mp.Vector3(CELL_WIDTH, CELL_LENGTH, 0),
        geometry=build_step_geometry(),
        sources=[
            mp.EigenModeSource(
                src=mp.GaussianSource(frequency=FREQUENCY, fwidth=0.1 * FREQUENCY),
                center=mp.Vector3(0, SOURCE_Y, 0),
                size=mp.Vector3(WIDTH_INPUT, 0, 0),
                eig_band=1,
                eig_kpoint=mp.Vector3(0, 1, 0),
                eig_match_freq=True,
                direction=mp.Y,
            )
        ],
        default_material=mp.Medium(index=N_CLAD),
        boundary_layers=[mp.PML(PML_THICKNESS)],
        resolution=resolution,
        dimensions=2,
    )

    input_monitor = sim.add_mode_monitor(
        FREQUENCY,
        0,
        1,
        mp.ModeRegion(center=mp.Vector3(0, -1.75, 0), size=mp.Vector3(CELL_WIDTH - 2 * PML_THICKNESS, 0, 0)),
    )
    output_monitor = sim.add_mode_monitor(
        FREQUENCY,
        0,
        1,
        mp.ModeRegion(center=mp.Vector3(0, 2.25, 0), size=mp.Vector3(CELL_WIDTH - 2 * PML_THICKNESS, 0, 0)),
    )

    print(f"Running width-step simulation at {resolution} pixels/um...")
    sim.run(until_after_sources=200)

    bands = [1, 2, 3, 4]
    input_coeff = sim.get_eigenmode_coefficients(
        input_monitor,
        bands=bands,
        direction=mp.NO_DIRECTION,
        kpoint_func=lambda f, n: mp.Vector3(0, 1, 0),
    )
    output_coeff = sim.get_eigenmode_coefficients(
        output_monitor,
        bands=bands,
        direction=mp.NO_DIRECTION,
        kpoint_func=lambda f, n: mp.Vector3(0, 1, 0),
    )

    alpha_in = input_coeff.alpha[:, 0, :]
    alpha_out = output_coeff.alpha[:, 0, :]
    incident = float(np.abs(alpha_in[0, 0]) ** 2)

    if incident <= 0:
        raise RuntimeError("Incident guided-mode power is zero; the source/monitor normalization failed.")

    per_mode: dict[str, dict[str, float]] = {}
    reflected_total = 0.0
    transmitted_total = 0.0
    for idx, band in enumerate(bands):
        label = f"band_{band}"
        coeff_ref = float(np.abs(alpha_in[idx, 1]) ** 2) / incident
        coeff_trans = float(np.abs(alpha_out[idx, 0]) ** 2) / incident
        per_mode[label] = {"reflected": coeff_ref, "transmitted": coeff_trans}
        reflected_total += coeff_ref
        transmitted_total += coeff_trans

    guided_total = reflected_total + transmitted_total
    unaccounted = max(0.0, 1.0 - guided_total)

    # Save raw coefficients for inspection.
    csv_path = RESULTS_DIR / f"{save_prefix}_modal_coefficients.csv"
    with csv_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "band",
            "mode_label",
            "reflected_coeff_real",
            "reflected_coeff_imag",
            "transmitted_coeff_real",
            "transmitted_coeff_imag",
            "reflected_power",
            "transmitted_power",
            "incident_power",
        ])
        for band in bands:
            label = f"band_{band}"
            coeff_ref = alpha_in[band - 1, 1]
            coeff_trans = alpha_out[band - 1, 0]
            writer.writerow([
                band,
                label,
                np.real(coeff_ref),
                np.imag(coeff_ref),
                np.real(coeff_trans),
                np.imag(coeff_trans),
                per_mode[label]["reflected"],
                per_mode[label]["transmitted"],
                incident,
            ])

    field = np.asarray(sim.get_array(component=mp.Ez, center=mp.Vector3(), size=mp.Vector3(CELL_WIDTH, CELL_LENGTH, 0)))
    np.save(RESULTS_DIR / f"{save_prefix}_field.npy", field)

    plt.figure(figsize=(7, 8))
    plt.imshow(
        field.T,
        origin="lower",
        extent=[-CELL_WIDTH / 2, CELL_WIDTH / 2, -CELL_LENGTH / 2, CELL_LENGTH / 2],
        aspect="auto",
        cmap="RdBu",
    )
    plt.xlabel("Transverse coordinate x (um)")
    plt.ylabel("Propagation coordinate z (um)")
    plt.title("Uniform Slab Width Step: FDTD TE Mode")
    plt.colorbar(label="Ez field (arbitrary units)")
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / f"{save_prefix}_field.png", dpi=200)
    plt.close()

    result = {
        "resolution": resolution,
        "incident_power": incident,
        "guided_reflected_total": reflected_total,
        "guided_transmitted_total": transmitted_total,
        "guiding_total": guided_total,
        "unaccounted_power": unaccounted,
        "per_mode": per_mode,
    }
    print("\nMode power summary (band order from Meep, not TE_m labels)")
    for label, values in per_mode.items():
        print(f"  {label}: reflected={values['reflected']:.6f}, transmitted={values['transmitted']:.6f}")
    print(f"  guided reflected total={reflected_total:.6f}")
    print(f"  guided transmitted total={transmitted_total:.6f}")
    print(f"  unaccounted power={unaccounted:.6f}")
    return result


if __name__ == "__main__":
    result = run_step_simulation()
    print(f"\nSaved modal-power CSV to {RESULTS_DIR / 'width_step_modal_coefficients.csv'}")
    print(f"Saved field map to {RESULTS_DIR / 'width_step_field.png'}")
