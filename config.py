
"""
Shared configuration for the Orange slab-waveguide FDTD project.
All lengths are in micrometres (um).
"""

from pathlib import Path

# --------------------------------------------------
# Output directories
# --------------------------------------------------
PROJECT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = PROJECT_DIR / "results"
FIELDS_DIR = RESULTS_DIR / "fields"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIELDS_DIR.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------
# Physical parameters
# --------------------------------------------------
WAVELENGTH = 1.55          # Free-space wavelength (um)
FREQUENCY = 1.0 / WAVELENGTH  # Meep frequency units (1/um)

N_CORE = 3.40
N_CLAD = 1.45

WIDTH_INPUT = 0.90         # Input waveguide width (um)
WIDTH_OUTPUT = 0.60        # Output waveguide width (um)

# --------------------------------------------------
# Computational parameters
# --------------------------------------------------
RESOLUTION = 20            # Grid points per um; default FDTD run value
PML_THICKNESS = 1.0        # PML thickness (um)

CELL_WIDTH = 4.0           # Transverse x dimension (um)
CELL_LENGTH = 12.0         # Propagation y dimension (um)

# Meep uses x-y coordinates here:
# x = physical transverse coordinate
# y = physical propagation coordinate z
# Ez = physical electric field Ey for TE polarization.

# --------------------------------------------------
# Source and monitor positions
# --------------------------------------------------
SOURCE_Y = -3.0
REFLECTION_MONITOR_Y = -1.5
TRANSMISSION_MONITOR_Y = 2.5

# --------------------------------------------------
# Simulation controls
# --------------------------------------------------
RUN_TIME = 50.0

# Analytical reference values supplied for validation.
# These are not FDTD outputs and should not be overwritten.
ANALYTICAL_EFFECTIVE_INDICES = {
    WIDTH_INPUT: {
        "TE0": 3.320758,
        "TE1": 3.074103,
        "TE2": 2.628060,
        "TE3": 1.911362,
    },
    WIDTH_OUTPUT: {
        "TE0": 3.244942,
        "TE1": 2.748945,
        "TE2": 1.808935,
    },
}

# The guided-mode truncation reference for the 0.90-to-0.60 um interface.
MODE_MATCHING_REFERENCE = {
    "guided_reflected_total": 0.002756,
    "guided_transmitted_total": 0.997244,
    "reflected_TE0": 0.000448,
    "reflected_TE2": 0.002308,
    "transmitted_TE0": 0.940719,
    "transmitted_TE2": 0.056525,
}
