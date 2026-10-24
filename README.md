# Mode Matching and FDTD Analysis of Optical Waveguide Discontinuities

A computational study of electromagnetic scattering at an abrupt width discontinuity in an optical slab waveguide, using **analytical mode matching** and **finite-difference time-domain (FDTD)** simulations with [Meep](https://meep.readthedocs.io/en/latest/).

The project investigates how an incident guided mode interacts with a change in waveguide width, with particular attention to modal reflection, transmission, higher-order mode coupling, and numerical convergence.

## Overview

Optical waveguide discontinuities can cause reflections and redistribute optical power among the guided modes supported by the structure. Understanding this behaviour is important when analysing integrated photonic devices and transitions between waveguides of different dimensions.

This project studies a two-dimensional, lossless optical slab waveguide with an abrupt transition from a wider input guide to a narrower output guide. Two complementary approaches are used:

- **Mode matching:** Calculate modal reflection and transmission coefficients by matching electromagnetic fields at the discontinuity.
- **FDTD simulation:** Numerically simulate electromagnetic field propagation through the structure using Meep.

The resulting modal coefficients and power distributions are compared to investigate the agreement between the analytical and numerical approaches.

## Physical Structure

The simulated structure consists of a symmetric dielectric slab waveguide with a width step at the interface.

| Parameter | Value |
|---|---:|
| Free-space wavelength | 1.55 µm |
| Core refractive index | 3.40 |
| Cladding refractive index | 1.45 |
| Input waveguide width | 0.90 µm |
| Output waveguide width | 0.60 µm |
| Incident mode | Fundamental TE₀ mode |
| Material model | Lossless, nondispersive |
| Simulation method | 2D FDTD |
| FDTD solver | Meep |

The fundamental TE₀ mode is incident from the wider waveguide toward the narrower waveguide. The discontinuity can produce reflected fields and redistribute transmitted power among the modes supported by the output guide.

## Methods

### 1. Analytical Mode Matching

The electromagnetic field on either side of the discontinuity is represented using the guided modes of the corresponding uniform waveguide sections.

The fields are matched at the interface to obtain a system of equations for the unknown modal amplitudes. Solving this system yields the reflection and transmission coefficients associated with the incident mode and the supported output and reflected modes.

The analysis includes:

- Calculation of modal propagation constants and effective refractive indices.
- Construction of modal overlap and matching equations.
- Solution for reflected and transmitted modal amplitudes.
- Evaluation of guided-mode reflection and transmission contributions.

### 2. Finite-Difference Time-Domain Simulation

Meep is used to simulate electromagnetic field propagation through the width-step geometry.

The simulation workflow includes:

- Defining the dielectric geometry and material properties.
- Exciting the fundamental guided mode in the input waveguide.
- Applying absorbing perfectly matched layer (PML) boundaries.
- Recording electromagnetic fields and evaluating modal or power quantities.
- Repeating simulations at different spatial resolutions to assess convergence.

### 3. Comparison and Numerical Convergence

The analytical and FDTD approaches are compared in terms of modal scattering and reflected and transmitted power.

A spatial-resolution study is used to investigate the sensitivity of the numerical solution to discretization. Differences between the two methods are examined with consideration of modal normalization, finite computational domains, absorbing boundaries, monitor placement, and radiation or other power contributions not captured by a guided-mode-only calculation.

## Repository Structure

```text
waveguide-mode-matching-fdtd/
├── config.py
├── mode_matching.py
├── modal_analysis.py
├── fdtd_step.py
├── fdtd_identity.py
├── convergence.py
├── plots.py
├── results/
├── requirements.txt
└── README.md
```

- `config.py` — Shared physical and simulation parameters.
- `mode_matching.py` — Analytical mode-matching implementation.
- `modal_analysis.py` — Waveguide modal analysis.
- `fdtd_step.py` — FDTD simulation of the width discontinuity.
- `fdtd_identity.py` — Reference simulation for checking the uniform-waveguide setup.
- `convergence.py` — Spatial-resolution and numerical convergence study.
- `plots.py` — Plotting and visualization utilities.
- `results/` — Selected simulation outputs, numerical data, and generated figures.

*The file descriptions above reflect the intended roles of the project scripts; update them if the implementation differs.*

## Software and Dependencies

The project uses:

- **Python** — Numerical implementation and analysis.
- **NumPy** — Array operations and matrix construction.
- **SciPy** — Linear algebra and numerical computation.
- **Meep** — Finite-difference time-domain electromagnetic simulation.
- **Matplotlib** — Field visualization and plotting.

Meep must be installed in a compatible Python environment. Refer to the [official Meep installation guide](https://meep.readthedocs.io/en/latest/Installation/) for platform-specific instructions.

## Reproducing the Simulations

Clone the repository:

```bash
git clone https://github.com/YOUR-USERNAME/waveguide-mode-matching-fdtd.git
cd waveguide-mode-matching-fdtd
```

Activate the Python environment in which Meep and the required numerical libraries are installed.

If the dependencies are listed in `requirements.txt`, install them using:

```bash
python -m pip install -r requirements.txt
```

Run the scripts individually according to the desired analysis:

```bash
python mode_matching.py
python modal_analysis.py
python fdtd_step.py
python convergence.py
```

The commands above assume the scripts are independently executable from the repository root. Adjust them if a script requires additional arguments or a specific execution order.

## Results

This section will document the computed results and their interpretation as the analytical and numerical comparisons are finalized.

Planned outputs include:

- Modal effective-index and propagation-constant data.
- Electromagnetic field distributions across the width discontinuity.
- Modal reflection and transmission coefficients.
- Guided reflected and transmitted power.
- Spatial-resolution convergence plots and numerical error analysis.
- Comparison between mode-matching predictions and FDTD results.

Numerical values and figures should be added here after their normalization, monitor dependence, and convergence have been checked.

## Scope and Limitations

The current study considers a two-dimensional, lossless dielectric slab waveguide with an abrupt width change. The analysis is focused on guided-mode scattering and comparison between analytical mode matching and numerical FDTD simulation.

Agreement between the two approaches depends on consistent field and power normalization, correct modal identification, and accounting for power outside the guided modes. Convergence of the FDTD solution and the interpretation of higher-order modal contributions are therefore important parts of the analysis.

The results should be interpreted within these modelling assumptions and the limitations of the implemented numerical methods.

## Acknowledgements

This project was developed as part of an academic study in photonics and electromagnetic wave propagation.

## License

No license has been assigned yet. Add an appropriate license if you intend to permit reuse or redistribution of the code.
