# Monte Carlo Photon Transport and DPF Research

MATLAB research code for Monte Carlo photon transport, near-infrared diffuse optics, differential pathlength factor (DPF) modeling, breast and tissue phantom studies, and localized absorber concentration recovery.

> **Repository status:** Research archive under active curation. The numbered files preserve the chronological development of the methods. Some studies require external simulation packages or datasets that are not distributed here.

## Research themes

| Directory | Scope |
|---|---|
| `01_monte_carlo_foundations` | Early photon-packet transport, Beer–Lambert studies, and MCX/MCmatlab validation |
| `02_optical_properties_and_kramers_kronig` | Optical-property datasets and Kramers–Kronig studies |
| `03_breast_tumor_geometry` | Breast/tumor geometry, depth, wavelength, and boundary studies |
| `04_compression_beam_and_dpf` | Compression, source geometry, two-dimensional measurements, and DPF analysis |
| `05_digital_breast_phantoms` | Voxelized breast phantoms and tumor/hematocrit simulations |
| `06_laboratory_phantom_models` | TiO2, India ink, agar/gelatin, and experimental phantom models |
| `07_dpf_parameter_studies` | Absorption/scattering sweeps, DPF models, thickness, gradients, and sensitivity |
| `08_oxygenation_and_proposal_studies` | Oxygenation and blood-volume studies |
| `09_concentration_recovery` | Localized absorbing inclusions and concentration-recovery simulations |
| `10_transport_models` | Advection–diffusion transport models |
| `11_validation_and_calibration` | Calibration and cross-dataset validation studies |

See [`docs/CODE_CATALOG.md`](docs/CODE_CATALOG.md) for the complete file-level catalog and [`docs/VERSIONING.md`](docs/VERSIONING.md) for the naming convention.

## Requirements

- MATLAB (recent versions recommended)
- [MCmatlab](https://github.com/ankrh/MCmatlab) for most Monte Carlo simulations
- MCXLAB/MCX for the early MCX example (`Photon_03.m`)
- MATLAB toolboxes depend on the individual study; plotting, statistics, and curve fitting functions may require additional products
- External experimental and simulation datasets for analysis files that call `load`, `readtable`, or import mesh files

Third-party packages and external research data are **not vendored** in this repository. Install them separately and confirm their licenses before use.

## Setup

1. Clone this repository.
2. Copy `config/project_config.example.m` to `config/project_config.m`.
3. Set the local paths to MCmatlab and the external data directory.
4. In MATLAB, run:

   ```matlab
   run('config/project_config.m')
   addpath(genpath('code'))
   ```

5. Open the relevant cluster and consult its version family in the catalog before running a file.

Many historical versions still contain the original machine-specific paths. They are flagged by `tools/check_repository.sh` and retained to preserve the exact research record. Replace those paths with values from `project_config.m` when promoting a workflow to a reproducible release.

## Data and results policy

Large generated `.mat` files, raw measurements, meshes, and most figures are excluded from version control. Put local input data under `data/` and generated outputs under `results/`. Only small, redistribution-safe examples should be committed.

## Reproducibility status

The repository distinguishes preservation from reproducibility:

- **Historical versions** preserve the original numerical work and its evolution.
- **Analysis variants** often require output from the matching primary simulation.
- **Reproducible workflows** should eventually receive a dedicated driver, example data, fixed configuration, and expected-output test.

The current curation preserves all 194 supplied MATLAB files. No scientific equations or parameter values were silently changed.

## Citation and license

Academic use should cite the associated publication when a script corresponds to published work. A software license has intentionally not been selected yet; see [`LICENSE-DECISION.md`](LICENSE-DECISION.md). Until a license is added, normal copyright applies and reuse requires permission.

## Author

Kaiser Niknam — computational electromagnetics, biomedical optics, and numerical modeling.
