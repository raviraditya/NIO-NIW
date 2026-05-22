# NIO-NIW

MATLAB analysis code for tropical cyclone near-inertial wave energetics in the North Indian Ocean.

This repository accompanies the manuscript:

Singh, R., & Behera, M. R. (2025). *Tropical Cyclone Generated Near-Inertial Energy and Its Vertical Redistribution in the North Indian Ocean*. JGR: Oceans (submitted).

## Contents

- `Main_combined.m` — Core diagnostic engine. Runs per-TC extraction of HYCOM ocean state, near-inertial velocity decomposition, NIKE budget, WKB normalisation, vertical group velocity, wind power input, resonance parameter, and right/left TC-relative partition. Set `STORM_ID` at the top to one of `KYARR`, `AMPHAN`, `FANI`, `TAUKTAE`.
- `Post_combined.m` — Post-processing. Reads the .mat/.xlsx output of `Main_combined.m` and computes penetration metrics, conversion efficiency, Richardson-number statistics, R/L asymmetry classification, and Langmuir-turbulence flags.
- `HYCOM_Argo_validation.m` — Validates the HYCOM upper-ocean state against quality-controlled Argo profiles co-located within 50 km and 3 h of each TC track. Produces the 5 m T/S biases, 50–500 m correlations, and MLD comparisons reported in Text S1, Figure S2, and Table S1 of the Supporting Information.

## Requirements

- MATLAB R2021a or newer (uses `datetime`, `readtable`, `ncread`)
- [GSW Gibbs Sea-Water Oceanographic Toolbox for MATLAB](https://www.teos-10.org/software.htm) on the path

## Data layout expected

Place input data under a `./data/` folder structured as:

```
data/
  IMD/          # cyclone best-track CSVs (Kyarr.csv, Amphan.csv, Fani.csv, Tauktae.csv)
  HYCOM/        # per-TC subdirs of 3-hourly HYCOM NetCDF snapshots
  Stress/       # ERA5 surface stress NetCDF
  Wind/         # ERA5 10-m wind NetCDF
  Currents/     # ocean-current NetCDF (if used)
  INCOIS_OHC/   # ocean heat-content NetCDF (optional)
  Argo/         # Argo delayed-mode profile NetCDFs (for validation only)
```

Sources: HYCOM GLBy0.08 (https://www.hycom.org/dataserver/gofs-3pt1/reanalysis), ERA5 (Copernicus Climate Data Store), IMD cyclone e-Atlas best tracks, Argo Global Data Assembly Centre (https://argo.ucsd.edu/).

## How to run

1. Edit `STORM_ID` at the top of `Main_combined.m`, then run it in MATLAB.
2. After `Main_combined.m` produces `*_V16_SUPERCHARGED.mat` and `.xlsx`, run `Post_combined.m`.
3. To reproduce the HYCOM–Argo validation (Supporting Information), run `HYCOM_Argo_validation.m` independently — it does not depend on the output of the main pipeline.

## License

MIT (see LICENSE).

## Citation

If you use this code, please cite the Zenodo archive (DOI listed under Releases) and the manuscript above.
