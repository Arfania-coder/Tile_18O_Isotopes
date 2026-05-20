# Tile_18O_Isotopes

**Data and R code accompanying:**

Arfania, H., Kayler, Z. E., Strawn, D. G., Brooks, E. S., & Laan, M. (in review). *Tillage × Flow Regime Interaction Controls Phosphorus Source Dynamics in Tile Drainage: Integrating Continuous Discharge Monitoring, Water Residence Time, and Phosphate Oxygen Isotopes.* Submitted to *Journal of Geophysical Research: Biogeosciences*.

**Study site:** R. J. Cook Agronomy Farm (CAF), USDA-ARS Long-Term Agroecosystem Research (LTAR) network, Palouse, Pullman, WA, USA.

---

## Overview

This repository contains all raw input data and R scripts required to reproduce every table, figure, and statistic reported in the manuscript and supplementary material. The pipeline ingests continuous discharge and water-temperature loggers, grab-sample chemistry, and dual water isotopes (δ¹⁸O–H₂O, δ²H–H₂O) and δ¹⁸O–PO₄, and produces:

- Discharge summary statistics and Q75 flow-regime thresholds
- Dissolved reactive phosphorus (DRP) summaries by treatment × flow regime
- Local Meteoric Water Line (LMWL) and tile-water evaporation line (Craig diagram with inset)
- Mean transit time / young water fraction (damping ratio and amplitude ratio methods)
- **Dynamic, sample-specific equilibrium δ¹⁸O–PO₄** using the Chang & Blake (2015) equation, with sensitivity comparison against the fixed-equilibrium assumption
- Two-way ANOVA on observed δ¹⁸O–PO₄ and on Δ_eq (dynamic and fixed)

---

## Repository structure
```
Tile_18O_Isotopes/
├── R/
│   ├── 00_run_all.R           # master runner — sources both scripts in order
│   ├── 01_main_analysis.R     # full pipeline: Tables 1–8, Figs 2–4, 6–7, S1–S2
│   └── 02_craig_diagram.R     # Craig diagram with magnified inset
├── data/                      # raw input CSVs (see Data dictionary below)
├── figures/                   # PNG (300 dpi), PDF, TIFF (600 dpi LZW) outputs
├── output/                    # CSV tables + analysis_summary.rds
├── docs/
│   └── data_dictionary.md
├── CITATION.cff
├── LICENSE
├── .gitignore
└── README.md
```

## Reproducing the analysis

### Requirements

- **R** ≥ 4.2.0
- The following CRAN packages (installed automatically on first run): `here`, `tidyverse`, `lubridate`, `ggpubr`, `car`, `effsize`, `minpack.lm`, `scales`, `cowplot`, `viridis`.

### Run

1. Clone the repository: `git clone https://github.com/Arfania-coder/Tile_18O_Isotopes.git`
2. Open the project from the repository root in RStudio (or set R's working directory to the repo root).
3. In the R console:
```r
   source("R/00_run_all.R")
```

All paths are resolved with the `here` package, so the scripts work identically on Windows, macOS, and Linux without modification. Outputs are written to `output/` and `figures/`.

To run only one part:
```r
source("R/01_main_analysis.R")   # everything except the Craig diagram
source("R/02_craig_diagram.R")   # Craig diagram with inset (standalone)
```

---

## Data dictionary

| File | Description |
|---|---|
| `18Ow_Prcep.csv` | Precipitation δ¹⁸O–H₂O and δ²H–H₂O (per mil VSMOW), event-basis. |
| `18Ow_Tile.csv` | Tile-drain δ¹⁸O–H₂O and δ²H–H₂O by treatment (CT, NT), grab samples. |
| `Tile_d18Op.csv` | Tile-drain δ¹⁸O–PO₄ (per mil VSMOW) by treatment, grab samples. |
| `DRP_2022_23.csv` | Dissolved reactive P, water year 2022–23 (paired CT/NT layout). |
| `DRP_2023_24.csv` | Dissolved reactive P, water year 2023–24 (stacked; tile samples only retained). |
| `CT_Discharge_2023.csv` / `_2024.csv` | Conventional-till tile discharge, 15-min logger. |
| `NT_Discharge_2023.csv` / `_2024.csv` | No-till tile discharge, 15-min logger. |
| `CT_water_tempreture_2023.csv` / `_2024.csv` | CT tile water temperature, 15-min logger. |
| `NT_water_tempreture_2023.csv` / `_2024.csv` | NT tile water temperature, 15-min logger. |
| `NT_P_Fraction.csv` | Hedley sequential P fractions (top vs. toe slope), reference data. |

Detailed column definitions, units, sampling design, and quality-control flags are documented in [`docs/data_dictionary.md`](docs/data_dictionary.md).

**Treatment codes:** CT = conventional tillage, NT = no-till.
**Sampling site code:** *_TL = tile-drain outlet (used for analysis); *_SR = surface-runoff samples (filtered out, not used in this paper).
**MDL for DRP:** 0.005 mg L⁻¹ (values below MDL set to MDL prior to flux calculations).
**Isotope reporting:** all δ values in per mil relative to VSMOW.

---

## Key methods and equations

### Chang & Blake (2015) equilibrium equation

δ¹⁸O–PO₄(eq) = exp[(14.43 / T_K) − 0.02654] × (δ¹⁸O–H₂O + 1000) − 1000

In this revision (script v2.0), δ¹⁸O–PO₄(eq) is computed *per sample* using the treatment-specific tile-drain δ¹⁸O–H₂O closest in time to each grab and the daily-mean tile water temperature derived from the 15-min logger record. Fixed-equilibrium (δ¹⁸O–H₂O = −14.1 ‰, T = 5 °C → 11.2 ‰) results are retained as a sensitivity check.

### Flow regime classification

Daily-mean tile discharge ≥ treatment-specific Q75 → "Event"; otherwise → "Baseflow".

### Mean transit time

- **Damping ratio:** τ = 1 / (ω × √((1/DR²) − 1)), where DR = σ_tile / σ_precip and ω = 2π / 365.25.
- **Amplitude ratio:** sinusoidal fit to δ¹⁸O–H₂O time series in precipitation and tile drain; Fyw (young water fraction) = A_tile / A_precip.

Fyw is reported as a tracer-derived index, not a literal residence-time estimate (Kirchner, 2016).

---

## Outputs

Running the full pipeline produces:

| File | Contents |
|---|---|
| `output/Table1_Discharge_Statistics.csv` | Continuous discharge summary by treatment |
| `output/Table2_DRP_Summary.csv` | DRP means, flow-weighted means, and loads by treatment × flow regime |
| `output/Table3_Water_Isotope_Summary.csv` | δ¹⁸O–H₂O summaries (precip vs. tile) |
| `output/Table4_Residence_Time.csv` | Mean transit time and Fyw |
| `output/Table5_PO4_Overall_Summary.csv` | δ¹⁸O–PO₄ summary by treatment |
| `output/Table6_PO4_by_Flow_Regime.csv` | δ¹⁸O–PO₄ and Δ_eq by treatment × flow regime |
| `output/Table7_ANOVA_Results.csv` | Two-way ANOVA, three response specifications |
| `output/Table8_Equilibrium_Parameters.csv` | Equilibrium calculation inputs and outputs |
| `output/po4_matched_dynamic.csv` | Each PO4 sample with matched H₂O, T, Q, and Δ_eq |
| `output/treatment_eq_reference.csv` | Treatment-mean dynamic equilibrium |
| `output/analysis_summary.rds` | Structured R list of all key results |
| `figures/Fig2_…` through `Fig7_…`, `FigS1_…`, `FigS2_…`, `Fig_Craig_…` | All figures in PNG/PDF/TIFF |

---

## Citation

If you use this code or data, please cite both the manuscript (citation will be updated upon acceptance) and this repository (see `CITATION.cff`).

## License

- **Code:** MIT License (see `LICENSE`).
- **Data:** Creative Commons Attribution 4.0 International (CC BY 4.0).
