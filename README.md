# Tile_18O_Isotopes

Data and analysis code for:

> Arfania, H., Kayler, Z. E., Strawn, D. G., Brooks, E. S., & Laan, M. (2026).
> **Flow-Dependent Phosphate Oxygen-Isotope (δ¹⁸O-PO₄) Signatures in Tile Drainage from Two Long-Term Tillage Systems.**
> *Journal of Geophysical Research: Biogeosciences* (manuscript 2026JG010070, in revision).

Continuous tile-drain discharge, water temperature, dissolved reactive phosphorus, dual water isotopes,
and phosphate oxygen isotopes from two adjacent long-term tillage systems (conventional tillage, CT, and
no-till, NT) at the R.J. Cook Agronomy Farm LTAR site near Pullman, Washington, USA, September 2022 – July 2024.

**Study design note.** One field was monitored per management system. The comparison is between two
long-term systems, not replicated treatments; see the paper for the interpretation this permits.

## Contents

```
data/                      raw monitoring data (unchanged from field/lab exports)
R/                         analysis pipeline as six modules + runner
CAF_tile_analysis.R        the same pipeline as a single file
output/tables/             Tables 1–7 and Data Sets S1–S2 as CSV; analysis_log.txt
output/figures/            Figures 1–5, S1–S3 (PNG, 200 dpi, review copies)
output/figures_hires/      the same figures as 600-dpi TIFF and vector PDF (journal files)
supporting_information/    Supporting Information PDF and Data Sets S1–S2
```

## Data files (`data/`)

| File | Contents | Notes |
|---|---|---|
| `Tile_d18Op.csv` | δ¹⁸O-PO₄ of tile-drain phosphate, ‰ VSMOW | 27 dates × CT/NT = 54 values |
| `18Ow_Tile.csv` | δ¹⁸O and δ²H of tile-drain water, ‰ VSMOW | 25 dates × CT/NT = 50 values; first column (date) is unlabelled in the raw export |
| `18Ow_Prcep.csv` | δ¹⁸O and δ²H of precipitation | 22 samples; 12 Apr 2023 is a d-excess outlier excluded from the LMWL |
| `Pullman_Climate_Daily.csv` | Daily station precipitation, Pullman WA, 1 Jan 2023 – 31 Dec 2024 | used for Figure 1; begins after discharge monitoring started |
| `Prceipitation.csv` | Earlier wet-day-only precipitation file | superseded by the station record; retained for provenance |
| `CT_Discharge_2023.csv`, `NT_Discharge_2023.csv` | 15-min tile discharge, L s⁻¹, Sep 2022 – Aug 2023 | NT file is a pivot-table export containing `(blank)` and `Grand Total` rows that must be dropped |
| `CT_Discharge_2024.csv`, `NT_Discharge_2024.csv` | Sub-daily tile discharge, **ft³ s⁻¹**, Nov 2023 – Jul 2024 | date-only timestamps; NT file has a 4-line TOA5 header |
| `*_water_tempreture_*.csv` | Sub-daily tile-water temperature, °C | sentinel −6999 present; NT 2024 column is headed `SR_Temp_C` |
| `DRP_2022_23.csv`, `DRP_2023_24.csv` | Dissolved reactive phosphorus, mg L⁻¹ | CT and NT stored side-by-side in separate column blocks |

Every one of these quirks is handled explicitly in `R/01_ingest_and_clean.R`, which logs what it removed.

## Reproducing the analysis

Requires R (≥ 4.3; the paper used 4.5.1) and `nlme`, which ships with R. No other packages.

```r
source("CAF_tile_analysis.R")        # or: source("R/00_run_all.R")
```

Run from the repository root, or set `PROJ` to the repository path first. Everything is written to
`output_revision/`; runtime is a few seconds. Type III sums of squares, HC3 robust standard errors and the
Breusch–Pagan test are implemented directly in the code so that no contributed packages are needed.

The pipeline regenerates every number, table and figure in the paper. `output/tables/analysis_log.txt`
is the log from the run used for the submitted manuscript, for cross-checking.

## Citing

Please cite the paper, and cite this archive as:

> Arfania, H., Kayler, Z. E., Strawn, D. G., Brooks, E. S., & Laan, M. (2026). *Tile_18O_Isotopes: data and
> code for flow-dependent phosphate oxygen-isotope signatures in tile drainage* (v1.0.0) [Data set and software].
> Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX

(The DOI is minted when the first GitHub release is archived by Zenodo — see `RELEASE_STEPS.md`.)

## Licence

Data: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Code: [MIT](LICENSE).
