# Data Dictionary

This document describes each input data file in the `data/` directory.

**Site:** R. J. Cook Agronomy Farm (CAF), USDA-ARS Long-Term Agroecosystem Research (LTAR) network, Palouse, Pullman, WA, USA.
**Treatments:** CT = conventional tillage; NT = no-till.
**Period:** Water years 2022–2024 (Oct 2022 – Sep 2024, varying by stream).
**Isotope notation:** δ values reported in per mil (‰) relative to VSMOW for water and PO₄.

---

## Hydrology

### `CT_Discharge_2023.csv`, `CT_Discharge_2024.csv`, `NT_Discharge_2023.csv`, `NT_Discharge_2024.csv`

Continuous tile-drain discharge logged at 15-minute intervals. **Note:** the four files use different headers (full instrument header, two-column "Date, cfs", row-labels layout, and Campbell Scientific TOA5). The loader in `R/01_main_analysis.R` (`read_discharge()`) auto-detects timestamp and discharge columns.

| Column (after import) | Units | Description |
|---|---|---|
| `treatment` | — | "CT" or "NT" |
| `datetime` | POSIXct | 15-min timestamp |
| `discharge_L_s` | L s⁻¹ | Tile discharge (converted to L s⁻¹ when source is cfs: × 28.317) |

### `CT_water_tempreture_2023.csv`, `CT_water_tempreture_2024.csv`, `NT_water_tempreture_2023.csv`, `NT_water_tempreture_2024.csv`

Continuous tile-drain water temperature, 15-minute intervals.

| Column | Units | Description |
|---|---|---|
| `Timestamp` | POSIXct | 15-min timestamp |
| `Temp_C` | °C | Water temperature; QC filter applied (−2 °C ≤ T ≤ 30 °C; sentinel value −6999 removed) |

---

## Water chemistry

### `DRP_2022_23.csv` (paired CT/NT layout)

| Column | Units | Description |
|---|---|---|
| `Sample_ID` (cols 1, 8) | — | Sample identifier (`CAF_<trt>_TL`) |
| `Date` (cols 2, 9) | date | Grab sample date |
| `DRP_mg/L` (cols 3, 10) | mg L⁻¹ | Dissolved reactive phosphorus |

The file is paired side-by-side: columns 1–3 are CT, columns 8–10 are NT.

### `DRP_2023_24.csv` (stacked layout)

| Column | Units | Description |
|---|---|---|
| `SampleID` | — | Sequential ID |
| `Location` | — | Sampling location code: `CAF_<trt>_TL` (tile, used) or `CAF_<trt>_SR` (surface runoff, **filtered out** in analysis) |
| `Date` | date | Grab sample date |
| `Time` | HHMM | Sample time |
| `mg/L DRP` | mg L⁻¹ | Dissolved reactive phosphorus |
| `pH` | — | Field pH |
| `EC (mS)` | mS cm⁻¹ | Electrical conductivity |
| `Turbidity (FTU)` | FTU | Turbidity |

**Method detection limit (DRP):** 0.005 mg L⁻¹. Values below MDL are set to MDL before flux calculations.

---

## Stable isotopes

### `18Ow_Prcep.csv` (precipitation)

| Column | Units | Description |
|---|---|---|
| `Date` | date | Precipitation event collection date |
| `18O_H2O` (loaded as `d18O_H2O`) | ‰ VSMOW | δ¹⁸O of precipitation H₂O |
| `Detrium_H2O` (loaded as `d2H_H2O`) | ‰ VSMOW | δ²H of precipitation H₂O |

Outlier exclusion: for the boxplot and main analysis, samples with |δ¹⁸O − mean| > 5 SD are removed. For the Craig diagram, the additional d-excess criterion |d − mean| > 3 SD is applied.

### `18Ow_Tile.csv` (tile-drain dual isotopes)

Side-by-side CT / NT layout.

| Column | Units | Description |
|---|---|---|
| (first column) | date | Sampling date |
| `NT_18OW`, `CT_18OW` | ‰ VSMOW | δ¹⁸O–H₂O for NT, CT |
| `NT_Deutrium`, `CT_Deutrium` | ‰ VSMOW | δ²H–H₂O for NT, CT |

### `Tile_d18Op.csv` (phosphate oxygen isotopes)

| Column | Units | Description |
|---|---|---|
| `Date` | date | Grab sample date |
| `CT` | ‰ VSMOW | δ¹⁸O–PO₄ for the CT tile |
| `NT` | ‰ VSMOW | δ¹⁸O–PO₄ for the NT tile |

---

## Auxiliary

### `NT_P_Fraction.csv`

Reference Hedley sequential P fractionation data, top-slope vs. toe-slope of the NT field. Used in the manuscript for context only — not consumed by the main analysis script.

| Column | Units | Description |
|---|---|---|
| `Date` | date | Sampling date |
| `DI water` | mg kg⁻¹ | Loosely bound P (DI water extract) |
| `NaHCO3` | mg kg⁻¹ | Labile P (NaHCO₃ extract) |
| `NaOH` | mg kg⁻¹ | Fe/Al-bound P (NaOH extract) |
| `NHO3` (sic) | mg kg⁻¹ | Ca-bound P (HNO₃ extract) |

Top-slope and toe-slope blocks are separated within the file.

---

## Notes on quality control

- All date parsers accept multiple formats (mdY, Ymd, dmY, with optional time).
- Discharge data are imported via a single robust loader; unit conversion from cfs to L s⁻¹ is applied where the source file is in cfs.
- Tile water temperature observations outside [−2, 30] °C and sentinel values (e.g., −6999) are excluded.
- Sample matching for the dynamic equilibrium calculation uses nearest-date within ±7 days for water isotopes and temperature and ±3 days for discharge. Lag days are recorded in `output/po4_matched_dynamic.csv`.
