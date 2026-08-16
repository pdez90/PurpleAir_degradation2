# PurpleAir sensor degradation analysis

Code for: *Long-term degradation, survivorship, and environmental drivers of instability in the global PurpleAir PM2.5 sensor network* (submitted to *Environmental Science & Technology*).

We analyze 452,143,148 hourly observations (2018–2025) from 36,328 outdoor PurpleAir sensors across CONUS, Africa, India, and South America to characterize sensor degradation, survivorship, and their environmental drivers.

## Main analysis scripts

Run in order — the numbering is a genuine dependency chain, not just a suggestion: `00_build_analysis_dataset.R` writes the analytic dataset every later script loads, and `02_climate_zone_assignment.R` writes the time-varying Köppen file that the climate-stratified models in `03`–`08` read. Each script otherwise loads its own inputs and libraries and can be run on its own once its upstream outputs exist.

| Script | Produces |
|---|---|
| `00_config.R` | Paths, GEE/GCS identifiers, key analysis constants — **edit this first** |
| `00_build_analysis_dataset.R` | Rebuilds and **verifies** the analytic dataset: **452,143,148 rows / 36,328 sensors** (applies the cf_1 > 1500 µg/m³ exclusion and all T/RH filters) |
| `01_descriptive_statistics.R` | Sensor movement/density (Figure S1); Tables S1–S2 and sample numbers; PA-vs-reference correlations; deployment timelines (Figures S4–S25) |
| `02_climate_zone_assignment.R` | Time-varying Köppen–Geiger join and climate-zone figures (Figures S2–S3) |
| `03_primary_age_degradation_models.R` | Primary age-dependent cf_1 degradation models — global, regional, climate-zone and region × climate → **Figure 2; Tables S3–S6** |
| `04_model_reference_divergence.R` | Age-dependent divergence vs GEOS-CF PM2.5 and ERA5-Land T/RH → **Figure 3; Tables S7–S9; Figures S36–S38** |
| `05_survivorship_attrition.R` | Kaplan–Meier, RMST, landmark Cox (10,000 h), pre-offline change, and effective lifetime → **Figure 4; Tables S10–S18; Figures S39–S43** |
| `06_high_pm_exposure.R` | High-PM exposure → future degradation (7/14/28-day windows × 35/55/100/150 µg/m³) → **Figure 5; Tables S19–S24; Figures S44–S49** |
| `07_spatial_analysis.R` | Nearby-pair similarity, Moran's I, semivariograms, spatial GAMs → **Section S4** |
| `08_persistent_degradation.R` | Stringent persistent / permanent degradation classification → **Section S5** (the **1,146 persistent / 206 permanent** result). Uses the sensor-specific **90th**-percentile relative-divergence flag rather than the 85th used in the primary analysis — documented in SI section S5.1 |
| `09_pa_epa_validation.R` | Join to EPA AQS monitors within 100 m, calibration models and calibration drift → **Section S3; Table S3.1; Figures S3.1–S3.4** |

Regional and climate-zone heterogeneity is not a separate script: those models are estimated alongside the global models inside `03_primary_age_degradation_models.R` and `04_model_reference_divergence.R`, as in the original code.

## `sensitivity/` — robustness and extension analyses

Run after the main scripts. Each is self-contained and requires `outdoor_final.RData` to exist (run `00_build_analysis_dataset.R` first).

| Script | Produces |
|---|---|
| `S01_age_period_cohort.R` | Operational age separated from calendar time and deployment cohort → **Section S7; Table S7.1; Figure S7.1** |
| `S02_ab_asymmetry.R` | Channel-asymmetry tests: temperature-dependent signed bias, channel-specific noise → **Section S8; Figures S8.1–S8.2** |
| `S03_prospective_flagging.R` | Prospectively applicable degradation-flag variants → **Section S6; Table S6.1** |
| `S04_wildfire_season.R` | CONUS exposure models stratified by fire season × west/east → **Section S9; Table S9.1** |
| `S05_trajectory_support.R` | Sensitivity of lifetime trajectories to truncating sparse oldest-age tails |
| `S06_exposure_windows.R` | Reported-day vs calendar-day windows × all-flagged outcome exclusion → **Table S9.2** |
| `S07_calibration_cv.R` | Sensor-grouped 5-fold cross-validation of the calibration models → **Table S3.2** |

## `data_prep/` — upstream data assembly

Only needed to rebuild the dataset from raw sources: PurpleAir download/merge, Google Earth Engine export of ERA5-Land and GEOS-CF, the reference join, and the original cleaning chunk. Requires data access and cloud credentials (see REPRODUCING.md).

## `extra/` — retained but not used in the paper

The indoor-sensor download branch (the paper analyzes outdoor sensors only) and two superseded passes: an earlier join pipeline (disabled in the source notebook) and an earlier persistent-degradation run using the 85th-percentile flag. Retained for completeness of the record. **Do not run these.**

## Verification

Run `Rscript verify_paper_numbers.R` to check the paper's reported values against the analysis outputs. It reads existing output files only (no model fitting) and completes in seconds: **34/34 checks pass**, covering the global, regional and climate-zone age models and all seven robustness analyses.

## Provenance

These scripts were extracted verbatim from the working notebook the analyses were run from; **all 37 analysis chunks were verified byte-identical**, and no analysis code was edited. The notebook itself is not distributed (it contained exploratory and superseded code). See **[PROVENANCE.md](PROVENANCE.md)** for the verification, the notebook's checksum, the full chunk-to-script map, and the three deliberate post-extraction changes.

## Reproducing

See **[REPRODUCING.md](REPRODUCING.md)** for data access, GEE/GCS setup, path configuration, run order, the row-count checkpoint, and a manuscript-to-code map for every table and SI section.

## Data

PurpleAir data were obtained through the U.S. EPA PurpleAir API and the Remote Sensing Information Gateway (RSIG) under MTA #1650-26. GEOS-CF and ERA5-Land are publicly available. Raw sensor data are not redistributed here; scripts use absolute paths rooted at `/Users/priyanka/Downloads/PurpleAir_EPA/` — see `00_config.R` for a one-line rewrite. No credentials are stored in this repository.

## Requirements

R (≥ 4.3) with: `data.table`, `lubridate`, `ggplot2`, `scales`, `dplyr`, `patchwork`, `sf`, `mgcv`, `fixest`, `survival`, `sandwich`, `lmtest`, `splines`, `zoo`, `terra`, `stringr`, `cowplot`, `rnaturalearth`, `rnaturalearthdata`, `flextable`, `officer`, `arrow`, `duckdb`, and `rgee` / `googleCloudStorageR` (data-prep stages only).
