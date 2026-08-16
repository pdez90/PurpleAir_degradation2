# Reproducing the analysis

This repository contains the complete analysis code for the paper. It is a research pipeline, not a one-click build: reproducing the results end-to-end requires (1) data access, (2) your own paths and cloud credentials, and (3) running the scripts in order.

If you only want to **check that the reported numbers match the code's outputs**, skip to [§5 Verifying](#5-verifying) — that takes seconds and needs no rerun.

## 1. What you need

**Data.** Raw data are not redistributed here. Hourly outdoor PurpleAir observations (cf_1 channels A and B, temperature, RH; 2018–2025) came from the U.S. EPA PurpleAir API and the Remote Sensing Information Gateway (RSIG) under MTA #1650-26 — request your own access from EPA/RSIG or PurpleAir. Regulatory PM2.5 for the calibration analysis comes from EPA AQS monitors within 100 m of a PurpleAir sensor.

**Cloud accounts** (only for `data_prep/`). A Google Earth Engine project and a Google Cloud Storage bucket. No credentials are stored in this repository; replace `YOUR-GEE-PROJECT-ID` and `YOUR-GCS-BUCKET` in `data_prep/02_gee_extract_era5.R` and `data_prep/03_gee_extract_geoscf.R`. Authentication is interactive via `rgee` / `googleCloudStorageR` — never commit a service-account JSON.

**Software.** R ≥ 4.3 with `data.table`, `lubridate`, `ggplot2`, `scales`, `dplyr`, `patchwork`, `sf`, `mgcv`, `fixest`, `survival`, `sandwich`, `lmtest`, `splines`, `zoo`, `terra`, `stringr`, `cowplot`, `rnaturalearth(data)`, `flextable`, `officer`, `arrow`, `duckdb`; plus `rgee` (needs a working Python earthengine-api) for `data_prep/`.

**Hardware.** The analytic dataset is 452 million hourly rows. The heaviest scripts hold it in memory alongside a second dataset and have peaked near 100 GB; 64 GB is a practical minimum and more is better. Model fits are capped at 2.5 M sampled rows (`MAX_MODEL_ROWS`) where noted.

**Paths.** Scripts use absolute paths rooted at `/Users/priyanka/Downloads/PurpleAir_EPA/`. Edit `00_config.R` — it documents every setting and gives a one-line `sed` command that rewrites all of them at once.

## 2. Run order

The numbering is a real dependency chain, not a suggestion:

| Step | Scripts | Notes |
|---|---|---|
| Raw assembly (optional) | `data_prep/01`–`06` | Only if rebuilding from source data; needs GEE/GCS |
| **Build dataset** | `00_build_analysis_dataset.R` | Writes `outdoor_final.RData`; **every later script loads it** |
| Climate join | `02_climate_zone_assignment.R` | Writes `outdoor_with_koppen_timevarying.RData`, read by `03`–`08` |
| Analyses | `01`, `03`–`09` | Otherwise independent; each loads its own inputs |
| Robustness | `sensitivity/S01`–`S07` | Require `outdoor_final.RData` to exist |

**Checkpoint.** `outdoor_final.RData` must contain **452,143,148 rows** from **36,328 sensors** (CONUS 439,768,685; Africa 3,134,219; India 2,798,734; South America 6,441,510). `00_build_analysis_dataset.R` verifies these counts and stops if they do not match.

The `sensitivity/` scripts deliberately do **not** rebuild the dataset — they stop and point at `00_build_analysis_dataset.R` — so there is no way to end up analysing a file that differs from the published one.

`extra/` holds code retained from the working notebook but not used in the paper (the indoor-sensor branch and two superseded passes). Do not run it.

## 3. Manuscript ↔ code map

| Paper element | Script |
|---|---|
| Analytic sample (452,143,148 rows / 36,328 sensors) | `00_build_analysis_dataset.R` |
| Tables S1–S2, Figures S1, S4–S25 | `01_descriptive_statistics.R` |
| Figures S2–S3 (climate zones) | `02_climate_zone_assignment.R` |
| **Figure 2; Tables S3–S6** | `03_primary_age_degradation_models.R` |
| **Figure 3; Tables S7–S9; Figures S36–S38** | `04_model_reference_divergence.R` |
| **Figure 4; Tables S10–S18** (S18 = effective lifetime) | `05_survivorship_attrition.R` |
| **Figure 5; Tables S19–S24; Figures S44–S49** | `06_high_pm_exposure.R` |
| Section S4 (spatial dependence) | `07_spatial_analysis.R` |
| Section S5 (persistent/permanent; 1,146 / 206) | `08_persistent_degradation.R` |
| Section S3, Table S3.1 | `09_pa_epa_validation.R` |
| Table S3.2 (calibration CV) | `sensitivity/S07_calibration_cv.R` |
| Section S6, Table S6.1 (prospective flagging) | `sensitivity/S03_prospective_flagging.R` |
| Section S7, Table S7.1 (age vs calendar/cohort) | `sensitivity/S01_age_period_cohort.R` |
| Section S8 (channel asymmetry) | `sensitivity/S02_ab_asymmetry.R` |
| Section S9, Table S9.1 (wildfire season) | `sensitivity/S04_wildfire_season.R` |
| Table S9.2 (exposure-window sensitivity) | `sensitivity/S06_exposure_windows.R` |
| Trajectory-tail sensitivity | `sensitivity/S05_trajectory_support.R` |

Two parameter choices differ from the primary analysis by design and are documented in the SI: `08_persistent_degradation.R` uses the sensor-specific **90th**-percentile relative-divergence flag (SI S5.1), and the early-life variant in `sensitivity/S03` leaves the 90-day estimation window unflagged so that flagging is strictly prospective (SI S6.1).

## 4. What has been verified

- **Code provenance.** Every analysis chunk of the original working notebook was extracted byte-for-byte into these scripts: 37 of 37 chunks, zero mismatches. See [PROVENANCE.md](PROVENANCE.md).
- **Syntax.** All scripts parse cleanly under R.
- **Dataset.** An independent rerun of `03_primary_age_degradation_models.R` reproduced 452,143,148 rows, 36,328 sensors and 22,581,395 global station-bins — the last matching the stored output exactly.
- **Reported numbers.** 34 of 34 checks pass (see below), covering the global, regional and climate-zone age models and all seven robustness analyses.

## 5. Verifying

`verify_paper_numbers.R` compares the values printed in the paper against the analysis output files. It reads existing outputs only — it fits no models, so it runs in seconds:

```bash
Rscript verify_paper_numbers.R
```

Each check prints PASS, FAIL, or SKIP (SKIP = that script's outputs are not present). Edit `ROOT` at the top if your data directory differs.

## 6. Known limits

- Raw PurpleAir and AQS extracts are not included (data-use terms).
- Package versions are not pinned; the analyses ran under R 4.5.x.
- Sampled fits use fixed seeds (`set.seed(123)`) and are reproducible on the same data.
- GEE extractions depend on the upstream ERA5-Land and GEOS-CF collections, which are versioned by Google/NASA/ECMWF and may be revised.
