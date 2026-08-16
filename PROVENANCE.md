# Provenance

The analysis scripts in this repository were not rewritten from scratch. They were
extracted **verbatim** from the working R Markdown notebook that the published
analyses were actually run from (`purpleair_degradation2.Rmd`, 24,215 lines,
38 code chunks). No analysis code was edited during extraction.

## Verification performed

Every code chunk in the notebook was compared byte-for-byte against the
corresponding section of the extracted scripts:

- **37 of 37** analysis chunks captured (chunk #1 is the `knitr` setup chunk and
  carries no analysis code).
- **0 byte-level mismatches.**
- Verification was re-run after every reorganisation of the files.

Each script retains a header comment naming the source chunk number and its line
range in the notebook, and each extracted block is preceded by a
`# ---- source chunk #N (Rmd lines A-B) ----` marker.

## Source notebook checksum

The notebook itself is not distributed here (it contained exploratory and
superseded code, an indoor-sensor branch not used in the paper, and cloud
identifiers). Its identity is recorded so the extraction can be re-verified by
anyone holding the archival copy:

```
file:    purpleair_degradation2.Rmd
lines:   24215
chunks:  38 (37 analysis + 1 setup)
sha256:  196d671eacbff5b339f2afd3b95db045b6cb82535dee329ac1edd822ce3cff45
```

That checksum is for the archival copy with the Google Earth Engine project ID
and Google Cloud Storage bucket replaced by the placeholders `YOUR-GEE-PROJECT-ID`
and `YOUR-GCS-BUCKET`; those substitutions are the only differences from the
notebook as run.

## Deliberate departures from the notebook

Three changes were made *after* extraction, each recorded in the git history:

1. **Removed an unsafe dataset-rebuild fallback** from `sensitivity/S01`–`S04`.
   The fallback reconstructed `outdoor_final.RData` if it was absent, but applied
   only the temperature and relative-humidity filters — it omitted the
   cf_1 > 1500 µg/m³ exclusion, so it could silently produce a dataset differing
   from the published one. Those scripts now stop and direct the user to
   `00_build_analysis_dataset.R`. This changes no published result.
2. **Made the early-life flag variant strictly prospective** in
   `sensitivity/S03_prospective_flagging.R`. The sensor-specific threshold is
   estimated on each sensor's first 90 operational days; those hours are now left
   unflagged rather than being flagged against a threshold derived from them.
   This analysis was re-run and SI section S6 / Table S6.1 updated accordingly.
   No other result is affected.
3. **Added `00_config.R`**, a new file documenting paths, cloud identifiers and
   key analysis constants. It contains no analysis logic.

## Chunk-to-script map

| Rmd chunk | Rmd lines | Script |
|---|---|---|
| #2 | 12-257 | `data_prep/01_download_purpleair_outdoor.R` |
| #3 | 260-505 | `extra/indoor_download_not_in_paper.R` |
| #4 | 508-752 | `data_prep/01_download_purpleair_outdoor.R` |
| #5 | 755-997 | `extra/indoor_download_not_in_paper.R` |
| #6 | 1000-1244 | `data_prep/01_download_purpleair_outdoor.R` |
| #7 | 1247-1489 | `extra/indoor_download_not_in_paper.R` |
| #8 | 1492-1736 | `data_prep/01_download_purpleair_outdoor.R` |
| #9 | 1739-1982 | `extra/indoor_download_not_in_paper.R` |
| #10 | 1985-2001 | `data_prep/01_download_purpleair_outdoor.R` |
| #11 | 2004-2172 | `data_prep/02_gee_extract_era5.R` |
| #12 | 2175-2360 | `data_prep/02_gee_extract_era5.R` |
| #13 | 2363-2548 | `data_prep/02_gee_extract_era5.R` |
| #14 | 2551-2737 | `data_prep/02_gee_extract_era5.R` |
| #15 | 2740-2931 | `data_prep/03_gee_extract_geoscf.R` |
| #16 | 2934-3012 | `data_prep/03_gee_extract_geoscf.R` |
| #17 | 3016-3810 | `extra/superseded_join_pipeline.R` |
| #18 | 3813-4649 | `data_prep/04_join_reference_parquet.R` |
| #19 | 4652-4922 | `data_prep/05_combine_joined.R` |
| #20 | 4925-5319 | `data_prep/06_clean_outdoor_final.R` |
| #21 | 5322-6334 | `01_descriptive_statistics.R` |
| #22 | 6337-6673 | `02_climate_zone_assignment.R` |
| #23 | 6676-7000 | `02_climate_zone_assignment.R` |
| #24 | 7003-7312 | `02_climate_zone_assignment.R` |
| #25 | 7315-7623 | `01_descriptive_statistics.R` |
| #26 | 7626-8539 | `01_descriptive_statistics.R` |
| #27 | 8542-9219 | `01_descriptive_statistics.R` |
| #28 | 9222-9987 | `01_descriptive_statistics.R` |
| #29 | 9991-11579 | `03_primary_age_degradation_models.R` |
| #30 | 11582-13472 | `04_model_reference_divergence.R` |
| #31 | 13476-15228 | `05_survivorship_attrition.R` |
| #32 | 15231-17396 | `06_high_pm_exposure.R` |
| #33 | 17399-18951 | `06_high_pm_exposure.R` |
| #34 | 18955-19630 | `09_pa_epa_validation.R` |
| #35 | 19633-20968 | `09_pa_epa_validation.R` |
| #36 | 20971-21808 | `07_spatial_analysis.R` |
| #37 | 21811-23034 | `extra/superseded_persistent_degradation_q85.R` |
| #38 | 23037-24214 | `08_persistent_degradation.R` |
