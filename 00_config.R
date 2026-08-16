# ==========================================================================
# 00_config.R
#
# Central configuration for the PurpleAir degradation analysis.
#
# The analysis scripts in this directory are VERBATIM extracts from the
# original working notebook (purpleair_degradation2.Rmd) and therefore still
# contain hard-coded absolute paths and cloud identifiers. Nothing was edited,
# so the scripts remain exactly the code that produced the published numbers.
#
# To run them on your own machine, set the values below and then rewrite the
# paths once (see "One-time path rewrite" at the bottom of this file).
# ==========================================================================

# --- Where your data lives -------------------------------------------------
# Root directory containing the PurpleAir downloads, joined files, and
# outdoor_final.RData. All scripts currently expect:
#   /Users/priyanka/Downloads/PurpleAir_EPA
DATA_DIR <- "/path/to/your/PurpleAir_EPA"

# --- Where outputs should be written --------------------------------------
# The scripts write figures, tables and CSVs into subdirectories of DATA_DIR.
OUTPUT_DIR <- file.path(DATA_DIR, "outputs")

# --- Google Earth Engine / Google Cloud Storage ----------------------------
# Needed only for the extraction scripts (02_gee_extract_era5.R,
# data_prep/03_gee_extract_geoscf.R). In those files replace the placeholders
# YOUR-GEE-PROJECT-ID and YOUR-GCS-BUCKET with these values.
GEE_PROJECT <- "your-gee-project-id"
GCS_BUCKET  <- "your-gcs-bucket"

# --- Key analysis constants (documented here for reference; do not change
#     unless you intend to depart from the published specification) ---------
FLAG_ABS_THRESHOLD <- 5     # ug/m3; |A-B| criterion for flagged degradation
FLAG_REL_Q         <- 0.85  # sensor-specific percentile of relative divergence
                            # NOTE: 08_persistent_degradation.R
                            # deliberately uses 0.90 (see SI section S5.1)
PM_MAX             <- 1500  # ug/m3; cf_1 plausibility ceiling
MAX_MODEL_ROWS     <- 2.5e6 # sampling cap for the largest model fits
LANDMARK_HOURS     <- 10000 # landmark for the Cox survivorship models
ROLLING_WINDOW_DAYS <- 180  # persistent-degradation rolling window

# ==========================================================================
# One-time path rewrite
#
# There are 179 hard-coded occurrences of the original data root across
# R/, R/extra/ and analysis_addons/. Rewrite them all in place with:
#
#   cd <repo root>
#   grep -rl '/Users/priyanka/Downloads/PurpleAir_EPA' .  \
#     | xargs sed -i.bak 's|/Users/priyanka/Downloads/PurpleAir_EPA|/your/data/dir|g'
#
# (On macOS use `sed -i ''` instead of `sed -i.bak` if you do not want
# backup files.) Verify afterwards with:
#
#   grep -rc '/Users/priyanka' . 
# ==========================================================================
