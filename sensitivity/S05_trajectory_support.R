# ============================================================
# 05_truncate_trajectory_tails.R
#
# TRUNCATE TRAJECTORY DISPLAYS WHERE SENSOR SUPPORT IS THIN
#
# The right ends of the Figure 2/3 trajectories (6+ years) rest on very
# few sensors. This drop-in helper truncates each plotted series at the
# first operational-age bin where fewer than MIN_SENSORS_DISPLAY sensors
# remain, so the figures stop where the data get thin.
#
# IMPORTANT: this is a DISPLAY rule only. Leave the trend models
# (Models 1-3) untouched — they are observation-weighted and
# sensor-clustered, so sparse tail bins carry almost no weight; the
# figures are where thin tails mislead.
#
# HOW TO USE in the cf1 degradation chunk (and the model-divergence
# equivalent):
#   1) source() this file, or paste the function in.
#   2) wherever a binned summary table feeds ggplot (the *_bin tables
#      built per group before plotting), apply:
#         plot_dt <- truncate_low_support(plot_dt,
#                        group_vars = c("region"),      # or climate_zone, etc.
#                        min_sensors = 100L)
#   3) add to the caption: "Trajectories are truncated at the first
#      operational-age bin with fewer than 100 contributing sensors."
#
# The function expects station-level rows (one row per station x bin,
# like the make_sens_bin() output). If your plotting table is already
# aggregated across stations but carries an n_sensors column, use
# truncate_low_support_agg() instead.
# ============================================================

suppressPackageStartupMessages(library(data.table))

MIN_SENSORS_DISPLAY <- 100L

# --- for station-level tables (station x op_bin rows) -------------
truncate_low_support <- function(dt, group_vars = character(),
                                 min_sensors = MIN_SENSORS_DISPLAY,
                                 bin_var = "op_bin",
                                 station_var = "station") {
  dt <- as.data.table(dt)
  support <- dt[, .(n_sensors_bin = uniqueN(get(station_var))),
                by = c(group_vars, bin_var)]
  setorderv(support, c(group_vars, bin_var))
  # first under-supported bin per group; keep everything strictly before it
  support[, under := n_sensors_bin < min_sensors]
  cutoffs <- support[, .(
    cutoff_bin = if (any(under)) min(get(bin_var)[under]) else Inf
  ), by = group_vars]
  if (length(group_vars)) {
    dt <- cutoffs[dt, on = group_vars]
  } else {
    dt[, cutoff_bin := cutoffs$cutoff_bin]
  }
  out <- dt[get(bin_var) < cutoff_bin]
  dropped <- nrow(dt) - nrow(out)
  message(sprintf(
    "truncate_low_support: dropped %s rows (%.1f%%) beyond the <%d-sensor cutoff",
    format(dropped, big.mark = ","), 100 * dropped / max(nrow(dt), 1),
    min_sensors))
  out[, cutoff_bin := NULL]
  out[]
}

# --- for pre-aggregated tables that carry n_sensors per bin -------
truncate_low_support_agg <- function(dt, group_vars = character(),
                                     min_sensors = MIN_SENSORS_DISPLAY,
                                     bin_var = "op_bin",
                                     nsensor_var = "n_sensors") {
  dt <- as.data.table(dt)
  setorderv(dt, c(group_vars, bin_var))
  dt[, under := get(nsensor_var) < min_sensors]
  dt[, cutoff_bin := if (any(under)) min(get(bin_var)[under]) else Inf,
     by = group_vars]
  out <- dt[get(bin_var) < cutoff_bin]
  out[, c("under", "cutoff_bin") := NULL]
  out[]
}

# --- report the truncation ages (for the SI / captions) -----------
report_cutoffs <- function(dt, group_vars = character(),
                           min_sensors = MIN_SENSORS_DISPLAY,
                           bin_var = "op_bin",
                           station_var = "station") {
  dt <- as.data.table(dt)
  support <- dt[, .(n_sensors_bin = uniqueN(get(station_var))),
                by = c(group_vars, bin_var)]
  setorderv(support, c(group_vars, bin_var))
  support[, under := n_sensors_bin < min_sensors]
  support[, .(
    cutoff_years = if (any(under)) min(get(bin_var)[under]) / (365.25 * 24)
                   else max(get(bin_var)) / (365.25 * 24),
    reached_cutoff = any(under)
  ), by = group_vars]
}

# ------------------------------------------------------------------
# Example (cf1 pipeline):
#
#   sens_bin_global <- make_sens_bin(dt_global)          # existing code
#   plot_global <- truncate_low_support(sens_bin_global,
#                      group_vars = character(),          # global panel
#                      min_sensors = 100L)
#   plot_region <- truncate_low_support(sens_bin_global,
#                      group_vars = "region",
#                      min_sensors = 100L)
#   report_cutoffs(sens_bin_global, "region")             # for captions
#
# Manuscript caption addition (Figures 2-3):
#   "Trajectories are truncated at the first operational-age bin with
#    fewer than 100 contributing sensors."
# ------------------------------------------------------------------
