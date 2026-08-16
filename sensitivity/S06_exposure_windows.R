# ============================================================
# 06_exposure_window_sensitivity.R
#
# SENSITIVITY OF THE HIGH-PM EXPOSURE MODELS TO:
#   (A) ROW-BASED vs CALENDAR-BASED WINDOWS  (review flag 2a)
#   (B) THE ALL-FLAGGED OUTCOME EXCLUSION    (review flag 4)
#
# The paper's exposure panel rolls windows over *reported* operational
# days, so "7-day" windows stretch across reporting gaps; and the panel
# filter `future_total_h > future_flagged_h` drops station-days whose
# entire future window is flagged. This script re-estimates the
# headline GLOBAL joint burden+spike models (threshold 35 µg/m³;
# 7/14/28-day windows; outcomes: future |A-B|, future flagged rate,
# future |PA-GEOS|) under a 2x2 design:
#
#   windows:   {reported-day rows (paper) , true calendar days}
#   exclusion: {drop all-flagged futures (paper) , keep them}
#
# If coefficients are stable across the four cells, one definitional
# sentence in the methods covers both concerns; the comparison table is
# SI-ready either way.
#
# INPUT:  outdoor_final.RData  (run 00_build_analysis_dataset.R first!)
# OUTPUT: exposure_window_sensitivity/ CSV + Word table
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
  library(splines)
  library(flextable)
  library(officer)
})

global_rdata <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_final.RData"
outdir <- "/Users/priyanka/Downloads/PurpleAir_EPA/exposure_window_sensitivity"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

PM_MAX <- 1500
THR <- 35
WINDOWS <- c(7L, 14L, 28L)
FUTURE_WINDOW <- 7L
DROP_FIRST_DAYS <- 7L
MIN_FUTURE_HOURS <- 72
FLAG_ABS_CUTOFF <- 5
FLAG_PCTL <- 0.85
CAL_COVERAGE <- 0.7          # calendar windows need >=70% of days observed

a_col_default <- "pm25_cf_1_a_hourly(ug/m3)_pm25_cf_1_a"
b_col_default <- "pm25_cf_1_b_hourly(ug/m3)_pm25_cf_1_b"

say <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")),
      paste0(..., collapse = ""), "\n")
}

# ------------------------------------------------------------
# 0) PANEL CACHE: skip the expensive build on reruns
# ------------------------------------------------------------
cache_file <- file.path(outdir, "sensitivity_panel_cache.rds")
if (file.exists(cache_file)) {
  say("Loading cached panel: ", cache_file)
  panel <- readRDS(cache_file)
  SKIP_BUILD <- TRUE
} else SKIP_BUILD <- FALSE

if (!SKIP_BUILD) {
# ------------------------------------------------------------
# 1) LOAD + HOURLY METRICS (paper-identical)
# ------------------------------------------------------------
say("Loading: ", global_rdata)
e <- new.env(); load(global_rdata, envir = e)
outdoor <- if (exists("outdoor_final", envir = e)) get("outdoor_final", envir = e) else get("outdoor", envir = e)
setDT(outdoor); rm(e)

a_col <- if (a_col_default %in% names(outdoor)) a_col_default else "cf1_a"
b_col <- if (b_col_default %in% names(outdoor)) b_col_default else "cf1_b"
dt <- outdoor[, c("station", "region", "datetime_utc", a_col, b_col), with = FALSE]
rm(outdoor); gc()
setnames(dt, c(a_col, b_col), c("cf1_a", "cf1_b"))
dt[, `:=`(station = as.character(station),
          ts = as.POSIXct(datetime_utc, tz = "UTC"),
          cf1_a = as.numeric(cf1_a), cf1_b = as.numeric(cf1_b))]
dt[, datetime_utc := NULL]
dt <- dt[is.finite(cf1_a) & is.finite(cf1_b) &
           cf1_a >= 0 & cf1_b >= 0 & cf1_a <= PM_MAX & cf1_b <= PM_MAX]

setorder(dt, station, ts)
dt[, first_ts := ts[1L], by = station]
dt[, op_hours := as.numeric(difftime(ts, first_ts, units = "hours"))]
dt <- dt[op_hours >= DROP_FIRST_DAYS * 24]
dt[, pm := (cf1_a + cf1_b) / 2]
dt[, absdiff := abs(cf1_a - cf1_b)]
dt[, pctdiff := fifelse((cf1_a + cf1_b) > 0, 2 * absdiff / (cf1_a + cf1_b), NA_real_)]
thr_dt <- dt[is.finite(pctdiff),
             .(thr85 = quantile(pctdiff, FLAG_PCTL, na.rm = TRUE, type = 7)), by = station]
dt <- thr_dt[dt, on = "station"]
dt[, flagged := as.integer(absdiff > FLAG_ABS_CUTOFF & pctdiff > thr85)]
rm(thr_dt); gc()

# ------------------------------------------------------------
# 2) DAILY PANEL (calendar-dated), THEN FULL CALENDAR GRID
# ------------------------------------------------------------
say("Building daily panel")
dt[, date := as.Date(ts)]
daily <- dt[, .(
  n_hours = .N,
  n_flag = sum(flagged, na.rm = TRUE),
  pm_day = mean(pm, na.rm = TRUE),
  pm_sd_day = sd(pm, na.rm = TRUE),
  absdiff_day = mean(absdiff, na.rm = TRUE),
  age_days = min(op_hours) / 24
), by = .(station, date)]
rm(dt); gc()

say("Expanding to calendar grid")
grid <- daily[, .(date = as.Date(
  seq(as.numeric(min(date)), as.numeric(max(date)), by = 1),
  origin = "1970-01-01")), by = station]
daily <- daily[grid, on = c("station", "date")]
setorder(daily, station, date)
daily[, obs := as.integer(is.finite(pm_day))]
daily[, excess := pmax(pm_day - THR, 0)]

# ------------------------------------------------------------
# 3) ROLLING EXPOSURES/OUTCOMES UNDER BOTH WINDOW DEFINITIONS
#    ROW  = roll over observed rows only (paper's construction)
#    CAL  = roll over the calendar grid with a coverage rule
# ------------------------------------------------------------
say("Computing rolling windows (both definitions)")

# CAL: on the full grid
for (pw in WINDOWS) {
  daily[, paste0("burden_cal_", pw) :=
          shift(frollsum(excess, pw, align = "right", na.rm = TRUE), 1L), by = station]
  daily[, paste0("spike_cal_", pw) := {
    s <- shift(frollapply(pm_day, pw,
          function(z) if (all(is.na(z))) NA_real_ else max(z, na.rm = TRUE),
          align = "right"), 1L)
    pmax(s - THR, 0) }, by = station]
  daily[, paste0("pastobs_cal_", pw) :=
          shift(frollsum(obs, pw, align = "right"), 1L), by = station]
}
daily[, fut_abs := shift(frollmean(absdiff_day, FUTURE_WINDOW, align = "left", na.rm = TRUE), -1L), by = station]
daily[, fut_flag_h := shift(frollsum(n_flag, FUTURE_WINDOW, align = "left", na.rm = TRUE), -1L), by = station]
daily[, fut_tot_h := shift(frollsum(n_hours, FUTURE_WINDOW, align = "left", na.rm = TRUE), -1L), by = station]
daily[, prev_flag := shift(n_flag, 1L), by = station]
daily[, prev_n := shift(n_hours, 1L), by = station]

# ROW: collapse to observed rows and roll there (mirrors the op-bin rows)
obsd <- daily[obs == 1L]
setorder(obsd, station, date)
for (pw in WINDOWS) {
  obsd[, paste0("burden_row_", pw) :=
         shift(frollsum(excess, pw, align = "right", na.rm = TRUE), 1L), by = station]
  obsd[, paste0("spike_row_", pw) := {
    s <- shift(frollapply(pm_day, pw, max, na.rm = TRUE, align = "right"), 1L)
    pmax(s - THR, 0) }, by = station]
}
obsd[, fut_abs_row := shift(frollmean(absdiff_day, FUTURE_WINDOW, align = "left"), -1L), by = station]
obsd[, fut_flag_h_row := shift(frollsum(n_flag, FUTURE_WINDOW, align = "left"), -1L), by = station]
obsd[, fut_tot_h_row := shift(frollsum(n_hours, FUTURE_WINDOW, align = "left"), -1L), by = station]

panel <- obsd  # every model row is an observed station-day
rm(daily, grid); gc()
saveRDS(panel, file.path(outdir, "sensitivity_panel_cache.rds"))
say("Panel cached to sensitivity_panel_cache.rds")

}  # end !SKIP_BUILD

panel[, p_prev := (prev_flag + 0.5) / (prev_n + 1.0)]
panel[, off_prev := qlogis(p_prev)]

# ------------------------------------------------------------
# 4) MODELS: 2 window defs x 2 exclusion rules x 3 windows x 3 outcomes
# ------------------------------------------------------------
fit_cell <- function(dd, burden_col, spike_col, fut_abs_col, fut_flag_col,
                     fut_tot_col, windows_lab, excl_lab, pw) {
  dd <- dd[is.finite(get(burden_col)) & is.finite(get(spike_col)) &
             is.finite(get(fut_abs_col)) & is.finite(get(fut_tot_col)) &
             get(fut_tot_col) >= MIN_FUTURE_HOURS &
             is.finite(pm_day) & is.finite(pm_sd_day) & is.finite(off_prev)]
  if (excl_lab == "drop all-flagged (paper)")
    dd <- dd[get(fut_tot_col) > get(fut_flag_col)]
  if (nrow(dd) < 10000L) return(data.table())
  dd[, b100 := get(burden_col) / 100]
  dd[, s10 := get(spike_col) / 10]

  out <- list()
  m1 <- tryCatch(feols(
    as.formula(paste0(fut_abs_col,
      " ~ b100 + s10 + ns(pm_day,3) + ns(pm_sd_day,3) + ns(age_days,3) | station")),
    data = dd, cluster = ~station), error = function(e) NULL)
  if (!is.null(m1)) {
    ct <- coeftable(m1)
    out[["abs"]] <- data.table(outcome = "Future |A-B|",
      est = ct["b100","Estimate"], se = ct["b100","Std. Error"], n = m1$nobs)
  }
  # fixest::feglm does not accept cbind() outcomes: use proportion + weights
  dd[, flag_prop := get(fut_flag_col) / pmax(get(fut_tot_col), 1)]
  dd[, flag_wt := pmax(get(fut_tot_col), 1)]
  m2 <- tryCatch(feglm(
    flag_prop ~ b100 + s10 + ns(pm_day,3) + ns(pm_sd_day,3) +
      ns(age_days,3) | station,
    data = dd, family = binomial(), weights = ~flag_wt, offset = ~off_prev,
    cluster = ~station), error = function(e) NULL)
  if (!is.null(m2)) {
    ct <- coeftable(m2)
    out[["flag"]] <- data.table(outcome = "Future flagged (log-odds)",
      est = ct["b100","Estimate"], se = ct["b100","Std. Error"], n = m2$nobs)
  }
  res <- rbindlist(out)
  if (nrow(res)) res[, `:=`(windows = windows_lab, exclusion = excl_lab,
                            past_window_days = pw,
                            n_sensors = uniqueN(dd$station))]
  res
}

results <- list()
for (pw in WINDOWS) {
  say("Window ", pw, " days")
  for (excl in c("drop all-flagged (paper)", "keep all rows")) {
    results[[length(results)+1]] <- fit_cell(panel,
      paste0("burden_row_", pw), paste0("spike_row_", pw),
      "fut_abs_row", "fut_flag_h_row", "fut_tot_h_row",
      "Reported-day rows (paper)", excl, pw)
    dd_cal <- panel[get(paste0("pastobs_cal_", pw)) >= ceiling(CAL_COVERAGE * pw)]
    results[[length(results)+1]] <- fit_cell(dd_cal,
      paste0("burden_cal_", pw), paste0("spike_cal_", pw),
      "fut_abs", "fut_flag_h", "fut_tot_h",
      "Calendar days (>=70% coverage)", excl, pw)
  }
}
res <- rbindlist(results)
res[, `:=`(ci_lwr = est - 1.96 * se, ci_upr = est + 1.96 * se)]
res[outcome == "Future flagged (log-odds)",
    `:=`(or = exp(est), or_lwr = exp(ci_lwr), or_upr = exp(ci_upr))]
fwrite(res, file.path(outdir, "exposure_window_sensitivity.csv"))
print(res[, .(past_window_days, windows, exclusion, outcome,
              est = round(est, 3), ci_lwr = round(ci_lwr, 3),
              ci_upr = round(ci_upr, 3), n_sensors)])

# ------------------------------------------------------------
# 5) WORD TABLE
# ------------------------------------------------------------
ft_dt <- res[, .(
  `Window (days)` = past_window_days,
  `Window definition` = windows, Exclusion = exclusion, Outcome = outcome,
  `Burden effect per +100 µg/m³-days (95% CI)` =
    sprintf("%.3f (%.3f, %.3f)", est, ci_lwr, ci_upr),
  `N sensors` = format(n_sensors, big.mark = ","))]
ft <- autofit(fontsize(flextable(ft_dt), size = 8, part = "all"))
doc <- read_docx()
doc <- body_add_par(doc, paste0(
  "Table SX: Sensitivity of global joint burden + spike models (threshold 35 µg/m³) ",
  "to the rolling-window definition (reported operational days, as in the main ",
  "analysis, versus true calendar days with >=70% coverage) and to the exclusion of ",
  "station-days whose entire future window is flagged. Sensor fixed effects; ",
  "sensor-clustered SEs."), style = "Normal")
doc <- body_add_flextable(doc, ft)
print(doc, target = file.path(outdir, "exposure_window_sensitivity.docx"))
say("Done. Outputs in: ", outdir)
