# ============================================================
# 04_wildfire_season_stratification.R
#
# IS THE CUMULATIVE-BURDEN EFFECT A WILDFIRE-SMOKE EFFECT?
#
# The CONUS network is concentrated in Csb/Csa (California / Pacific
# Northwest), so high-PM exposure there is plausibly dominated by smoke
# episodes. This script stratifies the exposure -> future degradation
# analysis by (a) fire season (Jun-Oct) vs non-fire season, and (b)
# western vs eastern CONUS (lon < -100°), and tests the interaction.
#
# Design mirrors the paper's daily fixed-effects panel, simplified to
# the headline specification: 7-day cumulative burden above 35 µg/m³ +
# spike exceedance (joint model), sensor FE, spline in age and current
# PM, sensor-clustered SEs. Outcomes: future 7-day |A-B| and future
# flagged rate.
#
# Interpretation guide:
#   - burden effect concentrated in fire-season x west  -> "the
#     cumulative-burden association in CONUS is driven largely by
#     wildfire-smoke episodes"
#   - similar across strata -> "the association is not specific to
#     wildfire smoke"
#
# NOTE: the rigorous upgrade is NOAA HMS smoke-plume polygons
# (https://www.ospo.noaa.gov/Products/land/hms.html) joined to
# station-days by location; the hook below marks where that flag would
# replace fire_season. Season x geography is the honest first cut and
# is defensible for an SI analysis.
#
# INPUT:  outdoor_final.RData  (CONUS subset)
# OUTPUT: wildfire_stratification/ CSVs + Word table
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
  library(splines)
  library(flextable)
  library(officer)
})

global_rdata <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_final.RData"
outdir <- "/Users/priyanka/Downloads/PurpleAir_EPA/wildfire_stratification"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

PM_MAX <- 1500
BURDEN_THRESHOLD <- 35
PAST_WINDOW <- 7L
FUTURE_WINDOW <- 7L
MIN_FUTURE_DAYS <- 3L
FIRE_MONTHS <- 6:10          # Jun-Oct
WEST_LON_CUT <- -100         # west of 100°W
FLAG_ABS_CUTOFF <- 5
FLAG_PCTL <- 0.85

a_col_default <- "pm25_cf_1_a_hourly(ug/m3)_pm25_cf_1_a"
b_col_default <- "pm25_cf_1_b_hourly(ug/m3)_pm25_cf_1_b"
lon_candidates <- c("LONGITUDE(deg)", "longitude", "lon", "LONGITUDE")

say <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")),
      paste0(..., collapse = ""), "\n")
}


load_outdoor_object <- function(infile) {
  # outdoor_final.RData must already exist and must be the paper's exact
  # analytic dataset (452,143,148 rows / 36,328 sensors). It is NOT rebuilt
  # here: an ad-hoc rebuild risks producing a file that differs from the
  # published dataset (e.g. omitting the cf_1 > 1500 ug/m3 exclusion).
  if (!file.exists(infile)) {
    stop("outdoor_final.RData not found at: ", infile,
         "\n  Run 00_build_analysis_dataset.R first; it reproduces and",
         "\n  verifies the paper's exact analytic dataset.")
  }
  e <- new.env(parent = emptyenv())
  load(infile, envir = e)
  if (exists("outdoor", envir = e)) x <- get("outdoor", envir = e)
  else if (exists("outdoor_final", envir = e)) x <- get("outdoor_final", envir = e)
  else stop("Could not find outdoor or outdoor_final in: ", infile)
  setDT(x)
  x[]
}

say("Loading: ", global_rdata)
outdoor <- load_outdoor_object(global_rdata)

a_col <- if (a_col_default %in% names(outdoor)) a_col_default else "cf1_a"
b_col <- if (b_col_default %in% names(outdoor)) b_col_default else "cf1_b"
lon_col <- intersect(lon_candidates, names(outdoor))[1]
if (is.na(lon_col)) stop("No longitude column found among: ",
                         paste(lon_candidates, collapse = ", "),
                         " — check names(outdoor) and add the right one.")

dt <- outdoor[region == "CONUS",
              c("station", "datetime_utc", a_col, b_col, lon_col),
              with = FALSE]
rm(outdoor); gc()

setnames(dt, c(a_col, b_col, lon_col), c("cf1_a", "cf1_b", "lon"))
dt[, `:=`(
  station = as.character(station),
  ts = as.POSIXct(datetime_utc, tz = "UTC"),
  cf1_a = as.numeric(cf1_a),
  cf1_b = as.numeric(cf1_b),
  lon = as.numeric(lon)
)]
dt[, datetime_utc := NULL]
dt <- dt[is.finite(cf1_a) & is.finite(cf1_b) &
           cf1_a >= 0 & cf1_b >= 0 & cf1_a <= PM_MAX & cf1_b <= PM_MAX]

setorder(dt, station, ts)
dt[, first_ts := ts[1L], by = station]
dt[, op_hours := as.numeric(difftime(ts, first_ts, units = "hours"))]
dt[, pm := (cf1_a + cf1_b) / 2]
dt[, `:=`(
  absdiff_cf1 = abs(cf1_a - cf1_b),
  pctdiff_cf1 = fifelse((cf1_a + cf1_b) > 0,
                        2 * abs(cf1_a - cf1_b) / (cf1_a + cf1_b), NA_real_)
)]
thr_dt <- dt[is.finite(pctdiff_cf1),
             .(thr85 = quantile(pctdiff_cf1, FLAG_PCTL, na.rm = TRUE, type = 7)),
             by = station]
dt <- thr_dt[dt, on = "station"]
dt[, flagged_cf1 := as.integer(absdiff_cf1 > FLAG_ABS_CUTOFF & pctdiff_cf1 > thr85)]
rm(thr_dt); gc()

# ------------------------------------------------------------
# 1) DAILY PANEL
# ------------------------------------------------------------
say("Building daily station panel")

dt[, date := as.Date(ts)]
daily <- dt[, .(
  n_hours = .N,
  pm_day = mean(pm, na.rm = TRUE),
  absdiff_day = mean(absdiff_cf1, na.rm = TRUE),
  n_flag = sum(flagged_cf1, na.rm = TRUE),
  age_days = min(op_hours) / 24,
  lon = first(lon)
), by = .(station, date)]
rm(dt); gc()

# --- complete each station's record to a FULL CALENDAR GRID so rolling
# --- windows are true calendar-day windows; data gaps become NA rows.
# --- (Rolling over observed rows only would silently stretch "7-day"
# --- windows across reporting gaps.)
say("Expanding to full calendar grid per station")
grid <- daily[, .(date = as.Date(
  seq(as.numeric(min(date)), as.numeric(max(date)), by = 1),
  origin = "1970-01-01")), by = station]
daily <- daily[grid, on = c("station", "date")]
setorder(daily, station, date)
daily[, lon := lon[!is.na(lon)][1L], by = station]

daily[, west := as.integer(lon < WEST_LON_CUT)]
daily[, fire_season := as.integer(month(date) %in% FIRE_MONTHS)]
# >>> HMS HOOK: replace fire_season with a station-day smoke-plume flag
# >>> from NOAA HMS polygons for the rigorous version.

# ------------------------------------------------------------
# 2) EXPOSURE (past 7 calendar days) AND OUTCOMES (future 7 calendar
#    days), rolling on the calendar grid; require minimum coverage
# ------------------------------------------------------------
say("Computing rolling exposure and future outcomes")

MIN_PAST_DAYS <- 5L   # >=5 of the 7 past calendar days observed

daily[, excess := pmax(pm_day - BURDEN_THRESHOLD, 0)]
daily[, obs_day := as.integer(is.finite(pm_day))]

roll_stats <- function(d) {
  # past window = the 7 calendar days before today (excluding today)
  burden <- shift(frollsum(d$excess, PAST_WINDOW, align = "right", na.rm = TRUE), 1L)
  spike_raw <- shift(
    frollapply(d$pm_day, PAST_WINDOW,
               function(z) if (all(is.na(z))) NA_real_ else max(z, na.rm = TRUE),
               align = "right"), 1L)
  spike <- pmax(spike_raw - BURDEN_THRESHOLD, 0)
  past_n <- shift(frollsum(d$obs_day, PAST_WINDOW, align = "right"), 1L)
  # fire-season share of the past 7 calendar days (defined for every day)
  fire_share <- shift(frollmean(d$fire_season, PAST_WINDOW, align = "right"), 1L)
  # future window = the 7 calendar days after today (excluding today)
  fut_abs <- shift(frollmean(d$absdiff_day, FUTURE_WINDOW, align = "left", na.rm = TRUE), -1L)
  fut_flag <- shift(frollsum(d$n_flag, FUTURE_WINDOW, align = "left", na.rm = TRUE), -1L)
  fut_hours <- shift(frollsum(d$n_hours, FUTURE_WINDOW, align = "left", na.rm = TRUE), -1L)
  fut_n <- shift(frollsum(d$obs_day, FUTURE_WINDOW, align = "left"), -1L)
  list(burden = burden, spike = spike, past_n = past_n, fire_share = fire_share,
       fut_abs = fut_abs, fut_flag = fut_flag, fut_hours = fut_hours, fut_n = fut_n)
}

daily[, c("burden", "spike", "past_n", "fire_share",
          "fut_abs", "fut_flag", "fut_hours", "fut_n") := roll_stats(.SD),
      by = station]

# keep rows where: today is observed, >=5/7 past calendar days observed,
# >=3/7 future calendar days observed
panel <- daily[
  is.finite(pm_day) & is.finite(age_days) &
    past_n >= MIN_PAST_DAYS & fut_n >= MIN_FUTURE_DAYS &
    is.finite(burden) & is.finite(spike) & is.finite(fut_abs)
]
# classify window by majority fire-season share
panel[, window_season := fifelse(fire_share >= 0.5, "Fire season (Jun-Oct)",
                                 "Non-fire season")]
panel[, region_ew := fifelse(west == 1, "Western CONUS", "Eastern CONUS")]
panel[, burden100 := burden / 100]
panel[, spike10 := spike / 10]
panel[, flag_rate_fut := fut_flag / pmax(fut_hours, 1)]

say("Panel rows: ", format(nrow(panel), big.mark = ","),
    " stations: ", format(uniqueN(panel$station), big.mark = ","))

# ------------------------------------------------------------
# 3) STRATIFIED JOINT MODELS + INTERACTION
# ------------------------------------------------------------
say("Fitting stratified models")

fit_stratum <- function(dd, outcome, label) {
  if (nrow(dd) < 5000L || uniqueN(dd$station) < 100L) return(data.table())
  fml <- as.formula(paste0(outcome,
    " ~ burden100 + spike10 + ns(pm_day, 3) + ns(age_days, 3) | station"))
  m <- tryCatch(
    feols(fml, data = dd, cluster = ~station),
    error = function(e) NULL)
  if (is.null(m)) return(data.table())
  ct <- coeftable(m)
  rbindlist(lapply(c("burden100", "spike10"), function(tm) {
    if (!(tm %in% rownames(ct))) return(data.table())
    data.table(stratum = label, outcome = outcome, term = tm,
               estimate = ct[tm, "Estimate"], se = ct[tm, "Std. Error"],
               ci_lwr = ct[tm, "Estimate"] - 1.96 * ct[tm, "Std. Error"],
               ci_upr = ct[tm, "Estimate"] + 1.96 * ct[tm, "Std. Error"],
               p = ct[tm, 4], n_obs = m$nobs,
               n_sensors = uniqueN(dd$station))
  }))
}

strata <- CJ(season = unique(panel$window_season),
             ew = unique(panel$region_ew))
res <- rbindlist(lapply(seq_len(nrow(strata)), function(i) {
  ss <- strata$season[i]; ee <- strata$ew[i]
  dd <- panel[window_season == ss & region_ew == ee]
  rbindlist(list(
    fit_stratum(dd, "fut_abs", paste(ee, "x", ss)),
    fit_stratum(dd, "flag_rate_fut", paste(ee, "x", ss))
  ))
}))

# formal interaction test (future |A-B|)
m_int <- feols(
  fut_abs ~ burden100 * fire_share + spike10 + ns(pm_day, 3) + ns(age_days, 3) |
    station,
  data = panel, cluster = ~station)
say("Interaction burden x fire-season share (future |A-B|):")
print(coeftable(m_int)[grep("burden", rownames(coeftable(m_int))), , drop = FALSE])

int_dt <- data.table(term = rownames(coeftable(m_int)), coeftable(m_int))
fwrite(res, file.path(outdir, "wildfire_stratified_models.csv"))
fwrite(int_dt, file.path(outdir, "wildfire_interaction_model.csv"))

# ------------------------------------------------------------
# 4) WORD TABLE
# ------------------------------------------------------------
ft_dt <- res[term == "burden100" & outcome == "fut_abs", .(
  Stratum = stratum,
  `Burden effect on future |A-B| (per +100 µg/m³-days, 95% CI)` =
    sprintf("%.2f (%.2f, %.2f)", estimate, ci_lwr, ci_upr),
  `N station-days` = format(n_obs, big.mark = ","),
  `N sensors` = format(n_sensors, big.mark = ",")
)]
ft <- autofit(fontsize(flextable(ft_dt), size = 9, part = "all"))
doc <- read_docx()
doc <- body_add_par(doc, paste0(
  "Table SX: Cumulative-burden associations with future 7-day |A-B| in ",
  "CONUS, stratified by fire season (Jun-Oct) and western (west of 100°W) ",
  "vs eastern CONUS. Joint burden + spike sensor fixed-effects models with ",
  "splines in current PM2.5 and operational age; sensor-clustered SEs."),
  style = "Normal")
doc <- body_add_flextable(doc, ft)
print(doc, target = file.path(outdir, "wildfire_stratification_table.docx"))
say("Done. Outputs in: ", outdir)
