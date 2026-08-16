# ==========================================================================
# 13_survivorship.R
#
# Survivorship and attrition: Kaplan-Meier, RMST, landmark Cox (10,000 h), pre-offline
# changes, and effective lifetime (cumulative flagged thresholds, >=100 h follow-up rule).
# => Figure 4; Tables S10-S18; Figures S39-S43.
#
# Extracted VERBATIM from purpleair_degradation2.Rmd (chunk #31 lines 13476-15228).
# Code is unmodified; see REPRODUCING.md for paths and run order.
# ==========================================================================

# ---- source chunk #31 (Rmd lines 13476-15228) ------------------------

# ============================================================
# SURVIVORSHIP / ATTRITION / EFFECTIVE LIFETIME ANALYSIS
# Robust landmark design
#
# Key updates:
#  1) Sensor lifetime is defined from PurpleAir observations only.
#  2) GEOS-CF is required only for GEOS-based metrics, not lifetime.
#  3) Cox models use a 10,000-hour landmark design.
#  4) Cox follow-up begins after 10,000 operational hours.
#  5) Cox models use robust Huber-White SEs.
#  6) Early-life metrics are measured before the survival risk period.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(mgcv)
  library(cowplot)
  library(scales)
  library(survival)
  library(patchwork)
})

# ------------------------------------------------------------
# 0) PATHS / SETTINGS
# ------------------------------------------------------------

global_rdata  <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_final.RData"
climate_rdata <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_with_koppen_timevarying.RData"
koppen_transition_csv <- "/Users/priyanka/Downloads/PurpleAir_EPA/station_koppen_transition_summary.csv"

outdir <- "/Users/priyanka/Downloads/PurpleAir_EPA/survivorship_landmark_UPDATED"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

early_cut <- 10000
bin_w <- 24L
censor_buffer_days <- 7L
pre_window_days <- 30L
baseline_start_days <- 180L
baseline_end_days <- 30L

min_stratum_events <- 10L
min_stratum_sensors <- 50L
min_at_risk_hazard <- 25L

hazard_bin_hours <- 24L * 30L
effective_thresholds <- c(0.05, 0.10)
min_post_hours_effective <- 100L

year_lines <- c(2, 4, 6) * 365.25 * 24

a_col_default <- "pm25_cf_1_a_hourly(ug/m3)_pm25_cf_1_a"
b_col_default <- "pm25_cf_1_b_hourly(ug/m3)_pm25_cf_1_b"
geos_col <- "geoscf_pm25"

fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  tolower(x)
}

say <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")),
      paste0(..., collapse = ""), "\n")
}

auto_ncol <- function(n_panels) {
  if (n_panels <= 6) return(n_panels)
  if (n_panels <= 12) return(4L)
  5L
}

auto_height <- function(n_panels, base = 3.2, per_row = 2.6) {
  nc <- auto_ncol(n_panels)
  nr <- ceiling(n_panels / nc)
  max(7, base + per_row * nr)
}

# ------------------------------------------------------------
# 1) LOAD HELPERS
# ------------------------------------------------------------

load_outdoor_object <- function(infile) {
  e <- new.env(parent = emptyenv())
  load(infile, envir = e)

  if (exists("outdoor", envir = e)) {
    x <- get("outdoor", envir = e)
  } else if (exists("outdoor_final", envir = e)) {
    x <- get("outdoor_final", envir = e)
  } else {
    stop("Could not find object named outdoor or outdoor_final in: ", infile)
  }

  setDT(x)
  x[]
}

get_single_zone_station_table <- function() {
  if (!file.exists(koppen_transition_csv)) {
    stop("Missing Köppen transition file: ", koppen_transition_csv)
  }

  kop <- fread(koppen_transition_csv)
  kop[, station := as.character(station)]

  kop <- kop[
    !is.na(first_koppen) &
      first_koppen != "" &
      n_koppen_codes == 1,
    .(
      station,
      climate_zone = as.character(first_koppen),
      climate_desc = as.character(first_koppen_desc)
    )
  ]

  kop[is.na(climate_desc) | climate_desc == "", climate_desc := climate_zone]
  kop[, climate_label := paste0(climate_zone, "\n", climate_desc)]
  unique(kop)
}

single_zone_stations <- get_single_zone_station_table()

# ------------------------------------------------------------
# 2) LOAD AND PREPARE HOURLY DATA
# ------------------------------------------------------------

make_surv_dt <- function(infile, require_climate = FALSE) {

  say("Loading: ", infile)

  outdoor <- load_outdoor_object(infile)

  a_col <- a_col_default
  b_col <- b_col_default

  if (!(a_col %in% names(outdoor)) && "cf1_a" %in% names(outdoor)) a_col <- "cf1_a"
  if (!(b_col %in% names(outdoor)) && "cf1_b" %in% names(outdoor)) b_col <- "cf1_b"

  needed <- c("station", "region", "datetime_utc", a_col, b_col)
  optional <- character()

  if (geos_col %in% names(outdoor)) optional <- c(optional, geos_col)

  miss <- setdiff(needed, names(outdoor))
  if (length(miss) > 0) {
    stop("Missing required columns in ", infile, ": ", paste(miss, collapse = ", "))
  }

  outdoor[, station := as.character(station)]

  if (require_climate) {
    outdoor <- outdoor[station %in% single_zone_stations$station]
    say("Rows after single-zone restriction: ", format(nrow(outdoor), big.mark = ","))
    say("Stations after single-zone restriction: ", format(uniqueN(outdoor$station), big.mark = ","))
  }

  dt <- outdoor[, c(needed, optional), with = FALSE]
  rm(outdoor)
  gc()

  setnames(dt, a_col, "cf1_a")
  setnames(dt, b_col, "cf1_b")

  dt[, `:=`(
    station = as.character(station),
    region = as.character(region),
    ts = as.POSIXct(datetime_utc, tz = "UTC"),
    cf1_a = as.numeric(cf1_a),
    cf1_b = as.numeric(cf1_b)
  )]

  if (geos_col %in% names(dt)) {
    dt[, geos_pm := as.numeric(get(geos_col))]
    dt[, (geos_col) := NULL]
  } else {
    dt[, geos_pm := NA_real_]
  }

  dt[, datetime_utc := NULL]

  # PurpleAir lifetime should not depend on GEOS availability.
  dt <- dt[
    !is.na(station) &
      !is.na(region) &
      !is.na(ts) &
      is.finite(cf1_a) &
      is.finite(cf1_b) &
      cf1_a >= 0 &
      cf1_b >= 0 &
      cf1_a <= 1500 &
      cf1_b <= 1500
  ]

  if (require_climate) {
    dt <- merge(dt, single_zone_stations, by = "station", all.x = FALSE)
  } else {
    dt[, `:=`(
      climate_zone = NA_character_,
      climate_desc = NA_character_,
      climate_label = NA_character_
    )]
  }

  setorder(dt, station, ts)

  dt[, first_ts := min(ts, na.rm = TRUE), by = station]
  dt[, op_hours := as.numeric(difftime(ts, first_ts, units = "hours"))]
  dt[, op_bin := floor(op_hours / bin_w) * bin_w]

  dt <- dt[is.finite(op_hours) & op_hours >= 0 & is.finite(op_bin)]

  dt[, `:=`(
    pa_cf1 = (cf1_a + cf1_b) / 2,
    ab_abs = abs(cf1_a - cf1_b),
    ab_bias = cf1_a - cf1_b,
    ab_rel = fifelse(
      (cf1_a + cf1_b) > 0,
      2 * abs(cf1_a - cf1_b) / (cf1_a + cf1_b),
      NA_real_
    )
  )]

  dt[, `:=`(
    pa_geos_signed = fifelse(is.finite(geos_pm), pa_cf1 - geos_pm, NA_real_),
    pa_geos_abs = fifelse(is.finite(geos_pm), abs(pa_cf1 - geos_pm), NA_real_)
  )]

  # Avoid using full-lifetime distribution to define early-life flag thresholds.
  # Here we use a fixed rule based on absolute and relative A-B disagreement.
  dt[, flagged_ab := as.integer(is.finite(ab_abs) & is.finite(ab_rel) & ab_abs > 5 & ab_rel > 0.20)]
  dt[is.na(flagged_ab), flagged_ab := 0L]

  say("Rows retained: ", format(nrow(dt), big.mark = ","))
  say("Stations retained: ", format(uniqueN(dt$station), big.mark = ","))

  dt[]
}

# ------------------------------------------------------------
# 3) METRIC MAP
# ------------------------------------------------------------

metric_map <- data.table(
  metric = c(
    "pa_geos_abs",
    "pa_geos_signed",
    "ab_abs",
    "ab_bias",
    "ab_rel",
    "flagged_ab"
  ),
  metric_label = c(
    "|PA CF1 − GEOS-CF|",
    "PA CF1 − GEOS-CF",
    "|A − B|",
    "A − B",
    "Relative A-B divergence",
    "Flagged rate"
  )
)

# ------------------------------------------------------------
# 4) STATION LIFE + EARLY METRICS + SENSOR-BIN DATA
# ------------------------------------------------------------

make_station_life <- function(dt, label = "dataset") {
  say("Building station lifetime table: ", label)

  study_end <- max(dt$ts, na.rm = TRUE)

  out <- dt[, .(
    region = region[1L],
    climate_zone = climate_zone[1L],
    climate_desc = climate_desc[1L],
    climate_label = climate_label[1L],
    first_ts = min(ts, na.rm = TRUE),
    last_ts = max(ts, na.rm = TRUE),
    max_op_hours = max(op_hours, na.rm = TRUE),
    max_op_years = max(op_hours, na.rm = TRUE) / (24 * 365.25),
    n_obs = .N
  ), by = station]

  out[, event_offline := fifelse(
    last_ts < study_end - censor_buffer_days * 24 * 3600,
    1L,
    0L
  )]

  out[, life_group := fifelse(
    max_op_hours < 10000, "<10k h",
    fifelse(
      max_op_hours < 30000, "10k-30k h",
      fifelse(max_op_hours < 60000, "30k-60k h", ">=60k h")
    )
  )]

  out[, life_group := factor(
    life_group,
    levels = c("<10k h", "10k-30k h", "30k-60k h", ">=60k h")
  )]

  out[]
}

add_early_metrics <- function(dt, station_life, label = "dataset") {
  say("Adding early-life metrics: ", label)

  early <- dt[
    op_hours <= early_cut,
    .(
      early_pa_geos_abs = mean(pa_geos_abs, na.rm = TRUE),
      early_pa_geos_signed = mean(pa_geos_signed, na.rm = TRUE),
      early_ab_abs = mean(ab_abs, na.rm = TRUE),
      early_ab_bias = mean(ab_bias, na.rm = TRUE),
      early_ab_rel = mean(ab_rel, na.rm = TRUE),
      early_flagged_ab = mean(flagged_ab, na.rm = TRUE),
      early_n_obs = .N
    ),
    by = station
  ]

  out <- merge(station_life, early, by = "station", all.x = TRUE)

  for (m in metric_map$metric) {
    v <- paste0("early_", m)
    vt <- paste0(v, "_trim")

    if (!v %in% names(out)) {
      out[, (vt) := NA_real_]
      next
    }

    cap <- quantile(out[[v]], 0.99, na.rm = TRUE, type = 7)
    if (!is.finite(cap)) cap <- Inf

    out[, (vt) := pmin(get(v), cap)]
  }

  # Landmark survival design:
  # include only sensors that survived beyond the early-life window.
  out[, landmark_eligible := is.finite(max_op_hours) & max_op_hours > early_cut]

  out[, `:=`(
    followup_after_landmark = pmax(max_op_hours - early_cut, 0),
    event_after_landmark = fifelse(landmark_eligible, event_offline, NA_integer_)
  )]

  out[]
}

make_sens_bin <- function(dt, station_life, label = "dataset") {
  say("Building station-bin summaries: ", label)

  out <- dt[, .(
    n_obs = .N,
    pa_geos_abs = mean(pa_geos_abs, na.rm = TRUE),
    pa_geos_signed = mean(pa_geos_signed, na.rm = TRUE),
    ab_abs = mean(ab_abs, na.rm = TRUE),
    ab_bias = mean(ab_bias, na.rm = TRUE),
    ab_rel = mean(ab_rel, na.rm = TRUE),
    flagged_ab = mean(flagged_ab, na.rm = TRUE),
    n_flagged_ab = sum(flagged_ab, na.rm = TRUE)
  ), by = .(
    station,
    region,
    climate_zone,
    climate_desc,
    climate_label,
    op_bin
  )]

  out <- merge(
    out,
    station_life[, .(station, life_group, max_op_hours, event_offline)],
    by = "station",
    all.x = TRUE
  )

  out[]
}

# ------------------------------------------------------------
# 5) CURVE HELPERS
# ------------------------------------------------------------

mean_ci_simple <- function(x) {
  x <- x[is.finite(x)]

  if (!length(x)) {
    return(list(mean = NA_real_, lwr = NA_real_, upr = NA_real_, n = 0L))
  }

  mu <- mean(x, na.rm = TRUE)
  se <- if (length(x) > 1) sd(x, na.rm = TRUE) / sqrt(length(x)) else 0

  list(mean = mu, lwr = mu - 1.96 * se, upr = mu + 1.96 * se, n = length(x))
}

build_curve <- function(d, metric, by_vars = NULL) {
  by_cols <- c(by_vars, "op_bin")

  d[
    is.finite(get(metric)),
    {
      z <- mean_ci_simple(get(metric))

      data.table(
        n_sensors = uniqueN(station),
        n_obs = sum(n_obs, na.rm = TRUE),
        mean = z$mean,
        lwr = z$lwr,
        upr = z$upr
      )
    },
    by = by_cols
  ][order(op_bin)]
}

gam_pred <- function(d) {
  dd <- d[is.finite(mean) & is.finite(op_bin) & n_sensors >= 5]

  if (nrow(dd) < 25 || uniqueN(dd$op_bin) < 20) return(data.table())

  fit <- tryCatch(
    gam(mean ~ s(op_bin, bs = "ps"), data = dd, weights = n_sensors, method = "REML"),
    error = function(e) NULL
  )

  if (is.null(fit)) return(data.table())

  nd <- data.table(op_bin = seq(min(dd$op_bin), max(dd$op_bin), by = 168))
  pr <- predict(fit, newdata = nd, se.fit = TRUE)

  nd[, `:=`(
    fit = pr$fit,
    lwr = pr$fit - 1.96 * pr$se.fit,
    upr = pr$fit + 1.96 * pr$se.fit
  )]

  nd[]
}

# ------------------------------------------------------------
# 6) LANDMARK COX MODEL HELPERS
# ------------------------------------------------------------

run_one_cox <- function(dd, metric, early_var, stratification, stratum,
                        adjust_region = FALSE) {

  dd <- as.data.table(dd)

  dd <- dd[
    landmark_eligible == TRUE &
      is.finite(followup_after_landmark) &
      followup_after_landmark > 0 &
      is.finite(get(early_var)) &
      !is.na(event_after_landmark)
  ]

  if (adjust_region) dd <- dd[!is.na(region)]

  n_here <- nrow(dd)
  ev_here <- sum(dd$event_after_landmark, na.rm = TRUE)

  if (n_here < min_stratum_sensors || ev_here < min_stratum_events) {
    return(data.table())
  }

  rhs <- if (adjust_region) paste0(early_var, " + region") else early_var

  fit <- tryCatch(
    coxph(
      as.formula(paste0("Surv(followup_after_landmark, event_after_landmark) ~ ", rhs)),
      data = dd,
      robust = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) return(data.table())

  sm <- summary(fit)
  if (!(early_var %in% rownames(sm$coefficients))) return(data.table())

  se_col <- if ("robust se" %in% colnames(sm$coefficients)) "robust se" else "se(coef)"

  data.table(
    stratification = stratification,
    stratum = stratum,
    metric = metric,
    term = early_var,
    beta = sm$coefficients[early_var, "coef"],
    se = sm$coefficients[early_var, se_col],
    hazard_ratio = sm$coefficients[early_var, "exp(coef)"],
    ci_lwr = sm$conf.int[early_var, "lower .95"],
    ci_upr = sm$conf.int[early_var, "upper .95"],
    p_value = sm$coefficients[early_var, "Pr(>|z|)"],
    n = sm$n,
    events = sm$nevent,
    landmark_hours = early_cut
  )
}

run_cox_metric_on_life <- function(station_life, metric, strat_vars = NULL,
                                   overall_adjust_region = FALSE,
                                   strat_label = NULL) {

  early_var <- paste0("early_", metric, "_trim")

  if (!(early_var %in% names(station_life))) {
    warning("Missing early variable: ", early_var)
    return(data.table())
  }

  if (is.null(strat_vars)) {
    dd <- station_life[is.finite(get(early_var))]
    return(run_one_cox(
      dd = dd,
      metric = metric,
      early_var = early_var,
      stratification = ifelse(overall_adjust_region, "Overall + region adjustment", "Overall"),
      stratum = "Overall",
      adjust_region = overall_adjust_region
    ))
  }

  dd0 <- station_life[is.finite(get(early_var))]
  dd0 <- dd0[complete.cases(dd0[, ..strat_vars])]

  if (nrow(dd0) == 0) return(data.table())

  dd0[, .group_id := do.call(paste, c(.SD, sep = " x ")), .SDcols = strat_vars]
  group_list <- split(dd0, by = ".group_id", keep.by = FALSE)

  rbindlist(lapply(names(group_list), function(g) {
    run_one_cox(
      dd = group_list[[g]],
      metric = metric,
      early_var = early_var,
      stratification = ifelse(is.null(strat_label), paste(strat_vars, collapse = " x "), strat_label),
      stratum = g,
      adjust_region = FALSE
    )
  }), fill = TRUE)
}

# ------------------------------------------------------------
# 7) PRE-FAILURE CHANGE
# ------------------------------------------------------------

compute_prefailure <- function(dt, station_life, label = "dataset") {

  say("Computing pre-failure changes: ", label)

  dt2 <- merge(
    dt,
    station_life[, .(station, last_ts_station = last_ts, event_offline)],
    by = "station",
    all.x = TRUE
  )

  dt2[, time_to_end_days := as.numeric(difftime(last_ts_station, ts, units = "days"))]

  dt2[, period_flag := fifelse(
    time_to_end_days <= pre_window_days,
    "last30",
    fifelse(
      time_to_end_days > baseline_end_days &
        time_to_end_days <= baseline_start_days,
      "baseline30_180",
      NA_character_
    )
  )]

  prefail_station <- dt2[
    event_offline == 1L & !is.na(period_flag),
    .(
      pa_geos_abs = mean(pa_geos_abs, na.rm = TRUE),
      pa_geos_signed = mean(pa_geos_signed, na.rm = TRUE),
      ab_abs = mean(ab_abs, na.rm = TRUE),
      ab_bias = mean(ab_bias, na.rm = TRUE),
      ab_rel = mean(ab_rel, na.rm = TRUE),
      flagged_ab = mean(flagged_ab, na.rm = TRUE),
      n_obs = .N
    ),
    by = .(station, region, climate_label, period_flag)
  ]

  prefail_long <- melt(
    prefail_station,
    id.vars = c("station", "region", "climate_label", "period_flag", "n_obs"),
    measure.vars = metric_map$metric,
    variable.name = "metric",
    value.name = "value"
  )

  prefail_long <- merge(
    prefail_long,
    metric_map[, .(metric, metric_label)],
    by = "metric",
    all.x = TRUE
  )

  prefail_wide <- dcast(
    prefail_long,
    station + region + climate_label + metric + metric_label ~ period_flag,
    value.var = "value"
  )

  if (!("last30" %in% names(prefail_wide))) prefail_wide[, last30 := NA_real_]
  if (!("baseline30_180" %in% names(prefail_wide))) prefail_wide[, baseline30_180 := NA_real_]

  prefail_wide[, change_last30_vs_baseline := last30 - baseline30_180]

  prefail_summary <- prefail_wide[
    is.finite(change_last30_vs_baseline),
    .(
      stations = .N,
      median_change = median(change_last30_vs_baseline, na.rm = TRUE),
      p25_change = quantile(change_last30_vs_baseline, 0.25, na.rm = TRUE),
      p75_change = quantile(change_last30_vs_baseline, 0.75, na.rm = TRUE),
      mean_change = mean(change_last30_vs_baseline, na.rm = TRUE)
    ),
    by = .(metric, metric_label, region)
  ]

  prefail_tests <- prefail_wide[
    is.finite(last30) & is.finite(baseline30_180),
    {
      wt <- tryCatch(wilcox.test(last30, baseline30_180, paired = TRUE), error = function(e) NULL)
      data.table(p_value = if (is.null(wt)) NA_real_ else wt$p.value)
    },
    by = .(metric, metric_label)
  ]

  prefail_tests[, p_value_fmt := fmt_p(p_value)]

  fwrite(prefail_wide, file.path(outdir, paste0("prefailure_station_changes_by_metric_", label, ".csv")))
  fwrite(prefail_summary, file.path(outdir, paste0("prefailure_summary_by_metric_region_", label, ".csv")))
  fwrite(prefail_tests, file.path(outdir, paste0("prefailure_wilcoxon_tests_by_metric_", label, ".csv")))

  rm(dt2, prefail_station, prefail_long)
  gc()

  list(wide = prefail_wide, summary = prefail_summary, tests = prefail_tests)
}

# ------------------------------------------------------------
# 8) RMST / HAZARD / EFFECTIVE LIFETIME HELPERS
# ------------------------------------------------------------

save_lifetime_summaries <- function(station_life, label, strat_var = "region") {

  say("Saving lifetime summaries: ", label, " by ", strat_var)

  dd <- station_life[
    !is.na(get(strat_var)) &
      is.finite(max_op_hours) &
      !is.na(event_offline)
  ]

  counts <- dd[, .(
    stations = .N,
    events = sum(event_offline, na.rm = TRUE),
    landmark_eligible = sum(landmark_eligible == TRUE, na.rm = TRUE)
  ), by = strat_var]

  fwrite(counts, file.path(outdir, paste0("lifetime_support_", label, "_by_", strat_var, ".csv")))

  keep <- counts[
    stations >= min_stratum_sensors &
      events >= min_stratum_events,
    get(strat_var)
  ]

  dd <- dd[get(strat_var) %in% keep]

  if (nrow(dd) < 50 || uniqueN(dd[[strat_var]]) < 2) return(invisible(NULL))

  form <- as.formula(paste0("Surv(max_op_hours, event_offline) ~ ", strat_var))
  km <- survfit(form, data = dd)

  km_tab <- as.data.table(summary(km)$table, keep.rownames = "stratum")
  km_tab[, group := sub(paste0("^", strat_var, "="), "", stratum)]
  fwrite(km_tab, file.path(outdir, paste0("km_summary_", label, "_by_", strat_var, ".csv")))

  tau <- quantile(dd$max_op_hours, 0.90, na.rm = TRUE, type = 7)

  rmst <- dd[, {
    fit <- tryCatch(survfit(Surv(max_op_hours, event_offline) ~ 1, data = .SD), error = function(e) NULL)
    if (is.null(fit)) {
      data.table()
    } else {
      tab <- summary(fit, rmean = tau)$table
      vals <- as.numeric(tab)
      nm <- names(tab)
      data.table(
        stations = .N,
        events = sum(event_offline, na.rm = TRUE),
        tau_hours = tau,
        rmst_hours = vals[grep("rmean", nm)[1]],
        rmst_se_hours = vals[grep("se\\(rmean\\)", nm)[1]],
        rmst_years = vals[grep("rmean", nm)[1]] / (24 * 365.25)
      )
    }
  }, by = strat_var]

  fwrite(rmst, file.path(outdir, paste0("rmst_", label, "_by_", strat_var, ".csv")))
  invisible(list(km = km_tab, rmst = rmst))
}

save_hazard_curves <- function(station_life, label, strat_var = "region") {

  say("Saving hazard curves: ", label, " by ", strat_var)

  dd <- station_life[
    !is.na(get(strat_var)) &
      is.finite(max_op_hours) &
      !is.na(event_offline)
  ]

  counts <- dd[, .(
    stations = .N,
    events = sum(event_offline, na.rm = TRUE)
  ), by = strat_var]

  keep <- counts[
    stations >= min_stratum_sensors &
      events >= min_stratum_events,
    get(strat_var)
  ]

  dd <- dd[get(strat_var) %in% keep]

  if (nrow(dd) < 50 || uniqueN(dd[[strat_var]]) < 1) return(invisible(NULL))

  max_follow <- ceiling(max(dd$max_op_hours, na.rm = TRUE) / hazard_bin_hours) * hazard_bin_hours
  breaks <- seq(0, max_follow, by = hazard_bin_hours)

  haz <- rbindlist(lapply(sort(unique(dd[[strat_var]])), function(g) {
    sub <- dd[get(strat_var) == g]

    rbindlist(lapply(seq_len(length(breaks) - 1L), function(i) {
      start <- breaks[i]
      end <- breaks[i + 1L]

      exposure <- pmax(0, pmin(sub$max_op_hours, end) - start)
      events <- sum(
        sub$event_offline == 1L &
          sub$max_op_hours >= start &
          sub$max_op_hours < end,
        na.rm = TRUE
      )
      at_risk <- sum(sub$max_op_hours >= start, na.rm = TRUE)

      data.table(
        group = g,
        age_start_hours = start,
        age_end_hours = end,
        age_mid_hours = (start + end) / 2,
        exposure_hours = sum(exposure, na.rm = TRUE),
        events = events,
        at_risk = at_risk
      )
    }))
  }), fill = TRUE)

  haz[, hazard_per_1000_sensor_hours := fifelse(
    exposure_hours > 0,
    1000 * events / exposure_hours,
    NA_real_
  )]

  haz[, age_mid_years := age_mid_hours / (24 * 365.25)]
  haz[at_risk < min_at_risk_hazard, hazard_per_1000_sensor_hours := NA_real_]

  fwrite(haz, file.path(outdir, paste0("hazard_", label, "_by_", strat_var, ".csv")))

  p <- ggplot(
    haz[is.finite(hazard_per_1000_sensor_hours)],
    aes(age_mid_years, hazard_per_1000_sensor_hours, color = group)
  ) +
    geom_line(linewidth = 0.8) +
    labs(
      title = paste0("Failure hazard over operational age: ", label),
      subtitle = paste0("Strata with ≥", min_stratum_sensors, " sensors and ≥", min_stratum_events, " events"),
      x = "Operational age years",
      y = "Hazard per 1,000 sensor-hours",
      color = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

  ggsave(
    file.path(outdir, paste0("hazard_", label, "_by_", strat_var, ".jpeg")),
    p,
    width = 11,
    height = 7,
    dpi = 500,
    bg = "white"
  )

  invisible(haz)
}

save_effective_lifetime <- function(sens_bin, station_life, label, strat_var = "region") {

  say("Saving effective lifetime: ", label, " by ", strat_var)

  flag_bin <- sens_bin[, .(
    flagged_n = sum(n_flagged_ab, na.rm = TRUE),
    total_n = sum(n_obs, na.rm = TRUE)
  ), by = .(station, op_bin)]

  setorder(flag_bin, station, op_bin)

  flag_bin[, `:=`(
    cum_flagged_n = cumsum(flagged_n),
    cum_total_n = cumsum(total_n),
    cum_flagged_prop = cumsum(flagged_n) / pmax(cumsum(total_n), 1)
  ), by = station]

  station_base <- station_life[, .(
    station,
    max_op_hours,
    event_offline,
    group_value = get(strat_var)
  )]

  out_all <- rbindlist(lapply(effective_thresholds, function(thr) {
    crossing <- flag_bin[
      cum_flagged_prop >= thr,
      .SD[1L],
      by = station
    ][, .(station, first_cross_hours = op_bin)]

    out <- merge(station_base, crossing, by = "station", all.x = TRUE)

    out[, enough_followup_after_cross := fifelse(
      is.finite(first_cross_hours),
      (max_op_hours - first_cross_hours) >= min_post_hours_effective,
      FALSE
    )]

    out[, `:=`(
      time_to_degrade_hours = fifelse(enough_followup_after_cross, first_cross_hours, max_op_hours),
      event_degrade = fifelse(enough_followup_after_cross, 1L, 0L),
      threshold = thr
    )]

    out
  }), fill = TRUE)

  fwrite(out_all, file.path(outdir, paste0("effective_lifetime_station_", label, "_by_", strat_var, ".csv")))

  km_eff <- out_all[
    !is.na(group_value),
    {
      fit <- tryCatch(
        survfit(Surv(time_to_degrade_hours, event_degrade) ~ 1, data = .SD),
        error = function(e) NULL
      )

      if (is.null(fit)) {
        data.table()
      } else {
        tab <- summary(fit)$table
        data.table(
          stations = .N,
          events = sum(event_degrade, na.rm = TRUE),
          median_effective_hours = tab["median"],
          median_effective_lcl_hours = tab["0.95LCL"],
          median_effective_ucl_hours = tab["0.95UCL"],
          median_effective_years = tab["median"] / (24 * 365.25)
        )
      }
    },
    by = .(threshold, group_value)
  ]

  fwrite(km_eff, file.path(outdir, paste0("effective_lifetime_km_", label, "_by_", strat_var, ".csv")))
  invisible(list(station = out_all, km = km_eff))
}

# ------------------------------------------------------------
# 9) PLOT HELPERS
# ------------------------------------------------------------

theme_fig <- theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8.5),
    strip.text = element_text(face = "bold", size = 8.5),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    legend.position = "bottom",
    plot.margin = margin(5, 5, 5, 5)
  )

save_km_stratified <- function(station_life, strat_var, outfile, title_text,
                               wrapped = TRUE) {

  dd <- station_life[
    !is.na(get(strat_var)) &
      is.finite(max_op_hours)
  ]

  counts <- dd[, .(
    stations = .N,
    events = sum(event_offline, na.rm = TRUE),
    landmark_eligible = sum(landmark_eligible == TRUE, na.rm = TRUE)
  ), by = strat_var]

  keep <- counts[
    stations >= min_stratum_sensors &
      events >= min_stratum_events,
    get(strat_var)
  ]

  dd <- dd[get(strat_var) %in% keep]

  fwrite(
    counts,
    file.path(outdir, paste0("km_support_", tools::file_path_sans_ext(outfile), ".csv"))
  )

  if (nrow(dd) < 50 || uniqueN(dd[[strat_var]]) < 2) {
    say("Skipping KM plot, insufficient strata: ", title_text)
    return(invisible(NULL))
  }

  form <- as.formula(paste0("Surv(max_op_hours, event_offline) ~ ", strat_var))
  km <- survfit(form, data = dd)

  km_dt <- data.table(
    time = km$time,
    surv = km$surv,
    strata = rep(names(km$strata), km$strata)
  )

  km_dt[, group := sub(paste0("^", strat_var, "="), "", strata)]
  km_dt[, time_years := time / (24 * 365.25)]

  if (wrapped) {
    n_groups <- uniqueN(km_dt$group)
    nc <- auto_ncol(n_groups)
    ht <- auto_height(n_groups, base = 3, per_row = 2.8)

    p <- ggplot(km_dt, aes(time_years, surv)) +
      geom_step(linewidth = 0.8, color = "black") +
      facet_wrap(~ group, ncol = nc) +
      labs(
        title = title_text,
        subtitle = paste0(
          "Showing strata with ≥", min_stratum_sensors,
          " sensors and ≥", min_stratum_events, " events"
        ),
        x = "Observed lifetime years",
        y = "Survival probability"
      ) +
      theme_bw(base_size = 10) +
      theme(plot.title = element_text(face = "bold"),
            strip.text = element_text(face = "bold", size = 7))

    ggsave(
      file.path(outdir, outfile),
      p,
      width = min(18, max(10, 3.2 * nc)),
      height = ht,
      dpi = 500,
      bg = "white"
    )

  } else {
    p <- ggplot(km_dt, aes(time_years, surv, color = group)) +
      geom_step(linewidth = 0.8) +
      labs(
        title = title_text,
        subtitle = paste0(
          "Showing strata with ≥", min_stratum_sensors,
          " sensors and ≥", min_stratum_events, " events"
        ),
        x = "Observed lifetime years",
        y = "Survival probability",
        color = NULL
      ) +
      theme_bw(base_size = 10) +
      theme(plot.title = element_text(face = "bold"),
            legend.position = "bottom")

    ggsave(file.path(outdir, outfile), p, width = 11, height = 7, dpi = 500, bg = "white")
  }

  invisible(p)
}

make_AG_figure <- function(metric_i, metric_label_i) {

  say("Making A-G figure for: ", metric_label_i)

  early_var <- paste0("early_", metric_i)
  early_trim_var <- paste0("early_", metric_i, "_trim")

  curve_all <- build_curve(sens_bin_global, metric_i, NULL)
  curve_all[, cohort := "All sensors"]

  curve_cohort <- build_curve(sens_bin_global, metric_i, "life_group")
  setnames(curve_cohort, "life_group", "cohort")

  curve_dt <- rbindlist(list(curve_all, curve_cohort), fill = TRUE)

  curve_dt[, cohort := factor(
    cohort,
    levels = c("All sensors", "<10k h", "10k-30k h", "30k-60k h", ">=60k h")
  )]

  gam_dt <- curve_dt[
    is.finite(mean),
    {
      gp <- gam_pred(.SD)
      if (nrow(gp) == 0) NULL else gp
    },
    by = cohort
  ]

  cap <- quantile(abs(curve_dt$mean), 0.995, na.rm = TRUE)
  if (!is.finite(cap) || cap <= 0) cap <- max(abs(curve_dt$mean), na.rm = TRUE)
  if (!is.finite(cap) || cap <= 0) cap <- 1

  curve_dt[, `:=`(
    mean_plot = pmax(pmin(mean, cap), -cap),
    lwr_plot = pmax(pmin(lwr, cap), -cap),
    upr_plot = pmax(pmin(upr, cap), -cap)
  )]

  if (nrow(gam_dt) > 0) gam_dt[, fit_plot := pmax(pmin(fit, cap), -cap)]

  pA <- ggplot(curve_dt, aes(op_bin, mean_plot, color = cohort, fill = cohort)) +
    geom_vline(xintercept = year_lines, linetype = 2, alpha = 0.35) +
    geom_ribbon(aes(ymin = lwr_plot, ymax = upr_plot), alpha = 0.10, color = NA) +
    geom_line(linewidth = 0.55) +
    geom_line(
      data = gam_dt,
      aes(op_bin, fit_plot, color = cohort),
      inherit.aes = FALSE,
      linewidth = 1
    ) +
    scale_x_continuous(labels = comma) +
    labs(
      title = paste0("A) ", metric_label_i, " by lifetime cohort"),
      subtitle = paste0("Y-axis capped at ", round(cap, 2)),
      x = "Operational time hours",
      y = metric_label_i,
      color = "Cohort",
      fill = "Cohort"
    ) +
    theme_fig

  composition <- sens_bin_global[, .(
    n_sensors = uniqueN(station)
  ), by = .(op_bin, life_group)]

  composition[, pct := 100 * n_sensors / sum(n_sensors), by = op_bin]

  pB <- ggplot(composition, aes(op_bin, pct, fill = life_group)) +
    geom_vline(xintercept = year_lines, linetype = 2, alpha = 0.35) +
    geom_area(position = "stack", alpha = 0.9) +
    scale_y_continuous(labels = percent_format(scale = 1)) +
    scale_x_continuous(labels = comma) +
    labs(
      title = "B) Active sensor composition over age",
      x = "Operational time hours",
      y = "Percent of active sensors",
      fill = "Eventual lifetime"
    ) +
    theme_fig

  trim_cap <- quantile(station_life_global[[early_var]], 0.99, na.rm = TRUE, type = 7)

  pC <- ggplot(
    station_life_global[is.finite(get(early_trim_var))],
    aes(x = life_group, y = get(early_trim_var))
  ) +
    geom_violin(fill = "grey88", color = "grey55", trim = TRUE, scale = "width") +
    geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white") +
    stat_summary(fun = median, geom = "point", size = 2, color = "red") +
    labs(
      title = "C) Early-life metric by eventual lifetime group",
      subtitle = paste0("Descriptive only; 99th percentile cap = ", round(trim_cap, 2)),
      x = "Eventual lifetime group",
      y = paste0("Early-life ", metric_label_i)
    ) +
    theme_fig

  scatter_dt <- station_life_global[
    landmark_eligible == TRUE &
      is.finite(get(early_var)) &
      is.finite(max_op_hours) &
      get(early_var) > 0
  ]

  pD <- ggplot(scatter_dt, aes(x = get(early_var), y = max_op_hours)) +
    geom_point(alpha = 0.08, size = 0.6) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "blue", se = TRUE) +
    scale_x_log10(labels = label_number()) +
    labs(
      title = "D) Early metric vs observed lifetime",
      subtitle = "Landmark-eligible sensors only; x-axis log scale",
      x = paste0("Early-life ", metric_label_i),
      y = "Observed lifetime hours"
    ) +
    theme_fig

  km_fit <- survfit(Surv(max_op_hours, event_offline) ~ life_group, data = station_life_global)

  km_dt <- data.table(
    time = km_fit$time,
    surv = km_fit$surv,
    strata = rep(names(km_fit$strata), km_fit$strata)
  )

  km_dt[, life_group := sub("^life_group=", "", strata)]

  pE <- ggplot(km_dt, aes(time / (24 * 365.25), surv, color = life_group)) +
    geom_step(linewidth = 0.9) +
    labs(
      title = "E) Kaplan-Meier survival by eventual lifetime cohort",
      subtitle = "Descriptive only; groups are defined by observed lifetime",
      x = "Observed lifetime years",
      y = "Survival probability",
      color = "Cohort"
    ) +
    theme_fig

  pf <- prefail_global$wide[
    metric == metric_i &
      is.finite(change_last30_vs_baseline)
  ]

  pF <- ggplot(pf, aes(region, change_last30_vs_baseline)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey45") +
    geom_violin(fill = "grey88", color = "grey55", trim = TRUE, scale = "width") +
    geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white") +
    stat_summary(fun = median, geom = "point", size = 2, color = "red") +
    labs(
      title = "F) Change before failure",
      subtitle = "Last 30 days minus 30–180 day baseline",
      x = "Region",
      y = paste0("Change in ", metric_label_i)
    ) +
    theme_fig

  cx <- cox_all[
    metric == metric_i &
      stratification == "Overall + region adjustment"
  ]

  pG <- ggplot(cx, aes(x = reorder(term, hazard_ratio), y = hazard_ratio)) +
    geom_hline(yintercept = 1, linetype = 2, color = "grey45") +
    geom_pointrange(aes(ymin = ci_lwr, ymax = ci_upr)) +
    coord_flip() +
    labs(
      title = "G) Landmark Cox hazard ratio",
      subtitle = "Failure risk after 10,000 h; early metric + region; robust SEs",
      x = NULL,
      y = "Hazard ratio 95% CI"
    ) +
    theme_fig

  row1 <- plot_grid(pA, pB, ncol = 2, align = "hv")
  row2 <- plot_grid(pC, pD, ncol = 2, align = "hv")
  row3 <- plot_grid(pE, pF, pG, ncol = 3, align = "hv")

  title_grob <- ggdraw() +
    draw_label(
      paste0("Survivorship and attrition analyses for ", metric_label_i),
      fontface = "bold",
      size = 16,
      x = 0,
      hjust = 0
    )

  fig <- plot_grid(
    title_grob,
    row1,
    row2,
    row3,
    ncol = 1,
    rel_heights = c(0.05, 1, 1, 1.05)
  )

  outfile <- file.path(
    outdir,
    paste0("survivorship_attrition_A_to_G_", safe_name(metric_label_i), ".jpeg")
  )

  ggsave(outfile, fig, width = 16, height = 17, dpi = 500, bg = "white")
  invisible(outfile)
}

# ------------------------------------------------------------
# 10) BUILD GLOBAL DATA
# ------------------------------------------------------------

say("Building global/regional dataset from outdoor_final.RData")

dt_global <- make_surv_dt(global_rdata, require_climate = FALSE)

station_life_global <- make_station_life(dt_global, "global_region")
station_life_global <- add_early_metrics(dt_global, station_life_global, "global_region")
sens_bin_global <- make_sens_bin(dt_global, station_life_global, "global_region")

prefail_global <- compute_prefailure(dt_global, station_life_global, "global_region")

rm(dt_global)
gc()

# ------------------------------------------------------------
# 11) BUILD CLIMATE DATA
# ------------------------------------------------------------

say("Building climate dataset from outdoor_with_koppen_timevarying.RData")

dt_climate <- make_surv_dt(climate_rdata, require_climate = TRUE)

station_life_climate <- make_station_life(dt_climate, "climate_single_zone")
station_life_climate <- add_early_metrics(dt_climate, station_life_climate, "climate_single_zone")
sens_bin_climate <- make_sens_bin(dt_climate, station_life_climate, "climate_single_zone")

prefail_climate <- compute_prefailure(dt_climate, station_life_climate, "climate_single_zone")

rm(dt_climate)
gc()

# ------------------------------------------------------------
# 12) COX MODELS
# ------------------------------------------------------------

say("Running landmark Cox models")

cox_all <- rbindlist(lapply(metric_map$metric, function(m) {
  rbindlist(list(
    run_cox_metric_on_life(
      station_life_global,
      m,
      strat_vars = NULL,
      overall_adjust_region = TRUE
    ),
    run_cox_metric_on_life(
      station_life_global,
      m,
      strat_vars = "region",
      strat_label = "region"
    ),
    run_cox_metric_on_life(
      station_life_climate,
      m,
      strat_vars = "climate_label",
      strat_label = "climate_label"
    ),
    run_cox_metric_on_life(
      station_life_climate,
      m,
      strat_vars = c("region", "climate_label"),
      strat_label = "region x climate_label"
    )
  ), fill = TRUE)
}), fill = TRUE)

cox_all <- merge(cox_all, metric_map, by = "metric", all.x = TRUE)
cox_all[, p_value_fmt := fmt_p(p_value)]
cox_all[, hr_95ci := sprintf("%.4f (%.4f, %.4f)", hazard_ratio, ci_lwr, ci_upr)]

fwrite(cox_all, file.path(outdir, "landmark_cox_models_global_region_climate_region_x_climate.csv"))

# ------------------------------------------------------------
# 13) A-G FIGURES
# ------------------------------------------------------------

say("Saving A-G survivorship figures")

for (i in seq_len(nrow(metric_map))) {
  make_AG_figure(metric_map$metric[i], metric_map$metric_label[i])
  gc()
}

# ------------------------------------------------------------
# 14) STRATIFIED KM CURVES
# ------------------------------------------------------------

station_life_climate[, region_climate := paste(region, climate_label, sep = " | ")]

save_km_stratified(
  station_life_global,
  "region",
  "km_survival_by_region_from_global_rdata.jpeg",
  "Kaplan-Meier survival curves by region",
  wrapped = FALSE
)

save_km_stratified(
  station_life_climate,
  "climate_label",
  "km_survival_by_climate_zone_from_climate_rdata_WRAPPED.jpeg",
  "Kaplan-Meier survival curves by climate zone",
  wrapped = TRUE
)

save_km_stratified(
  station_life_climate,
  "region_climate",
  "km_survival_by_region_x_climate_from_climate_rdata_WRAPPED.jpeg",
  "Kaplan-Meier survival curves by region × climate zone",
  wrapped = TRUE
)

for (rr in sort(unique(station_life_climate$region))) {
  save_km_stratified(
    station_life_climate[region == rr],
    "climate_label",
    paste0("km_survival_by_climate_zone_", safe_name(rr), "_WRAPPED.jpeg"),
    paste0("Kaplan-Meier survival by climate zone: ", rr),
    wrapped = TRUE
  )
}

# ------------------------------------------------------------
# 15) RMST / HAZARD / EFFECTIVE LIFETIME
# ------------------------------------------------------------

save_lifetime_summaries(station_life_global, "global_region", "region")
save_lifetime_summaries(station_life_climate, "single_zone_climate", "climate_label")
save_lifetime_summaries(station_life_climate, "region_x_climate", "region_climate")

save_hazard_curves(station_life_global, "global_region", "region")
save_hazard_curves(station_life_climate, "single_zone_climate", "climate_label")
save_hazard_curves(station_life_climate, "region_x_climate", "region_climate")

save_effective_lifetime(sens_bin_global, station_life_global, "global_region", "region")
save_effective_lifetime(sens_bin_climate, station_life_climate, "single_zone_climate", "climate_label")
save_effective_lifetime(sens_bin_climate, station_life_climate, "region_x_climate", "region_climate")

# ------------------------------------------------------------
# 16) SAVE FINAL TABLES
# ------------------------------------------------------------

fwrite(cox_all, file.path(outdir, "cox_all_metrics_landmark_global_region_climate.csv"))

fwrite(prefail_global$wide, file.path(outdir, "prefailure_station_metric_changes_global_region.csv"))
fwrite(prefail_global$summary, file.path(outdir, "prefailure_summary_global_region.csv"))
fwrite(prefail_global$tests, file.path(outdir, "prefailure_metric_tests_global_region.csv"))

fwrite(prefail_climate$wide, file.path(outdir, "prefailure_station_metric_changes_single_zone_climate.csv"))
fwrite(prefail_climate$summary, file.path(outdir, "prefailure_summary_single_zone_climate.csv"))
fwrite(prefail_climate$tests, file.path(outdir, "prefailure_metric_tests_single_zone_climate.csv"))

fwrite(station_life_global, file.path(outdir, "station_life_global_region_final.csv"))
fwrite(station_life_climate, file.path(outdir, "station_life_climate_final.csv"))

cat("\nDone. Outputs saved to:\n")
cat(outdir, "\n")

# ------------------------------------------------------------
# 17) WORD-READY TABLES FOR SURVIVORSHIP RESULTS
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(officer)
  library(flextable)
})

tabledir <- file.path(outdir, "word_ready_tables")
dir.create(tabledir, showWarnings = FALSE, recursive = TRUE)

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "", formatC(as.numeric(x), format = "f", digits = digits, big.mark = ","))
}

fmt_int <- function(x) {
  ifelse(is.na(x), "", formatC(as.integer(round(x)), format = "d", big.mark = ","))
}

fmt_hr <- function(hr, lwr, upr, digits = 2) {
  ifelse(
    is.na(hr) | is.na(lwr) | is.na(upr),
    "",
    paste0(
      formatC(hr, format = "f", digits = digits),
      " (",
      formatC(lwr, format = "f", digits = digits),
      ", ",
      formatC(upr, format = "f", digits = digits),
      ")"
    )
  )
}

clean_group <- function(x) {
  x <- gsub("\n", " ", x)
  x <- gsub("\\s+", " ", x)
  x
}

make_docx_table <- function(tab, title, footnote, outfile) {
  doc <- read_docx()
  doc <- body_add_par(doc, title, style = "heading 1")

  ft <- flextable(tab)
  ft <- autofit(ft)
  ft <- fontsize(ft, size = 9, part = "all")
  ft <- fontsize(ft, size = 8, part = "footer")
  ft <- bold(ft, part = "header")
  ft <- align(ft, align = "center", part = "all")
  ft <- align(ft, j = 1, align = "left", part = "body")
  ft <- theme_booktabs(ft)

  doc <- body_add_flextable(doc, ft)
  doc <- body_add_par(doc, footnote, style = "Normal")
  print(doc, target = outfile)
}

# ------------------------------------------------------------
# Table S. Survivorship status by region
# ------------------------------------------------------------

survival_status_region <- station_life_global[
  ,
  .(
    `Total sensors` = .N,
    `Went offline` = sum(event_offline == 1L, na.rm = TRUE),
    `Right-censored` = sum(event_offline == 0L, na.rm = TRUE),
    `Landmark eligible` = sum(landmark_eligible == TRUE, na.rm = TRUE),
    `Median lifetime (hours)` = median(max_op_hours, na.rm = TRUE),
    `Median lifetime (years)` = median(max_op_hours, na.rm = TRUE) / (24 * 365.25)
  ),
  by = .(Region = region)
][order(-`Total sensors`)]

survival_status_region[
  ,
  `:=`(
    `Went offline, %` = 100 * `Went offline` / `Total sensors`,
    `Right-censored, %` = 100 * `Right-censored` / `Total sensors`
  )
]

survival_status_region_word <- survival_status_region[
  ,
  .(
    Region,
    `Total sensors` = fmt_int(`Total sensors`),
    `Went offline` = paste0(fmt_int(`Went offline`), " (", fmt_num(`Went offline, %`, 1), "%)"),
    `Right-censored` = paste0(fmt_int(`Right-censored`), " (", fmt_num(`Right-censored, %`, 1), "%)"),
    `Landmark eligible` = fmt_int(`Landmark eligible`),
    `Median lifetime, h` = fmt_int(`Median lifetime (hours)`),
    `Median lifetime, y` = fmt_num(`Median lifetime (years)`, 2)
  )
]

fwrite(survival_status_region, file.path(tabledir, "survival_status_by_region.csv"))

make_docx_table(
  survival_status_region_word,
  "Table S. Survivorship status by region",
  "Offline status was defined using the final PurpleAir observation. Sensors whose last valid observation occurred within 7 days of the study end date were treated as right-censored. Landmark eligible sensors survived beyond 10,000 operational hours.",
  file.path(tabledir, "word_table_survival_status_by_region.docx")
)

# ------------------------------------------------------------
# Table S. Kaplan-Meier / RMST summaries by region
# ------------------------------------------------------------

km_region_file <- file.path(outdir, "km_summary_global_region_by_region.csv")
rmst_region_file <- file.path(outdir, "rmst_global_region_by_region.csv")

if (file.exists(km_region_file) & file.exists(rmst_region_file)) {

  km_region <- fread(km_region_file)
  rmst_region <- fread(rmst_region_file)

  km_region[, Region := clean_group(group)]
  rmst_region[, Region := clean_group(region)]

  km_rmst_region <- merge(
    km_region[, .(
      Region,
      records,
      events,
      rmean,
      `se(rmean)`,
      median,
      `0.95LCL`,
      `0.95UCL`
    )],
    rmst_region[, .(
      Region,
      tau_hours,
      rmst_hours,
      rmst_se_hours,
      rmst_years
    )],
    by = "Region",
    all = TRUE
  )

  km_rmst_region_word <- km_rmst_region[
    ,
    .(
      Region,
      `N sensors` = fmt_int(records),
      `Events` = fmt_int(events),
      `KM median, h` = fmt_int(median),
      `KM median 95% CI, h` = paste0(fmt_int(`0.95LCL`), ", ", fmt_int(`0.95UCL`)),
      `RMST, y` = fmt_num(rmst_years, 2),
      `RMST SE, h` = fmt_int(rmst_se_hours)
    )
  ]

  fwrite(km_rmst_region, file.path(tabledir, "km_rmst_summary_by_region.csv"))

  make_docx_table(
    km_rmst_region_word,
    "Table S. Kaplan-Meier and RMST survivorship summaries by region",
    "Kaplan-Meier medians and restricted mean survival times summarize observed PurpleAir operational lifetimes by region. RMST was calculated up to the analysis-specific follow-up horizon used in the survivorship code.",
    file.path(tabledir, "word_table_km_rmst_by_region.docx")
  )
}

# ------------------------------------------------------------
# Table S. Landmark Cox models, overall + region adjustment
# ------------------------------------------------------------

cox_overall_word <- cox_all[
  stratification == "Overall + region adjustment",
  .(
    Metric = metric_label,
    `Early-life predictor` = term,
    `Hazard ratio (95% CI)` = fmt_hr(hazard_ratio, ci_lwr, ci_upr, 3),
    `P value` = p_value_fmt,
    `N sensors` = fmt_int(n),
    `Events` = fmt_int(events),
    `Landmark, h` = fmt_int(landmark_hours)
  )
]

fwrite(
  cox_all[stratification == "Overall + region adjustment"],
  file.path(tabledir, "landmark_cox_overall_region_adjusted.csv")
)

make_docx_table(
  cox_overall_word,
  "Table S. Landmark Cox models for post-10,000-hour attrition",
  "Hazard ratios estimate the association between early-life degradation during the first 10,000 operational hours and subsequent sensor failure after the 10,000-hour landmark. Models adjust for region and use robust Huber-White standard errors.",
  file.path(tabledir, "word_table_landmark_cox_overall_region_adjusted.docx")
)

# ------------------------------------------------------------
# Table S. Landmark Cox models by region
# ------------------------------------------------------------

cox_region_word <- cox_all[
  stratification == "region",
  .(
    Region = clean_group(stratum),
    Metric = metric_label,
    `Hazard ratio (95% CI)` = fmt_hr(hazard_ratio, ci_lwr, ci_upr, 3),
    `P value` = p_value_fmt,
    `N sensors` = fmt_int(n),
    `Events` = fmt_int(events)
  )
][order(Metric, Region)]

fwrite(
  cox_all[stratification == "region"],
  file.path(tabledir, "landmark_cox_by_region.csv")
)

make_docx_table(
  cox_region_word,
  "Table S. Landmark Cox models by region",
  "Models were fit separately by region among sensors surviving beyond 10,000 operational hours. Hazard ratios represent post-landmark failure risk associated with each early-life degradation metric.",
  file.path(tabledir, "word_table_landmark_cox_by_region.docx")
)

# ------------------------------------------------------------
# Table S. Landmark Cox models by climate zone
# ------------------------------------------------------------

cox_climate_word <- cox_all[
  stratification == "climate_label",
  .(
    `Climate zone` = clean_group(stratum),
    Metric = metric_label,
    `Hazard ratio (95% CI)` = fmt_hr(hazard_ratio, ci_lwr, ci_upr, 3),
    `P value` = p_value_fmt,
    `N sensors` = fmt_int(n),
    `Events` = fmt_int(events)
  )
][order(Metric, `Climate zone`)]

fwrite(
  cox_all[stratification == "climate_label"],
  file.path(tabledir, "landmark_cox_by_climate_zone.csv")
)

make_docx_table(
  cox_climate_word,
  "Table S. Landmark Cox models by climate zone",
  "Models were fit separately by Köppen-Geiger climate zone among sensors surviving beyond 10,000 operational hours. Sparse climate zones should be interpreted cautiously.",
  file.path(tabledir, "word_table_landmark_cox_by_climate_zone.docx")
)

# ------------------------------------------------------------
# Table S. Lifetime summaries by climate zone
# ------------------------------------------------------------

km_climate_file <- file.path(outdir, "km_summary_single_zone_climate_by_climate_label.csv")

if (file.exists(km_climate_file)) {

  km_climate <- fread(km_climate_file)

  km_climate_word <- km_climate[
    ,
    .(
      `Climate zone` = clean_group(group),
      `N sensors` = fmt_int(records),
      `Events` = fmt_int(events),
      `KM median, h` = fmt_int(median),
      `KM median 95% CI, h` = paste0(fmt_int(`0.95LCL`), ", ", fmt_int(`0.95UCL`)),
      `RMST, h` = fmt_int(rmean),
      `RMST SE, h` = fmt_int(`se(rmean)`)
    )
  ][order(`Climate zone`)]

  fwrite(km_climate, file.path(tabledir, "km_lifetime_summary_by_climate_zone.csv"))

  make_docx_table(
    km_climate_word,
    "Table S. Kaplan-Meier survivorship summaries by climate zone",
    "Kaplan-Meier estimates summarize observed PurpleAir operational lifetimes among single-zone sensors. Only strata meeting the minimum support thresholds used in the survivorship code are shown.",
    file.path(tabledir, "word_table_km_by_climate_zone.docx")
  )
}

# ------------------------------------------------------------
# Table S. Lifetime summaries by region × climate
# ------------------------------------------------------------

km_region_climate_file <- file.path(outdir, "km_summary_region_x_climate_by_region_climate.csv")

if (file.exists(km_region_climate_file)) {

  km_region_climate <- fread(km_region_climate_file)

  km_region_climate_word <- km_region_climate[
    ,
    .(
      `Region × climate` = clean_group(group),
      `N sensors` = fmt_int(records),
      `Events` = fmt_int(events),
      `KM median, h` = fmt_int(median),
      `KM median 95% CI, h` = paste0(fmt_int(`0.95LCL`), ", ", fmt_int(`0.95UCL`)),
      `RMST, h` = fmt_int(rmean),
      `RMST SE, h` = fmt_int(`se(rmean)`)
    )
  ][order(`Region × climate`)]

  fwrite(km_region_climate, file.path(tabledir, "km_lifetime_summary_by_region_x_climate.csv"))

  make_docx_table(
    km_region_climate_word,
    "Table S. Kaplan-Meier survivorship summaries by region × climate zone",
    "Kaplan-Meier estimates summarize observed PurpleAir operational lifetimes by region and Köppen-Geiger climate zone. Only strata meeting minimum support thresholds are shown.",
    file.path(tabledir, "word_table_km_by_region_x_climate.docx")
  )
}

# ------------------------------------------------------------
# Table S. Pre-failure change summaries
# ------------------------------------------------------------

prefail_summary_word <- prefail_global$summary[
  ,
  .(
    Region = region,
    Metric = metric_label,
    `N sensors` = fmt_int(stations),
    `Median change` = fmt_num(median_change, 3),
    `IQR` = paste0(fmt_num(p25_change, 3), ", ", fmt_num(p75_change, 3)),
    `Mean change` = fmt_num(mean_change, 3)
  )
][order(Metric, Region)]

fwrite(prefail_global$summary, file.path(tabledir, "prefailure_summary_by_region.csv"))

make_docx_table(
  prefail_summary_word,
  "Table S. Pre-failure changes in degradation metrics by region",
  "Pre-failure change was defined as the mean degradation metric during the final 30 days of operation minus the mean during the 30–180 day pre-failure baseline period among sensors that went offline.",
  file.path(tabledir, "word_table_prefailure_changes_by_region.docx")
)

# ------------------------------------------------------------
# Table S. Effective lifetime based on cumulative flagged degradation
# ------------------------------------------------------------

eff_region_file <- file.path(outdir, "effective_lifetime_km_global_region_by_region.csv")

if (file.exists(eff_region_file)) {

  eff_region <- fread(eff_region_file)

  eff_region_word <- eff_region[
    ,
    .(
      Threshold = threshold,
      Region = clean_group(group_value),
      `N sensors` = fmt_int(stations),
      `Events` = fmt_int(events),
      `Median effective lifetime, h` = fmt_int(median_effective_hours),
      `Median effective lifetime, y` = fmt_num(median_effective_years, 2),
      `95% CI, h` = paste0(
        fmt_int(median_effective_lcl_hours),
        ", ",
        fmt_int(median_effective_ucl_hours)
      )
    )
  ][order(Threshold, Region)]

  fwrite(eff_region, file.path(tabledir, "effective_lifetime_by_region.csv"))

  make_docx_table(
    eff_region_word,
    "Table S. Effective lifetime based on cumulative flagged degradation by region",
    "Effective lifetime was defined as the operational age at which cumulative flagged degradation exceeded the specified threshold. Sensors not crossing the threshold were treated as censored.",
    file.path(tabledir, "word_table_effective_lifetime_by_region.docx")
  )
}

cat("\nWord-ready survivorship tables saved to:\n")
cat(tabledir, "\n")
