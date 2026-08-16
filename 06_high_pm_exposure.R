# ==========================================================================
# 14_highpm_exposure.R
#
# High-PM exposure and future degradation: 7/14/28-day windows x 35/55/100/150 ug/m3
# thresholds; frequency, cumulative burden and spike formulations; joint models.
# => Figure 5; Table S19; Figures S44-S49.
#
# Extracted VERBATIM from purpleair_degradation2.Rmd (chunks #32 lines 15231-17396, #33 lines 17399-18951).
# Code is unmodified; see REPRODUCING.md for paths and run order.
# ==========================================================================

# ---- source chunk #32 (Rmd lines 15231-17396) ------------------------

# ============================================================
# HIGH-PM EVENTS AND FUTURE SENSOR DEGRADATION
# MEMORY-EFFICIENT VERSION
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
  library(ggplot2)
  library(scales)
  library(splines)
})

# ------------------------------------------------------------
# 0) PATHS / SETTINGS
# ------------------------------------------------------------

global_rdata  <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_final.RData"
climate_rdata <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_with_koppen_timevarying.RData"
transition_csv <- "/Users/priyanka/Downloads/PurpleAir_EPA/station_koppen_transition_summary.csv"

OUT_DIR <- "/Users/priyanka/Downloads/PurpleAir_EPA/high_pm_future_degradation_combined"
MODEL_DIR <- file.path(OUT_DIR, "models")
FIG_DIR <- file.path(OUT_DIR, "figures")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(MODEL_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

BIN_W <- 24L
DROP_FIRST_DAYS <- 7L
DROP_H <- DROP_FIRST_DAYS * 24L

PM_THRESHOLDS <- c(35, 55, 100, 150)
PAST_WINDOWS <- c(7L, 14L, 28L)
FUTURE_WINDOW <- 7L

MIN_FUTURE_HOURS <- 24L * 3L
MIN_ROWS_MODEL <- 500L
MIN_STATIONS_MODEL <- 20L

# Critical memory fix: do not save sampled model data.
SAVE_SAMPLE_PER_COMBO <- 0L

# Critical disk/memory fix: do not save huge fixest model objects.
SAVE_MODEL_RDS <- FALSE

a_col_default <- "pm25_cf_1_a_hourly(ug/m3)_pm25_cf_1_a"
b_col_default <- "pm25_cf_1_b_hourly(ug/m3)_pm25_cf_1_b"
geos_col <- "geoscf_pm25"
koppen_col <- "koppen_code"

out_model_csv <- file.path(OUT_DIR, "high_pm_future_degradation_model_results.csv")
out_sample_csv <- file.path(OUT_DIR, "high_pm_future_degradation_analysis_SAMPLE.csv")
out_emp_csv <- file.path(OUT_DIR, "high_pm_future_degradation_empirical_summary.csv")

for (f in c(out_model_csv, out_sample_csv, out_emp_csv)) {
  if (file.exists(f)) file.remove(f)
}

say <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste0(..., collapse = ""), "\n")
}

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  tolower(x)
}

fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
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

get_single_zone_table <- function() {
  if (!file.exists(transition_csv)) {
    stop("Missing Köppen transition file: ", transition_csv)
  }

  z <- fread(transition_csv)
  z[, station := as.character(station)]

  z <- z[
    !is.na(first_koppen) &
      first_koppen != "" &
      n_koppen_codes == 1,
    .(
      station,
      climate_zone = as.character(first_koppen),
      climate_desc = as.character(first_koppen_desc)
    )
  ]

  z[is.na(climate_desc) | climate_desc == "", climate_desc := climate_zone]
  z[, climate_label := paste0(climate_zone, "\n", climate_desc)]

  unique(z)
}

single_zone <- get_single_zone_table()

# ------------------------------------------------------------
# 2) BUILD DAILY DATA
# ------------------------------------------------------------

make_daily_highpm <- function(infile, require_climate = FALSE) {

  say("Loading: ", infile)

  outdoor <- load_outdoor_object(infile)

  a_col <- a_col_default
  b_col <- b_col_default

  if (!(a_col %in% names(outdoor)) && "cf1_a" %in% names(outdoor)) a_col <- "cf1_a"
  if (!(b_col %in% names(outdoor)) && "cf1_b" %in% names(outdoor)) b_col <- "cf1_b"

  needed <- c("station", "region", "datetime_utc", a_col, b_col)

  has_geos <- geos_col %in% names(outdoor)
  if (has_geos) needed <- c(needed, geos_col)

  if (require_climate && koppen_col %in% names(outdoor)) {
    needed <- c(needed, koppen_col)
  }

  miss <- setdiff(needed, names(outdoor))
  if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))

  outdoor[, station := as.character(station)]

  if (require_climate) {
    outdoor <- outdoor[station %in% single_zone$station]
    say("Rows after single-zone restriction: ", format(nrow(outdoor), big.mark = ","))
    say("Stations after single-zone restriction: ", format(uniqueN(outdoor$station), big.mark = ","))
  }

  dt <- outdoor[, ..needed]
  rm(outdoor)
  gc()

  setnames(dt, a_col, "pm_a")
  setnames(dt, b_col, "pm_b")

  dt[, `:=`(
    station = as.character(station),
    region = as.character(region),
    ts = as.POSIXct(datetime_utc, tz = "UTC"),
    pm_a = as.numeric(pm_a),
    pm_b = as.numeric(pm_b)
  )]

  dt[, datetime_utc := NULL]

  if (has_geos) {
    dt[, geos_pm := as.numeric(get(geos_col))]
    if (geos_col %in% names(dt)) dt[, (geos_col) := NULL]
  } else {
    dt[, geos_pm := NA_real_]
  }

  dt <- dt[
    !is.na(station) &
      !is.na(region) &
      !is.na(ts) &
      is.finite(pm_a) &
      is.finite(pm_b)
  ]

  dt[, pm_cf1 := (pm_a + pm_b) / 2]

  dt[, abs_ab := abs(pm_a - pm_b)]
  dt[, ab_bias := pm_a - pm_b]
  dt[, ab_rel := fifelse(
    (pm_a + pm_b) > 0,
    2 * abs(pm_a - pm_b) / (pm_a + pm_b),
    NA_real_
  )]

  thr_dt <- dt[
    is.finite(ab_rel),
    .(thr85_ab_rel = quantile(ab_rel, 0.85, na.rm = TRUE, type = 7)),
    by = station
  ]

  dt <- thr_dt[dt, on = "station"]
  dt[, flagged := as.integer(abs_ab > 5 & ab_rel > thr85_ab_rel)]
  dt[is.na(flagged), flagged := 0L]

  dt[, `:=`(
    pa_geos_signed = fifelse(is.finite(geos_pm), pm_cf1 - geos_pm, NA_real_),
    pa_geos_abs = fifelse(is.finite(geos_pm), abs(pm_cf1 - geos_pm), NA_real_)
  )]

  setorder(dt, station, ts)
  dt[, first_ts := ts[1L], by = station]
  dt[, op_hours := as.numeric(difftime(ts, first_ts, units = "hours"))]
  dt[, op_bin := floor(op_hours / BIN_W) * BIN_W]

  dt <- dt[
    is.finite(op_bin) &
      op_bin >= DROP_H &
      is.finite(pm_cf1) &
      !is.na(flagged)
  ]

  daily <- dt[, .(
    region = region[1L],
    n = .N,
    n_flagged = sum(flagged, na.rm = TRUE),
    pm_mean = mean(pm_cf1, na.rm = TRUE),
    pm_max = max(pm_cf1, na.rm = TRUE),
    pm_sd = sd(pm_cf1, na.rm = TRUE),

    ab_abs = mean(abs_ab, na.rm = TRUE),
    ab_bias = mean(ab_bias, na.rm = TRUE),
    ab_rel = mean(ab_rel, na.rm = TRUE),
    flagged_rate = mean(flagged, na.rm = TRUE),

    pa_geos_abs = mean(pa_geos_abs, na.rm = TRUE),
    pa_geos_signed = mean(pa_geos_signed, na.rm = TRUE)
  ), by = .(station, op_bin)]

  rm(dt, thr_dt)
  gc()

  daily[, n_unflagged := n - n_flagged]

  if (require_climate) {
    daily <- merge(daily, single_zone, by = "station", all.x = FALSE)
  } else {
    daily[, `:=`(
      climate_zone = NA_character_,
      climate_desc = NA_character_,
      climate_label = NA_character_
    )]
  }

  setorder(daily, station, op_bin)

  say("Daily rows: ", format(nrow(daily), big.mark = ","))
  say("Daily stations: ", format(uniqueN(daily$station), big.mark = ","))

  daily[]
}

# ------------------------------------------------------------
# 3) ADD FUTURE OUTCOMES
# ------------------------------------------------------------

add_future_outcomes <- function(daily) {
  setorder(daily, station, op_bin)

  daily[, future_flagged_h := frollsum(
    shift(n_flagged, n = 1L, type = "lead"),
    n = FUTURE_WINDOW,
    align = "left",
    na.rm = TRUE
  ), by = station]

  daily[, future_total_h := frollsum(
    shift(n, n = 1L, type = "lead"),
    n = FUTURE_WINDOW,
    align = "left",
    na.rm = TRUE
  ), by = station]

  daily[, future_flagged_rate := fifelse(
    future_total_h > 0,
    future_flagged_h / future_total_h,
    NA_real_
  )]

  future_vars <- c(
    "ab_abs",
    "ab_bias",
    "ab_rel",
    "pa_geos_abs",
    "pa_geos_signed"
  )

  for (v in future_vars) {
    daily[, paste0("future_", v) := frollmean(
      shift(get(v), n = 1L, type = "lead"),
      n = FUTURE_WINDOW,
      align = "left",
      na.rm = TRUE
    ), by = station]
  }

  daily[]
}

# ------------------------------------------------------------
# 4) BUILD EXPOSURE DATA FOR ONE THRESHOLD/WINDOW
# ------------------------------------------------------------

safe_roll_max_lag <- function(x, n) {
  x_lag <- data.table::shift(x, n = 1L, type = "lag")

  if (length(x_lag) == 0L) return(numeric())

  if (requireNamespace("slider", quietly = TRUE)) {
    return(slider::slide_dbl(
      x_lag,
      .f = function(z) {
        z <- z[is.finite(z)]
        if (length(z) == 0L) NA_real_ else max(z)
      },
      .before = n - 1L,
      .complete = FALSE
    ))
  }

  out <- rep(NA_real_, length(x_lag))

  for (i in seq_along(x_lag)) {
    lo <- max(1L, i - n + 1L)
    z <- x_lag[lo:i]
    z <- z[is.finite(z)]
    out[i] <- if (length(z) == 0L) NA_real_ else max(z)
  }

  out
}

make_exposure_dt <- function(base_daily, thr, pw, analysis_group, group_var = "All") {

  dd <- copy(base_daily)
  setorder(dd, station, op_bin)

  dd[, high_pm_day := as.integer(pm_mean > thr)]
  dd[, excess_pm_day := pmax(pm_mean - thr, 0)]

  dd[, past_high_n := frollsum(
    shift(high_pm_day, n = 1L, type = "lag", fill = 0L),
    n = pw,
    align = "right",
    na.rm = TRUE
  ), by = station]

  dd[, past_obs_n := frollsum(
    shift(as.integer(is.finite(pm_mean)), n = 1L, type = "lag", fill = 0L),
    n = pw,
    align = "right",
    na.rm = TRUE
  ), by = station]

  dd[, past_high_prop := fifelse(
    past_obs_n > 0,
    past_high_n / past_obs_n,
    NA_real_
  )]

  dd[, past_excess_pm_sum := frollsum(
    shift(excess_pm_day, n = 1L, type = "lag", fill = 0),
    n = pw,
    align = "right",
    na.rm = TRUE
  ), by = station]

  dd[, past_excess_pm_mean := fifelse(
    past_obs_n > 0,
    past_excess_pm_sum / past_obs_n,
    NA_real_
  )]

  dd[, past_max_pm := safe_roll_max_lag(pm_max, pw), by = station]
  dd[, past_max_excess_pm := pmax(past_max_pm - thr, 0)]

  dd[, n_prev := shift(n, 1L), by = station]
  dd[, f_prev := shift(n_flagged, 1L), by = station]
  dd[, p_prev_smooth := (f_prev + 0.5) / (n_prev + 1.0)]
  dd[, offset_logit_prev := qlogis(p_prev_smooth)]

  dd[, step_h := op_bin - shift(op_bin), by = station]
  dd[, is_consecutive := is.finite(step_h) & step_h == BIN_W]

  out <- dd[
    is_consecutive &
      is.finite(offset_logit_prev) &
      past_obs_n >= max(3L, floor(pw / 2)) &
      is.finite(pm_mean) &
      is.finite(pm_sd) &
      is.finite(past_high_prop) &
      is.finite(past_excess_pm_sum) &
      is.finite(past_max_excess_pm) &
      is.finite(future_flagged_h) &
      is.finite(future_total_h) &
      future_total_h >= MIN_FUTURE_HOURS,
    .(
      station,
      region,
      climate_zone,
      climate_desc,
      climate_label,
      op_bin,
      future_flagged_h,
      future_total_h,
      future_unflagged_h = future_total_h - future_flagged_h,
      future_flagged_rate,
      future_ab_abs,
      future_ab_bias,
      future_ab_rel,
      future_pa_geos_abs,
      future_pa_geos_signed,
      past_high_prop,
      past_excess_pm_sum,
      past_excess_pm_mean,
      past_max_pm,
      past_max_excess_pm,
      pm_mean,
      pm_sd,
      offset_logit_prev,
      threshold = thr,
      past_window_days = pw,
      analysis_group = analysis_group,
      group_var = group_var
    )
  ][future_flagged_h >= 0 & future_total_h > future_flagged_h]

  rm(dd)
  gc()

  out[]
}

# ------------------------------------------------------------
# 5) MODEL HELPERS
# ------------------------------------------------------------

extract_fixest_term <- function(fit, term, scale_mult = 1, effect_type = "mean_difference") {
  s <- summary(fit, vcov = ~ station)
  cf <- s$coeftable

  if (!(term %in% rownames(cf))) return(NULL)

  beta <- cf[term, "Estimate"]
  se <- cf[term, "Std. Error"]
  pcol <- grep("^Pr", colnames(cf), value = TRUE)[1]
  pval <- cf[term, pcol]

  if (effect_type == "odds_ratio") {
    data.table(
      term = term,
      effect_type = effect_type,
      beta = beta,
      se = se,
      estimate = exp(beta * scale_mult),
      lwr = exp((beta - 1.96 * se) * scale_mult),
      upr = exp((beta + 1.96 * se) * scale_mult),
      p_value = pval
    )
  } else {
    data.table(
      term = term,
      effect_type = effect_type,
      beta = beta,
      se = se,
      estimate = beta * scale_mult,
      lwr = (beta - 1.96 * se) * scale_mult,
      upr = (beta + 1.96 * se) * scale_mult,
      p_value = pval
    )
  }
}

fit_one_model_set <- function(dt_fit, outcome_key, outcome_var, outcome_label) {

  if (outcome_key == "future_flagged") {
    keep <- is.finite(dt_fit$future_flagged_h) &
      is.finite(dt_fit$future_total_h) &
      dt_fit$future_total_h > dt_fit$future_flagged_h
  } else {
    keep <- is.finite(dt_fit[[outcome_var]])
  }

  if (!any(keep)) return(data.table())

  dd <- dt_fit[keep]

  if (nrow(dd) < MIN_ROWS_MODEL || uniqueN(dd$station) < MIN_STATIONS_MODEL) {
    rm(dd)
    gc()
    return(data.table())
  }

  common_rhs <- paste0(
    "ns(pm_mean, df = 4) + ",
    "ns(pm_sd, df = 3) + ",
    "ns(op_bin, df = 4) + ",
    "offset(offset_logit_prev)"
  )

  model_specs <- list(
    cumulative_only = list(
      rhs = paste0("past_excess_pm_sum + ", common_rhs),
      terms = "past_excess_pm_sum",
      scale = c(past_excess_pm_sum = 100),
      scale_label = c(past_excess_pm_sum = "+100 µg/m³-days")
    ),
    spike_only = list(
      rhs = paste0("past_max_excess_pm + ", common_rhs),
      terms = "past_max_excess_pm",
      scale = c(past_max_excess_pm = 10),
      scale_label = c(past_max_excess_pm = "+10 µg/m³ spike excess")
    ),
    frequency_only = list(
      rhs = paste0("past_high_prop + ", common_rhs),
      terms = "past_high_prop",
      scale = c(past_high_prop = 0.10),
      scale_label = c(past_high_prop = "10 percentage-point increase")
    ),
    joint_cumulative_spike = list(
      rhs = paste0("past_excess_pm_sum + past_max_excess_pm + ", common_rhs),
      terms = c("past_excess_pm_sum", "past_max_excess_pm"),
      scale = c(past_excess_pm_sum = 100, past_max_excess_pm = 10),
      scale_label = c(
        past_excess_pm_sum = "+100 µg/m³-days, adjusted for spike",
        past_max_excess_pm = "+10 µg/m³ spike excess, adjusted for cumulative burden"
      )
    )
  )

  out_rows <- list()
  idx <- 0L

  for (mn in names(model_specs)) {
    spec <- model_specs[[mn]]

    if (outcome_key == "future_flagged") {
      form <- as.formula(paste0(
        "future_flagged_h / future_total_h ~ ",
        spec$rhs,
        " | station"
      ))

      fit <- tryCatch(
        feglm(
          form,
          data = dd,
          family = "binomial",
          weights = ~ future_total_h,
          notes = TRUE
        ),
        error = function(e) NULL
      )

      effect_type <- "odds_ratio"

    } else {
      form <- as.formula(paste0(
        outcome_var,
        " ~ ",
        spec$rhs,
        " | station"
      ))

      fit <- tryCatch(
        feols(
          form,
          data = dd,
          weights = ~ future_total_h,
          notes = TRUE
        ),
        error = function(e) NULL
      )

      effect_type <- "mean_difference"
    }

    if (is.null(fit)) next

    for (tt in spec$terms) {
      res <- extract_fixest_term(
        fit,
        term = tt,
        scale_mult = spec$scale[[tt]],
        effect_type = effect_type
      )

      if (!is.null(res)) {
        idx <- idx + 1L
        out_rows[[idx]] <- cbind(
          data.table(
            model_name = mn,
            outcome_key = outcome_key,
            outcome_label = outcome_label,
            scale_unit = spec$scale_label[[tt]],
            n_rows = nrow(dd),
            n_stations = uniqueN(dd$station)
          ),
          res
        )
      }
    }

    rm(fit)
    gc()
  }

  rm(dd)
  gc()

  rbindlist(out_rows, use.names = TRUE, fill = TRUE)
}

# ------------------------------------------------------------
# 5B) EMPIRICAL SUMMARY HELPERS
# ------------------------------------------------------------

safe_empirical_cut <- function(x, probs = seq(0, 1, 0.1)) {
  x <- x[is.finite(x)]
  if (length(x) < 100L) return(NULL)

  qs <- unique(as.numeric(quantile(x, probs = probs, na.rm = TRUE, type = 7)))
  qs <- qs[is.finite(qs)]

  if (length(qs) < 3L) return(NULL)

  qs
}

save_empirical_summary <- function(dt_fit, exposure_var, exposure_label, outcome_var, outcome_label) {
  if (!(exposure_var %in% names(dt_fit)) || !(outcome_var %in% names(dt_fit))) {
    return(NULL)
  }

  dd <- dt_fit[
    is.finite(get(exposure_var)) &
      is.finite(get(outcome_var))
  ]

  if (nrow(dd) < 100L) return(NULL)

  qs <- safe_empirical_cut(dd[[exposure_var]])
  if (is.null(qs)) return(NULL)

  dd[, exposure_bin := cut(
    get(exposure_var),
    breaks = qs,
    include.lowest = TRUE,
    ordered_result = FALSE
  )]

  out <- dd[
    !is.na(exposure_bin),
    .(
      rows = .N,
      stations = uniqueN(station),
      mean_exposure = mean(get(exposure_var), na.rm = TRUE),
      outcome_mean = mean(get(outcome_var), na.rm = TRUE)
    ),
    by = exposure_bin
  ]

  if (nrow(out) == 0L) return(NULL)

  out[, `:=`(
    exposure_var = exposure_var,
    exposure_label = exposure_label,
    outcome_var = outcome_var,
    outcome_label = outcome_label,
    threshold = dt_fit$threshold[1],
    past_window_days = dt_fit$past_window_days[1],
    analysis_group = dt_fit$analysis_group[1],
    group_var = dt_fit$group_var[1]
  )]

  out[]
}

# ------------------------------------------------------------
# 6) RUN ANALYSES
# ------------------------------------------------------------

run_highpm_analysis <- function(daily, analysis_type) {

  say("Running analysis type: ", analysis_type)

  if (analysis_type == "global_region") {
    groups <- c("Global", sort(unique(na.omit(daily$region))))
    group_var_name <- "region"
  } else if (analysis_type == "climate") {
    groups <- sort(unique(na.omit(daily$climate_label)))
    group_var_name <- "climate_label"
  } else if (analysis_type == "region_climate") {
    daily[, region_climate := paste(region, climate_label, sep = " | ")]
    groups <- sort(unique(na.omit(daily$region_climate)))
    group_var_name <- "region_climate"
  } else {
    stop("Unknown analysis_type")
  }

  outcomes <- list(
    list(key = "future_flagged", var = "future_flagged_rate", label = "Future flagged probability"),
    list(key = "future_ab_abs", var = "future_ab_abs", label = "Future |A-B|"),
    list(key = "future_ab_bias", var = "future_ab_bias", label = "Future A-B"),
    list(key = "future_ab_rel", var = "future_ab_rel", label = "Future relative A-B"),
    list(key = "future_pa_geos_abs", var = "future_pa_geos_abs", label = "Future |PA-GEOS|"),
    list(key = "future_pa_geos_signed", var = "future_pa_geos_signed", label = "Future PA-GEOS")
  )

  for (grp in groups) {

    say("Group: ", grp)

    if (analysis_type == "global_region" && grp == "Global") {
      base_grp <- daily
    } else if (analysis_type == "global_region") {
      base_grp <- daily[region == grp]
    } else if (analysis_type == "climate") {
      base_grp <- daily[climate_label == grp]
    } else {
      base_grp <- daily[region_climate == grp]
    }

    if (nrow(base_grp) < MIN_ROWS_MODEL || uniqueN(base_grp$station) < MIN_STATIONS_MODEL) {
      say("  Skipping group: too few rows/stations.")
      next
    }

    setorder(base_grp, station, op_bin)

    for (thr in PM_THRESHOLDS) {
      for (pw in PAST_WINDOWS) {

        say("  Threshold=", thr, " window=", pw)

        dt_fit <- make_exposure_dt(
          base_daily = base_grp,
          thr = thr,
          pw = pw,
          analysis_group = grp,
          group_var = group_var_name
        )

        if (nrow(dt_fit) < MIN_ROWS_MODEL || uniqueN(dt_fit$station) < MIN_STATIONS_MODEL) {
          say("    Skipping threshold/window: too few fitted rows.")
          rm(dt_fit)
          gc()
          next
        }

        for (oo in outcomes) {
          if (!(oo$var %in% names(dt_fit))) next

          res <- fit_one_model_set(
            dt_fit = dt_fit,
            outcome_key = oo$key,
            outcome_var = oo$var,
            outcome_label = oo$label
          )

          if (nrow(res) > 0) {
            res[, `:=`(
              analysis_type = analysis_type,
              analysis_group = grp,
              group_var = group_var_name,
              threshold = thr,
              past_window_days = pw
            )]

            fwrite(
              res,
              out_model_csv,
              append = file.exists(out_model_csv)
            )
          }

          rm(res)
          gc()
        }

        emp_outcomes <- c(
          future_flagged_rate = "Future flagged probability",
          future_ab_abs = "Future |A-B|",
          future_ab_rel = "Future relative A-B"
        )

        for (ev in c("past_excess_pm_sum", "past_max_excess_pm", "past_high_prop")) {
          exp_lab <- switch(
            ev,
            past_excess_pm_sum = "Cumulative excess PM",
            past_max_excess_pm = "Sharp spike excess PM",
            past_high_prop = "Proportion of days above threshold"
          )

          for (ov in names(emp_outcomes)) {
            emp <- save_empirical_summary(
              dt_fit,
              exposure_var = ev,
              exposure_label = exp_lab,
              outcome_var = ov,
              outcome_label = emp_outcomes[[ov]]
            )

            if (!is.null(emp) && nrow(emp) > 0) {
              fwrite(
                emp,
                out_emp_csv,
                append = file.exists(out_emp_csv)
              )
            }

            rm(emp)
            gc()
          }
        }

        rm(dt_fit)
        gc()
      }
    }

    if (!(analysis_type == "global_region" && grp == "Global")) {
      rm(base_grp)
      gc()
    }
  }

  invisible(TRUE)
}

# ------------------------------------------------------------
# 7) BUILD DATASETS AND RUN
# ------------------------------------------------------------

say("Building global/regional daily dataset from outdoor_final.RData")
daily_global <- make_daily_highpm(global_rdata, require_climate = FALSE)
daily_global <- add_future_outcomes(daily_global)

run_highpm_analysis(daily_global, "global_region")

rm(daily_global)
gc()

say("Building climate daily dataset from outdoor_with_koppen_timevarying.RData")
daily_climate <- make_daily_highpm(climate_rdata, require_climate = TRUE)
daily_climate <- add_future_outcomes(daily_climate)

run_highpm_analysis(daily_climate, "climate")
run_highpm_analysis(daily_climate, "region_climate")

rm(daily_climate)
gc()

# ------------------------------------------------------------
# 8) CLEAN RESULTS
# ------------------------------------------------------------

model_dt <- fread(out_model_csv)

model_dt[, threshold_label := paste0(">", threshold, " µg/m³")]
model_dt[, window_label := paste0(past_window_days, "-day prior window")]
model_dt[, p_value_fmt := fmt_p(p_value)]

model_dt[, estimate_95ci := sprintf(
  "%.4f (%.4f, %.4f)",
  estimate,
  lwr,
  upr
)]

model_dt[, estimate_95ci_p := paste0(
  estimate_95ci,
  "; p = ",
  p_value_fmt
)]

fwrite(model_dt, file.path(OUT_DIR, "high_pm_future_degradation_model_results_clean.csv"))

# ------------------------------------------------------------
# 9) JOINT MODEL INTERPRETATION TABLE
# ------------------------------------------------------------

joint_parse <- model_dt[
  model_name == "joint_cumulative_spike" &
    term %in% c("past_excess_pm_sum", "past_max_excess_pm"),
  .(
    analysis_type,
    analysis_group,
    group_var,
    outcome_label,
    threshold,
    past_window_days,
    term,
    scale_unit,
    effect_type,
    estimate,
    lwr,
    upr,
    p_value,
    estimate_95ci_p,
    n_rows,
    n_stations
  )
]

joint_wide <- dcast(
  joint_parse,
  analysis_type + analysis_group + group_var + outcome_label + threshold + past_window_days ~ term,
  value.var = c("estimate", "p_value", "estimate_95ci_p")
)

fwrite(joint_parse, file.path(OUT_DIR, "joint_models_cumulative_vs_spike_long.csv"))
fwrite(joint_wide, file.path(OUT_DIR, "joint_models_cumulative_vs_spike_wide.csv"))

# ------------------------------------------------------------
# 10) MAIN FIGURES
# ------------------------------------------------------------

plot_dt <- copy(model_dt)

plot_dt[, exposure_class := fifelse(
  term == "past_excess_pm_sum", "Cumulative burden",
  fifelse(term == "past_max_excess_pm", "Sharp spike",
          fifelse(term == "past_high_prop", "Frequency above threshold", term))
)]

plot_dt[, model_name := factor(
  model_name,
  levels = c("frequency_only", "cumulative_only", "spike_only", "joint_cumulative_spike")
)]

plot_dt[, exposure_class := factor(
  exposure_class,
  levels = c("Frequency above threshold", "Cumulative burden", "Sharp spike")
)]

p_global <- ggplot(
  plot_dt[
    analysis_type == "global_region" &
      analysis_group == "Global" &
      outcome_key %in% c("future_flagged", "future_ab_abs", "future_ab_rel")
  ],
  aes(x = threshold_label, y = estimate, ymin = lwr, ymax = upr, color = exposure_class)
) +
  geom_hline(
    aes(yintercept = ifelse(effect_type == "odds_ratio", 1, 0)),
    linetype = 2,
    color = "grey50"
  ) +
  geom_pointrange(position = position_dodge(width = 0.55)) +
  facet_grid(outcome_label ~ model_name + window_label, scales = "free_y") +
  theme_bw(base_size = 10) +
  labs(
    title = "Global associations: high-PM exposure and future degradation",
    subtitle = "Separate and joint models parse cumulative exposure from sharp spikes",
    x = "PM2.5 threshold",
    y = "Effect estimate",
    color = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y = element_text(angle = 0)
  )

ggsave(
  file.path(FIG_DIR, "fig_global_cumulative_vs_spike.jpeg"),
  p_global,
  width = 18,
  height = 10,
  dpi = 500
)

p_region <- ggplot(
  plot_dt[
    analysis_type == "global_region" &
      analysis_group != "Global" &
      outcome_key %in% c("future_flagged", "future_ab_abs", "future_ab_rel") &
      model_name == "joint_cumulative_spike" &
      past_window_days == 28
  ],
  aes(x = threshold_label, y = estimate, ymin = lwr, ymax = upr, color = exposure_class)
) +
  geom_hline(
    aes(yintercept = ifelse(effect_type == "odds_ratio", 1, 0)),
    linetype = 2,
    color = "grey50"
  ) +
  geom_pointrange(position = position_dodge(width = 0.55)) +
  facet_grid(outcome_label ~ analysis_group, scales = "free_y") +
  theme_bw(base_size = 10) +
  labs(
    title = "Regional heterogeneity: cumulative burden vs sharp PM spikes",
    subtitle = "Joint models, 28-day exposure window",
    x = "PM2.5 threshold",
    y = "Effect estimate",
    color = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y = element_text(angle = 0)
  )

ggsave(
  file.path(FIG_DIR, "fig_region_joint_cumulative_vs_spike_28day.jpeg"),
  p_region,
  width = 16,
  height = 9,
  dpi = 500
)

p_climate <- ggplot(
  plot_dt[
    analysis_type == "climate" &
      outcome_key %in% c("future_flagged", "future_ab_abs", "future_ab_rel") &
      model_name == "joint_cumulative_spike" &
      past_window_days == 28
  ],
  aes(x = threshold_label, y = estimate, ymin = lwr, ymax = upr, color = exposure_class)
) +
  geom_hline(
    aes(yintercept = ifelse(effect_type == "odds_ratio", 1, 0)),
    linetype = 2,
    color = "grey50"
  ) +
  geom_pointrange(position = position_dodge(width = 0.55)) +
  facet_grid(outcome_label ~ analysis_group, scales = "free_y") +
  theme_bw(base_size = 9) +
  labs(
    title = "Climate-zone heterogeneity: cumulative burden vs sharp PM spikes",
    subtitle = "Joint models, 28-day exposure window; valid single-zone Köppen stations only",
    x = "PM2.5 threshold",
    y = "Effect estimate",
    color = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(size = 6),
    strip.text.y = element_text(angle = 0)
  )

ggsave(
  file.path(FIG_DIR, "fig_climate_joint_cumulative_vs_spike_28day.jpeg"),
  p_climate,
  width = 22,
  height = 10,
  dpi = 500,
  limitsize = FALSE
)

# ------------------------------------------------------------
# 11) IMPROVED EMPIRICAL CURVES
#     Includes internal degradation + PA-GEOS external divergence
# ------------------------------------------------------------

if (file.exists(out_emp_csv)) {
  emp_dt <- fread(out_emp_csv)

  emp_dt[, threshold_label := paste0(">", threshold, " µg/m³")]
  emp_dt[, window_label := paste0(past_window_days, "-day prior window")]

  MIN_EMP_ROWS <- 500L
  MIN_EMP_STATIONS <- 20L

  emp_dt[, outcome_domain := fifelse(
    outcome_var %in% c("future_pa_geos_abs", "future_pa_geos_signed"),
    "External-reference disagreement",
    "Internal sensor disagreement"
  )]

  emp_dt[, outcome_label := fifelse(
    outcome_var == "future_pa_geos_abs",
    "Future |PA-GEOS|",
    fifelse(
      outcome_var == "future_pa_geos_signed",
      "Future PA-GEOS",
      outcome_label
    )
  )]

  emp_plot <- emp_dt[
    analysis_group == "Global" &
      past_window_days == 28 &
      threshold %in% c(35, 100, 150) &
      rows >= MIN_EMP_ROWS &
      stations >= MIN_EMP_STATIONS
  ]

  if (nrow(emp_plot) > 0L) {
    emp_plot[, exposure_label := factor(
      exposure_label,
      levels = c(
        "Cumulative excess PM",
        "Sharp spike excess PM",
        "Proportion of days above threshold"
      )
    )]

    p_emp_clean <- ggplot(
      emp_plot,
      aes(x = mean_exposure, y = outcome_mean)
    ) +
      geom_point(aes(size = stations), alpha = 0.8) +
      geom_line(linewidth = 0.6) +
      facet_grid(
        outcome_label ~ exposure_label + threshold_label,
        scales = "free"
      ) +
      scale_size_continuous(range = c(1.2, 3.2)) +
      theme_bw(base_size = 10) +
      labs(
        title = "Empirical global future degradation by prior high-PM exposure",
        subtitle = paste0(
          "28-day prior window; bins retained only if ≥",
          MIN_EMP_ROWS, " rows and ≥", MIN_EMP_STATIONS, " sensors"
        ),
        x = "Mean exposure in bin",
        y = "Future outcome mean",
        size = "Sensors"
      ) +
      theme(
        legend.position = "bottom",
        strip.text.x = element_text(size = 8),
        strip.text.y = element_text(size = 8)
      )

    ggsave(
      file.path(FIG_DIR, "fig_empirical_global_28day_clean_with_pa_geos.jpeg"),
      p_emp_clean,
      width = 16,
      height = 11,
      dpi = 500
    )
  } else {
    warning("No empirical bins met support thresholds for the cleaned global empirical plot.")
  }

  fwrite(
    emp_plot,
    file.path(OUT_DIR, "high_pm_future_degradation_empirical_summary_clean_global_28day_with_pa_geos.csv")
  )
}

# ------------------------------------------------------------
# 12) CLEAN / FILTERED MODEL TABLES
# ------------------------------------------------------------

model_dt <- fread(file.path(OUT_DIR, "high_pm_future_degradation_model_results_clean.csv"))

model_dt[, threshold_label := paste0(">", threshold, " µg/m³")]
model_dt[, window_label := paste0(past_window_days, "-day prior window")]

model_dt[, outcome_label := fifelse(
  outcome_key == "future_pa_geos_abs",
  "Future |PA-GEOS|",
  fifelse(
    outcome_key == "future_pa_geos_signed",
    "Future PA-GEOS",
    outcome_label
  )
)]

model_dt[, degradation_domain := fifelse(
  outcome_key %in% c("future_pa_geos_abs", "future_pa_geos_signed"),
  "External-reference disagreement",
  "Internal sensor disagreement"
)]

model_dt[, exposure_class := fifelse(
  term == "past_excess_pm_sum", "Cumulative burden",
  fifelse(
    term == "past_max_excess_pm", "Sharp spike",
    fifelse(term == "past_high_prop", "Frequency above threshold", term)
  )
)]

model_dt[, exposure_class := factor(
  exposure_class,
  levels = c("Frequency above threshold", "Cumulative burden", "Sharp spike")
)]

model_dt[, model_name := factor(
  model_name,
  levels = c(
    "frequency_only",
    "cumulative_only",
    "spike_only",
    "joint_cumulative_spike"
  )
)]

MIN_PLOT_ROWS <- 1000L
MIN_PLOT_STATIONS <- 20L

outcomes_all <- c(
  "future_flagged",
  "future_ab_abs",
  "future_ab_rel",
  "future_pa_geos_abs",
  "future_pa_geos_signed"
)

plot_dt_clean <- model_dt[
  n_rows >= MIN_PLOT_ROWS &
    n_stations >= MIN_PLOT_STATIONS
]

fwrite(
  plot_dt_clean,
  file.path(OUT_DIR, "high_pm_future_degradation_model_results_clean_supported_with_pa_geos.csv")
)

# ------------------------------------------------------------
# 13) CLEAN GLOBAL FIGURE
# ------------------------------------------------------------

p_global_clean <- ggplot(
  plot_dt_clean[
    analysis_type == "global_region" &
      analysis_group == "Global" &
      outcome_key %in% outcomes_all
  ],
  aes(
    x = threshold_label,
    y = estimate,
    ymin = lwr,
    ymax = upr,
    color = exposure_class
  )
) +
  geom_hline(
    aes(yintercept = ifelse(effect_type == "odds_ratio", 1, 0)),
    linetype = 2,
    color = "grey50"
  ) +
  geom_pointrange(position = position_dodge(width = 0.55)) +
  facet_grid(outcome_label ~ model_name + window_label, scales = "free_y") +
  theme_bw(base_size = 9.5) +
  labs(
    title = "Global associations: high-PM exposure and future sensor degradation",
    subtitle = paste0(
      "Internal A/B metrics plus external PA-GEOS disagreement; supported panels only: ≥",
      MIN_PLOT_ROWS, " rows and ≥", MIN_PLOT_STATIONS, " sensors"
    ),
    x = "PM2.5 threshold",
    y = "Effect estimate",
    color = NULL
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y = element_text(angle = 0)
  )

ggsave(
  file.path(FIG_DIR, "fig_global_cumulative_vs_spike_supported_with_pa_geos.jpeg"),
  p_global_clean,
  width = 18,
  height = 13,
  dpi = 500
)

# ------------------------------------------------------------
# 14) CLEAN REGIONAL JOINT MODEL FIGURE
# ------------------------------------------------------------

p_region_clean <- ggplot(
  plot_dt_clean[
    analysis_type == "global_region" &
      analysis_group != "Global" &
      outcome_key %in% outcomes_all &
      model_name == "joint_cumulative_spike" &
      past_window_days == 28 &
      term %in% c("past_excess_pm_sum", "past_max_excess_pm")
  ],
  aes(
    x = threshold_label,
    y = estimate,
    ymin = lwr,
    ymax = upr,
    color = exposure_class
  )
) +
  geom_hline(
    aes(yintercept = ifelse(effect_type == "odds_ratio", 1, 0)),
    linetype = 2,
    color = "grey50"
  ) +
  geom_pointrange(position = position_dodge(width = 0.55)) +
  facet_grid(outcome_label ~ analysis_group, scales = "free_y") +
  theme_bw(base_size = 9.5) +
  labs(
    title = "Regional heterogeneity: cumulative burden vs sharp PM spikes",
    subtitle = "Joint models, 28-day exposure window; includes PA-GEOS external disagreement",
    x = "PM2.5 threshold",
    y = "Effect estimate",
    color = NULL
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y = element_text(angle = 0)
  )

ggsave(
  file.path(FIG_DIR, "fig_region_joint_cumulative_vs_spike_28day_supported_with_pa_geos.jpeg"),
  p_region_clean,
  width = 16,
  height = 12,
  dpi = 500
)

# ------------------------------------------------------------
# 15) CLEAN CLIMATE JOINT MODEL FIGURE
# ------------------------------------------------------------

climate_support <- plot_dt_clean[
  analysis_type == "climate" &
    model_name == "joint_cumulative_spike" &
    past_window_days == 28,
  .(
    min_rows = min(n_rows, na.rm = TRUE),
    min_stations = min(n_stations, na.rm = TRUE),
    n_terms = .N
  ),
  by = analysis_group
][
  min_rows >= MIN_PLOT_ROWS &
    min_stations >= MIN_PLOT_STATIONS
]

climate_keep <- climate_support$analysis_group

p_climate_clean <- ggplot(
  plot_dt_clean[
    analysis_type == "climate" &
      analysis_group %in% climate_keep &
      outcome_key %in% outcomes_all &
      model_name == "joint_cumulative_spike" &
      past_window_days == 28 &
      term %in% c("past_excess_pm_sum", "past_max_excess_pm")
  ],
  aes(
    x = threshold_label,
    y = estimate,
    ymin = lwr,
    ymax = upr,
    color = exposure_class
  )
) +
  geom_hline(
    aes(yintercept = ifelse(effect_type == "odds_ratio", 1, 0)),
    linetype = 2,
    color = "grey50"
  ) +
  geom_pointrange(position = position_dodge(width = 0.55)) +
  facet_grid(outcome_label ~ analysis_group, scales = "free_y") +
  theme_bw(base_size = 8.2) +
  labs(
    title = "Climate-zone heterogeneity: cumulative burden vs sharp PM spikes",
    subtitle = paste0(
      "Joint models, 28-day exposure window; includes PA-GEOS; panels require ≥",
      MIN_PLOT_ROWS, " rows and ≥", MIN_PLOT_STATIONS, " sensors"
    ),
    x = "PM2.5 threshold",
    y = "Effect estimate",
    color = NULL
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5.8),
    strip.text.x = element_text(size = 5.2),
    strip.text.y = element_text(angle = 0, size = 7)
  )

ggsave(
  file.path(FIG_DIR, "fig_climate_joint_cumulative_vs_spike_28day_supported_with_pa_geos.jpeg"),
  p_climate_clean,
  width = max(16, 0.9 * length(climate_keep)),
  height = 12,
  dpi = 500,
  limitsize = FALSE
)

fwrite(
  climate_support,
  file.path(OUT_DIR, "climate_panels_retained_for_supported_plot_with_pa_geos.csv")
)

# ------------------------------------------------------------
# 16) SUMMARY TABLE: CUMULATIVE VS SPIKE IN JOINT MODELS
# ------------------------------------------------------------

joint_summary <- model_dt[
  model_name == "joint_cumulative_spike" &
    term %in% c("past_excess_pm_sum", "past_max_excess_pm"),
  .(
    analysis_type,
    analysis_group,
    group_var,
    degradation_domain,
    outcome_key,
    outcome_label,
    threshold,
    past_window_days,
    exposure_class,
    term,
    effect_type,
    estimate,
    lwr,
    upr,
    p_value,
    n_rows,
    n_stations
  )
]

joint_summary[, significant := p_value < 0.05]

joint_direction <- joint_summary[
  ,
  .(
    n_models = .N,
    n_significant = sum(significant, na.rm = TRUE),
    pct_significant = mean(significant, na.rm = TRUE),
    median_estimate = median(estimate, na.rm = TRUE),
    median_lwr = median(lwr, na.rm = TRUE),
    median_upr = median(upr, na.rm = TRUE)
  ),
  by = .(
    analysis_type,
    degradation_domain,
    outcome_label,
    past_window_days,
    exposure_class,
    effect_type
  )
]

fwrite(
  joint_summary,
  file.path(OUT_DIR, "joint_models_cumulative_vs_spike_summary_long_with_pa_geos.csv")
)

fwrite(
  joint_direction,
  file.path(OUT_DIR, "joint_models_cumulative_vs_spike_direction_summary_with_pa_geos.csv")
)

# ------------------------------------------------------------
# 17) HEATMAP OF JOINT MODEL SIGNALS
# ------------------------------------------------------------

heat_dt <- joint_summary[
  analysis_type %in% c("global_region", "climate") &
    outcome_key %in% outcomes_all
]

heat_dt[, signed_logp := -log10(pmax(p_value, 1e-300)) * sign(estimate)]

p_heat <- ggplot(
  heat_dt[
    past_window_days == 28 &
      threshold %in% c(35, 55, 100, 150)
  ],
  aes(
    x = paste0(">", threshold),
    y = analysis_group,
    fill = signed_logp
  )
) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_grid(outcome_label ~ exposure_class, scales = "free_y", space = "free_y") +
  scale_fill_gradient2(
    low = "blue3",
    mid = "white",
    high = "red3",
    midpoint = 0,
    name = "signed -log10(p)"
  ) +
  labs(
    title = "Joint model signal strength: cumulative burden vs sharp spikes",
    subtitle = "Internal A/B metrics plus external PA-GEOS disagreement; 28-day prior window",
    x = "PM2.5 threshold",
    y = NULL
  ) +
  theme_bw(base_size = 8.5) +
  theme(
    axis.text.y = element_text(size = 4.8),
    strip.text.y = element_text(angle = 0),
    legend.position = "bottom"
  )

ggsave(
  file.path(FIG_DIR, "fig_joint_model_signal_heatmap_28day_with_pa_geos.jpeg"),
  p_heat,
  width = 12,
  height = 15,
  dpi = 500,
  limitsize = FALSE
)

# ------------------------------------------------------------
# 18) PRINT INTERPRETATION SUMMARY
# ------------------------------------------------------------

cat("\n=== JOINT MODEL DIRECTION SUMMARY WITH PA-GEOS ===\n")
print(joint_direction[order(analysis_type, degradation_domain, outcome_label, past_window_days, exposure_class)])

cat("\n=== Supported climate panels retained ===\n")
print(climate_support)

cat("\nUpdated figures saved to:\n")
cat(FIG_DIR, "\n")

# ------------------------------------------------------------
# 19) CLEAN GLOBAL SUMMARY FIGURES
#     Outdoor PurpleAir sensors; global and regional models only
#     Includes internal A/B metrics + PA-GEOS Metric 4
# ------------------------------------------------------------

std_dt <- copy(model_dt)

std_dt[, exposure_class := fifelse(
  term == "past_excess_pm_sum", "Cumulative burden",
  fifelse(
    term == "past_max_excess_pm", "Sharp spike",
    fifelse(term == "past_high_prop", "Frequency above threshold", term)
  )
)]

std_dt[, exposure_class := factor(
  exposure_class,
  levels = c("Frequency above threshold", "Cumulative burden", "Sharp spike")
)]

std_dt[, degradation_domain := fifelse(
  outcome_key %in% c("future_pa_geos_abs", "future_pa_geos_signed"),
  "External-reference disagreement",
  "Internal sensor disagreement"
)]

std_dt <- std_dt[
  analysis_type == "global_region" &
    n_rows >= MIN_PLOT_ROWS &
    n_stations >= MIN_PLOT_STATIONS &
    outcome_key %in% outcomes_all
]

std_dt[, ci_half_width := (upr - lwr) / 2]

std_dt[, standardized_signal_abs := fifelse(
  is.finite(ci_half_width) & ci_half_width > 0,
  abs(estimate - fifelse(effect_type == "odds_ratio", 1, 0)) / ci_half_width,
  NA_real_
)]

std_dt <- std_dt[is.finite(standardized_signal_abs)]

std_summary <- std_dt[
  ,
  .(
    n_models = .N,
    median_signal = median(standardized_signal_abs, na.rm = TRUE),
    mean_signal = mean(standardized_signal_abs, na.rm = TRUE),
    pct_significant = mean(p_value < 0.05, na.rm = TRUE),
    median_n_rows = median(n_rows, na.rm = TRUE),
    median_n_stations = median(n_stations, na.rm = TRUE)
  ),
  by = .(
    degradation_domain,
    outcome_label,
    past_window_days,
    exposure_class
  )
]

fwrite(
  std_summary,
  file.path(OUT_DIR, "standardized_exposure_signal_summary_global_region_only_with_pa_geos.csv")
)

p_std <- ggplot(
  std_summary[
    past_window_days %in% c(7, 14, 28)
  ],
  aes(
    x = factor(past_window_days),
    y = median_signal,
    fill = exposure_class
  )
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.68
  ) +
  facet_grid(
    outcome_label ~ .,
    scales = "free_y"
  ) +
  theme_bw(base_size = 11) +
  labs(
    title = "Comparable high-PM exposure signals for future degradation",
    subtitle = paste0(
      "Outdoor PurpleAir sensors; global and regional models only. ",
      "Includes internal A/B metrics and external PA-GEOS disagreement."
    ),
    x = "Prior exposure window, days",
    y = "Median standardized signal strength",
    fill = NULL
  ) +
  theme(
    legend.position = "bottom",
    strip.text.y = element_text(angle = 0),
    axis.text.x = element_text(size = 10)
  )

ggsave(
  file.path(FIG_DIR, "fig_standardized_exposure_signal_summary_global_only_with_pa_geos.jpeg"),
  p_std,
  width = 10,
  height = 11,
  dpi = 500
)

p_sig <- ggplot(
  std_summary[
    past_window_days %in% c(7, 14, 28)
  ],
  aes(
    x = factor(past_window_days),
    y = pct_significant,
    fill = exposure_class
  )
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.68
  ) +
  facet_grid(
    outcome_label ~ .,
    scales = "free_y"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  theme_bw(base_size = 11) +
  labs(
    title = "Share of supported models with statistically detectable high-PM effects",
    subtitle = "Outdoor PurpleAir sensors; global and regional models only; p < 0.05",
    x = "Prior exposure window, days",
    y = "Percent of models significant",
    fill = NULL
  ) +
  theme(
    legend.position = "bottom",
    strip.text.y = element_text(angle = 0),
    axis.text.x = element_text(size = 10)
  )

ggsave(
  file.path(FIG_DIR, "fig_percent_significant_exposure_summary_global_only_with_pa_geos.jpeg"),
  p_sig,
  width = 10,
  height = 11,
  dpi = 500
)

compare_summary <- std_summary[
  ,
  .(
    median_signal = median(median_signal, na.rm = TRUE),
    pct_significant = mean(pct_significant, na.rm = TRUE)
  ),
  by = .(degradation_domain, exposure_class)
]

compare_summary[, exposure_class := factor(
  exposure_class,
  levels = c(
    "Frequency above threshold",
    "Cumulative burden",
    "Sharp spike"
  )
)]

p_compare <- ggplot(
  compare_summary,
  aes(
    exposure_class,
    median_signal,
    fill = exposure_class
  )
) +
  geom_col(width = 0.7) +
  facet_wrap(~ degradation_domain, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(
    title = "Overall degradation signal strength by exposure type",
    subtitle = "Outdoor PurpleAir sensors; global and regional models only",
    x = NULL,
    y = "Median standardized signal strength",
    fill = NULL
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 12, hjust = 1)
  )

ggsave(
  file.path(FIG_DIR, "fig_overall_exposure_type_comparison_with_pa_geos.jpeg"),
  p_compare,
  width = 10,
  height = 5,
  dpi = 500
)

cat("\n=== OVERALL EXPOSURE-TYPE COMPARISON WITH PA-GEOS ===\n")
print(compare_summary)

cat("\nUpdated summary figures saved to:\n")
cat(file.path(FIG_DIR, "fig_standardized_exposure_signal_summary_global_only_with_pa_geos.jpeg"), "\n")
cat(file.path(FIG_DIR, "fig_percent_significant_exposure_summary_global_only_with_pa_geos.jpeg"), "\n")
cat(file.path(FIG_DIR, "fig_overall_exposure_type_comparison_with_pa_geos.jpeg"), "\n")

# ------------------------------------------------------------
# 20) WORD-READY TABLES FOR HIGH-PM FUTURE DEGRADATION ANALYSIS
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(officer)
  library(flextable)
})

TABLE_DIR <- file.path(OUT_DIR, "word_ready_tables")
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)

model_file <- file.path(OUT_DIR, "high_pm_future_degradation_model_results_clean.csv")

if (!exists("model_dt")) {
  model_dt <- fread(model_file)
}

# -----------------------------
# Formatting helpers
# -----------------------------

fmt_num <- function(x, digits = 3) {
  ifelse(
    is.na(x),
    "",
    formatC(as.numeric(x), format = "f", digits = digits, big.mark = ",")
  )
}

fmt_int <- function(x) {
  ifelse(
    is.na(x),
    "",
    formatC(as.integer(round(x)), format = "d", big.mark = ",")
  )
}

fmt_p <- function(p) {
  ifelse(
    is.na(p), "",
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
}

fmt_est_ci <- function(est, lwr, upr, digits = 3) {
  ifelse(
    is.na(est) | is.na(lwr) | is.na(upr),
    "",
    paste0(
      fmt_num(est, digits),
      " (",
      fmt_num(lwr, digits),
      ", ",
      fmt_num(upr, digits),
      ")"
    )
  )
}

clean_label <- function(x) {
  x <- gsub("\n", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

make_docx_table <- function(tab, title, note, outfile) {

  if (is.null(tab) || nrow(tab) == 0 || ncol(tab) == 0) {
    warning("Skipping empty table: ", title)
    return(invisible(NULL))
  }

  doc <- read_docx()
  doc <- body_add_par(doc, title, style = "heading 1")

  ft <- flextable(tab)
  ft <- theme_booktabs(ft)
  ft <- autofit(ft)
  ft <- fontsize(ft, size = 8, part = "all")
  ft <- bold(ft, part = "header")
  ft <- align(ft, align = "center", part = "all")
  ft <- align(ft, j = 1, align = "left", part = "body")

  doc <- body_add_flextable(doc, ft)

  if (!is.null(note) && nzchar(note)) {
    doc <- body_add_par(doc, note, style = "Normal")
  }

  print(doc, target = outfile)
  invisible(outfile)
}

# -----------------------------
# Standardize labels
# -----------------------------

model_dt[, threshold_label := paste0(">", threshold, " µg/m³")]
model_dt[, window_label := paste0(past_window_days, "-day")]
model_dt[, group_clean := clean_label(analysis_group)]

model_dt[, outcome_label := fifelse(
  outcome_key == "future_pa_geos_abs",
  "Future |PA-GEOS|",
  fifelse(
    outcome_key == "future_pa_geos_signed",
    "Future PA-GEOS",
    outcome_label
  )
)]

model_dt[, degradation_domain := fifelse(
  outcome_key %in% c("future_pa_geos_abs", "future_pa_geos_signed"),
  "External-reference disagreement",
  "Internal sensor disagreement"
)]

model_dt[, exposure_class := fifelse(
  term == "past_excess_pm_sum", "Cumulative burden",
  fifelse(
    term == "past_max_excess_pm", "Sharp spike",
    fifelse(term == "past_high_prop", "Frequency above threshold", term)
  )
)]

model_dt[, model_label := fifelse(
  model_name == "frequency_only", "Frequency only",
  fifelse(
    model_name == "cumulative_only", "Cumulative only",
    fifelse(
      model_name == "spike_only", "Spike only",
      "Joint cumulative + spike"
    )
  )
)]

model_dt[, effect_label := fifelse(
  effect_type == "odds_ratio",
  "Odds ratio",
  "Mean difference"
)]

model_dt[, estimate_95ci := fmt_est_ci(estimate, lwr, upr, 3)]
model_dt[, p_value_fmt := fmt_p(p_value)]

supported_dt <- model_dt[
  n_rows >= MIN_PLOT_ROWS &
    n_stations >= MIN_PLOT_STATIONS
]

fwrite(
  supported_dt,
  file.path(TABLE_DIR, "high_pm_supported_model_results_for_tables.csv")
)

# ------------------------------------------------------------
# Table S19. Global high-PM exposure associations
# ------------------------------------------------------------

tab_s19 <- supported_dt[
  analysis_type == "global_region" &
    analysis_group == "Global" &
    outcome_key %in% c(
      "future_flagged",
      "future_ab_abs",
      "future_ab_rel",
      "future_pa_geos_abs",
      "future_pa_geos_signed"
    ),
  .(
    Outcome = outcome_label,
    `Exposure model` = model_label,
    `Exposure metric` = exposure_class,
    Threshold = threshold_label,
    `Prior window` = window_label,
    `Effect scale` = scale_unit,
    `Effect type` = effect_label,
    `Estimate (95% CI)` = estimate_95ci,
    `P value` = p_value_fmt,
    `N rows` = fmt_int(n_rows),
    `N sensors` = fmt_int(n_stations)
  )
][order(Outcome, `Prior window`, Threshold, `Exposure model`, `Exposure metric`)]

fwrite(tab_s19, file.path(TABLE_DIR, "Table_S19_global_high_pm_associations.csv"))

make_docx_table(
  tab_s19,
  "Table S19. Global associations between prior high-PM exposure and future degradation",
  "Models include sensor fixed effects, spline terms for contemporaneous mean PM₂.₅, PM₂.₅ variability, and operational age, and an offset for previous flagged degradation. Standard errors are clustered at the sensor level. Odds ratios are reported for future flagged degradation; mean differences are reported for continuous outcomes.",
  file.path(TABLE_DIR, "Table_S19_global_high_pm_associations.docx")
)

# ------------------------------------------------------------
# Table S20. Regional joint cumulative + spike models, 28-day window
# ------------------------------------------------------------

tab_s20 <- supported_dt[
  analysis_type == "global_region" &
    analysis_group != "Global" &
    model_name == "joint_cumulative_spike" &
    past_window_days == 28 &
    term %in% c("past_excess_pm_sum", "past_max_excess_pm") &
    outcome_key %in% c(
      "future_flagged",
      "future_ab_abs",
      "future_ab_rel",
      "future_pa_geos_abs",
      "future_pa_geos_signed"
    ),
  .(
    Region = group_clean,
    Outcome = outcome_label,
    `Exposure metric` = exposure_class,
    Threshold = threshold_label,
    `Effect scale` = scale_unit,
    `Effect type` = effect_label,
    `Estimate (95% CI)` = estimate_95ci,
    `P value` = p_value_fmt,
    `N rows` = fmt_int(n_rows),
    `N sensors` = fmt_int(n_stations)
  )
][order(Region, Outcome, Threshold, `Exposure metric`)]

fwrite(tab_s20, file.path(TABLE_DIR, "Table_S20_regional_joint_28day_models.csv"))

make_docx_table(
  tab_s20,
  "Table S20. Regional joint models separating cumulative PM burden and sharp PM spikes",
  "Results are from joint cumulative-burden and spike-exposure models using the 28-day prior exposure window. Models include sensor fixed effects, spline covariates, and sensor-clustered standard errors.",
  file.path(TABLE_DIR, "Table_S20_regional_joint_28day_models.docx")
)

# ------------------------------------------------------------
# Table S21. Climate-zone joint cumulative + spike models, 28-day window
# ------------------------------------------------------------

tab_s21 <- supported_dt[
  analysis_type == "climate" &
    model_name == "joint_cumulative_spike" &
    past_window_days == 28 &
    term %in% c("past_excess_pm_sum", "past_max_excess_pm") &
    outcome_key %in% c(
      "future_flagged",
      "future_ab_abs",
      "future_ab_rel",
      "future_pa_geos_abs",
      "future_pa_geos_signed"
    ),
  .(
    `Climate zone` = group_clean,
    Outcome = outcome_label,
    `Exposure metric` = exposure_class,
    Threshold = threshold_label,
    `Effect scale` = scale_unit,
    `Effect type` = effect_label,
    `Estimate (95% CI)` = estimate_95ci,
    `P value` = p_value_fmt,
    `N rows` = fmt_int(n_rows),
    `N sensors` = fmt_int(n_stations)
  )
][order(`Climate zone`, Outcome, Threshold, `Exposure metric`)]

fwrite(tab_s21, file.path(TABLE_DIR, "Table_S21_climate_joint_28day_models.csv"))

make_docx_table(
  tab_s21,
  "Table S21. Climate-zone joint models separating cumulative PM burden and sharp PM spikes",
  "Results are restricted to sensors assigned to a single Köppen–Geiger climate zone throughout the observation period. Models use the 28-day prior exposure window and include sensor fixed effects and sensor-clustered standard errors.",
  file.path(TABLE_DIR, "Table_S21_climate_joint_28day_models.docx")
)

# ------------------------------------------------------------
# Table S22. Region × climate joint cumulative + spike models, 28-day window
# ------------------------------------------------------------

tab_s22 <- supported_dt[
  analysis_type == "region_climate" &
    model_name == "joint_cumulative_spike" &
    past_window_days == 28 &
    term %in% c("past_excess_pm_sum", "past_max_excess_pm") &
    outcome_key %in% c(
      "future_flagged",
      "future_ab_abs",
      "future_ab_rel",
      "future_pa_geos_abs",
      "future_pa_geos_signed"
    ),
  .(
    `Region × climate` = group_clean,
    Outcome = outcome_label,
    `Exposure metric` = exposure_class,
    Threshold = threshold_label,
    `Effect scale` = scale_unit,
    `Effect type` = effect_label,
    `Estimate (95% CI)` = estimate_95ci,
    `P value` = p_value_fmt,
    `N rows` = fmt_int(n_rows),
    `N sensors` = fmt_int(n_stations)
  )
][order(`Region × climate`, Outcome, Threshold, `Exposure metric`)]

fwrite(tab_s22, file.path(TABLE_DIR, "Table_S22_region_x_climate_joint_28day_models.csv"))

make_docx_table(
  tab_s22,
  "Table S22. Region × climate joint models separating cumulative PM burden and sharp PM spikes",
  "Results are shown for supported region × Köppen–Geiger climate strata. Models use the 28-day prior exposure window and include sensor fixed effects and sensor-clustered standard errors.",
  file.path(TABLE_DIR, "Table_S22_region_x_climate_joint_28day_models.docx")
)

# ------------------------------------------------------------
# Table S23. Summary of joint model signal direction
# ------------------------------------------------------------

joint_summary <- supported_dt[
  model_name == "joint_cumulative_spike" &
    term %in% c("past_excess_pm_sum", "past_max_excess_pm")
]

joint_summary[, significant := p_value < 0.05]

tab_s23_raw <- joint_summary[
  ,
  .(
    n_models = .N,
    n_significant = sum(significant, na.rm = TRUE),
    pct_significant = 100 * mean(significant, na.rm = TRUE),
    median_estimate = median(estimate, na.rm = TRUE),
    median_lwr = median(lwr, na.rm = TRUE),
    median_upr = median(upr, na.rm = TRUE),
    median_n_rows = median(n_rows, na.rm = TRUE),
    median_n_sensors = median(n_stations, na.rm = TRUE)
  ),
  by = .(
    `Analysis type` = analysis_type,
    Domain = degradation_domain,
    Outcome = outcome_label,
    `Prior window` = window_label,
    `Exposure metric` = exposure_class,
    `Effect type` = effect_label
  )
]

tab_s23 <- tab_s23_raw[
  ,
  .(
    `Analysis type`,
    Domain,
    Outcome,
    `Prior window`,
    `Exposure metric`,
    `Effect type`,
    `Models` = fmt_int(n_models),
    `Significant models` = fmt_int(n_significant),
    `% significant` = fmt_num(pct_significant, 1),
    `Median estimate (median CI)` = fmt_est_ci(median_estimate, median_lwr, median_upr, 3),
    `Median rows` = fmt_int(median_n_rows),
    `Median sensors` = fmt_int(median_n_sensors)
  )
][order(`Analysis type`, Domain, Outcome, `Prior window`, `Exposure metric`)]

fwrite(tab_s23_raw, file.path(TABLE_DIR, "Table_S23_joint_model_signal_summary_raw.csv"))
fwrite(tab_s23, file.path(TABLE_DIR, "Table_S23_joint_model_signal_summary.csv"))

make_docx_table(
  tab_s23,
  "Table S23. Summary of joint cumulative-burden and spike-exposure model signals",
  "Summary table aggregates supported joint models across global, regional, climate-zone, and region × climate analyses. Values summarize the number of fitted models, the share with p < 0.05, and the median effect estimate across models.",
  file.path(TABLE_DIR, "Table_S23_joint_model_signal_summary.docx")
)

# ------------------------------------------------------------
# Table S24. Empirical exposure-response summaries
# ------------------------------------------------------------

emp_file <- file.path(OUT_DIR, "high_pm_future_degradation_empirical_summary.csv")

if (file.exists(emp_file)) {

  emp_dt <- fread(emp_file)

  emp_dt[, group_clean := clean_label(analysis_group)]
  emp_dt[, threshold_label := paste0(">", threshold, " µg/m³")]
  emp_dt[, window_label := paste0(past_window_days, "-day")]

  emp_keep <- emp_dt[
    rows >= MIN_EMP_ROWS &
      stations >= MIN_EMP_STATIONS &
      analysis_group == "Global" &
      past_window_days == 28
  ]

  tab_s24 <- emp_keep[
    ,
    .(
      Group = group_clean,
      Outcome = outcome_label,
      `Exposure metric` = exposure_label,
      Threshold = threshold_label,
      `Exposure bin` = as.character(exposure_bin),
      `Mean exposure` = fmt_num(mean_exposure, 3),
      `Future outcome mean` = fmt_num(outcome_mean, 3),
      `Rows` = fmt_int(rows),
      `Sensors` = fmt_int(stations)
    )
  ][order(Outcome, `Exposure metric`, Threshold, `Exposure bin`)]

  fwrite(tab_s24, file.path(TABLE_DIR, "Table_S24_empirical_exposure_response_global_28day.csv"))

  make_docx_table(
    tab_s24,
    "Table S24. Empirical exposure-response summaries for global 28-day high-PM exposure",
    "Empirical summaries group observations into exposure bins and report mean future degradation outcomes within each bin. Rows are restricted to bins with adequate support.",
    file.path(TABLE_DIR, "Table_S24_empirical_exposure_response_global_28day.docx")
  )
}

cat("\nWord-ready high-PM degradation tables saved to:\n")
cat(TABLE_DIR, "\n")

# ---- source chunk #33 (Rmd lines 17399-18951) ------------------------

# ============================================================
# COMBINED SENSOR LIFECYCLE + HIGH-PM FUTURE DEGRADATION ANALYSIS
#
# GLOBAL + REGION:
#   uses outdoor_final.RData
#
# CLIMATE + REGION x CLIMATE:
#   uses outdoor_with_koppen_timevarying.RData
#   restricted to valid single-zone Köppen stations
#
# Produces:
#   1) Station lifetime / KM / RMST / hazard by region
#   2) Effective lifetime until cumulative flagged degradation threshold
#   3) Pre-failure instability analysis
#   4) Drift coupling: PM drift vs temp/RH drift
#   5) High-PM future degradation models:
#        cumulative burden vs sharp spike vs frequency above threshold
#   6) Global, regional, climate, and region x climate model outputs
#   7) Empirical SI exposure-response plots
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(survival)
  library(cowplot)
  library(fixest)
  library(splines)
})

# ------------------------------------------------------------
# 0) PATHS / SETTINGS
# ------------------------------------------------------------

global_rdata  <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_final.RData"
climate_rdata <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_with_koppen_timevarying.RData"
transition_csv <- "/Users/priyanka/Downloads/PurpleAir_EPA/station_koppen_transition_summary.csv"

OUTDIR <- "/Users/priyanka/Downloads/PurpleAir_EPA/combined_lifecycle_highpm_degradation"
FIGDIR <- file.path(OUTDIR, "figures")
MODELDIR <- file.path(OUTDIR, "models")

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(MODELDIR, showWarnings = FALSE, recursive = TRUE)

BIN_W <- 24L
DROP_FIRST_DAYS <- 7L
DROP_H <- DROP_FIRST_DAYS * 24L

PM_THRESHOLDS <- c(35, 55, 100, 150)
PAST_WINDOWS <- c(7L, 14L, 28L)
FUTURE_WINDOW <- 7L
MIN_FUTURE_HOURS <- 24L * 3L

MIN_ROWS_MODEL <- 500L
MIN_STATIONS_MODEL <- 20L
SAVE_SAMPLE_PER_COMBO <- 25000L

FLAG_ABS_DIFF_CUTOFF <- 5
FLAG_PERCENTILE <- 0.85

CENSOR_BUFFER_DAYS <- 7L
HAZARD_BIN_HOURS <- 24L * 30L
MIN_AT_RISK_HAZARD <- 25L

EFFECTIVE_THRESHOLDS <- c(0.05, 0.10)
MIN_POST_HOURS_EFFECTIVE <- 100L

PM_MAX <- 1500

a_col_default <- "pm25_cf_1_a_hourly(ug/m3)_pm25_cf_1_a"
b_col_default <- "pm25_cf_1_b_hourly(ug/m3)_pm25_cf_1_b"
geos_col <- "geoscf_pm25"
temp_col <- "temperature_hourly(C)_temperature"
rh_col <- "humidity_hourly(%)_humidity"
era5_temp_col <- "era5_t2m_c"
era5_rh_col <- "era5_rh_pct"
koppen_col <- "koppen_code"

out_model_csv <- file.path(OUTDIR, "high_pm_future_degradation_model_results.csv")
out_emp_csv <- file.path(OUTDIR, "high_pm_future_degradation_empirical_summary.csv")
out_sample_csv <- file.path(OUTDIR, "high_pm_future_degradation_analysis_SAMPLE.csv")

for (f in c(out_model_csv, out_emp_csv, out_sample_csv)) {
  if (file.exists(f)) file.remove(f)
}

say <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste0(..., collapse = ""), "\n")
}

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  tolower(x)
}

fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

safe_quantile <- function(x, probs) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(quantile(x, probs = probs, na.rm = TRUE, type = 7))
}

calc_slope <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 10L) return(NA_real_)
  vx <- var(x)
  if (!is.finite(vx) || vx <= 0) return(NA_real_)
  cov(x, y) / vx
}

safe_empirical_cut <- function(x, probs = seq(0, 1, 0.1)) {
  qs <- unique(quantile(x, probs = probs, na.rm = TRUE))
  if (length(qs) < 3L) return(NULL)
  qs
}

theme_pub <- theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

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
    stop("Could not find outdoor or outdoor_final in: ", infile)
  }

  setDT(x)
  x[]
}

get_single_zone_table <- function() {
  if (!file.exists(transition_csv)) {
    stop("Missing Köppen transition file: ", transition_csv)
  }

  z <- fread(transition_csv)
  z[, station := as.character(station)]

  z <- z[
    !is.na(first_koppen) &
      first_koppen != "" &
      n_koppen_codes == 1,
    .(
      station,
      climate_zone = as.character(first_koppen),
      climate_desc = as.character(first_koppen_desc)
    )
  ]

  z[is.na(climate_desc) | climate_desc == "", climate_desc := climate_zone]
  z[, climate_label := paste0(climate_zone, "\n", climate_desc)]
  unique(z)
}

single_zone <- get_single_zone_table()

# ------------------------------------------------------------
# 2) BUILD COMPACT DAILY DATASET
# ------------------------------------------------------------

make_daily_dataset <- function(infile, require_climate = FALSE) {

  say("Loading: ", infile)

  outdoor <- load_outdoor_object(infile)

  a_col <- a_col_default
  b_col <- b_col_default

  if (!(a_col %in% names(outdoor)) && "cf1_a" %in% names(outdoor)) a_col <- "cf1_a"
  if (!(b_col %in% names(outdoor)) && "cf1_b" %in% names(outdoor)) b_col <- "cf1_b"

  needed <- c("station", "region", "datetime_utc", a_col, b_col)

  has_geos <- geos_col %in% names(outdoor)
  has_temp <- temp_col %in% names(outdoor)
  has_rh <- rh_col %in% names(outdoor)
  has_era5_temp <- era5_temp_col %in% names(outdoor)
  has_era5_rh <- era5_rh_col %in% names(outdoor)

  if (has_geos) needed <- c(needed, geos_col)
  if (has_temp) needed <- c(needed, temp_col)
  if (has_rh) needed <- c(needed, rh_col)
  if (has_era5_temp) needed <- c(needed, era5_temp_col)
  if (has_era5_rh) needed <- c(needed, era5_rh_col)
  if (require_climate && koppen_col %in% names(outdoor)) needed <- c(needed, koppen_col)

  miss <- setdiff(needed, names(outdoor))
  if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))

  outdoor[, station := as.character(station)]

  if (require_climate) {
    outdoor <- outdoor[station %in% single_zone$station]
    say("Rows after single-zone restriction: ", format(nrow(outdoor), big.mark = ","))
    say("Stations after single-zone restriction: ", format(uniqueN(outdoor$station), big.mark = ","))
  }

  dt <- outdoor[, ..needed]
  rm(outdoor)
  gc()

  setnames(dt, a_col, "pm_a")
  setnames(dt, b_col, "pm_b")

  dt[, `:=`(
    station = as.character(station),
    region = as.character(region),
    ts = as.POSIXct(datetime_utc, tz = "UTC"),
    pm_a = as.numeric(pm_a),
    pm_b = as.numeric(pm_b)
  )]

  dt[, datetime_utc := NULL]

  dt[, pm_cf1 := (pm_a + pm_b) / 2]
  dt <- dt[
    !is.na(station) &
      !is.na(region) &
      !is.na(ts) &
      is.finite(pm_a) &
      is.finite(pm_b) &
      is.finite(pm_cf1) &
      pm_cf1 >= 0 &
      pm_cf1 <= PM_MAX
  ]

  dt[, abs_ab := abs(pm_a - pm_b)]
  dt[, ab_bias := pm_a - pm_b]
  dt[, ab_rel := fifelse(
    (pm_a + pm_b) > 0,
    2 * abs(pm_a - pm_b) / (pm_a + pm_b),
    NA_real_
  )]

  station_cut <- dt[
    is.finite(ab_rel),
    .(thr85_ab_rel = quantile(ab_rel, FLAG_PERCENTILE, na.rm = TRUE, type = 7)),
    by = station
  ]

  dt <- station_cut[dt, on = "station"]

  dt[, flagged := as.integer(
    abs_ab > FLAG_ABS_DIFF_CUTOFF &
      ab_rel > thr85_ab_rel
  )]

  dt[is.na(flagged), flagged := 0L]

  if (has_geos) {
    dt[, geos_pm := as.numeric(get(geos_col))]
  } else {
    dt[, geos_pm := NA_real_]
  }

  if (has_temp) {
    dt[, pa_temp := as.numeric(get(temp_col))]
  } else {
    dt[, pa_temp := NA_real_]
  }

  if (has_rh) {
    dt[, pa_rh := as.numeric(get(rh_col))]
  } else {
    dt[, pa_rh := NA_real_]
  }

  if (has_era5_temp) {
    dt[, era5_temp := as.numeric(get(era5_temp_col))]
  } else {
    dt[, era5_temp := NA_real_]
  }

  if (has_era5_rh) {
    dt[, era5_rh := as.numeric(get(era5_rh_col))]
  } else {
    dt[, era5_rh := NA_real_]
  }

  dt[, `:=`(
    pm_geos_signed = pm_cf1 - geos_pm,
    pm_geos_abs = abs(pm_cf1 - geos_pm),
    temp_div = pa_temp - era5_temp,
    rh_div = pa_rh - era5_rh
  )]

  setorder(dt, station, ts)
  dt[, first_ts := ts[1L], by = station]
  dt[, op_hours := as.numeric(difftime(ts, first_ts, units = "hours"))]
  dt[, op_bin := floor(op_hours / BIN_W) * BIN_W]
  dt[, date := as.IDate(ts)]

  dt <- dt[is.finite(op_bin) & op_bin >= DROP_H]

  daily <- dt[, .(
    region = region[1L],
    n = .N,
    n_flagged = sum(flagged, na.rm = TRUE),
    pm_mean = mean(pm_cf1, na.rm = TRUE),
    pm_max = max(pm_cf1, na.rm = TRUE),
    pm_sd = sd(pm_cf1, na.rm = TRUE),

    ab_abs = mean(abs_ab, na.rm = TRUE),
    ab_bias = mean(ab_bias, na.rm = TRUE),
    ab_rel = mean(ab_rel, na.rm = TRUE),
    flagged_rate = mean(flagged, na.rm = TRUE),

    pa_geos_abs = mean(pm_geos_abs, na.rm = TRUE),
    pa_geos_signed = mean(pm_geos_signed, na.rm = TRUE),
    temp_div = mean(temp_div, na.rm = TRUE),
    rh_div = mean(rh_div, na.rm = TRUE)
  ), by = .(station, op_bin)]

  station_life <- dt[, .(
    region = region[1L],
    first_ts = min(ts, na.rm = TRUE),
    last_ts = max(ts, na.rm = TRUE),
    max_op_hours = max(op_hours, na.rm = TRUE),
    max_op_years = max(op_hours, na.rm = TRUE) / (24 * 365.25),
    n_obs = .N
  ), by = station]

  rm(dt, station_cut)
  gc()

  daily[, n_unflagged := n - n_flagged]

  if (require_climate) {
    daily <- merge(daily, single_zone, by = "station", all.x = FALSE)
    station_life <- merge(station_life, single_zone, by = "station", all.x = FALSE)
  } else {
    daily[, `:=`(
      climate_zone = NA_character_,
      climate_desc = NA_character_,
      climate_label = NA_character_
    )]
    station_life[, `:=`(
      climate_zone = NA_character_,
      climate_desc = NA_character_,
      climate_label = NA_character_
    )]
  }

  setorder(daily, station, op_bin)

  say("Daily rows: ", format(nrow(daily), big.mark = ","))
  say("Stations: ", format(uniqueN(daily$station), big.mark = ","))

  list(daily = daily, station_life = station_life)
}

# ------------------------------------------------------------
# 3) FUTURE OUTCOMES + EXPOSURES
# ------------------------------------------------------------

add_future_outcomes <- function(daily) {
  setorder(daily, station, op_bin)

  daily[, future_flagged_h := frollsum(
    shift(n_flagged, n = 1L, type = "lead"),
    n = FUTURE_WINDOW,
    align = "left",
    na.rm = TRUE
  ), by = station]

  daily[, future_total_h := frollsum(
    shift(n, n = 1L, type = "lead"),
    n = FUTURE_WINDOW,
    align = "left",
    na.rm = TRUE
  ), by = station]

  daily[, future_flagged_rate := fifelse(
    future_total_h > 0,
    future_flagged_h / future_total_h,
    NA_real_
  )]

  for (v in c("ab_abs", "ab_bias", "ab_rel", "pa_geos_abs", "pa_geos_signed", "temp_div", "rh_div")) {
    daily[, paste0("future_", v) := frollmean(
      shift(get(v), n = 1L, type = "lead"),
      n = FUTURE_WINDOW,
      align = "left",
      na.rm = TRUE
    ), by = station]
  }

  daily[]
}

make_exposure_dt <- function(base_daily, thr, pw, analysis_group, group_var) {
  dd <- copy(base_daily)

  dd[, high_pm_day := as.integer(pm_mean > thr)]
  dd[, excess_pm_day := pmax(pm_mean - thr, 0)]
  dd[, max_excess_pm_day := pmax(pm_max - thr, 0)]

  dd[, past_high_n := frollsum(
    shift(high_pm_day, n = 1L, type = "lag", fill = 0L),
    n = pw,
    align = "right",
    na.rm = TRUE
  ), by = station]

  dd[, past_obs_n := frollsum(
    shift(as.integer(is.finite(pm_mean)), n = 1L, type = "lag", fill = 0L),
    n = pw,
    align = "right",
    na.rm = TRUE
  ), by = station]

  dd[, past_high_prop := fifelse(past_obs_n > 0, past_high_n / past_obs_n, NA_real_)]

  dd[, past_excess_pm_sum := frollsum(
    shift(excess_pm_day, n = 1L, type = "lag", fill = 0),
    n = pw,
    align = "right",
    na.rm = TRUE
  ), by = station]

  dd[, past_excess_pm_mean := fifelse(past_obs_n > 0, past_excess_pm_sum / past_obs_n, NA_real_)]

  dd[, past_max_pm := frollapply(
    shift(pm_max, n = 1L, type = "lag"),
    n = pw,
    FUN = max,
    align = "right",
    na.rm = TRUE
  ), by = station]

  dd[, past_max_excess_pm := pmax(past_max_pm - thr, 0)]

  dd[, n_prev := shift(n, 1L), by = station]
  dd[, f_prev := shift(n_flagged, 1L), by = station]
  dd[, p_prev_smooth := (f_prev + 0.5) / (n_prev + 1.0)]
  dd[, offset_logit_prev := qlogis(p_prev_smooth)]

  dd[, step_h := op_bin - shift(op_bin), by = station]
  dd[, is_consecutive := is.finite(step_h) & step_h == BIN_W]

  out <- dd[
    is_consecutive &
      is.finite(offset_logit_prev) &
      past_obs_n >= max(3L, floor(pw / 2)) &
      is.finite(pm_mean) &
      is.finite(pm_sd) &
      is.finite(past_high_prop) &
      is.finite(past_excess_pm_sum) &
      is.finite(past_max_excess_pm) &
      is.finite(future_flagged_h) &
      is.finite(future_total_h) &
      future_total_h >= MIN_FUTURE_HOURS,
    .(
      station,
      region,
      climate_zone,
      climate_desc,
      climate_label,
      op_bin,
      future_flagged_h,
      future_total_h,
      future_unflagged_h = future_total_h - future_flagged_h,
      future_flagged_rate,
      future_ab_abs,
      future_ab_bias,
      future_ab_rel,
      future_pa_geos_abs,
      future_pa_geos_signed,
      future_temp_div,
      future_rh_div,
      past_high_prop,
      past_excess_pm_sum,
      past_excess_pm_mean,
      past_max_pm,
      past_max_excess_pm,
      pm_mean,
      pm_sd,
      offset_logit_prev,
      threshold = thr,
      past_window_days = pw,
      analysis_group = analysis_group,
      group_var = group_var
    )
  ][future_flagged_h >= 0 & future_total_h > future_flagged_h]

  rm(dd)
  gc()

  out[]
}

# ------------------------------------------------------------
# 4) MODEL HELPERS
# ------------------------------------------------------------

extract_fixest_term <- function(fit, term, scale_mult = 1, effect_type = "mean_difference") {
  s <- summary(fit, vcov = ~ station)
  cf <- s$coeftable

  if (!(term %in% rownames(cf))) return(NULL)

  beta <- cf[term, "Estimate"]
  se <- cf[term, "Std. Error"]
  pcol <- grep("^Pr", colnames(cf), value = TRUE)[1]
  pval <- cf[term, pcol]

  if (effect_type == "odds_ratio") {
    data.table(
      term = term,
      effect_type = effect_type,
      beta = beta,
      se = se,
      estimate = exp(beta * scale_mult),
      lwr = exp((beta - 1.96 * se) * scale_mult),
      upr = exp((beta + 1.96 * se) * scale_mult),
      p_value = pval
    )
  } else {
    data.table(
      term = term,
      effect_type = effect_type,
      beta = beta,
      se = se,
      estimate = beta * scale_mult,
      lwr = (beta - 1.96 * se) * scale_mult,
      upr = (beta + 1.96 * se) * scale_mult,
      p_value = pval
    )
  }
}

fit_one_model_set <- function(dt_fit, outcome_key, outcome_var, outcome_label) {
  out_rows <- list()

  if (outcome_key == "future_flagged") {
    dt_fit <- dt_fit[
      is.finite(future_flagged_h) &
        is.finite(future_total_h) &
        future_total_h > future_flagged_h
    ]
  } else {
    dt_fit <- dt_fit[is.finite(get(outcome_var))]
  }

  if (nrow(dt_fit) < MIN_ROWS_MODEL || uniqueN(dt_fit$station) < MIN_STATIONS_MODEL) {
    return(data.table())
  }

  common_rhs <- paste0(
    "ns(pm_mean, df = 4) + ",
    "ns(pm_sd, df = 3) + ",
    "ns(op_bin, df = 4) + ",
    "offset(offset_logit_prev)"
  )

  model_specs <- list(
    frequency_only = list(
      rhs = paste0("past_high_prop + ", common_rhs),
      terms = c("past_high_prop"),
      scale = c(past_high_prop = 0.10),
      scale_label = c(past_high_prop = "10 percentage-point increase")
    ),
    cumulative_only = list(
      rhs = paste0("past_excess_pm_sum + ", common_rhs),
      terms = c("past_excess_pm_sum"),
      scale = c(past_excess_pm_sum = 100),
      scale_label = c(past_excess_pm_sum = "+100 µg/m³-days")
    ),
    spike_only = list(
      rhs = paste0("past_max_excess_pm + ", common_rhs),
      terms = c("past_max_excess_pm"),
      scale = c(past_max_excess_pm = 10),
      scale_label = c(past_max_excess_pm = "+10 µg/m³ spike excess")
    ),
    joint_cumulative_spike = list(
      rhs = paste0("past_excess_pm_sum + past_max_excess_pm + ", common_rhs),
      terms = c("past_excess_pm_sum", "past_max_excess_pm"),
      scale = c(past_excess_pm_sum = 100, past_max_excess_pm = 10),
      scale_label = c(
        past_excess_pm_sum = "+100 µg/m³-days, adjusted for spike",
        past_max_excess_pm = "+10 µg/m³ spike excess, adjusted for cumulative burden"
      )
    )
  )

  for (mn in names(model_specs)) {
    spec <- model_specs[[mn]]

    if (outcome_key == "future_flagged") {
      form <- as.formula(paste0(
        "future_flagged_h / future_total_h ~ ",
        spec$rhs,
        " | station"
      ))

      fit <- tryCatch(
        feglm(
          form,
          data = dt_fit,
          family = "binomial",
          weights = ~ future_total_h
        ),
        error = function(e) NULL
      )

      effect_type <- "odds_ratio"
    } else {
      form <- as.formula(paste0(
        outcome_var,
        " ~ ",
        spec$rhs,
        " | station"
      ))

      fit <- tryCatch(
        feols(
          form,
          data = dt_fit,
          weights = ~ future_total_h
        ),
        error = function(e) NULL
      )

      effect_type <- "mean_difference"
    }

    if (is.null(fit)) next

    for (tt in spec$terms) {
      res <- extract_fixest_term(
        fit,
        term = tt,
        scale_mult = spec$scale[[tt]],
        effect_type = effect_type
      )

      if (!is.null(res)) {
        out_rows[[length(out_rows) + 1L]] <- cbind(
          data.table(
            model_name = mn,
            outcome_key = outcome_key,
            outcome_label = outcome_label,
            scale_unit = spec$scale_label[[tt]],
            n_rows = nrow(dt_fit),
            n_stations = uniqueN(dt_fit$station)
          ),
          res
        )
      }
    }

    saveRDS(
      fit,
      file.path(
        MODELDIR,
        paste0(
          "model_",
          safe_name(dt_fit$analysis_group[1]),
          "_thr", dt_fit$threshold[1],
          "_pw", dt_fit$past_window_days[1],
          "_", outcome_key,
          "_", mn,
          ".rds"
        )
      )
    )

    rm(fit)
    gc()
  }

  rbindlist(out_rows, use.names = TRUE, fill = TRUE)
}

save_empirical_summary <- function(dt_fit, exposure_var, exposure_label, outcome_var, outcome_label) {
  if (!(exposure_var %in% names(dt_fit)) || !(outcome_var %in% names(dt_fit))) return(NULL)

  dd <- dt_fit[is.finite(get(exposure_var)) & is.finite(get(outcome_var))]
  if (nrow(dd) < 100L) return(NULL)

  qs <- safe_empirical_cut(dd[[exposure_var]])
  if (is.null(qs)) return(NULL)

  dd[, exposure_bin := cut(
    get(exposure_var),
    breaks = qs,
    include.lowest = TRUE,
    ordered_result = FALSE
  )]

  out <- dd[, .(
    rows = .N,
    stations = uniqueN(station),
    mean_exposure = mean(get(exposure_var), na.rm = TRUE),
    outcome_mean = mean(get(outcome_var), na.rm = TRUE)
  ), by = exposure_bin]

  out[, `:=`(
    exposure_var = exposure_var,
    exposure_label = exposure_label,
    outcome_var = outcome_var,
    outcome_label = outcome_label,
    threshold = dt_fit$threshold[1],
    past_window_days = dt_fit$past_window_days[1],
    analysis_group = dt_fit$analysis_group[1],
    group_var = dt_fit$group_var[1]
  )]

  out[]
}

# ------------------------------------------------------------
# 5) RUN HIGH-PM ANALYSIS
# ------------------------------------------------------------

run_highpm_analysis <- function(daily, analysis_type) {
  say("Running high-PM analysis: ", analysis_type)

  if (analysis_type == "global_region") {
    groups <- c("Global", sort(unique(na.omit(daily$region))))
    get_group_dt <- function(g) if (g == "Global") daily else daily[region == g]
    group_var_name <- "region"
  } else if (analysis_type == "climate") {
    groups <- sort(unique(na.omit(daily$climate_label)))
    get_group_dt <- function(g) daily[climate_label == g]
    group_var_name <- "climate_label"
  } else if (analysis_type == "region_climate") {
    daily[, region_climate := paste(region, climate_label, sep = " | ")]
    groups <- sort(unique(na.omit(daily$region_climate)))
    get_group_dt <- function(g) daily[region_climate == g]
    group_var_name <- "region_climate"
  } else {
    stop("Unknown analysis_type")
  }

  outcomes <- list(
    list(key = "future_flagged", var = "future_flagged_rate", label = "Future flagged probability"),
    list(key = "future_ab_abs", var = "future_ab_abs", label = "Future |A-B|"),
    list(key = "future_ab_bias", var = "future_ab_bias", label = "Future A-B"),
    list(key = "future_ab_rel", var = "future_ab_rel", label = "Future relative A-B"),
    list(key = "future_pa_geos_abs", var = "future_pa_geos_abs", label = "Future |PA-GEOS|"),
    list(key = "future_pa_geos_signed", var = "future_pa_geos_signed", label = "Future PA-GEOS"),
    list(key = "future_temp_div", var = "future_temp_div", label = "Future PA-ERA5 temperature"),
    list(key = "future_rh_div", var = "future_rh_div", label = "Future PA-ERA5 RH")
  )

  for (grp in groups) {
    say("Group: ", grp)

    base_grp <- get_group_dt(grp)

    if (nrow(base_grp) < MIN_ROWS_MODEL || uniqueN(base_grp$station) < MIN_STATIONS_MODEL) {
      say("  Skipping: too few rows/stations.")
      next
    }

    setorder(base_grp, station, op_bin)

    for (thr in PM_THRESHOLDS) {
      for (pw in PAST_WINDOWS) {
        say("  Threshold=", thr, " window=", pw)

        dt_fit <- make_exposure_dt(
          base_daily = base_grp,
          thr = thr,
          pw = pw,
          analysis_group = grp,
          group_var = group_var_name
        )

        if (nrow(dt_fit) < MIN_ROWS_MODEL || uniqueN(dt_fit$station) < MIN_STATIONS_MODEL) {
          say("    Skipping threshold/window.")
          rm(dt_fit)
          gc()
          next
        }

        samp_n <- min(nrow(dt_fit), SAVE_SAMPLE_PER_COMBO)

        fwrite(
          dt_fit[sample.int(nrow(dt_fit), samp_n)],
          out_sample_csv,
          append = file.exists(out_sample_csv)
        )

        model_rows <- list()

        for (oo in outcomes) {
          if (!(oo$var %in% names(dt_fit))) next

          res <- fit_one_model_set(
            dt_fit = dt_fit,
            outcome_key = oo$key,
            outcome_var = oo$var,
            outcome_label = oo$label
          )

          if (nrow(res) > 0) {
            res[, `:=`(
              analysis_type = analysis_type,
              analysis_group = grp,
              group_var = group_var_name,
              threshold = thr,
              past_window_days = pw
            )]

            model_rows[[length(model_rows) + 1L]] <- res
          }
        }

        model_out <- rbindlist(model_rows, use.names = TRUE, fill = TRUE)

        if (nrow(model_out) > 0) {
          fwrite(model_out, out_model_csv, append = file.exists(out_model_csv))
        }

        emp_rows <- list()

        emp_outcomes <- c(
          future_flagged_rate = "Future flagged probability",
          future_ab_abs = "Future |A-B|",
          future_ab_rel = "Future relative A-B"
        )

        for (ev in c("past_high_prop", "past_excess_pm_sum", "past_max_excess_pm")) {
          exp_lab <- switch(
            ev,
            past_high_prop = "Frequency above threshold",
            past_excess_pm_sum = "Cumulative excess PM",
            past_max_excess_pm = "Sharp spike excess PM"
          )

          for (ov in names(emp_outcomes)) {
            emp <- save_empirical_summary(
              dt_fit,
              exposure_var = ev,
              exposure_label = exp_lab,
              outcome_var = ov,
              outcome_label = emp_outcomes[[ov]]
            )

            if (!is.null(emp)) emp_rows[[length(emp_rows) + 1L]] <- emp
          }
        }

        emp_out <- rbindlist(emp_rows, use.names = TRUE, fill = TRUE)

        if (nrow(emp_out) > 0) {
          fwrite(emp_out, out_emp_csv, append = file.exists(out_emp_csv))
        }

        rm(dt_fit, model_rows, model_out, emp_rows, emp_out)
        gc()
      }
    }

    rm(base_grp)
    gc()
  }
}

# ------------------------------------------------------------
# 6) LIFECYCLE ANALYSES
# ------------------------------------------------------------

run_lifecycle_analysis <- function(daily, station_life, label = "global_region") {

  say("Running lifecycle analysis: ", label)

  study_end <- max(station_life$last_ts, na.rm = TRUE)

  station_life[, event_offline := fifelse(
    last_ts < (study_end - CENSOR_BUFFER_DAYS * 24 * 3600),
    1L,
    0L
  )]

  station_life[, life_group := fifelse(
    max_op_hours < 10000, "<10k h",
    fifelse(
      max_op_hours < 30000, "10k-30k h",
      fifelse(max_op_hours < 60000, "30k-60k h", ">=60k h")
    )
  )]

  naive <- station_life[, .(
    stations = .N,
    observed_failures = sum(event_offline, na.rm = TRUE),
    censored_active = sum(event_offline == 0, na.rm = TRUE),
    mean_lifetime_years = mean(max_op_years, na.rm = TRUE),
    median_lifetime_years = median(max_op_years, na.rm = TRUE)
  ), by = region]

  fwrite(naive, file.path(OUTDIR, paste0("naive_lifetime_", label, ".csv")))

  km <- survfit(Surv(max_op_hours, event_offline) ~ region, data = station_life)
  km_tab <- summary(km)$table

  km_region <- data.table(
    region = sub("^region=", "", rownames(km_tab)),
    records = km_tab[, "records"],
    events = km_tab[, "events"],
    median_survival_hours = km_tab[, "median"],
    median_survival_lcl_hours = km_tab[, "0.95LCL"],
    median_survival_ucl_hours = km_tab[, "0.95UCL"]
  )

  km_region[, `:=`(
    median_survival_years = median_survival_hours / (24 * 365.25),
    median_survival_lcl_years = median_survival_lcl_hours / (24 * 365.25),
    median_survival_ucl_years = median_survival_ucl_hours / (24 * 365.25)
  )]

  fwrite(km_region, file.path(OUTDIR, paste0("km_lifetime_", label, ".csv")))

  tau <- safe_quantile(station_life$max_op_hours, 0.90)

  rmst_one <- function(dd) {
    fit <- survfit(Surv(max_op_hours, event_offline) ~ 1, data = dd)
    tab <- summary(fit, rmean = tau)$table
    vals <- as.numeric(tab)
    nm <- names(tab)

    data.table(
      stations = nrow(dd),
      observed_failures = sum(dd$event_offline, na.rm = TRUE),
      rmst_hours = vals[grep("rmean", nm)[1]],
      rmst_se_hours = vals[grep("se\\(rmean\\)", nm)[1]],
      tau_hours = tau
    )
  }

  rmst_region <- rbindlist(lapply(split(station_life, by = "region", keep.by = TRUE), function(dd) {
    out <- rmst_one(dd)
    out[, region := dd$region[1]]
    out
  }), fill = TRUE)

  rmst_region[, `:=`(
    rmst_years = rmst_hours / (24 * 365.25),
    rmst_se_years = rmst_se_hours / (24 * 365.25),
    tau_years = tau_hours / (24 * 365.25)
  )]

  fwrite(rmst_region, file.path(OUTDIR, paste0("rmst_lifetime_", label, ".csv")))

  # Hazard by age
  max_follow <- ceiling(max(station_life$max_op_hours, na.rm = TRUE) / HAZARD_BIN_HOURS) * HAZARD_BIN_HOURS
  haz_breaks <- seq(0, max_follow, by = HAZARD_BIN_HOURS)

  haz_region <- rbindlist(lapply(sort(unique(station_life$region)), function(rg) {
    dd <- station_life[region == rg]

    rbindlist(lapply(seq_len(length(haz_breaks) - 1), function(i) {
      start <- haz_breaks[i]
      end <- haz_breaks[i + 1]

      exposure <- pmax(0, pmin(dd$max_op_hours, end) - start)
      events <- sum(
        dd$event_offline == 1 &
          dd$max_op_hours >= start &
          dd$max_op_hours < end,
        na.rm = TRUE
      )
      at_risk <- sum(dd$max_op_hours >= start, na.rm = TRUE)

      data.table(
        region = rg,
        age_start_hours = start,
        age_end_hours = end,
        age_mid_hours = (start + end) / 2,
        exposure_hours = sum(exposure, na.rm = TRUE),
        events = events,
        at_risk = at_risk
      )
    }))
  }), fill = TRUE)

  haz_region[, hazard_per_1000_sensor_hours := fifelse(
    exposure_hours > 0,
    1000 * events / exposure_hours,
    NA_real_
  )]

  haz_region[, age_mid_years := age_mid_hours / (24 * 365.25)]
  fwrite(haz_region, file.path(OUTDIR, paste0("hazard_by_region_", label, ".csv")))

  # Effective lifetime
  flag_bin <- daily[, .(
    flagged_n = sum(n_flagged, na.rm = TRUE),
    total_n = sum(n, na.rm = TRUE)
  ), by = .(station, region, op_bin)]

  setorder(flag_bin, station, op_bin)

  flag_bin[, `:=`(
    cum_flagged_n = cumsum(flagged_n),
    cum_total_n = cumsum(total_n),
    cum_flagged_prop = cumsum(flagged_n) / cumsum(total_n)
  ), by = station]

  station_max <- station_life[, .(station, region, max_op_hours, event_offline)]

  effective_thresholds <- rbindlist(lapply(EFFECTIVE_THRESHOLDS, function(thr) {
    crossing <- flag_bin[
      cum_flagged_prop >= thr,
      .SD[1L],
      by = station
    ][, .(station, region, first_cross_hours = op_bin)]

    out <- merge(station_max, crossing, by = c("station", "region"), all.x = TRUE)

    out[, enough_followup_after_cross := fifelse(
      is.finite(first_cross_hours),
      (max_op_hours - first_cross_hours) >= MIN_POST_HOURS_EFFECTIVE,
      FALSE
    )]

    out[, `:=`(
      time_to_degrade_hours = fifelse(enough_followup_after_cross, first_cross_hours, max_op_hours),
      event_degrade = fifelse(enough_followup_after_cross, 1L, 0L),
      threshold = thr
    )]

    out
  }), fill = TRUE)

  effective_km <- rbindlist(lapply(EFFECTIVE_THRESHOLDS, function(thr) {
    dd <- effective_thresholds[threshold == thr]
    fit <- survfit(Surv(time_to_degrade_hours, event_degrade) ~ region, data = dd)
    tab <- summary(fit)$table

    out <- data.table(
      region = sub("^region=", "", rownames(tab)),
      threshold = thr,
      median_effective_hours = tab[, "median"],
      median_effective_lcl_hours = tab[, "0.95LCL"],
      median_effective_ucl_hours = tab[, "0.95UCL"]
    )

    out[, `:=`(
      median_effective_years = median_effective_hours / (24 * 365.25),
      median_effective_lcl_years = median_effective_lcl_hours / (24 * 365.25),
      median_effective_ucl_years = median_effective_ucl_hours / (24 * 365.25)
    )]

    out
  }), fill = TRUE)

  fwrite(effective_thresholds, file.path(OUTDIR, paste0("effective_lifetime_station_", label, ".csv")))
  fwrite(effective_km, file.path(OUTDIR, paste0("effective_lifetime_km_", label, ".csv")))

  # Drift coupling
  drift <- daily[, .(
    n_days = .N,
    slope_pm_geos_abs_per_hour = calc_slope(op_bin, pa_geos_abs),
    slope_temp_div_per_hour = calc_slope(op_bin, temp_div),
    slope_rh_div_per_hour = calc_slope(op_bin, rh_div),
    slope_flagged_per_hour = calc_slope(op_bin, flagged_rate)
  ), by = .(station, region)]

  drift[, `:=`(
    slope_pm_geos_abs_per_year = slope_pm_geos_abs_per_hour * 24 * 365.25,
    slope_temp_div_per_year = slope_temp_div_per_hour * 24 * 365.25,
    slope_rh_div_per_year = slope_rh_div_per_hour * 24 * 365.25,
    slope_flagged_per_year = slope_flagged_per_hour * 24 * 365.25
  )]

  fwrite(drift, file.path(OUTDIR, paste0("drift_by_sensor_", label, ".csv")))

  drift_fit <- tryCatch(
    lm(
      slope_pm_geos_abs_per_year ~ slope_temp_div_per_year + slope_rh_div_per_year + region,
      data = drift[
        is.finite(slope_pm_geos_abs_per_year) &
          is.finite(slope_temp_div_per_year) &
          is.finite(slope_rh_div_per_year)
      ]
    ),
    error = function(e) NULL
  )

  if (!is.null(drift_fit)) {
    capture.output(
      print(summary(drift_fit)),
      file = file.path(OUTDIR, paste0("drift_link_lm_", label, ".txt"))
    )
  }

  # Pre-failure instability
  prefail_daily <- merge(
    daily,
    station_life[, .(station, last_ts, event_offline, max_op_hours)],
    by = "station",
    all.x = TRUE
  )

  prefail_daily[, time_to_end_days := (max_op_hours - op_bin) / 24]

  prefail_daily[, period_flag := fifelse(
    time_to_end_days <= 30,
    "last30",
    fifelse(time_to_end_days > 30 & time_to_end_days <= 180, "baseline30_180", NA_character_)
  )]

  prefail_station <- prefail_daily[
    event_offline == 1L & !is.na(period_flag),
    .(
      mean_abs_ab = mean(ab_abs, na.rm = TRUE),
      mean_flagged_rate = mean(flagged_rate, na.rm = TRUE),
      mean_pa_geos_abs = mean(pa_geos_abs, na.rm = TRUE)
    ),
    by = .(station, region, period_flag)
  ]

  prefail_wide <- dcast(
    prefail_station,
    station + region ~ period_flag,
    value.var = c("mean_abs_ab", "mean_flagged_rate", "mean_pa_geos_abs")
  )

  prefail_wide[, `:=`(
    change_abs_ab = mean_abs_ab_last30 - mean_abs_ab_baseline30_180,
    change_flagged_rate = mean_flagged_rate_last30 - mean_flagged_rate_baseline30_180,
    change_pa_geos_abs = mean_pa_geos_abs_last30 - mean_pa_geos_abs_baseline30_180
  )]

  fwrite(prefail_wide, file.path(OUTDIR, paste0("prefailure_instability_", label, ".csv")))

  # Figures
  km_sum <- summary(km)
  km_plot_dt <- data.table(
    time_hours = km_sum$time,
    surv = km_sum$surv,
    lower = km_sum$lower,
    upper = km_sum$upper,
    strata = km_sum$strata
  )
  km_plot_dt[, region := sub("^region=", "", strata)]

  p_km <- ggplot(km_plot_dt, aes(time_hours / (24 * 365.25), surv, color = region, fill = region)) +
    geom_step(linewidth = 0.9) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, color = NA) +
    theme_pub +
    labs(
      title = paste0("Kaplan-Meier sensor survival by region: ", label),
      x = "Observed lifetime years",
      y = "Survival probability",
      color = "Region",
      fill = "Region"
    )

  ggsave(
    file.path(FIGDIR, paste0("km_survival_", label, ".jpeg")),
    p_km,
    width = 9,
    height = 6,
    dpi = 500
  )

  haz_plot <- copy(haz_region)
  haz_plot[at_risk < MIN_AT_RISK_HAZARD, hazard_per_1000_sensor_hours := NA_real_]

  p_haz <- ggplot(
    haz_plot[is.finite(hazard_per_1000_sensor_hours)],
    aes(age_mid_years, hazard_per_1000_sensor_hours, color = region)
  ) +
    geom_line(linewidth = 0.9) +
    theme_pub +
    labs(
      title = paste0("Failure hazard over operational age: ", label),
      subtitle = paste0("Bins with fewer than ", MIN_AT_RISK_HAZARD, " sensors at risk omitted"),
      x = "Operational age years",
      y = "Hazard per 1,000 sensor-hours",
      color = "Region"
    )

  ggsave(
    file.path(FIGDIR, paste0("hazard_by_region_", label, ".jpeg")),
    p_haz,
    width = 9,
    height = 6,
    dpi = 500
  )

  invisible(list(
    station_life = station_life,
    km_region = km_region,
    rmst_region = rmst_region,
    haz_region = haz_region,
    effective_km = effective_km,
    drift = drift,
    prefail = prefail_wide
  ))
}

# ------------------------------------------------------------
# 7) BUILD DATA + RUN ALL ANALYSES
# ------------------------------------------------------------

say("Building global/regional dataset")
global_obj <- make_daily_dataset(global_rdata, require_climate = FALSE)
daily_global <- add_future_outcomes(global_obj$daily)

lifecycle_global <- run_lifecycle_analysis(
  daily = daily_global,
  station_life = global_obj$station_life,
  label = "global_region"
)

run_highpm_analysis(daily_global, "global_region")

rm(daily_global, global_obj)
gc()

say("Building climate dataset")
climate_obj <- make_daily_dataset(climate_rdata, require_climate = TRUE)
daily_climate <- add_future_outcomes(climate_obj$daily)

lifecycle_climate <- run_lifecycle_analysis(
  daily = daily_climate,
  station_life = climate_obj$station_life,
  label = "single_zone_climate"
)

run_highpm_analysis(daily_climate, "climate")
run_highpm_analysis(daily_climate, "region_climate")

rm(daily_climate, climate_obj)
gc()

# ------------------------------------------------------------
# 8) CLEAN HIGH-PM RESULTS
# ------------------------------------------------------------

model_dt <- fread(out_model_csv)

model_dt[, threshold_label := paste0(">", threshold, " µg/m³")]
model_dt[, window_label := paste0(past_window_days, "-day prior window")]
model_dt[, p_value_fmt := fmt_p(p_value)]

model_dt[, estimate_95ci := sprintf(
  "%.4f (%.4f, %.4f)",
  estimate,
  lwr,
  upr
)]

model_dt[, estimate_95ci_p := paste0(
  estimate_95ci,
  "; p = ",
  p_value_fmt
)]

model_dt[, exposure_class := fifelse(
  term == "past_excess_pm_sum",
  "Cumulative burden",
  fifelse(
    term == "past_max_excess_pm",
    "Sharp spike",
    fifelse(term == "past_high_prop", "Frequency above threshold", term)
  )
)]

fwrite(model_dt, file.path(OUTDIR, "high_pm_future_degradation_model_results_clean.csv"))

joint_parse <- model_dt[
  model_name == "joint_cumulative_spike" &
    term %in% c("past_excess_pm_sum", "past_max_excess_pm"),
  .(
    analysis_type,
    analysis_group,
    group_var,
    outcome_key,
    outcome_label,
    threshold,
    past_window_days,
    term,
    exposure_class,
    scale_unit,
    effect_type,
    estimate,
    lwr,
    upr,
    p_value,
    estimate_95ci_p,
    n_rows,
    n_stations
  )
]

joint_wide <- dcast(
  joint_parse,
  analysis_type + analysis_group + group_var + outcome_label + threshold + past_window_days ~ term,
  value.var = c("estimate", "p_value", "estimate_95ci_p")
)

fwrite(joint_parse, file.path(OUTDIR, "joint_models_cumulative_vs_spike_long.csv"))
fwrite(joint_wide, file.path(OUTDIR, "joint_models_cumulative_vs_spike_wide.csv"))

# ------------------------------------------------------------
# 9) HIGH-PM FIGURES
# ------------------------------------------------------------

model_dt[, model_name := factor(
  model_name,
  levels = c("frequency_only", "cumulative_only", "spike_only", "joint_cumulative_spike")
)]

model_dt[, exposure_class := factor(
  exposure_class,
  levels = c("Frequency above threshold", "Cumulative burden", "Sharp spike")
)]

p_global <- ggplot(
  model_dt[
    analysis_type == "global_region" &
      analysis_group == "Global" &
      outcome_key %in% c("future_flagged", "future_ab_abs", "future_ab_rel")
  ],
  aes(threshold_label, estimate, ymin = lwr, ymax = upr, color = exposure_class)
) +
  geom_hline(
    aes(yintercept = ifelse(effect_type == "odds_ratio", 1, 0)),
    linetype = 2,
    color = "grey50"
  ) +
  geom_pointrange(position = position_dodge(width = 0.55)) +
  facet_grid(outcome_label ~ model_name + window_label, scales = "free_y") +
  theme_pub +
  labs(
    title = "Global high-PM exposure and future degradation",
    subtitle = "Separate and joint models parse cumulative burden from sharp spikes",
    x = "PM2.5 threshold",
    y = "Effect estimate",
    color = NULL
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(FIGDIR, "fig_global_cumulative_vs_spike.jpeg"),
  p_global,
  width = 18,
  height = 10,
  dpi = 500
)

p_region <- ggplot(
  model_dt[
    analysis_type == "global_region" &
      analysis_group != "Global" &
      outcome_key %in% c("future_flagged", "future_ab_abs", "future_ab_rel") &
      model_name == "joint_cumulative_spike" &
      past_window_days == 28
  ],
  aes(threshold_label, estimate, ymin = lwr, ymax = upr, color = exposure_class)
) +
  geom_hline(
    aes(yintercept = ifelse(effect_type == "odds_ratio", 1, 0)),
    linetype = 2,
    color = "grey50"
  ) +
  geom_pointrange(position = position_dodge(width = 0.55)) +
  facet_grid(outcome_label ~ analysis_group, scales = "free_y") +
  theme_pub +
  labs(
    title = "Regional cumulative burden vs sharp PM spikes",
    subtitle = "Joint models, 28-day exposure window",
    x = "PM2.5 threshold",
    y = "Effect estimate",
    color = NULL
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(FIGDIR, "fig_region_joint_cumulative_vs_spike_28day.jpeg"),
  p_region,
  width = 16,
  height = 9,
  dpi = 500
)

p_climate <- ggplot(
  model_dt[
    analysis_type == "climate" &
      outcome_key %in% c("future_flagged", "future_ab_abs", "future_ab_rel") &
      model_name == "joint_cumulative_spike" &
      past_window_days == 28
  ],
  aes(threshold_label, estimate, ymin = lwr, ymax = upr, color = exposure_class)
) +
  geom_hline(
    aes(yintercept = ifelse(effect_type == "odds_ratio", 1, 0)),
    linetype = 2,
    color = "grey50"
  ) +
  geom_pointrange(position = position_dodge(width = 0.55)) +
  facet_grid(outcome_label ~ analysis_group, scales = "free_y") +
  theme_pub +
  labs(
    title = "Climate-zone cumulative burden vs sharp PM spikes",
    subtitle = "Joint models, 28-day exposure window; valid single-zone Köppen stations only",
    x = "PM2.5 threshold",
    y = "Effect estimate",
    color = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(size = 6)
  )

ggsave(
  file.path(FIGDIR, "fig_climate_joint_cumulative_vs_spike_28day.jpeg"),
  p_climate,
  width = 22,
  height = 10,
  dpi = 500,
  limitsize = FALSE
)

# ------------------------------------------------------------
# 10) EMPIRICAL SI FIGURES
# ------------------------------------------------------------

if (file.exists(out_emp_csv)) {
  emp_dt <- fread(out_emp_csv)

  emp_dt[, threshold_label := factor(
    paste0(">", threshold, " µg/m³"),
    levels = c(">35 µg/m³", ">55 µg/m³", ">100 µg/m³", ">150 µg/m³")
  )]

  emp_dt[, window_label := factor(
    paste0("Past ", past_window_days, " days"),
    levels = c("Past 7 days", "Past 14 days", "Past 28 days")
  )]

  combo_counts <- emp_dt[, .(
    n_bins = .N,
    total_rows = sum(rows, na.rm = TRUE),
    total_stations = sum(stations, na.rm = TRUE)
  ), by = .(analysis_group, group_var, exposure_label, threshold, threshold_label, window_label)]

  fwrite(combo_counts, file.path(OUTDIR, "empirical_support_table.csv"))

  emp_plot <- emp_dt[
    combo_counts[n_bins >= 3],
    on = .(analysis_group, group_var, exposure_label, threshold, threshold_label, window_label)
  ]

  p_emp_global <- ggplot(
    emp_plot[
      analysis_group == "Global" &
        past_window_days == 28 &
        threshold %in% c(35, 55, 100, 150)
    ],
    aes(mean_exposure, outcome_mean)
  ) +
    geom_line() +
    geom_point(aes(size = stations), alpha = 0.85) +
    facet_grid(outcome_label ~ exposure_label + threshold_label, scales = "free") +
    theme_pub +
    labs(
      title = "Global empirical future degradation by prior high-PM exposure",
      subtitle = "28-day prior window; points are exposure quantile bins",
      x = "Mean exposure in bin",
      y = "Future outcome mean",
      size = "Stations"
    )

  ggsave(
    file.path(FIGDIR, "SI_empirical_global_28day_all_exposures.jpeg"),
    p_emp_global,
    width = 16,
    height = 9,
    dpi = 500
  )

  p_emp_region <- ggplot(
    emp_plot[
      group_var == "region" &
        analysis_group != "Global" &
        past_window_days == 28 &
        threshold %in% c(35, 55)
    ],
    aes(mean_exposure, outcome_mean, color = analysis_group)
  ) +
    geom_line() +
    geom_point(aes(size = stations), alpha = 0.85) +
    facet_grid(outcome_label ~ exposure_label + threshold_label, scales = "free") +
    theme_pub +
    labs(
      title = "Regional empirical future degradation by prior high-PM exposure",
      subtitle = "28-day prior window; thresholds >35 and >55 shown for readability",
      x = "Mean exposure in bin",
      y = "Future outcome mean",
      color = "Region",
      size = "Stations"
    )

  ggsave(
    file.path(FIGDIR, "SI_empirical_region_28day_threshold35_55.jpeg"),
    p_emp_region,
    width = 16,
    height = 9,
    dpi = 500
  )
}

# ------------------------------------------------------------
# 11) PRINT KEY RESULTS
# ------------------------------------------------------------

cat("\nGLOBAL JOINT MODEL RESULTS, 28-day window:\n")
print(
  joint_parse[
    analysis_type == "global_region" &
      analysis_group == "Global" &
      past_window_days == 28
  ][order(outcome_label, threshold, term)]
)

cat("\nOutputs saved to:\n")
cat(OUTDIR, "\n")
