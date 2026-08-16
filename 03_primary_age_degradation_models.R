# ==========================================================================
# 11_age_trends_cf1.R
#
# Age-dependent cf_1 channel degradation: global, region, climate zone, region x climate.
# => Figure 2; Tables S3-S6.
#
# Extracted VERBATIM from purpleair_degradation2.Rmd (chunk #29 lines 9991-11579).
# Code is unmodified; see REPRODUCING.md for paths and run order.
# ==========================================================================

# ---- source chunk #29 (Rmd lines 9991-11579) ------------------------

# ============================================================
# MEMORY-EFFICIENT CF1 DEGRADATION ANALYSIS
#
# GLOBAL/REGION ANALYSES FROM:
#   outdoor_final.RData
#
# CLIMATE-ZONE ANALYSES FROM:
#   outdoor_with_koppen_timevarying.RData
#
# Important:
#   - Global/region analyses are NOT restricted by Köppen zone.
#   - Climate analyses are restricted to stations with exactly ONE valid Köppen code.
#   - GAM lines are shown as smooth trend lines only.
#   - 95% CI ribbons are empirical CIs from all PurpleAir measurements
#     in each operational-age bin, not GAM prediction intervals.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(mgcv)
  library(scales)
  library(patchwork)
})

# ------------------------------------------------------------
# 1) PATHS / SETTINGS
# ------------------------------------------------------------

global_rdata <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_final.RData"
climate_rdata <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_with_koppen_timevarying.RData"
in_transition_csv <- "/Users/priyanka/Downloads/PurpleAir_EPA/station_koppen_transition_summary.csv"

outdir <- "/Users/priyanka/Downloads/PurpleAir_EPA/cf1_degradation_split_inputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

bin_w <- 24L
year_lines <- c(2, 4, 6) * 365.25 * 24

min_sensors_per_bin <- 5L
min_rows_for_trend <- 20L
min_stations_group <- 10L
min_total_sensor_bins_region_climate <- 200L
min_unique_bins_region_climate <- 25L
ncol_climate <- 4L

plot_clip_probs <- c(0.02, 0.98)

say <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")),
      paste0(..., collapse = ""), "\n")
}

safe_file <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

# ------------------------------------------------------------
# 2) COLUMN NAMES
# ------------------------------------------------------------

a_col_default <- "pm25_cf_1_a_hourly(ug/m3)_pm25_cf_1_a"
b_col_default <- "pm25_cf_1_b_hourly(ug/m3)_pm25_cf_1_b"

koppen_col <- "koppen_code"

# ------------------------------------------------------------
# 3) LOAD + REDUCE HELPERS
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

make_cf_dt <- function(infile, require_koppen = FALSE, restrict_single_zone = FALSE) {

  say("Loading: ", infile)
  outdoor <- load_outdoor_object(infile)

  a_col <- a_col_default
  b_col <- b_col_default

  if (!(a_col %in% names(outdoor)) && "cf1_a" %in% names(outdoor)) a_col <- "cf1_a"
  if (!(b_col %in% names(outdoor)) && "cf1_b" %in% names(outdoor)) b_col <- "cf1_b"

  needed <- c("station", "region", "datetime_utc", a_col, b_col)
  if (require_koppen) needed <- c(needed, koppen_col)

  miss <- setdiff(needed, names(outdoor))
  if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))

  outdoor[, station := as.character(station)]

  if (restrict_single_zone) {
    if (!file.exists(in_transition_csv)) {
      stop("Missing station Köppen transition file: ", in_transition_csv)
    }

    station_change <- fread(in_transition_csv)
    station_change[, station := as.character(station)]

    keep_stations <- station_change[
      !is.na(first_koppen) & n_koppen_codes == 1,
      station
    ]

    outdoor <- outdoor[station %in% keep_stations]

    say("Restricted to valid single-zone stations.")
    say("Rows after single-zone restriction: ", format(nrow(outdoor), big.mark = ","))
    say("Stations after single-zone restriction: ", format(uniqueN(outdoor$station), big.mark = ","))

    rm(station_change, keep_stations)
    gc()
  }

  if (require_koppen) {
    outdoor <- outdoor[
      !is.na(get(koppen_col)) &
        get(koppen_col) != ""
    ]
  }

  keep_cols <- c("station", "region", "datetime_utc", a_col, b_col)
  if (require_koppen) keep_cols <- c(keep_cols, koppen_col)

  dt <- outdoor[, ..keep_cols]

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

  dt[, datetime_utc := NULL]

  if (require_koppen) {
    setnames(dt, koppen_col, "climate_zone")
    dt[, climate_zone := as.character(climate_zone)]
  } else {
    dt[, climate_zone := NA_character_]
  }

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

  say("Rows after CF1 filtering: ", format(nrow(dt), big.mark = ","))
  say("Stations after CF1 filtering: ", format(uniqueN(dt$station), big.mark = ","))

  setorder(dt, station, ts)

  dt[, first_ts := ts[1L], by = station]
  dt[, op_hours := as.numeric(difftime(ts, first_ts, units = "hours"))]
  dt[, op_bin := floor(op_hours / bin_w) * bin_w]
  dt[, op_years := op_bin / (365.25 * 24)]

  dt <- dt[is.finite(op_bin) & op_bin >= 0]

  dt[, `:=`(
    absdiff_cf1 = abs(cf1_a - cf1_b),
    bias_cf1 = cf1_a - cf1_b,
    pctdiff_cf1 = fifelse(
      (cf1_a + cf1_b) > 0,
      2 * abs(cf1_a - cf1_b) / (cf1_a + cf1_b),
      NA_real_
    )
  )]

  thr_dt <- dt[
    is.finite(pctdiff_cf1),
    .(thr85 = quantile(pctdiff_cf1, 0.85, na.rm = TRUE, type = 7)),
    by = station
  ]

  dt <- thr_dt[dt, on = "station"]
  dt[, flagged_cf1 := as.integer(absdiff_cf1 > 5 & pctdiff_cf1 > thr85)]

  rm(thr_dt)
  gc()

  dt[]
}

# ------------------------------------------------------------
# 4) MEMORY-EFFICIENT GLOBAL/REGION PROCESSING
# ------------------------------------------------------------

make_sens_bin <- function(dt, by_extra = character()) {

  by_vars <- c("station", "region", by_extra, "op_bin", "op_years")

  dt[, .(
    n_obs = .N,
    n_flagged = sum(flagged_cf1, na.rm = TRUE),
    mean_a_s = mean(cf1_a, na.rm = TRUE),
    mean_b_s = mean(cf1_b, na.rm = TRUE),
    mean_abs_s = mean(absdiff_cf1, na.rm = TRUE),
    mean_bias_s = mean(bias_cf1, na.rm = TRUE),
    mean_pct_s = mean(pctdiff_cf1, na.rm = TRUE)
  ), by = by_vars]
}

say("Building global/region dataset from outdoor_final.RData")

dt_global <- make_cf_dt(
  infile = global_rdata,
  require_koppen = FALSE,
  restrict_single_zone = FALSE
)

say("Aggregating global station-bin summaries")

sens_global <- make_sens_bin(dt_global)

say("Global station-bin rows: ", format(nrow(sens_global), big.mark = ","))

rm(dt_global)
gc()

# ------------------------------------------------------------
# 5) MEMORY-EFFICIENT CLIMATE PROCESSING
# ------------------------------------------------------------

say("Building climate dataset from outdoor_with_koppen_timevarying.RData")

dt_climate <- make_cf_dt(
  infile = climate_rdata,
  require_koppen = TRUE,
  restrict_single_zone = TRUE
)

if (file.exists(in_transition_csv)) {
  kop_desc <- fread(
    in_transition_csv,
    select = c("station", "first_koppen", "first_koppen_desc", "n_koppen_codes")
  )

  kop_desc[, station := as.character(station)]

  kop_desc <- kop_desc[
    !is.na(first_koppen) & n_koppen_codes == 1,
    .(
      station,
      climate_desc = as.character(first_koppen_desc)
    )
  ]

  dt_climate <- merge(dt_climate, kop_desc, by = "station", all.x = TRUE)

  rm(kop_desc)
  gc()
} else {
  dt_climate[, climate_desc := NA_character_]
}

dt_climate[
  is.na(climate_desc) | climate_desc == "",
  climate_desc := climate_zone
]

dt_climate[, climate_label := paste0(climate_zone, "\n", climate_desc)]

say("Aggregating climate station-bin summaries")

sens_climate <- make_sens_bin(
  dt_climate,
  by_extra = c("climate_zone", "climate_desc", "climate_label")
)

say("Climate station-bin rows: ", format(nrow(sens_climate), big.mark = ","))

rm(dt_climate)
gc()

# ------------------------------------------------------------
# 6) TREND TEST HELPERS + CLUSTER-ROBUST WORD-READY TABLES
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(officer)
  library(flextable)
  library(sandwich)
  library(lmtest)
})

fmt_p <- function(p) {
  ifelse(
    is.na(p), NA_character_,
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), NA_character_, sprintf(paste0("%.", digits, "f"), x))
}

fmt_ci <- function(est, lwr, upr, digits = 3) {
  paste0(
    fmt_num(est, digits), " (",
    fmt_num(lwr, digits), ", ",
    fmt_num(upr, digits), ")"
  )
}

extract_term <- function(ct, term = "op_years") {
  if (is.null(ct)) return(NULL)
  if (!(term %in% rownames(ct))) return(NULL)

  p_col <- grep("^Pr\\(", colnames(ct), value = TRUE)

  if (length(p_col) == 0L) {
    p_val <- NA_real_
  } else {
    p_val <- as.numeric(ct[term, p_col[1]])
  }

  list(
    estimate = as.numeric(ct[term, "Estimate"]),
    se = as.numeric(ct[term, "Std. Error"]),
    p_value = p_val
  )
}

get_cluster_test <- function(fit, dd) {
  vc <- sandwich::vcovCL(
    fit,
    cluster = dd$station,
    type = "HC1"
  )

  lmtest::coeftest(fit, vcov. = vc)
}

tidy_lm_trend <- function(dd, yvar, group_name, outcome_label) {
  dd <- dd[
    is.finite(get(yvar)) &
      is.finite(op_years) &
      is.finite(n_obs) &
      n_obs > 0 &
      !is.na(station)
  ]

  if (nrow(dd) < min_rows_for_trend || uniqueN(dd$station) < min_stations_group) {
    return(data.table())
  }

  fit <- tryCatch(
    lm(
      as.formula(paste0(yvar, " ~ op_years")),
      data = dd,
      weights = n_obs
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) return(data.table())

  ct <- tryCatch(
    get_cluster_test(fit, dd),
    error = function(e) NULL
  )

  ext <- extract_term(ct, "op_years")
  if (is.null(ext)) return(data.table())

  beta <- ext$estimate
  se <- ext$se

  data.table(
    group = group_name,
    outcome = outcome_label,
    model = "Weighted linear regression with sensor-clustered SEs",
    effect_type = "Beta per operational year",
    n_sensor_bins = nrow(dd),
    n_sensors = uniqueN(dd$station),
    min_year = min(dd$op_years, na.rm = TRUE),
    max_year = max(dd$op_years, na.rm = TRUE),
    estimate = beta,
    ci_lwr = beta - 1.96 * se,
    ci_upr = beta + 1.96 * se,
    se_cluster = se,
    p_value = ext$p_value
  )
}

tidy_flag_trend <- function(dd, group_name) {
  dd <- dd[
    is.finite(op_years) &
      is.finite(n_flagged) &
      is.finite(n_obs) &
      n_obs > 0 &
      n_flagged >= 0 &
      n_flagged <= n_obs &
      !is.na(station)
  ]

  if (nrow(dd) < min_rows_for_trend || uniqueN(dd$station) < min_stations_group) {
    return(data.table())
  }

  fit <- tryCatch(
    glm(
      cbind(n_flagged, n_obs - n_flagged) ~ op_years,
      data = dd,
      family = quasibinomial()
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) return(data.table())

  ct <- tryCatch(
    get_cluster_test(fit, dd),
    error = function(e) NULL
  )

  ext <- extract_term(ct, "op_years")
  if (is.null(ext)) return(data.table())

  beta <- ext$estimate
  se <- ext$se

  data.table(
    group = group_name,
    outcome = "Flagged degradation",
    model = "Quasibinomial regression with sensor-clustered SEs",
    effect_type = "Odds ratio per operational year",
    n_sensor_bins = nrow(dd),
    n_sensors = uniqueN(dd$station),
    min_year = min(dd$op_years, na.rm = TRUE),
    max_year = max(dd$op_years, na.rm = TRUE),
    estimate = exp(beta),
    ci_lwr = exp(beta - 1.96 * se),
    ci_upr = exp(beta + 1.96 * se),
    log_or = beta,
    se_cluster = se,
    p_value = ext$p_value
  )
}

run_trend_set <- function(dd, group_label) {
  rbindlist(
    list(
      tidy_lm_trend(dd, "mean_abs_s",  group_label, "Absolute channel divergence |A-B|"),
      tidy_lm_trend(dd, "mean_bias_s", group_label, "Signed channel bias A-B"),
      tidy_lm_trend(dd, "mean_pct_s",  group_label, "Relative channel divergence"),
      tidy_lm_trend(dd, "mean_a_s",    group_label, "Channel A concentration"),
      tidy_lm_trend(dd, "mean_b_s",    group_label, "Channel B concentration"),
      tidy_flag_trend(dd, group_label)
    ),
    use.names = TRUE,
    fill = TRUE
  )
}

make_word_table <- function(dt, title, outfile) {
  tab <- copy(dt)

  tab[, `:=`(
    `Estimate (95% CI)` = fmt_ci(estimate, ci_lwr, ci_upr, digits = 3),
    `P value` = fmt_p(p_value),
    `Operational years` = paste0(round(min_year, 2), "-", round(max_year, 2))
  )]

  tab <- tab[
    , .(
      Group = group,
      Outcome = outcome,
      Model = model,
      `Effect reported` = effect_type,
      `Estimate (95% CI)`,
      `P value`,
      `N sensors` = n_sensors,
      `N sensor-bins` = n_sensor_bins,
      `Operational years`
    )
  ]

  ft <- flextable(tab)
  ft <- autofit(ft)
  ft <- fontsize(ft, size = 8, part = "all")
  ft <- bold(ft, part = "header")
  ft <- align(ft, align = "center", part = "header")
  ft <- align(
    ft,
    j = c("Estimate (95% CI)", "P value", "N sensors", "N sensor-bins", "Operational years"),
    align = "center",
    part = "body"
  )
  ft <- set_table_properties(ft, width = 1, layout = "autofit")

  doc <- read_docx()
  doc <- body_add_par(doc, title, style = "heading 1")
  doc <- body_add_par(
    doc,
    paste0(
      "Continuous degradation outcomes are reported as beta coefficients per operational year. ",
      "Flagged degradation is reported as an odds ratio per operational year from quasibinomial regression. ",
      "Confidence intervals are 95% intervals calculated using sensor-clustered robust standard errors."
    ),
    style = "Normal"
  )
  doc <- body_add_flextable(doc, ft)

  print(doc, target = outfile)
}

# ------------------------------------------------------------
# 7) RUN CLUSTER-ROBUST TREND TESTS
# ------------------------------------------------------------

say("Running cluster-robust trend tests for associations between degradation metrics and operational age")

trend_global <- run_trend_set(sens_global, "Global")

trend_region <- sens_global[
  , run_trend_set(.SD, paste0("Region: ", unique(region))),
  by = region
][, region := NULL]

trend_climate <- sens_climate[
  , run_trend_set(.SD, paste0("Climate: ", unique(climate_zone), " - ", unique(climate_desc))),
  by = .(climate_zone, climate_desc)
][, c("climate_zone", "climate_desc") := NULL]

trend_region_climate <- sens_climate[
  , run_trend_set(
    .SD,
    paste0("Region x climate: ", unique(region), " / ", unique(climate_zone), " - ", unique(climate_desc))
  ),
  by = .(region, climate_zone, climate_desc)
][, c("region", "climate_zone", "climate_desc") := NULL]

trend_all <- rbindlist(
  list(trend_global, trend_region, trend_climate, trend_region_climate),
  use.names = TRUE,
  fill = TRUE
)

trend_table_clean <- copy(trend_all)[
  , `:=`(
    estimate_95ci = fmt_ci(estimate, ci_lwr, ci_upr, digits = 3),
    p_value_fmt = fmt_p(p_value),
    estimate_95ci_p = paste0(fmt_ci(estimate, ci_lwr, ci_upr, digits = 3), "; p=", fmt_p(p_value))
  )
]

fwrite(trend_all, file.path(outdir, "degradation_time_associations_all_groups_cluster_robust_raw.csv"))
fwrite(trend_table_clean, file.path(outdir, "degradation_time_associations_all_groups_cluster_robust_clean.csv"))
fwrite(trend_global, file.path(outdir, "degradation_time_associations_global_cluster_robust.csv"))
fwrite(trend_region, file.path(outdir, "degradation_time_associations_by_region_cluster_robust.csv"))
fwrite(trend_climate, file.path(outdir, "degradation_time_associations_by_climate_zone_cluster_robust.csv"))
fwrite(trend_region_climate, file.path(outdir, "degradation_time_associations_by_region_x_climate_zone_cluster_robust.csv"))



# ------------------------------------------------------------
# 8) SENSOR-LEVEL 95% CI CURVES
# ------------------------------------------------------------

mean_ci_sensor <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)

  if (n == 0L) {
    return(c(mean = NA_real_, lwr = NA_real_, upr = NA_real_, n = 0))
  }

  mu <- mean(x)
  se <- if (n > 1L) sd(x) / sqrt(n) else 0

  c(
    mean = mu,
    lwr = mu - 1.96 * se,
    upr = mu + 1.96 * se,
    n = n
  )
}

prop_ci_sensor <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)

  if (n == 0L) {
    return(c(mean = NA_real_, lwr = NA_real_, upr = NA_real_, n = 0))
  }

  mu <- mean(x)
  se <- if (n > 1L) sd(x) / sqrt(n) else 0

  c(
    mean = mu,
    lwr = max(0, mu - 1.96 * se),
    upr = min(1, mu + 1.96 * se),
    n = n
  )
}

summarize_curve_sensor_ci <- function(sens_dt, by_vars) {
  sens_dt[, {
    a  <- mean_ci_sensor(mean_a_s)
    b  <- mean_ci_sensor(mean_b_s)
    ab <- mean_ci_sensor(mean_abs_s)
    bi <- mean_ci_sensor(mean_bias_s)
    pc <- mean_ci_sensor(mean_pct_s)
    fl <- prop_ci_sensor(n_flagged / pmax(n_obs, 1))

    .(
      n_sensors = as.numeric(ab["n"]),
      n_observations = sum(n_obs, na.rm = TRUE),

      mean_a = a["mean"], lwr_a = a["lwr"], upr_a = a["upr"],
      mean_b = b["mean"], lwr_b = b["lwr"], upr_b = b["upr"],

      mean_abs = ab["mean"], lwr_abs = ab["lwr"], upr_abs = ab["upr"],
      mean_bias = bi["mean"], lwr_bias = bi["lwr"], upr_bias = bi["upr"],
      mean_pct = pc["mean"], lwr_pct = pc["lwr"], upr_pct = pc["upr"],
      pct_flagged = fl["mean"], lwr_flag = fl["lwr"], upr_flag = fl["upr"]
    )
  }, by = by_vars]
}

say("Building sensor-level 95% CI curves")

curve_global <- summarize_curve_sensor_ci(sens_global, c("op_bin", "op_years"))
curve_global[, region := "Global"]
setcolorder(curve_global, c("region", setdiff(names(curve_global), "region")))

curve_region <- summarize_curve_sensor_ci(
  sens_global,
  c("region", "op_bin", "op_years")
)

curve_region_plus <- rbindlist(
  list(curve_global, curve_region),
  use.names = TRUE,
  fill = TRUE
)

curve_region_plus <- curve_region_plus[n_sensors >= min_sensors_per_bin]

fwrite(
  curve_region_plus,
  file.path(outdir, "cf1_lifetime_curves_global_and_region_sensor_ci.csv")
)

rm(curve_global, curve_region, sens_global)
gc()

curve_climate <- summarize_curve_sensor_ci(
  sens_climate,
  c("climate_zone", "climate_desc", "climate_label", "op_bin", "op_years")
)

curve_region_climate <- summarize_curve_sensor_ci(
  sens_climate,
  c("region", "climate_zone", "climate_desc", "climate_label", "op_bin", "op_years")
)

curve_climate <- curve_climate[n_sensors >= min_sensors_per_bin]
curve_region_climate <- curve_region_climate[n_sensors >= min_sensors_per_bin]

clim_order <- curve_climate[
  , .(total_sensor_bins = sum(n_sensors, na.rm = TRUE)),
  by = climate_label
][order(-total_sensor_bins)]$climate_label

curve_climate[, climate_label := factor(climate_label, levels = clim_order)]
curve_region_climate[, climate_label := factor(climate_label, levels = clim_order)]

fwrite(
  curve_climate,
  file.path(outdir, "cf1_lifetime_curves_by_climate_zone_sensor_ci.csv")
)

fwrite(
  curve_region_climate,
  file.path(outdir, "cf1_lifetime_curves_by_region_x_climate_zone_sensor_ci.csv")
)

rm(sens_climate)
gc()

# ------------------------------------------------------------
# 9) REGION × CLIMATE FILTER
# ------------------------------------------------------------

valid_region_climate <- curve_region_climate[
  , .(
    total_sensor_bins = sum(n_sensors, na.rm = TRUE),
    n_bins = uniqueN(op_bin)
  ),
  by = .(region, climate_label)
][
  total_sensor_bins >= min_total_sensor_bins_region_climate &
    n_bins >= min_unique_bins_region_climate
]

curve_region_climate_plot <- merge(
  curve_region_climate,
  valid_region_climate[, .(region, climate_label)],
  by = c("region", "climate_label"),
  all = FALSE
)

fwrite(
  curve_region_climate_plot,
  file.path(outdir, "cf1_lifetime_curves_region_x_climate_PLOTTED_sensor_ci.csv")
)

say("Region x climate panels retained: ", nrow(valid_region_climate))

rm(valid_region_climate)
gc()

# ------------------------------------------------------------
# 10) GAM SMOOTH HELPER
# ------------------------------------------------------------

make_gam_predictions <- function(d, yvar, by_vars, step = 168) {
  d <- copy(d)
  d <- d[
    is.finite(get(yvar)) &
      is.finite(op_bin) &
      is.finite(n_sensors) &
      n_sensors > 0
  ]

  if (nrow(d) == 0L) return(data.table())

  d[, group_key := do.call(paste, c(.SD, sep = "___")), .SDcols = by_vars]

  out <- d[, {
    dd <- .SD[
      is.finite(get(yvar)) &
        is.finite(op_bin) &
        is.finite(n_sensors) &
        n_sensors > 0
    ]

    if (.N < 25L || uniqueN(op_bin) < 25L) {
      data.table()
    } else {
      fit <- tryCatch(
        gam(
          as.formula(sprintf("%s ~ s(op_bin, bs = 'ps')", yvar)),
          data = dd,
          weights = n_sensors,
          method = "REML"
        ),
        error = function(e) NULL
      )

      if (is.null(fit)) {
        data.table()
      } else {
        newd <- data.table(op_bin = seq(min(dd$op_bin), max(dd$op_bin), by = step))
        pr <- predict(fit, newdata = newd, se.fit = FALSE)

        newd[, `:=`(
          fit = as.numeric(pr),
          variable = yvar
        )]

        for (bv in by_vars) newd[, (bv) := dd[[bv]][1]]
        newd
      }
    }
  }, by = group_key]

  out[, group_key := NULL]
  out[]
}

# ------------------------------------------------------------
# 11) GAM SMOOTH HELPER
#     GAM is trend line only; no GAM CI is plotted.
# ------------------------------------------------------------

make_gam_predictions <- function(d, yvar, by_vars, step = 168) {
  d <- copy(d)
  d <- d[
    is.finite(get(yvar)) &
      is.finite(op_bin) &
      is.finite(n_sensors) &
      n_sensors > 0
  ]

  if (nrow(d) == 0L) return(data.table())

  d[, group_key := do.call(paste, c(.SD, sep = "___")), .SDcols = by_vars]

  out <- d[, {
    dd <- .SD[
      is.finite(get(yvar)) &
        is.finite(op_bin) &
        is.finite(n_sensors) &
        n_sensors > 0
    ]

    if (.N < 25L || uniqueN(op_bin) < 25L) {
      data.table()
    } else {
      fit <- tryCatch(
        gam(
          as.formula(sprintf("%s ~ s(op_bin, bs = 'ps')", yvar)),
          data = dd,
          weights = n_sensors,
          method = "REML"
        ),
        error = function(e) NULL
      )

      if (is.null(fit)) {
        data.table()
      } else {
        newd <- data.table(op_bin = seq(min(dd$op_bin), max(dd$op_bin), by = step))
        pr <- predict(fit, newdata = newd, se.fit = FALSE)

        newd[, `:=`(
          fit = as.numeric(pr),
          variable = yvar
        )]

        for (bv in by_vars) newd[, (bv) := dd[[bv]][1]]
        newd
      }
    }
  }, by = group_key]

  out[, group_key := NULL]
  out[]
}

# ------------------------------------------------------------
# 12) PLOTTING HELPERS
#     Taller wrapped grids + purple n_sensors
# ------------------------------------------------------------

make_panel_labels <- function(n) {
  paste0(
    vapply(seq_len(n), function(i) {
      if (i <= 26) LETTERS[i] else paste0(
        LETTERS[((i - 1) %/% 26)],
        LETTERS[((i - 1) %% 26) + 1]
      )
    }, character(1)),
    ")"
  )
}

auto_ncol <- function(n_panels) {
  if (n_panels <= 6) return(n_panels)
  if (n_panels <= 12) return(4L)
  return(5L)
}

theme_pub <- theme_minimal(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 9),
    strip.text = element_text(face = "bold", size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey87", linewidth = 0.25),
    legend.position = "bottom"
  )

theme_n <- theme_pub +
  theme(
    strip.text = element_blank(),
    axis.title.y = element_text(color = "purple4", face = "bold"),
    axis.text.y = element_text(color = "purple4"),
    axis.ticks.y = element_line(color = "purple4")
  )

add_lifetime_x <- function(p) {
  p +
    geom_vline(
      xintercept = year_lines,
      linetype = 2,
      linewidth = 0.25,
      alpha = 0.45
    ) +
    scale_x_continuous(labels = comma) +
    labs(x = "Operational hour")
}

clip_metric_for_plot <- function(plot_dt, gam_dt, by_vars, y, lwr, upr,
                                 metric_type = c("positive", "signed", "proportion", "free"),
                                 probs = plot_clip_probs) {
  metric_type <- match.arg(metric_type)

  d <- copy(plot_dt)
  g <- copy(gam_dt)

  d[, n_sensors := fifelse(!is.finite(n_sensors) | n_sensors < 0, 0, n_sensors)]

  lims <- d[
    is.finite(get(y)) & is.finite(get(lwr)) & is.finite(get(upr)),
    {
      vals <- c(get(y), get(lwr), get(upr))

      if (metric_type == "signed") {
        hi <- as.numeric(quantile(abs(vals), probs[2], na.rm = TRUE, type = 7))
        if (!is.finite(hi) || hi == 0) hi <- max(abs(vals), na.rm = TRUE)
        if (!is.finite(hi) || hi == 0) hi <- 1
        .(ymin_clip = -hi, ymax_clip = hi)
      } else if (metric_type == "proportion") {
        hi <- as.numeric(quantile(vals, probs[2], na.rm = TRUE, type = 7))
        if (!is.finite(hi) || hi <= 0) hi <- max(vals, na.rm = TRUE)
        if (!is.finite(hi) || hi <= 0) hi <- 0.1
        .(ymin_clip = 0, ymax_clip = min(1, hi * 1.08))
      } else if (metric_type == "positive") {
        hi <- as.numeric(quantile(vals, probs[2], na.rm = TRUE, type = 7))
        if (!is.finite(hi) || hi <= 0) hi <- max(vals, na.rm = TRUE)
        if (!is.finite(hi) || hi <= 0) hi <- 1
        .(ymin_clip = 0, ymax_clip = hi * 1.08)
      } else {
        qs <- quantile(vals, probs, na.rm = TRUE, type = 7)
        .(ymin_clip = qs[1], ymax_clip = qs[2])
      }
    },
    by = by_vars
  ]

  d <- merge(d, lims, by = by_vars, all.x = TRUE)

  d[, y_plot := pmin(pmax(get(y), ymin_clip), ymax_clip)]
  d[, lwr_plot := pmin(pmax(get(lwr), ymin_clip), ymax_clip)]
  d[, upr_plot := pmin(pmax(get(upr), ymin_clip), ymax_clip)]

  if (nrow(g) > 0) {
    g <- merge(g, lims, by = by_vars, all.x = TRUE)
    g[, fit_plot := pmin(pmax(fit, ymin_clip), ymax_clip)]
  } else {
    g[, fit_plot := numeric()]
  }

  list(plot_dt = d, gam_dt = g)
}

make_metric_plus_n_plot <- function(
    plot_dt,
    gam_dt,
    y,
    lwr,
    upr,
    ylab,
    title,
    subtitle,
    facet_obj,
    by_vars,
    outfile,
    width,
    height,
    metric_type = "positive"
) {
  clipped <- clip_metric_for_plot(
    plot_dt = plot_dt,
    gam_dt = gam_dt,
    by_vars = by_vars,
    y = y,
    lwr = lwr,
    upr = upr,
    metric_type = metric_type
  )

  plot_dt <- clipped$plot_dt
  gam_dt <- clipped$gam_dt

  p_metric <- ggplot(plot_dt, aes(x = op_bin, y = y_plot)) +
    geom_ribbon(
      aes(ymin = lwr_plot, ymax = upr_plot),
      fill = "grey45",
      alpha = 0.42,
      color = "grey25",
      linewidth = 0.12
    ) +
    geom_line(color = "black", linewidth = 0.44, alpha = 0.98) +
    geom_line(
      data = gam_dt,
      aes(x = op_bin, y = fit_plot),
      inherit.aes = FALSE,
      color = "blue4",
      linewidth = 0.95
    ) +
    facet_obj +
    scale_y_continuous(name = ylab) +
    labs(
      title = title,
      subtitle = paste0(
        subtitle,
        " Grey ribbon = sensor-level 95% CI across sensors in each age bin; blue = GAM smooth only."
      ),
      x = NULL
    ) +
    theme_pub +
    theme(axis.text.x = element_blank()) +
    geom_vline(
      xintercept = year_lines,
      linetype = 2,
      linewidth = 0.25,
      alpha = 0.45
    )

  p_n <- ggplot(plot_dt, aes(x = op_bin, y = n_sensors)) +
    geom_line(
      color = "purple4",
      linewidth = 0.55,
      alpha = 0.95,
      linetype = "dotdash"
    ) +
    facet_obj +
    scale_y_continuous(
      name = "n sensors",
      labels = comma,
      limits = c(0, NA)
    ) +
    labs(x = "Operational hour") +
    theme_n

  p_n <- add_lifetime_x(p_n)

  p_combined <- p_metric / p_n + plot_layout(heights = c(4, 1))

  ggsave(
    file.path(outdir, outfile),
    p_combined,
    width = width,
    height = height,
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )

  rm(p_metric, p_n, p_combined, plot_dt, gam_dt, clipped)
  gc()

  invisible(NULL)
}

# Helper: use this in Sections 13, 15, and 16 to avoid ultra-wide plots.
get_facet_ncol <- function(dt, facet_var) {
  auto_ncol(uniqueN(dt[[facet_var]]))
}

get_wrapped_height <- function(dt, facet_var, base = 4, per_row = 3.2) {
  n <- uniqueN(dt[[facet_var]])
  nc <- auto_ncol(n)
  nr <- ceiling(n / nc)
  max(8, base + per_row * nr)
}

# ------------------------------------------------------------
# 13) CLIMATE-ZONE FIGURES
# ------------------------------------------------------------

say("Saving climate-zone figures")

save_climate_metric_plot <- function(y, lwr, upr, ylab, title, outfile, gam_var, metric_type) {
  gd <- make_gam_predictions(curve_climate, gam_var, "climate_label")

  make_metric_plus_n_plot(
    plot_dt = curve_climate,
    gam_dt = gd,
    y = y,
    lwr = lwr,
    upr = upr,
    ylab = ylab,
    title = title,
    subtitle = "Climate-zone analyses use only valid single-zone Köppen stations.",
    facet_obj = facet_wrap(~ climate_label, scales = "free_y", ncol = ncol_climate),
    by_vars = "climate_label",
    outfile = outfile,
    width = 16,
    height = 14,
    metric_type = metric_type
  )

  rm(gd)
  gc()
}

save_climate_metric_plot("mean_abs", "lwr_abs", "upr_abs", "|A − B| (µg/m³)",
                         "Lifetime absolute divergence by climate zone (CF1)",
                         "cf1_absolute_divergence_all_climate_zones_pub.jpeg",
                         "mean_abs", "positive")

save_climate_metric_plot("mean_bias", "lwr_bias", "upr_bias", "A − B (µg/m³)",
                         "Lifetime signed bias by climate zone (CF1)",
                         "cf1_signed_bias_all_climate_zones_pub.jpeg",
                         "mean_bias", "signed")

save_climate_metric_plot("mean_pct", "lwr_pct", "upr_pct", "Relative divergence",
                         "Lifetime relative divergence by climate zone (CF1)",
                         "cf1_relative_divergence_all_climate_zones_pub.jpeg",
                         "mean_pct", "positive")

save_climate_metric_plot(
  "pct_flagged", "lwr_flag", "upr_flag",
  "Proportion flagged",
  "Lifetime flagged rate by climate zone (CF1)",
  "cf1_flagged_rate_all_climate_zones_pub.jpeg",
  "pct_flagged",
  "proportion"
)

# ------------------------------------------------------------
# 14) REGION + GLOBAL FIGURES
# ------------------------------------------------------------

say("Saving region/global figures")

save_region_metric_plot <- function(y, lwr, upr, ylab, title, outfile, gam_var, metric_type) {
  gd <- make_gam_predictions(curve_region_plus, gam_var, "region")

  make_metric_plus_n_plot(
    plot_dt = curve_region_plus,
    gam_dt = gd,
    y = y,
    lwr = lwr,
    upr = upr,
    ylab = ylab,
    title = title,
    subtitle = "Global and regional analyses use outdoor_final.RData without Köppen restriction.",
    facet_obj = facet_wrap(~ region, scales = "free_y", nrow = 1),
    by_vars = "region",
    outfile = outfile,
    width = 16,
    height = 7.5,
    metric_type = metric_type
  )

  rm(gd)
  gc()
}

save_region_metric_plot("mean_abs", "lwr_abs", "upr_abs", "|A − B| (µg/m³)",
                        "Lifetime absolute divergence by region and globally (CF1)",
                        "cf1_absolute_divergence_region_global_pub.jpeg",
                        "mean_abs", "positive")

save_region_metric_plot("mean_bias", "lwr_bias", "upr_bias", "A − B (µg/m³)",
                        "Lifetime signed bias by region and globally (CF1)",
                        "cf1_signed_bias_region_global_pub.jpeg",
                        "mean_bias", "signed")

save_region_metric_plot("mean_pct", "lwr_pct", "upr_pct", "Relative divergence",
                        "Lifetime relative divergence by region and globally (CF1)",
                        "cf1_relative_divergence_region_global_pub.jpeg",
                        "mean_pct", "positive")

save_region_metric_plot(
  "pct_flagged", "lwr_flag", "upr_flag",
  "Proportion flagged",
  "Lifetime flagged rate by region and globally (CF1)",
  "cf1_flagged_rate_region_global_pub.jpeg",
  "pct_flagged",
  "proportion"
)

# ------------------------------------------------------------
# 15) REGION × CLIMATE FIGURES
# ------------------------------------------------------------

say("Saving region x climate figures")

save_region_climate_metric_plot <- function(
    y, lwr, upr, ylab, title_stub, outfile_stub, gam_var, metric_type
) {
  for (rr in sort(unique(curve_region_climate_plot$region))) {
    plot_dt <- curve_region_climate_plot[region == rr]
    if (nrow(plot_dt) == 0L) next

    gd <- make_gam_predictions(plot_dt, gam_var, c("region", "climate_label"))

    make_metric_plus_n_plot(
      plot_dt = plot_dt,
      gam_dt = gd,
      y = y,
      lwr = lwr,
      upr = upr,
      ylab = ylab,
      title = paste0(title_stub, ": ", rr),
      subtitle = paste0(
        "Region × climate analyses use valid single-zone Köppen stations only. Panels require ≥",
        min_total_sensor_bins_region_climate,
        " station-bins and ≥", min_unique_bins_region_climate,
        " age bins."
      ),
      facet_obj = facet_wrap(~ climate_label, scales = "free_y", ncol = 4),
      by_vars = c("region", "climate_label"),
      outfile = paste0(outfile_stub, "_", safe_file(rr), ".jpeg"),
      width = 16,
      height = 12,
      metric_type = metric_type
    )

    rm(plot_dt, gd)
    gc()
  }
}

save_region_climate_metric_plot("mean_abs", "lwr_abs", "upr_abs", "|A − B| (µg/m³)",
                                "Lifetime absolute divergence by region and climate zone (CF1)",
                                "cf1_absolute_divergence_region_x_climate_pub",
                                "mean_abs", "positive")

save_region_climate_metric_plot("mean_bias", "lwr_bias", "upr_bias", "A − B (µg/m³)",
                                "Lifetime signed bias by region and climate zone (CF1)",
                                "cf1_signed_bias_region_x_climate_pub",
                                "mean_bias", "signed")

save_region_climate_metric_plot("mean_pct", "lwr_pct", "upr_pct", "Relative divergence",
                                "Lifetime relative divergence by region and climate zone (CF1)",
                                "cf1_relative_divergence_region_x_climate_pub",
                                "mean_pct", "positive")

save_region_climate_metric_plot(
  "pct_flagged", "lwr_flag", "upr_flag",
  "Proportion flagged",
  "Lifetime flagged rate by region and climate zone (CF1)",
  "cf1_flagged_rate_region_x_climate_pub",
  "pct_flagged",
  "proportion"
)
# ------------------------------------------------------------
# 16) CHANNEL A/B FIGURES WITH SENSOR-LEVEL CI RIBBONS
# ------------------------------------------------------------

say("Saving Channel A/B figures")

make_ab_plot <- function(curve_dt, facet_formula, outfile, title, width, height, by_vars) {

  curve_ab <- melt(
    curve_dt,
    id.vars = intersect(
      c("region", "climate_zone", "climate_label", "op_bin", "n_sensors", "n_observations"),
      names(curve_dt)
    ),
    measure.vars = list(
      mean = c("mean_a", "mean_b"),
      lwr = c("lwr_a", "lwr_b"),
      upr = c("upr_a", "upr_b")
    ),
    variable.name = "channel"
  )

  curve_ab[, channel := factor(channel, levels = c(1, 2), labels = c("Channel A", "Channel B"))]

  gd_a <- make_gam_predictions(curve_dt, "mean_a", by_vars)
  gd_b <- make_gam_predictions(curve_dt, "mean_b", by_vars)

  gam_ab <- rbindlist(list(gd_a, gd_b), use.names = TRUE, fill = TRUE)

  if (nrow(gam_ab) > 0) {
    gam_ab[, channel := fifelse(variable == "mean_a", "Channel A", "Channel B")]
    gam_ab[, channel := factor(channel, levels = c("Channel A", "Channel B"))]
  }

  p <- ggplot(curve_ab, aes(x = op_bin, y = mean, color = channel, fill = channel)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.35, color = NA) +
    geom_line(aes(y = lwr), linewidth = 0.16, alpha = 0.6) +
    geom_line(aes(y = upr), linewidth = 0.16, alpha = 0.6) +
    geom_line(linewidth = 0.35, alpha = 0.95) +
    geom_line(
      data = gam_ab,
      aes(x = op_bin, y = fit, color = channel),
      inherit.aes = FALSE,
      linewidth = 0.9
    ) +
    facet_formula +
    scale_color_manual(values = c("Channel A" = "#1f78b4", "Channel B" = "#e31a1c")) +
    scale_fill_manual(values = c("Channel A" = "#1f78b4", "Channel B" = "#e31a1c")) +
    labs(
      title = title,
      subtitle = "Ribbons are sensor-level 95% CIs across sensors in each operational-age bin; smooth lines are GAM fits only.",
      y = expression("Mean CF1 PM"[2.5]~"("*mu*g/m^3*")"),
      x = "Operational hour",
      color = NULL,
      fill = NULL
    ) +
    theme_pub

  p <- add_lifetime_x(p)

  ggsave(file.path(outdir, outfile), p, width = width, height = height, dpi = 300, bg = "white")

  rm(curve_ab, gd_a, gd_b, gam_ab, p)
  gc()
}

make_ab_plot(
  curve_region_plus,
  facet_wrap(~ region, scales = "free_y", nrow = 1),
  "cf1_channel_ab_concentrations_region_global_pub.jpeg",
  "Lifetime trajectories of CF1 channel concentrations by region and globally",
  width = 16,
  height = 5.8,
  by_vars = "region"
)

make_ab_plot(
  curve_climate,
  facet_wrap(~ climate_label, scales = "free_y", ncol = ncol_climate),
  "cf1_channel_ab_concentrations_climate_pub.jpeg",
  "Lifetime trajectories of CF1 channel concentrations by climate zone",
  width = 16,
  height = 11,
  by_vars = "climate_label"
)

for (rr in sort(unique(curve_region_climate_plot$region))) {
  plot_dt <- curve_region_climate_plot[region == rr]

  make_ab_plot(
    plot_dt,
    facet_wrap(~ climate_label, scales = "free_y", ncol = ncol_climate),
    paste0("cf1_channel_ab_concentrations_region_x_climate_pub_", safe_file(rr), ".jpeg"),
    paste0("Lifetime trajectories of CF1 channel concentrations by climate zone: ", rr),
    width = 16,
    height = 9.5,
    by_vars = c("region", "climate_label")
  )

  rm(plot_dt)
  gc()
}

# ------------------------------------------------------------
# 17) WORD-READY TABLE EXPORTS WITH SENSOR-CLUSTERED SEs
# ------------------------------------------------------------

say("Saving Word-ready cluster-robust tables")

make_word_table(
  trend_global,
  "Table S: Global associations between degradation metrics and operational age",
  file.path(outdir, "word_table_global_degradation_time_associations_cluster_robust.docx")
)

make_word_table(
  trend_region,
  "Table S: Regional associations between degradation metrics and operational age",
  file.path(outdir, "word_table_region_degradation_time_associations_cluster_robust.docx")
)

make_word_table(
  trend_climate,
  "Table S: Climate-zone associations between degradation metrics and operational age",
  file.path(outdir, "word_table_climate_degradation_time_associations_cluster_robust.docx")
)

make_word_table(
  trend_region_climate,
  "Table S: Region-by-climate associations between degradation metrics and operational age",
  file.path(outdir, "word_table_region_x_climate_degradation_time_associations_cluster_robust.docx")
)

make_word_table(
  trend_all,
  "Table S: Associations between degradation metrics and operational age across all groups",
  file.path(outdir, "word_table_all_degradation_time_associations_cluster_robust.docx")
)

cat("\nCluster-robust Word-ready tables saved to:\n")
cat(outdir, "\n")

# ------------------------------------------------------------
# 18) 4 x 5 SUMMARY GRID WITH SENSOR-LEVEL CI RIBBONS
#     Adds visible blue GAM smooths like Metric 4
#     n_sensors = purple dot-dash line on right axis
# ------------------------------------------------------------

say("Saving 4 x 5 summary grid with sensor-level CIs, GAMs, and purple n sensors")

setDT(curve_region_plus)

curve_region_plus_plot <- curve_region_plus[
  !is.na(region) &
    is.finite(op_bin) &
    is.finite(n_sensors) &
    n_sensors >= min_sensors_per_bin
]

if (nrow(curve_region_plus_plot) == 0L) {
  stop("No rows available for 4x5 plot.")
}

region_order <- c("Global", sort(setdiff(unique(curve_region_plus_plot$region), "Global")))

metric_specs <- list(
  list(y = "mean_abs", lwr = "lwr_abs", upr = "upr_abs", ylab = "|A − B|", metric_type = "positive"),
  list(y = "mean_bias", lwr = "lwr_bias", upr = "upr_bias", ylab = "A − B", metric_type = "signed"),
  list(y = "mean_pct", lwr = "lwr_pct", upr = "upr_pct", ylab = "Relative divergence", metric_type = "positive"),
  list(y = "pct_flagged", lwr = "lwr_flag", upr = "upr_flag", ylab = "Proportion flagged", metric_type = "proportion")
)

say("Building GAMs for 4 x 5 summary grid")

gam_grid <- rbindlist(
  lapply(metric_specs, function(m) {
    make_gam_predictions(
      d = curve_region_plus_plot,
      yvar = m$y,
      by_vars = "region"
    )
  }),
  use.names = TRUE,
  fill = TRUE
)

make_metric_n_panel <- function(d, gam_d, region_name, y, lwr, upr, ylab, metric_type = "positive") {

  dd <- copy(d[region == region_name])

  dd <- dd[
    is.finite(op_bin) &
      is.finite(get(y)) &
      is.finite(get(lwr)) &
      is.finite(get(upr)) &
      is.finite(n_sensors)
  ]

  if (nrow(dd) == 0L) return(NULL)

  if (metric_type == "signed") {
    lim_vec <- c(dd[[lwr]], dd[[upr]], dd[[y]])
    ymax <- as.numeric(quantile(abs(lim_vec), 0.98, na.rm = TRUE, type = 7))
    if (!is.finite(ymax) || ymax <= 0) ymax <- max(abs(lim_vec), na.rm = TRUE)
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1
    ymin <- -ymax
  } else if (metric_type == "proportion") {
    ymin <- 0
    ymax <- as.numeric(quantile(c(dd[[upr]], dd[[y]]), 0.98, na.rm = TRUE, type = 7)) * 1.08
    ymax <- min(1, ymax)
    if (!is.finite(ymax) || ymax <= 0) ymax <- 0.1
  } else {
    ymin <- 0
    ymax <- as.numeric(quantile(c(dd[[upr]], dd[[y]]), 0.98, na.rm = TRUE, type = 7)) * 1.08
    if (!is.finite(ymax) || ymax <= 0) ymax <- max(c(dd[[upr]], dd[[y]]), na.rm = TRUE)
    if (!is.finite(ymax) || ymax <= 0) ymax <- 1
  }

  dd[, y_plot := pmin(pmax(get(y), ymin), ymax)]
  dd[, lwr_plot := pmin(pmax(get(lwr), ymin), ymax)]
  dd[, upr_plot := pmin(pmax(get(upr), ymin), ymax)]

  gd <- copy(gam_d[region == region_name & variable == y])
  if (nrow(gd) > 0L) {
    gd[, fit_plot := pmin(pmax(fit, ymin), ymax)]
  }

  n_max <- max(dd$n_sensors, na.rm = TRUE)
  if (!is.finite(n_max) || n_max <= 0) n_max <- 1

  dd[, n_scaled := ymin + (n_sensors / n_max) * (ymax - ymin)]

  ggplot(dd, aes(x = op_bin)) +
    geom_ribbon(
      aes(ymin = lwr_plot, ymax = upr_plot),
      fill = "grey45",
      alpha = 0.42,
      color = "grey25",
      linewidth = 0.12
    ) +
    geom_line(aes(y = y_plot), color = "black", linewidth = 0.44, alpha = 0.98) +
    geom_line(
      data = gd,
      aes(x = op_bin, y = fit_plot),
      inherit.aes = FALSE,
      color = "blue4",
      linewidth = 1.05,
      alpha = 0.98
    ) +
    geom_line(
      aes(y = n_scaled),
      color = "purple4",
      linewidth = 0.45,
      alpha = 0.95,
      linetype = "dotdash"
    ) +
    geom_vline(
      xintercept = year_lines,
      linetype = 2,
      linewidth = 0.2,
      alpha = 0.35
    ) +
    scale_x_continuous(labels = comma) +
    scale_y_continuous(
      name = ylab,
      limits = c(ymin, ymax),
      sec.axis = sec_axis(
        trans = ~ (. - ymin) / (ymax - ymin) * n_max,
        name = "n sensors",
        labels = comma
      )
    ) +
    labs(title = region_name, x = "Operational hour") +
    theme_minimal(base_size = 8) +
    theme(
      plot.title = element_text(face = "bold", size = 9),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      axis.title.y.left = element_text(size = 7),
      axis.title.y.right = element_text(size = 7, color = "purple4"),
      axis.text.y.right = element_text(size = 6, color = "purple4"),
      axis.text.y.left = element_text(size = 6),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.2)
    )
}

panel_list <- list()

for (m in metric_specs) {
  for (rr in region_order) {
    p_tmp <- make_metric_n_panel(
      d = curve_region_plus_plot,
      gam_d = gam_grid,
      region_name = rr,
      y = m$y,
      lwr = m$lwr,
      upr = m$upr,
      ylab = m$ylab,
      metric_type = m$metric_type
    )

    if (!is.null(p_tmp)) panel_list[[length(panel_list) + 1L]] <- p_tmp
  }
}

if (length(panel_list) == 0L) {
  stop("No panels were created for 4x5 plot.")
}

p_4x5 <- wrap_plots(panel_list, ncol = length(region_order), byrow = TRUE) +
  plot_annotation(
    title = "CF1 divergence trajectories by region and globally",
    subtitle = "Black = binned sensor-level mean; grey ribbon = sensor-level 95% CI across sensors; blue = GAM smooth; purple dot-dash = n sensors on right y-axis.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10)
    )
  )

ggsave(
  file.path(outdir, "cf1_divergence_4x5_region_global_with_n_sensors_sensor_ci_purple_BLUE_GAM.jpeg"),
  p_4x5,
  width = 20,
  height = 13,
  dpi = 300,
  bg = "white"
)

rm(panel_list, p_4x5, curve_region_plus_plot, gam_grid)
gc()

say("Saved 4 x 5 summary grid with visible blue GAM lines.")
