# ==========================================================================
# 12_divergence_reference.R
#
# Age-dependent divergence vs GEOS-CF PM2.5 and ERA5-Land T/RH.
# => Figure 3; Tables S7-S9; Figures S36-S38.
#
# Extracted VERBATIM from purpleair_degradation2.Rmd (chunk #30 lines 11582-13472).
# Code is unmodified; see REPRODUCING.md for paths and run order.
# ==========================================================================

# ---- source chunk #30 (Rmd lines 11582-13472) ------------------------

# ============================================================
# MEMORY-EFFICIENT METRIC 4
# AGE-DEPENDENT DIVERGENCE + PAIRED REFERENCE TRAJECTORIES
#
# GLOBAL/REGION:
#   Uses outdoor_final.RData
#
# CLIMATE + REGION x CLIMATE:
#   Uses outdoor_with_koppen_timevarying.RData
#   Restricted to valid single-zone Köppen stations only
#
# Metrics:
#   Divergence:
#     PM2.5: PA CF1 average - GEOS-CF PM2.5
#     Temperature: PA temperature - ERA5 temperature
#     RH: PA RH - ERA5 RH
#
#   Paired trajectories:
#     PA PM2.5 and GEOS-CF PM2.5
#     PA temperature and ERA5 temperature
#     PA RH and ERA5 RH
#
# Data filtering:
#   No physical range filtering.
#   Keeps all finite paired observations.
#
# CI:
#   Sensor-level 95% CI across sensors in each operational-age bin.
#   Each sensor contributes one station-bin mean.
#
# GAM:
#   Smooth line only. No GAM CI.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(mgcv)
  library(scales)
})

# ============================================================
# ADD TO PACKAGE LIST AT TOP
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(mgcv)
  library(scales)
  library(sandwich)
  library(lmtest)
  library(officer)
  library(flextable)
})

# ------------------------------------------------------------
# 0) PATHS / SETTINGS
# ------------------------------------------------------------

global_rdata  <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_final.RData"
climate_rdata <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_with_koppen_timevarying.RData"
CHANGE_FILE   <- "/Users/priyanka/Downloads/PurpleAir_EPA/station_koppen_transition_summary.csv"

OUTDIR <- "/Users/priyanka/Downloads/PurpleAir_EPA/metric4_memory_efficient_sensor_ci_no_range_filter"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

BIN_W <- 24L

MIN_SENSORS_PER_BIN_GLOBAL  <- 25L
MIN_SENSORS_PER_BIN_REGION  <- 10L
MIN_SENSORS_PER_BIN_CLIMATE <- 10L
MIN_TOTAL_SENSOR_BINS_REGION_CLIMATE <- 200L
MIN_UNIQUE_BINS_REGION_CLIMATE <- 25L

year_lines <- c(2, 4, 6) * 365.25 * 24

say <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste0(..., collapse = ""), "\n")
}

safe_file <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

# ------------------------------------------------------------
# 1) COLUMN NAMES + LABELS
# ------------------------------------------------------------

a_col <- "pm25_cf_1_a_hourly(ug/m3)_pm25_cf_1_a"
b_col <- "pm25_cf_1_b_hourly(ug/m3)_pm25_cf_1_b"

geos_col      <- "geoscf_pm25"
pa_temp_col   <- "temperature_hourly(C)_temperature"
pa_rh_col     <- "humidity_hourly(%)_humidity"
era5_temp_col <- "era5_t2m_c"
era5_rh_col   <- "era5_rh_pct"
koppen_col    <- "koppen_code"

metric_order <- c(
  "PM2.5: PA_cf1 - GEOS",
  "T: PA - ERA5",
  "RH: PA - ERA5"
)

metric_labels <- list(
  "PM2.5: PA_cf1 - GEOS" = "PM2.5 divergence\n(PA CF1 − GEOS-CF)",
  "T: PA - ERA5"         = "Temperature divergence\n(PA − ERA5, °C)",
  "RH: PA - ERA5"        = "RH divergence\n(PA − ERA5, %)"
)

paired_metric_order <- c("PM2.5", "Temperature", "RH")

paired_metric_labels <- list(
  "PM2.5"       = "PM2.5\n(µg/m³)",
  "Temperature" = "Temperature\n(°C)",
  "RH"          = "Relative humidity\n(%)"
)

# ------------------------------------------------------------
# 2) LOAD HELPERS
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

restrict_single_zone <- function(outdoor) {
  if (!file.exists(CHANGE_FILE)) {
    stop("Missing station Köppen transition file: ", CHANGE_FILE)
  }

  station_change <- fread(CHANGE_FILE)
  station_change[, station := as.character(station)]

  keep_stations <- station_change[
    !is.na(first_koppen) & n_koppen_codes == 1,
    station
  ]

  outdoor[, station := as.character(station)]
  outdoor <- outdoor[station %in% keep_stations]

  say("Single-zone stations retained: ", format(length(keep_stations), big.mark = ","))
  say("Rows after single-zone restriction: ", format(nrow(outdoor), big.mark = ","))
  say("Sensors after single-zone restriction: ", format(uniqueN(outdoor$station), big.mark = ","))

  rm(station_change, keep_stations)
  gc()

  outdoor[]
}

make_koppen_labels <- function() {
  if (!file.exists(CHANGE_FILE)) {
    warning("Köppen transition file missing; climate labels will use codes only.")
    return(data.table(koppen_code = character(), climate_desc = character()))
  }

  out <- fread(CHANGE_FILE)[
    !is.na(first_koppen),
    .(
      koppen_code = as.character(first_koppen),
      climate_desc = as.character(first_koppen_desc)
    )
  ]

  out <- unique(out)
  out[is.na(climate_desc) | climate_desc == "", climate_desc := koppen_code]
  out[]
}

koppen_labels <- make_koppen_labels()

add_climate_labels <- function(dt) {
  if (!("koppen_code" %in% names(dt))) return(dt)

  dt[, koppen_code := as.character(koppen_code)]

  dt <- merge(
    dt,
    koppen_labels,
    by = "koppen_code",
    all.x = TRUE
  )

  dt[is.na(climate_desc) | climate_desc == "", climate_desc := koppen_code]
  dt[, climate_label := paste0(koppen_code, "\n", climate_desc)]
  dt[]
}

# ------------------------------------------------------------
# 3) SENSOR-BIN + CURVE HELPERS
# ------------------------------------------------------------

sensor_ci_curve <- function(sensor_bins, by_vars, min_sensors_keep) {
  out <- sensor_bins[
    is.finite(sensor_mean),
    .(
      n_sensors = uniqueN(station),
      n_sensor_bins = .N,
      n_obs = sum(n_obs, na.rm = TRUE),
      mean = mean(sensor_mean, na.rm = TRUE),
      sd_sensor = sd(sensor_mean, na.rm = TRUE)
    ),
    by = c(by_vars, "metric", "op_bin", "op_years")
  ]

  out[, se := fifelse(n_sensors > 1 & is.finite(sd_sensor), sd_sensor / sqrt(n_sensors), NA_real_)]
  out[, lwr := mean - 1.96 * se]
  out[, upr := mean + 1.96 * se]
  out <- out[n_sensors >= min_sensors_keep]

  setorder(out, metric, op_bin)
  out[]
}

sensor_ci_curve_paired <- function(sensor_bins, by_vars, min_sensors_keep) {
  out <- sensor_bins[
    is.finite(sensor_mean),
    .(
      n_sensors = uniqueN(station),
      n_sensor_bins = .N,
      n_obs = sum(n_obs, na.rm = TRUE),
      mean = mean(sensor_mean, na.rm = TRUE),
      sd_sensor = sd(sensor_mean, na.rm = TRUE)
    ),
    by = c(by_vars, "metric", "source", "op_bin", "op_years")
  ]

  out[, se := fifelse(n_sensors > 1 & is.finite(sd_sensor), sd_sensor / sqrt(n_sensors), NA_real_)]
  out[, lwr := mean - 1.96 * se]
  out[, upr := mean + 1.96 * se]
  out <- out[n_sensors >= min_sensors_keep]

  setorder(out, metric, source, op_bin)
  out[]
}

build_metric_sensor_bins <- function(infile, metric_name, require_koppen = FALSE, single_zone = FALSE) {
  say("Building sensor-bins for: ", metric_name)

  outdoor <- load_outdoor_object(infile)

  if (single_zone) {
    outdoor <- restrict_single_zone(outdoor)
  }

  outdoor[, station := as.character(station)]

  if (metric_name == "PM2.5: PA_cf1 - GEOS") {
    needed <- c("station", "region", "datetime_utc", a_col, b_col, geos_col)
    pair_metric <- "PM2.5"
    pa_source <- "PurpleAir PM2.5"
    ref_source <- "GEOS-CF PM2.5"

    if (require_koppen) needed <- c(needed, koppen_col)
    miss <- setdiff(needed, names(outdoor))
    if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))

    dt <- outdoor[, ..needed]
    rm(outdoor)
    gc()

    dt[, `:=`(
      station = as.character(station),
      region = as.character(region),
      ts = as.POSIXct(datetime_utc, tz = "UTC"),
      pa_value = (as.numeric(get(a_col)) + as.numeric(get(b_col))) / 2,
      ref_value = as.numeric(get(geos_col))
    )]
  }

  if (metric_name == "T: PA - ERA5") {
    needed <- c("station", "region", "datetime_utc", pa_temp_col, era5_temp_col)
    pair_metric <- "Temperature"
    pa_source <- "PurpleAir temperature"
    ref_source <- "ERA5 temperature"

    if (require_koppen) needed <- c(needed, koppen_col)
    miss <- setdiff(needed, names(outdoor))
    if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))

    dt <- outdoor[, ..needed]
    rm(outdoor)
    gc()

    dt[, `:=`(
      station = as.character(station),
      region = as.character(region),
      ts = as.POSIXct(datetime_utc, tz = "UTC"),
      pa_value = as.numeric(get(pa_temp_col)),
      ref_value = as.numeric(get(era5_temp_col))
    )]
  }

  if (metric_name == "RH: PA - ERA5") {
    needed <- c("station", "region", "datetime_utc", pa_rh_col, era5_rh_col)
    pair_metric <- "RH"
    pa_source <- "PurpleAir RH"
    ref_source <- "ERA5 RH"

    if (require_koppen) needed <- c(needed, koppen_col)
    miss <- setdiff(needed, names(outdoor))
    if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))

    dt <- outdoor[, ..needed]
    rm(outdoor)
    gc()

    dt[, `:=`(
      station = as.character(station),
      region = as.character(region),
      ts = as.POSIXct(datetime_utc, tz = "UTC"),
      pa_value = as.numeric(get(pa_rh_col)),
      ref_value = as.numeric(get(era5_rh_col))
    )]
  }

  dt[, datetime_utc := NULL]

  if (require_koppen) {
    dt[, koppen_code := as.character(get(koppen_col))]
    dt <- dt[!is.na(koppen_code) & koppen_code != ""]
  } else {
    dt[, koppen_code := NA_character_]
  }

  # No physical-range filtering. Only finite paired observations.
  dt <- dt[
    !is.na(station) &
      !is.na(region) &
      !is.na(ts) &
      is.finite(pa_value) &
      is.finite(ref_value)
  ]

  dt[, div_value := pa_value - ref_value]
  dt <- dt[is.finite(div_value)]

  setorder(dt, station, ts)
  dt[, first_ts := ts[1L], by = station]
  dt[, op_hours := as.numeric(difftime(ts, first_ts, units = "hours"))]
  dt <- dt[is.finite(op_hours) & op_hours >= 0]
  dt[, op_bin := floor(op_hours / BIN_W) * BIN_W]
  dt[, op_years := op_bin / (365.25 * 24)]

  by_vars <- if (require_koppen) {
    c("station", "region", "koppen_code", "op_bin", "op_years")
  } else {
    c("station", "region", "op_bin", "op_years")
  }

  sens_div <- dt[, .(
    n_obs = .N,
    sensor_mean = mean(div_value, na.rm = TRUE)
  ), by = by_vars]

  sens_div[, metric := metric_name]

  sens_pa <- dt[, .(
    n_obs = .N,
    sensor_mean = mean(pa_value, na.rm = TRUE)
  ), by = by_vars]

  sens_pa[, `:=`(
    metric = pair_metric,
    source = pa_source
  )]

  sens_ref <- dt[, .(
    n_obs = .N,
    sensor_mean = mean(ref_value, na.rm = TRUE)
  ), by = by_vars]

  sens_ref[, `:=`(
    metric = pair_metric,
    source = ref_source
  )]

  sens_pair <- rbindlist(list(sens_pa, sens_ref), use.names = TRUE, fill = TRUE)

  rm(dt, sens_pa, sens_ref)
  gc()

  list(div = sens_div, pair = sens_pair)
}

# ------------------------------------------------------------
# 3B) CLUSTER-ROBUST TREND HELPERS + WORD TABLE HELPERS
# ------------------------------------------------------------

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

  p_val <- if (length(p_col) == 0L) {
    NA_real_
  } else {
    as.numeric(ct[term, p_col[1]])
  }

  list(
    estimate = as.numeric(ct[term, "Estimate"]),
    se = as.numeric(ct[term, "Std. Error"]),
    p_value = p_val
  )
}

fit_cluster_trend <- function(dt, group_vars, table_type = "divergence") {
  dd <- copy(dt)

  dd <- dd[
    is.finite(sensor_mean) &
      is.finite(op_years) &
      is.finite(n_obs) &
      n_obs > 0 &
      !is.na(station)
  ]

  if (nrow(dd) == 0L) return(data.table())

  dd[, {
    dsub <- .SD[
      is.finite(sensor_mean) &
        is.finite(op_years) &
        is.finite(n_obs) &
        n_obs > 0 &
        !is.na(station)
    ]

    if (nrow(dsub) < 20L || uniqueN(dsub$op_years) < 5L || uniqueN(dsub$station) < 10L) {
      data.table()
    } else {
      fit <- tryCatch(
        lm(sensor_mean ~ op_years, data = dsub, weights = n_obs),
        error = function(e) NULL
      )

      if (is.null(fit)) {
        data.table()
      } else {
        ct <- tryCatch(
          {
            vc <- sandwich::vcovCL(
              fit,
              cluster = dsub$station,
              type = "HC1"
            )
            lmtest::coeftest(fit, vcov. = vc)
          },
          error = function(e) NULL
        )

        ext <- extract_term(ct, "op_years")

        if (is.null(ext)) {
          data.table()
        } else {
          beta <- ext$estimate
          se <- ext$se

          data.table(
            table_type = table_type,
            n_sensor_bins = nrow(dsub),
            n_sensors = uniqueN(dsub$station),
            min_year = min(dsub$op_years, na.rm = TRUE),
            max_year = max(dsub$op_years, na.rm = TRUE),
            estimate = beta,
            ci_lwr = beta - 1.96 * se,
            ci_upr = beta + 1.96 * se,
            se_cluster = se,
            p_value = ext$p_value
          )
        }
      }
    }
  }, by = group_vars]
}

format_trend_results <- function(dt) {
  if (nrow(dt) == 0L) return(dt)

  dt[, `:=`(
    estimate_95ci = fmt_ci(estimate, ci_lwr, ci_upr, digits = 3),
    p_value_fmt = fmt_p(p_value),
    estimate_95ci_p = paste0(
      fmt_ci(estimate, ci_lwr, ci_upr, digits = 3),
      "; p=",
      fmt_p(p_value)
    )
  )]

  dt[]
}

make_metric4_word_table <- function(dt, title, outfile) {
  if (nrow(dt) == 0L) return(invisible(NULL))

  tab <- copy(dt)

  if (!("source" %in% names(tab))) tab[, source := NA_character_]
  if (!("region" %in% names(tab))) tab[, region := NA_character_]
  if (!("group" %in% names(tab))) tab[, group := NA_character_]

  tab[, `:=`(
    `Estimate (95% CI)` = fmt_ci(estimate, ci_lwr, ci_upr, digits = 3),
    `P value` = fmt_p(p_value),
    `Operational years` = paste0(round(min_year, 2), "-", round(max_year, 2))
  )]

  tab_out <- tab[
    , .(
      Group = group,
      Region = region,
      Metric = as.character(metric),
      Source = source,
      `Effect reported` = "Beta per operational year",
      `Estimate (95% CI)`,
      `P value`,
      `N sensors` = n_sensors,
      `N sensor-bins` = n_sensor_bins,
      `Operational years`
    )
  ]

  tab_out <- tab_out[
    ,
    lapply(.SD, function(x) {
      if (is.character(x)) fifelse(is.na(x) | x == "", "-", x) else x
    })
  ]

  ft <- flextable(tab_out)
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
      "Metric 4 trend estimates are reported as beta coefficients per operational year. ",
      "Models were fit to sensor-by-operational-age-bin means and weighted by the number of contributing hourly observations. ",
      "Confidence intervals are 95% intervals calculated using sensor-clustered robust standard errors."
    ),
    style = "Normal"
  )
  doc <- body_add_flextable(doc, ft)

  print(doc, target = outfile)

  invisible(NULL)
}

# ------------------------------------------------------------
# 4) BUILD CURVES METRIC-BY-METRIC
# ------------------------------------------------------------

curves_global <- list()
curves_region <- list()
curves_climate <- list()
curves_region_climate <- list()

paired_global <- list()
paired_region <- list()
paired_climate <- list()
paired_region_climate <- list()

# ------------------------------------------------------------
# 4A) TREND RESULT LISTS
# ------------------------------------------------------------

trend_m4_global_list <- list()
trend_m4_region_list <- list()
trend_m4_climate_list <- list()
trend_m4_region_climate_list <- list()

trend_pair_m4_global_list <- list()
trend_pair_m4_region_list <- list()
trend_pair_m4_climate_list <- list()
trend_pair_m4_region_climate_list <- list()

for (met in metric_order) {

  say("Processing global/region metric: ", met)

  sens_gr <- build_metric_sensor_bins(
    infile = global_rdata,
    metric_name = met,
    require_koppen = FALSE,
    single_zone = FALSE
  )

    # Cluster-robust Metric 4 trend tests: global and region
  tg <- fit_cluster_trend(
    sens_gr$div,
    group_vars = "metric",
    table_type = "divergence"
  )
  tg[, group := "Global"]

  tr <- fit_cluster_trend(
    sens_gr$div,
    group_vars = c("region", "metric"),
    table_type = "divergence"
  )
  tr[, group := region]

  tpg <- fit_cluster_trend(
    sens_gr$pair,
    group_vars = c("metric", "source"),
    table_type = "paired_value"
  )
  tpg[, group := "Global"]

  tpr <- fit_cluster_trend(
    sens_gr$pair,
    group_vars = c("region", "metric", "source"),
    table_type = "paired_value"
  )
  tpr[, group := region]

  trend_m4_global_list[[met]] <- tg
  trend_m4_region_list[[met]] <- tr
  trend_pair_m4_global_list[[met]] <- tpg
  trend_pair_m4_region_list[[met]] <- tpr

  rm(tg, tr, tpg, tpr)
  gc()
  
  curve_g <- sensor_ci_curve(
    sens_gr$div,
    by_vars = character(),
    min_sensors_keep = MIN_SENSORS_PER_BIN_GLOBAL
  )
  curve_g[, group := "Global"]

  curve_r <- sensor_ci_curve(
    sens_gr$div,
    by_vars = "region",
    min_sensors_keep = MIN_SENSORS_PER_BIN_REGION
  )
  curve_r[, group := region]

  pair_g <- sensor_ci_curve_paired(
    sens_gr$pair,
    by_vars = character(),
    min_sensors_keep = MIN_SENSORS_PER_BIN_GLOBAL
  )
  pair_g[, group := "Global"]

  pair_r <- sensor_ci_curve_paired(
    sens_gr$pair,
    by_vars = "region",
    min_sensors_keep = MIN_SENSORS_PER_BIN_REGION
  )
  pair_r[, group := region]

  curves_global[[met]] <- curve_g
  curves_region[[met]] <- curve_r
  paired_global[[met]] <- pair_g
  paired_region[[met]] <- pair_r

  rm(sens_gr, curve_g, curve_r, pair_g, pair_r)
  gc()

  say("Processing climate and region x climate metric: ", met)

  sens_clim <- build_metric_sensor_bins(
    infile = climate_rdata,
    metric_name = met,
    require_koppen = TRUE,
    single_zone = TRUE
  )

    # Cluster-robust Metric 4 trend tests: climate and region x climate
  tc <- fit_cluster_trend(
    sens_clim$div,
    group_vars = c("koppen_code", "metric"),
    table_type = "divergence"
  )
  tc <- add_climate_labels(tc)
  tc[, group := climate_label]

  trc <- fit_cluster_trend(
    sens_clim$div,
    group_vars = c("region", "koppen_code", "metric"),
    table_type = "divergence"
  )
  trc <- add_climate_labels(trc)
  trc[, group := climate_label]

  tpc <- fit_cluster_trend(
    sens_clim$pair,
    group_vars = c("koppen_code", "metric", "source"),
    table_type = "paired_value"
  )
  tpc <- add_climate_labels(tpc)
  tpc[, group := climate_label]

  tprc <- fit_cluster_trend(
    sens_clim$pair,
    group_vars = c("region", "koppen_code", "metric", "source"),
    table_type = "paired_value"
  )
  tprc <- add_climate_labels(tprc)
  tprc[, group := climate_label]

  trend_m4_climate_list[[met]] <- tc
  trend_m4_region_climate_list[[met]] <- trc
  trend_pair_m4_climate_list[[met]] <- tpc
  trend_pair_m4_region_climate_list[[met]] <- tprc

  rm(tc, trc, tpc, tprc)
  gc()
  
  curve_c <- sensor_ci_curve(
    sens_clim$div,
    by_vars = "koppen_code",
    min_sensors_keep = MIN_SENSORS_PER_BIN_CLIMATE
  )
  curve_c <- add_climate_labels(curve_c)
  curve_c[, group := climate_label]

  curve_rc <- sensor_ci_curve(
    sens_clim$div,
    by_vars = c("region", "koppen_code"),
    min_sensors_keep = MIN_SENSORS_PER_BIN_CLIMATE
  )
  curve_rc <- add_climate_labels(curve_rc)
  curve_rc[, group := climate_label]

  pair_c <- sensor_ci_curve_paired(
    sens_clim$pair,
    by_vars = "koppen_code",
    min_sensors_keep = MIN_SENSORS_PER_BIN_CLIMATE
  )
  pair_c <- add_climate_labels(pair_c)
  pair_c[, group := climate_label]

  pair_rc <- sensor_ci_curve_paired(
    sens_clim$pair,
    by_vars = c("region", "koppen_code"),
    min_sensors_keep = MIN_SENSORS_PER_BIN_CLIMATE
  )
  pair_rc <- add_climate_labels(pair_rc)
  pair_rc[, group := climate_label]

  curves_climate[[met]] <- curve_c
  curves_region_climate[[met]] <- curve_rc
  paired_climate[[met]] <- pair_c
  paired_region_climate[[met]] <- pair_rc

  rm(sens_clim, curve_c, curve_rc, pair_c, pair_rc)
  gc()
}

curve_m4_global <- rbindlist(curves_global, use.names = TRUE, fill = TRUE)
curve_m4_region <- rbindlist(curves_region, use.names = TRUE, fill = TRUE)
curve_m4_climate <- rbindlist(curves_climate, use.names = TRUE, fill = TRUE)
curve_m4_region_climate <- rbindlist(curves_region_climate, use.names = TRUE, fill = TRUE)

pair_m4_global <- rbindlist(paired_global, use.names = TRUE, fill = TRUE)
pair_m4_region <- rbindlist(paired_region, use.names = TRUE, fill = TRUE)
pair_m4_climate <- rbindlist(paired_climate, use.names = TRUE, fill = TRUE)
pair_m4_region_climate <- rbindlist(paired_region_climate, use.names = TRUE, fill = TRUE)

curve_m4_global[, metric := factor(metric, levels = metric_order)]
curve_m4_region[, metric := factor(metric, levels = metric_order)]
curve_m4_climate[, metric := factor(metric, levels = metric_order)]
curve_m4_region_climate[, metric := factor(metric, levels = metric_order)]

pair_m4_global[, metric := factor(metric, levels = paired_metric_order)]
pair_m4_region[, metric := factor(metric, levels = paired_metric_order)]
pair_m4_climate[, metric := factor(metric, levels = paired_metric_order)]
pair_m4_region_climate[, metric := factor(metric, levels = paired_metric_order)]

curve_m4_global_region <- rbindlist(
  list(
    curve_m4_global[, .(group, metric, op_bin, op_years, n_sensors, n_sensor_bins, n_obs, mean, lwr, upr)],
    curve_m4_region[, .(group, metric, op_bin, op_years, n_sensors, n_sensor_bins, n_obs, mean, lwr, upr)]
  ),
  use.names = TRUE,
  fill = TRUE
)

pair_m4_global_region <- rbindlist(
  list(
    pair_m4_global[, .(group, metric, source, op_bin, op_years, n_sensors, n_sensor_bins, n_obs, mean, lwr, upr)],
    pair_m4_region[, .(group, metric, source, op_bin, op_years, n_sensors, n_sensor_bins, n_obs, mean, lwr, upr)]
  ),
  use.names = TRUE,
  fill = TRUE
)

valid_region_climate_panels <- curve_m4_region_climate[
  ,
  .(
    total_sensor_bins = sum(n_sensor_bins, na.rm = TRUE),
    n_bins = uniqueN(op_bin)
  ),
  by = .(region, group)
][
  total_sensor_bins >= MIN_TOTAL_SENSOR_BINS_REGION_CLIMATE &
    n_bins >= MIN_UNIQUE_BINS_REGION_CLIMATE
]

curve_m4_region_climate_plot <- merge(
  curve_m4_region_climate,
  valid_region_climate_panels[, .(region, group)],
  by = c("region", "group"),
  all = FALSE
)

pair_m4_region_climate_plot <- merge(
  pair_m4_region_climate,
  valid_region_climate_panels[, .(region, group)],
  by = c("region", "group"),
  all = FALSE
)

# ------------------------------------------------------------
# 5) WRITE CURVES
# ------------------------------------------------------------

fwrite(curve_m4_global, file.path(OUTDIR, "curve_m4_global_divergence_sensor_ci.csv"))
fwrite(curve_m4_region, file.path(OUTDIR, "curve_m4_region_divergence_sensor_ci.csv"))
fwrite(curve_m4_global_region, file.path(OUTDIR, "curve_m4_global_region_divergence_sensor_ci.csv"))
fwrite(curve_m4_climate, file.path(OUTDIR, "curve_m4_climate_divergence_sensor_ci.csv"))
fwrite(curve_m4_region_climate, file.path(OUTDIR, "curve_m4_region_x_climate_divergence_sensor_ci.csv"))
fwrite(curve_m4_region_climate_plot, file.path(OUTDIR, "curve_m4_region_x_climate_PLOTTED_divergence_sensor_ci.csv"))

fwrite(pair_m4_global, file.path(OUTDIR, "curve_m4_global_paired_values_sensor_ci.csv"))
fwrite(pair_m4_region, file.path(OUTDIR, "curve_m4_region_paired_values_sensor_ci.csv"))
fwrite(pair_m4_global_region, file.path(OUTDIR, "curve_m4_global_region_paired_values_sensor_ci.csv"))
fwrite(pair_m4_climate, file.path(OUTDIR, "curve_m4_climate_paired_values_sensor_ci.csv"))
fwrite(pair_m4_region_climate, file.path(OUTDIR, "curve_m4_region_x_climate_paired_values_sensor_ci.csv"))
fwrite(pair_m4_region_climate_plot, file.path(OUTDIR, "curve_m4_region_x_climate_PLOTTED_paired_values_sensor_ci.csv"))

# ------------------------------------------------------------
# 6) GAM SMOOTH
#     GAM is trend line only.
#     For region x climate, suppress GAM only when late-life support is sparse.
#     This does NOT suppress climate-zone panels.
# ------------------------------------------------------------

# ------------------------------------------------------------
# 6) GAM SMOOTH
#     GAM is trend line only.
#     For region x climate, suppress GAM when late-life support is sparse.
#     Also truncate GAM predictions to well-supported bins.
#     This does NOT suppress climate-zone panels.
# ------------------------------------------------------------

MIN_GAM_BINS <- 30L
MIN_GAM_TOTAL_SENSOR_BINS <- 250L
MIN_GAM_LATE_BINS <- 8L
MIN_GAM_LATE_MEDIAN_SENSORS <- 10L
MIN_GAM_SENSORS_PER_BIN <- 10L
MAX_GAM_GAP_HOURS <- 24L * 90L

get_supported_gam_range <- function(dd) {
  dd <- copy(dd)[
    is.finite(mean) &
      is.finite(op_bin) &
      is.finite(n_sensors) &
      n_sensors >= MIN_GAM_SENSORS_PER_BIN
  ]

  if (nrow(dd) == 0L) return(NULL)

  setorder(dd, op_bin)

  dd[, gap_prev := op_bin - shift(op_bin)]
  dd[, segment_id := cumsum(is.na(gap_prev) | gap_prev > MAX_GAM_GAP_HOURS)]

  seg <- dd[
    ,
    .(
      n_bins = uniqueN(op_bin),
      total_sensors = sum(n_sensors, na.rm = TRUE),
      min_bin = min(op_bin, na.rm = TRUE),
      max_bin = max(op_bin, na.rm = TRUE)
    ),
    by = segment_id
  ][
    n_bins >= MIN_GAM_LATE_BINS &
      total_sensors >= MIN_GAM_TOTAL_SENSOR_BINS * 0.10
  ]

  if (nrow(seg) == 0L) {
    return(list(
      min_bin = min(dd$op_bin, na.rm = TRUE),
      max_bin = max(dd$op_bin, na.rm = TRUE)
    ))
  }

  list(
    min_bin = min(seg$min_bin, na.rm = TRUE),
    max_bin = max(seg$max_bin, na.rm = TRUE)
  )
}

gam_quality_ok <- function(dd) {
  dd <- copy(dd)[
    is.finite(mean) &
      is.finite(op_bin) &
      is.finite(n_sensors) &
      n_sensors > 0
  ]

  if (nrow(dd) < MIN_GAM_BINS) return(FALSE)
  if (uniqueN(dd$op_bin) < MIN_GAM_BINS) return(FALSE)
  if (sum(dd$n_sensors, na.rm = TRUE) < MIN_GAM_TOTAL_SENSOR_BINS) return(FALSE)

  support <- get_supported_gam_range(dd)
  if (is.null(support)) return(FALSE)

  dd <- dd[op_bin >= support$min_bin & op_bin <= support$max_bin]

  mid_bin <- median(range(dd$op_bin, na.rm = TRUE), na.rm = TRUE)
  late <- dd[op_bin >= mid_bin & n_sensors >= MIN_GAM_SENSORS_PER_BIN]

  if (nrow(late) < MIN_GAM_LATE_BINS) return(FALSE)
  if (uniqueN(late$op_bin) < MIN_GAM_LATE_BINS) return(FALSE)
  if (median(late$n_sensors, na.rm = TRUE) < MIN_GAM_LATE_MEDIAN_SENSORS) return(FALSE)

  TRUE
}

fit_gam_pred_curve <- function(dt, by_step = 168, sparse_guard = FALSE) {
  dd <- copy(dt)[
    is.finite(mean) &
      is.finite(op_bin) &
      is.finite(n_sensors) &
      n_sensors > 0
  ]

  if (sparse_guard && !gam_quality_ok(dd)) return(NULL)

  support <- get_supported_gam_range(dd)
  if (is.null(support)) return(NULL)

  dd_fit <- dd[op_bin >= support$min_bin & op_bin <= support$max_bin]

  if (nrow(dd_fit) < 80L || uniqueN(dd_fit$op_bin) < 30L) return(NULL)

  g <- tryCatch(
    gam(
      mean ~ s(op_bin, bs = "ps"),
      data = dd_fit,
      weights = n_sensors,
      method = "REML"
    ),
    error = function(e) NULL
  )

  if (is.null(g)) return(NULL)

  xg <- data.table(
    op_bin = seq(
      min(dd_fit$op_bin, na.rm = TRUE),
      max(dd_fit$op_bin, na.rm = TRUE),
      by = by_step
    )
  )

  xg[, fit := as.numeric(predict(g, newdata = xg, se.fit = FALSE))]
  xg[]
}

make_gam_by_group_metric <- function(curve_dt, group_vars = c("group", "metric"), sparse_guard = FALSE) {
  out <- curve_dt[, {
    gp <- fit_gam_pred_curve(.SD, sparse_guard = sparse_guard)
    if (is.null(gp)) NULL else gp
  }, by = group_vars]

  if (nrow(out) > 0) out[] else data.table()
}

make_gam_by_group_metric_source <- function(curve_dt, sparse_guard = FALSE) {
  out <- curve_dt[, {
    gp <- fit_gam_pred_curve(.SD, sparse_guard = sparse_guard)
    if (is.null(gp)) NULL else gp
  }, by = .(group, metric, source)]

  if (nrow(out) > 0) out[] else data.table()
}

say("Building GAM smooths")

gam_global_region <- make_gam_by_group_metric(curve_m4_global_region)
gam_climate <- make_gam_by_group_metric(curve_m4_climate)

gam_pair_global_region <- make_gam_by_group_metric_source(pair_m4_global_region)
gam_pair_climate <- make_gam_by_group_metric_source(pair_m4_climate)

gam_region_climate <- make_gam_by_group_metric(
  curve_m4_region_climate_plot,
  sparse_guard = TRUE
)

gam_pair_region_climate <- make_gam_by_group_metric_source(
  pair_m4_region_climate_plot,
  sparse_guard = TRUE
)

# ------------------------------------------------------------
# 6B) CLUSTER-ROBUST STATISTICAL TESTING FOR LINEAR TIME TRENDS
#     Uses sensor-bin-level data.
#     Model: sensor_mean ~ op_years, weighted by n_obs.
#     SEs clustered at sensor level.
# ------------------------------------------------------------

say("Combining cluster-robust Metric 4 trend-test tables")

trend_m4_global_region <- rbindlist(
  list(
    rbindlist(trend_m4_global_list, use.names = TRUE, fill = TRUE),
    rbindlist(trend_m4_region_list, use.names = TRUE, fill = TRUE)
  ),
  use.names = TRUE,
  fill = TRUE
)

trend_m4_climate <- rbindlist(
  trend_m4_climate_list,
  use.names = TRUE,
  fill = TRUE
)

trend_m4_region_climate <- rbindlist(
  trend_m4_region_climate_list,
  use.names = TRUE,
  fill = TRUE
)

trend_pair_m4_global_region <- rbindlist(
  list(
    rbindlist(trend_pair_m4_global_list, use.names = TRUE, fill = TRUE),
    rbindlist(trend_pair_m4_region_list, use.names = TRUE, fill = TRUE)
  ),
  use.names = TRUE,
  fill = TRUE
)

trend_pair_m4_climate <- rbindlist(
  trend_pair_m4_climate_list,
  use.names = TRUE,
  fill = TRUE
)

trend_pair_m4_region_climate <- rbindlist(
  trend_pair_m4_region_climate_list,
  use.names = TRUE,
  fill = TRUE
)

trend_m4_global_region <- format_trend_results(trend_m4_global_region)
trend_m4_climate <- format_trend_results(trend_m4_climate)
trend_m4_region_climate <- format_trend_results(trend_m4_region_climate)

trend_pair_m4_global_region <- format_trend_results(trend_pair_m4_global_region)
trend_pair_m4_climate <- format_trend_results(trend_pair_m4_climate)
trend_pair_m4_region_climate <- format_trend_results(trend_pair_m4_region_climate)

trend_m4_all <- rbindlist(
  list(
    trend_m4_global_region,
    trend_m4_climate,
    trend_m4_region_climate
  ),
  use.names = TRUE,
  fill = TRUE
)

trend_pair_m4_all <- rbindlist(
  list(
    trend_pair_m4_global_region,
    trend_pair_m4_climate,
    trend_pair_m4_region_climate
  ),
  use.names = TRUE,
  fill = TRUE
)

# Save CSVs
fwrite(
  trend_m4_global_region,
  file.path(OUTDIR, "trend_tests_metric4_divergence_global_region_cluster_robust.csv")
)

fwrite(
  trend_m4_climate,
  file.path(OUTDIR, "trend_tests_metric4_divergence_climate_cluster_robust.csv")
)

fwrite(
  trend_m4_region_climate,
  file.path(OUTDIR, "trend_tests_metric4_divergence_region_x_climate_cluster_robust.csv")
)

fwrite(
  trend_m4_all,
  file.path(OUTDIR, "trend_tests_metric4_divergence_all_groups_cluster_robust.csv")
)

fwrite(
  trend_pair_m4_global_region,
  file.path(OUTDIR, "trend_tests_metric4_paired_values_global_region_cluster_robust.csv")
)

fwrite(
  trend_pair_m4_climate,
  file.path(OUTDIR, "trend_tests_metric4_paired_values_climate_cluster_robust.csv")
)

fwrite(
  trend_pair_m4_region_climate,
  file.path(OUTDIR, "trend_tests_metric4_paired_values_region_x_climate_cluster_robust.csv")
)

fwrite(
  trend_pair_m4_all,
  file.path(OUTDIR, "trend_tests_metric4_paired_values_all_groups_cluster_robust.csv")
)

# Save Word-ready tables
say("Saving Metric 4 Word-ready cluster-robust tables")

make_metric4_word_table(
  trend_m4_global_region,
  "Table S: Metric 4 divergence trends by global and regional groups",
  file.path(OUTDIR, "word_table_metric4_divergence_global_region_cluster_robust.docx")
)

make_metric4_word_table(
  trend_m4_climate,
  "Table S: Metric 4 divergence trends by climate zone",
  file.path(OUTDIR, "word_table_metric4_divergence_climate_cluster_robust.docx")
)

make_metric4_word_table(
  trend_m4_region_climate,
  "Table S: Metric 4 divergence trends by region and climate zone",
  file.path(OUTDIR, "word_table_metric4_divergence_region_x_climate_cluster_robust.docx")
)

make_metric4_word_table(
  trend_m4_all,
  "Table S: Metric 4 divergence trends across all groups",
  file.path(OUTDIR, "word_table_metric4_divergence_all_groups_cluster_robust.docx")
)

make_metric4_word_table(
  trend_pair_m4_global_region,
  "Table S: Metric 4 paired PurpleAir and reference trends by global and regional groups",
  file.path(OUTDIR, "word_table_metric4_paired_values_global_region_cluster_robust.docx")
)

make_metric4_word_table(
  trend_pair_m4_climate,
  "Table S: Metric 4 paired PurpleAir and reference trends by climate zone",
  file.path(OUTDIR, "word_table_metric4_paired_values_climate_cluster_robust.docx")
)

make_metric4_word_table(
  trend_pair_m4_region_climate,
  "Table S: Metric 4 paired PurpleAir and reference trends by region and climate zone",
  file.path(OUTDIR, "word_table_metric4_paired_values_region_x_climate_cluster_robust.docx")
)

make_metric4_word_table(
  trend_pair_m4_all,
  "Table S: Metric 4 paired PurpleAir and reference trends across all groups",
  file.path(OUTDIR, "word_table_metric4_paired_values_all_groups_cluster_robust.docx")
)

cat("\nMetric 4 cluster-robust trend tables saved to:\n")
cat(OUTDIR, "\n")

cat("\nMetric 4 divergence trend tests: global + region\n")
print(trend_m4_global_region)

cat("\nMetric 4 paired-value trend tests: global + region\n")
print(trend_pair_m4_global_region)

# ------------------------------------------------------------
# 7) PLOT HELPERS — WRAPPED CLIMATE TITLES
#     n_sensors shown as purple dot-dash line/right axis
# ------------------------------------------------------------
# ------------------------------------------------------------
# 7) PLOT HELPERS — WRAPPED PANEL GRID
#     Fixes panel labels beyond Z and makes climate grids taller
# ------------------------------------------------------------

make_panel_labels <- function(n) {
  paste0(
    vapply(seq_len(n), function(i) {
      if (i <= 26) {
        LETTERS[i]
      } else {
        paste0(
          LETTERS[((i - 1) %/% 26)],
          LETTERS[((i - 1) %% 26) + 1]
        )
      }
    }, character(1)),
    ")"
  )
}

wrap_group_title <- function(x, width = 18) {
  vapply(
    as.character(x),
    function(z) paste(strwrap(z, width = width), collapse = "\n"),
    character(1)
  )
}

auto_panel_ncol <- function(n_groups) {
  if (n_groups <= 6) return(n_groups)
  if (n_groups <= 12) return(4L)
  return(5L)
}

auto_fig_width <- function(n_groups, panel_ncol = auto_panel_ncol(n_groups)) {
  max(10, min(22, 3.2 * panel_ncol))
}

auto_fig_height <- function(n_metric_rows, n_groups, panel_ncol = auto_panel_ncol(n_groups)) {
  panel_nrow <- ceiling(n_groups / panel_ncol)
  max(9, 2.6 * n_metric_rows * panel_nrow + 1.4)
}

make_metric_strip <- function(label, size = 10) {
  ggdraw() +
    draw_label(label, angle = 90, fontface = "bold", size = size, x = 0.5, y = 0.5)
}

add_year_lines <- function(p) {
  p +
    geom_vline(
      xintercept = year_lines,
      linetype = 2,
      alpha = 0.35,
      linewidth = 0.45
    )
}

add_nsensors_layer <- function(p, dt,
                               ypad = 0.05,
                               alpha = 0.75,
                               linewidth = 0.75,
                               color = "purple4",
                               show_right_axis = TRUE) {
  y_min <- suppressWarnings(min(c(dt$lwr, dt$mean), na.rm = TRUE))
  y_max <- suppressWarnings(max(c(dt$upr, dt$mean), na.rm = TRUE))

  if (!is.finite(y_min) || !is.finite(y_max) || y_min == y_max) return(p)

  nmax <- suppressWarnings(max(dt$n_sensors, na.rm = TRUE))
  if (!is.finite(nmax) || nmax <= 0) return(p)

  yr <- y_max - y_min
  y0 <- y_min + ypad * yr
  y1 <- y_max - ypad * yr

  dt_ns <- copy(dt)
  dt_ns[, y_ns := y0 + (n_sensors / nmax) * (y1 - y0)]

  inv_fun <- function(yval) (yval - y0) / (y1 - y0) * nmax

  p <- p +
    geom_line(
      data = dt_ns,
      aes(x = op_bin, y = y_ns),
      inherit.aes = FALSE,
      alpha = alpha,
      linewidth = linewidth,
      color = color,
      linetype = "dotdash"
    )

  if (show_right_axis) {
    p <- p +
      scale_y_continuous(
        sec.axis = sec_axis(~ inv_fun(.), name = "n sensors")
      )
  } else {
    p <- p + scale_y_continuous()
  }

  p
}

th_grid <- theme_minimal(base_size = 9) +
  theme(
    plot.title = element_text(
      size = 8.4,
      face = "bold",
      hjust = 0.5,
      lineheight = 0.9,
      margin = margin(b = 8)
    ),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8, color = "black"),
    axis.title.y.right = element_text(size = 8, color = "purple4"),
    axis.text.y.right = element_text(size = 7, color = "purple4"),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey92", linewidth = 0.25),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.35),
    plot.margin = margin(10, 4, 6, 4)
  )

make_cell_m4 <- function(dt_curve, dt_gam,
                         ylab_left = NULL,
                         title = NULL,
                         panel_label = NULL,
                         show_x = TRUE,
                         show_left_y = TRUE,
                         show_right_y = TRUE,
                         gam_color = "blue4") {
  dd <- copy(dt_curve)[order(op_bin)]
  if (nrow(dd) == 0L) return(ggplot() + theme_void())

  p <- ggplot(dd, aes(x = op_bin, y = mean)) +
    geom_ribbon(
      aes(ymin = lwr, ymax = upr),
      fill = "grey45",
      alpha = 0.45,
      color = "grey25",
      linewidth = 0.12
    ) +
    geom_line(linewidth = 0.65, color = "black") +
    labs(title = title, x = NULL, y = ylab_left) +
    th_grid

  if (!is.null(dt_gam) && nrow(dt_gam) > 0L) {
    p <- p +
      geom_line(
        data = dt_gam,
        aes(x = op_bin, y = fit),
        inherit.aes = FALSE,
        linewidth = 1.15,
        color = gam_color
      )
  }

  p <- add_year_lines(p)

  p <- add_nsensors_layer(
    p,
    dd,
    color = "purple4",
    alpha = 0.75,
    linewidth = 0.75,
    show_right_axis = show_right_y
  )

  if (!show_x) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  } else {
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  }

  if (!show_left_y) {
    p <- p + theme(
      axis.title.y.left = element_blank(),
      axis.text.y.left = element_blank(),
      axis.ticks.y.left = element_blank()
    )
  }

  if (!show_right_y) {
    p <- p + theme(
      axis.title.y.right = element_blank(),
      axis.text.y.right = element_blank(),
      axis.ticks.y.right = element_blank()
    )
  }

  if (!is.null(panel_label)) {
    p <- p +
      annotate(
        "text",
        x = -Inf, y = Inf,
        label = panel_label,
        hjust = -0.12,
        vjust = 1.22,
        fontface = "bold",
        size = 3.3
      )
  }

  p
}

make_cell_pair <- function(dt_curve, dt_gam,
                           ylab_left = NULL,
                           title = NULL,
                           panel_label = NULL,
                           show_x = TRUE,
                           show_left_y = TRUE,
                           show_right_y = TRUE) {
  dd <- copy(dt_curve)[order(source, op_bin)]
  if (nrow(dd) == 0L) return(ggplot() + theme_void())

  p <- ggplot(dd, aes(x = op_bin, y = mean, color = source, fill = source)) +
    geom_ribbon(
      aes(ymin = lwr, ymax = upr),
      alpha = 0.25,
      color = NA
    ) +
    geom_line(linewidth = 0.55, alpha = 0.92) +
    labs(title = title, x = NULL, y = ylab_left, color = NULL, fill = NULL) +
    th_grid +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 6),
      legend.key.size = unit(0.35, "lines")
    )

  if (!is.null(dt_gam) && nrow(dt_gam) > 0L) {
    p <- p +
      geom_line(
        data = dt_gam,
        aes(x = op_bin, y = fit, color = source),
        inherit.aes = FALSE,
        linewidth = 1.05
      )
  }

  p <- add_year_lines(p)

  p <- add_nsensors_layer(
    p,
    dd,
    color = "purple4",
    alpha = 0.75,
    linewidth = 0.75,
    show_right_axis = show_right_y
  )

  if (!show_x) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  } else {
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  }

  if (!show_left_y) {
    p <- p + theme(
      axis.title.y.left = element_blank(),
      axis.text.y.left = element_blank(),
      axis.ticks.y.left = element_blank()
    )
  }

  if (!show_right_y) {
    p <- p + theme(
      axis.title.y.right = element_blank(),
      axis.text.y.right = element_blank(),
      axis.ticks.y.right = element_blank()
    )
  }

  if (!is.null(panel_label)) {
    p <- p +
      annotate(
        "text",
        x = -Inf,
        y = Inf,
        label = panel_label,
        hjust = -0.12,
        vjust = 1.22,
        fontface = "bold",
        size = 3.3
      )
  }

  p
}

make_metric_grid <- function(curve_dt, gam_dt, group_order, outfile, title_text, note_text,
                             width = NULL, height = NULL, panel_ncol = NULL) {
  curve_dt <- copy(curve_dt)
  gam_dt <- copy(gam_dt)

  n_groups <- length(group_order)
  if (is.null(panel_ncol)) panel_ncol <- auto_panel_ncol(n_groups)

  if (is.null(width)) width <- auto_fig_width(n_groups, panel_ncol)
  if (is.null(height)) height <- auto_fig_height(length(metric_order), n_groups, panel_ncol)

  title_wrap_width <- ifelse(panel_ncol <= 4, 20, 16)
  title_size <- ifelse(panel_ncol <= 4, 8.2, 7.2)
  strip_size <- 9.5

  curve_dt[, group := factor(group, levels = group_order)]
  curve_dt[, metric := factor(metric, levels = metric_order)]

  if (nrow(gam_dt) > 0) {
    gam_dt[, group := factor(group, levels = group_order)]
    gam_dt[, metric := factor(metric, levels = metric_order)]
  }

  panel_labels <- make_panel_labels(length(group_order) * length(metric_order))

  metric_blocks <- list()
  k <- 0L

  for (i in seq_along(metric_order)) {
    met <- metric_order[i]

    panel_plots <- list()

    for (j in seq_along(group_order)) {
      grp <- group_order[j]
      k <- k + 1L

      dt_curve <- curve_dt[group == grp & metric == met][order(op_bin)]
      dt_gam <- gam_dt[group == grp & metric == met][order(op_bin)]

      panel_plots[[j]] <- make_cell_m4(
        dt_curve = dt_curve,
        dt_gam = dt_gam,
        ylab_left = if ((j - 1L) %% panel_ncol == 0L) "Mean divergence" else NULL,
        title = wrap_group_title(as.character(grp), width = title_wrap_width),
        panel_label = panel_labels[k],
        show_x = TRUE,
        show_left_y = ((j - 1L) %% panel_ncol == 0L),
        show_right_y = (j %% panel_ncol == 0L || j == length(group_order))
      ) +
        theme(
          plot.title = element_text(
            size = title_size,
            face = "bold",
            hjust = 0.5,
            lineheight = 0.9,
            margin = margin(b = 8)
          )
        )
    }

    panel_grid <- plot_grid(
      plotlist = panel_plots,
      ncol = panel_ncol,
      align = "hv",
      axis = "tblr"
    )

    metric_blocks[[i]] <- plot_grid(
      make_metric_strip(metric_labels[[met]], size = strip_size),
      panel_grid,
      nrow = 1,
      rel_widths = c(0.055, 1)
    )
  }

  final_m4 <- plot_grid(
    ggdraw() + draw_label(title_text, fontface = "bold", size = 15, x = 0, hjust = 0),
    plot_grid(plotlist = metric_blocks, ncol = 1, align = "v", axis = "lr"),
    ggdraw() + draw_label("Operational time (hours since first observation)", size = 11, x = 0.5),
    ggdraw() + draw_label(note_text, size = 8.5, x = 0, hjust = 0),
    ncol = 1,
    rel_heights = c(0.045, 0.90, 0.025, 0.03)
  )

  ggsave(outfile, final_m4, width = width, height = height, dpi = 600, limitsize = FALSE)

  rm(curve_dt, gam_dt, final_m4, metric_blocks)
  gc()

  say("Wrote figure to: ", outfile)
}

make_pair_grid <- function(curve_dt, gam_dt, group_order, outfile, title_text, note_text,
                           width = NULL, height = NULL, panel_ncol = NULL) {
  curve_dt <- copy(curve_dt)
  gam_dt <- copy(gam_dt)

  n_groups <- length(group_order)
  if (is.null(panel_ncol)) panel_ncol <- auto_panel_ncol(n_groups)

  if (is.null(width)) width <- auto_fig_width(n_groups, panel_ncol)
  if (is.null(height)) height <- auto_fig_height(length(paired_metric_order), n_groups, panel_ncol)

  title_wrap_width <- ifelse(panel_ncol <= 4, 20, 16)
  title_size <- ifelse(panel_ncol <= 4, 8.2, 7.2)
  strip_size <- 9.5

  curve_dt[, group := factor(group, levels = group_order)]
  curve_dt[, metric := factor(metric, levels = paired_metric_order)]

  if (nrow(gam_dt) > 0) {
    gam_dt[, group := factor(group, levels = group_order)]
    gam_dt[, metric := factor(metric, levels = paired_metric_order)]
  }

  panel_labels <- make_panel_labels(length(group_order) * length(paired_metric_order))

  metric_blocks <- list()
  k <- 0L

  for (i in seq_along(paired_metric_order)) {
    met <- paired_metric_order[i]

    panel_plots <- list()

    for (j in seq_along(group_order)) {
      grp <- group_order[j]
      k <- k + 1L

      dt_curve <- curve_dt[group == grp & metric == met][order(source, op_bin)]
      dt_gam <- gam_dt[group == grp & metric == met][order(source, op_bin)]

      panel_plots[[j]] <- make_cell_pair(
        dt_curve = dt_curve,
        dt_gam = dt_gam,
        ylab_left = if ((j - 1L) %% panel_ncol == 0L) "Mean value" else NULL,
        title = wrap_group_title(as.character(grp), width = title_wrap_width),
        panel_label = panel_labels[k],
        show_x = TRUE,
        show_left_y = ((j - 1L) %% panel_ncol == 0L),
        show_right_y = (j %% panel_ncol == 0L || j == length(group_order))
      ) +
        theme(
          plot.title = element_text(
            size = title_size,
            face = "bold",
            hjust = 0.5,
            lineheight = 0.9,
            margin = margin(b = 8)
          )
        )
    }

    panel_grid <- plot_grid(
      plotlist = panel_plots,
      ncol = panel_ncol,
      align = "hv",
      axis = "tblr"
    )

    metric_blocks[[i]] <- plot_grid(
      make_metric_strip(paired_metric_labels[[met]], size = strip_size),
      panel_grid,
      nrow = 1,
      rel_widths = c(0.055, 1)
    )
  }

  final_pair <- plot_grid(
    ggdraw() + draw_label(title_text, fontface = "bold", size = 15, x = 0, hjust = 0),
    plot_grid(plotlist = metric_blocks, ncol = 1, align = "v", axis = "lr"),
    ggdraw() + draw_label("Operational time (hours since first observation)", size = 11, x = 0.5),
    ggdraw() + draw_label(note_text, size = 8.5, x = 0, hjust = 0),
    ncol = 1,
    rel_heights = c(0.045, 0.90, 0.025, 0.03)
  )

  ggsave(outfile, final_pair, width = width, height = height, dpi = 600, limitsize = FALSE)

  rm(curve_dt, gam_dt, final_pair, metric_blocks)
  gc()

  say("Wrote figure to: ", outfile)
}

# ------------------------------------------------------------
# 8) GLOBAL + REGION FIGURES
# ------------------------------------------------------------

region_order <- c("Global", sort(setdiff(unique(curve_m4_global_region$group), "Global")))

make_metric_grid(
  curve_dt = curve_m4_global_region,
  gam_dt = gam_global_region,
  group_order = region_order,
  outfile = file.path(OUTDIR, "metric4_divergence_global_region_sensor_ci.jpeg"),
  title_text = "Metric 4: Age-dependent divergence versus external references by region",
  note_text = "Black line = binned sensor-level mean; grey band = sensor-level 95% CI across sensors; blue line = GAM smooth only; purple dot-dashed line = number of sensors. Global and regional panels use outdoor_final.RData. No physical-range filtering applied.",
  width = 13,
  height = 9
)

make_pair_grid(
  curve_dt = pair_m4_global_region,
  gam_dt = gam_pair_global_region,
  group_order = region_order,
  outfile = file.path(OUTDIR, "metric4_paired_values_global_region_sensor_ci.jpeg"),
  title_text = "Metric 4: Paired PurpleAir and reference trajectories by region",
  note_text = "Lines = binned sensor-level means; ribbons = sensor-level 95% CI across sensors; smooth lines = GAM fits only; red line = number of sensors. Global and regional panels use outdoor_final.RData. No physical-range filtering applied.",
  width = 13,
  height = 9
)

# ------------------------------------------------------------
# 9) GLOBAL + ALL CLIMATE-ZONE FIGURES
# ------------------------------------------------------------

climate_order <- as.character(
  curve_m4_climate[
    ,
    .(total_sensor_bins = sum(n_sensor_bins, na.rm = TRUE)),
    by = group
  ][order(-total_sensor_bins), group]
)

curve_m4_global_for_climate <- copy(
  curve_m4_global[, .(group, metric, op_bin, op_years, n_sensors, n_sensor_bins, n_obs, mean, lwr, upr)]
)

curve_m4_global_plus_climate <- rbindlist(
  list(
    curve_m4_global_for_climate,
    curve_m4_climate[, .(group, metric, op_bin, op_years, n_sensors, n_sensor_bins, n_obs, mean, lwr, upr)]
  ),
  use.names = TRUE,
  fill = TRUE
)

gam_global_for_climate <- make_gam_by_group_metric(curve_m4_global_for_climate)

gam_global_plus_climate <- rbindlist(
  list(gam_global_for_climate, gam_climate),
  use.names = TRUE,
  fill = TRUE
)

pair_global_for_climate <- copy(
  pair_m4_global[, .(group, metric, source, op_bin, op_years, n_sensors, n_sensor_bins, n_obs, mean, lwr, upr)]
)

pair_global_plus_climate <- rbindlist(
  list(
    pair_global_for_climate,
    pair_m4_climate[, .(group, metric, source, op_bin, op_years, n_sensors, n_sensor_bins, n_obs, mean, lwr, upr)]
  ),
  use.names = TRUE,
  fill = TRUE
)

gam_pair_global_for_climate <- make_gam_by_group_metric_source(pair_global_for_climate)

gam_pair_global_plus_climate <- rbindlist(
  list(gam_pair_global_for_climate, gam_pair_climate),
  use.names = TRUE,
  fill = TRUE
)

make_metric_grid(
  curve_dt = curve_m4_global_plus_climate,
  gam_dt = gam_global_plus_climate,
  group_order = c("Global", climate_order),
  outfile = file.path(OUTDIR, "metric4_divergence_global_all_climate_sensor_ci.jpeg"),
  title_text = "Metric 4: Age-dependent divergence versus external references by climate zone",
  note_text = "Black line = binned sensor-level mean; grey band = sensor-level 95% CI across sensors; blue line = GAM smooth only; red line = number of sensors. Global panel uses outdoor_final.RData; climate panels use valid single-zone stations from outdoor_with_koppen_timevarying.RData. No physical-range filtering applied."
)

make_pair_grid(
  curve_dt = pair_global_plus_climate,
  gam_dt = gam_pair_global_plus_climate,
  group_order = c("Global", climate_order),
  outfile = file.path(OUTDIR, "metric4_paired_values_global_all_climate_sensor_ci.jpeg"),
  title_text = "Metric 4: Paired PurpleAir and reference trajectories by climate zone",
  note_text = "Lines = binned sensor-level means; ribbons = sensor-level 95% CI across sensors; smooth lines = GAM fits only; red line = number of sensors. Global panel uses outdoor_final.RData; climate panels use valid single-zone stations from outdoor_with_koppen_timevarying.RData. No physical-range filtering applied."
)

# ------------------------------------------------------------
# 10) REGION x CLIMATE FIGURES
# ------------------------------------------------------------

for (rr in sort(unique(curve_m4_region_climate_plot$region))) {

  curve_rr <- curve_m4_region_climate_plot[region == rr]
  pair_rr <- pair_m4_region_climate_plot[region == rr]

  group_order_rr <- as.character(
    curve_rr[
      ,
      .(total_sensor_bins = sum(n_sensor_bins, na.rm = TRUE)),
      by = group
    ][order(-total_sensor_bins), group]
  )

  if (length(group_order_rr) == 0L) next

  gam_rr <- gam_region_climate[group %in% group_order_rr]
  gam_pair_rr <- gam_pair_region_climate[group %in% group_order_rr]

  make_metric_grid(
    curve_dt = curve_rr,
    gam_dt = gam_rr,
    group_order = group_order_rr,
    outfile = file.path(
      OUTDIR,
      paste0("metric4_divergence_region_x_climate_sensor_ci_", safe_file(rr), ".jpeg")
    ),
    title_text = paste0("Metric 4: Age-dependent divergence versus external references by climate zone: ", rr),
    note_text = paste0(
      "Black line = binned sensor-level mean; grey band = sensor-level 95% CI across sensors; blue line = GAM smooth only; red line = number of sensors. ",
      "Region × climate panels use valid single-zone stations from outdoor_with_koppen_timevarying.RData. No physical-range filtering applied. Panels shown require ≥",
      MIN_TOTAL_SENSOR_BINS_REGION_CLIMATE,
      " sensor-bins and ≥",
      MIN_UNIQUE_BINS_REGION_CLIMATE,
      " operational-age bins."
    )
  )

  make_pair_grid(
    curve_dt = pair_rr,
    gam_dt = gam_pair_rr,
    group_order = group_order_rr,
    outfile = file.path(
      OUTDIR,
      paste0("metric4_paired_values_region_x_climate_sensor_ci_", safe_file(rr), ".jpeg")
    ),
    title_text = paste0("Metric 4: Paired PurpleAir and reference trajectories by climate zone: ", rr),
    note_text = paste0(
      "Lines = binned sensor-level means; ribbons = sensor-level 95% CI across sensors; smooth lines = GAM fits only; red line = number of sensors. ",
      "Region × climate panels use valid single-zone stations from outdoor_with_koppen_timevarying.RData. No physical-range filtering applied. Panels shown require ≥",
      MIN_TOTAL_SENSOR_BINS_REGION_CLIMATE,
      " sensor-bins and ≥",
      MIN_UNIQUE_BINS_REGION_CLIMATE,
      " operational-age bins."
    )
  )
}
# ------------------------------------------------------------
# 11) PREVIEW
# ------------------------------------------------------------

cat("\ncurve_m4_global divergence preview:\n")
print(curve_m4_global[1:12])

cat("\npaired global preview:\n")
print(pair_m4_global[1:12])

cat("\ncurve_m4_climate divergence preview:\n")
print(curve_m4_climate[1:12])

cat("\npaired climate preview:\n")
print(pair_m4_climate[1:12])
