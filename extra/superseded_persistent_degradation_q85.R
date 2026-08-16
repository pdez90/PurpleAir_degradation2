# ==========================================================================
# superseded_persistent_degradation_q85.R
#
# SUPERSEDED - do not run. Earlier persistent/permanent degradation pass using
# FLAG_REL_Q = 0.85. The paper reports the stringent (0.90) version in
# R/18_persistent_degradation_stringent.R.
#
# Extracted VERBATIM from purpleair_degradation2.Rmd (chunk #37 lines 21811-23034).
# Code is unmodified; see REPRODUCING.md for paths and run order.
# ==========================================================================

# ---- source chunk #37 (Rmd lines 21811-23034) ------------------------

# ============================================================
# IDENTIFY POTENTIALLY PERMANENTLY DEGRADED PURPLEAIR SENSORS
# Hourly-record version; memory-conscious outputs
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(zoo)
})

global_rdata  <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_final.RData"
climate_rdata <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_with_koppen_timevarying.RData"
transition_csv <- "/Users/priyanka/Downloads/PurpleAir_EPA/station_koppen_transition_summary.csv"

OUTDIR <- "/Users/priyanka/Downloads/PurpleAir_EPA/permanent_degradation_hourly"
FIGDIR <- file.path(OUTDIR, "figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

ROLL_DAYS <- 90L
ROLL_HOURS <- ROLL_DAYS * 24L

PERSIST_DAYS <- 90L
PERSIST_HOURS <- PERSIST_DAYS * 24L

RECOVERY_DAYS <- 90L
RECOVERY_HOURS <- RECOVERY_DAYS * 24L

MIN_HOURS_SENSOR <- 24L * 180L

FLAG_ABS_DIFF <- 5
FLAG_REL_Q <- 0.85

DEGRADED_FLAG_RATE <- 0.25
RECOVERED_FLAG_RATE <- 0.10

DEGRADED_REL_AB <- 0.30
RECOVERED_REL_AB <- 0.15

DEGRADED_ABS_AB <- 10

HIGH_PM_THRESHOLDS <- c(35, 55, 100, 150)

CHG_MIN_SEG_HOURS <- 24L * 30L
CHG_MEAN_SHIFT_MIN <- 0.10

SAVE_HOURLY_PROCESSED <- FALSE

a_col_default <- "pm25_cf_1_a_hourly(ug/m3)_pm25_cf_1_a"
b_col_default <- "pm25_cf_1_b_hourly(ug/m3)_pm25_cf_1_b"
geos_col <- "geoscf_pm25"

say <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste0(..., collapse = ""), "\n")
}

find_col <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) == 0L) return(NA_character_)
  hit[1]
}

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

empty_episode_table <- function() {
  data.table(
    station = character(),
    metric = character(),
    episode_id = integer(),
    start_index = integer(),
    end_index = integer(),
    start_hour = numeric(),
    end_hour = numeric(),
    duration_hours = numeric()
  )
}

empty_sensor_summary <- function() {
  data.table(
    station = character(),
    region = character(),
    climate_zone = character(),
    climate_desc = character(),
    climate_label = character(),
    lat = numeric(),
    lon = numeric(),
    n_hours = integer(),
    max_hour = numeric(),
    max_year = numeric(),
    mean_flagged = numeric(),
    mean_rel_ab = numeric(),
    mean_abs_ab = numeric(),
    mean_pa_geos_abs = numeric(),
    max_roll_flag_rate = numeric(),
    max_roll_rel_ab = numeric(),
    max_roll_abs_ab = numeric(),
    has_persistent_degradation = logical(),
    permanent_degradation = logical(),
    persistent_onset_hour = numeric(),
    persistent_onset_year = numeric(),
    recovery_hour = numeric(),
    recovery_year = numeric(),
    persistence_score = numeric(),
    effective_lifetime_hour = numeric(),
    effective_lifetime_year = numeric()
  )
}

make_hourly_degradation_dt <- function(infile, require_climate = FALSE, label = "global") {

  say("Loading hourly data: ", label)
  outdoor <- load_outdoor_object(infile)

  a_col <- find_col(names(outdoor), c(
    a_col_default,
    "pm25_cf_1_a_hourly(ug/m3)_pm25_cf_1_a",
    "PM2.5_CF_1_ug/m3_A",
    "cf1_a",
    "pm_a"
  ))

  b_col <- find_col(names(outdoor), c(
    b_col_default,
    "pm25_cf_1_b_hourly(ug/m3)_pm25_cf_1_b",
    "PM2.5_CF_1_ug/m3_B",
    "cf1_b",
    "pm_b"
  ))

  lat_col <- find_col(names(outdoor), c(
    "LATITUDE(deg)",
    "lat", "Lat", "latitude", "Latitude",
    "gee_latitude", "sensor_latitude"
  ))

  lon_col <- find_col(names(outdoor), c(
    "LONGITUDE(deg)",
    "lon", "Lon", "lng", "longitude", "Longitude",
    "gee_longitude", "sensor_longitude"
  ))

  time_col <- find_col(names(outdoor), c(
    "datetime_utc",
    "ts_hr",
    "date",
    "datetime",
    "timestamp",
    "time"
  ))

  needed <- c("station", "region", time_col, lat_col, lon_col, a_col, b_col)

  has_geos <- geos_col %in% names(outdoor)
  if (has_geos) needed <- c(needed, geos_col)

  if (any(is.na(needed))) {
    stop(
      "Could not identify required columns.\n\n",
      "Identified:\n",
      "time_col = ", time_col, "\n",
      "lat_col = ", lat_col, "\n",
      "lon_col = ", lon_col, "\n",
      "a_col = ", a_col, "\n",
      "b_col = ", b_col, "\n",
      "geos_col_present = ", has_geos, "\n\n",
      "Available columns:\n",
      paste(names(outdoor), collapse = ", ")
    )
  }

  miss <- setdiff(needed, names(outdoor))
  if (length(miss) > 0L) {
    stop("Missing required columns: ", paste(miss, collapse = ", "))
  }

  say("Using time column: ", time_col)
  say("Using coordinate columns: ", lat_col, ", ", lon_col)
  say("Using PM columns: ", a_col, ", ", b_col)

  outdoor[, station := as.character(station)]

  if (require_climate) {
    outdoor <- outdoor[station %in% single_zone$station]
    say("Rows after single-zone restriction: ", format(nrow(outdoor), big.mark = ","))
  }

  dt <- outdoor[, ..needed]
  rm(outdoor)
  gc()

  setnames(dt, time_col, "datetime_utc")
  setnames(dt, a_col, "cf1_a")
  setnames(dt, b_col, "cf1_b")
  setnames(dt, lat_col, "lat")
  setnames(dt, lon_col, "lon")

  dt[, `:=`(
    station = as.character(station),
    region = as.character(region),
    ts = as.POSIXct(datetime_utc, tz = "UTC"),
    lat = as.numeric(lat),
    lon = as.numeric(lon),
    cf1_a = as.numeric(cf1_a),
    cf1_b = as.numeric(cf1_b)
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
      is.finite(lat) &
      is.finite(lon) &
      is.finite(cf1_a) &
      is.finite(cf1_b)
  ]

  dt[, pm_cf1 := (cf1_a + cf1_b) / 2]
  dt[, abs_ab := abs(cf1_a - cf1_b)]
  dt[, rel_ab := fifelse((cf1_a + cf1_b) > 0, 2 * abs_ab / (cf1_a + cf1_b), NA_real_)]
  dt[, bias_ab := cf1_a - cf1_b]
  dt[, pa_geos_abs := fifelse(is.finite(geos_pm), abs(pm_cf1 - geos_pm), NA_real_)]
  dt[, pa_geos_signed := fifelse(is.finite(geos_pm), pm_cf1 - geos_pm, NA_real_)]

  setorder(dt, station, ts)

  thr <- dt[
    is.finite(rel_ab),
    .(rel_thr = quantile(rel_ab, FLAG_REL_Q, na.rm = TRUE, type = 7)),
    by = station
  ]

  dt <- thr[dt, on = "station"]

  dt[, flagged := as.integer(abs_ab > FLAG_ABS_DIFF & rel_ab > rel_thr)]
  dt[is.na(flagged), flagged := 0L]

  dt[, first_ts := ts[1L], by = station]
  dt[, hour := as.numeric(difftime(ts, first_ts, units = "hours"))]
  dt <- dt[is.finite(hour) & hour >= 0]

  if (require_climate) {
    dt <- merge(dt, single_zone, by = "station", all.x = FALSE)
  } else {
    dt[, `:=`(
      climate_zone = NA_character_,
      climate_desc = NA_character_,
      climate_label = NA_character_
    )]
  }

  keep <- c(
    "station", "region", "climate_zone", "climate_desc", "climate_label",
    "ts", "hour", "lat", "lon", "pm_cf1", "abs_ab", "bias_ab",
    "rel_ab", "flagged", "pa_geos_abs", "pa_geos_signed"
  )

  dt <- dt[, ..keep]

  say("Hourly rows retained: ", format(nrow(dt), big.mark = ","))
  say("Stations retained: ", format(uniqueN(dt$station), big.mark = ","))

  dt[]
}

add_rolling_metrics <- function(dt) {

  say("Adding rolling degradation metrics")
  setorder(dt, station, hour)

  dt[, roll_flag_rate := frollmean(flagged, n = ROLL_HOURS, align = "right", na.rm = TRUE), by = station]
  dt[, roll_rel_ab := frollmean(rel_ab, n = ROLL_HOURS, align = "right", na.rm = TRUE), by = station]
  dt[, roll_abs_ab := frollmean(abs_ab, n = ROLL_HOURS, align = "right", na.rm = TRUE), by = station]
  dt[, roll_pa_geos_abs := frollmean(pa_geos_abs, n = ROLL_HOURS, align = "right", na.rm = TRUE), by = station]
  dt[, roll_pm := frollmean(pm_cf1, n = ROLL_HOURS, align = "right", na.rm = TRUE), by = station]

  for (thr in HIGH_PM_THRESHOLDS) {
    nm <- paste0("highpm_", thr)
    cn <- paste0("cum_highpm_", thr)
    en <- paste0("cum_excesspm_", thr)

    dt[, (nm) := as.integer(pm_cf1 > thr)]
    dt[, (cn) := cumsum(get(nm)), by = station]
    dt[, (en) := cumsum(pmax(pm_cf1 - thr, 0)), by = station]
  }

  dt[]
}

episode_summary <- function(x, hour, station_id, metric_name) {
  x <- as.logical(x)
  x[is.na(x)] <- FALSE
  hour <- as.numeric(hour)

  if (length(x) == 0L || length(hour) == 0L || !any(x)) {
    return(empty_episode_table())
  }

  r <- rle(x)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L
  idx <- which(r$values)

  if (length(idx) == 0L) return(empty_episode_table())

  out <- data.table(
    station = as.character(station_id),
    metric = metric_name,
    episode_id = seq_along(idx),
    start_index = starts[idx],
    end_index = ends[idx]
  )

  out[, `:=`(
    start_hour = hour[start_index],
    end_hour = hour[end_index]
  )]

  out[, duration_hours := end_hour - start_hour + 1]

  out <- out[
    is.finite(start_hour) &
      is.finite(end_hour) &
      is.finite(duration_hours)
  ]

  if (nrow(out) == 0L) return(empty_episode_table())

  out[]
}

classify_one_sensor <- function(d) {

  st <- as.character(d$station[1])
  setorder(d, hour)

  if (nrow(d) < MIN_HOURS_SENSOR) return(NULL)

  degraded_flag <- d$roll_flag_rate > DEGRADED_FLAG_RATE
  degraded_rel  <- d$roll_rel_ab > DEGRADED_REL_AB
  degraded_abs  <- d$roll_abs_ab > DEGRADED_ABS_AB

  degraded_any <- degraded_flag | degraded_rel | degraded_abs
  degraded_any[is.na(degraded_any)] <- FALSE

  ep <- rbindlist(
    list(
      episode_summary(degraded_flag, d$hour, st, "rolling_flag_rate"),
      episode_summary(degraded_rel,  d$hour, st, "rolling_relative_AB"),
      episode_summary(degraded_abs,  d$hour, st, "rolling_absolute_AB"),
      episode_summary(degraded_any,  d$hour, st, "any_degradation")
    ),
    use.names = TRUE,
    fill = TRUE
  )

  if (nrow(ep) > 0L) {
    ep[, `:=`(
      region = d$region[1],
      climate_zone = d$climate_zone[1],
      climate_desc = d$climate_desc[1],
      climate_label = d$climate_label[1]
    )]
  }

  persistent_ep <- ep[
    metric == "any_degradation" &
      is.finite(duration_hours) &
      duration_hours >= PERSIST_HOURS
  ]

  has_persistent <- nrow(persistent_ep) > 0L
  onset_hour <- if (has_persistent) persistent_ep$start_hour[1] else NA_real_

  permanent <- FALSE
  recovery_hour <- NA_real_

  if (has_persistent && is.finite(onset_hour)) {
    future <- d[hour >= onset_hour]

    recovered_vec <- future$roll_flag_rate < RECOVERED_FLAG_RATE &
      future$roll_rel_ab < RECOVERED_REL_AB
    recovered_vec[is.na(recovered_vec)] <- FALSE

    rec_ep <- episode_summary(recovered_vec, future$hour, st, "recovery")
    rec_ep <- rec_ep[is.finite(duration_hours) & duration_hours >= RECOVERY_HOURS]

    if (nrow(rec_ep) > 0L) {
      recovery_hour <- rec_ep$start_hour[1]
      permanent <- FALSE
    } else {
      permanent <- TRUE
    }
  }

  persistence_score <- NA_real_
  if (has_persistent && is.finite(onset_hour)) {
    persistence_score <- mean(degraded_any[d$hour >= onset_hour], na.rm = TRUE)
  }

  effective_lifetime_hour <- if (has_persistent && is.finite(onset_hour)) {
    onset_hour
  } else {
    max(d$hour, na.rm = TRUE)
  }

  before_onset <- if (has_persistent && is.finite(onset_hour)) {
    d[hour <= onset_hour]
  } else {
    d
  }

  summary <- data.table(
    station = st,
    region = d$region[1],
    climate_zone = d$climate_zone[1],
    climate_desc = d$climate_desc[1],
    climate_label = d$climate_label[1],
    lat = median(d$lat, na.rm = TRUE),
    lon = median(d$lon, na.rm = TRUE),
    n_hours = nrow(d),
    max_hour = max(d$hour, na.rm = TRUE),
    max_year = max(d$hour, na.rm = TRUE) / (24 * 365.25),
    mean_flagged = mean(d$flagged, na.rm = TRUE),
    mean_rel_ab = mean(d$rel_ab, na.rm = TRUE),
    mean_abs_ab = mean(d$abs_ab, na.rm = TRUE),
    mean_pa_geos_abs = mean(d$pa_geos_abs, na.rm = TRUE),
    max_roll_flag_rate = max(d$roll_flag_rate, na.rm = TRUE),
    max_roll_rel_ab = max(d$roll_rel_ab, na.rm = TRUE),
    max_roll_abs_ab = max(d$roll_abs_ab, na.rm = TRUE),
    has_persistent_degradation = has_persistent,
    permanent_degradation = permanent,
    persistent_onset_hour = onset_hour,
    persistent_onset_year = onset_hour / (24 * 365.25),
    recovery_hour = recovery_hour,
    recovery_year = recovery_hour / (24 * 365.25),
    persistence_score = persistence_score,
    effective_lifetime_hour = effective_lifetime_hour,
    effective_lifetime_year = effective_lifetime_hour / (24 * 365.25)
  )

  for (thr in HIGH_PM_THRESHOLDS) {
    high_col <- paste0("highpm_", thr)
    excess_col <- paste0("cum_excesspm_", thr)

    summary[, paste0("pre_onset_n_highpm_", thr) :=
              sum(before_onset[[high_col]], na.rm = TRUE)]

    summary[, paste0("pre_onset_excesspm_", thr) :=
              max(before_onset[[excess_col]], na.rm = TRUE)]
  }

  list(summary = summary, episodes = ep)
}

classify_all_sensors <- function(dt, label = "global_region") {

  say("Classifying persistent/permanent degradation: ", label)

  setorder(dt, station, hour)

  stations <- unique(dt$station)
  n_st <- length(stations)

  summary_list <- vector("list", n_st)
  episode_list <- vector("list", n_st)

  failed <- 0L
  pb <- txtProgressBar(min = 0, max = n_st, style = 3)

  for (ii in seq_along(stations)) {
    st <- stations[ii]
    d_st <- dt[station == st]

    res <- tryCatch(
      classify_one_sensor(d_st),
      error = function(e) {
        failed <<- failed + 1L
        if (failed <= 10L) {
          message("Classification failed for station ", st, ": ", conditionMessage(e))
        }
        NULL
      }
    )

    if (!is.null(res)) {
      summary_list[[ii]] <- res$summary
      episode_list[[ii]] <- res$episodes
    }

    if (ii %% 100 == 0L || ii == n_st) setTxtProgressBar(pb, ii)
    if (ii %% 500 == 0L) gc()
  }

  close(pb)

  sensor_summary <- rbindlist(summary_list, use.names = TRUE, fill = TRUE)
  episode_summary_dt <- rbindlist(episode_list, use.names = TRUE, fill = TRUE)

  if (nrow(sensor_summary) == 0L) sensor_summary <- empty_sensor_summary()
  if (nrow(episode_summary_dt) == 0L) episode_summary_dt <- empty_episode_table()

  sensor_summary[, analysis_label := label]
  episode_summary_dt[, analysis_label := label]

  fwrite(sensor_summary, file.path(OUTDIR, paste0("permanent_degradation_sensor_summary_", label, ".csv")))
  fwrite(episode_summary_dt, file.path(OUTDIR, paste0("permanent_degradation_episode_summary_", label, ".csv")))

  say("Classified sensors: ", format(nrow(sensor_summary), big.mark = ","))
  say("Episodes detected: ", format(nrow(episode_summary_dt), big.mark = ","))
  say("Classification failures: ", format(failed, big.mark = ","))

  list(summary = sensor_summary, episodes = episode_summary_dt)
}

detect_changepoint_one_sensor <- function(d) {

  st <- d$station[1]

  d <- d[is.finite(roll_flag_rate) & is.finite(hour)]
  if (nrow(d) < MIN_HOURS_SENSOR) return(NULL)

  d[, day_bin := floor(hour / 24) * 24]

  dd <- d[
    ,
    .(
      roll_flag_rate = mean(roll_flag_rate, na.rm = TRUE),
      roll_rel_ab = mean(roll_rel_ab, na.rm = TRUE),
      roll_abs_ab = mean(roll_abs_ab, na.rm = TRUE)
    ),
    by = day_bin
  ][order(day_bin)]

  if (nrow(dd) < 180L) return(NULL)

  y <- dd$roll_flag_rate
  n <- length(y)
  min_seg <- max(30L, floor(CHG_MIN_SEG_HOURS / 24))

  if (n < 2 * min_seg) return(NULL)

  candidates <- seq(min_seg, n - min_seg)

  cs <- cumsum(y)
  cs2 <- cumsum(y^2)

  total_sum <- cs[n]
  total_sum2 <- cs2[n]

  best_sse <- Inf
  best_k <- NA_integer_
  best_before <- NA_real_
  best_after <- NA_real_

  for (k in candidates) {
    n1 <- k
    n2 <- n - k

    s1 <- cs[k]
    s2 <- total_sum - cs[k]

    ss1 <- cs2[k]
    ss2 <- total_sum2 - cs2[k]

    m1 <- s1 / n1
    m2 <- s2 / n2

    sse <- (ss1 - n1 * m1^2) + (ss2 - n2 * m2^2)

    if (is.finite(sse) && sse < best_sse) {
      best_sse <- sse
      best_k <- k
      best_before <- m1
      best_after <- m2
    }
  }

  if (!is.finite(best_sse) || is.na(best_k)) return(NULL)

  shift <- best_after - best_before

  data.table(
    station = st,
    region = d$region[1],
    climate_zone = d$climate_zone[1],
    climate_desc = d$climate_desc[1],
    climate_label = d$climate_label[1],
    changepoint_hour = dd$day_bin[best_k],
    changepoint_year = dd$day_bin[best_k] / (24 * 365.25),
    mean_before = best_before,
    mean_after = best_after,
    mean_shift = shift,
    likely_degradation_changepoint =
      is.finite(shift) &&
      shift >= CHG_MEAN_SHIFT_MIN &&
      best_after >= DEGRADED_FLAG_RATE
  )
}

detect_changepoints_all <- function(dt, label) {

  say("Running change-point detection: ", label)

  ids <- unique(dt$station)
  out <- vector("list", length(ids))

  for (i in seq_along(ids)) {
    if (i %% 500 == 0) say("  changepoints ", i, " / ", length(ids))

    d <- dt[station == ids[i]]
    out[[i]] <- detect_changepoint_one_sensor(d)

    rm(d)
    if (i %% 500 == 0) gc()
  }

  out <- rbindlist(out, fill = TRUE)

  if (nrow(out) == 0L) {
    out <- data.table(
      station = character(),
      region = character(),
      climate_zone = character(),
      climate_desc = character(),
      climate_label = character(),
      changepoint_hour = numeric(),
      changepoint_year = numeric(),
      mean_before = numeric(),
      mean_after = numeric(),
      mean_shift = numeric(),
      likely_degradation_changepoint = logical()
    )
  }

  fwrite(out, file.path(OUTDIR, paste0("changepoint_results_", label, ".csv")))

  out[]
}

assign_rule_states <- function(dt) {
  dt[, state_rule := fifelse(
    roll_flag_rate >= DEGRADED_FLAG_RATE | roll_rel_ab >= DEGRADED_REL_AB,
    "Degraded",
    fifelse(
      roll_flag_rate >= RECOVERED_FLAG_RATE | roll_rel_ab >= RECOVERED_REL_AB,
      "Noisy",
      "Healthy"
    )
  )]

  dt[]
}

fit_hmm_if_available <- function(dt, label, max_rows = 250000L) {

  if (!requireNamespace("depmixS4", quietly = TRUE)) {
    say("depmixS4 not installed; using rule-based HMM-style states.")
    return(NULL)
  }

  say("Fitting sampled 3-state HMM using depmixS4: ", label)

  d <- dt[
    is.finite(roll_flag_rate) &
      is.finite(roll_rel_ab) &
      is.finite(roll_abs_ab)
  ]

  if (nrow(d) > max_rows) {
    set.seed(123)
    d <- d[sample.int(nrow(d), max_rows)]
  }

  d[, station_num := as.integer(factor(station))]
  setorder(d, station_num, hour)

  mod <- depmixS4::depmix(
    list(
      roll_flag_rate ~ 1,
      roll_rel_ab ~ 1
    ),
    data = d,
    nstates = 3,
    family = list(gaussian(), gaussian()),
    ntimes = as.numeric(d[, .N, by = station_num]$N)
  )

  fit <- tryCatch(
    depmixS4::fit(mod, verbose = FALSE),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NULL)

  capture.output(
    summary(fit),
    file = file.path(OUTDIR, paste0("hmm_depmixS4_summary_", label, ".txt"))
  )

  invisible(fit)
}

summarize_states <- function(dt, label) {

  say("Summarizing states: ", label)

  state_summary <- dt[
    !is.na(state_rule),
    .(
      n_hours = .N,
      healthy_hours = sum(state_rule == "Healthy", na.rm = TRUE),
      noisy_hours = sum(state_rule == "Noisy", na.rm = TRUE),
      degraded_hours = sum(state_rule == "Degraded", na.rm = TRUE),
      pct_degraded = mean(state_rule == "Degraded", na.rm = TRUE),
      pct_noisy_or_degraded = mean(state_rule %in% c("Noisy", "Degraded"), na.rm = TRUE),
      ever_degraded = any(state_rule == "Degraded"),
      final_state = state_rule[.N],
      final_90d_degraded_share = mean(tail(state_rule, ROLL_HOURS) == "Degraded", na.rm = TRUE)
    ),
    by = .(station, region, climate_zone, climate_desc, climate_label)
  ]

  state_summary[, hmm_like_permanent_degraded :=
                  final_90d_degraded_share >= 0.80 &
                  pct_degraded >= 0.20]

  fwrite(state_summary, file.path(OUTDIR, paste0("hmm_like_state_summary_", label, ".csv")))

  state_summary[]
}

make_trajectory_classes <- function(sensor_summary, changepoints, state_summary, label) {

  say("Creating joint trajectory degradation classes: ", label)

  if (!("station" %in% names(sensor_summary))) sensor_summary <- empty_sensor_summary()

  if (!("station" %in% names(changepoints))) {
    changepoints <- data.table(
      station = character(),
      likely_degradation_changepoint = logical(),
      changepoint_hour = numeric(),
      changepoint_year = numeric(),
      mean_shift = numeric()
    )
  }

  if (!("station" %in% names(state_summary))) {
    state_summary <- data.table(
      station = character(),
      pct_degraded = numeric(),
      pct_noisy_or_degraded = numeric(),
      final_90d_degraded_share = numeric(),
      hmm_like_permanent_degraded = logical()
    )
  }

  x <- merge(
    sensor_summary,
    changepoints[, .(
      station,
      likely_degradation_changepoint,
      changepoint_hour,
      changepoint_year,
      mean_shift
    )],
    by = "station",
    all.x = TRUE
  )

  x <- merge(
    x,
    state_summary[, .(
      station,
      pct_degraded,
      pct_noisy_or_degraded,
      final_90d_degraded_share,
      hmm_like_permanent_degraded
    )],
    by = "station",
    all.x = TRUE
  )

  x[is.na(likely_degradation_changepoint), likely_degradation_changepoint := FALSE]
  x[is.na(hmm_like_permanent_degraded), hmm_like_permanent_degraded := FALSE]

  x[, degradation_class := fifelse(
    permanent_degradation == TRUE & likely_degradation_changepoint == TRUE,
    "Likely permanently degraded with changepoint",
    fifelse(
      permanent_degradation == TRUE,
      "Likely permanently degraded",
      fifelse(
        has_persistent_degradation == TRUE & !permanent_degradation,
        "Persistent but recovered",
        fifelse(
          hmm_like_permanent_degraded == TRUE,
          "HMM-like final degraded",
          fifelse(
            max_roll_flag_rate >= DEGRADED_FLAG_RATE |
              max_roll_rel_ab >= DEGRADED_REL_AB,
            "Intermittently noisy/degraded",
            "No strong degradation evidence"
          )
        )
      )
    )
  )]

  fwrite(x, file.path(OUTDIR, paste0("degraded_sensor_IDs_final_classes_", label, ".csv")))

  x[]
}

make_recovery_table <- function(sensor_summary, label) {

  if (nrow(sensor_summary) == 0L) {
    rec <- data.table(
      region = character(),
      n_persistent = integer(),
      n_permanent = integer(),
      n_recovered = integer(),
      recovery_probability = numeric(),
      median_onset_year = numeric(),
      median_effective_lifetime_year = numeric()
    )
  } else {
    rec <- sensor_summary[
      has_persistent_degradation == TRUE,
      .(
        n_persistent = .N,
        n_permanent = sum(permanent_degradation, na.rm = TRUE),
        n_recovered = sum(!permanent_degradation, na.rm = TRUE),
        recovery_probability = mean(!permanent_degradation, na.rm = TRUE),
        median_onset_year = median(persistent_onset_year, na.rm = TRUE),
        median_effective_lifetime_year = median(effective_lifetime_year, na.rm = TRUE)
      ),
      by = region
    ]
  }

  fwrite(rec, file.path(OUTDIR, paste0("recovery_summary_by_region_", label, ".csv")))

  rec[]
}

make_summary_plots <- function(final_classes, label) {

  if (nrow(final_classes) == 0L) {
    say("No final classes to plot for: ", label)
    return(invisible(NULL))
  }

  p1 <- ggplot(final_classes, aes(degradation_class)) +
    geom_bar(fill = "grey35") +
    coord_flip() +
    scale_y_continuous(labels = comma) +
    labs(
      title = paste0("Sensor degradation classes: ", label),
      x = NULL,
      y = "Number of sensors"
    ) +
    theme_bw(base_size = 12)

  ggsave(
    file.path(FIGDIR, paste0("degradation_classes_", label, ".jpeg")),
    p1,
    width = 10,
    height = 6,
    dpi = 500
  )

  p2 <- ggplot(
    final_classes[is.finite(effective_lifetime_year)],
    aes(effective_lifetime_year, fill = permanent_degradation)
  ) +
    geom_histogram(bins = 50, alpha = 0.75, position = "identity") +
    scale_x_continuous(labels = number_format(accuracy = 0.5)) +
    scale_y_continuous(labels = comma) +
    labs(
      title = paste0("Effective usable lifetime: ", label),
      x = "Effective lifetime until persistent degradation, or censoring year",
      y = "Number of sensors",
      fill = "Permanent degradation"
    ) +
    theme_bw(base_size = 12)

  ggsave(
    file.path(FIGDIR, paste0("effective_lifetime_histogram_", label, ".jpeg")),
    p2,
    width = 10,
    height = 6,
    dpi = 500
  )

  invisible(NULL)
}

run_permanent_degradation_pipeline <- function(infile, label, require_climate = FALSE) {

  dt <- make_hourly_degradation_dt(infile, require_climate = require_climate, label = label)
  dt <- add_rolling_metrics(dt)
  dt <- assign_rule_states(dt)

  if (isTRUE(SAVE_HOURLY_PROCESSED)) {
    fwrite(dt, file.path(OUTDIR, paste0("hourly_processed_", label, ".csv")))
  }

  hmm_fit <- fit_hmm_if_available(dt, label)

  res_class <- classify_all_sensors(dt, label)
  changepoints <- detect_changepoints_all(dt, label)
  states <- summarize_states(dt, label)

  final_classes <- make_trajectory_classes(
    sensor_summary = res_class$summary,
    changepoints = changepoints,
    state_summary = states,
    label = label
  )

  recovery <- make_recovery_table(res_class$summary, label)

  make_summary_plots(final_classes, label)

  save(
    final_classes,
    changepoints,
    states,
    recovery,
    file = file.path(OUTDIR, paste0("permanent_degradation_compact_outputs_", label, ".RData"))
  )

  rm(dt, res_class, changepoints, states, final_classes, recovery, hmm_fit)
  gc()

  invisible(TRUE)
}

run_permanent_degradation_pipeline(
  infile = global_rdata,
  label = "global_region",
  require_climate = FALSE
)

run_permanent_degradation_pipeline(
  infile = climate_rdata,
  label = "single_zone_climate",
  require_climate = TRUE
)

cat("\nDone. Outputs saved to:\n")
cat(OUTDIR, "\n")

# ============================================================
# GLOBAL MAPS + TABLES OF PERMANENTLY DEGRADED SENSORS
# 1) Likely permanently degraded sensors
# 2) Likely permanently degraded sensors with changepoint
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(scales)
})

OUTDIR <- "/Users/priyanka/Downloads/PurpleAir_EPA/permanent_degradation_hourly"
FIGDIR <- file.path(OUTDIR, "figures")
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

# ---------------------------
# Load final degradation classes
# ---------------------------
classes_file <- file.path(
  OUTDIR,
  "degraded_sensor_IDs_final_classes_global_region.csv"
)

degraded_sensors <- fread(classes_file)
degraded_sensors[, station := as.character(station)]

# ---------------------------
# Filter two permanent-degradation definitions
# ---------------------------

perm_degraded <- degraded_sensors[
  degradation_class %in% c(
    "Likely permanently degraded",
    "Likely permanently degraded with changepoint"
  ) &
    is.finite(lat) &
    is.finite(lon)
]

perm_degraded_changepoint <- degraded_sensors[
  degradation_class == "Likely permanently degraded with changepoint" &
    is.finite(lat) &
    is.finite(lon)
]

# ---------------------------
# Save two separate CSV tables
# ---------------------------

keep_cols <- intersect(
  c(
    "station", "region", "lat", "lon",
    "degradation_class",
    "persistent_onset_hour", "persistent_onset_year",
    "effective_lifetime_hour", "effective_lifetime_year",
    "max_hour", "max_year",
    "mean_flagged", "mean_rel_ab", "mean_abs_ab",
    "max_roll_flag_rate", "max_roll_rel_ab", "max_roll_abs_ab",
    "persistence_score",
    "likely_degradation_changepoint",
    "changepoint_hour", "changepoint_year",
    "mean_before", "mean_after", "mean_shift"
  ),
  names(degraded_sensors)
)

fwrite(
  perm_degraded[, ..keep_cols],
  file.path(OUTDIR, "table_likely_permanently_degraded_sensors.csv")
)

fwrite(
  perm_degraded_changepoint[, ..keep_cols],
  file.path(OUTDIR, "table_likely_permanently_degraded_sensors_with_changepoint.csv")
)

# ---------------------------
# Map boundaries
# ---------------------------

world <- rnaturalearth::ne_countries(
  scale = "medium",
  returnclass = "sf"
)

make_bbox_sf <- function(xmin, xmax, ymin, ymax, region_name, crs = 4326) {
  bb <- st_bbox(
    c(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax),
    crs = st_crs(crs)
  )

  st_sf(
    region = region_name,
    geometry = st_as_sfc(bb)
  )
}

region_boxes <- rbind(
  make_bbox_sf(-125, -66,  24,  50, "CONUS"),
  make_bbox_sf(  68,  90,   6,  36, "India"),
  make_bbox_sf( -82, -35, -56,  13, "South America"),
  make_bbox_sf( -20,  55, -35,  38, "Africa")
)

region_labels <- data.frame(
  label = c("CONUS", "India", "South America", "Africa"),
  x = c(-95.5, 79, -58.5, 17.5),
  y = c(37, 21, -21.5, 1.5)
)

# ---------------------------
# Plot helper
# ---------------------------

make_perm_map <- function(dt, title_text, subtitle_text, outfile_stub) {

  p <- ggplot() +
    geom_sf(
      data = world,
      fill = "grey95",
      color = "grey55",
      linewidth = 0.20
    ) +
    geom_sf(
      data = region_boxes,
      fill = NA,
      color = "black",
      linewidth = 0.45
    ) +
    geom_text(
      data = region_labels,
      aes(x = x, y = y, label = label),
      size = 3,
      fontface = "bold"
    ) +
    geom_point(
      data = dt,
      aes(
        x = lon,
        y = lat,
        color = effective_lifetime_year
      ),
      alpha = 0.75,
      size = 1.3
    ) +
    scale_color_viridis_c(
      option = "plasma",
      name = "Effective lifetime\nuntil persistent\ndegradation (years)",
      labels = label_number(accuracy = 0.1),
      na.value = "grey50"
    ) +
    coord_sf(
      xlim = c(-140, 100),
      ylim = c(-60, 60),
      expand = FALSE
    ) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      plot.title = element_text(face = "bold")
    )

  ggsave(
    file.path(FIGDIR, paste0(outfile_stub, ".png")),
    p,
    width = 13,
    height = 7.5,
    dpi = 500,
    bg = "white"
  )

  ggsave(
    file.path(FIGDIR, paste0(outfile_stub, ".jpeg")),
    p,
    width = 13,
    height = 7.5,
    dpi = 500,
    bg = "white"
  )

  p
}

# ---------------------------
# Produce maps
# ---------------------------

p_perm <- make_perm_map(
  dt = perm_degraded,
  title_text = "Likely permanently degraded PurpleAir sensors",
  subtitle_text = "Sensors classified as likely permanently degraded; color indicates effective lifetime before persistent degradation",
  outfile_stub = "global_map_likely_permanently_degraded_sensors"
)

p_perm_cp <- make_perm_map(
  dt = perm_degraded_changepoint,
  title_text = "Likely permanently degraded PurpleAir sensors with changepoint",
  subtitle_text = "Subset with both persistent permanent degradation and detected degradation changepoint",
  outfile_stub = "global_map_likely_permanently_degraded_sensors_with_changepoint"
)

cat("\nSaved CSV tables:\n")
cat(file.path(OUTDIR, "table_likely_permanently_degraded_sensors.csv"), "\n")
cat(file.path(OUTDIR, "table_likely_permanently_degraded_sensors_with_changepoint.csv"), "\n")

cat("\nSaved maps:\n")
cat(file.path(FIGDIR, "global_map_likely_permanently_degraded_sensors.png"), "\n")
cat(file.path(FIGDIR, "global_map_likely_permanently_degraded_sensors_with_changepoint.png"), "\n")
