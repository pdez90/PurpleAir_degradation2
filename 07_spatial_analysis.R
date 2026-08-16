# ==========================================================================
# 17_spatial_dependence.R
#
# Spatial dependence of sensor degradation: nearby-pair similarity, Moran's I,
# semivariograms, spatial GAMs. => Section S4.
#
# Extracted VERBATIM from purpleair_degradation2.Rmd (chunk #36 lines 20971-21808).
# Code is unmodified; see REPRODUCING.md for paths and run order.
# ==========================================================================

# ---- source chunk #36 (Rmd lines 20971-21808) ------------------------

# ============================================================
# MEMORY-EFFICIENT SPATIAL DEPENDENCE OF SENSOR DEGRADATION
# Corrected for LATITUDE(deg) / LONGITUDE(deg)
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(ggplot2)
  library(mgcv)
  library(scales)
})

global_rdata  <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_final.RData"
climate_rdata <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_with_koppen_timevarying.RData"
transition_csv <- "/Users/priyanka/Downloads/PurpleAir_EPA/station_koppen_transition_summary.csv"

OUTDIR <- "/Users/priyanka/Downloads/PurpleAir_EPA/spatial_degradation_dependence_memory_efficient"
FIGDIR <- file.path(OUTDIR, "figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

BIN_W <- 24L
DROP_H <- 7L * 24L

NEAR_DIST_M <- 5000
K_NEIGHBORS <- 8L
PERM_N <- 300L

MIN_OBS_SENSOR <- 500L
MIN_BINS_SENSOR <- 20L
MIN_STATIONS_GROUP <- 25L

MAX_NEAR_PAIRS_PER_VAR <- 250000L
MAX_RANDOM_PAIRS <- 250000L
MAX_SEMIVARIOGRAM_PAIRS <- 300000L
SAVE_PAIR_SAMPLE <- TRUE

a_col_default <- "pm25_cf_1_a_hourly(ug/m3)_pm25_cf_1_a"
b_col_default <- "pm25_cf_1_b_hourly(ug/m3)_pm25_cf_1_b"
geos_col <- "geoscf_pm25"

say <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste0(..., collapse = ""), "\n")
}

safe_file <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  tolower(x)
}

find_col <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) == 0L) return(NA_character_)
  hit[1]
}

calc_slope <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < MIN_BINS_SENSOR) return(NA_real_)
  vx <- var(x)
  if (!is.finite(vx) || vx <= 0) return(NA_real_)
  cov(x, y) / vx
}

load_outdoor_object <- function(infile) {
  e <- new.env(parent = emptyenv())
  load(infile, envir = e)

  if (exists("outdoor", envir = e)) {
    x <- get("outdoor", envir = e)
  } else if (exists("outdoor_final", envir = e)) {
    x <- get("outdoor_final", envir = e)
  } else {
    stop("Could not find outdoor/outdoor_final in: ", infile)
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

# ------------------------------------------------------------
# 1) Build compact sensor-bin data
# ------------------------------------------------------------

make_sensor_bin_dt <- function(infile, require_climate = FALSE) {
  say("Loading: ", infile)

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
    "gee_latitude",
    "sensor_latitude"
  ))

  lon_col <- find_col(names(outdoor), c(
    "LONGITUDE(deg)",
    "lon", "Lon", "lng", "longitude", "Longitude",
    "gee_longitude",
    "sensor_longitude"
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
    stop("Missing columns: ", paste(miss, collapse = ", "))
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
  dt[, bias_ab := cf1_a - cf1_b]
  dt[, rel_ab := fifelse(
    (cf1_a + cf1_b) > 0,
    2 * abs_ab / (cf1_a + cf1_b),
    NA_real_
  )]

  thr <- dt[
    is.finite(rel_ab),
    .(thr85_rel = quantile(rel_ab, 0.85, na.rm = TRUE, type = 7)),
    by = station
  ]

  dt <- thr[dt, on = "station"]
  dt[, flagged_ab := as.integer(abs_ab > 5 & rel_ab > thr85_rel)]
  dt[is.na(flagged_ab), flagged_ab := 0L]

  dt[, pa_geos_abs := fifelse(is.finite(geos_pm), abs(pm_cf1 - geos_pm), NA_real_)]
  dt[, pa_geos_signed := fifelse(is.finite(geos_pm), pm_cf1 - geos_pm, NA_real_)]

  setorder(dt, station, ts)
  dt[, first_ts := ts[1L], by = station]
  dt[, op_hours := as.numeric(difftime(ts, first_ts, units = "hours"))]
  dt[, op_bin := floor(op_hours / BIN_W) * BIN_W]

  dt <- dt[is.finite(op_bin) & op_bin >= DROP_H]

  sens_bin <- dt[, .(
    region = region[1L],
    lat = median(lat, na.rm = TRUE),
    lon = median(lon, na.rm = TRUE),
    n_obs = .N,
    mean_abs_ab = mean(abs_ab, na.rm = TRUE),
    mean_bias_ab = mean(bias_ab, na.rm = TRUE),
    mean_rel_ab = mean(rel_ab, na.rm = TRUE),
    flagged_rate = mean(flagged_ab, na.rm = TRUE),
    mean_pa_geos_abs = mean(pa_geos_abs, na.rm = TRUE),
    mean_pa_geos_signed = mean(pa_geos_signed, na.rm = TRUE),
    mean_pm = mean(pm_cf1, na.rm = TRUE)
  ), by = .(station, op_bin)]

  rm(dt, thr)
  gc()

  if (require_climate) {
    sens_bin <- merge(sens_bin, single_zone, by = "station", all.x = FALSE)
  } else {
    sens_bin[, `:=`(
      climate_zone = NA_character_,
      climate_desc = NA_character_,
      climate_label = NA_character_
    )]
  }

  sens_bin[]
}

# ------------------------------------------------------------
# 2) Sensor-level slopes
# ------------------------------------------------------------

make_sensor_slopes <- function(sens_bin, label) {
  say("Estimating slopes: ", label)

  slopes <- sens_bin[
    ,
    .(
      region = region[1L],
      climate_zone = climate_zone[1L],
      climate_desc = climate_desc[1L],
      climate_label = climate_label[1L],
      lat = median(lat, na.rm = TRUE),
      lon = median(lon, na.rm = TRUE),
      n_bins = .N,
      n_obs = sum(n_obs, na.rm = TRUE),
      max_op_years = max(op_bin, na.rm = TRUE) / (24 * 365.25),

      slope_abs_ab_per_year = calc_slope(op_bin, mean_abs_ab) * 24 * 365.25,
      slope_bias_ab_per_year = calc_slope(op_bin, mean_bias_ab) * 24 * 365.25,
      slope_rel_ab_per_year = calc_slope(op_bin, mean_rel_ab) * 24 * 365.25,
      slope_flagged_per_year = calc_slope(op_bin, flagged_rate) * 24 * 365.25,
      slope_pa_geos_abs_per_year = calc_slope(op_bin, mean_pa_geos_abs) * 24 * 365.25,
      slope_pa_geos_signed_per_year = calc_slope(op_bin, mean_pa_geos_signed) * 24 * 365.25,
      slope_pm_per_year = calc_slope(op_bin, mean_pm) * 24 * 365.25
    ),
    by = station
  ]

  slopes <- slopes[
    n_obs >= MIN_OBS_SENSOR &
      n_bins >= MIN_BINS_SENSOR &
      is.finite(lat) &
      is.finite(lon)
  ]

  fwrite(slopes, file.path(OUTDIR, paste0("sensor_degradation_slopes_", label, ".csv")))
  slopes[]
}

# ------------------------------------------------------------
# 3) Spatial helper
# ------------------------------------------------------------

prep_spatial <- function(slopes, slope_var) {
  dd <- slopes[is.finite(get(slope_var)) & is.finite(lat) & is.finite(lon)]
  if (nrow(dd) < MIN_STATIONS_GROUP) return(NULL)

  sf_pts <- st_as_sf(dd, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
  sf_m <- st_transform(sf_pts, 5070)
  xy <- st_coordinates(sf_m)

  list(dd = dd, xy = xy, sf_m = sf_m)
}

# ------------------------------------------------------------
# 4) Nearby-pair summary
# ------------------------------------------------------------

nearby_pair_summary <- function(slopes, slope_var, label) {
  obj <- prep_spatial(slopes, slope_var)
  if (is.null(obj)) return(data.table())

  dd <- obj$dd
  xy <- obj$xy
  sf_m <- obj$sf_m

  nb <- st_is_within_distance(sf_m, sf_m, dist = NEAR_DIST_M)

  pair_i <- integer()
  pair_j <- integer()

  for (i in seq_along(nb)) {
    js <- nb[[i]]
    js <- js[js > i]
    if (length(js)) {
      pair_i <- c(pair_i, rep(i, length(js)))
      pair_j <- c(pair_j, js)
    }
  }

  n_pairs_total <- length(pair_i)
  if (n_pairs_total == 0L) return(data.table())

  if (n_pairs_total > MAX_NEAR_PAIRS_PER_VAR) {
    set.seed(123)
    keep <- sample.int(n_pairs_total, MAX_NEAR_PAIRS_PER_VAR)
    pair_i <- pair_i[keep]
    pair_j <- pair_j[keep]
  }

  dist_m <- sqrt(
    (xy[pair_i, 1] - xy[pair_j, 1])^2 +
      (xy[pair_i, 2] - xy[pair_j, 2])^2
  )

  absdiff <- abs(dd[[slope_var]][pair_i] - dd[[slope_var]][pair_j])

  out <- data.table(
    label = label,
    slope_var = slope_var,
    n_stations = nrow(dd),
    n_near_pairs_total = n_pairs_total,
    n_near_pairs_used = length(absdiff),
    mean_near_abs_diff = mean(absdiff, na.rm = TRUE),
    median_near_abs_diff = median(absdiff, na.rm = TRUE),
    mean_near_distance_m = mean(dist_m, na.rm = TRUE)
  )

  set.seed(123)
  n_rand <- min(MAX_RANDOM_PAIRS, length(absdiff), nrow(dd) * 20L)

  random_means <- replicate(PERM_N, {
    ii <- sample.int(nrow(dd), n_rand, replace = TRUE)
    jj <- sample.int(nrow(dd), n_rand, replace = TRUE)
    mean(abs(dd[[slope_var]][ii] - dd[[slope_var]][jj]), na.rm = TRUE)
  })

  out[, `:=`(
    random_mean_abs_diff = mean(random_means, na.rm = TRUE),
    random_sd_abs_diff = sd(random_means, na.rm = TRUE),
    similarity_ratio = mean_near_abs_diff / mean(random_means, na.rm = TRUE),
    p_near_more_similar = mean(random_means <= mean_near_abs_diff, na.rm = TRUE)
  )]

  fwrite(out, file.path(OUTDIR, paste0("nearby_similarity_", slope_var, "_", label, ".csv")))

  if (isTRUE(SAVE_PAIR_SAMPLE)) {
    fwrite(
      data.table(
        station_i = dd$station[pair_i],
        station_j = dd$station[pair_j],
        distance_m = dist_m,
        abs_slope_diff = absdiff,
        slope_var = slope_var,
        label = label
      ),
      file.path(OUTDIR, paste0("nearby_pair_SAMPLE_", slope_var, "_", label, ".csv"))
    )
  }

  rm(nb, pair_i, pair_j, dist_m, absdiff, random_means, obj)
  gc()

  out[]
}

# ------------------------------------------------------------
# 5) Moran's I using k-nearest neighbors
# ------------------------------------------------------------

moran_knn <- function(slopes, slope_var, label, k = K_NEIGHBORS) {
  obj <- prep_spatial(slopes, slope_var)
  if (is.null(obj)) return(data.table())

  dd <- obj$dd
  xy <- obj$xy
  n <- nrow(dd)

  if (n > 8000L) {
    set.seed(123)
    idx <- sample.int(n, 8000L)
    dd <- dd[idx]
    xy <- xy[idx, ]
    n <- nrow(dd)
  }

  dmat <- as.matrix(dist(xy))
  diag(dmat) <- Inf

  nb <- t(apply(dmat, 1, function(x) order(x)[seq_len(min(k, length(x)))]))

  x <- dd[[slope_var]]
  z <- x - mean(x, na.rm = TRUE)
  denominator <- sum(z^2, na.rm = TRUE)
  wij_sum <- length(nb)

  numerator <- 0
  for (i in seq_len(n)) {
    js <- nb[i, ]
    numerator <- numerator + sum(z[i] * z[js], na.rm = TRUE)
  }

  I <- (n / wij_sum) * numerator / denominator

  set.seed(123)
  perm_I <- replicate(PERM_N, {
    zp <- sample(z)
    num <- 0
    for (i in seq_len(n)) {
      js <- nb[i, ]
      num <- num + sum(zp[i] * zp[js], na.rm = TRUE)
    }
    (n / wij_sum) * num / denominator
  })

  out <- data.table(
    label = label,
    slope_var = slope_var,
    k = k,
    n_stations = n,
    moran_I = I,
    perm_mean = mean(perm_I, na.rm = TRUE),
    perm_sd = sd(perm_I, na.rm = TRUE),
    p_greater = mean(perm_I >= I, na.rm = TRUE),
    p_two_sided = mean(abs(perm_I) >= abs(I), na.rm = TRUE)
  )

  fwrite(out, file.path(OUTDIR, paste0("morans_I_", slope_var, "_", label, ".csv")))

  rm(dmat, nb, perm_I, obj)
  gc()

  out[]
}

# ------------------------------------------------------------
# 6) Semivariogram
# ------------------------------------------------------------

make_semivariogram <- function(slopes, slope_var, label,
                               max_dist_km = 500, bin_km = 10) {
  obj <- prep_spatial(slopes, slope_var)
  if (is.null(obj)) return(data.table())

  dd <- obj$dd
  xy <- obj$xy
  n <- nrow(dd)

  set.seed(123)
  np <- min(MAX_SEMIVARIOGRAM_PAIRS, n * 50L)

  i <- sample.int(n, np, replace = TRUE)
  j <- sample.int(n, np, replace = TRUE)

  keep <- i < j
  i <- i[keep]
  j <- j[keep]

  dist_km <- sqrt(
    (xy[i, 1] - xy[j, 1])^2 +
      (xy[i, 2] - xy[j, 2])^2
  ) / 1000

  keep2 <- is.finite(dist_km) & dist_km <= max_dist_km
  i <- i[keep2]
  j <- j[keep2]
  dist_km <- dist_km[keep2]

  if (!length(dist_km)) return(data.table())

  semivar <- 0.5 * (dd[[slope_var]][i] - dd[[slope_var]][j])^2

  semi <- data.table(
    dist_km = dist_km,
    semivariance = semivar
  )

  semi[, dist_bin_km := floor(dist_km / bin_km) * bin_km]

  semi_sum <- semi[, .(
    n_pairs = .N,
    mean_dist_km = mean(dist_km, na.rm = TRUE),
    semivariance = mean(semivariance, na.rm = TRUE),
    median_semivariance = median(semivariance, na.rm = TRUE)
  ), by = dist_bin_km][order(dist_bin_km)]

  semi_sum[, `:=`(
    label = label,
    slope_var = slope_var
  )]

  fwrite(semi_sum, file.path(OUTDIR, paste0("semivariogram_", slope_var, "_", label, ".csv")))

  p <- ggplot(semi_sum[n_pairs >= 20], aes(mean_dist_km, semivariance)) +
    geom_point(size = 1.8) +
    geom_line(linewidth = 0.7) +
    labs(
      title = paste0("Semivariogram of degradation slopes: ", slope_var),
      subtitle = label,
      x = "Distance between sensors (km)",
      y = "Semivariance"
    ) +
    theme_bw(base_size = 12)

  ggsave(
    file.path(FIGDIR, paste0("semivariogram_", slope_var, "_", label, ".jpeg")),
    p,
    width = 8,
    height = 5.5,
    dpi = 500
  )

  rm(i, j, dist_km, semivar, semi, p, obj)
  gc()

  semi_sum[]
}

# ------------------------------------------------------------
# 7) Spatial GAM map
# ------------------------------------------------------------

fit_spatial_gam_slopes <- function(slopes, slope_var, label) {
  dd <- slopes[is.finite(get(slope_var)) & is.finite(lon) & is.finite(lat)]

  if (nrow(dd) < 100L) return(NULL)

  if (nrow(dd) > 10000L) {
    set.seed(123)
    dd <- dd[sample.int(nrow(dd), 10000L)]
  }

  fit <- tryCatch(
    gam(
      as.formula(paste0(slope_var, " ~ s(lon, lat, bs = 'tp')")),
      data = dd,
      method = "REML"
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NULL)

  capture.output(
    summary(fit),
    file = file.path(OUTDIR, paste0("spatial_gam_summary_", slope_var, "_", label, ".txt"))
  )

  p <- ggplot(dd, aes(lon, lat, color = get(slope_var))) +
    borders("world", colour = "grey80", fill = NA, linewidth = 0.2) +
    geom_point(size = 1.2, alpha = 0.75) +
    coord_quickmap() +
    scale_color_viridis_c(option = "plasma") +
    labs(
      title = paste0("Sensor degradation slope map: ", slope_var),
      subtitle = label,
      x = "Longitude",
      y = "Latitude",
      color = "Slope/year"
    ) +
    theme_bw(base_size = 11)

  ggsave(
    file.path(FIGDIR, paste0("map_sensor_slopes_", slope_var, "_", label, ".jpeg")),
    p,
    width = 9,
    height = 6,
    dpi = 500
  )

  rm(p, dd)
  gc()

  invisible(fit)
}

# ------------------------------------------------------------
# 8) Run suite one slope at a time
# ------------------------------------------------------------

run_spatial_suite <- function(slopes, label) {
  slope_vars <- c(
    "slope_abs_ab_per_year",
    "slope_rel_ab_per_year",
    "slope_flagged_per_year",
    "slope_pa_geos_abs_per_year"
  )

  summary_all <- list()

  for (sv in slope_vars) {
    say("Running spatial suite: ", label, " | ", sv)

    near <- nearby_pair_summary(slopes, sv, label)
    moran <- moran_knn(slopes, sv, label)
    semi <- make_semivariogram(slopes, sv, label)
    fit_spatial_gam_slopes(slopes, sv, label)

    summary_all[[sv]] <- rbindlist(list(
      if (nrow(near)) near[, .(
        label, slope_var,
        metric = "nearby_similarity_ratio",
        estimate = similarity_ratio,
        p_value = p_near_more_similar,
        n = n_near_pairs_used
      )],
      if (nrow(moran)) moran[, .(
        label, slope_var,
        metric = "morans_I",
        estimate = moran_I,
        p_value = p_two_sided,
        n = n_stations
      )]
    ), fill = TRUE)

    rm(near, moran, semi)
    gc()
  }

  summary_all <- rbindlist(summary_all, fill = TRUE)

  fwrite(
    summary_all,
    file.path(OUTDIR, paste0("spatial_dependence_summary_", label, ".csv"))
  )

  if (nrow(summary_all) > 0) {
    p <- ggplot(summary_all, aes(slope_var, estimate, color = metric)) +
      geom_hline(yintercept = 1, linetype = 2, color = "grey60") +
      geom_point(size = 2.8, position = position_dodge(width = 0.5)) +
      coord_flip() +
      labs(
        title = paste0("Spatial dependence summary: ", label),
        subtitle = "Similarity ratio < 1 means nearby sensors are more similar than random pairs",
        x = NULL,
        y = "Estimate",
        color = NULL
      ) +
      theme_bw(base_size = 12)

    ggsave(
      file.path(FIGDIR, paste0("spatial_dependence_summary_", label, ".jpeg")),
      p,
      width = 9,
      height = 5.5,
      dpi = 500
    )
  }

  gc()
  invisible(summary_all)
}

# ------------------------------------------------------------
# 9) Global / regional analyses from global_rdata
# ------------------------------------------------------------

say("GLOBAL / REGION ANALYSES")

sens_global <- make_sensor_bin_dt(global_rdata, require_climate = FALSE)
slopes_global <- make_sensor_slopes(sens_global, "global")

rm(sens_global)
gc()

run_spatial_suite(slopes_global, "global")

for (rr in sort(unique(na.omit(slopes_global$region)))) {
  sub <- slopes_global[region == rr]

  if (nrow(sub) >= MIN_STATIONS_GROUP) {
    lab <- paste0("region_", safe_file(rr))
    run_spatial_suite(sub, lab)
  }

  rm(sub)
  gc()
}

# ------------------------------------------------------------
# 10) Climate analyses from climate_rdata
# ------------------------------------------------------------

say("CLIMATE ANALYSES")

sens_climate <- make_sensor_bin_dt(climate_rdata, require_climate = TRUE)
slopes_climate <- make_sensor_slopes(sens_climate, "single_zone_climate")

rm(sens_climate)
gc()

run_spatial_suite(slopes_climate, "single_zone_climate")

for (cc in sort(unique(na.omit(slopes_climate$climate_label)))) {
  sub <- slopes_climate[climate_label == cc]

  if (nrow(sub) >= MIN_STATIONS_GROUP) {
    lab <- paste0("climate_", safe_file(cc))
    run_spatial_suite(sub, lab)
  }

  rm(sub)
  gc()
}

slopes_climate[, region_climate := paste(region, climate_label, sep = " | ")]

for (rc in sort(unique(na.omit(slopes_climate$region_climate)))) {
  sub <- slopes_climate[region_climate == rc]

  if (nrow(sub) >= MIN_STATIONS_GROUP) {
    lab <- paste0("region_x_climate_", safe_file(rc))
    run_spatial_suite(sub, lab)
  }

  rm(sub)
  gc()
}

# ------------------------------------------------------------
# 11) Spatial dependence beyond climate
# ------------------------------------------------------------

say("Testing Moran's I after climate adjustment")

resid_results <- list()

for (sv in c(
  "slope_abs_ab_per_year",
  "slope_rel_ab_per_year",
  "slope_flagged_per_year",
  "slope_pa_geos_abs_per_year"
)) {
  dd <- slopes_climate[
    is.finite(get(sv)) &
      !is.na(climate_label) &
      is.finite(lat) &
      is.finite(lon)
  ]

  if (nrow(dd) >= 100L && uniqueN(dd$climate_label) >= 2L) {
    fit <- lm(as.formula(paste0(sv, " ~ climate_label")), data = dd)
    dd[, resid_slope := resid(fit)]
    setnames(dd, "resid_slope", paste0(sv, "_climate_resid"))

    resid_var <- paste0(sv, "_climate_resid")

    resid_results[[sv]] <- moran_knn(
      dd,
      resid_var,
      paste0("climate_adjusted_", sv)
    )

    rm(fit, dd)
    gc()
  }
}

resid_results <- rbindlist(resid_results, fill = TRUE)
fwrite(resid_results, file.path(OUTDIR, "morans_I_after_climate_adjustment.csv"))

cat("\nDone. Outputs saved to:\n")
cat(OUTDIR, "\n")
