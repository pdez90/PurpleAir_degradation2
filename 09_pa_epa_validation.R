# ==========================================================================
# 09_pa_epa_validation.R
#
# PurpleAir vs EPA regulatory monitors: join to AQS monitors within 100 m (CONUS), then
# calibration models, age-dependent calibration drift, and climate heterogeneity.
# => Section S3; Table S3.1; Figures S3.1-S3.4.
# (Cross-validated performance: sensitivity/S07_calibration_cv.R => Table S3.2.)
#
# Extracted VERBATIM from purpleair_degradation2.Rmd (chunks #34 lines 18955-19630, #35 lines 19633-20968).
# Code is unmodified; see REPRODUCING.md for paths and run order.
# ==========================================================================

# ---- source chunk #34 (Rmd lines 18955-19630) ------------------------

# ============================================================
# JOIN PURPLEAIR outdoor_final.RData TO EPA AQS PM2.5 MONITORS
# within 100 m, USA/CONUS only
#
# Updates:
#  1) Retains PurpleAir cf_1 channel A and B when available.
#  2) Retains PurpleAir temperature and RH.
#  3) Uses exact outdoor_final.RData T/RH columns:
#       temperature_hourly(C)_temperature
#       humidity_hourly(%)_humidity
#  4) Calculates operational age from each sensor's earliest
#     PurpleAir observation BEFORE EPA matching.
#  5) Avoids unstable relative-difference values when EPA or PA <= 0.
#  6) Removes duplicated summary section.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(lubridate)
})

# ------------------------------------------------------------
# 0) PATHS
# ------------------------------------------------------------

pa_rdata <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_final.RData"
aqs_dir  <- "/Users/priyanka/Downloads/AQS"

outdir <- "/Users/priyanka/Downloads/PurpleAir_EPA/pa_epa_100m_join"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

aqs_monitor_file <- file.path(aqs_dir, "aqs_monitors.csv")
aqs_years <- 2018:2025
join_radius_m <- 100

say <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste0(..., collapse = ""), "\n")
}

# ------------------------------------------------------------
# 1) HELPERS
# ------------------------------------------------------------

load_pa <- function(infile) {
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

find_col <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) == 0L) return(NA_character_)
  hit[1]
}

clean_names_simple <- function(x) {
  x2 <- tolower(x)
  x2 <- gsub("[^a-z0-9]+", "_", x2)
  x2 <- gsub("^_|_$", "", x2)
  x2
}

pick_col <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) == 0L) return(NA_character_)
  hit[1]
}

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(NA_real_)
  cor(x[ok], y[ok])
}

# ------------------------------------------------------------
# 2) LOAD PURPLEAIR DATA
# ------------------------------------------------------------

say("Loading PurpleAir data")
pa_raw <- load_pa(pa_rdata)

lat_col <- find_col(names(pa_raw), c(
  "LATITUDE(deg)", "lat", "Lat", "latitude", "Latitude",
  "gee_latitude", "sensor_latitude"
))

lon_col <- find_col(names(pa_raw), c(
  "LONGITUDE(deg)", "lon", "Lon", "lng", "longitude", "Longitude",
  "gee_longitude", "sensor_longitude"
))

time_col <- find_col(names(pa_raw), c(
  "datetime_utc", "ts_hr", "date", "datetime", "timestamp", "time"
))

pm_avg_col <- find_col(names(pa_raw), c(
  "pm25_cf1_avg",
  "PM2.5_CF_1_ug/m3",
  "pa_pm25_cf1",
  "pm_cf1"
))

pm_a_col <- find_col(names(pa_raw), c(
  "pm25_cf_1_a_hourly(ug/m3)_pm25_cf_1_a",
  "PM2.5_CF_1_ug/m3_A",
  "cf1_a",
  "pm_a"
))

pm_b_col <- find_col(names(pa_raw), c(
  "pm25_cf_1_b_hourly(ug/m3)_pm25_cf_1_b",
  "PM2.5_CF_1_ug/m3_B",
  "cf1_b",
  "pm_b"
))

temp_col <- "temperature_hourly(C)_temperature"
rh_col   <- "humidity_hourly(%)_humidity"

needed_core <- c("station", "region", time_col, lat_col, lon_col)

if (any(is.na(needed_core)) || (is.na(pm_avg_col) && (is.na(pm_a_col) || is.na(pm_b_col)))) {
  stop(
    "Could not identify required PA columns.\n\n",
    "Identified:\n",
    "time_col = ", time_col, "\n",
    "lat_col = ", lat_col, "\n",
    "lon_col = ", lon_col, "\n",
    "pm_avg_col = ", pm_avg_col, "\n",
    "pm_a_col = ", pm_a_col, "\n",
    "pm_b_col = ", pm_b_col, "\n\n",
    "Available columns:\n",
    paste(names(pa_raw), collapse = ", ")
  )
}

if (!temp_col %in% names(pa_raw)) {
  stop("Missing PurpleAir temperature column: ", temp_col)
}

if (!rh_col %in% names(pa_raw)) {
  stop("Missing PurpleAir RH column: ", rh_col)
}

say("Using coordinate columns: ", lat_col, ", ", lon_col)
say("Using time column: ", time_col)
say("Using PurpleAir temperature column: ", temp_col)
say("Using PurpleAir RH column: ", rh_col)

if (!is.na(pm_avg_col)) say("Using PurpleAir CF1 average column: ", pm_avg_col)
if (!is.na(pm_a_col))   say("Retaining PurpleAir channel A column: ", pm_a_col)
if (!is.na(pm_b_col))   say("Retaining PurpleAir channel B column: ", pm_b_col)

# ------------------------------------------------------------
# 3) STANDARDIZE PURPLEAIR DATA
#    Important: retain A/B even when pm25_cf1_avg exists.
# ------------------------------------------------------------

pa <- pa_raw[, .(
  station = as.character(station),
  region = as.character(region),
  datetime_utc = as.POSIXct(get(time_col), tz = "UTC"),
  lat = as.numeric(get(lat_col)),
  lon = as.numeric(get(lon_col)),

  pa_pm25_a = if (!is.na(pm_a_col)) as.numeric(get(pm_a_col)) else NA_real_,
  pa_pm25_b = if (!is.na(pm_b_col)) as.numeric(get(pm_b_col)) else NA_real_,

  pa_pm25_cf1 = if (!is.na(pm_avg_col)) {
    as.numeric(get(pm_avg_col))
  } else {
    NA_real_
  },

  pa_temperature = as.numeric(get(temp_col)),
  pa_rh = as.numeric(get(rh_col))
)]

# If no precomputed average exists, calculate it from A/B.
pa[
  !is.finite(pa_pm25_cf1) & is.finite(pa_pm25_a) & is.finite(pa_pm25_b),
  pa_pm25_cf1 := (pa_pm25_a + pa_pm25_b) / 2
]

rm(pa_raw)
gc()

# ------------------------------------------------------------
# 4) SUBSET USA / CONUS
# ------------------------------------------------------------

pa <- pa[
  region %in% c("USA", "US", "CONUS", "United States", "North America") |
    grepl("usa|conus|united", region, ignore.case = TRUE)
]

pa <- pa[
  !is.na(station) &
    is.finite(lat) &
    is.finite(lon) &
    !is.na(datetime_utc) &
    is.finite(pa_pm25_cf1)
]

pa[, datetime_utc := floor_date(datetime_utc, unit = "hour")]

say("Rows after USA/CONUS subset: ", format(nrow(pa), big.mark = ","))
say("Stations after USA/CONUS subset: ", format(uniqueN(pa$station), big.mark = ","))
say("Non-missing PA channel A rows: ", format(sum(is.finite(pa$pa_pm25_a)), big.mark = ","))
say("Non-missing PA channel B rows: ", format(sum(is.finite(pa$pa_pm25_b)), big.mark = ","))
say("Non-missing PA temperature rows: ", format(sum(is.finite(pa$pa_temperature)), big.mark = ","))
say("Non-missing PA RH rows: ", format(sum(is.finite(pa$pa_rh)), big.mark = ","))

# Operational start from earliest PurpleAir observation BEFORE EPA matching
pa_sensor_start <- pa[, .(
  first_pa_observation = min(datetime_utc, na.rm = TRUE),
  last_pa_observation  = max(datetime_utc, na.rm = TRUE),
  n_total_pa_hours = .N
), by = station]

pa <- merge(pa, pa_sensor_start, by = "station", all.x = TRUE)

pa[, `:=`(
  op_hours_since_first_pa = as.numeric(difftime(datetime_utc, first_pa_observation, units = "hours")),
  op_days_since_first_pa  = as.numeric(difftime(datetime_utc, first_pa_observation, units = "days")),
  op_years_since_first_pa = as.numeric(difftime(datetime_utc, first_pa_observation, units = "days")) / 365.25
)]

pa_station <- pa[, .(
  lat = median(lat, na.rm = TRUE),
  lon = median(lon, na.rm = TRUE),
  region = region[1],
  n_pa_hours = .N,
  first_pa_observation = first_pa_observation[1],
  last_pa_observation = last_pa_observation[1],
  n_total_pa_hours = n_total_pa_hours[1]
), by = station]

# ------------------------------------------------------------
# 5) LOAD EPA AQS MONITOR METADATA
# ------------------------------------------------------------

say("Loading AQS monitor metadata")
aqs_mon <- fread(aqs_monitor_file)

aqs_mon <- aqs_mon[`Parameter Name` == "PM2.5 - Local Conditions"]

aqs_mon[, `:=`(
  state_code = sprintf("%02d", as.integer(`State Code`)),
  county_code = sprintf("%03d", as.integer(`County Code`)),
  site_num = sprintf("%04d", as.integer(`Site Number`)),
  aqs_lat = as.numeric(Latitude),
  aqs_lon = as.numeric(Longitude)
)]

aqs_mon <- aqs_mon[is.finite(aqs_lat) & is.finite(aqs_lon)]
aqs_mon[, aqs_id := paste(state_code, county_code, site_num, sep = "-")]

aqs_mon_small <- unique(aqs_mon[, .(
  aqs_id,
  state_code,
  county_code,
  site_num,
  aqs_lat,
  aqs_lon,
  state_name = `State Name`,
  county_name = `County Name`,
  parameter_name = `Parameter Name`
)])

say("AQS PM2.5 monitors: ", format(nrow(aqs_mon_small), big.mark = ","))

# ------------------------------------------------------------
# 6) SPATIAL MATCH: PA STATIONS TO NEAREST AQS MONITOR WITHIN 100 m
# ------------------------------------------------------------

say("Spatially matching PA stations to AQS monitors within ", join_radius_m, " m")

pa_sf <- st_as_sf(pa_station, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
aqs_sf <- st_as_sf(aqs_mon_small, coords = c("aqs_lon", "aqs_lat"), crs = 4326, remove = FALSE)

pa_sf_m <- st_transform(pa_sf, 5070)
aqs_sf_m <- st_transform(aqs_sf, 5070)

nearest_idx <- st_nearest_feature(pa_sf_m, aqs_sf_m)

nearest_dist <- as.numeric(st_distance(
  pa_sf_m,
  aqs_sf_m[nearest_idx, ],
  by_element = TRUE
))

station_match <- cbind(
  st_drop_geometry(pa_sf_m),
  st_drop_geometry(aqs_sf_m[nearest_idx, ])
)

station_match <- as.data.table(station_match)
station_match[, distance_m := nearest_dist]
station_match <- station_match[distance_m <= join_radius_m]

station_match <- station_match[, .(
  station,
  region,
  pa_lat = lat,
  pa_lon = lon,
  n_pa_hours,
  first_pa_observation,
  last_pa_observation,
  n_total_pa_hours,
  aqs_id,
  state_code,
  county_code,
  site_num,
  state_name,
  county_name,
  aqs_lat,
  aqs_lon,
  distance_m
)]

fwrite(station_match, file.path(outdir, "pa_station_aqs_monitor_matches_within_100m.csv"))

say("Matched PA stations: ", format(uniqueN(station_match$station), big.mark = ","))
say("Matched AQS monitors: ", format(uniqueN(station_match$aqs_id), big.mark = ","))

if (nrow(station_match) == 0L) {
  stop("No PurpleAir stations matched to AQS PM2.5 monitors within 100 m.")
}

# ------------------------------------------------------------
# 7) KEEP MATCHED PA HOURLY DATA
# ------------------------------------------------------------

matched_stations <- unique(station_match$station)
pa_matched <- pa[station %in% matched_stations]

pa_matched <- merge(
  pa_matched,
  station_match[, .(
    station,
    aqs_id,
    state_code,
    county_code,
    site_num,
    state_name,
    county_name,
    aqs_lat,
    aqs_lon,
    distance_m
  )],
  by = "station",
  all.x = FALSE,
  allow.cartesian = TRUE
)

setkey(pa_matched, aqs_id, datetime_utc)

rm(pa, pa_station, pa_sf, aqs_sf, pa_sf_m, aqs_sf_m)
gc()

# ------------------------------------------------------------
# 8) LOAD AQS HOURLY PM2.5 YEAR-BY-YEAR AND JOIN
# ------------------------------------------------------------

read_aqs_year <- function(year) {
  f <- file.path(aqs_dir, paste0("hourly_88101_", year, ".csv"))

  if (!file.exists(f)) {
    warning("Missing AQS file: ", f)
    return(data.table())
  }

  say("Reading AQS hourly file: ", basename(f))

  x <- fread(f)
  setnames(x, clean_names_simple(names(x)))

  state_col  <- pick_col(names(x), c("state_code", "state"))
  county_col <- pick_col(names(x), c("county_code", "county"))
  site_col   <- pick_col(names(x), c("site_num", "site_number", "site_code"))
  param_col  <- pick_col(names(x), c("parameter_name", "parameter"))
  date_col   <- pick_col(names(x), c("date_gmt", "date_local"))
  time_col   <- pick_col(names(x), c("time_gmt", "time_local"))
  meas_col   <- pick_col(names(x), c("sample_measurement", "sample_value", "arithmetic_mean"))
  method_col <- pick_col(names(x), c("method_name", "method"))
  poc_col    <- pick_col(names(x), c("poc"))

  required <- c(state_col, county_col, site_col, date_col, time_col, meas_col)

  if (any(is.na(required))) {
    stop(
      "Could not identify required AQS columns for year ", year, ".\n",
      "Available standardized columns:\n",
      paste(names(x), collapse = ", ")
    )
  }

  if (!is.na(param_col)) {
    x <- x[
      get(param_col) == "PM2.5 - Local Conditions" |
        is.na(get(param_col))
    ]
  }

  x[, `:=`(
    state_code = sprintf("%02d", as.integer(get(state_col))),
    county_code = sprintf("%03d", as.integer(get(county_col))),
    site_num = sprintf("%04d", as.integer(get(site_col))),
    epa_pm25 = as.numeric(get(meas_col)),
    datetime_utc = ymd_hm(paste(get(date_col), get(time_col)), tz = "UTC")
  )]

  x[, datetime_utc := floor_date(datetime_utc, unit = "hour")]
  x[, aqs_id := paste(state_code, county_code, site_num, sep = "-")]

  x <- x[
    aqs_id %in% unique(station_match$aqs_id) &
      !is.na(datetime_utc) &
      is.finite(epa_pm25)
  ]

  x <- x[, .(
    epa_pm25 = mean(epa_pm25, na.rm = TRUE),
    epa_n_records = .N,
    epa_method = if (!is.na(method_col)) {
      paste(unique(na.omit(get(method_col))), collapse = " | ")
    } else {
      NA_character_
    },
    epa_poc = if (!is.na(poc_col)) {
      paste(unique(na.omit(get(poc_col))), collapse = " | ")
    } else {
      NA_character_
    }
  ), by = .(aqs_id, datetime_utc)]

  setkey(x, aqs_id, datetime_utc)
  x[]
}

out_join_csv <- file.path(outdir, "pa_epa_hourly_joined_within_100m_2018_2025.csv")
if (file.exists(out_join_csv)) file.remove(out_join_csv)

join_counts <- list()

for (yy in aqs_years) {
  epa_y <- read_aqs_year(yy)

  if (nrow(epa_y) == 0L) next

  pa_y <- pa_matched[
    year(datetime_utc) == yy,
    .(
      station,
      region,
      datetime_utc,

      first_pa_observation,
      last_pa_observation,
      n_total_pa_hours,
      op_hours_since_first_pa,
      op_days_since_first_pa,
      op_years_since_first_pa,

      pa_lat = lat,
      pa_lon = lon,

      pa_pm25_a,
      pa_pm25_b,
      pa_pm25_cf1,
      pa_temperature,
      pa_rh,

      aqs_id,
      state_code,
      county_code,
      site_num,
      state_name,
      county_name,
      aqs_lat,
      aqs_lon,
      distance_m
    )
  ]

  setkey(pa_y, aqs_id, datetime_utc)

  joined_y <- epa_y[pa_y, on = .(aqs_id, datetime_utc), nomatch = 0L]

  if (nrow(joined_y) > 0L) {
    joined_y[, `:=`(
      year = yy,
      pa_epa_diff = pa_pm25_cf1 - epa_pm25,
      pa_epa_absdiff = abs(pa_pm25_cf1 - epa_pm25),
      pa_epa_rel_diff = fifelse(
        is.finite(pa_pm25_cf1) & is.finite(epa_pm25) &
          pa_pm25_cf1 > 0 & epa_pm25 > 0,
        2 * abs(pa_pm25_cf1 - epa_pm25) / (pa_pm25_cf1 + epa_pm25),
        NA_real_
      )
    )]

    fwrite(joined_y, out_join_csv, append = file.exists(out_join_csv))
  }

  join_counts[[as.character(yy)]] <- data.table(
    year = yy,
    n_joined_rows = nrow(joined_y),
    n_pa_stations = uniqueN(joined_y$station),
    n_aqs_monitors = uniqueN(joined_y$aqs_id),
    n_nonmissing_pa_pm25_a = sum(is.finite(joined_y$pa_pm25_a)),
    n_nonmissing_pa_pm25_b = sum(is.finite(joined_y$pa_pm25_b)),
    n_nonmissing_pa_temperature = sum(is.finite(joined_y$pa_temperature)),
    n_nonmissing_pa_rh = sum(is.finite(joined_y$pa_rh))
  )

  say(
    "Year ", yy,
    ": joined rows = ", format(nrow(joined_y), big.mark = ","),
    "; PA stations = ", format(uniqueN(joined_y$station), big.mark = ","),
    "; AQS monitors = ", format(uniqueN(joined_y$aqs_id), big.mark = ","),
    "; non-missing A = ", format(sum(is.finite(joined_y$pa_pm25_a)), big.mark = ","),
    "; non-missing B = ", format(sum(is.finite(joined_y$pa_pm25_b)), big.mark = ","),
    "; non-missing T = ", format(sum(is.finite(joined_y$pa_temperature)), big.mark = ","),
    "; non-missing RH = ", format(sum(is.finite(joined_y$pa_rh)), big.mark = ",")
  )

  rm(epa_y, pa_y, joined_y)
  gc()
}

join_counts <- rbindlist(join_counts, fill = TRUE)
fwrite(join_counts, file.path(outdir, "pa_epa_join_counts_by_year.csv"))

# ------------------------------------------------------------
# 9) SAVE SUMMARY + RDATA
# ------------------------------------------------------------

if (file.exists(out_join_csv)) {
  joined <- fread(out_join_csv)

  summary_overall <- joined[, .(
    n_rows = .N,
    n_pa_stations = uniqueN(station),
    n_aqs_monitors = uniqueN(aqs_id),

    min_datetime_utc = min(datetime_utc, na.rm = TRUE),
    max_datetime_utc = max(datetime_utc, na.rm = TRUE),
    min_op_years_since_first_pa = min(op_years_since_first_pa, na.rm = TRUE),
    max_op_years_since_first_pa = max(op_years_since_first_pa, na.rm = TRUE),

    mean_pa = mean(pa_pm25_cf1, na.rm = TRUE),
    mean_epa = mean(epa_pm25, na.rm = TRUE),
    mean_bias_pa_minus_epa = mean(pa_epa_diff, na.rm = TRUE),
    median_bias_pa_minus_epa = median(pa_epa_diff, na.rm = TRUE),
    mean_absdiff = mean(pa_epa_absdiff, na.rm = TRUE),
    mean_relative_difference = mean(pa_epa_rel_diff, na.rm = TRUE),
    correlation_pa_epa = safe_cor(pa_pm25_cf1, epa_pm25),

    mean_pa_pm25_a = mean(pa_pm25_a, na.rm = TRUE),
    mean_pa_pm25_b = mean(pa_pm25_b, na.rm = TRUE),
    n_nonmissing_pa_pm25_a = sum(is.finite(pa_pm25_a)),
    n_nonmissing_pa_pm25_b = sum(is.finite(pa_pm25_b)),

    mean_pa_temperature = mean(pa_temperature, na.rm = TRUE),
    median_pa_temperature = median(pa_temperature, na.rm = TRUE),
    n_nonmissing_pa_temperature = sum(is.finite(pa_temperature)),

    mean_pa_rh = mean(pa_rh, na.rm = TRUE),
    median_pa_rh = median(pa_rh, na.rm = TRUE),
    n_nonmissing_pa_rh = sum(is.finite(pa_rh))
  )]

  summary_by_station <- joined[, .(
    n_hours = .N,
    aqs_id = aqs_id[1],
    state_name = state_name[1],
    county_name = county_name[1],
    distance_m = distance_m[1],

    first_pa_observation = first_pa_observation[1],
    last_pa_observation = last_pa_observation[1],
    first_joined_observation = min(datetime_utc, na.rm = TRUE),
    last_joined_observation = max(datetime_utc, na.rm = TRUE),
    min_op_years_since_first_pa = min(op_years_since_first_pa, na.rm = TRUE),
    max_op_years_since_first_pa = max(op_years_since_first_pa, na.rm = TRUE),

    mean_pa = mean(pa_pm25_cf1, na.rm = TRUE),
    mean_epa = mean(epa_pm25, na.rm = TRUE),
    mean_bias_pa_minus_epa = mean(pa_epa_diff, na.rm = TRUE),
    mean_absdiff = mean(pa_epa_absdiff, na.rm = TRUE),
    mean_relative_difference = mean(pa_epa_rel_diff, na.rm = TRUE),
    correlation_pa_epa = safe_cor(pa_pm25_cf1, epa_pm25),

    mean_pa_pm25_a = mean(pa_pm25_a, na.rm = TRUE),
    mean_pa_pm25_b = mean(pa_pm25_b, na.rm = TRUE),
    n_nonmissing_pa_pm25_a = sum(is.finite(pa_pm25_a)),
    n_nonmissing_pa_pm25_b = sum(is.finite(pa_pm25_b)),

    mean_pa_temperature = mean(pa_temperature, na.rm = TRUE),
    median_pa_temperature = median(pa_temperature, na.rm = TRUE),
    n_nonmissing_pa_temperature = sum(is.finite(pa_temperature)),

    mean_pa_rh = mean(pa_rh, na.rm = TRUE),
    median_pa_rh = median(pa_rh, na.rm = TRUE),
    n_nonmissing_pa_rh = sum(is.finite(pa_rh))
  ), by = station]

  summary_by_monitor <- joined[, .(
    n_hours = .N,
    n_pa_stations = uniqueN(station),

    mean_pa = mean(pa_pm25_cf1, na.rm = TRUE),
    mean_epa = mean(epa_pm25, na.rm = TRUE),
    mean_bias_pa_minus_epa = mean(pa_epa_diff, na.rm = TRUE),
    mean_absdiff = mean(pa_epa_absdiff, na.rm = TRUE),
    mean_relative_difference = mean(pa_epa_rel_diff, na.rm = TRUE),
    correlation_pa_epa = safe_cor(pa_pm25_cf1, epa_pm25),

    mean_pa_pm25_a = mean(pa_pm25_a, na.rm = TRUE),
    mean_pa_pm25_b = mean(pa_pm25_b, na.rm = TRUE),
    n_nonmissing_pa_pm25_a = sum(is.finite(pa_pm25_a)),
    n_nonmissing_pa_pm25_b = sum(is.finite(pa_pm25_b)),

    mean_pa_temperature = mean(pa_temperature, na.rm = TRUE),
    median_pa_temperature = median(pa_temperature, na.rm = TRUE),
    n_nonmissing_pa_temperature = sum(is.finite(pa_temperature)),

    mean_pa_rh = mean(pa_rh, na.rm = TRUE),
    median_pa_rh = median(pa_rh, na.rm = TRUE),
    n_nonmissing_pa_rh = sum(is.finite(pa_rh))
  ), by = .(aqs_id, state_name, county_name)]

  fwrite(summary_overall, file.path(outdir, "pa_epa_summary_overall.csv"))
  fwrite(summary_by_station, file.path(outdir, "pa_epa_summary_by_station.csv"))
  fwrite(summary_by_monitor, file.path(outdir, "pa_epa_summary_by_monitor.csv"))

  save(
    joined,
    station_match,
    join_counts,
    summary_overall,
    summary_by_station,
    summary_by_monitor,
    file = file.path(outdir, "pa_epa_joined_within_100m_2018_2025.RData")
  )

  say("Overall correlation PA CF1 vs EPA PM2.5: ", round(summary_overall$correlation_pa_epa, 3))
  say("Joined rows: ", format(summary_overall$n_rows, big.mark = ","))
  say("Joined PA stations: ", format(summary_overall$n_pa_stations, big.mark = ","))
  say("Joined AQS monitors: ", format(summary_overall$n_aqs_monitors, big.mark = ","))
  say("Non-missing joined PA channel A: ", format(summary_overall$n_nonmissing_pa_pm25_a, big.mark = ","))
  say("Non-missing joined PA channel B: ", format(summary_overall$n_nonmissing_pa_pm25_b, big.mark = ","))
  say("Non-missing joined PA temperature: ", format(summary_overall$n_nonmissing_pa_temperature, big.mark = ","))
  say("Non-missing joined PA RH: ", format(summary_overall$n_nonmissing_pa_rh, big.mark = ","))
} else {
  warning("No joined PA-EPA CSV was created. Check station matches and AQS files.")
}

cat("\nDone. Outputs saved to:\n")
cat(outdir, "\n")

# ---- source chunk #35 (Rmd lines 19633-20968) ------------------------

# ============================================================
# PA–EPA DEGRADATION + CALIBRATION + KÖPPEN HETEROGENEITY
# Cluster-robust SEs + Word-ready tables
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(mgcv)
  library(scales)
  library(cowplot)
  library(fixest)
  library(splines)
  library(officer)
  library(flextable)
})

# ------------------------------------------------------------
# 0) PATHS
# ------------------------------------------------------------

global_rdata  <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_final.RData"
climate_rdata <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_with_koppen_timevarying.RData"

join_rdata <- "/Users/priyanka/Downloads/PurpleAir_EPA/pa_epa_100m_join/pa_epa_joined_within_100m_2018_2025.RData"
transition_csv <- "/Users/priyanka/Downloads/PurpleAir_EPA/station_koppen_transition_summary.csv"

outdir <- "/Users/priyanka/Downloads/PurpleAir_EPA/pa_epa_degradation_analysis_with_koppen"
figdir <- file.path(outdir, "figures")
tabledir <- file.path(outdir, "word_tables")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
dir.create(tabledir, showWarnings = FALSE, recursive = TRUE)

MIN_N_PANEL <- 100
MIN_STATIONS_PANEL <- 5
AGE_BIN_HOURS <- 168
MAX_MODEL_ROWS <- 2500000L

say <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")),
      paste0(..., collapse = ""), "\n")
}

fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), NA_character_, sprintf(paste0("%.", digits, "f"), x))
}

fmt_ci <- function(est, lwr, upr, digits = 3) {
  paste0(fmt_num(est, digits), " (", fmt_num(lwr, digits), ", ", fmt_num(upr, digits), ")")
}

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(NA_real_)
  cor(x[ok], y[ok])
}

calc_nerr <- function(pred, ref) {
  fifelse(
    is.finite(pred) & is.finite(ref) & ref > 0,
    abs(pred - ref) / pmax(ref, 1),
    NA_real_
  )
}

calc_err <- function(pred, ref) abs(pred - ref)

first_nonmissing <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA)
  x[1]
}

# ------------------------------------------------------------
# ROBUST WORD TABLE HELPER
# ------------------------------------------------------------

make_word_table <- function(dt, title, outfile, note = NULL) {
  if (is.null(dt) || nrow(dt) == 0L) {
    warning("No rows available for Word table: ", title)
    return(invisible(NULL))
  }

  tab <- copy(dt)

  if (all(c("beta_per_year", "ci_lwr_per_year", "ci_upr_per_year") %in% names(tab))) {
    tab[, `Estimate (95% CI)` := fmt_ci(
      beta_per_year,
      ci_lwr_per_year,
      ci_upr_per_year,
      4
    )]
  }

  if ("p_value" %in% names(tab)) {
    tab[, `P value` := fmt_p(p_value)]
  }

  preferred <- c(
    "group_type", "group", "model", "outcome",
    "Estimate (95% CI)", "P value",
    "n_stations", "n_sensor_bins", "n_rows", "n_obs",
    "min_year", "max_year",
    "state_name", "climate_zone", "climate_desc", "climate_label",
    "n_sensors", "n_aqs_monitors",
    "median_hours", "median_duration_days", "median_operational_year",
    "median_age_te_norm_error",
    "correlation", "rmse", "mean_abs_error", "mean_norm_error",
    "pred"
  )

  keep <- intersect(preferred, names(tab))
  if (length(keep) == 0L) keep <- names(tab)

  tab <- tab[, ..keep]

  rename_map <- c(
    group_type = "Group type",
    group = "Group",
    model = "Model",
    outcome = "Outcome",
    n_stations = "N stations",
    n_sensor_bins = "N sensor-bins",
    n_rows = "N rows",
    n_obs = "N observations",
    min_year = "Minimum operational year",
    max_year = "Maximum operational year",
    state_name = "State",
    climate_zone = "Climate zone",
    climate_desc = "Climate description",
    climate_label = "Climate label",
    n_sensors = "N sensors",
    n_aqs_monitors = "N AQS monitors",
    median_hours = "Median matched hours",
    median_duration_days = "Median duration days",
    median_operational_year = "Median operational years",
    median_age_te_norm_error = "Median age-tensor normalized error",
    correlation = "Correlation",
    rmse = "RMSE",
    mean_abs_error = "Mean absolute error",
    mean_norm_error = "Mean normalized error",
    pred = "Prediction variable"
  )

  old <- intersect(names(rename_map), names(tab))
  if (length(old) > 0L) {
    setnames(tab, old = old, new = unname(rename_map[old]))
  }

  ft <- flextable(tab)
  ft <- autofit(ft)
  ft <- fontsize(ft, size = 8, part = "all")
  ft <- bold(ft, part = "header")
  ft <- align(ft, align = "center", part = "header")
  ft <- set_table_properties(ft, width = 1, layout = "autofit")

  doc <- read_docx()
  doc <- body_add_par(doc, title, style = "heading 1")

  if (!is.null(note)) {
    doc <- body_add_par(doc, note, style = "Normal")
  }

  doc <- body_add_flextable(doc, ft)
  print(doc, target = outfile)

  invisible(NULL)
}

# ------------------------------------------------------------
# 1) LOAD JOINED PA–EPA DATA
# ------------------------------------------------------------

say("Loading joined PA–EPA dataset")

load(join_rdata)
setDT(joined)

joined <- joined[
  is.finite(pa_pm25_cf1) &
    is.finite(epa_pm25) &
    epa_pm25 >= 0 &
    pa_pm25_cf1 >= 0 &
    pa_pm25_cf1 <= 1500
]

joined[, station := as.character(station)]
joined[, datetime_utc := as.POSIXct(datetime_utc, tz = "UTC")]

if ("op_hours_since_first_pa" %in% names(joined)) {
  joined[, hour := as.numeric(op_hours_since_first_pa)]
  joined[, year_op := as.numeric(op_years_since_first_pa)]
} else {
  setorder(joined, station, datetime_utc)
  joined[, first_time := min(datetime_utc, na.rm = TRUE), by = station]
  joined[, hour := as.numeric(difftime(datetime_utc, first_time, units = "hours"))]
  joined[, year_op := hour / (24 * 365.25)]
}

joined[, op_bin := floor(hour / AGE_BIN_HOURS) * AGE_BIN_HOURS]
joined[, op_year_bin := op_bin / (24 * 365.25)]

# ------------------------------------------------------------
# 2) ADD KÖPPEN CLIMATE ZONE
# ------------------------------------------------------------

say("Adding single-zone Köppen climate labels")

get_single_zone_table <- function() {
  if (!file.exists(transition_csv)) stop("Missing transition file: ", transition_csv)

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

load_climate_stations <- function(infile) {
  e <- new.env(parent = emptyenv())
  load(infile, envir = e)

  if (exists("outdoor", envir = e)) {
    x <- get("outdoor", envir = e)
  } else if (exists("outdoor_final", envir = e)) {
    x <- get("outdoor_final", envir = e)
  } else {
    stop("Could not find outdoor/outdoor_final in climate_rdata")
  }

  setDT(x)
  x[, station := as.character(station)]
  unique(x[, .(station)])
}

single_zone <- get_single_zone_table()
climate_stations <- load_climate_stations(climate_rdata)

single_zone <- merge(single_zone, climate_stations, by = "station", all = FALSE)

joined <- merge(joined, single_zone, by = "station", all.x = TRUE)
joined[, has_single_zone_climate := !is.na(climate_zone)]

say("Joined PA–EPA rows: ", format(nrow(joined), big.mark = ","))
say("PA–EPA stations: ", format(uniqueN(joined$station), big.mark = ","))
say("Rows with valid single-zone climate: ", format(nrow(joined[has_single_zone_climate == TRUE]), big.mark = ","))
say("Stations with valid single-zone climate: ", format(uniqueN(joined[has_single_zone_climate == TRUE]$station), big.mark = ","))

fwrite(
  unique(joined[, .(
    station, aqs_id, state_name, county_name,
    distance_m, climate_zone, climate_desc, climate_label,
    has_single_zone_climate
  )]),
  file.path(outdir, "pa_epa_station_climate_linkage.csv")
)

# ------------------------------------------------------------
# 3) RAW ERRORS
# ------------------------------------------------------------

joined[, `:=`(
  error_raw = abs(pa_pm25_cf1 - epa_pm25),
  signed_error_raw = pa_pm25_cf1 - epa_pm25,
  nerror_raw = calc_nerr(pa_pm25_cf1, epa_pm25)
)]

# ------------------------------------------------------------
# 4) CALIBRATION MODELS
# ------------------------------------------------------------

say("Fitting calibration models")

c1 <- lm(epa_pm25 ~ pa_pm25_cf1, data = joined)

c2 <- lm(
  epa_pm25 ~ pa_pm25_cf1 + I(pa_pm25_cf1^2),
  data = joined
)

c3 <- lm(
  epa_pm25 ~ splines::ns(pa_pm25_cf1, df = 4),
  data = joined
)

c4 <- mgcv::bam(
  epa_pm25 ~ s(pa_pm25_cf1, bs = "ps"),
  data = joined,
  family = gaussian(),
  method = "fREML",
  discrete = TRUE
)

joined[, `:=`(
  pm_c1 = as.numeric(predict(c1, joined)),
  pm_c2 = as.numeric(predict(c2, joined)),
  pm_c3 = as.numeric(predict(c3, joined)),
  pm_c4 = as.numeric(predict(c4, joined))
)]

# ------------------------------------------------------------
# 5) CALIBRATED ERRORS
# ------------------------------------------------------------

joined[, `:=`(
  error_c1 = calc_err(pm_c1, epa_pm25),
  error_c2 = calc_err(pm_c2, epa_pm25),
  error_c3 = calc_err(pm_c3, epa_pm25),
  error_c4 = calc_err(pm_c4, epa_pm25),

  signed_error_c1 = pm_c1 - epa_pm25,
  signed_error_c2 = pm_c2 - epa_pm25,
  signed_error_c3 = pm_c3 - epa_pm25,
  signed_error_c4 = pm_c4 - epa_pm25,

  nerror_c1 = calc_nerr(pm_c1, epa_pm25),
  nerror_c2 = calc_nerr(pm_c2, epa_pm25),
  nerror_c3 = calc_nerr(pm_c3, epa_pm25),
  nerror_c4 = calc_nerr(pm_c4, epa_pm25)
)]

# ------------------------------------------------------------
# 6) CALIBRATION PERFORMANCE TABLE
# ------------------------------------------------------------

rmse <- function(x, y) sqrt(mean((x - y)^2, na.rm = TRUE))

perf <- rbindlist(list(
  data.table(model = "Raw",       pred = "pa_pm25_cf1"),
  data.table(model = "Linear",    pred = "pm_c1"),
  data.table(model = "Quadratic", pred = "pm_c2"),
  data.table(model = "Spline",    pred = "pm_c3"),
  data.table(model = "GAM",       pred = "pm_c4")
))[
  ,
  {
    p <- joined[[pred]]

    .(
      correlation = safe_cor(p, joined$epa_pm25),
      rmse = rmse(p, joined$epa_pm25),
      mean_abs_error = mean(abs(p - joined$epa_pm25), na.rm = TRUE),
      mean_norm_error = mean(calc_nerr(p, joined$epa_pm25), na.rm = TRUE),
      n_rows = sum(is.finite(p) & is.finite(joined$epa_pm25)),
      n_stations = uniqueN(joined$station)
    )
  },
  by = .(model, pred)
]

fwrite(perf, file.path(outdir, "calibration_performance.csv"))

make_word_table(
  perf,
  "Table S. PA–EPA calibration performance",
  file.path(tabledir, "word_table_calibration_performance.docx"),
  note = "Performance metrics compare raw and calibrated PurpleAir PM2.5 against EPA PM2.5 for sensors within 100 m of EPA monitors."
)

# ------------------------------------------------------------
# 7) SENSOR-BIN DATASET FOR CLUSTER-ROBUST TREND TESTS
# ------------------------------------------------------------

say("Creating station-by-operational-age-bin dataset for trend tests")

trend_vars <- c(
  "nerror_raw", "nerror_c1", "nerror_c2", "nerror_c3", "nerror_c4",
  "error_raw", "error_c4",
  "signed_error_raw", "signed_error_c4"
)

sensor_bin <- joined[
  is.finite(op_bin),
  c(
    .(
      n_obs = .N,
      op_years = mean(year_op, na.rm = TRUE),
      op_bin = op_bin[1],
      state_name = as.character(first_nonmissing(state_name)),
      climate_zone = as.character(first_nonmissing(climate_zone)),
      climate_desc = as.character(first_nonmissing(climate_desc)),
      climate_label = as.character(first_nonmissing(climate_label)),
      has_single_zone_climate = any(has_single_zone_climate == TRUE, na.rm = TRUE)
    ),
    lapply(.SD, mean, na.rm = TRUE)
  ),
  by = .(station, op_bin),
  .SDcols = trend_vars
]

sensor_bin <- sensor_bin[is.finite(op_years)]

# ------------------------------------------------------------
# 8) CLUSTER-ROBUST TREND TESTS
# ------------------------------------------------------------

say("Running cluster-robust trend tests on station-bin data")

empty_trend_table <- function() {
  data.table(
    group_type = character(),
    group = character(),
    model = character(),
    outcome = character(),
    n_sensor_bins = integer(),
    n_stations = integer(),
    min_year = numeric(),
    max_year = numeric(),
    beta_per_year = numeric(),
    ci_lwr_per_year = numeric(),
    ci_upr_per_year = numeric(),
    se_cluster = numeric(),
    p_value = numeric(),
    p_value_fmt = character(),
    beta_95ci_per_year = character(),
    beta_95ci_p = character()
  )
}

trend_fit_lm <- function(d, yvar, group_name, group_type, model_label) {
  if (!(yvar %in% names(d))) return(empty_trend_table())

  d <- d[
    is.finite(get(yvar)) &
      is.finite(op_years) &
      is.finite(n_obs) &
      n_obs > 0 &
      !is.na(station)
  ]

  if (nrow(d) < MIN_N_PANEL || uniqueN(d$station) < MIN_STATIONS_PANEL) {
    return(empty_trend_table())
  }

  fit <- tryCatch(
    feols(
      as.formula(paste0(yvar, " ~ op_years | station")),
      data = d,
      weights = ~ n_obs,
      vcov = ~ station,
      notes = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) return(empty_trend_table())

  ct <- tryCatch(coeftable(fit), error = function(e) NULL)
  if (is.null(ct) || !("op_years" %in% rownames(ct))) return(empty_trend_table())

  beta <- as.numeric(ct["op_years", "Estimate"])
  se <- as.numeric(ct["op_years", "Std. Error"])
  pcol <- grep("^Pr", colnames(ct), value = TRUE)[1]
  pval <- if (!is.na(pcol)) as.numeric(ct["op_years", pcol]) else NA_real_

  out <- data.table(
    group_type = group_type,
    group = as.character(group_name)[1],
    model = model_label,
    outcome = yvar,
    n_sensor_bins = nrow(d),
    n_stations = uniqueN(d$station),
    min_year = min(d$op_years, na.rm = TRUE),
    max_year = max(d$op_years, na.rm = TRUE),
    beta_per_year = beta,
    ci_lwr_per_year = beta - 1.96 * se,
    ci_upr_per_year = beta + 1.96 * se,
    se_cluster = se,
    p_value = pval
  )

  out[, p_value_fmt := fmt_p(p_value)]
  out[, beta_95ci_per_year := fmt_ci(beta_per_year, ci_lwr_per_year, ci_upr_per_year, 5)]
  out[, beta_95ci_p := paste0(beta_95ci_per_year, "; p = ", p_value_fmt)]

  out[]
}

run_trend_set <- function(d, group_name, group_type) {
  rbindlist(list(
    trend_fit_lm(d, "nerror_raw", group_name, group_type, "Raw normalized error"),
    trend_fit_lm(d, "nerror_c1",  group_name, group_type, "Linear-calibrated normalized error"),
    trend_fit_lm(d, "nerror_c2",  group_name, group_type, "Quadratic-calibrated normalized error"),
    trend_fit_lm(d, "nerror_c3",  group_name, group_type, "Spline-calibrated normalized error"),
    trend_fit_lm(d, "nerror_c4",  group_name, group_type, "GAM-calibrated normalized error"),
    trend_fit_lm(d, "error_raw",  group_name, group_type, "Raw absolute error"),
    trend_fit_lm(d, "error_c4",   group_name, group_type, "GAM-calibrated absolute error"),
    trend_fit_lm(d, "signed_error_raw", group_name, group_type, "Raw signed error"),
    trend_fit_lm(d, "signed_error_c4",  group_name, group_type, "GAM-calibrated signed error")
  ), use.names = TRUE, fill = TRUE)
}

trend_global <- run_trend_set(sensor_bin, "Global", "Global")

trend_state <- sensor_bin[
  !is.na(state_name),
  run_trend_set(.SD, state_name[1], "State"),
  by = state_name
][, state_name := NULL]

trend_climate <- sensor_bin[
  has_single_zone_climate == TRUE & !is.na(climate_label),
  run_trend_set(.SD, climate_label[1], "Köppen climate zone"),
  by = climate_label
][, climate_label := NULL]

trend_all <- rbindlist(
  list(trend_global, trend_state, trend_climate),
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  trend_all,
  file.path(outdir, "pa_epa_degradation_trend_tests_global_state_climate_cluster_robust.csv")
)

make_word_table(
  trend_all,
  "Table S. PA–EPA degradation trends by global, state, and climate groups",
  file.path(tabledir, "word_table_pa_epa_degradation_trend_tests_cluster_robust.docx"),
  note = paste0(
    "Trend estimates are beta coefficients per operational year. ",
    "Models were fit to station-by-operational-age-bin means, weighted by the number of contributing hourly observations, ",
    "with station fixed effects and sensor-clustered robust standard errors."
  )
)

cat("\nTrend models estimated:\n")
print(trend_all[, .N, by = .(group_type, model)][order(group_type, model)])

# ------------------------------------------------------------
# 9) BINNED AGE CURVES
# ------------------------------------------------------------

make_curve <- function(d, by_vars = NULL) {
  byv <- c(by_vars, "op_bin")

  d[
    is.finite(op_bin),
    .(
      n = .N,
      n_stations = uniqueN(station),
      mean_raw = mean(nerror_raw, na.rm = TRUE),
      mean_c1 = mean(nerror_c1, na.rm = TRUE),
      mean_c2 = mean(nerror_c2, na.rm = TRUE),
      mean_c3 = mean(nerror_c3, na.rm = TRUE),
      mean_c4 = mean(nerror_c4, na.rm = TRUE),
      mean_abs_raw = mean(error_raw, na.rm = TRUE),
      mean_abs_c4 = mean(error_c4, na.rm = TRUE),
      mean_signed_raw = mean(signed_error_raw, na.rm = TRUE),
      mean_signed_c4 = mean(signed_error_c4, na.rm = TRUE)
    ),
    by = byv
  ][order(op_bin)]
}

curve_global <- make_curve(joined)
curve_global[, group := "Global"]
curve_global[, group_type := "Global"]

curve_state <- make_curve(joined[!is.na(state_name)], "state_name")
setnames(curve_state, "state_name", "group")
curve_state[, group_type := "State"]

curve_climate <- make_curve(
  joined[has_single_zone_climate == TRUE],
  c("climate_zone", "climate_desc", "climate_label")
)
setnames(curve_climate, "climate_label", "group")
curve_climate[, group_type := "Köppen climate zone"]

fwrite(curve_global, file.path(outdir, "curves_global.csv"))
fwrite(curve_state, file.path(outdir, "curves_by_state.csv"))
fwrite(curve_climate, file.path(outdir, "curves_by_koppen_climate_zone.csv"))

# ------------------------------------------------------------
# Remaining plotting, high-PM interaction, age-aware GAM,
# sensor summaries, maps, and final save can stay as in your current chunk,
# now using this robust make_word_table() and trend_all table.
# ------------------------------------------------------------

# ------------------------------------------------------------
# 10) PLOTS
# ------------------------------------------------------------

curve_long <- melt(
  curve_global,
  id.vars = c("op_bin", "n", "n_stations", "group", "group_type"),
  measure.vars = c("mean_raw", "mean_c1", "mean_c2", "mean_c3", "mean_c4"),
  variable.name = "model",
  value.name = "nerror"
)

curve_long[, model := factor(
  model,
  levels = c("mean_raw", "mean_c1", "mean_c2", "mean_c3", "mean_c4"),
  labels = c("Raw", "Linear", "Quadratic", "Spline", "GAM")
)]

p_curve <- ggplot(curve_long, aes(op_bin, nerror, color = model)) +
  geom_line(linewidth = 0.8) +
  scale_x_continuous(labels = comma) +
  labs(
    title = "PA–EPA normalized error over sensor age",
    subtitle = "Raw and calibrated PurpleAir PM2.5; PA–EPA matches within 100 m",
    x = "Operational hours since first PurpleAir observation",
    y = "Normalized absolute error",
    color = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(
  file.path(figdir, "global_binned_normalized_error_raw_vs_calibrated.jpeg"),
  p_curve,
  width = 11,
  height = 7,
  dpi = 500
)

say("Fitting global GAM degradation smooths")

pd <- data.table(
  hour = seq(0, quantile(joined$hour, 0.995, na.rm = TRUE), by = 24)
)

fit_raw <- bam(
  nerror_raw ~ s(hour, bs = "ps"),
  joined,
  family = gaussian(),
  method = "fREML",
  discrete = TRUE
)

fit_cal <- bam(
  nerror_c4 ~ s(hour, bs = "ps"),
  joined,
  family = gaussian(),
  method = "fREML",
  discrete = TRUE
)

pr_raw <- predict(fit_raw, newdata = pd, se.fit = TRUE)
pr_cal <- predict(fit_cal, newdata = pd, se.fit = TRUE)

pd[, `:=`(
  raw_fit = pr_raw$fit,
  raw_low = pr_raw$fit - 1.96 * pr_raw$se.fit,
  raw_high = pr_raw$fit + 1.96 * pr_raw$se.fit,
  cal_fit = pr_cal$fit,
  cal_low = pr_cal$fit - 1.96 * pr_cal$se.fit,
  cal_high = pr_cal$fit + 1.96 * pr_cal$se.fit
)]

p_gam <- ggplot(pd, aes(hour)) +
  geom_ribbon(aes(ymin = raw_low, ymax = raw_high), fill = "grey60", alpha = 0.35) +
  geom_line(aes(y = raw_fit), color = "black", linewidth = 1) +
  geom_ribbon(aes(ymin = cal_low, ymax = cal_high), fill = "steelblue", alpha = 0.25) +
  geom_line(aes(y = cal_fit), color = "blue4", linewidth = 1) +
  scale_x_continuous(labels = comma) +
  labs(
    title = "Age-dependent PA–EPA error",
    subtitle = "Black = raw PurpleAir; blue = GAM-calibrated PurpleAir; ribbons are approximate GAM confidence intervals",
    x = "Operational hours since first PurpleAir observation",
    y = "Normalized absolute error"
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(figdir, "global_gam_age_degradation_raw_vs_calibrated.jpeg"),
  p_gam,
  width = 10,
  height = 6,
  dpi = 500
)

# ------------------------------------------------------------
# 11) STATE AND CLIMATE HETEROGENEITY PLOTS
# ------------------------------------------------------------

state_plot <- curve_state[n >= MIN_N_PANEL & n_stations >= MIN_STATIONS_PANEL]

p_state <- ggplot(state_plot, aes(op_bin, mean_c4)) +
  geom_line(linewidth = 0.65) +
  facet_wrap(~ group, scales = "free_y") +
  scale_x_continuous(labels = comma) +
  labs(
    title = "Age-dependent calibrated PA–EPA disagreement by state",
    subtitle = "GAM-calibrated PurpleAir; panels require minimum observation and station support",
    x = "Operational hours since first PurpleAir observation",
    y = "Normalized absolute error"
  ) +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(face = "bold", size = 7))

ggsave(
  file.path(figdir, "state_level_calibrated_degradation_curves.jpeg"),
  p_state,
  width = 15,
  height = 10,
  dpi = 500,
  limitsize = FALSE
)

climate_plot <- curve_climate[n >= MIN_N_PANEL & n_stations >= MIN_STATIONS_PANEL]

climate_order <- climate_plot[
  ,
  .(total_n = sum(n, na.rm = TRUE), total_stations = max(n_stations, na.rm = TRUE)),
  by = group
][order(-total_n)]$group

climate_plot[, group := factor(group, levels = climate_order)]

p_climate <- ggplot(climate_plot, aes(op_bin, mean_c4)) +
  geom_line(linewidth = 0.7, color = "blue4") +
  facet_wrap(~ group, scales = "free_y", ncol = 4) +
  scale_x_continuous(labels = comma) +
  labs(
    title = "Age-dependent calibrated PA–EPA disagreement by Köppen climate zone",
    subtitle = "Single-zone PurpleAir stations only",
    x = "Operational hours since first PurpleAir observation",
    y = "Normalized absolute error"
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.text = element_text(face = "bold", size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(figdir, "koppen_climate_calibrated_degradation_curves.jpeg"),
  p_climate,
  width = 14,
  height = max(8, 2.6 * ceiling(length(climate_order) / 4)),
  dpi = 500,
  limitsize = FALSE
)

climate_long <- melt(
  climate_plot,
  id.vars = c("climate_zone", "climate_desc", "group", "group_type", "op_bin", "n", "n_stations"),
  measure.vars = c("mean_raw", "mean_c4"),
  variable.name = "model",
  value.name = "nerror"
)

climate_long[, model := factor(
  model,
  levels = c("mean_raw", "mean_c4"),
  labels = c("Raw", "GAM calibrated")
)]

p_climate_raw_cal <- ggplot(climate_long, aes(op_bin, nerror, color = model)) +
  geom_line(linewidth = 0.65) +
  facet_wrap(~ group, scales = "free_y", ncol = 4) +
  scale_x_continuous(labels = comma) +
  labs(
    title = "Raw vs calibrated PA–EPA disagreement by Köppen climate zone",
    subtitle = "Single-zone stations only",
    x = "Operational hours since first PurpleAir observation",
    y = "Normalized absolute error",
    color = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 7),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(figdir, "koppen_climate_raw_vs_calibrated_degradation_curves.jpeg"),
  p_climate_raw_cal,
  width = 14,
  height = max(8, 2.6 * ceiling(length(climate_order) / 4)),
  dpi = 500,
  limitsize = FALSE
)

# ------------------------------------------------------------
# 12) HIGH-PM EXPOSURE AND AGE INTERACTION
# ------------------------------------------------------------

say("Evaluating cumulative high-PM exposure")

setorder(joined, station, datetime_utc)

joined[, `:=`(
  pm50 = as.integer(pa_pm25_cf1 > 50),
  pm100 = as.integer(pa_pm25_cf1 > 100)
)]

joined[, `:=`(
  cumpm50 = cumsum(pm50),
  cumpm100 = cumsum(pm100)
), by = station]

fit_pm50 <- feols(
  nerror_c4 ~ year_op * cumpm50 | station,
  data = joined[is.finite(nerror_c4)],
  vcov = ~ station,
  notes = FALSE
)

fit_pm100 <- feols(
  nerror_c4 ~ year_op * cumpm100 | station,
  data = joined[is.finite(nerror_c4)],
  vcov = ~ station,
  notes = FALSE
)

capture.output(summary(fit_pm50), file = file.path(outdir, "feols_nerror_c4_yearop_x_cumpm50_cluster_robust.txt"))
capture.output(summary(fit_pm100), file = file.path(outdir, "feols_nerror_c4_yearop_x_cumpm100_cluster_robust.txt"))

if (nrow(joined[has_single_zone_climate == TRUE]) > 0) {
  fit_climate_interaction <- feols(
    nerror_c4 ~ year_op * climate_label | station,
    data = joined[has_single_zone_climate == TRUE & is.finite(nerror_c4)],
    vcov = ~ station,
    notes = FALSE
  )

  capture.output(
    summary(fit_climate_interaction),
    file = file.path(outdir, "feols_nerror_c4_yearop_x_koppen_climate_cluster_robust.txt")
  )
}

# ------------------------------------------------------------
# 13) AGE-AWARE PA–EPA CALIBRATION
# ------------------------------------------------------------

say("Running age-aware PA–EPA calibration/degradation models")

mod_dt <- joined[
  is.finite(pa_pm25_cf1) &
    is.finite(epa_pm25) &
    is.finite(hour) &
    !is.na(station),
  .(
    station,
    datetime_utc,
    hour,
    year_op,
    op_bin,
    op_year_bin,
    pa_pm25_cf1,
    epa_pm25,
    state_name,
    region,
    climate_zone,
    climate_desc,
    climate_label,
    has_single_zone_climate,
    lat = pa_lat,
    lon = pa_lon,
    aqs_id,
    distance_m,
    pm_c4,
    nerror_c4,
    error_c4
  )
]

MAX_MODEL_ROWS <- 2500000L

if (nrow(mod_dt) > MAX_MODEL_ROWS) {
  say("Sampling rows for age-aware GAMs: ", format(MAX_MODEL_ROWS, big.mark = ","))
  set.seed(123)
  mod_fit <- mod_dt[sample.int(.N, MAX_MODEL_ROWS)]
} else {
  mod_fit <- copy(mod_dt)
}

mod_fit[, hour_k := hour / 1000]
mod_dt[, hour_k := hour / 1000]

say("Fitting additive age-aware calibration: EPA ~ s(PA) + s(age)")

gam_add_age <- mgcv::bam(
  epa_pm25 ~ s(pa_pm25_cf1, bs = "ps") + s(hour_k, bs = "ps"),
  data = mod_fit,
  family = gaussian(),
  method = "fREML",
  discrete = TRUE
)

say("Fitting tensor age-aware calibration: EPA ~ te(PA, age)")

gam_te_age <- mgcv::bam(
  epa_pm25 ~ te(pa_pm25_cf1, hour_k, bs = c("ps", "ps")),
  data = mod_fit,
  family = gaussian(),
  method = "fREML",
  discrete = TRUE
)

capture.output(summary(gam_add_age), file = file.path(outdir, "gam_additive_age_calibration_summary.txt"))
capture.output(summary(gam_te_age), file = file.path(outdir, "gam_tensor_age_calibration_summary.txt"))

say("Predicting age-aware calibrated PM2.5")

mod_dt[, `:=`(
  pm_age_add = as.numeric(predict(gam_add_age, newdata = mod_dt)),
  pm_age_te  = as.numeric(predict(gam_te_age,  newdata = mod_dt))
)]

mod_dt[, `:=`(
  error_age_add = abs(pm_age_add - epa_pm25),
  error_age_te  = abs(pm_age_te  - epa_pm25),
  nerror_age_add = calc_nerr(pm_age_add, epa_pm25),
  nerror_age_te  = calc_nerr(pm_age_te,  epa_pm25),
  signed_error_age_add = pm_age_add - epa_pm25,
  signed_error_age_te  = pm_age_te  - epa_pm25
)]

perf_age <- rbindlist(list(
  data.table(model = "Raw", pred = "pa_pm25_cf1"),
  data.table(model = "Original GAM: EPA ~ s(PA)", pred = "pm_c4"),
  data.table(model = "Age additive GAM: EPA ~ s(PA)+s(age)", pred = "pm_age_add"),
  data.table(model = "Age tensor GAM: EPA ~ te(PA,age)", pred = "pm_age_te")
), fill = TRUE)[
  ,
  {
    p <- mod_dt[[pred]]
    .(
      correlation = safe_cor(p, mod_dt$epa_pm25),
      rmse = sqrt(mean((p - mod_dt$epa_pm25)^2, na.rm = TRUE)),
      mean_abs_error = mean(abs(p - mod_dt$epa_pm25), na.rm = TRUE),
      mean_norm_error = mean(calc_nerr(p, mod_dt$epa_pm25), na.rm = TRUE),
      n_rows = sum(is.finite(p) & is.finite(mod_dt$epa_pm25)),
      n_stations = uniqueN(mod_dt$station)
    )
  },
  by = .(model, pred)
]

fwrite(perf_age, file.path(outdir, "calibration_performance_age_aware_models.csv"))

make_word_table(
  perf_age,
  "Table S. Age-aware PA–EPA calibration performance",
  file.path(tabledir, "word_table_age_aware_calibration_performance.docx"),
  note = "Performance metrics compare raw, original GAM-calibrated, and age-aware calibrated PurpleAir PM2.5 against EPA PM2.5."
)

# ------------------------------------------------------------
# 14) CLUSTER-ROBUST TRENDS FOR AGE-AWARE RESIDUALS
# ------------------------------------------------------------

say("Testing age trends in age-aware calibration residuals")

age_trend_vars <- c(
  "nerror_age_add", "nerror_age_te",
  "error_age_add", "error_age_te",
  "signed_error_age_add", "signed_error_age_te"
)

age_sensor_bin <- mod_dt[
  is.finite(op_bin),
  c(
    .(
      n_obs = .N,
      op_years = mean(year_op, na.rm = TRUE),
      op_bin = op_bin[1]
    ),
    lapply(.SD, mean, na.rm = TRUE)
  ),
  by = .(station, op_bin),
  .SDcols = age_trend_vars
]

trend_age_models <- function(d, yvar, label) {
  trend_fit_lm(d, yvar, "Global", "Global", label)
}

trend_age <- rbindlist(list(
  trend_age_models(age_sensor_bin, "nerror_age_add", "Age additive GAM normalized error"),
  trend_age_models(age_sensor_bin, "nerror_age_te",  "Age tensor GAM normalized error"),
  trend_age_models(age_sensor_bin, "error_age_add",  "Age additive GAM absolute error"),
  trend_age_models(age_sensor_bin, "error_age_te",   "Age tensor GAM absolute error"),
  trend_age_models(age_sensor_bin, "signed_error_age_add", "Age additive GAM signed error"),
  trend_age_models(age_sensor_bin, "signed_error_age_te",  "Age tensor GAM signed error")
), fill = TRUE)

fwrite(
  trend_age,
  file.path(outdir, "age_aware_calibration_degradation_trend_tests_cluster_robust.csv")
)

make_word_table(
  trend_age,
  "Table S. Age-aware calibration degradation trends",
  file.path(tabledir, "word_table_age_aware_degradation_tests_cluster_robust.docx"),
  note = paste0(
    "Trend estimates are beta coefficients per operational year. ",
    "Models were fit to station-by-operational-age-bin means, weighted by the number of contributing hourly observations, ",
    "with station fixed effects and sensor-clustered robust standard errors."
  )
)

# ------------------------------------------------------------
# REPLACE SECTION 15 WITH THIS ROBUST VERSION
# ------------------------------------------------------------

say("Creating EPA-colocated sensor summary table")

first_nonmissing <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0L) return(NA)
  x[1]
}

sensor_summary_epa <- mod_dt[
  ,
  .(
    n_hours = as.integer(.N),
    first_obs = min(datetime_utc, na.rm = TRUE),
    last_obs = max(datetime_utc, na.rm = TRUE),
    duration_days = as.numeric(
      difftime(
        max(datetime_utc, na.rm = TRUE),
        min(datetime_utc, na.rm = TRUE),
        units = "days"
      )
    ),
    max_operational_hour = as.numeric(max(hour, na.rm = TRUE)),
    max_operational_year = as.numeric(max(hour, na.rm = TRUE) / (24 * 365.25)),
    mean_pa = as.numeric(mean(pa_pm25_cf1, na.rm = TRUE)),
    mean_epa = as.numeric(mean(epa_pm25, na.rm = TRUE)),
    mean_raw_abs_error = as.numeric(mean(abs(pa_pm25_cf1 - epa_pm25), na.rm = TRUE)),
    mean_raw_norm_error = as.numeric(mean(calc_nerr(pa_pm25_cf1, epa_pm25), na.rm = TRUE)),
    mean_age_add_norm_error = as.numeric(mean(nerror_age_add, na.rm = TRUE)),
    mean_age_te_norm_error = as.numeric(mean(nerror_age_te, na.rm = TRUE)),
    median_distance_m = as.numeric(median(distance_m, na.rm = TRUE)),
    lat = as.numeric(median(lat, na.rm = TRUE)),
    lon = as.numeric(median(lon, na.rm = TRUE)),
    state_name = as.character(first_nonmissing(state_name)),
    region = as.character(first_nonmissing(region)),
    climate_zone = as.character(first_nonmissing(climate_zone)),
    climate_desc = as.character(first_nonmissing(climate_desc)),
    climate_label = as.character(first_nonmissing(climate_label))
  ),
  by = .(station, aqs_id)
]

sensor_summary_epa[, hours_per_day_observed := as.numeric(n_hours) / pmax(duration_days, 1)]

fwrite(
  sensor_summary_epa,
  file.path(outdir, "epa_colocated_sensor_summary_stats.csv")
)

epa_summary_state <- sensor_summary_epa[
  !is.na(state_name),
  .(
    n_sensors = as.integer(uniqueN(station)),
    n_aqs_monitors = as.integer(uniqueN(aqs_id)),
    median_hours = as.numeric(median(as.numeric(n_hours), na.rm = TRUE)),
    median_duration_days = as.numeric(median(as.numeric(duration_days), na.rm = TRUE)),
    median_operational_year = as.numeric(median(as.numeric(max_operational_year), na.rm = TRUE)),
    median_age_te_norm_error = as.numeric(median(as.numeric(mean_age_te_norm_error), na.rm = TRUE))
  ),
  by = state_name
][order(-n_sensors)]

epa_summary_climate <- sensor_summary_epa[
  !is.na(climate_label),
  .(
    n_sensors = as.integer(uniqueN(station)),
    n_aqs_monitors = as.integer(uniqueN(aqs_id)),
    median_hours = as.numeric(median(as.numeric(n_hours), na.rm = TRUE)),
    median_duration_days = as.numeric(median(as.numeric(duration_days), na.rm = TRUE)),
    median_operational_year = as.numeric(median(as.numeric(max_operational_year), na.rm = TRUE)),
    median_age_te_norm_error = as.numeric(median(as.numeric(mean_age_te_norm_error), na.rm = TRUE))
  ),
  by = .(climate_zone, climate_desc, climate_label)
][order(-n_sensors)]

fwrite(
  epa_summary_state,
  file.path(outdir, "epa_colocated_summary_by_state.csv")
)

fwrite(
  epa_summary_climate,
  file.path(outdir, "epa_colocated_summary_by_climate.csv")
)

make_word_table(
  epa_summary_state,
  "Table S. EPA-colocated PurpleAir sensor summary by state",
  file.path(tabledir, "word_table_epa_colocated_summary_by_state.docx")
)

make_word_table(
  epa_summary_climate,
  "Table S. EPA-colocated PurpleAir sensor summary by Köppen climate zone",
  file.path(tabledir, "word_table_epa_colocated_summary_by_climate.docx")
)

# ------------------------------------------------------------
# 16) MAP AND COVERAGE FIGURES
# ------------------------------------------------------------

say("Drawing EPA-colocated sensor map")

map_dt <- sensor_summary_epa[is.finite(lat) & is.finite(lon)]

p_map <- ggplot(map_dt, aes(lon, lat)) +
  borders("state", colour = "grey75", fill = NA, linewidth = 0.25) +
  geom_point(aes(size = n_hours, color = max_operational_year), alpha = 0.75) +
  scale_size_continuous(labels = comma, range = c(1.5, 7), name = "Matched hours") +
  scale_color_viridis_c(option = "plasma", name = "Operational years") +
  coord_quickmap(xlim = c(-126, -66), ylim = c(24, 50)) +
  labs(
    title = "PurpleAir sensors within 100 m of EPA PM2.5 monitors",
    subtitle = "Point size = matched PA–EPA hourly observations; color = observed operational lifetime",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

ggsave(
  file.path(figdir, "map_epa_colocated_purpleair_sensors.jpeg"),
  p_map,
  width = 11,
  height = 7,
  dpi = 500
)

say("Drawing EPA-colocated sensor Gantt chart")

gantt_dt <- copy(sensor_summary_epa)
gantt_dt <- gantt_dt[order(first_obs)]
gantt_dt[, sensor_order := frank(first_obs, ties.method = "dense")]

p_gantt <- ggplot(
  gantt_dt,
  aes(
    x = first_obs,
    xend = last_obs,
    y = reorder(station, first_obs),
    yend = reorder(station, first_obs),
    color = n_hours
  )
) +
  geom_segment(linewidth = 0.45, alpha = 0.85) +
  scale_color_viridis_c(option = "plasma", trans = "log10", labels = comma, name = "Matched hours\n(log scale)") +
  labs(
    title = "Observation windows for EPA-colocated PurpleAir sensors",
    subtitle = "Sensors within 100 m of EPA PM2.5 monitors",
    x = "Calendar time",
    y = "PurpleAir sensor"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  file.path(figdir, "gantt_epa_colocated_sensor_coverage.jpeg"),
  p_gantt,
  width = 12,
  height = 9,
  dpi = 500
)

# ------------------------------------------------------------
# 17) MAIN AGE-AWARE DEGRADATION FIGURE
# ------------------------------------------------------------

say("Creating main PA–EPA degradation figure")

curve_age <- mod_dt[
  is.finite(op_bin),
  .(
    n = .N,
    n_stations = uniqueN(station),
    raw = mean(calc_nerr(pa_pm25_cf1, epa_pm25), na.rm = TRUE),
    original_gam = mean(nerror_c4, na.rm = TRUE),
    age_additive_gam = mean(nerror_age_add, na.rm = TRUE),
    age_tensor_gam = mean(nerror_age_te, na.rm = TRUE),
    raw_abs = mean(abs(pa_pm25_cf1 - epa_pm25), na.rm = TRUE),
    age_tensor_abs = mean(error_age_te, na.rm = TRUE)
  ),
  by = op_bin
][
  n >= 1000 & n_stations >= 20
][order(op_bin)]

curve_age_long <- melt(
  curve_age,
  id.vars = c("op_bin", "n", "n_stations"),
  measure.vars = c("raw", "original_gam", "age_additive_gam", "age_tensor_gam"),
  variable.name = "model",
  value.name = "normalized_error"
)

curve_age_long[, model := factor(
  model,
  levels = c("raw", "original_gam", "age_additive_gam", "age_tensor_gam"),
  labels = c(
    "Raw PurpleAir",
    "Original GAM: EPA ~ s(PA)",
    "Age additive GAM: EPA ~ s(PA)+s(age)",
    "Age tensor GAM: EPA ~ te(PA, age)"
  )
)]

p_main <- ggplot(curve_age_long, aes(op_bin, normalized_error, color = model)) +
  geom_line(linewidth = 0.9) +
  scale_x_continuous(labels = comma) +
  labs(
    title = "PA–EPA disagreement over PurpleAir sensor age",
    subtitle = "EPA-colocated sensors within 100 m; operational age is measured from first PurpleAir observation",
    x = "Operational hours since first PurpleAir observation",
    y = "Mean normalized absolute PA–EPA error",
    color = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(
  file.path(figdir, "main_pa_epa_degradation_age_aware_calibration.jpeg"),
  p_main,
  width = 12,
  height = 7,
  dpi = 500
)

# ------------------------------------------------------------
# 18) AGE EFFECT SURFACE FROM TENSOR MODEL
# ------------------------------------------------------------

say("Creating tensor GAM prediction surface")

pa_seq <- quantile(mod_dt$pa_pm25_cf1, probs = seq(0.02, 0.98, length.out = 80), na.rm = TRUE)
hour_seq <- seq(
  quantile(mod_dt$hour, 0.01, na.rm = TRUE),
  quantile(mod_dt$hour, 0.99, na.rm = TRUE),
  length.out = 80
)

surface_dt <- CJ(
  pa_pm25_cf1 = as.numeric(pa_seq),
  hour = as.numeric(hour_seq)
)

surface_dt[, hour_k := hour / 1000]
surface_dt[, pred_epa := as.numeric(predict(gam_te_age, newdata = surface_dt))]
surface_dt[, bias := pa_pm25_cf1 - pred_epa]

p_surface <- ggplot(surface_dt, aes(hour, pa_pm25_cf1, fill = bias)) +
  geom_tile() +
  scale_x_continuous(labels = comma) +
  scale_fill_gradient2(
    low = "blue3",
    mid = "white",
    high = "red3",
    midpoint = 0,
    name = "PA − predicted EPA"
  ) +
  labs(
    title = "Age-aware calibration surface",
    subtitle = "Tensor GAM: EPA PM2.5 ~ te(PurpleAir PM2.5, sensor age)",
    x = "Operational hours since first PurpleAir observation",
    y = "PurpleAir PM2.5"
  ) +
  theme_bw(base_size = 12)

ggsave(
  file.path(figdir, "tensor_gam_age_pa_bias_surface.jpeg"),
  p_surface,
  width = 10,
  height = 7,
  dpi = 500
)

# ------------------------------------------------------------
# 19) FINAL SAVE
# ------------------------------------------------------------

save(
  joined,
  perf,
  trend_all,
  curve_global,
  curve_state,
  curve_climate,
  pd,
  c1, c2, c3, c4,
  fit_raw,
  fit_cal,
  gam_add_age,
  gam_te_age,
  perf_age,
  trend_age,
  sensor_summary_epa,
  epa_summary_state,
  epa_summary_climate,
  curve_age,
  file = file.path(outdir, "pa_epa_calibration_degradation_with_koppen_cluster_robust.RData")
)

cat("\nDone. Outputs saved to:\n")
cat(outdir, "\n")
cat("\nWord-ready tables saved to:\n")
cat(tabledir, "\n")
