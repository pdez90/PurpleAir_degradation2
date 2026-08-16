# ==========================================================================
# 06_clean_outdoor_final.R
#
# QC and cleaning -> the analytic dataset outdoor_final.RData.
# CHECKPOINT: 452,143,148 rows / 36,328 sensors.
# (analysis_addons/00_build_analysis_dataset.R rebuilds and verifies these counts.)
#
# Extracted VERBATIM from purpleair_degradation2.Rmd (chunk #20 lines 4925-5319).
# Code is unmodified; see REPRODUCING.md for paths and run order.
# ==========================================================================

# ---- source chunk #20 (Rmd lines 4925-5319) ------------------------

# ============================================================
# QC + SUMMARY STATS + PM EXTREMES + CORRELATIONS
# ============================================================

library(data.table)
outdoor<-readRDS("/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_all_joined.rds")
# ---- 0) Read + data.table ----
outdoor <- readRDS("/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_all_joined.rds")
setDT(outdoor)

# ---- 1) Column names ----
a_col   <- "pm25_cf_1_a_hourly(ug/m3)_pm25_cf_1_a"
b_col   <- "pm25_cf_1_b_hourly(ug/m3)_pm25_cf_1_b"
tcol    <- "temperature_hourly(C)_temperature"
rhcol   <- "humidity_hourly(%)_humidity"
era5_tcol  <- "temperature_2m"
era5_tdcol <- "dewpoint_temperature_2m"
geoscf_col <- "geoscf_pm25"

# ---- 2) Detect coordinate columns safely ----
pa_lon <- intersect(c("LONGITUDE(deg)", "longitude"), names(outdoor))
pa_lat <- intersect(c("LATITUDE(deg)",  "latitude"),  names(outdoor))

gee_lon <- intersect(c("gee_longitude"), names(outdoor))
gee_lat <- intersect(c("gee_latitude"),  names(outdoor))

# ---- 3) Build KEEP list (robust) ----
keep <- unique(intersect(
  c(
    "datetime_utc","station","region",
    pa_lon, pa_lat,
    gee_lon, gee_lat,
    "era_pid","pm_pid","gee_pid",
    a_col,b_col,tcol,rhcol,
    era5_tcol,era5_tdcol,geoscf_col
  ),
  names(outdoor)
))

outdoor <- outdoor[, ..keep]
#498135221

gc()
# ---- 2) Require BOTH PM cf1 A & B ----
outdoor <- outdoor[!is.na(get(a_col)) & !is.na(get(b_col))]
#493787385
gc()

# create hourly timestamp
outdoor[, ts_hr := lubridate::floor_date(as.POSIXct(datetime_utc, tz = "UTC"), "hour")]

# ---- 3) Fast hourly timestamp ----
outdoor[, ts_hr :=
  as.POSIXct((as.numeric(as.POSIXct(datetime_utc, tz="UTC")) %/% 3600L) * 3600L,
             origin="1970-01-01", tz="UTC")]
gc()

# ============================================================
# A) SENSOR COUNTS + MISSING HOURS (PM-filtered)
# ============================================================

sensor_stats <- outdoor[
  , .(
    n_obs_hours = uniqueN(ts_hr),
    first_ts = min(ts_hr),
    last_ts  = max(ts_hr)
  ),
  by = .(region, station)
]

sensor_stats[, span_hours :=
  (as.numeric(last_ts) %/% 3600L) -
  (as.numeric(first_ts) %/% 3600L) + 1L
]
sensor_stats[, missing_hours := pmax(0L, span_hours - n_obs_hours)]

sensor_counts_overall <- sensor_stats[, .(n_sensors = uniqueN(station))]
sensor_counts_region  <- sensor_stats[, .(n_sensors = uniqueN(station)), by=region]

overall_missing <- sensor_stats[
  , .(
    total_obs_hours = sum(n_obs_hours),
    total_possible_hours = sum(span_hours),
    total_missing_hours = sum(missing_hours),
    pct_missing_hours = 100 * sum(missing_hours) / sum(span_hours)
  )
]

region_missing <- sensor_stats[
  , .(
    total_obs_hours = sum(n_obs_hours),
    total_possible_hours = sum(span_hours),
    total_missing_hours = sum(missing_hours),
    pct_missing_hours = 100 * sum(missing_hours) / sum(span_hours)
  ),
  by = region
]

# ============================================================
# B) REMOVE NA RH / T
# ============================================================

missing_rh_t_overall <- outdoor[
  , .(
    n_total = .N,
    pct_miss_T  = 100 * sum(is.na(get(tcol))) / .N,
    pct_miss_RH = 100 * sum(is.na(get(rhcol))) / .N
  )
]

missing_rh_t_region <- outdoor[
  , .(
    n_total = .N,
    pct_miss_T   = 100 * mean(is.na(get(tcol))),
    pct_miss_RH  = 100 * mean(is.na(get(rhcol))),
    pct_miss_T_RH = 100 * mean(is.na(get(tcol)) | is.na(get(rhcol)))
  ),
  by = region
]

gc()
outdoor <- outdoor[!is.na(get(tcol)) & !is.na(get(rhcol))]
gc()

print(nrow(outdoor))
print(table(outdoor$region))
# ============================================================
# C) REMOVE EXTREME T / RH
# ============================================================

gc()
outdoor[, bad_T  := !is.na(get(tcol))  & (get(tcol) <= -50 | get(tcol) >= 100)]
outdoor[, bad_RH := !is.na(get(rhcol)) & (get(rhcol) > 99)]
outdoor[, bad_env := bad_T | bad_RH]
gc()

extremes_env_overall <- outdoor[
  , .(
    n_total = .N,
    n_bad   = sum(bad_env, na.rm = TRUE),
    pct_bad = 100 * mean(bad_env, na.rm = TRUE)
  )
]

extremes_env_region <- outdoor[
  , .(
    n_total = .N,
    n_bad   = sum(bad_env, na.rm = TRUE),
    pct_bad = 100 * mean(bad_env, na.rm = TRUE)
  ),
  by = region
]

gc()
outdoor <- outdoor[bad_env == FALSE][, c("bad_T","bad_RH","bad_env") := NULL]
gc()

print(nrow(outdoor))
#464848642
print(table(outdoor$region))
#      Africa        CONUS        India SouthAmerica 
#     3255433    452053816      2865850      6673543 
# ============================================================
# C2) REMOVE EXTREME PM (cf1 A OR B > 1500)
# ============================================================
gc()
print(nrow(outdoor))
print(table(outdoor$region))
#464848642
outdoor[, bad_pm := get(a_col) > 1500 | get(b_col) > 1500]
print(nrow(outdoor))
print(table(outdoor$region))
#464848642
#      Africa        CONUS        India SouthAmerica 
#     3255433    452053816      2865850      6673543 
gc()

pm_extreme_overall <- outdoor[
  , .(
    n_before = .N,
    n_removed = sum(bad_pm),
    pct_removed = 100 * sum(bad_pm) / .N
  )
]

pm_extreme_region <- outdoor[
  , .(
    n_before = .N,
    n_removed = sum(bad_pm),
    pct_removed = 100 * sum(bad_pm) / .N
  ),
  by = region
]

gc()
outdoor_final <- outdoor[bad_pm == FALSE][, bad_pm := NULL]
gc()
print(nrow(outdoor_final))
print(table(outdoor_final$region))
#452143148
#      Africa        CONUS        India SouthAmerica 
#     3134219    439768685      2798734      6441510 

# Remaining measurements + sensors
remaining_meas_overall <- data.table(n_remaining = nrow(outdoor_final))
remaining_meas_region  <- outdoor_final[, .(n_remaining = .N), by = region]

remaining_sensors_overall <- outdoor_final[, .(n_sensors = uniqueN(station))]
remaining_sensors_region  <- outdoor_final[, .(n_sensors = uniqueN(station)), by = region]

# ============================================================
# D) CLEANED PM + DERIVE ERA5 RH + CORRELATIONS
#   - pm25_cf1_avg = mean(cf1_a, cf1_b)
#   - Convert ERA5 T2m & Td2m to °C if needed
#   - Derive ERA5 RH (%) from T and Td
#   - Correlate pm25_cf1_avg with: GEOS-CF PM2.5, ERA5 T (°C), ERA5 RH (%)
# ============================================================

# 1) Cleaned PM
outdoor_final[, pm25_cf1_avg := (get(a_col) + get(b_col)) / 2]

# 2) Convert ERA5 units (K -> C) using scalar logic (fast + correct)
have_t  <- era5_tcol  %in% names(outdoor_final)
have_td <- era5_tdcol %in% names(outdoor_final)

if (have_t) {
  med_t <- median(outdoor_final[[era5_tcol]], na.rm = TRUE)
  if (is.finite(med_t) && med_t > 150) {
    outdoor_final[, era5_t2m_c := get(era5_tcol) - 273.15]
  } else {
    outdoor_final[, era5_t2m_c := get(era5_tcol)]
  }
}

if (have_td) {
  med_td <- median(outdoor_final[[era5_tdcol]], na.rm = TRUE)
  if (is.finite(med_td) && med_td > 150) {
    outdoor_final[, era5_td2m_c := get(era5_tdcol) - 273.15]
  } else {
    outdoor_final[, era5_td2m_c := get(era5_tdcol)]
  }
}

# 3) Convert dewpoint to RH (%) using Magnus formula
# RH = 100 * exp( (a*Td)/(b+Td) - (a*T)/(b+T) )
# Valid for typical near-surface conditions; uses °C.
if (all(c("era5_t2m_c", "era5_td2m_c") %in% names(outdoor_final))) {
  a <- 17.625
  b <- 243.04
  outdoor_final[, era5_rh_pct := {
    T  <- era5_t2m_c
    Td <- era5_td2m_c
    rh <- 100 * exp((a * Td) / (b + Td) - (a * T) / (b + T))
    pmin(100, pmax(0, rh))  # clamp to [0, 100]
  }]
}

save(outdoor_final, file="/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_final.RData")
#452143148
#      Africa        CONUS        India SouthAmerica 
#     3134219    439768685      2798734      6441510 
# ============================================================
# OUTPUT OBJECTS
# ============================================================

sensor_counts_overall
sensor_counts_region
overall_missing
region_missing

missing_rh_t_overall
missing_rh_t_region

extremes_env_overall
extremes_env_region

pm_extreme_overall
pm_extreme_region

remaining_meas_overall
remaining_meas_region
remaining_sensors_overall
remaining_sensors_region


# Cleaned dataset: outdoor_final
load("/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_final.RData")
outdoor<-outdoor_final
rm(outdoor_final)
gc()
setDT(outdoor)
cat("Rows loaded:", nrow(outdoor), "\n")
#452143148

#Correlations
library(data.table)

setDT(outdoor)

# --- Define columns you want to correlate ---
pa_pm   <- "pm25_cf1_avg"
pa_t    <- "temperature_hourly(C)_temperature"
pa_rh   <- "humidity_hourly(%)_humidity"

era_t   <- "era5_t2m_c"
era_td  <- "era5_td2m_c"
era_rh  <- "era5_rh_pct"

geos_pm <- "geoscf_pm25"

vars <- c(pa_pm, pa_t, pa_rh, era_t, era_td, era_rh, geos_pm)
vars <- vars[vars %in% names(outdoor)]

# Keep only what we need for correlation to reduce memory pressure
outdoor_cor <- outdoor[, c("region", vars), with = FALSE]
rm(outdoor); gc()

# --- Helper: correlation from sufficient statistics (scalar) ---
cor_from_sums <- function(n, sx, sy, sxx, syy, sxy) {
  num <- n * sxy - sx * sy
  den <- sqrt((n * sxx - sx * sx) * (n * syy - sy * sy))
  ifelse(n >= 3 & is.finite(den) & den > 0, num / den, NA_real_)
}

# --- Build list of variable pairs to evaluate ---
# (edit this list if you want fewer pairs)
pairs <- CJ(var1 = vars, var2 = vars)[var1 < var2]  # all unique pairs
# If you only want cross-source pairs, replace above with a custom list.

# ============================================================
# 1) OVERALL correlations (sums; no big ok vectors)
# ============================================================

cor_overall <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  v1 <- pairs$var1[i]
  v2 <- pairs$var2[i]

  # Aggregate sums over finite pairs only
  tmp <- outdoor_cor[
    is.finite(get(v1)) & is.finite(get(v2)),
    .(
      n   = .N,
      sx  = sum(get(v1)),
      sy  = sum(get(v2)),
      sxx = sum(get(v1)^2),
      syy = sum(get(v2)^2),
      sxy = sum(get(v1) * get(v2))
    )
  ]

  data.table(
    var1 = v1,
    var2 = v2,
    n_pair = tmp$n,
    cor = cor_from_sums(tmp$n, tmp$sx, tmp$sy, tmp$sxx, tmp$syy, tmp$sxy)
  )
}), use.names = TRUE)

# ============================================================
# 2) BY-REGION correlations (same approach; grouped)
# ============================================================

cor_region <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  v1 <- pairs$var1[i]
  v2 <- pairs$var2[i]

  tmp <- outdoor_cor[
    is.finite(get(v1)) & is.finite(get(v2)),
    .(
      n   = .N,
      sx  = sum(get(v1)),
      sy  = sum(get(v2)),
      sxx = sum(get(v1)^2),
      syy = sum(get(v2)^2),
      sxy = sum(get(v1) * get(v2))
    ),
    by = region
  ]

  tmp[, `:=`(
    var1 = v1,
    var2 = v2,
    cor  = cor_from_sums(n, sx, sy, sxx, syy, sxy)
  )]

  tmp[, .(region, var1, var2, n_pair = n, cor)]
}), use.names = TRUE, fill = TRUE)

setorder(cor_region, var1, var2, region)

# --- Results ---
cor_overall
cor_region
