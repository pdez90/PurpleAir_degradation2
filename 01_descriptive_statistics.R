# ==========================================================================
# 01_descriptive_statistics.R
#
# Descriptive statistics for the analytic sample: sensor movement/relocation and density
# (Figure S1); Tables S1-S2 and the manuscript sample numbers; full-data summaries and
# PurpleAir-vs-reference correlations; sensor distributions, Metric 3 thresholds and
# deployment timelines (Figures S4-S25).
#
# Extracted VERBATIM from purpleair_degradation2.Rmd (chunks #21 lines 5322-6334, #25 lines 7315-7623, #26 lines 7626-8539, #27 lines 8542-9219, #28 lines 9222-9987).
# Code is unmodified; see REPRODUCING.md for paths and run order.
# ==========================================================================

# ---- source chunk #21 (Rmd lines 5322-6334) ------------------------

# ============================================================
# EVALUATE + MAP + ANIMATE PURPLEAIR SENSOR MOVEMENT / DENSITY
# + REPEAT CORE DESCRIPTIVE DEPLOYMENT ANALYSES
# + COMBINED DEPLOYMENT FIGURE (PANEL A + B)
# + ONE SUMMARY TABLE BY REGION
#
# What this does:
# 1) loads outdoor_final
# 2) chooses time-varying coordinates (PurpleAir first, then GEE fallback)
# 3) builds monthly station positions
# 4) evaluates whether a station ever moved at all
# 5) distinguishes:
#      - moved
#      - did not move
#      - only 1 observed month
# 6) creates:
#      - monthly animation
#      - global sensor count-density map (bin2d)
#      - regional sensor count-density maps
#      - map of stations colored by months active
#      - map of stations colored by movement status
#      - histogram of months active, global + by region
#      - histogram of approximate movement distance
#      - combined deployment figure with panels A and B
# 7) writes summary CSVs
# 8) repeats core descriptive deployment analyses:
#      - deployment tenure from first to last observation
#      - longevity groups: <1, 1–2, 2–3, >3 years
#      - global + regional summaries
#      - region-by-unit counts if a unit column exists
#      - exemplar-unit PM summaries by longevity
#      - optional urban/rural summaries if urban_rural exists
# 9) creates one region summary table with:
#      - number of sensors
#      - lifetime summary stats
#      - movement summary
# ============================================================

load("/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_final.RData")
outdoor <- outdoor_final
rm(outdoor_final)
gc()

suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
  library(ggplot2)
  library(gganimate)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(av)
  library(viridis)
  library(scales)
  library(patchwork)
})

setDT(outdoor)

# ------------------------------------------------------------
# 0) HELPERS
# ------------------------------------------------------------
rename_if_present <- function(dt, old, new) {
  if (old %in% names(dt)) setnames(dt, old, new)
}

q_num <- function(x, p) {
  as.numeric(quantile(x, p, names = FALSE, type = 7, na.rm = TRUE))
}

safe_div_pct <- function(num, den) {
  ifelse(den > 0, round(100 * num / den, 2), NA_real_)
}

mode_value <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_character_)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# ------------------------------------------------------------
# 1) STANDARDIZE / CHOOSE COORDINATES
# ------------------------------------------------------------
rename_if_present(outdoor, "LONGITUDE(deg)", "pa_longitude")
rename_if_present(outdoor, "LATITUDE(deg)",  "pa_latitude")

outdoor[, station := as.character(station)]

if (!inherits(outdoor$datetime_utc, "POSIXct")) {
  outdoor[, datetime_utc := as.POSIXct(datetime_utc, tz = "UTC")]
} else {
  attr(outdoor$datetime_utc, "tzone") <- "UTC"
}

outdoor[, longitude := NA_real_]
outdoor[, latitude  := NA_real_]

if ("pa_longitude" %in% names(outdoor)) {
  outdoor[, longitude := pa_longitude]
}
if ("pa_latitude" %in% names(outdoor)) {
  outdoor[, latitude := pa_latitude]
}
if ("gee_longitude" %in% names(outdoor)) {
  outdoor[is.na(longitude), longitude := gee_longitude]
}
if ("gee_latitude" %in% names(outdoor)) {
  outdoor[is.na(latitude), latitude := gee_latitude]
}

candidate_unit_cols <- c(
  "country", "Country", "country_name",
  "state", "State", "admin1", "admin_name",
  "subregion", "city"
)
unit_col <- intersect(candidate_unit_cols, names(outdoor))[1]
if (length(unit_col) == 0) unit_col <- NA_character_

base_keep <- c("station", "datetime_utc", "region", "longitude", "latitude", "pm25_cf1_avg")
opt_keep <- c()
if (!is.na(unit_col)) opt_keep <- c(opt_keep, unit_col)
if ("urban_rural" %in% names(outdoor)) opt_keep <- c(opt_keep, "urban_rural")

keep_cols <- unique(c(base_keep, opt_keep))

outdoor_geo <- outdoor[
  !is.na(station) &
    !is.na(datetime_utc) &
    !is.na(region) &
    is.finite(longitude) &
    is.finite(latitude),
  ..keep_cols
]

if (!is.na(unit_col)) {
  setnames(outdoor_geo, unit_col, "unit_name")
} else {
  outdoor_geo[, unit_name := NA_character_]
}

if (!"urban_rural" %in% names(outdoor_geo)) {
  outdoor_geo[, urban_rural := NA_character_]
}

cat("Rows with valid coords:", format(nrow(outdoor_geo), big.mark = ","), "\n")
cat("Stations with valid coords:", format(uniqueN(outdoor_geo$station), big.mark = ","), "\n")
#452,143,148 
cat("Unit column used:", ifelse(is.na(unit_col), "NONE", unit_col), "\n")
#36,328
rm(outdoor)
gc()

# ------------------------------------------------------------
# 2) BUILD MONTHLY STATION POSITIONS
# ------------------------------------------------------------
outdoor_geo[, ym := floor_date(datetime_utc, "month")]

sensor_monthly <- outdoor_geo[
  ,
  .(
    lon = median(longitude, na.rm = TRUE),
    lat = median(latitude,  na.rm = TRUE),
    n_rows_month = .N,
    first_seen_month = min(datetime_utc),
    last_seen_month  = max(datetime_utc)
  ),
  by = .(station, region, ym)
]

setorder(sensor_monthly, region, station, ym)

cat("Station-month rows:", format(nrow(sensor_monthly), big.mark = ","), "\n")
#929,765 
cat("Unique stations in station-month table:", format(uniqueN(sensor_monthly$station), big.mark = ","), "\n")
#36,328 

# ------------------------------------------------------------
# 3) STATION ACTIVE MONTHS
# ------------------------------------------------------------
station_months <- sensor_monthly[
  ,
  .(
    n_months = uniqueN(ym),
    first_month = min(ym),
    last_month  = max(ym)
  ),
  by = .(station, region)
]

months_summary_global <- station_months[, .(
  stations = .N,
  min_months = as.numeric(min(n_months)),
  p25 = q_num(n_months, 0.25),
  median_months = as.numeric(median(n_months)),
  mean_months = as.numeric(mean(n_months)),
  p75 = q_num(n_months, 0.75),
  p90 = q_num(n_months, 0.90),
  p95 = q_num(n_months, 0.95),
  max_months = as.numeric(max(n_months))
)]

months_summary_region <- station_months[, .(
  stations = .N,
  min_months = as.numeric(min(n_months)),
  p25 = q_num(n_months, 0.25),
  median_months = as.numeric(median(n_months)),
  mean_months = as.numeric(mean(n_months)),
  p75 = q_num(n_months, 0.75),
  p90 = q_num(n_months, 0.90),
  p95 = q_num(n_months, 0.95),
  max_months = as.numeric(max(n_months))
), by = region][order(-stations)]

cat("\nGLOBAL MONTHS ACTIVE SUMMARY:\n")
print(months_summary_global)
#   stations min_months   p25 median_months mean_months   p75   p90   p95

#1:    36346          1     7            20    25.58094    41    58    65
#   max_months
#1:         92
cat("\nREGIONAL MONTHS ACTIVE SUMMARY:\n")
print(months_summary_region)
#         region stations min_months   p25 median_months mean_months   p75
#1:        CONUS    33926          1     8            21    26.47332    43
#2: SouthAmerica     1008          1     4            11    14.74603    20
#3:        India      783          1     3             5    10.13665    15
#4:       Africa      629          1     4             9    14.03816    20
#     p90   p95 max_months
#1:    59    66         92
#2:    32    46         84
#3:    27    36         92
#4:    34    44         82
cat("\nLIFETIME BREAKDOWN (% of stations):\n")
print(
  station_months[, .(
    pct_1_month = 100 * mean(n_months == 1),
    pct_lt_6    = 100 * mean(n_months < 6),
    pct_ge_12   = 100 * mean(n_months >= 12),
    pct_ge_24   = 100 * mean(n_months >= 24)
  )]
)
#   pct_1_month pct_lt_6 pct_ge_12 pct_ge_24
#    5.742035 21.49618  64.63985  44.84675
   
# ------------------------------------------------------------
# 4) EVALUATE WHETHER A STATION EVER MOVED AT ALL
# ------------------------------------------------------------
move_thr <- 0.01  # ~1.1 km latitude

station_move <- sensor_monthly[
  ,
  .(
    n_months = .N,
    first_month = min(ym),
    last_month  = max(ym),
    lon_min = min(lon, na.rm = TRUE),
    lon_max = max(lon, na.rm = TRUE),
    lat_min = min(lat, na.rm = TRUE),
    lat_max = max(lat, na.rm = TRUE),
    max_lon_shift = max(lon, na.rm = TRUE) - min(lon, na.rm = TRUE),
    max_lat_shift = max(lat, na.rm = TRUE) - min(lat, na.rm = TRUE),
    first_lon = lon[1],
    first_lat = lat[1],
    last_lon  = lon[.N],
    last_lat  = lat[.N]
  ),
  by = .(station, region)
]

station_move[, max_shift_km_approx :=
  sqrt(
    (max_lat_shift * 111)^2 +
      (max_lon_shift * 111 * cospi((lat_min + lat_max) / 360))^2
  )
]

station_move[, movement_assessable := n_months >= 2]
station_move[, moved_any := movement_assessable &
               ((abs(max_lon_shift) > move_thr) | (abs(max_lat_shift) > move_thr))]

station_move[, movement_status := fifelse(
  n_months == 1, "Only 1 month observed",
  fifelse(moved_any, "Moved", "Did not move")
)]

station_move[, movement_status := factor(
  movement_status,
  levels = c("Moved", "Did not move", "Only 1 month observed")
)]

cat("\nMOVEMENT COUNTS:\n")
cat("Stations total:", format(nrow(station_move), big.mark = ","), "\n")
#36,346 
cat("Stations assessable (>=2 months):", format(sum(station_move$movement_assessable), big.mark = ","), "\n")
#Stations assessable (>=2 months): 34,259
cat("Stations moved at all:", format(sum(station_move$moved_any), big.mark = ","), "\n")
#Stations moved at all: 1,812 
cat("Stations with only 1 month:", format(sum(station_move$n_months == 1), big.mark = ","), "\n")
#Stations with only 1 month: 2,087 

movement_summary_overall <- station_move[, .(
  stations_total = .N,
  stations_assessable = as.numeric(sum(movement_assessable)),
  moved_stations = as.numeric(sum(moved_any)),
  pct_moved_among_all = round(100 * mean(moved_any), 2),
  pct_moved_among_assessable = safe_div_pct(sum(moved_any), sum(movement_assessable)),
  stations_one_month_only = as.numeric(sum(n_months == 1)),
  pct_one_month_only = round(100 * mean(n_months == 1), 2),
  median_months_observed = as.numeric(median(n_months)),
  median_shift_km = as.numeric(median(max_shift_km_approx, na.rm = TRUE)),
  p95_shift_km = q_num(max_shift_km_approx, 0.95)
)]

movement_summary_region <- station_move[, .(
  stations_total = .N,
  stations_assessable = as.numeric(sum(movement_assessable)),
  moved_stations = as.numeric(sum(moved_any)),
  pct_moved_among_all = round(100 * mean(moved_any), 2),
  pct_moved_among_assessable = safe_div_pct(sum(moved_any), sum(movement_assessable)),
  stations_one_month_only = as.numeric(sum(n_months == 1)),
  pct_one_month_only = round(100 * mean(n_months == 1), 2),
  median_months_observed = as.numeric(median(n_months)),
  median_shift_km = as.numeric(median(max_shift_km_approx, na.rm = TRUE)),
  p95_shift_km = q_num(max_shift_km_approx, 0.95)
), by = region][order(-stations_total)]

cat("\nMOVEMENT SUMMARY OVERALL:\n")
print(movement_summary_overall)
#   stations_total stations_assessable moved_stations pct_moved_among_all
#          36346               34259           1812                4.99
#   pct_moved_among_assessable stations_one_month_only pct_one_month_only
#                     5.29                    2087               5.74
#   median_months_observed median_shift_km p95_shift_km
#                    20               0     1.026174
cat("\nMOVEMENT SUMMARY BY REGION:\n")
print(movement_summary_region)
#         region stations_total stations_assessable moved_stations
#        CONUS          33926               32025           1557
# SouthAmerica           1008                 934            115
#        India            783                 734             79
#      Africa            629                 566             61
#   pct_moved_among_all pct_moved_among_assessable stations_one_month_only
#                4.59                       4.86                    1901
#               11.41                      12.31                      74
#               10.09                      10.76                      49
#                9.70                      10.78                      63
#   pct_one_month_only median_months_observed median_shift_km p95_shift_km
#               5.60                     21               0    0.5280846
#               7.34                     11               0   13.2716777
#              6.26                      5               0   17.0144356
#             10.02                      9               0    8.6821101

# ------------------------------------------------------------
# 4A) COMBINED REGION SUMMARY TABLE
# ------------------------------------------------------------
region_summary_table <- merge(
  months_summary_region[
    , .(
      region,
      n_sensors = stations,
      months_min = min_months,
      months_p25 = p25,
      months_median = median_months,
      months_mean = round(mean_months, 2),
      months_p75 = p75,
      months_p95 = p95,
      months_max = max_months
    )
  ],
  movement_summary_region[
    , .(
      region,
      sensors_assessable_for_movement = stations_assessable,
      moved_sensors = moved_stations,
      pct_moved_among_all,
      pct_moved_among_assessable,
      one_month_only_sensors = stations_one_month_only,
      pct_one_month_only,
      movement_median_shift_km = round(median_shift_km, 3),
      movement_p95_shift_km = round(p95_shift_km, 3)
    )
  ],
  by = "region",
  all = TRUE
)

setorder(region_summary_table, -n_sensors)

cat("\nREGION SUMMARY TABLE:\n")
print(region_summary_table)

# region n_sensors months_min months_p25 months_median months_mean
#        CONUS     33926          1          8            21       26.47
# SouthAmerica      1008          1          4            11       14.75
#        India       783          1          3             5       10.14
#      Africa       629          1          4             9       14.04
 #  months_p75 months_p95 months_max sensors_assessable_for_movement
#         43         66         92                           32025
#         20         46         84                             934
#         15         36         92                             734
#         20         44         82                             566
#   moved_sensors pct_moved_among_all pct_moved_among_assessable
#         1557                4.59                       4.86
#           115               11.41                      12.31
#            79               10.09                      10.76
#           61                9.70                      10.78
#   one_month_only_sensors pct_one_month_only movement_median_shift_km
#                   1901               5.60                        0
#                     74               7.34                        0
#                     49               6.26                        0
#                    63              10.02                        0
#   movement_p95_shift_km
#                 0.528
#                13.272
#                17.014
#                8.682

# ------------------------------------------------------------
# 5) STATION-LEVEL POINT TABLES FOR STATIC MAPS
# ------------------------------------------------------------
station_points_latest <- sensor_monthly[
  order(region, station, ym),
  .SD[.N],
  by = .(station, region)
]

station_points_latest <- merge(
  station_points_latest,
  station_move[, .(
    station, region, n_months, moved_any, movement_assessable,
    movement_status, max_shift_km_approx
  )],
  by = c("station", "region"),
  all.x = TRUE
)

station_points_latest <- station_points_latest[!is.na(movement_status)]

# ------------------------------------------------------------
# 6) ANIMATION DATA
# ------------------------------------------------------------
anim_dt <- merge(
  sensor_monthly,
  station_move[, .(
    station, region, n_months, moved_any, movement_assessable,
    movement_status, max_shift_km_approx
  )],
  by = c("station", "region"),
  all.x = TRUE
)

anim_dt <- anim_dt[!is.na(movement_status)]

cat("Animation rows:", format(nrow(anim_dt), big.mark = ","), "\n")
#929,765 
cat("Animation stations:", format(uniqueN(anim_dt$station), big.mark = ","), "\n")
#36,328

# ------------------------------------------------------------
# 7) BASEMAP + EXTENT
# ------------------------------------------------------------
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

x_rng <- range(sensor_monthly$lon, na.rm = TRUE)
y_rng <- range(sensor_monthly$lat, na.rm = TRUE)

x_pad <- diff(x_rng) * 0.05
y_pad <- diff(y_rng) * 0.05

if (!is.finite(x_pad) || x_pad == 0) x_pad <- 1
if (!is.finite(y_pad) || y_pad == 0) y_pad <- 1

# ------------------------------------------------------------
# 8) MONTHLY ANIMATION
# ------------------------------------------------------------
p_anim <- ggplot() +
  geom_sf(data = world, linewidth = 0.15, fill = "grey98", color = "grey70") +
  geom_point(
    data = anim_dt,
    aes(x = lon, y = lat, color = movement_status),
    alpha = 0.65,
    size = 0.8
  ) +
  scale_color_manual(
    values = c(
      "Moved" = "red",
      "Did not move" = "steelblue",
      "Only 1 month observed" = "goldenrod"
    ),
    drop = FALSE,
    name = "Movement status"
  ) +
  coord_sf(
    xlim = c(x_rng[1] - x_pad, x_rng[2] + x_pad),
    ylim = c(y_rng[1] - y_pad, y_rng[2] + y_pad),
    expand = FALSE
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title = "PurpleAir sensor locations over time",
    subtitle = "Month: {frame_time}",
    x = NULL,
    y = NULL,
    caption = paste0("Moved = any station with > ", move_thr, " degree change across observed months")
  ) +
  transition_time(ym) +
  ease_aes("linear")

out_mp4 <- "/Users/priyanka/Downloads/PurpleAir_EPA/station_moves_monthly.mp4"

animate(
  p_anim,
  fps = 8,
  width = 1200,
  height = 650,
  renderer = av_renderer(out_mp4)
)

cat("Wrote animation:", out_mp4, "\n")

# ------------------------------------------------------------
# 9) STATIC MAP: MOVEMENT STATUS
# ------------------------------------------------------------
p_move <- ggplot() +
  geom_sf(data = world, linewidth = 0.15, fill = "grey98", color = "grey70") +
  geom_point(
    data = station_points_latest,
    aes(x = lon, y = lat, color = movement_status),
    alpha = 0.70,
    size = 0.8
  ) +
  scale_color_manual(
    values = c(
      "Moved" = "red",
      "Did not move" = "steelblue",
      "Only 1 month observed" = "goldenrod"
    ),
    drop = FALSE,
    name = "Movement status"
  ) +
  coord_sf(
    xlim = c(x_rng[1] - x_pad, x_rng[2] + x_pad),
    ylim = c(y_rng[1] - y_pad, y_rng[2] + y_pad),
    expand = FALSE
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title = "PurpleAir stations colored by movement status",
    subtitle = "Using latest observed location per station",
    x = NULL, y = NULL
  )

ggsave(
  "/Users/priyanka/Downloads/PurpleAir_EPA/station_movement_status_map.png",
  p_move, width = 12, height = 6.5, dpi = 300
)

# ------------------------------------------------------------
# 10) STATIC MAP: MONTHS ACTIVE
# ------------------------------------------------------------
p_months <- ggplot() +
  geom_sf(data = world, linewidth = 0.15, fill = "grey98", color = "grey70") +
  geom_point(
    data = station_points_latest,
    aes(x = lon, y = lat, color = n_months),
    alpha = 0.75,
    size = 0.8
  ) +
  scale_color_viridis_c(
    option = "plasma",
    trans = "sqrt",
    name = "Months active"
  ) +
  coord_sf(
    xlim = c(x_rng[1] - x_pad, x_rng[2] + x_pad),
    ylim = c(y_rng[1] - y_pad, y_rng[2] + y_pad),
    expand = FALSE
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title = "PurpleAir stations colored by months active",
    subtitle = "Using latest observed location per station",
    x = NULL, y = NULL
  )

ggsave(
  "/Users/priyanka/Downloads/PurpleAir_EPA/station_months_active_map.png",
  p_months, width = 12, height = 6.5, dpi = 300
)

# ------------------------------------------------------------
# 11) SENSOR COUNT-DENSITY MAP
# ------------------------------------------------------------
p_density <- ggplot() +
  geom_sf(data = world, linewidth = 0.15, fill = "grey98", color = "grey70") +
  stat_bin_2d(
    data = station_points_latest,
    aes(x = lon, y = lat, fill = after_stat(count)),
    bins = 180,
    alpha = 0.85
  ) +
  scale_fill_viridis_c(
    option = "magma",
    trans = "log10",
    name = "Station count"
  ) +
  coord_sf(
    xlim = c(x_rng[1] - x_pad, x_rng[2] + x_pad),
    ylim = c(y_rng[1] - y_pad, y_rng[2] + y_pad),
    expand = FALSE
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title = "PurpleAir station density",
    subtitle = "2D binned count of latest observed station locations",
    x = NULL, y = NULL
  )

ggsave(
  "/Users/priyanka/Downloads/PurpleAir_EPA/station_density_map.png",
  p_density, width = 12, height = 6.5, dpi = 300
)

# ------------------------------------------------------------
# 12) REGIONAL SENSOR COUNT-DENSITY MAPS
# ------------------------------------------------------------
p_density_region <- ggplot() +
  geom_sf(data = world, linewidth = 0.15, fill = "grey98", color = "grey70") +
  stat_bin_2d(
    data = station_points_latest,
    aes(x = lon, y = lat, fill = after_stat(count)),
    bins = 60,
    alpha = 0.85
  ) +
  scale_fill_viridis_c(
    option = "magma",
    trans = "log10",
    name = "Station count"
  ) +
  coord_sf(expand = FALSE) +
  facet_wrap(~region) +
  theme_minimal(base_size = 11) +
  labs(
    title = "PurpleAir station density by region",
    subtitle = "2D binned count of latest observed station locations",
    x = NULL, y = NULL
  )

ggsave(
  "/Users/priyanka/Downloads/PurpleAir_EPA/station_density_map_by_region.png",
  p_density_region, width = 14, height = 8, dpi = 300
)

# ------------------------------------------------------------
# 13) HISTOGRAM: MONTHS ACTIVE
# ------------------------------------------------------------
p_hist_months <- ggplot(station_move, aes(x = n_months)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white", boundary = 0) +
  scale_x_continuous(
    breaks = pretty_breaks(n = 12),
    labels = label_number(accuracy = 1)
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Distribution of months active per station",
    x = "Months active",
    y = "Number of stations"
  )

ggsave(
  "/Users/priyanka/Downloads/PurpleAir_EPA/station_months_active_histogram.png",
  p_hist_months, width = 9, height = 5.5, dpi = 300
)

p_hist_months_region <- ggplot(station_move, aes(x = n_months)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white", boundary = 0) +
  scale_x_continuous(
    breaks = pretty_breaks(n = 8),
    labels = label_number(accuracy = 1)
  ) +
  facet_wrap(~region, scales = "free_y") +
  theme_minimal(base_size = 11) +
  labs(
    title = "Distribution of months active per station by region",
    x = "Months active",
    y = "Number of stations"
  )

ggsave(
  "/Users/priyanka/Downloads/PurpleAir_EPA/station_months_active_histogram_by_region.png",
  p_hist_months_region, width = 13, height = 8, dpi = 300
)

# ------------------------------------------------------------
# 14) HISTOGRAM: MOVEMENT DISTANCE
# ------------------------------------------------------------
shift_hist_dt <- station_move[
  movement_assessable == TRUE & is.finite(max_shift_km_approx)
]

shift_positive_dt <- shift_hist_dt[max_shift_km_approx > 0]

if (nrow(shift_positive_dt) > 0) {
  x_upper <- q_num(shift_positive_dt$max_shift_km_approx, 0.995)
  shift_plot_dt <- shift_positive_dt[max_shift_km_approx <= x_upper]

  move_breaks <- c(0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100)

  p_hist_shift <- ggplot(shift_plot_dt, aes(x = max_shift_km_approx)) +
    geom_histogram(bins = 60, fill = "firebrick", color = "white") +
    scale_x_log10(
      breaks = move_breaks,
      labels = label_number(accuracy = 0.001, trim = TRUE)
    ) +
    scale_y_continuous(labels = comma_format()) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Distribution of approximate station movement distance",
      subtitle = paste0(
        "Stations observed in at least 2 months and with non-zero shift\n",
        "Zero-shift assessable stations excluded: ",
        comma(sum(shift_hist_dt$max_shift_km_approx == 0, na.rm = TRUE)),
        " | Plot capped at 99.5th percentile = ",
        round(x_upper, 2), " km"
      ),
      x = "Approximate max shift (km, log scale)",
      y = "Number of stations"
    )

  ggsave(
    "/Users/priyanka/Downloads/PurpleAir_EPA/station_movement_distance_histogram.png",
    p_hist_shift, width = 10, height = 6, dpi = 300
  )
}

# ------------------------------------------------------------
# 15) PAPER-STYLE DEPLOYMENT ANALYSES
# ------------------------------------------------------------
station_deploy <- outdoor_geo[
  ,
  .(
    first_seen = min(datetime_utc, na.rm = TRUE),
    last_seen  = max(datetime_utc, na.rm = TRUE),
    n_obs = .N,
    n_obs_months = uniqueN(ym),
    pm25_mean = mean(pm25_cf1_avg, na.rm = TRUE),
    pm25_min  = min(pm25_cf1_avg, na.rm = TRUE),
    pm25_max  = max(pm25_cf1_avg, na.rm = TRUE),
    unit_name = mode_value(unit_name),
    urban_rural = mode_value(urban_rural)
  ),
  by = .(station, region)
]

station_deploy[, deploy_days := as.numeric(difftime(last_seen, first_seen, units = "days"))]
station_deploy[, deploy_years := deploy_days / 365.25]

station_deploy[, longevity_group := fifelse(
  deploy_years < 1, "<1 year",
  fifelse(
    deploy_years < 2, "1-2 years",
    fifelse(
      deploy_years < 3, "2-3 years",
      ">3 years"
    )
  )
)]

station_deploy[, longevity_group := factor(
  longevity_group,
  levels = c("<1 year", "1-2 years", "2-3 years", ">3 years")
)]

cat("\nDEPLOYMENT TABLE STATIONS:", format(nrow(station_deploy), big.mark = ","), "\n")

deployment_global <- station_deploy[, .N, by = longevity_group][order(longevity_group)]
deployment_global[, pct := round(100 * N / sum(N), 2)]

deployment_region <- station_deploy[, .N, by = .(region, longevity_group)][order(region, longevity_group)]
deployment_region[, pct_within_region := round(100 * N / sum(N), 2), by = region]

cat("\nGLOBAL DEPLOYMENT LONGEVITY:\n")
print(deployment_global)

cat("\nREGIONAL DEPLOYMENT LONGEVITY:\n")
print(deployment_region)

# panel A: counts
p_deploy_region_counts <- ggplot(
  deployment_region,
  aes(x = region, y = N, fill = longevity_group)
) +
  geom_col() +
  theme_minimal(base_size = 12) +
  labs(
    title = "A) Sensor deployment longevity by region",
    subtitle = "Counts by tenure category",
    x = "Region",
    y = "Sensor count",
    fill = "Longevity group"
  )

# panel B: within-region proportions
p_deploy_region_prop <- ggplot(
  deployment_region,
  aes(x = region, y = N, fill = longevity_group)
) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_minimal(base_size = 12) +
  labs(
    title = "B) Sensor deployment longevity by region",
    subtitle = "Within-region distribution of station tenure",
    x = "Region",
    y = "Percent of sensors",
    fill = "Longevity group"
  )

# combined panel figure
library(patchwork)

p_deploy_region_counts2 <- p_deploy_region_counts +
  theme(legend.position = "right")

p_deploy_region_prop2 <- p_deploy_region_prop +
  theme(legend.position = "right")

p_deploy_combined <- (
  p_deploy_region_counts2 /
  p_deploy_region_prop2
) + plot_layout(guides = "collect", heights = c(1, 1))

ggsave(
  "/Users/priyanka/Downloads/PurpleAir_EPA/deployment_longevity_combined.png",
  p_deploy_combined,
  width = 9,
  height = 10,
  dpi = 300
)

# optional separate files too
ggsave(
  "/Users/priyanka/Downloads/PurpleAir_EPA/deployment_longevity_by_region.png",
  p_deploy_region_counts, width = 9, height = 5.5, dpi = 300
)

ggsave(
  "/Users/priyanka/Downloads/PurpleAir_EPA/deployment_longevity_by_region_proportion.png",
  p_deploy_region_prop, width = 9, height = 5.5, dpi = 300
)

# ------------------------------------------------------------
# 15A) REGION-BY-UNIT COUNTS
# ------------------------------------------------------------
if (any(!is.na(station_deploy$unit_name))) {

  deployment_unit <- station_deploy[
    !is.na(unit_name),
    .N,
    by = .(region, unit_name, longevity_group)
  ][order(region, unit_name, longevity_group)]

  deployment_unit_wide <- dcast(
    deployment_unit,
    region + unit_name ~ longevity_group,
    value.var = "N",
    fill = 0
  )

  for (nm in c("<1 year", "1-2 years", "2-3 years", ">3 years")) {
    if (!nm %in% names(deployment_unit_wide)) {
      deployment_unit_wide[, (nm) := 0]
    }
  }

  deployment_unit_wide[
    ,
    total := get("<1 year") + get("1-2 years") + get("2-3 years") + get(">3 years")
  ]

  setorder(deployment_unit_wide, region, -total)

  fwrite(
    deployment_unit_wide,
    "/Users/priyanka/Downloads/PurpleAir_EPA/deployment_longevity_by_region_unit.csv"
  )

  exemplar_units <- deployment_unit_wide[, .SD[1], by = region]

  exemplar_pm_summary <- station_deploy[
    exemplar_units[, .(region, unit_name)],
    on = .(region, unit_name)
  ][
    ,
    .(
      stations = .N,
      pm25_min = min(pm25_mean, na.rm = TRUE),
      pm25_mean = mean(pm25_mean, na.rm = TRUE),
      pm25_max = max(pm25_mean, na.rm = TRUE)
    ),
    by = .(region, unit_name, longevity_group)
  ][order(region, longevity_group)]

  cat("\nEXEMPLAR UNIT PM SUMMARY BY LONGEVITY:\n")
  print(exemplar_pm_summary)

  fwrite(
    exemplar_pm_summary,
    "/Users/priyanka/Downloads/PurpleAir_EPA/exemplar_unit_pm_summary_by_longevity.csv"
  )

  if (any(!is.na(station_deploy$urban_rural))) {
    exemplar_urban_pm <- station_deploy[
      exemplar_units[, .(region, unit_name)],
      on = .(region, unit_name)
    ][
      !is.na(urban_rural),
      .(
        stations = .N,
        pm25_min = min(pm25_mean, na.rm = TRUE),
        pm25_mean = mean(pm25_mean, na.rm = TRUE),
        pm25_max = max(pm25_mean, na.rm = TRUE)
      ),
      by = .(region, unit_name, urban_rural)
    ][order(region, urban_rural)]

    cat("\nEXEMPLAR UNIT PM SUMMARY BY URBAN/RURAL:\n")
    print(exemplar_urban_pm)

    fwrite(
      exemplar_urban_pm,
      "/Users/priyanka/Downloads/PurpleAir_EPA/exemplar_unit_pm_summary_by_urban_rural.csv"
    )
  }
}

# ------------------------------------------------------------
# 16) CSV OUTPUTS
# ------------------------------------------------------------
fwrite(
  station_move[order(-moved_any, -max_shift_km_approx)],
  "/Users/priyanka/Downloads/PurpleAir_EPA/station_movement_summary.csv"
)

fwrite(
  station_points_latest[order(-n_months)],
  "/Users/priyanka/Downloads/PurpleAir_EPA/station_latest_locations_with_months.csv"
)

fwrite(
  months_summary_global,
  "/Users/priyanka/Downloads/PurpleAir_EPA/station_months_active_summary_global.csv"
)

fwrite(
  months_summary_region,
  "/Users/priyanka/Downloads/PurpleAir_EPA/station_months_active_summary_by_region.csv"
)

fwrite(
  movement_summary_overall,
  "/Users/priyanka/Downloads/PurpleAir_EPA/station_movement_summary_global.csv"
)

fwrite(
  movement_summary_region,
  "/Users/priyanka/Downloads/PurpleAir_EPA/station_movement_summary_by_region.csv"
)

fwrite(
  deployment_global,
  "/Users/priyanka/Downloads/PurpleAir_EPA/deployment_longevity_global.csv"
)

fwrite(
  deployment_region,
  "/Users/priyanka/Downloads/PurpleAir_EPA/deployment_longevity_by_region.csv"
)

fwrite(
  region_summary_table,
  "/Users/priyanka/Downloads/PurpleAir_EPA/region_sensor_lifetime_movement_summary.csv"
)

# ------------------------------------------------------------
# 17) CONSOLE SUMMARY
# ------------------------------------------------------------
cat("\nFILES WRITTEN:\n")
cat("  Animation: /Users/priyanka/Downloads/PurpleAir_EPA/station_moves_monthly.mp4\n")
cat("  Map:       /Users/priyanka/Downloads/PurpleAir_EPA/station_movement_status_map.png\n")
cat("  Map:       /Users/priyanka/Downloads/PurpleAir_EPA/station_months_active_map.png\n")
cat("  Map:       /Users/priyanka/Downloads/PurpleAir_EPA/station_density_map.png\n")
cat("  Map:       /Users/priyanka/Downloads/PurpleAir_EPA/station_density_map_by_region.png\n")
cat("  Figure:    /Users/priyanka/Downloads/PurpleAir_EPA/deployment_longevity_by_region_combined.png\n")
cat("  Hist:      /Users/priyanka/Downloads/PurpleAir_EPA/station_months_active_histogram.png\n")
cat("  Hist:      /Users/priyanka/Downloads/PurpleAir_EPA/station_months_active_histogram_by_region.png\n")
cat("  Hist:      /Users/priyanka/Downloads/PurpleAir_EPA/station_movement_distance_histogram.png\n")
cat("  CSV:       /Users/priyanka/Downloads/PurpleAir_EPA/station_movement_summary.csv\n")
cat("  CSV:       /Users/priyanka/Downloads/PurpleAir_EPA/station_latest_locations_with_months.csv\n")
cat("  CSV:       /Users/priyanka/Downloads/PurpleAir_EPA/station_months_active_summary_global.csv\n")
cat("  CSV:       /Users/priyanka/Downloads/PurpleAir_EPA/station_months_active_summary_by_region.csv\n")
cat("  CSV:       /Users/priyanka/Downloads/PurpleAir_EPA/station_movement_summary_global.csv\n")
cat("  CSV:       /Users/priyanka/Downloads/PurpleAir_EPA/station_movement_summary_by_region.csv\n")
cat("  CSV:       /Users/priyanka/Downloads/PurpleAir_EPA/deployment_longevity_global.csv\n")
cat("  CSV:       /Users/priyanka/Downloads/PurpleAir_EPA/deployment_longevity_by_region.csv\n")
cat("  CSV:       /Users/priyanka/Downloads/PurpleAir_EPA/region_sensor_lifetime_movement_summary.csv\n")
cat("  CSV:       /Users/priyanka/Downloads/PurpleAir_EPA/deployment_longevity_by_region_unit.csv (if unit_name available)\n")
cat("  CSV:       /Users/priyanka/Downloads/PurpleAir_EPA/exemplar_unit_pm_summary_by_longevity.csv (if unit_name available)\n")
cat("  CSV:       /Users/priyanka/Downloads/PurpleAir_EPA/exemplar_unit_pm_summary_by_urban_rural.csv (if urban_rural available)\n")

# ---- source chunk #25 (Rmd lines 7315-7623) ------------------------

# ============================================================
# TABLE S1 + TABLE S2 + MANUSCRIPT NUMBERS
# Paste at end of summary/climate workflow
# ============================================================

library(data.table)

OUTDIR <- "/Users/priyanka/Downloads/PurpleAir_EPA"
dir.create(file.path(OUTDIR, "tables"), showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Load saved outputs if needed
# ----------------------------

if (!exists("outdoor_final")) {
  load(file.path(OUTDIR, "outdoor_final.RData"))
}
setDT(outdoor_final)

transition_csv <- file.path(OUTDIR, "station_koppen_transition_summary.csv")
station_koppen_change <- fread(transition_csv)
station_koppen_change[, station := as.character(station)]

# ----------------------------
# Column names
# ----------------------------

a_col <- "pm25_cf_1_a_hourly(ug/m3)_pm25_cf_1_a"
b_col <- "pm25_cf_1_b_hourly(ug/m3)_pm25_cf_1_b"
tcol  <- "temperature_hourly(C)_temperature"
rhcol <- "humidity_hourly(%)_humidity"

# ----------------------------
# Known initial values from raw workflow
# ----------------------------

initial_total_rows <- 493787385
initial_region_rows <- data.table(
  region = c("Africa", "CONUS", "India", "SouthAmerica"),
  initial_rows = c(NA_integer_, NA_integer_, NA_integer_, NA_integer_)
)

initial_sensor_counts <- data.table(
  region = c("CONUS", "SouthAmerica", "India", "Africa"),
  initial_sensors = c(34071, 1013, 785, 639)
)

initial_total_sensors <- 36489

# ----------------------------
# Final observation and sensor counts
# ----------------------------

final_rows_region <- outdoor_final[, .(
  final_rows = .N,
  final_sensors = uniqueN(station)
), by = region][order(region)]

final_rows_overall <- outdoor_final[, .(
  final_rows = .N,
  final_sensors = uniqueN(station)
)]
final_rows_overall[, region := "Global"]

# ----------------------------
# Table S1: deployment/lifetime/movement
# Uses existing region_summary_table if available;
# otherwise rebuilds from saved CSV.
# ----------------------------

region_summary_path <- file.path(OUTDIR, "region_sensor_lifetime_movement_summary.csv")

if (!exists("region_summary_table")) {
  if (file.exists(region_summary_path)) {
    region_summary_table <- fread(region_summary_path)
  } else {
    stop("region_summary_table not found. Run the movement/lifetime chunk first.")
  }
}

table_s1_region <- copy(region_summary_table)

# Add global row if station_move exists
if (exists("station_move")) {
  table_s1_global <- station_move[, .(
    region = "Global",
    n_sensors = .N,
    months_min = min(n_months, na.rm = TRUE),
    months_p25 = as.numeric(quantile(n_months, 0.25, na.rm = TRUE)),
    months_median = median(n_months, na.rm = TRUE),
    months_mean = round(mean(n_months, na.rm = TRUE), 2),
    months_p75 = as.numeric(quantile(n_months, 0.75, na.rm = TRUE)),
    months_p95 = as.numeric(quantile(n_months, 0.95, na.rm = TRUE)),
    months_max = max(n_months, na.rm = TRUE),
    sensors_assessable_for_movement = sum(movement_assessable, na.rm = TRUE),
    moved_sensors = sum(moved_any, na.rm = TRUE),
    pct_moved_among_all = round(100 * mean(moved_any, na.rm = TRUE), 2),
    pct_moved_among_assessable =
      round(100 * sum(moved_any, na.rm = TRUE) / sum(movement_assessable, na.rm = TRUE), 2),
    one_month_only_sensors = sum(n_months == 1, na.rm = TRUE),
    pct_one_month_only = round(100 * mean(n_months == 1, na.rm = TRUE), 2),
    movement_median_shift_km = round(median(max_shift_km_approx, na.rm = TRUE), 3),
    movement_p95_shift_km = round(as.numeric(quantile(max_shift_km_approx, 0.95, na.rm = TRUE)), 3)
  )]
  
  table_s1 <- rbindlist(list(table_s1_global, table_s1_region), fill = TRUE)
} else {
  table_s1 <- table_s1_region
}

fwrite(table_s1, file.path(OUTDIR, "tables", "Table_S1_sensor_lifetime_movement_summary.csv"))

# ----------------------------
# Table S2: Köppen climate zone distribution
# ----------------------------

station_koppen_single <- station_koppen_change[
  !is.na(first_koppen) & n_koppen_codes == 1,
  .(
    station,
    region,
    koppen_code = first_koppen,
    koppen_desc = first_koppen_desc
  )
]

table_s2_counts <- station_koppen_single[
  ,
  .(n_sensors = .N),
  by = .(koppen_code, koppen_desc, region)
]

table_s2_wide <- dcast(
  table_s2_counts,
  koppen_code + koppen_desc ~ region,
  value.var = "n_sensors",
  fill = 0
)

region_cols <- intersect(c("Africa", "CONUS", "India", "SouthAmerica"), names(table_s2_wide))

table_s2_wide[, Global := rowSums(.SD), .SDcols = region_cols]
setcolorder(table_s2_wide, c("koppen_code", "koppen_desc", "Africa", "CONUS", "Global", "India", "SouthAmerica"))

excluded_by_region <- station_koppen_change[
  changed_koppen_zone == TRUE,
  .N,
  by = region
]

excluded_row <- data.table(
  koppen_code = "Excluded sensors",
  koppen_desc = "Stations changing Köppen zone during study period; excluded from climate-zone summaries",
  Africa = excluded_by_region[region == "Africa", N],
  CONUS = excluded_by_region[region == "CONUS", N],
  Global = sum(excluded_by_region$N),
  India = excluded_by_region[region == "India", N],
  SouthAmerica = excluded_by_region[region == "SouthAmerica", N]
)

excluded_row[is.na(excluded_row)] <- 0

table_s2 <- rbindlist(list(table_s2_wide, excluded_row), fill = TRUE)
setorder(table_s2, koppen_code)

fwrite(table_s2, file.path(OUTDIR, "tables", "Table_S2_koppen_climate_zone_distribution.csv"))

# ----------------------------
# Climate assignment numbers
# ----------------------------

n_final_sensors <- uniqueN(outdoor_final$station)

n_assigned_any <- station_koppen_change[!is.na(first_koppen), uniqueN(station)]
n_unassigned <- n_final_sensors - n_assigned_any

n_changed <- station_koppen_change[changed_koppen_zone == TRUE, uniqueN(station)]
pct_changed <- 100 * n_changed / n_final_sensors

changed_region <- station_koppen_change[
  changed_koppen_zone == TRUE,
  .(changed_sensors = uniqueN(station)),
  by = region
]

final_sensors_region <- outdoor_final[, .(final_sensors = uniqueN(station)), by = region]

changed_region <- merge(final_sensors_region, changed_region, by = "region", all.x = TRUE)
changed_region[is.na(changed_sensors), changed_sensors := 0]
changed_region[, pct_changed := 100 * changed_sensors / final_sensors]

# ----------------------------
# PM >1500 removal numbers
# If outdoor before PM filtering exists, use it.
# Otherwise use known values from your prior printed output.
# ----------------------------

pm_before_rows_region <- data.table(
  region = c("Africa", "CONUS", "India", "SouthAmerica"),
  before_pm_rows = c(3255433, 452053816, 2865850, 6673543)
)

pm_after_rows_region <- final_rows_region[, .(
  region,
  after_pm_rows = final_rows
)]

pm_filter_region <- merge(pm_before_rows_region, pm_after_rows_region, by = "region")
pm_filter_region[, n_removed_pm := before_pm_rows - after_pm_rows]
pm_filter_region[, pct_removed_pm := 100 * n_removed_pm / before_pm_rows]

pm_filter_overall <- pm_filter_region[, .(
  before_pm_rows = sum(before_pm_rows),
  after_pm_rows = sum(after_pm_rows),
  n_removed_pm = sum(n_removed_pm)
)]
pm_filter_overall[, pct_removed_pm := 100 * n_removed_pm / before_pm_rows]
pm_filter_overall[, region := "Global"]

pm_filter_table <- rbindlist(list(pm_filter_overall, pm_filter_region), fill = TRUE)
setcolorder(pm_filter_table, c("region", "before_pm_rows", "n_removed_pm", "pct_removed_pm", "after_pm_rows"))

fwrite(pm_filter_table, file.path(OUTDIR, "tables", "pm1500_filtering_numbers.csv"))

# ----------------------------
# Manuscript numbers summary
# ----------------------------

manuscript_numbers <- list(
  initial_total_rows = initial_total_rows,
  initial_total_sensors = initial_total_sensors,
  final_rows = nrow(outdoor_final),
  final_sensors = n_final_sensors,
  koppen_assigned_any = n_assigned_any,
  koppen_assigned_any_pct = 100 * n_assigned_any / n_final_sensors,
  koppen_unassigned = n_unassigned,
  koppen_changed_zone = n_changed,
  koppen_changed_zone_pct = pct_changed
)

manuscript_numbers_dt <- rbindlist(lapply(names(manuscript_numbers), function(x) {
  data.table(metric = x, value = manuscript_numbers[[x]])
}))

fwrite(manuscript_numbers_dt, file.path(OUTDIR, "tables", "manuscript_key_numbers_overall.csv"))
fwrite(changed_region, file.path(OUTDIR, "tables", "koppen_changed_zone_by_region.csv"))
fwrite(final_rows_region, file.path(OUTDIR, "tables", "final_rows_sensors_by_region.csv"))

# ----------------------------
# Create paste-ready manuscript text
# ----------------------------

fmt_int <- function(x) format(as.integer(round(x)), big.mark = ",", scientific = FALSE)
fmt_pct <- function(x, digits = 1) sprintf(paste0("%.", digits, "f"), x)

pm_global <- pm_filter_table[region == "Global"]
pm_conus  <- pm_filter_table[region == "CONUS"]
pm_india  <- pm_filter_table[region == "India"]
pm_africa <- pm_filter_table[region == "Africa"]
pm_sa     <- pm_filter_table[region == "SouthAmerica"]

final_conus  <- final_rows_region[region == "CONUS"]
final_india  <- final_rows_region[region == "India"]
final_africa <- final_rows_region[region == "Africa"]
final_sa     <- final_rows_region[region == "SouthAmerica"]

chg_conus  <- changed_region[region == "CONUS"]
chg_india  <- changed_region[region == "India"]
chg_africa <- changed_region[region == "Africa"]
chg_sa     <- changed_region[region == "SouthAmerica"]

paste_ready <- paste0(
"Hourly PM2.5 observations were additionally removed when either cf_1 channel A or B exceeded 1,500 µg/m3. ",
"This filtering removed ", fmt_int(pm_global$n_removed_pm), " observations (",
fmt_pct(pm_global$pct_removed_pm, 1), "% overall; ",
fmt_pct(pm_conus$pct_removed_pm, 1), "% in CONUS, ",
fmt_pct(pm_india$pct_removed_pm, 1), "% in India, ",
fmt_pct(pm_africa$pct_removed_pm, 1), "% in Africa, and ",
fmt_pct(pm_sa$pct_removed_pm, 1), "% in South America), leaving a final dataset comprising ",
fmt_int(nrow(outdoor_final)), " hourly observations from ", fmt_int(n_final_sensors), " sensors (",
fmt_int(final_conus$final_sensors), " in CONUS, ",
fmt_int(final_india$final_sensors), " in India, ",
fmt_int(final_africa$final_sensors), " in Africa, and ",
fmt_int(final_sa$final_sensors), " in South America).\n\n",
"Sensors were assigned Köppen-Geiger climate classifications using the 1991–2020 global climate maps. Of the ",
fmt_int(n_final_sensors), " sensors, ", fmt_int(n_assigned_any), " (",
fmt_pct(100 * n_assigned_any / n_final_sensors, 2), "%) were successfully assigned at least one Köppen-Geiger climate zone. ",
"A total of ", fmt_int(n_changed), " sensors (", fmt_pct(pct_changed, 2),
"%) were associated with multiple climate zones because of location changes over time, including ",
fmt_int(chg_conus$changed_sensors), " (", fmt_pct(chg_conus$pct_changed, 2), "%) in CONUS, ",
fmt_int(chg_sa$changed_sensors), " (", fmt_pct(chg_sa$pct_changed, 2), "%) in South America, ",
fmt_int(chg_india$changed_sensors), " (", fmt_pct(chg_india$pct_changed, 2), "%) in India, and ",
fmt_int(chg_africa$changed_sensors), " (", fmt_pct(chg_africa$pct_changed, 2), "%) in Africa. ",
"These relocated sensors were excluded from climate-zone-specific analyses."
)

writeLines(paste_ready, file.path(OUTDIR, "tables", "paste_ready_methods_numbers.txt"))

cat("\nSaved:\n")
cat(file.path(OUTDIR, "tables", "Table_S1_sensor_lifetime_movement_summary.csv"), "\n")
cat(file.path(OUTDIR, "tables", "Table_S2_koppen_climate_zone_distribution.csv"), "\n")
cat(file.path(OUTDIR, "tables", "pm1500_filtering_numbers.csv"), "\n")
cat(file.path(OUTDIR, "tables", "manuscript_key_numbers_overall.csv"), "\n")
cat(file.path(OUTDIR, "tables", "paste_ready_methods_numbers.txt"), "\n\n")

cat("Paste-ready text:\n\n")
cat(paste_ready, "\n")

# ---- source chunk #26 (Rmd lines 7626-8539) ------------------------

# ============================================================
# MEMORY-SAFE FULL-DATA SUMMARY + CORRELATION + PLOTTING
#
# GLOBAL/REGION ANALYSES FROM:
#   /Users/priyanka/Downloads/PurpleAir_EPA/outdoor_final.RData
#
# CLIMATE-ZONE ANALYSES FROM:
#   /Users/priyanka/Downloads/PurpleAir_EPA/outdoor_with_koppen_timevarying.RData
#
# Global/region analyses:
#   - no Köppen restriction
#   - no single-zone restriction
#
# Climate-zone analyses:
#   - valid Köppen code only
#   - restricted to stations with exactly one Köppen code
#
# Plotting:
#   - no lumping/capping into edge bins
#   - plot windows are implemented by filtering plotted points/bins
#   - summaries/correlations use all available finite data
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(viridisLite)
})

# ------------------------------------------------------------
# 0) PATHS / SETTINGS
# ------------------------------------------------------------

GLOBAL_REGION_INFILE <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_final.RData"
CLIMATE_INFILE <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_with_koppen_timevarying.RData"
KOPPEN_CHANGE_FILE <- "/Users/priyanka/Downloads/PurpleAir_EPA/station_koppen_transition_summary.csv"

OUTDIR <- "/Users/priyanka/Downloads/PurpleAir_EPA/full_data_summary_figures_split_inputs"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

USE_SINGLE_ZONE_ONLY_FOR_CLIMATE <- TRUE

# Plot windows only — these do NOT cap/lump values
PM_PLOT_MAX <- 1500

TEMP_RESID_MIN <- -40
TEMP_RESID_MAX <- 40

RH_RESID_MIN <- -100
RH_RESID_MAX <- 100

PM_RESID_MIN <- -1500
PM_RESID_MAX <- 1500

# 1D bin widths
PM_BIN_WIDTH <- 5
TEMP_BIN_WIDTH <- 0.25
RH_BIN_WIDTH <- 1

# 2D bin widths
PM_X_BIN <- 5
PM_Y_BIN <- 5

TEMP_X_BIN <- 0.5
TEMP_Y_BIN <- 0.25

RH_X_BIN <- 1
RH_Y_BIN <- 1

say <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste0(..., collapse = ""), "\n")
}

# ------------------------------------------------------------
# 1) COLUMN MAPPING
# ------------------------------------------------------------

a_col <- "pm25_cf_1_a_hourly(ug/m3)_pm25_cf_1_a"
b_col <- "pm25_cf_1_b_hourly(ug/m3)_pm25_cf_1_b"

pa_pm_avg_col <- "pm25_cf1_avg"

geos_col <- "geoscf_pm25"

pa_temp_col <- "temperature_hourly(C)_temperature"
pa_rh_col   <- "humidity_hourly(%)_humidity"

era5_temp_col <- "era5_t2m_c"
era5_rh_col   <- "era5_rh_pct"

koppen_col <- "koppen_code"

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
    stop("Could not find object named 'outdoor' or 'outdoor_final' in: ", infile)
  }

  setDT(x)
  x[]
}

check_required_cols <- function(outdoor, require_koppen = FALSE) {

  needed <- c(
    "station", "region",
    geos_col,
    pa_temp_col, pa_rh_col,
    era5_temp_col, era5_rh_col
  )

  if (!(pa_pm_avg_col %in% names(outdoor))) {
    needed <- c(needed, a_col, b_col)
  }

  if (require_koppen) {
    needed <- c(needed, koppen_col)
  }

  missing <- setdiff(needed, names(outdoor))

  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }

  invisible(TRUE)
}

make_analysis_dt <- function(outdoor, require_koppen = FALSE) {

  check_required_cols(outdoor, require_koppen = require_koppen)

  keep_cols <- c(
    "station", "region",
    geos_col,
    pa_temp_col, pa_rh_col,
    era5_temp_col, era5_rh_col
  )

  if (pa_pm_avg_col %in% names(outdoor)) {
    keep_cols <- c(keep_cols, pa_pm_avg_col)
  } else {
    keep_cols <- c(keep_cols, a_col, b_col)
  }

  if (koppen_col %in% names(outdoor)) {
    keep_cols <- c(keep_cols, koppen_col)
  }

  keep_cols <- unique(keep_cols)

  dt <- outdoor[, ..keep_cols]

  rm(outdoor)
  gc()

  dt[, station := as.character(station)]
  dt[, region := as.character(region)]

  if (koppen_col %in% names(dt)) {
    dt[, koppen_code := as.character(get(koppen_col))]
    if (koppen_col != "koppen_code") {
      dt[, (koppen_col) := NULL]
    }
  } else {
    dt[, koppen_code := NA_character_]
  }

  if (pa_pm_avg_col %in% names(dt)) {
    setnames(dt, pa_pm_avg_col, "pa_pm25")
    dt[, pa_pm25 := as.numeric(pa_pm25)]
  } else {
    dt[, pa_pm25 := (as.numeric(get(a_col)) + as.numeric(get(b_col))) / 2]
    dt[, c(a_col, b_col) := NULL]
  }

  setnames(dt, geos_col, "geos_pm")
  setnames(dt, pa_temp_col, "pa_temp")
  setnames(dt, pa_rh_col, "pa_rh")
  setnames(dt, era5_temp_col, "era5_temp")
  setnames(dt, era5_rh_col, "era5_rh")

  dt[, `:=`(
    geos_pm = as.numeric(geos_pm),
    pa_temp = as.numeric(pa_temp),
    pa_rh = as.numeric(pa_rh),
    era5_temp = as.numeric(era5_temp),
    era5_rh = as.numeric(era5_rh)
  )]

  dt[, pm_resid := pa_pm25 - geos_pm]
  dt[, temp_resid := pa_temp - era5_temp]
  dt[, rh_resid := pa_rh - era5_rh]

  # remove rows missing region only because region grouping cannot use NA
  dt <- dt[!is.na(region)]

  gc()
  dt[]
}

# ------------------------------------------------------------
# 3) SUMMARY + CORRELATION HELPERS
# ------------------------------------------------------------

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(NA_real_)
  cor(x[ok], y[ok], use = "complete.obs")
}

make_summary <- function(dt, group_var) {
  dt[, .(
    n_rows = .N,
    n_stations = uniqueN(station),

    mean_pa_pm25 = mean(pa_pm25, na.rm = TRUE),
    median_pa_pm25 = median(pa_pm25, na.rm = TRUE),
    sd_pa_pm25 = sd(pa_pm25, na.rm = TRUE),

    mean_geos_pm = mean(geos_pm, na.rm = TRUE),
    mean_pm_resid = mean(pm_resid, na.rm = TRUE),
    mean_abs_pm_resid = mean(abs(pm_resid), na.rm = TRUE),

    mean_pa_temp = mean(pa_temp, na.rm = TRUE),
    mean_era5_temp = mean(era5_temp, na.rm = TRUE),
    mean_temp_resid = mean(temp_resid, na.rm = TRUE),
    mean_abs_temp_resid = mean(abs(temp_resid), na.rm = TRUE),

    mean_pa_rh = mean(pa_rh, na.rm = TRUE),
    mean_era5_rh = mean(era5_rh, na.rm = TRUE),
    mean_rh_resid = mean(rh_resid, na.rm = TRUE),
    mean_abs_rh_resid = mean(abs(rh_resid), na.rm = TRUE)
  ), by = group_var]
}

make_cor <- function(dt, group_var) {
  dt[, .(
    n_rows = .N,
    n_stations = uniqueN(station),
    cor_pa_pm25_geos_pm = safe_cor(pa_pm25, geos_pm),
    cor_pa_temp_era5_temp = safe_cor(pa_temp, era5_temp),
    cor_pa_rh_era5_rh = safe_cor(pa_rh, era5_rh),
    cor_pa_temp_pa_rh = safe_cor(pa_temp, pa_rh),
    cor_era5_temp_era5_rh = safe_cor(era5_temp, era5_rh)
  ), by = group_var]
}

make_qa <- function(dt, label) {
  qa <- dt[, .(
    dataset = label,
    n_rows = .N,
    n_stations = uniqueN(station),
    pa_pm25_finite = sum(is.finite(pa_pm25)),
    geos_pm_finite = sum(is.finite(geos_pm)),
    pm_resid_finite = sum(is.finite(pm_resid)),
    pa_temp_finite = sum(is.finite(pa_temp)),
    era5_temp_finite = sum(is.finite(era5_temp)),
    temp_resid_finite = sum(is.finite(temp_resid)),
    pa_rh_finite = sum(is.finite(pa_rh)),
    era5_rh_finite = sum(is.finite(era5_rh)),
    rh_resid_finite = sum(is.finite(rh_resid))
  )]

  print(qa)
  qa[]
}

# ------------------------------------------------------------
# 4) BINNING HELPERS
# ------------------------------------------------------------

bin1d_counts <- function(x, by_group, width, x_min = NULL, x_max = NULL) {

  ok <- is.finite(x) & !is.na(by_group)

  if (!is.null(x_min)) ok <- ok & x >= x_min
  if (!is.null(x_max)) ok <- ok & x <= x_max

  if (!any(ok)) {
    return(data.table(group = character(), mid = numeric(), N = integer()))
  }

  x <- x[ok]
  by_group <- by_group[ok]

  if (is.null(x_min)) x_min <- floor(min(x, na.rm = TRUE) / width) * width
  if (is.null(x_max)) x_max <- ceiling(max(x, na.rm = TRUE) / width) * width

  bins <- floor((x - x_min) / width)
  mids <- x_min + bins * width + width / 2

  data.table(group = by_group, mid = mids)[
    mid >= x_min & mid <= x_max,
    .N,
    by = .(group, mid)
  ][]
}

bin2d_count <- function(x, y, by_group, x_width, y_width,
                        x_min = NULL, x_max = NULL,
                        y_min = NULL, y_max = NULL) {

  ok <- is.finite(x) & is.finite(y) & !is.na(by_group)

  if (!is.null(x_min)) ok <- ok & x >= x_min
  if (!is.null(x_max)) ok <- ok & x <= x_max
  if (!is.null(y_min)) ok <- ok & y >= y_min
  if (!is.null(y_max)) ok <- ok & y <= y_max

  if (!any(ok)) {
    return(data.table(group = character(), xmid = numeric(), ymid = numeric(), value = numeric()))
  }

  x <- x[ok]
  y <- y[ok]
  by_group <- by_group[ok]

  if (is.null(x_min)) x_min <- floor(min(x, na.rm = TRUE) / x_width) * x_width
  if (is.null(x_max)) x_max <- ceiling(max(x, na.rm = TRUE) / x_width) * x_width
  if (is.null(y_min)) y_min <- floor(min(y, na.rm = TRUE) / y_width) * y_width
  if (is.null(y_max)) y_max <- ceiling(max(y, na.rm = TRUE) / y_width) * y_width

  xb <- x_min + floor((x - x_min) / x_width) * x_width + x_width / 2
  yb <- y_min + floor((y - y_min) / y_width) * y_width + y_width / 2

  out <- data.table(group = by_group, xmid = xb, ymid = yb)[
    xmid >= x_min & xmid <= x_max &
      ymid >= y_min & ymid <= y_max,
    .N,
    by = .(group, xmid, ymid)
  ]

  setnames(out, "N", "value")
  out[]
}

# ------------------------------------------------------------
# 5) PLOT HELPERS
# ------------------------------------------------------------

save_density_plot <- function(bin_dt, title, xlab, outfile) {
  if (is.null(bin_dt) || nrow(bin_dt) == 0L) {
    message("Skipping empty density plot: ", outfile)
    return(invisible(NULL))
  }

  p <- ggplot(bin_dt, aes(mid, N, color = group, group = group)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~group, scales = "free_y") +
    theme_bw() +
    labs(title = title, x = xlab, y = "Count")

  ggsave(outfile, p, width = 12, height = 7, dpi = 300, bg = "white")
  invisible(p)
}

save_binned_heatmap <- function(dt2, title, xlab, ylab, outfile,
                                fill_lab = "Count",
                                palette = "viridis",
                                trans = "log10",
                                fixed_scales = FALSE,
                                add_identity = FALSE,
                                add_zero = FALSE) {
  if (is.null(dt2) || nrow(dt2) == 0L) {
    message("Skipping empty heatmap: ", outfile)
    return(invisible(NULL))
  }

  dt2 <- dt2[is.finite(xmid) & is.finite(ymid) & is.finite(value) & !is.na(group)]

  if (nrow(dt2) == 0L) {
    message("Skipping filtered-empty heatmap: ", outfile)
    return(invisible(NULL))
  }

  facet_scales <- if (fixed_scales) "fixed" else "free"

  p <- ggplot(dt2, aes(xmid, ymid, fill = value)) +
    geom_tile() +
    facet_wrap(~group, scales = facet_scales) +
    scale_fill_viridis_c(
      option = palette,
      trans = trans,
      na.value = "grey97"
    ) +
    theme_minimal(base_size = 12) +
    labs(title = title, x = xlab, y = ylab, fill = fill_lab) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey92", color = "grey55"),
      strip.text = element_text(face = "bold"),
      panel.border = element_rect(color = "grey65", fill = NA),
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )

  if (add_identity) {
    p <- p + geom_abline(intercept = 0, slope = 1, color = "white", linetype = 2, linewidth = 0.4)
  }

  if (add_zero) {
    p <- p + geom_hline(yintercept = 0, color = "white", linetype = 2, linewidth = 0.4)
  }

  ggsave(outfile, p, width = 12, height = 7, dpi = 300, bg = "white")
  invisible(p)
}

save_correlation_figures <- function(cor_dt, group_col, group_label, file_suffix, top_n = NULL) {

  setDT(cor_dt)

  corr_vars <- c(
    "cor_pa_pm25_geos_pm",
    "cor_pa_temp_era5_temp",
    "cor_pa_rh_era5_rh",
    "cor_pa_temp_pa_rh",
    "cor_era5_temp_era5_rh"
  )

  metric_labels <- c(
    cor_pa_pm25_geos_pm    = "PA PM2.5 vs GEOS-CF PM2.5",
    cor_pa_temp_era5_temp  = "PA temperature vs ERA5 temperature",
    cor_pa_rh_era5_rh      = "PA RH vs ERA5 RH",
    cor_pa_temp_pa_rh      = "PA temperature vs PA RH",
    cor_era5_temp_era5_rh  = "ERA5 temperature vs ERA5 RH"
  )

  if (!is.null(top_n)) {
    keep_groups <- cor_dt[order(-n_rows)][1:min(top_n, .N), get(group_col)]
    cor_dt <- cor_dt[get(group_col) %in% keep_groups]
  }

  corr_long <- melt(
    cor_dt,
    id.vars = c(group_col, "n_rows", "n_stations"),
    measure.vars = corr_vars,
    variable.name = "metric",
    value.name = "correlation"
  )

  setnames(corr_long, group_col, "group")
  corr_long[, metric_label := metric_labels[metric]]
  corr_long[, metric_label := factor(metric_label, levels = unname(metric_labels))]

  p_heat <- ggplot(
    corr_long,
    aes(x = metric_label, y = reorder(group, correlation, FUN = mean), fill = correlation)
  ) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", correlation)), size = 3.0) +
    scale_fill_gradient2(
      low = "#b2182b",
      mid = "white",
      high = "#2166ac",
      midpoint = 0,
      limits = c(-1, 1),
      oob = squish
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold")
    ) +
    labs(
      title = paste0("Correlations by ", group_label),
      x = NULL,
      y = NULL,
      fill = "r"
    )

  ggsave(
    file.path(OUTDIR, paste0("correlations_", file_suffix, "_heatmap.png")),
    p_heat,
    width = 10,
    height = ifelse(file_suffix == "by_region", 4.5, 7),
    dpi = 300,
    bg = "white"
  )

  p_dot <- ggplot(
    corr_long,
    aes(x = correlation, y = reorder(group, correlation, FUN = mean), color = metric_label)
  ) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
    geom_point(size = 2.7, alpha = 0.9, position = position_dodge(width = 0.5)) +
    scale_x_continuous(limits = c(-1, 1)) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    labs(
      title = paste0("Correlation coefficients by ", group_label),
      x = "Correlation coefficient (r)",
      y = NULL,
      color = NULL
    )

  ggsave(
    file.path(OUTDIR, paste0("correlations_", file_suffix, "_dotplot.png")),
    p_dot,
    width = 10,
    height = ifelse(file_suffix == "by_region", 5, 7),
    dpi = 300,
    bg = "white"
  )
}

# ------------------------------------------------------------
# 6) GLOBAL/REGION DATASET
# ------------------------------------------------------------

say("Loading global/region data from outdoor_final.RData ...")

outdoor_global <- load_outdoor_object(GLOBAL_REGION_INFILE)

say("Global/region rows loaded: ", format(nrow(outdoor_global), big.mark = ","))
say("Global/region stations loaded: ", format(uniqueN(as.character(outdoor_global$station)), big.mark = ","))

dt_global <- make_analysis_dt(outdoor_global, require_koppen = FALSE)

rm(outdoor_global)
gc()

qa_global <- make_qa(dt_global, "outdoor_final_global_region")
fwrite(qa_global, file.path(OUTDIR, "qa_global_region.csv"))

# ------------------------------------------------------------
# 7) GLOBAL/REGION SUMMARY + CORRELATION
# ------------------------------------------------------------

say("Writing global/region summary and correlation tables ...")

summary_by_region <- make_summary(dt_global, "region")[order(region)]

cor_by_region <- make_cor(dt_global, "region")[order(region)]

fwrite(summary_by_region, file.path(OUTDIR, "summary_stats_by_region.csv"))
fwrite(cor_by_region, file.path(OUTDIR, "correlations_by_region.csv"))

save_correlation_figures(
  cor_dt = cor_by_region,
  group_col = "region",
  group_label = "region",
  file_suffix = "by_region"
)

# ------------------------------------------------------------
# 8) REGION PLOTS FROM GLOBAL DATASET
# ------------------------------------------------------------

say("Building region plots from outdoor_final.RData ...")

pm_den_region <- bin1d_counts(
  x = dt_global$pa_pm25,
  by_group = dt_global$region,
  width = PM_BIN_WIDTH,
  x_min = 0,
  x_max = PM_PLOT_MAX
)

temp_den_region <- bin1d_counts(
  x = dt_global$temp_resid,
  by_group = dt_global$region,
  width = TEMP_BIN_WIDTH,
  x_min = TEMP_RESID_MIN,
  x_max = TEMP_RESID_MAX
)

rh_den_region <- bin1d_counts(
  x = dt_global$rh_resid,
  by_group = dt_global$region,
  width = RH_BIN_WIDTH,
  x_min = RH_RESID_MIN,
  x_max = RH_RESID_MAX
)

save_density_plot(
  pm_den_region,
  "PM2.5 distribution by region",
  "PurpleAir PM2.5 (µg/m³)",
  file.path(OUTDIR, "density_pm_by_region.png")
)

save_density_plot(
  temp_den_region,
  "Temperature residual distribution by region",
  "PurpleAir temperature - ERA5 temperature (°C)",
  file.path(OUTDIR, "density_temp_resid_by_region.png")
)

save_density_plot(
  rh_den_region,
  "RH residual distribution by region",
  "PurpleAir RH - ERA5 RH (%)",
  file.path(OUTDIR, "density_rh_resid_by_region.png")
)

rm(pm_den_region, temp_den_region, rh_den_region)
gc()

pm2d_region <- bin2d_count(
  x = dt_global$geos_pm,
  y = dt_global$pa_pm25,
  by_group = dt_global$region,
  x_width = PM_X_BIN,
  y_width = PM_Y_BIN,
  x_min = 0, x_max = PM_PLOT_MAX,
  y_min = 0, y_max = PM_PLOT_MAX
)

temp2d_region <- bin2d_count(
  x = dt_global$era5_temp,
  y = dt_global$pa_temp,
  by_group = dt_global$region,
  x_width = TEMP_X_BIN,
  y_width = TEMP_Y_BIN
)

rh2d_region <- bin2d_count(
  x = dt_global$era5_rh,
  y = dt_global$pa_rh,
  by_group = dt_global$region,
  x_width = RH_X_BIN,
  y_width = RH_Y_BIN,
  x_min = 0, x_max = 100,
  y_min = 0, y_max = 100
)

save_binned_heatmap(
  pm2d_region,
  "PurpleAir PM2.5 versus GEOS-CF PM2.5 by region",
  "GEOS-CF PM2.5 (µg/m³)",
  "PurpleAir PM2.5 (µg/m³)",
  file.path(OUTDIR, "binned_pm_vs_geos_by_region.png"),
  palette = "magma",
  trans = "log10",
  fixed_scales = TRUE,
  add_identity = TRUE
)

save_binned_heatmap(
  temp2d_region,
  "PurpleAir temperature versus ERA5 temperature by region",
  "ERA5 temperature (°C)",
  "PurpleAir temperature (°C)",
  file.path(OUTDIR, "binned_temp_vs_era5_by_region.png"),
  palette = "viridis",
  trans = "log10",
  fixed_scales = FALSE,
  add_identity = TRUE
)

save_binned_heatmap(
  rh2d_region,
  "PurpleAir RH versus ERA5 RH by region",
  "ERA5 RH (%)",
  "PurpleAir RH (%)",
  file.path(OUTDIR, "binned_rh_vs_era5_by_region.png"),
  palette = "cividis",
  trans = "log10",
  fixed_scales = TRUE,
  add_identity = TRUE
)

rm(pm2d_region, temp2d_region, rh2d_region)
gc()

pm_resid_region <- bin2d_count(
  x = dt_global$geos_pm,
  y = dt_global$pm_resid,
  by_group = dt_global$region,
  x_width = PM_X_BIN,
  y_width = PM_Y_BIN,
  x_min = 0, x_max = PM_PLOT_MAX,
  y_min = PM_RESID_MIN, y_max = PM_RESID_MAX
)

temp_resid_region <- bin2d_count(
  x = dt_global$era5_temp,
  y = dt_global$temp_resid,
  by_group = dt_global$region,
  x_width = TEMP_X_BIN,
  y_width = TEMP_Y_BIN,
  y_min = TEMP_RESID_MIN,
  y_max = TEMP_RESID_MAX
)

rh_resid_region <- bin2d_count(
  x = dt_global$era5_rh,
  y = dt_global$rh_resid,
  by_group = dt_global$region,
  x_width = RH_X_BIN,
  y_width = RH_Y_BIN,
  x_min = 0, x_max = 100,
  y_min = RH_RESID_MIN,
  y_max = RH_RESID_MAX
)

save_binned_heatmap(
  pm_resid_region,
  "PM residual versus GEOS-CF PM2.5 by region",
  "GEOS-CF PM2.5 (µg/m³)",
  "PurpleAir PM2.5 - GEOS-CF PM2.5 (µg/m³)",
  file.path(OUTDIR, "residual_pm_vs_pm_by_region.png"),
  palette = "inferno",
  trans = "log10",
  fixed_scales = TRUE,
  add_zero = TRUE
)

save_binned_heatmap(
  temp_resid_region,
  "Temperature residual versus ERA5 temperature by region",
  "ERA5 temperature (°C)",
  "PurpleAir temperature - ERA5 temperature (°C)",
  file.path(OUTDIR, "residual_temp_vs_temp_by_region.png"),
  palette = "plasma",
  trans = "log10",
  fixed_scales = FALSE,
  add_zero = TRUE
)

save_binned_heatmap(
  rh_resid_region,
  "RH residual versus ERA5 RH by region",
  "ERA5 RH (%)",
  "PurpleAir RH - ERA5 RH (%)",
  file.path(OUTDIR, "residual_rh_vs_rh_by_region.png"),
  palette = "cividis",
  trans = "log10",
  fixed_scales = TRUE,
  add_zero = TRUE
)

rm(pm_resid_region, temp_resid_region, rh_resid_region)
gc()

rm(dt_global)
gc()

# ------------------------------------------------------------
# 9) CLIMATE DATASET
# ------------------------------------------------------------

say("Loading climate-zone data from outdoor_with_koppen_timevarying.RData ...")

outdoor_climate <- load_outdoor_object(CLIMATE_INFILE)
setDT(outdoor_climate)
outdoor_climate[, station := as.character(station)]

say("Climate rows loaded before restriction: ", format(nrow(outdoor_climate), big.mark = ","))
say("Climate stations loaded before restriction: ", format(uniqueN(outdoor_climate$station), big.mark = ","))

if (USE_SINGLE_ZONE_ONLY_FOR_CLIMATE) {

  if (!file.exists(KOPPEN_CHANGE_FILE)) {
    stop("Required Köppen transition file not found: ", KOPPEN_CHANGE_FILE)
  }

  station_change <- fread(KOPPEN_CHANGE_FILE)
  station_change[, station := as.character(station)]

  keep_stations <- station_change[
    !is.na(first_koppen) & n_koppen_codes == 1,
    station
  ]

  outdoor_climate <- outdoor_climate[station %in% keep_stations]

  say("Restricted climate data to valid single-zone stations only.")
  say("Climate rows after single-zone restriction: ", format(nrow(outdoor_climate), big.mark = ","))
  say("Climate stations after single-zone restriction: ", format(uniqueN(outdoor_climate$station), big.mark = ","))

  rm(station_change, keep_stations)
  gc()
}

outdoor_climate <- outdoor_climate[
  !is.na(get(koppen_col)) &
    get(koppen_col) != ""
]

say("Climate rows after valid Köppen-code restriction: ", format(nrow(outdoor_climate), big.mark = ","))
say("Climate stations after valid Köppen-code restriction: ", format(uniqueN(outdoor_climate$station), big.mark = ","))

dt_climate <- make_analysis_dt(outdoor_climate, require_koppen = TRUE)

rm(outdoor_climate)
gc()

qa_climate <- make_qa(dt_climate, "outdoor_with_koppen_single_zone_climate")
fwrite(qa_climate, file.path(OUTDIR, "qa_climate_single_zone.csv"))

# ------------------------------------------------------------
# 10) CLIMATE SUMMARY + CORRELATION
# ------------------------------------------------------------

say("Writing climate-zone summary and correlation tables ...")

summary_by_climate <- make_summary(dt_climate, "koppen_code")[order(-n_rows)]

cor_by_climate <- make_cor(dt_climate, "koppen_code")[order(-n_rows)]

fwrite(summary_by_climate, file.path(OUTDIR, "summary_stats_by_climate_zone.csv"))
fwrite(cor_by_climate, file.path(OUTDIR, "correlations_by_climate_zone.csv"))

save_correlation_figures(
  cor_dt = cor_by_climate,
  group_col = "koppen_code",
  group_label = "Köppen climate zone",
  file_suffix = "by_climate",
  top_n = 12
)

# ------------------------------------------------------------
# 11) CLIMATE-ZONE PLOTS
# ------------------------------------------------------------

say("Building climate-zone plots from single-zone climate dataset ...")

top_climates <- summary_by_climate[order(-n_rows)][1:min(12, .N), koppen_code]

dt_climate_plot <- dt_climate[koppen_code %in% top_climates]

say("Climate zones plotted: ", paste(sort(unique(dt_climate_plot$koppen_code)), collapse = ", "))

pm_resid_climate <- bin2d_count(
  x = dt_climate_plot$geos_pm,
  y = dt_climate_plot$pm_resid,
  by_group = dt_climate_plot$koppen_code,
  x_width = PM_X_BIN,
  y_width = PM_Y_BIN,
  x_min = 0, x_max = PM_PLOT_MAX,
  y_min = PM_RESID_MIN, y_max = PM_RESID_MAX
)

save_binned_heatmap(
  pm_resid_climate,
  "PM residual versus GEOS-CF PM2.5 by climate zone",
  "GEOS-CF PM2.5 (µg/m³)",
  "PurpleAir PM2.5 - GEOS-CF PM2.5 (µg/m³)",
  file.path(OUTDIR, "pm_residual_surface_by_climate.png"),
  palette = "inferno",
  trans = "log10",
  fixed_scales = FALSE,
  add_zero = TRUE
)

rm(pm_resid_climate)
gc()

# Optional climate-zone agreement plots
pm2d_climate <- bin2d_count(
  x = dt_climate_plot$geos_pm,
  y = dt_climate_plot$pa_pm25,
  by_group = dt_climate_plot$koppen_code,
  x_width = PM_X_BIN,
  y_width = PM_Y_BIN,
  x_min = 0, x_max = PM_PLOT_MAX,
  y_min = 0, y_max = PM_PLOT_MAX
)

save_binned_heatmap(
  pm2d_climate,
  "PurpleAir PM2.5 versus GEOS-CF PM2.5 by climate zone",
  "GEOS-CF PM2.5 (µg/m³)",
  "PurpleAir PM2.5 (µg/m³)",
  file.path(OUTDIR, "binned_pm_vs_geos_by_climate.png"),
  palette = "magma",
  trans = "log10",
  fixed_scales = FALSE,
  add_identity = TRUE
)

rm(pm2d_climate, dt_climate_plot)
gc()

# ------------------------------------------------------------
# 12) SAVE SMALL ANALYSIS SNAPSHOTS
# ------------------------------------------------------------

say("Saving small analysis snapshots ...")

set.seed(123)

global_snapshot <- summary_by_region
climate_snapshot <- summary_by_climate

fwrite(global_snapshot, file.path(OUTDIR, "analysis_snapshot_region_summary.csv"))
fwrite(climate_snapshot, file.path(OUTDIR, "analysis_snapshot_climate_summary.csv"))

rm(dt_climate)
gc()

say("Done.")
say("Wrote outputs to: ", OUTDIR)

# ---- source chunk #27 (Rmd lines 8542-9219) ------------------------

# ============================================================
# MEMORY-EFFICIENT STAND-ALONE REPLACEMENT
#
# GLOBAL/REGION PLOTS FROM:
#   outdoor_final.RData
#
# CLIMATE-ZONE PLOTS FROM:
#   outdoor_with_koppen_timevarying.RData
#
# IMPORTANT:
#   - PM values are NOT capped.
#   - Values outside the PM plotting window are simply not displayed.
#   - This prevents artificial piles/clusters at PM_PLOT_MAX.
#   - Climate-zone plots are restricted to valid single-zone stations.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# ------------------------------------------------------------
# 0) PATHS / SETTINGS
# ------------------------------------------------------------

GLOBAL_REGION_INFILE <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_final.RData"
CLIMATE_INFILE <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_with_koppen_timevarying.RData"
KOPPEN_CHANGE_FILE <- "/Users/priyanka/Downloads/PurpleAir_EPA/station_koppen_transition_summary.csv"

OUTDIR <- "/Users/priyanka/Downloads/PurpleAir_EPA/agreement_residual_figures_no_pm_lumping"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

USE_SINGLE_ZONE_ONLY_FOR_CLIMATE <- TRUE

# Display window only; PM values are NOT capped/lumped.
PM_PLOT_MAX <- 1500

PM_RESID_MIN <- -1500
PM_RESID_MAX <- 1500

TEMP_RESID_MIN <- -40
TEMP_RESID_MAX <- 40

RH_RESID_MIN <- -100
RH_RESID_MAX <- 100

PM_X_BIN <- 5
PM_Y_BIN <- 5
PM_RESID_BIN <- 5

TEMP_X_BIN <- 0.5
TEMP_Y_BIN <- 0.5
TEMP_RESID_BIN <- 0.25

RH_X_BIN <- 1
RH_Y_BIN <- 1
RH_RESID_BIN <- 1

TOP_N_CLIMATE_ZONES <- Inf

MAKE_LOG_PM_PLOTS <- TRUE

say <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste0(..., collapse = ""), "\n")
}

cap <- function(x, lo = -Inf, hi = Inf) {
  pmax(pmin(x, hi), lo)
}

# ------------------------------------------------------------
# 1) COLUMN MAPPING
# ------------------------------------------------------------

a_col <- "pm25_cf_1_a_hourly(ug/m3)_pm25_cf_1_a"
b_col <- "pm25_cf_1_b_hourly(ug/m3)_pm25_cf_1_b"
pa_pm_avg_col <- "pm25_cf1_avg"

geos_col <- "geoscf_pm25"

pa_temp_col <- "temperature_hourly(C)_temperature"
pa_rh_col   <- "humidity_hourly(%)_humidity"

era5_temp_col <- "era5_t2m_c"
era5_rh_col   <- "era5_rh_pct"

koppen_col <- "koppen_code"

# ------------------------------------------------------------
# 2) LOAD HELPER
# ------------------------------------------------------------

load_outdoor_object <- function(infile) {
  e <- new.env(parent = globalenv())
  load(infile, envir = e)

  if (exists("outdoor", envir = e)) {
    x <- get("outdoor", envir = e)
  } else if (exists("outdoor_final", envir = e)) {
    x <- get("outdoor_final", envir = e)
  } else {
    stop("Could not find object named 'outdoor' or 'outdoor_final' in: ", infile)
  }

  rm(e)
  setDT(x)
  x[]
}

check_required_cols <- function(dt, require_koppen = FALSE) {
  required <- c(
    "station", "region",
    geos_col,
    pa_temp_col, pa_rh_col,
    era5_temp_col, era5_rh_col
  )

  if (!(pa_pm_avg_col %in% names(dt))) {
    required <- c(required, a_col, b_col)
  }

  if (require_koppen) {
    required <- c(required, koppen_col)
  }

  missing <- setdiff(required, names(dt))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }

  invisible(TRUE)
}

# ------------------------------------------------------------
# 3) BUILD MINIMAL DATASET
# ------------------------------------------------------------

make_minimal_dt <- function(dt, require_koppen = FALSE) {
  check_required_cols(dt, require_koppen = require_koppen)

  keep_cols <- c(
    "station", "region",
    geos_col,
    pa_temp_col, pa_rh_col,
    era5_temp_col, era5_rh_col
  )

  if (pa_pm_avg_col %in% names(dt)) {
    keep_cols <- c(keep_cols, pa_pm_avg_col)
  } else {
    keep_cols <- c(keep_cols, a_col, b_col)
  }

  if (koppen_col %in% names(dt)) {
    keep_cols <- c(keep_cols, koppen_col)
  }

  keep_cols <- unique(keep_cols)
  dt <- dt[, ..keep_cols]

  dt[, station := as.character(station)]
  dt[, region := as.character(region)]

  if (koppen_col %in% names(dt)) {
    dt[, koppen_code := as.character(get(koppen_col))]
    if (koppen_col != "koppen_code") {
      dt[, (koppen_col) := NULL]
    }
  } else {
    dt[, koppen_code := NA_character_]
  }

  if (pa_pm_avg_col %in% names(dt)) {
    setnames(dt, pa_pm_avg_col, "pa_pm25")
    dt[, pa_pm25 := as.numeric(pa_pm25)]
  } else {
    dt[, pa_pm25 := (as.numeric(get(a_col)) + as.numeric(get(b_col))) / 2]
    dt[, c(a_col, b_col) := NULL]
  }

  setnames(dt, geos_col, "geos_pm")
  setnames(dt, pa_temp_col, "pa_temp")
  setnames(dt, pa_rh_col, "pa_rh")
  setnames(dt, era5_temp_col, "era5_temp")
  setnames(dt, era5_rh_col, "era5_rh")

  dt[, `:=`(
    geos_pm = as.numeric(geos_pm),
    pa_temp = as.numeric(pa_temp),
    pa_rh = as.numeric(pa_rh),
    era5_temp = as.numeric(era5_temp),
    era5_rh = as.numeric(era5_rh)
  )]

  dt <- dt[!is.na(region)]

  gc()
  dt[]
}

# ------------------------------------------------------------
# 4) MEMORY-SAFE BINNING
# ------------------------------------------------------------

bin_agreement <- function(dt, x_expr, y_expr, group_var,
                          x_width, y_width,
                          x_min, x_max, y_min, y_max) {

  x_call <- substitute(x_expr)
  y_call <- substitute(y_expr)

  tmp <- data.table(
    group = dt[[group_var]],
    x = eval(x_call, envir = dt, enclos = parent.frame()),
    y = eval(y_call, envir = dt, enclos = parent.frame())
  )

  tmp <- tmp[
    !is.na(group) &
      is.finite(x) & is.finite(y) &
      x >= x_min & x <= x_max &
      y >= y_min & y <= y_max
  ]

  tmp[, `:=`(
    xmid = x_min + floor((x - x_min) / x_width) * x_width + x_width / 2,
    ymid = y_min + floor((y - y_min) / y_width) * y_width + y_width / 2
  )]

  out <- tmp[, .N, by = .(group, xmid, ymid)]
  setnames(out, "N", "count")

  rm(tmp)
  gc()

  out[]
}

bin_residual <- function(dt, x_expr, resid_expr, group_var,
                         x_width, y_width,
                         x_min, x_max, y_min, y_max) {

  x_call <- substitute(x_expr)
  y_call <- substitute(resid_expr)

  tmp <- data.table(
    group = dt[[group_var]],
    x = eval(x_call, envir = dt, enclos = parent.frame()),
    y = eval(y_call, envir = dt, enclos = parent.frame())
  )

  tmp <- tmp[
    !is.na(group) &
      is.finite(x) & is.finite(y) &
      x >= x_min & x <= x_max &
      y >= y_min & y <= y_max
  ]

  tmp[, `:=`(
    xmid = x_min + floor((x - x_min) / x_width) * x_width + x_width / 2,
    ymid = y_min + floor((y - y_min) / y_width) * y_width + y_width / 2
  )]

  h <- tmp[, .N, by = .(group, xmid, ymid)]
  setnames(h, "N", "count")

  s <- tmp[, .(
    mean_resid = mean(y, na.rm = TRUE),
    p25_resid = quantile(y, 0.25, na.rm = TRUE),
    p75_resid = quantile(y, 0.75, na.rm = TRUE),
    n = .N
  ), by = .(group, xmid)]

  rm(tmp)
  gc()

  list(h = h[], s = s[])
}

# ------------------------------------------------------------
# 5) PLOTTING HELPERS
# ------------------------------------------------------------

save_agreement_plot <- function(h, title, xlab, ylab, outfile,
                                palette = "viridis", fixed = TRUE,
                                one_to_one = TRUE) {

  p <- ggplot(h, aes(xmid, ymid, fill = count)) +
    geom_tile() +
    scale_fill_viridis_c(option = palette, trans = "log10", name = "Count") +
    facet_wrap(~group, scales = ifelse(fixed, "fixed", "free")) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey92", color = "grey60"),
      strip.text = element_text(face = "bold"),
      panel.border = element_rect(color = "grey70", fill = NA),
      plot.title = element_text(face = "bold")
    ) +
    labs(title = title, x = xlab, y = ylab)

  if (one_to_one) {
    p <- p + geom_abline(intercept = 0, slope = 1, linetype = 2,
                         linewidth = 0.4, color = "white")
  }

  ggsave(file.path(OUTDIR, outfile), p, width = 12, height = 7, dpi = 300, bg = "white")

  rm(p)
  gc()
}

save_residual_plot <- function(h, s, title, xlab, ylab, outfile,
                               palette = "viridis", fixed = TRUE) {

  p <- ggplot() +
    geom_tile(data = h, aes(xmid, ymid, fill = count)) +
    scale_fill_viridis_c(option = palette, trans = "log10", name = "Count") +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4, color = "white") +
    geom_ribbon(
      data = s,
      aes(x = xmid, ymin = p25_resid, ymax = p75_resid),
      inherit.aes = FALSE,
      fill = "cyan",
      alpha = 0.2
    ) +
    geom_line(
      data = s,
      aes(x = xmid, y = mean_resid),
      inherit.aes = FALSE,
      color = "cyan",
      linewidth = 0.8
    ) +
    facet_wrap(~group, scales = ifelse(fixed, "fixed", "free")) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey92", color = "grey60"),
      strip.text = element_text(face = "bold"),
      panel.border = element_rect(color = "grey70", fill = NA),
      plot.title = element_text(face = "bold")
    ) +
    labs(title = title, x = xlab, y = ylab)

  ggsave(file.path(OUTDIR, outfile), p, width = 12, height = 7, dpi = 300, bg = "white")

  rm(p)
  gc()
}

# ------------------------------------------------------------
# 6) SAVE ONE SET OF FIGURES
# ------------------------------------------------------------

save_all_plot_set <- function(dt, group_var, label_suffix, file_suffix) {

  say("Saving ", label_suffix, " figures...")

  if (group_var == "overall") {
    dt[, overall := "All observations"]
    group_var_use <- "overall"
  } else {
    group_var_use <- group_var
  }

  temp_x_min <- floor(min(dt$era5_temp, na.rm = TRUE))
  temp_x_max <- ceiling(max(dt$era5_temp, na.rm = TRUE))
  temp_y_min <- floor(min(dt$pa_temp, na.rm = TRUE))
  temp_y_max <- ceiling(max(dt$pa_temp, na.rm = TRUE))

  # ---------------- PM agreement: linear window, NO capping ----------------
  h <- bin_agreement(
    dt,
    x_expr = geos_pm,
    y_expr = pa_pm25,
    group_var = group_var_use,
    x_width = PM_X_BIN,
    y_width = PM_Y_BIN,
    x_min = 0,
    x_max = PM_PLOT_MAX,
    y_min = 0,
    y_max = PM_PLOT_MAX
  )

  save_agreement_plot(
    h,
    paste0("PurpleAir PM2.5 versus GEOS-CF PM2.5 ", label_suffix),
    "GEOS-CF PM2.5 (µg/m³)",
    "PurpleAir PM2.5 (µg/m³)",
    paste0("agreement_pm_", file_suffix, ".png"),
    "magma",
    TRUE,
    TRUE
  )

  rm(h)
  gc()

  # ---------------- Optional PM agreement: log10 scale, NO capping ----------------
  if (MAKE_LOG_PM_PLOTS) {

    h <- bin_agreement(
      dt,
      x_expr = log10(geos_pm + 1),
      y_expr = log10(pa_pm25 + 1),
      group_var = group_var_use,
      x_width = 0.025,
      y_width = 0.025,
      x_min = 0,
      x_max = log10(PM_PLOT_MAX + 1),
      y_min = 0,
      y_max = log10(PM_PLOT_MAX + 1)
    )

    save_agreement_plot(
      h,
      paste0("PurpleAir PM2.5 versus GEOS-CF PM2.5 ", label_suffix, " — log10 scale"),
      expression(log[10]("GEOS-CF PM2.5 + 1")),
      expression(log[10]("PurpleAir PM2.5 + 1")),
      paste0("agreement_pm_log10_", file_suffix, ".png"),
      "magma",
      TRUE,
      TRUE
    )

    rm(h)
    gc()
  }

  # ---------------- Temperature agreement ----------------
  h <- bin_agreement(
    dt,
    x_expr = era5_temp,
    y_expr = pa_temp,
    group_var = group_var_use,
    x_width = TEMP_X_BIN,
    y_width = TEMP_Y_BIN,
    x_min = temp_x_min,
    x_max = temp_x_max,
    y_min = temp_y_min,
    y_max = temp_y_max
  )

  save_agreement_plot(
    h,
    paste0("PurpleAir temperature versus ERA5 temperature ", label_suffix),
    "ERA5 temperature (°C)",
    "PurpleAir temperature (°C)",
    paste0("agreement_temp_", file_suffix, ".png"),
    "viridis",
    FALSE,
    TRUE
  )

  rm(h)
  gc()

  # ---------------- RH agreement ----------------
  h <- bin_agreement(
    dt,
    x_expr = cap(era5_rh, 0, 100),
    y_expr = cap(pa_rh, 0, 100),
    group_var = group_var_use,
    x_width = RH_X_BIN,
    y_width = RH_Y_BIN,
    x_min = 0,
    x_max = 100,
    y_min = 0,
    y_max = 100
  )

  save_agreement_plot(
    h,
    paste0("PurpleAir RH versus ERA5 RH ", label_suffix),
    "ERA5 RH (%)",
    "PurpleAir RH (%)",
    paste0("agreement_rh_", file_suffix, ".png"),
    "cividis",
    TRUE,
    TRUE
  )

  rm(h)
  gc()

  # ---------------- PM residual: NO capping ----------------
  res <- bin_residual(
    dt,
    x_expr = geos_pm,
    resid_expr = pa_pm25 - geos_pm,
    group_var = group_var_use,
    x_width = PM_X_BIN,
    y_width = PM_RESID_BIN,
    x_min = 0,
    x_max = PM_PLOT_MAX,
    y_min = PM_RESID_MIN,
    y_max = PM_RESID_MAX
  )

  save_residual_plot(
    res$h, res$s,
    paste0("PM2.5 residual versus GEOS-CF PM2.5 ", label_suffix),
    "GEOS-CF PM2.5 (µg/m³)",
    "PurpleAir PM2.5 − GEOS-CF PM2.5 (µg/m³)",
    paste0("residual_pm_", file_suffix, ".png"),
    "inferno",
    TRUE
  )

  rm(res)
  gc()

  # ---------------- Temperature residual ----------------
  res <- bin_residual(
    dt,
    x_expr = era5_temp,
    resid_expr = pa_temp - era5_temp,
    group_var = group_var_use,
    x_width = TEMP_X_BIN,
    y_width = TEMP_RESID_BIN,
    x_min = temp_x_min,
    x_max = temp_x_max,
    y_min = TEMP_RESID_MIN,
    y_max = TEMP_RESID_MAX
  )

  save_residual_plot(
    res$h, res$s,
    paste0("Temperature residual versus ERA5 temperature ", label_suffix),
    "ERA5 temperature (°C)",
    "PurpleAir temperature − ERA5 temperature (°C)",
    paste0("residual_temp_", file_suffix, ".png"),
    "plasma",
    FALSE
  )

  rm(res)
  gc()

  # ---------------- RH residual ----------------
  res <- bin_residual(
    dt,
    x_expr = cap(era5_rh, 0, 100),
    resid_expr = pa_rh - era5_rh,
    group_var = group_var_use,
    x_width = RH_X_BIN,
    y_width = RH_RESID_BIN,
    x_min = 0,
    x_max = 100,
    y_min = RH_RESID_MIN,
    y_max = RH_RESID_MAX
  )

  save_residual_plot(
    res$h, res$s,
    paste0("RH residual versus ERA5 RH ", label_suffix),
    "ERA5 RH (%)",
    "PurpleAir RH − ERA5 RH (%)",
    paste0("residual_rh_", file_suffix, ".png"),
    "cividis",
    TRUE
  )

  rm(res)
  gc()

  if (group_var == "overall") {
    dt[, overall := NULL]
    gc()
  }

  invisible(TRUE)
}

# ------------------------------------------------------------
# 7) GLOBAL + REGION
# ------------------------------------------------------------

say("Loading global/region input: ", GLOBAL_REGION_INFILE)
outdoor_gr <- load_outdoor_object(GLOBAL_REGION_INFILE)

say("Global/region rows loaded: ", format(nrow(outdoor_gr), big.mark = ","))
say("Global/region stations loaded: ", format(uniqueN(as.character(outdoor_gr$station)), big.mark = ","))

dt_gr <- make_minimal_dt(outdoor_gr, require_koppen = FALSE)
rm(outdoor_gr)
gc()

save_all_plot_set(
  dt = dt_gr,
  group_var = "overall",
  label_suffix = "(overall)",
  file_suffix = "overall"
)

save_all_plot_set(
  dt = dt_gr,
  group_var = "region",
  label_suffix = "by region",
  file_suffix = "by_region"
)

rm(dt_gr)
gc()

# ------------------------------------------------------------
# 8) CLIMATE ZONE
# ------------------------------------------------------------

say("Loading climate-zone input: ", CLIMATE_INFILE)
outdoor_clim <- load_outdoor_object(CLIMATE_INFILE)
setDT(outdoor_clim)
outdoor_clim[, station := as.character(station)]

say("Climate input rows loaded before restriction: ", format(nrow(outdoor_clim), big.mark = ","))
say("Climate input stations before restriction: ", format(uniqueN(outdoor_clim$station), big.mark = ","))

if (USE_SINGLE_ZONE_ONLY_FOR_CLIMATE) {

  if (!file.exists(KOPPEN_CHANGE_FILE)) {
    stop("Required station Köppen transition file not found: ", KOPPEN_CHANGE_FILE)
  }

  station_change <- fread(KOPPEN_CHANGE_FILE)
  station_change[, station := as.character(station)]

  keep_stations <- station_change[
    !is.na(first_koppen) & n_koppen_codes == 1,
    station
  ]

  outdoor_clim <- outdoor_clim[station %in% keep_stations]

  say("Restricted climate data to valid single-zone stations only.")
  say("Climate rows after single-zone restriction: ", format(nrow(outdoor_clim), big.mark = ","))
  say("Climate stations after single-zone restriction: ", format(uniqueN(outdoor_clim$station), big.mark = ","))

  rm(station_change, keep_stations)
  gc()
}

outdoor_clim <- outdoor_clim[
  !is.na(get(koppen_col)) &
    get(koppen_col) != ""
]

say("Climate rows after valid Köppen-code restriction: ", format(nrow(outdoor_clim), big.mark = ","))
say("Climate stations after valid Köppen-code restriction: ", format(uniqueN(outdoor_clim$station), big.mark = ","))

dt_clim <- make_minimal_dt(outdoor_clim, require_koppen = TRUE)
rm(outdoor_clim)
gc()

top_climates <- dt_clim[
  !is.na(koppen_code) & koppen_code != "",
  .N,
  by = koppen_code
][order(-N)][1:min(TOP_N_CLIMATE_ZONES, .N), koppen_code]

dt_clim <- dt_clim[koppen_code %in% top_climates]

say("Climate zones included: ", paste(sort(unique(dt_clim$koppen_code)), collapse = ", "))

save_all_plot_set(
  dt = dt_clim,
  group_var = "koppen_code",
  label_suffix = "by Köppen climate zone",
  file_suffix = "by_climate_zone"
)

rm(dt_clim)
gc()

say("Done. Figures saved to: ", OUTDIR)

# ---- source chunk #28 (Rmd lines 9222-9987) ------------------------

# ============================================================
# SENSOR DISTRIBUTIONS, METRIC 3 THRESHOLDS, AND DEPLOYMENT TIMELINES
#
# Non-climate figures:
#   outdoor_final.RData
#
# Climate figures:
#   outdoor_with_koppen_timevarying.RData
#   restricted to single-zone Köppen stations
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# ------------------------------------------------------------
# 0) PATHS
# ------------------------------------------------------------

global_rdata  <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_final.RData"
climate_rdata <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_with_koppen_timevarying.RData"

CHANGE_FILE <- "/Users/priyanka/Downloads/PurpleAir_EPA/station_koppen_transition_summary.csv"

OUTDIR <- "/Users/priyanka/Downloads/PurpleAir_EPA/sensor_distribution_timeline_figures"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Metric 3 threshold file, if already created by earlier code.
# If this file does not exist, the code will try to compute Metric 3 thresholds from paired A/B channels.
metric3_threshold_file <- file.path(OUTDIR, "metric3_thresholds_by_station.csv")

# ------------------------------------------------------------
# 1) SETTINGS
# ------------------------------------------------------------

a_col <- "pm25_cf_1_a_hourly(ug/m3)_pm25_cf_1_a"
b_col <- "pm25_cf_1_b_hourly(ug/m3)_pm25_cf_1_b"

datetime_col <- "datetime_utc"
station_col  <- "station"
region_col   <- "region"
koppen_col   <- "koppen_code"

metric3_prob <- 0.85

region_order <- c("Africa", "CONUS", "India", "SouthAmerica")

koppen_major_colors <- c(
  "Af"  = "#1b9e77", "Am"  = "#2ca25f", "Aw"  = "#78c679",
  "BSh" = "#e6550d", "BSk" = "#fdae6b", "BWh" = "#8c2d04", "BWk" = "#d94801",
  "Cfa" = "#08519c", "Cfb" = "#4292c6", "Csa" = "#084594", "Csb" = "#3182bd",
  "Cwa" = "#2171b5", "Cwb" = "#6baed6", "Cwc" = "#9ecae1",
  "Dfa" = "#88419d", "Dfb" = "#8c6bb1", "Dfc" = "#bcbddc",
  "Dsa" = "#54278f", "Dsb" = "#756bb1", "Dsc" = "#9e9ac8",
  "Dwa" = "#6a51a3", "Dwb" = "#807dba", "Dwd" = "#dadaeb",
  "ET"  = "#969696"
)

say <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste0(..., collapse = ""), "\n")
}

safe_file <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

# ------------------------------------------------------------
# 2) LOAD HELPERS
# ------------------------------------------------------------

load_outdoor_object <- function(infile, keep_cols = NULL) {
  e <- new.env(parent = emptyenv())
  load(infile, envir = e)

  obj_name <- intersect(c("outdoor", "outdoor_final"), ls(e))[1]
  if (is.na(obj_name)) stop("Could not find object named outdoor or outdoor_final in: ", infile)

  x <- get(obj_name, envir = e)
  setDT(x)

  if (!is.null(keep_cols)) {
    keep_cols <- intersect(keep_cols, names(x))
    x <- x[, ..keep_cols]
  }

  x[]
}

restrict_single_zone <- function(dt) {
  station_change <- fread(CHANGE_FILE)
  station_change[, station := as.character(station)]

  keep_stations <- station_change[
    !is.na(first_koppen) & n_koppen_codes == 1,
    .(
      station,
      koppen_code = as.character(first_koppen),
      koppen_desc = as.character(first_koppen_desc)
    )
  ]

  dt[, station := as.character(station)]

  if (koppen_col %in% names(dt)) dt[, (koppen_col) := NULL]

  dt <- merge(dt, keep_stations, by = "station", all = FALSE)

  dt[is.na(koppen_desc) | koppen_desc == "", koppen_desc := koppen_code]
  dt[, climate_label := paste0(koppen_code, "\n", koppen_desc)]

  rm(station_change, keep_stations)
  gc()

  dt[]
}

# ------------------------------------------------------------
# 3) STATION-LEVEL SUMMARIES
# ------------------------------------------------------------

make_station_summary <- function(infile, climate = FALSE) {

  keep_cols <- c(station_col, region_col, datetime_col)
  if (climate) keep_cols <- c(keep_cols, koppen_col)

  dt <- load_outdoor_object(infile, keep_cols = keep_cols)

  setnames(dt, old = c(station_col, region_col, datetime_col),
           new = c("station", "region", "datetime_utc"),
           skip_absent = TRUE)

  dt[, station := as.character(station)]
  dt[, region := as.character(region)]
  dt[, datetime_utc := as.POSIXct(datetime_utc, tz = "UTC")]

  dt <- dt[!is.na(station) & !is.na(region) & !is.na(datetime_utc)]

  if (climate) {
    dt <- restrict_single_zone(dt)
  }

  by_vars <- if (climate) c("station", "region", "koppen_code", "koppen_desc", "climate_label") else c("station", "region")

  out <- dt[
    ,
    .(
      start_date = as.Date(min(datetime_utc, na.rm = TRUE)),
      end_date   = as.Date(max(datetime_utc, na.rm = TRUE)),
      n_obs      = .N
    ),
    by = by_vars
  ]

  out[, months_active := as.numeric(difftime(end_date, start_date, units = "days")) / 30.4375]
  out[, years_active := months_active / 12]

  rm(dt)
  gc()

  out[]
}

say("Creating non-climate station summary")
station_global <- make_station_summary(global_rdata, climate = FALSE)

say("Creating climate station summary")
station_climate <- make_station_summary(climate_rdata, climate = TRUE)

fwrite(station_global, file.path(OUTDIR, "station_summary_global_region.csv"))
fwrite(station_climate, file.path(OUTDIR, "station_summary_single_zone_climate.csv"))

# ------------------------------------------------------------
# 4) FIGURE: GLOBAL + REGIONAL CLIMATE-ZONE DISTRIBUTION
# ------------------------------------------------------------

climate_counts_global <- station_climate[
  ,
  .N,
  by = koppen_code
][
  order(-N)
]

climate_counts_global[, percent := 100 * N / sum(N)]
climate_order <- climate_counts_global$koppen_code

climate_counts_region <- station_climate[
  ,
  .N,
  by = .(region, koppen_code)
]

region_totals <- station_climate[, .(region_total = .N), by = region]

climate_counts_region <- merge(climate_counts_region, region_totals, by = "region", all.x = TRUE)
climate_counts_region[, percent := 100 * N / region_total]

climate_counts_global[, koppen_code := factor(koppen_code, levels = climate_order)]
climate_counts_region[, koppen_code := factor(koppen_code, levels = climate_order)]
climate_counts_region[, region := factor(region, levels = region_order)]

p_clim_global <- ggplot(climate_counts_global, aes(koppen_code, percent, fill = koppen_code)) +
  geom_col(width = 0.85) +
  scale_fill_manual(values = koppen_major_colors, drop = FALSE) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "A) Global distribution of sensors by climate zone",
    x = "Köppen-Geiger climate zone",
    y = "Percent of sensors"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0, face = "plain")
  )

p_clim_region <- ggplot(climate_counts_region, aes(koppen_code, percent, fill = koppen_code)) +
  geom_col(width = 0.85) +
  facet_wrap(~ region, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = koppen_major_colors, drop = FALSE) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "B) Regional distribution of sensors by climate zone",
    x = "Köppen-Geiger climate zone",
    y = "Percent of sensors",
    fill = "Köppen class"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0, face = "plain")
  )

fig_climate_distribution <- p_clim_global / p_clim_region +
  plot_layout(heights = c(0.9, 1.55), guides = "collect")

ggsave(
  file.path(OUTDIR, "fig_climate_zone_distribution_global_region.png"),
  fig_climate_distribution,
  width = 14,
  height = 10,
  dpi = 600
)

ggsave(
  file.path(OUTDIR, "fig_climate_zone_distribution_global_region.jpeg"),
  fig_climate_distribution,
  width = 14,
  height = 10,
  dpi = 600
)

# ------------------------------------------------------------
# 5) METRIC 3 THRESHOLD BY REGION
# ------------------------------------------------------------
# If metric3_threshold_file already exists, this section uses it.
# Otherwise, it computes station-level 85th percentile percent difference:
#   abs(A - B) / rowMeans(A, B)
# capped visually at 2 in the plot only.

make_metric3_thresholds <- function() {

  if (file.exists(metric3_threshold_file)) {
    say("Reading existing Metric 3 threshold file")
    m3 <- fread(metric3_threshold_file)
    return(m3[])
  }

  say("Computing Metric 3 thresholds from outdoor_final.RData")

  needed <- c(station_col, region_col, a_col, b_col)

  dt <- load_outdoor_object(global_rdata, keep_cols = needed)

  setnames(dt, old = c(station_col, region_col), new = c("station", "region"), skip_absent = TRUE)

  dt[, station := as.character(station)]
  dt[, region := as.character(region)]

  dt[, a := as.numeric(get(a_col))]
  dt[, b := as.numeric(get(b_col))]

  dt <- dt[!is.na(station) & !is.na(region) & is.finite(a) & is.finite(b)]

  dt[, denom := rowMeans(cbind(a, b), na.rm = TRUE)]
  dt <- dt[is.finite(denom) & denom > 0]

  dt[, percent_difference := abs(a - b) / denom]
  dt <- dt[is.finite(percent_difference)]

  m3 <- dt[
    ,
    .(
      metric3_threshold_85 = as.numeric(quantile(percent_difference, probs = metric3_prob, na.rm = TRUE, names = FALSE)),
      n_obs = .N
    ),
    by = .(station, region)
  ]

  fwrite(m3, metric3_threshold_file)

  rm(dt)
  gc()

  m3[]
}

metric3 <- make_metric3_thresholds()
metric3 <- metric3[is.finite(metric3_threshold_85)]
metric3[, region := factor(region, levels = region_order)]

p_metric3_region <- ggplot(metric3, aes(x = metric3_threshold_85, y = region)) +
  geom_boxplot(outlier.shape = NA, width = 0.75, fill = "white", color = "grey20", linewidth = 0.8) +
  geom_jitter(height = 0.05, width = 0, alpha = 0.18, size = 1.1, color = "grey35") +
  coord_cartesian(xlim = c(0, 2)) +
  scale_x_continuous(breaks = seq(0, 2, 0.5)) +
  labs(
    title = "Distribution of Metric 3 threshold (85th percentile) across regions",
    x = "85th percentile threshold (percent difference)",
    y = "Region"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(OUTDIR, "fig_metric3_threshold_85_by_region.png"),
  p_metric3_region,
  width = 10,
  height = 7,
  dpi = 600
)

ggsave(
  file.path(OUTDIR, "fig_metric3_threshold_85_by_region.jpeg"),
  p_metric3_region,
  width = 10,
  height = 7,
  dpi = 600
)

# ------------------------------------------------------------
# 6) DEPLOYMENT TIMELINE + LONGEVITY HELPERS
# ------------------------------------------------------------

make_month_seq <- function(start_date, end_date) {
  start_month <- as.Date(format(as.Date(start_date), "%Y-%m-01"))
  end_month   <- as.Date(format(as.Date(end_date), "%Y-%m-01"))

  if (is.na(start_month) || is.na(end_month) || end_month < start_month) {
    return(as.Date(character()))
  }

  seq.Date(start_month, end_month, by = "month")
}

make_pretty_gantt_data <- function(dt, group_var = NULL, max_sensors_per_group = 500) {
  dd <- copy(dt)

  if (is.null(group_var)) {
    dd[, group := "Global"]
  } else {
    dd[, group := as.character(get(group_var))]
  }

  dd <- dd[!is.na(group) & !is.na(start_date) & !is.na(end_date)]
  dd[, duration_days := as.numeric(end_date - start_date)]
  setorder(dd, group, start_date, end_date)

  dd <- dd[
    ,
    {
      if (.N <= max_sensors_per_group) {
        .SD
      } else {
        idx <- unique(round(seq(1, .N, length.out = max_sensors_per_group)))
        .SD[idx]
      }
    },
    by = group
  ]

  dd[, sensor_plot_id := seq_len(.N), by = group]
  dd[]
}

plot_pretty_gantt <- function(dt, group_var = NULL, title, outfile,
                              max_sensors_per_group = 500,
                              width = 13, height = 8,
                              ncol = 1) {
  dd <- make_pretty_gantt_data(
    dt = dt,
    group_var = group_var,
    max_sensors_per_group = max_sensors_per_group
  )

  p <- ggplot(dd) +
    geom_segment(
      aes(
        x = start_date,
        xend = end_date,
        y = sensor_plot_id,
        yend = sensor_plot_id,
        color = duration_days
      ),
      linewidth = 0.42,
      alpha = 0.72,
      lineend = "round"
    ) +
    facet_wrap(~ group, scales = "free_y", ncol = ncol) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    scale_color_viridis_c(
      option = "plasma",
      labels = comma,
      name = "Deployment\nlength (days)"
    ) +
    labs(
      title = title,
      x = "Calendar date",
      y = "Sensors sorted by deployment start"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 16),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    )

  ggsave(outfile, p, width = width, height = height, dpi = 600, limitsize = FALSE)
  invisible(p)
}

plot_global_region_gantt <- function(dt, outfile) {
  global_dt <- copy(dt)
  global_dt[, panel_group := "Global"]

  region_dt <- copy(dt)
  region_dt[, panel_group := as.character(region)]

  both <- rbindlist(
    list(
      global_dt[, .(station, region, start_date, end_date, n_obs, months_active, years_active, panel_group)],
      region_dt[, .(station, region, start_date, end_date, n_obs, months_active, years_active, panel_group)]
    ),
    use.names = TRUE,
    fill = TRUE
  )

  panel_order <- c("Global", region_order)
  both[, panel_group := factor(panel_group, levels = panel_order)]

  plot_pretty_gantt(
    both,
    group_var = "panel_group",
    title = "PurpleAir sensor deployment timelines: Global and region",
    outfile = outfile,
    max_sensors_per_group = 900,
    width = 15,
    height = 12,
    ncol = 1
  )
}

plot_longevity_histograms <- function(dt, outfile) {
  global_dt <- copy(dt)
  global_dt[, panel_group := "Global"]

  region_dt <- copy(dt)
  region_dt[, panel_group := as.character(region)]

  both <- rbindlist(
    list(
      global_dt[, .(station, region, years_active, months_active, panel_group)],
      region_dt[, .(station, region, years_active, months_active, panel_group)]
    ),
    use.names = TRUE,
    fill = TRUE
  )

  panel_order <- c("Global", region_order)
  both[, panel_group := factor(panel_group, levels = panel_order)]

  p <- ggplot(both[is.finite(years_active)], aes(x = years_active)) +
    geom_histogram(
      bins = 45,
      fill = "grey35",
      color = "white",
      linewidth = 0.15
    ) +
    facet_wrap(~ panel_group, scales = "free_y", ncol = 2) +
    scale_x_continuous(labels = number_format(accuracy = 0.5)) +
    scale_y_continuous(labels = comma) +
    labs(
      title = "Distribution of sensor longevity: Global, CONUS, and regions",
      subtitle = "Longevity is end date minus start date for each sensor.",
      x = "Sensor longevity (years)",
      y = "Number of sensors"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 17),
      plot.subtitle = element_text(size = 11),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  ggsave(outfile, p, width = 13, height = 9, dpi = 600, limitsize = FALSE)
  invisible(p)
}

# ------------------------------------------------------------
# 7) ACTIVE SENSOR HEATMAP + COHORT HEATMAP HELPERS
#     Log-scaled color axes reveal low-count variation.
# ------------------------------------------------------------

make_sensor_month_heatmap_data <- function(dt, group_var = NULL) {
  dd <- copy(dt)

  if (is.null(group_var)) {
    dd[, group := "Global"]
  } else {
    dd[, group := as.character(get(group_var))]
  }

  dd <- dd[!is.na(group) & !is.na(start_date) & !is.na(end_date)]

  active_months <- dd[
    ,
    {
      mseq <- make_month_seq(start_date[1], end_date[1])
      .(month = mseq)
    },
    by = .(station, group)
  ]

  hm <- active_months[
    ,
    .(n_active_sensors = uniqueN(station)),
    by = .(group, month)
  ]

  hm[]
}

plot_active_sensor_heatmap <- function(dt, group_var = NULL, title, outfile,
                                       width = 13, height = 7) {
  hm <- make_sensor_month_heatmap_data(dt, group_var = group_var)

  p <- ggplot(hm, aes(x = month, y = group, fill = n_active_sensors)) +
    geom_tile(color = "white", linewidth = 0.08) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0, 0)) +
    scale_fill_viridis_c(
      option = "magma",
      trans = "log10",
      breaks = c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000),
      labels = comma,
      name = "Active\nsensors\n(log10)"
    ) +
    labs(
      title = title,
      x = "Calendar month",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 15),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      legend.position = "right"
    )

  ggsave(outfile, p, width = width, height = height, dpi = 600, limitsize = FALSE)
  invisible(p)
}

make_cohort_heatmap_data <- function(dt, group_var = NULL) {
  dd <- copy(dt)

  if (is.null(group_var)) {
    dd[, group := "Global"]
  } else {
    dd[, group := as.character(get(group_var))]
  }

  dd <- dd[!is.na(group) & !is.na(start_date) & !is.na(end_date)]

  dd[, start_month := as.Date(format(as.Date(start_date), "%Y-%m-01"))]
  dd[, end_month   := as.Date(format(as.Date(end_date), "%Y-%m-01"))]

  dd <- dd[!is.na(start_month) & !is.na(end_month) & end_month >= start_month]

  cohort_months <- dd[
    ,
    {
      active_seq <- seq.Date(start_month[1], end_month[1], by = "month")
      .(active_month = active_seq)
    },
    by = .(station, group, start_month)
  ]

  cohort <- cohort_months[
    ,
    .(n_active_sensors = uniqueN(station)),
    by = .(group, start_month, active_month)
  ]

  cohort[]
}

plot_cohort_heatmap <- function(dt, group_var = NULL, title, outfile,
                                width = 12, height = 9,
                                ncol = 1) {
  cohort <- make_cohort_heatmap_data(dt, group_var = group_var)

  p <- ggplot(cohort, aes(x = active_month, y = start_month, fill = n_active_sensors)) +
    geom_tile(color = "grey90", linewidth = 0.04) +
    facet_wrap(~ group, scales = "free", ncol = ncol) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0, 0)) +
    scale_y_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0, 0)) +
    scale_fill_viridis_c(
      option = "inferno",
      trans = "log10",
      breaks = c(1, 3, 10, 30, 100, 300, 1000, 3000),
      labels = comma,
      name = "Sensors\nactive\n(log10)"
    ) +
    labs(
      title = title,
      x = "Calendar month active",
      y = "Sensor start cohort"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 15),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    )

  ggsave(outfile, p, width = width, height = height, dpi = 600, limitsize = FALSE)
  invisible(p)
}

plot_global_region_cohort_heatmap <- function(dt, outfile) {
  global_dt <- copy(dt)
  global_dt[, panel_group := "Global"]

  region_dt <- copy(dt)
  region_dt[, panel_group := as.character(region)]

  both <- rbindlist(
    list(
      global_dt[, .(station, region, start_date, end_date, panel_group)],
      region_dt[, .(station, region, start_date, end_date, panel_group)]
    ),
    use.names = TRUE,
    fill = TRUE
  )

  panel_order <- c("Global", region_order)
  both[, panel_group := factor(panel_group, levels = panel_order)]

  plot_cohort_heatmap(
    both,
    group_var = "panel_group",
    title = "Sensor deployment cohorts and persistence: Global and region",
    outfile = outfile,
    width = 15,
    height = 13,
    ncol = 1
  )
}

# ------------------------------------------------------------
# 8) WRITE TIMELINE, LONGEVITY, AND HEATMAP FIGURES
# ------------------------------------------------------------

say("Writing combined global + regional Gantt figure")

plot_global_region_gantt(
  station_global,
  outfile = file.path(OUTDIR, "fig_pretty_gantt_global_and_region.png")
)

say("Writing sensor longevity histograms")

plot_longevity_histograms(
  station_global,
  outfile = file.path(OUTDIR, "fig_sensor_longevity_histograms_global_CONUS_region.png")
)

say("Writing climate-zone Gantt figure")

plot_pretty_gantt(
  station_climate,
  group_var = "koppen_code",
  title = "PurpleAir sensor deployment timelines by Köppen-Geiger climate zone",
  outfile = file.path(OUTDIR, "fig_pretty_gantt_by_climate_zone.png"),
  max_sensors_per_group = 300,
  width = 15,
  height = 13,
  ncol = 1
)

say("Writing active-sensor heatmaps with log color scales")

plot_active_sensor_heatmap(
  station_global,
  group_var = NULL,
  title = "Active PurpleAir sensors over time: Global",
  outfile = file.path(OUTDIR, "fig_active_sensor_heatmap_global_logscale.png"),
  width = 13,
  height = 4
)

plot_active_sensor_heatmap(
  station_global,
  group_var = "region",
  title = "Active PurpleAir sensors over time by region",
  outfile = file.path(OUTDIR, "fig_active_sensor_heatmap_by_region_logscale.png"),
  width = 13,
  height = 6
)

plot_active_sensor_heatmap(
  station_climate,
  group_var = "koppen_code",
  title = "Active PurpleAir sensors over time by Köppen-Geiger climate zone",
  outfile = file.path(OUTDIR, "fig_active_sensor_heatmap_by_climate_zone_logscale.png"),
  width = 15,
  height = 9
)

say("Writing cohort heatmaps with log color scales")

plot_global_region_cohort_heatmap(
  station_global,
  outfile = file.path(OUTDIR, "fig_cohort_heatmap_global_and_region_logscale.png")
)

plot_cohort_heatmap(
  station_climate,
  group_var = "koppen_code",
  title = "Sensor deployment cohorts and persistence by Köppen-Geiger climate zone",
  outfile = file.path(OUTDIR, "fig_cohort_heatmap_by_climate_zone_logscale.png"),
  width = 15,
  height = 14,
  ncol = 1
)
