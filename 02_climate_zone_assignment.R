# ==========================================================================
# 08_koppen_climate.R
#
# Time-varying Koppen-Geiger climate-zone join, climate-zone change analysis, and the
# climate-zone figure (Figure S2, Figure S3).
#
# Extracted VERBATIM from purpleair_degradation2.Rmd (chunks #22 lines 6337-6673, #23 lines 6676-7000, #24 lines 7003-7312).
# Code is unmodified; see REPRODUCING.md for paths and run order.
# ==========================================================================

# ---- source chunk #22 (Rmd lines 6337-6673) ------------------------

# ============================================================
# TIME-VARYING Köppen–Geiger join + station climate-zone change analysis
#
# What this does:
# 1) loads outdoor_final.RData
# 2) chooses time-varying coordinates (PurpleAir first, then GEE fallback)
# 3) collapses to one location per station-month
# 4) extracts Köppen class for each station-month
#    - simple nearest-cell extraction
#    - bilinear fallback
#    - buffered modal fallback
# 5) identifies stations that move between climate zones over time
# 6) joins station-month Köppen class back to outdoor via station + ym
# 7) writes outputs
#
# Outputs:
#   - station_month_koppen.csv
#   - station_koppen_transition_summary.csv
#   - outdoor_with_koppen_timevarying.RData
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(terra)
  library(lubridate)
})

# ------------------------------------------------------------
# 0) PATHS
# ------------------------------------------------------------
tif_path    <- "/Users/priyanka/Downloads/koppen_geiger_tif/1991_2020/koppen_geiger_0p00833333.tif"
legend_path <- "/Users/priyanka/Downloads/koppen_geiger_tif/legend.txt"

out_rdata <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_with_koppen_timevarying.RData"
out_station_month_csv <- "/Users/priyanka/Downloads/PurpleAir_EPA/station_month_koppen.csv"
out_transition_csv    <- "/Users/priyanka/Downloads/PurpleAir_EPA/station_koppen_transition_summary.csv"

# ------------------------------------------------------------
# 1) HELPERS
# ------------------------------------------------------------
rename_if_present <- function(dt, old, new) {
  if (old %in% names(dt)) setnames(dt, old, new)
}

mode_value <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_character_)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

say <- function(...) cat(sprintf(...), "\n")

# ------------------------------------------------------------
# 2) LOAD outdoor_final
# ------------------------------------------------------------
load("/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_final.RData")
outdoor <- outdoor_final
rm(outdoor_final)
gc()

setDT(outdoor)
outdoor[, station := as.character(station)]

if (!inherits(outdoor$datetime_utc, "POSIXct")) {
  outdoor[, datetime_utc := as.POSIXct(datetime_utc, tz = "UTC")]
} else {
  attr(outdoor$datetime_utc, "tzone") <- "UTC"
}

say("Outdoor rows loaded: %s", format(nrow(outdoor), big.mark = ","))
#452,143,148
say("Unique stations: %s", format(uniqueN(outdoor$station), big.mark = ","))
#36,328 

# ------------------------------------------------------------
# 3) CHOOSE TIME-VARYING COORDINATES
#    Priority:
#      1) PurpleAir lon/lat if present
#      2) fallback to GEE lon/lat
# ------------------------------------------------------------
rename_if_present(outdoor, "LONGITUDE(deg)", "pa_longitude")
rename_if_present(outdoor, "LATITUDE(deg)",  "pa_latitude")

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

# keep only what we need
keep_cols <- intersect(
  c("station", "datetime_utc", "region", "longitude", "latitude"),
  names(outdoor)
)
outdoor_geo <- outdoor[
  !is.na(station) &
    !is.na(datetime_utc) &
    is.finite(longitude) &
    is.finite(latitude),
  ..keep_cols
]

say("Rows with valid coords: %s", format(nrow(outdoor_geo), big.mark = ","))
#452,143,148
say("Stations with valid coords: %s", format(uniqueN(outdoor_geo$station), big.mark = ","))
#36,328

# ------------------------------------------------------------
# 4) BUILD STATION-MONTH LOCATION TABLE
#    One representative location per station-month
# ------------------------------------------------------------
outdoor_geo[, ym := floor_date(datetime_utc, "month")]

station_month <- outdoor_geo[
  ,
  .(
    region = mode_value(region),
    longitude = median(longitude, na.rm = TRUE),
    latitude  = median(latitude,  na.rm = TRUE),
    n_rows_month = .N,
    first_seen_month = min(datetime_utc),
    last_seen_month  = max(datetime_utc)
  ),
  by = .(station, ym)
]

setorder(station_month, station, ym)

say("Station-month rows: %s", format(nrow(station_month), big.mark = ","))
#929,758
say("Unique stations in station-month table: %s", format(uniqueN(station_month$station), big.mark = ","))
#36,328

# ------------------------------------------------------------
# 5) READ KÖPPEN RASTER + LEGEND
# ------------------------------------------------------------
r <- rast(tif_path)

leg_lines <- readLines(legend_path, warn = FALSE)
leg_lines <- leg_lines[grepl("^\\s*\\d+\\s*:\\s*", leg_lines)]

m <- regexec("^\\s*(\\d+)\\s*:\\s*([A-Za-z]{2,3})\\s+([^\\[]+)", leg_lines)
parts <- regmatches(leg_lines, m)

legend_dt <- rbindlist(lapply(parts, function(x) {
  if (length(x) < 4) return(NULL)
  data.table(
    kg_id = as.integer(x[2]),
    koppen_code = x[3],
    koppen_desc = trimws(x[4])
  )
}), use.names = TRUE, fill = TRUE)

stopifnot(nrow(legend_dt) > 0)
setkey(legend_dt, kg_id)

# ------------------------------------------------------------
# 6) GUARDRAILS AGAINST EDGE NAs
# ------------------------------------------------------------
exr <- ext(r)
eps <- 1e-6

station_month[, `:=`(
  longitude = pmin(pmax(longitude, exr[1] + eps), exr[2] - eps),
  latitude  = pmin(pmax(latitude,  exr[3] + eps), exr[4] - eps)
)]

pts <- vect(
  station_month[, .(longitude, latitude)],
  geom = c("longitude", "latitude"),
  crs = "EPSG:4326"
)

if (!same.crs(pts, r)) pts <- project(pts, crs(r))

# ------------------------------------------------------------
# 7) EXTRACT KÖPPEN CLASS FOR EACH STATION-MONTH
#    simple -> bilinear -> buffer fallback
# ------------------------------------------------------------
ex1 <- terra::extract(r, pts, method = "simple")
val_col <- setdiff(names(ex1), "ID")[1]
kg_id <- as.integer(ex1[[val_col]])

na_idx <- which(is.na(kg_id))
say("NA after simple extract: %s of %s", format(length(na_idx), big.mark = ","), format(length(kg_id), big.mark = ","))
#NA after simple extract: 5,490 of 929,758 

if (length(na_idx) > 0) {
  ex2 <- terra::extract(r, pts[na_idx], method = "bilinear")
  val2 <- setdiff(names(ex2), "ID")[1]
  kg2 <- suppressWarnings(as.integer(round(ex2[[val2]])))
  kg_id[na_idx] <- kg2
}

na_idx2 <- which(is.na(kg_id))
say("NA after bilinear fallback: %s of %s", format(length(na_idx2), big.mark = ","), format(length(kg_id), big.mark = ","))
#NA after bilinear fallback: 312 of 929,758

if (length(na_idx2) > 0) {
  pts_m <- project(pts, "EPSG:3857")
  r_m   <- project(r, "EPSG:3857", method = "near")

  buf_m <- 2000
  ex3 <- terra::extract(
    r_m,
    buffer(pts_m[na_idx2], width = buf_m),
    fun = "modal",
    na.rm = TRUE
  )
  val3 <- setdiff(names(ex3), "ID")[1]
  kg3 <- as.integer(ex3[[val3]])
  kg_id[na_idx2] <- kg3
}

say("Final NA after buffer fallback: %s of %s",
    format(sum(is.na(kg_id)), big.mark = ","),
    format(length(kg_id), big.mark = ","))
#Final NA after buffer fallback: 141 of 929,758 

station_month[, kg_id := kg_id]

station_month[legend_dt, on = .(kg_id), `:=`(
  koppen_code = i.koppen_code,
  koppen_desc = i.koppen_desc
)]

# ------------------------------------------------------------
# 8) CHECK WHETHER STATIONS CHANGE CLIMATE ZONE OVER TIME
# ------------------------------------------------------------
station_koppen_change <- station_month[
  ,
  .(
    region = mode_value(region),
    n_months = .N,
    n_koppen_ids = uniqueN(kg_id, na.rm = TRUE),
    n_koppen_codes = uniqueN(koppen_code, na.rm = TRUE),
    first_month = min(ym, na.rm = TRUE),
    last_month  = max(ym, na.rm = TRUE),
    first_koppen = koppen_code[which.min(ym)],
    last_koppen  = koppen_code[which.max(ym)],
    first_koppen_desc = koppen_desc[which.min(ym)],
    last_koppen_desc  = koppen_desc[which.max(ym)]
  ),
  by = station
]

station_koppen_change[, changed_koppen_zone :=
  n_koppen_codes > 1
]

# add a string of all observed classes
koppen_path <- station_month[
  !is.na(koppen_code),
  .(koppen_sequence = paste(unique(koppen_code), collapse = " -> ")),
  by = station
]

station_koppen_change <- koppen_path[station_koppen_change, on = "station"]

say("Stations with >=1 Köppen class: %s",
    format(sum(!is.na(station_koppen_change$first_koppen)), big.mark = ","))
#Stations with >=1 Köppen class: 36,318 

say("Stations changing Köppen zone over time: %s",
    format(sum(station_koppen_change$changed_koppen_zone, na.rm = TRUE), big.mark = ","))
#Stations changing Köppen zone over time: 568 

say("Percent changing Köppen zone: %.2f%%",
    100 * mean(station_koppen_change$changed_koppen_zone, na.rm = TRUE))
#Percent changing Köppen zone: 1.56%

# optional regional summary
koppen_change_region <- station_koppen_change[
  ,
  .(
    stations = .N,
    changed_koppen_zone = sum(changed_koppen_zone, na.rm = TRUE),
    pct_changed = round(100 * mean(changed_koppen_zone, na.rm = TRUE), 2)
  ),
  by = region
][order(-stations)]

say("Regional summary of climate-zone changes:")
print(koppen_change_region)
#         region stations changed_koppen_zone pct_changed
#        CONUS    33915                 497        1.47
# SouthAmerica     1007                  38        3.77
#        India      779                  18        2.31
#       Africa      627                  15        2.39

# ------------------------------------------------------------
# 9) JOIN STATION-MONTH KÖPPEN BACK TO outdoor
#    join by station + ym, not by station only
# ------------------------------------------------------------
outdoor[, ym := floor_date(datetime_utc, "month")]

station_month_keep <- station_month[, .(
  station, ym, kg_id, koppen_code, koppen_desc
)]

setkey(station_month_keep, station, ym)
setkey(outdoor, station, ym)

outdoor <- station_month_keep[outdoor]

say("Share of outdoor rows with Köppen code: %.2f%%",
    100 * mean(!is.na(outdoor$koppen_code)))
#Share of outdoor rows with Köppen code: 99.99% 

# ------------------------------------------------------------
# 10) SAVE OUTPUTS
# ------------------------------------------------------------
fwrite(station_month, out_station_month_csv)
fwrite(station_koppen_change[order(-changed_koppen_zone, station)], out_transition_csv)
save(outdoor, file = out_rdata)

say("Wrote station-month Köppen table to:")
say("%s", out_station_month_csv)

say("Wrote station Köppen transition summary to:")
say("%s", out_transition_csv)

say("Wrote outdoor table with time-varying Köppen join to:")
say("%s", out_rdata)

# ---- source chunk #23 (Rmd lines 6676-7000) ------------------------

library(data.table)
library(terra)
library(lubridate)
say <- function(...) cat(sprintf(...), "\n")

station_koppen_change <- fread(
  "/Users/priyanka/Downloads/PurpleAir_EPA/station_koppen_transition_summary.csv"
)

# ensure types
station_koppen_change[, station := as.character(station)]
# ------------------------------------------------------------
# 11) CLIMATE-ZONE SUMMARY FOR STATIONS THAT REMAIN IN ONE ZONE
#     Excludes stations that changed Köppen zone over time
# ------------------------------------------------------------

out_climate_summary_csv <- "/Users/priyanka/Downloads/PurpleAir_EPA/station_climate_zone_summary_single_zone_only.csv"
out_climate_summary_wide_csv <- "/Users/priyanka/Downloads/PurpleAir_EPA/station_climate_zone_summary_single_zone_only_wide.csv"

# identify stations that stayed in a single Köppen zone
station_koppen_single <- station_koppen_change[
  !is.na(first_koppen) & n_koppen_codes == 1,
  .(
    station,
    region,
    koppen_code = first_koppen,
    koppen_desc = first_koppen_desc,
    n_months,
    first_month,
    last_month
  )
]

# stations excluded because they changed climate zone
station_koppen_changed <- station_koppen_change[
  changed_koppen_zone == TRUE,
  .(
    station,
    region,
    n_months,
    n_koppen_codes,
    koppen_sequence,
    first_koppen,
    last_koppen
  )
]

say("Stations retained for climate-zone summaries (single zone only): %s",
    format(nrow(station_koppen_single), big.mark = ","))
say("Stations excluded for changing Köppen zone: %s",
    format(nrow(station_koppen_changed), big.mark = ","))

# ------------------------------------------------------------
# 11A) GLOBAL SUMMARY BY CLIMATE ZONE
# ------------------------------------------------------------
climate_summary_global <- station_koppen_single[
  ,
  .(n_sensors = .N),
  by = .(koppen_code, koppen_desc)
][order(-n_sensors, koppen_code)]

climate_summary_global[, pct_within_global_single_zone :=
  round(100 * n_sensors / sum(n_sensors), 2)
]
climate_summary_global[, region := "Global"]

# ------------------------------------------------------------
# 11B) REGIONAL SUMMARY BY CLIMATE ZONE
# ------------------------------------------------------------
climate_summary_region <- station_koppen_single[
  ,
  .(n_sensors = .N),
  by = .(region, koppen_code, koppen_desc)
][order(region, -n_sensors, koppen_code)]

climate_summary_region[, pct_within_region_single_zone :=
  round(100 * n_sensors / sum(n_sensors), 2),
  by = region
]

# ------------------------------------------------------------
# 11C) EXCLUSION ROWS
#     Add explicit row for stations excluded due to zone change
# ------------------------------------------------------------
excluded_row_global <- data.table(
  region = "Global",
  koppen_code = "EXCLUDED_CHANGED_ZONE",
  koppen_desc = "Stations changing Köppen zone during study period; excluded from climate-zone summaries",
  n_sensors = nrow(station_koppen_changed),
  pct_within_global_single_zone = NA_real_
)

excluded_row_region <- station_koppen_changed[
  ,
  .(n_sensors = .N),
  by = region
][
  ,
  `:=`(
    koppen_code = "EXCLUDED_CHANGED_ZONE",
    koppen_desc = "Stations changing Köppen zone during study period; excluded from climate-zone summaries",
    pct_within_region_single_zone = NA_real_
  )
]

# ------------------------------------------------------------
# 11D) COMBINED LONG TABLE
# ------------------------------------------------------------
climate_summary_all <- rbindlist(
  list(
    climate_summary_global[, .(
      region, koppen_code, koppen_desc, n_sensors, pct = pct_within_global_single_zone
    )],
    excluded_row_global[, .(
      region, koppen_code, koppen_desc, n_sensors, pct = pct_within_global_single_zone
    )],
    climate_summary_region[, .(
      region, koppen_code, koppen_desc, n_sensors, pct = pct_within_region_single_zone
    )],
    excluded_row_region[, .(
      region, koppen_code, koppen_desc, n_sensors, pct = pct_within_region_single_zone
    )]
  ),
  use.names = TRUE,
  fill = TRUE
)

# order regions explicitly
region_levels <- c("Global", sort(unique(station_koppen_single$region)))

climate_summary_all[, region_order := match(region, region_levels)]

setorder(climate_summary_all, region_order, -n_sensors, koppen_code)

climate_summary_all[, region_order := NULL]

# ------------------------------------------------------------
# 11E) WIDE TABLE FOR WORD / MANUSCRIPT
# ------------------------------------------------------------
# counts only
climate_summary_wide <- dcast(
  climate_summary_all,
  koppen_code + koppen_desc ~ region,
  value.var = "n_sensors",
  fill = 0
)

# make excluded row appear at bottom
climate_summary_wide[
  koppen_code == "EXCLUDED_CHANGED_ZONE",
  sort_order := 9999
]
climate_summary_wide[
  koppen_code != "EXCLUDED_CHANGED_ZONE",
  sort_order := 1
]
setorder(climate_summary_wide, sort_order, koppen_code)
climate_summary_wide[, sort_order := NULL]

cat("\nWIDE CLIMATE-ZONE TABLE:\n")
print(climate_summary_wide)

# ------------------------------------------------------------
# 11F) OPTIONAL REGION x CLIMATE TABLE INCLUDING RETAINED TOTALS
# ------------------------------------------------------------
single_zone_totals <- station_koppen_single[, .(n_single_zone_sensors = .N), by = region]
changed_zone_totals <- station_koppen_changed[, .(n_changed_zone_sensors = .N), by = region]

summary_region_totals <- merge(
  single_zone_totals,
  changed_zone_totals,
  by = "region",
  all = TRUE
)

summary_region_totals[is.na(n_single_zone_sensors), n_single_zone_sensors := 0L]
summary_region_totals[is.na(n_changed_zone_sensors), n_changed_zone_sensors := 0L]
summary_region_totals[, n_total_with_koppen := n_single_zone_sensors + n_changed_zone_sensors]
summary_region_totals[, pct_changed_zone :=
  round(100 * n_changed_zone_sensors / n_total_with_koppen, 2)
]

global_region_totals <- data.table(
  region = "Global",
  n_single_zone_sensors = nrow(station_koppen_single),
  n_changed_zone_sensors = nrow(station_koppen_changed),
  n_total_with_koppen = nrow(station_koppen_single) + nrow(station_koppen_changed),
  pct_changed_zone = round(
    100 * nrow(station_koppen_changed) /
      (nrow(station_koppen_single) + nrow(station_koppen_changed)),
    2
  )
)

summary_region_totals <- rbindlist(
  list(global_region_totals, summary_region_totals),
  use.names = TRUE,
  fill = TRUE
)

cat("\nREGIONAL TOTALS FOR CLIMATE-ZONE ANALYSIS:\n")
print(summary_region_totals)

# ------------------------------------------------------------
# 11G) SAVE OUTPUTS
# ------------------------------------------------------------
fwrite(climate_summary_all, out_climate_summary_csv)
fwrite(climate_summary_wide, out_climate_summary_wide_csv)

say("Wrote single-zone climate summary (long) to:")
say("%s", out_climate_summary_csv)

say("Wrote single-zone climate summary (wide) to:")
say("%s", out_climate_summary_wide_csv)

# ------------------------------------------------------------
# 12) FIGURE: CLIMATE ZONE DISTRIBUTION (MATCHED COLORS)
# ------------------------------------------------------------
library(ggplot2)
library(scales)
library(patchwork)

# ------------------------------------------------------------
# PREP DATA
# ------------------------------------------------------------
plot_dt <- copy(climate_summary_all)

plot_dt <- plot_dt[!is.na(koppen_code) & koppen_code != ""]
plot_dt[, is_excluded := koppen_code == "EXCLUDED_CHANGED_ZONE"]

# keep only valid zones
plot_dt_pct <- plot_dt[is_excluded == FALSE & !is.na(pct)]

# order by global importance
global_order <- plot_dt_pct[
  region == "Global",
  .(n_sensors = sum(n_sensors)),
  by = koppen_code
][order(-n_sensors)]$koppen_code

plot_dt_pct[, koppen_code := factor(koppen_code, levels = global_order)]

# ------------------------------------------------------------
# COLOR PALETTE (USE SAME AS MAP)
# ------------------------------------------------------------
# assumes koppen_colors already exists earlier in your script

# ensure palette matches only present levels
plot_colors <- koppen_colors[names(koppen_colors) %in% levels(plot_dt_pct$koppen_code)]

# safety fallback
missing_cols <- setdiff(levels(plot_dt_pct$koppen_code), names(plot_colors))
if (length(missing_cols) > 0) {
  plot_colors <- c(plot_colors, setNames(rep("#999999", length(missing_cols)), missing_cols))
}

# ------------------------------------------------------------
# 12A) GLOBAL FIGURE
# ------------------------------------------------------------
p_global <- ggplot(
  plot_dt_pct[region == "Global"],
  aes(x = koppen_code, y = pct, fill = koppen_code)
) +
  geom_col(width = 0.85) +
  scale_fill_manual(values = plot_colors, guide = "none") +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "A) Global distribution of sensors by climate zone",
    x = "Köppen–Geiger climate zone",
    y = "Percent of sensors"
  )

# ------------------------------------------------------------
# 12B) REGIONAL FIGURE
# ------------------------------------------------------------
p_region <- ggplot(
  plot_dt_pct[region != "Global"],
  aes(x = koppen_code, y = pct, fill = koppen_code)
) +
  geom_col(width = 0.85) +
  scale_fill_manual(values = plot_colors, name = "Köppen class") +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  facet_wrap(~region, scales = "free_y") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  ) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(
    title = "B) Regional distribution of sensors by climate zone",
    x = "Köppen–Geiger climate zone",
    y = "Percent of sensors"
  )

# ------------------------------------------------------------
# 12C) COMBINE
# ------------------------------------------------------------
p_climate_combined <- p_global / p_region +
  plot_layout(heights = c(1, 1.6))

# ------------------------------------------------------------
# 12D) SAVE
# ------------------------------------------------------------
out_fig_climate <- "/Users/priyanka/Downloads/PurpleAir_EPA/fig_climate_zone_distribution.png"

ggsave(
  out_fig_climate,
  p_climate_combined,
  width = 13,
  height = 10,
  dpi = 300,
  bg = "white"
)

cat("Wrote climate-zone distribution figure:\n", out_fig_climate, "\n")

# ---- source chunk #24 (Rmd lines 7003-7312) ------------------------

# ============================================================
# PAPER FIGURE: KÖPPEN–GEIGER CLIMATE ZONES FOR STUDY REGIONS
#
# Produces:
# 1) global overview map with study-region boxes
# 2) 4-panel regional climate-zone map:
#    - CONUS
#    - India
#    - South America
#    - Africa
#
# Inputs:
#   koppen_geiger_0p00833333.tif
#   legend.txt
#
# Output:
#   /Users/priyanka/Downloads/PurpleAir_EPA/fig_climate_zones_study_regions.png
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(terra)
  library(sf)
  library(ggplot2)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(patchwork)
})


koppen_colors <- c(
  Af  = "#006d2c", Am  = "#31a354", Aw  = "#74c476",
  BWh = "#8c2d04", BWk = "#d94801", BSh = "#f16913", BSk = "#fdae6b",
  Csa = "#08519c", Csb = "#3182bd", Csc = "#6baed6",
  Cwa = "#2171b5", Cwb = "#4292c6", Cwc = "#9ecae1",
  Cfa = "#084594", Cfb = "#4a90d9", Cfc = "#bdd7e7",
  Dsa = "#54278f", Dsb = "#756bb1", Dsc = "#9e9ac8", Dsd = "#cbc9e2",
  Dwa = "#6a51a3", Dwb = "#807dba", Dwc = "#9e9ac8", Dwd = "#dadaeb",
  Dfa = "#88419d", Dfb = "#8c6bb1", Dfc = "#bcbddc", Dfd = "#efedf5",
  ET  = "#969696", EF  = "#d9d9d9"
)


# ------------------------------------------------------------
# 0) PATHS
# ------------------------------------------------------------
tif_path <- "/Users/priyanka/Downloads/koppen_geiger_tif/1991_2020/koppen_geiger_0p00833333.tif"
legend_path <- "/Users/priyanka/Downloads/koppen_geiger_tif/legend.txt"
out_png <- "/Users/priyanka/Downloads/PurpleAir_EPA/fig_climate_zones_study_regions.png"

# ------------------------------------------------------------
# 1) HELPERS
# ------------------------------------------------------------
say <- function(...) cat(sprintf(...), "\n")

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

make_region_plot <- function(title_text, xlim, ylim, show_legend = FALSE) {
  dt <- koppen_df[
    x >= xlim[1] & x <= xlim[2] &
      y >= ylim[1] & y <= ylim[2]
  ]

  bb <- st_as_sfc(
    st_bbox(
      c(xmin = xlim[1], ymin = ylim[1], xmax = xlim[2], ymax = ylim[2]),
      crs = st_crs(4326)
    )
  )

  world_crop <- suppressWarnings(st_intersection(world, bb))

  ggplot() +
    geom_tile(
      data = dt,
      aes(x = x, y = y, fill = factor(koppen_code, levels = koppen_levels)),
      width = tile_w,
      height = tile_h
    ) +
    geom_sf(data = world_crop, fill = NA, color = "grey25", linewidth = 0.2) +
    scale_fill_manual(
      values = koppen_colors,
      drop = FALSE,
      name = "Köppen class"
    ) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      legend.position = if (show_legend) "right" else "none",
      plot.title = element_text(face = "bold", size = 11)
    ) +
    labs(title = title_text)
}

# ------------------------------------------------------------
# 2) READ KÖPPEN RASTER + LEGEND
# ------------------------------------------------------------
r <- rast(tif_path)

leg_lines <- readLines(legend_path, warn = FALSE)
leg_lines <- leg_lines[grepl("^\\s*\\d+\\s*:\\s*", leg_lines)]

m <- regexec("^\\s*(\\d+)\\s*:\\s*([A-Za-z]{2,3})\\s+([^\\[]+)", leg_lines)
parts <- regmatches(leg_lines, m)

legend_dt <- rbindlist(lapply(parts, function(x) {
  if (length(x) < 4) return(NULL)
  data.table(
    kg_id = as.integer(x[2]),
    koppen_code = x[3],
    koppen_desc = trimws(x[4])
  )
}), use.names = TRUE, fill = TRUE)

stopifnot(nrow(legend_dt) > 0)
say("Legend rows: %d", nrow(legend_dt))

# raster -> data.frame for ggplot
koppen_df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
val_col <- setdiff(names(koppen_df), c("x", "y"))[1]
setDT(koppen_df)
setnames(koppen_df, val_col, "kg_id")
koppen_df[, kg_id := as.integer(kg_id)]

koppen_df <- merge(koppen_df, legend_dt, by = "kg_id", all.x = TRUE)

say("Raster rows kept for plotting: %s", format(nrow(koppen_df), big.mark = ","))

# estimate tile size for geom_tile
ux <- sort(unique(koppen_df$x))
uy <- sort(unique(koppen_df$y))
tile_w <- if (length(ux) > 1) median(diff(ux)) else 0.00833333
tile_h <- if (length(uy) > 1) median(diff(uy)) else 0.00833333

# ------------------------------------------------------------
# 3) DEFINE STUDY REGIONS
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# 4) WORLDMAP
# ------------------------------------------------------------
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# ------------------------------------------------------------
# 5) COLOR PALETTE
# ------------------------------------------------------------
koppen_levels <- legend_dt[order(kg_id)]$koppen_code

koppen_colors <- c(
  Af  = "#006d2c", Am  = "#31a354", Aw  = "#74c476",
  BWh = "#8c2d04", BWk = "#d94801", BSh = "#f16913", BSk = "#fdae6b",
  Csa = "#08519c", Csb = "#3182bd", Csc = "#6baed6",
  Cwa = "#2171b5", Cwb = "#4292c6", Cwc = "#9ecae1",
  Cfa = "#084594", Cfb = "#4a90d9", Cfc = "#bdd7e7",
  Dsa = "#54278f", Dsb = "#756bb1", Dsc = "#9e9ac8", Dsd = "#cbc9e2",
  Dwa = "#6a51a3", Dwb = "#807dba", Dwc = "#9e9ac8", Dwd = "#dadaeb",
  Dfa = "#88419d", Dfb = "#8c6bb1", Dfc = "#bcbddc", Dfd = "#efedf5",
  ET  = "#969696", EF  = "#d9d9d9"
)

koppen_colors <- koppen_colors[names(koppen_colors) %in% koppen_levels]

missing_codes <- setdiff(koppen_levels, names(koppen_colors))
if (length(missing_codes) > 0) {
  koppen_colors <- c(
    koppen_colors,
    setNames(rep("#000000", length(missing_codes)), missing_codes)
  )
}
koppen_colors <- koppen_colors[koppen_levels]

# ------------------------------------------------------------
# 6) GLOBAL OVERVIEW PANEL
# ------------------------------------------------------------
p_overview <- ggplot() +
  geom_tile(
    data = koppen_df,
    aes(x = x, y = y, fill = factor(koppen_code, levels = koppen_levels)),
    width = tile_w,
    height = tile_h
  ) +
  geom_sf(data = world, fill = NA, color = "grey35", linewidth = 0.2) +
  geom_sf(data = region_boxes, fill = NA, color = "black", linewidth = 0.7) +
  geom_text(
    data = region_labels,
    aes(x = x, y = y, label = label),
    size = 3,
    fontface = "bold"
  ) +
  scale_fill_manual(values = koppen_colors, drop = FALSE, name = "Köppen class") +
  coord_sf(xlim = c(-140, 100), ylim = c(-60, 60), expand = FALSE) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank()
  ) +
  labs(
    title = "A) Study regions in the Köppen–Geiger climate classification"
  )

# ------------------------------------------------------------
# 7) REGIONAL PANELS
#    Keep legend only on Africa panel to avoid patchwork issues
# ------------------------------------------------------------
p_conus <- make_region_plot(
  title_text = "CONUS",
  xlim = c(-125, -66),
  ylim = c(24, 50),
  show_legend = FALSE
)

p_india <- make_region_plot(
  title_text = "India",
  xlim = c(68, 90),
  ylim = c(6, 36),
  show_legend = FALSE
)

p_sa <- make_region_plot(
  title_text = "South America",
  xlim = c(-82, -35),
  ylim = c(-56, 13),
  show_legend = FALSE
)

p_africa <- make_region_plot(
  title_text = "Africa",
  xlim = c(-20, 55),
  ylim = c(-35, 38),
  show_legend = TRUE
)

p_regions <- (p_conus | p_india) / (p_sa | p_africa)

# ------------------------------------------------------------
# 8) COMBINE + SAVE
# ------------------------------------------------------------
final_plot <- p_overview / p_regions +
  plot_layout(heights = c(0.9, 1.6))

ggsave(
  filename = out_png,
  plot = final_plot,
  width = 14,
  height = 11,
  dpi = 300,
  bg = "white"
)

say("Wrote figure:")
say("%s", out_png)

# ------------------------------------------------------------
# 9) OPTIONAL: SAVE THE TWO FIGURES SEPARATELY
# ------------------------------------------------------------
out_png_overview <- "/Users/priyanka/Downloads/PurpleAir_EPA/fig_climate_zones_overview.png"
out_png_regions  <- "/Users/priyanka/Downloads/PurpleAir_EPA/fig_climate_zones_regions.png"

ggsave(
  filename = out_png_overview,
  plot = p_overview,
  width = 13,
  height = 7.5,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = out_png_regions,
  plot = p_regions,
  width = 14,
  height = 10,
  dpi = 300,
  bg = "white"
)

say("Wrote separate overview figure:")
say("%s", out_png_overview)

say("Wrote separate regional figure:")
say("%s", out_png_regions)
