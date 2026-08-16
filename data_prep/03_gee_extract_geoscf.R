# ==========================================================================
# 03_gee_extract_geoscf.R
#
# Google Earth Engine (rgee): export hourly GEOS-CF PM2.5 at each sensor location to GCS.
#
# Extracted VERBATIM from purpleair_degradation2.Rmd (chunks #15 lines 2740-2931, #16 lines 2934-3012).
# Code is unmodified; see REPRODUCING.md for paths and run order.
# ==========================================================================

# ---- source chunk #15 (Rmd lines 2740-2931) ------------------------

# ============================================================
# GEMINI PM2.5 (GEOS-CF) — PREFIXED EXPORT + KEEP LAT/LON
# - New files go to: run=20260206_locid/GEOSCF_PM25_YYYYMM_batchXX.csv
# - Fast existence check using prefix= (no full-bucket listing)
# - Keeps longitude/latitude as columns in the exported CSV
# ============================================================

library(rgee)
library(dplyr)
library(lubridate)
library(googleCloudStorageR)

my_gcs_bucket <- "YOUR-GCS-BUCKET"
run_tag <- "run=20260206_locid"   # <-- change per rerun
ee_Initialize(project = "YOUR-GEE-PROJECT-ID", gcs = TRUE)

# -----------------------------
# 1) Load + prep points (KEEP LON/LAT AS PROPERTIES)
# -----------------------------
load("/Users/priyanka/Downloads/PurpleAir_EPA/annual_data_all_outdoor_sensorlocations.RData")
colnames(outdoor) <- c("longitude", "latitude", "station")

outdoor <- outdoor %>%
  mutate(
    pid       = row_number(),
    longitude = as.numeric(longitude),
    latitude  = as.numeric(latitude),
    station   = as.character(station)
  ) %>%
  select(pid, station, longitude, latitude)

# -----------------------------
# 2) Config
# -----------------------------
points_per_batch <- 10000
max_active_tasks <- 25
sleep_sec        <- 120

start_date    <- as.Date("2018-01-01")
end_date      <- as.Date("2025-10-01")
window_starts <- seq(start_date, end_date, by = "1 month")

# -----------------------------
# 3) Dataset helpers (PM2.5 only)
# -----------------------------
get_dataset <- function(name) {
  if (name == "GEOSCF_PM25") {
    ee$ImageCollection("NASA/GEOS-CF/v1/rpl/tavg1hr")$select("PM25_RH35_GCC")
  } else stop("Unknown dataset")
}

# -----------------------------
# 4) Native Scale Helper (Cached)
# -----------------------------
scale_cache <- new.env(parent = emptyenv())
get_native_scale <- function(dataset_name) {
  if (exists(dataset_name, envir = scale_cache)) return(get(dataset_name, envir = scale_cache))
  ic <- get_dataset(dataset_name)
  sc <- as.numeric(ic$first()$projection()$nominalScale()$getInfo())
  assign(dataset_name, sc, envir = scale_cache)
  sc
}

# -----------------------------
# 5) Batching + FeatureCollection
# -----------------------------
points_batches <- split(outdoor, ceiling(seq_len(nrow(outdoor)) / points_per_batch))
points_batches <- lapply(seq_along(points_batches), function(i) {
  x <- points_batches[[i]]; attr(x, "batch_id") <- i; x
})

# IMPORTANT: include lon/lat as properties so they appear in CSV
to_ee_points_manual <- function(df_batch) {
  feats <- lapply(seq_len(nrow(df_batch)), function(i) {
    ee$Feature(
      ee$Geometry$Point(list(df_batch$longitude[i], df_batch$latitude[i])),
      list(
        pid     = as.integer(df_batch$pid[i]),
        station = as.character(df_batch$station[i]),
        lon     = as.numeric(df_batch$longitude[i]),
        lat     = as.numeric(df_batch$latitude[i])
      )
    )
  })
  ee$FeatureCollection(feats)
}

# -----------------------------
# 6) Robust Task Management
# -----------------------------
task_state_by_desc <- function(desc) {
  tasks <- ee_manage_task()
  if (length(tasks) == 0 || (is.data.frame(tasks) && nrow(tasks) == 0)) return(NA_character_)
  if (is.data.frame(tasks)) {
    match <- tasks[tasks$description == desc, ]
    return(if (nrow(match) > 0) as.character(match$state[1]) else NA_character_)
  }
  NA_character_
}

wait_for_slot <- function(max_active, sleep_sec) {
  repeat {
    tasks <- ee_manage_task()
    active <- if (is.data.frame(tasks)) sum(tasks$state %in% c("READY", "RUNNING")) else 0
    cat(sprintf("[%s] Active: %d\n", format(Sys.time(), "%H:%M:%S"), active))
    if (active < max_active) break
    Sys.sleep(sleep_sec)
  }
}

# -----------------------------
# 7) RESUMABLE Export Function (Option 1 prefix folder)
# -----------------------------
export_monthly_mapped_gcs <- function(dataset_name, df_batch, t_start, t_end, bucket, run_tag) {

  b_id <- attr(df_batch, "batch_id")
  desc <- sprintf("%s_%s_batch%02d", dataset_name, format(t_start, "%Y%m"), b_id)

  # NEW: object lives under run_tag prefix
  object_name <- paste0(run_tag, "/", desc, ".csv")
  prefix_key  <- paste0(run_tag, "/", desc)  # fast existence check

  # CHECK 1: does this exact object already exist?
  bucket_files <- gcs_list_objects(bucket, prefix = prefix_key)
  if (nrow(bucket_files) > 0 && any(bucket_files$name == object_name)) {
    message("✅ Skipping (File already in GCS): ", object_name)
    return(invisible(NULL))
  }

  # CHECK 2: is a task already active with this description?
  st <- task_state_by_desc(desc)
  if (!is.na(st) && st %in% c("READY", "RUNNING")) {
    message("⏳ Skipping (Task currently active): ", desc)
    return(invisible(NULL))
  }

  dataset   <- get_dataset(dataset_name)
  scale     <- get_native_scale(dataset_name)
  ee_points <- to_ee_points_manual(df_batch)

  img_col <- dataset$filterDate(as.character(t_start), as.character(t_end))

  extract_fn <- function(img) {
    t_str <- img$date()$format("YYYY-MM-dd HH:mm:ss")

    img$reduceRegions(
      collection = ee_points,
      reducer    = ee$Reducer$first(),
      scale      = scale
    )$
      # keep properties (pid/station/lon/lat) and add datetime
      map(function(f) f$set("datetime", t_str)$setGeometry(NULL))
  }

  samples <- img_col$map(extract_fn)$flatten()

  message("🚀 Queueing: ", desc, " -> ", run_tag, "/")
  wait_for_slot(max_active_tasks, sleep_sec)

  task <- ee_table_to_gcs(
    collection     = samples,
    description    = desc,
    bucket         = bucket,
    fileNamePrefix = paste0(run_tag, "/", desc),  # <-- Option 1
    fileFormat     = "CSV"
  )
  task$start()
}

# -----------------------------
# SESSION: PM2.5 ONLY (GEOS-CF)
# -----------------------------
for (i in seq_along(window_starts)) {
  t_start <- window_starts[i]
  t_end   <- seq(t_start, length = 2, by = "1 month")[2]
  if (t_end > end_date) t_end <- end_date
  if (t_start >= end_date) next

  for (b in seq_along(points_batches)) {
    cat(sprintf("\nPM2.5 | %s | Batch %d\n", t_start, b))
    export_monthly_mapped_gcs("GEOSCF_PM25", points_batches[[b]], t_start, t_end,
                              my_gcs_bucket, run_tag)
    Sys.sleep(1)
  }
}

# Optional: list only NEW files
# objs_new <- gcs_list_objects(my_gcs_bucket, prefix = paste0(run_tag, "/"))
# tail(objs_new$name)

# ---- source chunk #16 (Rmd lines 2934-3012) ------------------------

library(rgee)
library(dplyr)
library(lubridate)
library(googleCloudStorageR)

my_gcs_bucket <- "YOUR-GCS-BUCKET"
run_tag <- "run=20260206_locid"   # <-- change per rerun
ee_Initialize(project = "YOUR-GEE-PROJECT-ID", gcs = TRUE)

# ----------------------------
# 0.5) DOWNLOAD FROM GCS (run_tag prefix) -> local ERA/PM folders
# ----------------------------
suppressPackageStartupMessages({
  library(googleCloudStorageR)
  library(stringr)
})

my_gcs_bucket <- "YOUR-GCS-BUCKET"
run_tag <- "run=20260206_locid"   # <-- MUST match your export run_tag
gcs_prefix <- paste0(run_tag, "/")

# Make sure you are authenticated (use whichever you already use)
# If you already have googleCloudStorageR working, you can skip re-auth.
# gcs_auth()   # or gcs_auth("path/to/your-service-account.json")

# Local dirs you already use
base_dir     <- "/Users/priyanka/Downloads/PurpleAir_EPA/GCS_exports"
era_csv_dir  <- file.path(base_dir, "ERA")
pm_csv_dir   <- file.path(base_dir, "PM")
dir.create(era_csv_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pm_csv_dir,  recursive = TRUE, showWarnings = FALSE)

# List ALL CSV objects under the run prefix (may take a bit if huge)
objs <- gcs_list_objects(my_gcs_bucket, prefix = gcs_prefix)

if (is.null(objs) || nrow(objs) == 0) stop("No objects found in bucket under prefix: ", gcs_prefix)

# Keep CSVs only
objs <- objs[grepl("\\.csv$", objs$name), ]
if (nrow(objs) == 0) stop("No CSVs found under prefix: ", gcs_prefix)

# Route to ERA vs PM based on filename
fname <- basename(objs$name)
is_era <- str_detect(fname, "^ERA5_")
is_pm  <- str_detect(fname, "^GEOSCF_PM25_")

objs_keep <- objs[is_era | is_pm, ]
if (nrow(objs_keep) == 0) stop("No ERA5_ or GEOSCF_PM25_ CSVs found under prefix: ", gcs_prefix)

download_one <- function(object_name) {
  f <- basename(object_name)

  out_dir <- if (str_detect(f, "^ERA5_")) era_csv_dir else pm_csv_dir
  out_path <- file.path(out_dir, f)

  if (file.exists(out_path) && file.info(out_path)$size > 0) {
    message("✅ already local: ", f)
    return(invisible(TRUE))
  }

  message("⬇️ downloading: ", object_name, " -> ", out_path)
  # This writes the object to disk
  googleCloudStorageR::gcs_get_object(
    object_name,
    bucket = my_gcs_bucket,
    saveToDisk = out_path,
    overwrite = TRUE
  )
  TRUE
}

# Download loop (serial; safest)
for (k in seq_len(nrow(objs_keep))) {
  download_one(objs_keep$name[k])
}
cat("✅ GCS download complete.\n")
