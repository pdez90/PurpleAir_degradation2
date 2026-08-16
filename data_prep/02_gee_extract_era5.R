# ==========================================================================
# 02_gee_extract_era5.R
#
# Google Earth Engine (rgee): export hourly ERA5-Land 2 m temperature and dewpoint at each
# sensor location to Google Cloud Storage, in monthly batches. Restartable.
#
# Extracted VERBATIM from purpleair_degradation2.Rmd (chunks #11 lines 2004-2172, #12 lines 2175-2360, #13 lines 2363-2548, #14 lines 2551-2737).
# Code is unmodified; see REPRODUCING.md for paths and run order.
# ==========================================================================

# ---- source chunk #11 (Rmd lines 2004-2172) ------------------------

# ============================================================
# FINAL ROBUST: GCS EXPORT + MAPPED HOURLY EXTRACTION (rgee)
# Fixed: "Computed value too large" by avoiding toBands()
# ============================================================

library(rgee)
library(dplyr)
library(lubridate)

my_gcs_bucket <- "YOUR-GCS-BUCKET" 
ee_Initialize(project = "YOUR-GEE-PROJECT-ID", gcs = TRUE)

# -----------------------------
# 1) Load + prep points
# -----------------------------
load("/Users/priyanka/Downloads/PurpleAir_EPA/annual_data_all_outdoor_sensorlocations.RData")
colnames(outdoor) <- c("longitude", "latitude", "station")

outdoor <- outdoor %>%
  dplyr::mutate(
    pid       = dplyr::row_number(),
    longitude = as.numeric(longitude),
    latitude  = as.numeric(latitude),
    station   = as.character(station)
  ) %>%
  dplyr::select(pid, station, longitude, latitude)

# -----------------------------
# 2) Config (Monthly Windows, Balanced Batches)
# -----------------------------
points_per_batch     <- 10000  # Dropped to 10k for memory stability
max_active_tasks     <- 25     # Increase parallelism
sleep_sec            <- 120    

datasets_met <- c("ERA5", "NLDAS")
datasets_pm  <- c("GEOSCF_PM25")

start_date   <- as.Date("2018-01-01")
end_date     <- as.Date("2025-10-01")
window_starts <- seq(start_date, end_date, by = "1 month")

# -----------------------------
# 3) Dataset helpers
# -----------------------------
get_dataset <- function(name) {
  if (name == "ERA5") {
    ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY")$
      select(c("temperature_2m", "surface_pressure", "dewpoint_temperature_2m"))
  } else if (name == "NLDAS") {
    ee$ImageCollection("NASA/NLDAS/FORA0125_H002")$
      select(c("temperature", "pressure", "specific_humidity"))
  } else if (name == "GEOSCF_PM25") {
    ee$ImageCollection("NASA/GEOS-CF/v1/rpl/tavg1hr")$
      select("PM25_RH35_GCC")
  } else { stop("Unknown dataset") }
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

to_ee_points_manual <- function(df_batch) {
  feats <- lapply(seq_len(nrow(df_batch)), function(i) {
    ee$Feature(ee$Geometry$Point(list(df_batch$longitude[i], df_batch$latitude[i])), 
               list(pid = as.integer(df_batch$pid[i]), station = df_batch$station[i]))
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
    return(if(nrow(match) > 0) as.character(match$state[1]) else NA_character_)
  }
  return(NA_character_)
}

wait_for_slot <- function(max_active, sleep_sec) {
  repeat {
    tasks <- ee_manage_task()
    active <- if(is.data.frame(tasks)) sum(tasks$state %in% c("READY", "RUNNING")) else 0
    cat(sprintf("[%s] Active: %d\n", format(Sys.time(), "%H:%M:%S"), active))
    if (active < max_active) break
    Sys.sleep(sleep_sec)
  }
}

# -----------------------------
# 7) Final Export Function: Mapped Hourly Collection
# -----------------------------
export_monthly_mapped_gcs <- function(dataset_name, df_batch, t_start, t_end, bucket) {
  
  b_id <- attr(df_batch, "batch_id")
  desc <- sprintf("%s_%s_batch%02d", dataset_name, format(t_start, "%Y%m"), b_id)
  
  st <- task_state_by_desc(desc)
  if (!is.na(st) && st %in% c("COMPLETED", "READY", "RUNNING")) return(NULL)

  dataset   <- get_dataset(dataset_name)
  scale     <- get_native_scale(dataset_name)
  ee_points <- to_ee_points_manual(df_batch)

  # MAPPING: Processes images one-by-one (Memory Efficient)
  img_col <- dataset$filterDate(as.character(t_start), as.character(t_end))
  
  extract_fn <- function(img) {
    t_str <- img$date()$format("YYYY-MM-dd HH:mm:ss")
    img$reduceRegions(
      collection = ee_points,
      reducer = ee$Reducer$first(),
      scale = scale
    )$map(function(f) f$set("datetime", t_str)$setGeometry(NULL))
  }

  samples <- img_col$map(extract_fn)$flatten()

  message("Queueing: ", desc)
  wait_for_slot(max_active_tasks, sleep_sec)

  task <- ee_table_to_gcs(
    collection = samples,
    description = desc,
    bucket = bucket,
    fileNamePrefix = desc,
    fileFormat = "CSV"
  )
  task$start()
}

# -----------------------------
# 8) EXECUTION
# -----------------------------
for (ds in c(datasets_met, datasets_pm)) {
  for (i in seq_along(window_starts)) {
    t_start <- window_starts[i]
    t_end   <- seq(t_start, length = 2, by = "1 month")[2]
    if (t_end > end_date) t_end <- end_date
    if (t_start >= end_date) next

    for (b in seq_along(points_batches)) {
      cat(sprintf("\n%s | %s | Batch %d\n", ds, t_start, b))
      export_monthly_mapped_gcs(ds, points_batches[[b]], t_start, t_end, my_gcs_bucket)
      Sys.sleep(2) 
    }
  }
}

# ---- source chunk #12 (Rmd lines 2175-2360) ------------------------

# ============================================================
# FINAL ROBUST: GCS EXPORT + MAPPED HOURLY EXTRACTION (rgee)
# Fixed: "Computed value too large" by avoiding toBands()
# ============================================================

library(rgee)
library(dplyr)
library(lubridate)

my_gcs_bucket <- "YOUR-GCS-BUCKET" 
ee_Initialize(project = "YOUR-GEE-PROJECT-ID", gcs = TRUE)

# -----------------------------
# 1) Load + prep points
# -----------------------------
load("/Users/priyanka/Downloads/PurpleAir_EPA/annual_data_all_outdoor_sensorlocations.RData")
colnames(outdoor) <- c("longitude", "latitude", "station")

outdoor <- outdoor %>%
  dplyr::mutate(
    pid       = dplyr::row_number(),
    longitude = as.numeric(longitude),
    latitude  = as.numeric(latitude),
    station   = as.character(station)
  ) %>%
  dplyr::select(pid, station, longitude, latitude)

# -----------------------------
# 2) Config (Monthly Windows, Balanced Batches)
# -----------------------------
points_per_batch     <- 10000  # Dropped to 10k for memory stability
max_active_tasks     <- 25     # Increase parallelism
sleep_sec            <- 120    

datasets_met <- c("ERA5", "NLDAS")
datasets_pm  <- c("GEOSCF_PM25")

start_date   <- as.Date("2018-01-01")
end_date     <- as.Date("2025-10-01")
window_starts <- seq(start_date, end_date, by = "1 month")

# -----------------------------
# 3) Dataset helpers
# -----------------------------
get_dataset <- function(name) {
  if (name == "ERA5") {
    ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY")$
      select(c("temperature_2m", "surface_pressure", "dewpoint_temperature_2m"))
  } else if (name == "NLDAS") {
    ee$ImageCollection("NASA/NLDAS/FORA0125_H002")$
      select(c("temperature", "pressure", "specific_humidity"))
  } else if (name == "GEOSCF_PM25") {
    ee$ImageCollection("NASA/GEOS-CF/v1/rpl/tavg1hr")$
      select("PM25_RH35_GCC")
  } else { stop("Unknown dataset") }
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

to_ee_points_manual <- function(df_batch) {
  feats <- lapply(seq_len(nrow(df_batch)), function(i) {
    ee$Feature(ee$Geometry$Point(list(df_batch$longitude[i], df_batch$latitude[i])), 
               list(pid = as.integer(df_batch$pid[i]), station = df_batch$station[i]))
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
    return(if(nrow(match) > 0) as.character(match$state[1]) else NA_character_)
  }
  return(NA_character_)
}

wait_for_slot <- function(max_active, sleep_sec) {
  repeat {
    tasks <- ee_manage_task()
    active <- if(is.data.frame(tasks)) sum(tasks$state %in% c("READY", "RUNNING")) else 0
    cat(sprintf("[%s] Active: %d\n", format(Sys.time(), "%H:%M:%S"), active))
    if (active < max_active) break
    Sys.sleep(sleep_sec)
  }
}

# -----------------------------
# 7) Final Export Function: Mapped Hourly Collection
# -----------------------------
library(googleCloudStorageR)

# -----------------------------
# 7) RESUMABLE Export Function: GCS File Check + Task History
# -----------------------------
export_monthly_mapped_gcs <- function(dataset_name, df_batch, t_start, t_end, bucket) {
  
  b_id <- attr(df_batch, "batch_id")
  desc <- sprintf("%s_%s_batch%02d", dataset_name, format(t_start, "%Y%m"), b_id)
  file_name <- paste0(desc, ".csv")

  # --- CHECK 1: Does the file exist in the GCS bucket? ---
  # This works even years later after GEE history is deleted.
  bucket_files <- gcs_list_objects(bucket)
  if (nrow(bucket_files) > 0 && file_name %in% bucket_files$name) {
    message("✅ Skipping (File already in GCS): ", file_name)
    return(NULL)
  }

  # --- CHECK 2: Is the task currently Running or Ready? ---
  st <- task_state_by_desc(desc)
  if (!is.na(st) && st %in% c("READY", "RUNNING")) {
    message("⏳ Skipping (Task currently active): ", desc)
    return(NULL)
  }

  # --- EXECUTION: Only runs if File NOT in GCS and Task NOT active ---
  dataset   <- get_dataset(dataset_name)
  scale     <- get_native_scale(dataset_name)
  ee_points <- to_ee_points_manual(df_batch)

  img_col <- dataset$filterDate(as.character(t_start), as.character(t_end))
  
  extract_fn <- function(img) {
    t_str <- img$date()$format("YYYY-MM-dd HH:mm:ss")
    img$reduceRegions(
      collection = ee_points,
      reducer = ee$Reducer$first(),
      scale = scale
    )$map(function(f) f$set("datetime", t_str)$setGeometry(NULL))
  }

  samples <- img_col$map(extract_fn)$flatten()

  message("🚀 Queueing: ", desc)
  wait_for_slot(max_active_tasks, sleep_sec)

  task <- ee_table_to_gcs(
    collection = samples,
    description = desc,
    bucket = bucket,
    fileNamePrefix = desc,
    fileFormat = "CSV"
  )
  task$start()
}

# -----------------------------
# SESSION 1: ERA5 ONLY
# -----------------------------
target_ds <- "ERA5"
for (i in seq_along(window_starts)) {
  t_start <- window_starts[i]
  t_end   <- seq(t_start, length = 2, by = "1 month")[2]
  if (t_end > end_date) t_end <- end_date
  if (t_start >= end_date) next

  for (b in seq_along(points_batches)) {
    cat(sprintf("\nERA5 | %s | Batch %d\n", t_start, b))
    export_monthly_mapped_gcs("ERA5", points_batches[[b]], t_start, t_end, my_gcs_bucket)
    Sys.sleep(1) # Faster submission for single dataset
  }
}

# ---- source chunk #13 (Rmd lines 2363-2548) ------------------------

# ============================================================
# FINAL ROBUST: GCS EXPORT + MAPPED HOURLY EXTRACTION (rgee)
# Fixed: "Computed value too large" by avoiding toBands()
# ============================================================

library(rgee)
library(dplyr)
library(lubridate)

my_gcs_bucket <- "YOUR-GCS-BUCKET" 
ee_Initialize(project = "YOUR-GEE-PROJECT-ID", gcs = TRUE)

# -----------------------------
# 1) Load + prep points
# -----------------------------
load("/Users/priyanka/Downloads/PurpleAir_EPA/annual_data_all_outdoor_sensorlocations.RData")
colnames(outdoor) <- c("longitude", "latitude", "station")

outdoor <- outdoor %>%
  dplyr::mutate(
    pid       = dplyr::row_number(),
    longitude = as.numeric(longitude),
    latitude  = as.numeric(latitude),
    station   = as.character(station)
  ) %>%
  dplyr::select(pid, station, longitude, latitude)

# -----------------------------
# 2) Config (Monthly Windows, Balanced Batches)
# -----------------------------
points_per_batch     <- 10000  # Dropped to 10k for memory stability
max_active_tasks     <- 25     # Increase parallelism
sleep_sec            <- 120    

datasets_met <- c("ERA5", "NLDAS")
datasets_pm  <- c("GEOSCF_PM25")

start_date   <- as.Date("2018-01-01")
end_date     <- as.Date("2025-10-01")
window_starts <- seq(start_date, end_date, by = "1 month")

# -----------------------------
# 3) Dataset helpers
# -----------------------------
get_dataset <- function(name) {
  if (name == "ERA5") {
    ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY")$
      select(c("temperature_2m", "surface_pressure", "dewpoint_temperature_2m"))
  } else if (name == "NLDAS") {
    ee$ImageCollection("NASA/NLDAS/FORA0125_H002")$
      select(c("temperature", "pressure", "specific_humidity"))
  } else if (name == "GEOSCF_PM25") {
    ee$ImageCollection("NASA/GEOS-CF/v1/rpl/tavg1hr")$
      select("PM25_RH35_GCC")
  } else { stop("Unknown dataset") }
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

to_ee_points_manual <- function(df_batch) {
  feats <- lapply(seq_len(nrow(df_batch)), function(i) {
    ee$Feature(ee$Geometry$Point(list(df_batch$longitude[i], df_batch$latitude[i])), 
               list(pid = as.integer(df_batch$pid[i]), station = df_batch$station[i]))
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
    return(if(nrow(match) > 0) as.character(match$state[1]) else NA_character_)
  }
  return(NA_character_)
}

wait_for_slot <- function(max_active, sleep_sec) {
  repeat {
    tasks <- ee_manage_task()
    active <- if(is.data.frame(tasks)) sum(tasks$state %in% c("READY", "RUNNING")) else 0
    cat(sprintf("[%s] Active: %d\n", format(Sys.time(), "%H:%M:%S"), active))
    if (active < max_active) break
    Sys.sleep(sleep_sec)
  }
}

# -----------------------------
# 7) Final Export Function: Mapped Hourly Collection
# -----------------------------
library(googleCloudStorageR)

# -----------------------------
# 7) RESUMABLE Export Function: GCS File Check + Task History
# -----------------------------
export_monthly_mapped_gcs <- function(dataset_name, df_batch, t_start, t_end, bucket) {
  
  b_id <- attr(df_batch, "batch_id")
  desc <- sprintf("%s_%s_batch%02d", dataset_name, format(t_start, "%Y%m"), b_id)
  file_name <- paste0(desc, ".csv")

  # --- CHECK 1: Does the file exist in the GCS bucket? ---
  # This works even years later after GEE history is deleted.
  bucket_files <- gcs_list_objects(bucket)
  if (nrow(bucket_files) > 0 && file_name %in% bucket_files$name) {
    message("✅ Skipping (File already in GCS): ", file_name)
    return(NULL)
  }

  # --- CHECK 2: Is the task currently Running or Ready? ---
  st <- task_state_by_desc(desc)
  if (!is.na(st) && st %in% c("READY", "RUNNING")) {
    message("⏳ Skipping (Task currently active): ", desc)
    return(NULL)
  }

  # --- EXECUTION: Only runs if File NOT in GCS and Task NOT active ---
  dataset   <- get_dataset(dataset_name)
  scale     <- get_native_scale(dataset_name)
  ee_points <- to_ee_points_manual(df_batch)

  img_col <- dataset$filterDate(as.character(t_start), as.character(t_end))
  
  extract_fn <- function(img) {
    t_str <- img$date()$format("YYYY-MM-dd HH:mm:ss")
    img$reduceRegions(
      collection = ee_points,
      reducer = ee$Reducer$first(),
      scale = scale
    )$map(function(f) f$set("datetime", t_str)$setGeometry(NULL))
  }

  samples <- img_col$map(extract_fn)$flatten()

  message("🚀 Queueing: ", desc)
  wait_for_slot(max_active_tasks, sleep_sec)

  task <- ee_table_to_gcs(
    collection = samples,
    description = desc,
    bucket = bucket,
    fileNamePrefix = desc,
    fileFormat = "CSV"
  )
  task$start()
}

# -----------------------------
# SESSION 3: PM2.5 ONLY
# -----------------------------
target_ds <- "GEOSCF_PM25"
for (i in seq_along(window_starts)) {
  t_start <- window_starts[i]
  t_end   <- seq(t_start, length = 2, by = "1 month")[2]
  if (t_end > end_date) t_end <- end_date
  if (t_start >= end_date) next

  for (b in seq_along(points_batches)) {
    cat(sprintf("\nPM2.5 | %s | Batch %d\n", t_start, b))
    export_monthly_mapped_gcs("GEOSCF_PM25", points_batches[[b]], t_start, t_end, my_gcs_bucket)
    Sys.sleep(1)
  }
}

# ---- source chunk #14 (Rmd lines 2551-2737) ------------------------

# ============================================================
# FINAL ROBUST: GCS EXPORT + MAPPED HOURLY EXTRACTION (rgee)
# UPDATED (Option 1): write NEW exports under a GCS prefix folder
#  - All new files go to:  run=20260206_locid/<desc>.csv
#  - Existence check uses prefix= (fast; does NOT list whole bucket)
# ============================================================

library(rgee)
library(dplyr)
library(lubridate)
library(googleCloudStorageR)

my_gcs_bucket <- "YOUR-GCS-BUCKET"
run_tag <- "run=20260206_locid"   # <-- change this when you start a new run
ee_Initialize(project = "YOUR-GEE-PROJECT-ID", gcs = TRUE)

# -----------------------------
# 1) Load + prep points
# -----------------------------
load("/Users/priyanka/Downloads/PurpleAir_EPA/annual_data_all_outdoor_sensorlocations.RData")
colnames(outdoor) <- c("longitude", "latitude", "station")

outdoor <- outdoor %>%
  dplyr::mutate(
    pid       = dplyr::row_number(),
    longitude = as.numeric(longitude),
    latitude  = as.numeric(latitude),
    station   = as.character(station)
  ) %>%
  dplyr::select(pid, station, longitude, latitude)

# -----------------------------
# 2) Config (Monthly Windows, Balanced Batches)
# -----------------------------
points_per_batch <- 10000
max_active_tasks <- 25
sleep_sec        <- 120

start_date    <- as.Date("2018-01-01")
end_date      <- as.Date("2025-10-01")
window_starts <- seq(start_date, end_date, by = "1 month")

# -----------------------------
# 3) Dataset helpers
# -----------------------------
get_dataset <- function(name) {
  if (name == "ERA5") {
    ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY")$
      select(c("temperature_2m", "surface_pressure", "dewpoint_temperature_2m"))
  } else if (name == "NLDAS") {
    ee$ImageCollection("NASA/NLDAS/FORA0125_H002")$
      select(c("temperature", "pressure", "specific_humidity"))
  } else if (name == "GEOSCF_PM25") {
    ee$ImageCollection("NASA/GEOS-CF/v1/rpl/tavg1hr")$
      select("PM25_RH35_GCC")
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

to_ee_points_manual <- function(df_batch) {
  feats <- lapply(seq_len(nrow(df_batch)), function(i) {
    ee$Feature(
      ee$Geometry$Point(list(df_batch$longitude[i], df_batch$latitude[i])),
      list(pid = as.integer(df_batch$pid[i]), station = as.character(df_batch$station[i]))
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

  # NEW: write into prefix folder
  object_name <- paste0(run_tag, "/", desc, ".csv")
  prefix_key  <- paste0(run_tag, "/", desc)  # used for fast existence check

  # --- CHECK 1: Does the NEW object already exist? (fast prefix listing) ---
  bucket_files <- gcs_list_objects(bucket, prefix = prefix_key)
  if (nrow(bucket_files) > 0 && any(bucket_files$name == object_name)) {
    message("✅ Skipping (File already in GCS): ", object_name)
    return(invisible(NULL))
  }

  # --- CHECK 2: Is the task currently Running or Ready? ---
  st <- task_state_by_desc(desc)
  if (!is.na(st) && st %in% c("READY", "RUNNING")) {
    message("⏳ Skipping (Task currently active): ", desc)
    return(invisible(NULL))
  }

  # --- EXECUTION ---
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
    )$map(function(f) f$set("datetime", t_str)$setGeometry(NULL))
  }

  samples <- img_col$map(extract_fn)$flatten()

  message("🚀 Queueing: ", desc, " -> ", run_tag, "/")
  wait_for_slot(max_active_tasks, sleep_sec)

  task <- ee_table_to_gcs(
    collection     = samples,
    description    = desc,                 # keep short for task mgmt
    bucket         = bucket,
    fileNamePrefix = paste0(run_tag, "/", desc),  # <-- Option 1
    fileFormat     = "CSV"
  )
  task$start()
}

# -----------------------------
# SESSION 1: ERA5 ONLY
# -----------------------------
for (i in seq_along(window_starts)) {
  t_start <- window_starts[i]
  t_end   <- seq(t_start, length = 2, by = "1 month")[2]
  if (t_end > end_date) t_end <- end_date
  if (t_start >= end_date) next

  for (b in seq_along(points_batches)) {
    cat(sprintf("\nERA5 | %s | Batch %d\n", t_start, b))
    export_monthly_mapped_gcs("ERA5", points_batches[[b]], t_start, t_end,
                              my_gcs_bucket, run_tag)
    Sys.sleep(1)
  }
}
