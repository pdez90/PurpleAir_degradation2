# ==========================================================================
# 04_join_reference_parquet.R
#
# Post-GCS-download pipeline: local ERA5/GEOS-CF CSVs -> parquet -> hourly reference dataset
# -> joined to outdoor PurpleAir observations (monthly joined RDS files).
#
# Extracted VERBATIM from purpleair_degradation2.Rmd (chunk #18 lines 3813-4649).
# Code is unmodified; see REPRODUCING.md for paths and run order.
# ==========================================================================

# ---- source chunk #18 (Rmd lines 3813-4649) ------------------------

# ============================================================
# ONE-STOP PIPELINE AFTER GCS DOWNLOAD (UPDATED VERSION: PID FIX)
#   Local ERA/PM CSVs -> Parquet -> gee_hourly_ds -> join to outdoor monthly
#
# CRITICAL FIX:
# - ERA and PM are merged by: pid + datetime_hr
#   NOT by station + datetime_hr
# - This avoids many-to-many row explosion when multiple pids share a station.
#
# YOUR RUN:
# - ERA5 export kept: pid, station, datetime
# - PM2.5 export kept: pid, station, datetime, lon, lat
#
# This pipeline keeps:
#   era_pid, pm_pid, gee_pid
#
# Produces:
#   .../GCS_exports/PARQUET/gee_hourly_ds/yyyymm=YYYYMM/part.parquet
#   .../outdoor_with_gee_hourly_UTC_ds/yyyymm=YYYYMM/part.parquet
#   .../outdoor_joined_month_rds/outdoor_joined_YYYYMM.rds
# ============================================================
# ============================================================
# CLEAN OLD OUTPUTS BEFORE REBUILD
# ============================================================

base_dir <- "/Users/priyanka/Downloads/PurpleAir_EPA/GCS_exports"

dirs_to_remove <- c(
  file.path(base_dir, "PARQUET", "gee_hourly_ds"),
  "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_with_gee_hourly_UTC_ds",
  "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_joined_month_rds"
)

files_to_remove <- c(
  "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_all_joined.rds"
)

say <- function(...) cat(sprintf(...), "\n")

say("Removing old directories...")
for (p in dirs_to_remove) {
  if (dir.exists(p)) {
    say("  deleting dir: %s", p)
    unlink(p, recursive = TRUE, force = TRUE)
  } else {
    say("  dir not found, skipping: %s", p)
  }
}

say("Removing old files...")
for (p in files_to_remove) {
  if (file.exists(p)) {
    say("  deleting file: %s", p)
    unlink(p, force = TRUE)
  } else {
    say("  file not found, skipping: %s", p)
  }
}

# recreate empty output dirs so the next chunk writes cleanly
dir.create(file.path(base_dir, "PARQUET", "gee_hourly_ds"), recursive = TRUE, showWarnings = FALSE)
dir.create("/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_with_gee_hourly_UTC_ds", recursive = TRUE, showWarnings = FALSE)
dir.create("/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_joined_month_rds", recursive = TRUE, showWarnings = FALSE)

say("Cleanup complete.")

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
  library(data.table)
  library(stringr)
  library(lubridate)
  library(arrow)
})

# ----------------------------
# PATHS
# ----------------------------
base_dir     <- "/Users/priyanka/Downloads/PurpleAir_EPA/GCS_exports"
era_csv_dir  <- file.path(base_dir, "ERA")
pm_csv_dir   <- file.path(base_dir, "PM")

parq_root    <- file.path(base_dir, "PARQUET")
era_parq_dir <- file.path(parq_root, "ERA_PARQUET")
pm_parq_dir  <- file.path(parq_root, "PM_PARQUET")
gee_hourly_parq <- file.path(parq_root, "gee_hourly_ds")

outdoor_rdata <- "/Users/priyanka/Downloads/PurpleAir_EPA/annual_data_all_outdoor.RData"
out_dir   <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_with_gee_hourly_UTC_ds"
out_store <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_joined_month_rds"

tmp_dir <- "/Users/priyanka/Downloads/PurpleAir_EPA/duckdb_tmp"
db_path_merge <- file.path(base_dir, "gee_merge.duckdb")
db_path_join  <- file.path(base_dir, "join_outdoor_gee.duckdb")

dir.create(era_parq_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pm_parq_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(gee_hourly_parq, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(out_store, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(dir.exists(era_csv_dir), dir.exists(pm_csv_dir))
stopifnot(file.exists(outdoor_rdata))

# ----------------------------
# UTIL
# ----------------------------
quiet_db <- function(con, sql) invisible(DBI::dbExecute(con, sql))
say <- function(...) cat(sprintf(...), "\n")

yyyymm_safe <- function(d) {
  d <- as.Date(d)
  strftime(as.POSIXct(d, tz = "UTC"), "%Y%m", tz = "UTC")
}

parq_exists <- function(p) {
  file.exists(p) && is.finite(file.info(p)$size) && file.info(p)$size > 0
}

fail_log <- file.path(parq_root, "pipeline_failures.txt")
if (file.exists(fail_log)) file.remove(fail_log)

log_fail <- function(step, f, msg) {
  cat(sprintf("[%s]\t%s\t%s\n", step, f, msg), file = fail_log, append = TRUE)
}

sniff_delim <- function(f) {
  hdr <- readLines(f, n = 1, warn = FALSE)
  if (!length(hdr)) return(",")
  n_tab <- length(strsplit(hdr, "\t", fixed = TRUE)[[1]])
  n_com <- length(strsplit(hdr, ",",  fixed = TRUE)[[1]])
  if (n_tab >= n_com) "\t" else ","
}

sniff_cols <- function(f, delim) {
  hdr <- readLines(f, n = 1, warn = FALSE)
  cols <- strsplit(hdr, delim, fixed = TRUE)[[1]]
  trimws(gsub('^"|"$', "", cols))
}

sniff_dt_format <- function(f, delim) {
  x <- readLines(f, n = 6, warn = FALSE)
  if (length(x) < 2) return("%Y-%m-%d %H:%M:%S")

  cols <- trimws(gsub('^"|"$', "", strsplit(x[1], delim, fixed = TRUE)[[1]]))
  j <- match("datetime", cols)
  if (is.na(j)) return("%Y-%m-%d %H:%M:%S")

  row <- strsplit(x[2], delim, fixed = TRUE)[[1]]
  dt  <- gsub('^"|"$', "", row[j])

  if (grepl("T", dt) && grepl("Z$", dt) && grepl("\\.\\d+", dt)) return("%Y-%m-%dT%H:%M:%S.%fZ")
  if (grepl("T", dt) && grepl("Z$", dt)) return("%Y-%m-%dT%H:%M:%SZ")
  if (grepl("T", dt)) return("%Y-%m-%dT%H:%M:%S")
  "%Y-%m-%d %H:%M:%S"
}

pick_latest_exports <- function(files, prefix_regex) {
  dt <- data.table(path = files, file = basename(files))

  mA <- str_match(
    dt$file,
    paste0("^", prefix_regex, "_(\\d{6})_batch(\\d{2})_(\\d{4}_\\d{2}_\\d{2}_\\d{2}_\\d{2}_\\d{2})\\.csv$")
  )
  okA <- !is.na(mA[, 1])

  mB <- str_match(
    dt$file,
    paste0("^", prefix_regex, "_(\\d{6})_batch(\\d{2})\\.csv$")
  )
  okB <- !is.na(mB[, 1])

  if (!any(okA | okB)) stop("No files matched expected pattern for ", prefix_regex)

  out <- rbindlist(list(
    dt[okA, .(
      path,
      yyyymm = mA[okA, 2],
      batch  = mA[okA, 3],
      ts     = as.POSIXct(mA[okA, 4], format = "%Y_%m_%d_%H_%M_%S", tz = "UTC")
    )],
    dt[okB, .(
      path,
      yyyymm = mB[okB, 2],
      batch  = mB[okB, 3],
      ts     = as.POSIXct(file.info(path)$mtime, tz = "UTC")
    )]
  ), use.names = TRUE, fill = TRUE)

  setorder(out, yyyymm, batch, ts)
  latest <- out[, .SD[.N], by = .(yyyymm, batch)]
  latest$path
}

file_yyyymm <- function(f, prefix_regex) {
  str_match(basename(f), paste0("^", prefix_regex, "_(\\d{6})_batch\\d{2}"))[, 2]
}

list_months <- function(root_dir) {
  d <- list.dirs(root_dir, full.names = FALSE, recursive = FALSE)
  d <- d[grepl("^yyyymm=", d)]
  sort(sub("^yyyymm=", "", d))
}

glob_parq <- function(root, yyyymm) file.path(root, paste0("yyyymm=", yyyymm), "*.parquet")

# ----------------------------
# 1) OPEN DUCKDB FOR CSV -> PARQUET
# ----------------------------
say("Opening DuckDB for CSV -> parquet conversion...")
con <- dbConnect(duckdb::duckdb(), db_path_merge)
quiet_db(con, "PRAGMA threads=16;")
quiet_db(con, "PRAGMA memory_limit='80GB';")
quiet_db(con, "PRAGMA preserve_insertion_order=false;")
quiet_db(con, "PRAGMA enable_progress_bar=false;")
quiet_db(con, sprintf("PRAGMA temp_directory='%s';", gsub("'", "''", tmp_dir)))

# ----------------------------
# 2) DEDUPE FILES
# ----------------------------
era_all <- list.files(era_csv_dir, pattern = "\\.csv$", full.names = TRUE)
pm_all  <- list.files(pm_csv_dir,  pattern = "\\.csv$", full.names = TRUE)
stopifnot(length(era_all) > 0, length(pm_all) > 0)

era_files <- pick_latest_exports(era_all, "ERA5")
pm_files  <- pick_latest_exports(pm_all,  "GEOSCF_PM25")

say("ERA raw files: %d | kept latest: %d", length(era_all), length(era_files))
say("PM  raw files: %d | kept latest: %d", length(pm_all),  length(pm_files))

# ----------------------------
# 3) BUILD SELECT LISTS
# ----------------------------
sel_for_era <- function(cols, dt_fmt) {
  has <- function(x) x %in% cols

  pid     <- if (has("pid")) "try_cast(pid AS INTEGER) AS pid" else "NULL::INTEGER AS pid"
  station <- if (has("station")) "try_cast(station AS BIGINT) AS station" else "NULL::BIGINT AS station"
  lon     <- if (has("lon")) "try_cast(lon AS DOUBLE) AS lon" else "NULL::DOUBLE AS lon"
  lat     <- if (has("lat")) "try_cast(lat AS DOUBLE) AS lat" else "NULL::DOUBLE AS lat"

  t2m <- if (has("temperature_2m")) "try_cast(temperature_2m AS DOUBLE) AS temperature_2m" else "NULL::DOUBLE AS temperature_2m"
  sp  <- if (has("surface_pressure")) "try_cast(surface_pressure AS DOUBLE) AS surface_pressure" else "NULL::DOUBLE AS surface_pressure"
  d2m <- if (has("dewpoint_temperature_2m")) "try_cast(dewpoint_temperature_2m AS DOUBLE) AS dewpoint_temperature_2m" else "NULL::DOUBLE AS dewpoint_temperature_2m"

  paste(
    pid,
    station,
    sprintf("date_trunc('hour', strptime(datetime, '%s')) AS datetime_hr", dt_fmt),
    t2m, sp, d2m,
    lon, lat,
    sep = ",\n        "
  )
}

sel_for_pm <- function(cols, dt_fmt) {
  has <- function(x) x %in% cols

  pid     <- if (has("pid")) "try_cast(pid AS INTEGER) AS pid" else "NULL::INTEGER AS pid"
  station <- if (has("station")) "try_cast(station AS BIGINT) AS station" else "NULL::BIGINT AS station"
  lon     <- if (has("lon")) "try_cast(lon AS DOUBLE) AS lon" else "NULL::DOUBLE AS lon"
  lat     <- if (has("lat")) "try_cast(lat AS DOUBLE) AS lat" else "NULL::DOUBLE AS lat"

  pm <- if (has("PM25_RH35_GCC")) {
    "try_cast(PM25_RH35_GCC AS DOUBLE) AS geoscf_pm25"
  } else if (has("first")) {
    "try_cast(first AS DOUBLE) AS geoscf_pm25"
  } else {
    "NULL::DOUBLE AS geoscf_pm25"
  }

  paste(
    pid,
    station,
    sprintf("date_trunc('hour', strptime(datetime, '%s')) AS datetime_hr", dt_fmt),
    pm,
    lon, lat,
    sep = ",\n        "
  )
}

# ----------------------------
# 4) CONVERT CSV -> PARQUET
# ----------------------------
convert_one <- function(f, kind = c("ERA", "PM")) {
  kind <- match.arg(kind)

  if (!file.exists(f) || file.info(f)$size < 500) {
    stop("CSV too small/empty: ", f)
  }

  prefix <- if (kind == "ERA") "ERA5" else "GEOSCF_PM25"
  yyyymm <- file_yyyymm(f, prefix)
  if (is.na(yyyymm)) stop("Cannot parse yyyymm from: ", basename(f))

  out_base <- if (kind == "ERA") era_parq_dir else pm_parq_dir
  out_part <- file.path(out_base, paste0("yyyymm=", yyyymm))
  dir.create(out_part, recursive = TRUE, showWarnings = FALSE)

  out <- file.path(out_part, paste0(basename(f), ".parquet"))
  if (parq_exists(out)) return(TRUE)

  delim  <- sniff_delim(f)
  cols   <- sniff_cols(f, delim)
  dt_fmt <- sniff_dt_format(f, delim)

  if (!all(c("pid", "station", "datetime") %in% cols)) {
    stop("Missing required columns pid/station/datetime. Found: ", paste(cols, collapse = ", "))
  }

  sel <- if (kind == "ERA") sel_for_era(cols, dt_fmt) else sel_for_pm(cols, dt_fmt)

  say("Converting %s | %s | delim='%s' | has pid=%s station=%s lon=%s lat=%s",
      kind, basename(f), delim,
      "pid" %in% cols, "station" %in% cols, "lon" %in% cols, "lat" %in% cols)

  sql <- sprintf("
    COPY (
      SELECT
        %s
      FROM read_csv('%s',
        delim='%s',
        header=true,
        strict_mode=false,
        ignore_errors=true,
        null_padding=true,
        all_varchar=true,
        sample_size=-1
      )
      WHERE pid IS NOT NULL
        AND station IS NOT NULL
        AND datetime_hr IS NOT NULL
    ) TO '%s' (FORMAT PARQUET, COMPRESSION ZSTD);
  ", sel, gsub("'", "''", f), delim, gsub("'", "''", out))

  dbExecute(con, sql)
  TRUE
}

say("Converting ERA -> parquet...")
for (i in seq_along(era_files)) {
  if (i %% 25 == 0 || i == length(era_files)) say("ERA progress: %d / %d", i, length(era_files))
  f <- era_files[i]
  ok <- tryCatch(convert_one(f, "ERA"),
                 error = function(e) {
                   log_fail("ERA_CONVERT", basename(f), conditionMessage(e))
                   FALSE
                 })
  if (!ok) say("FAILED ERA: %s", basename(f))
}

say("Converting PM -> parquet...")
for (i in seq_along(pm_files)) {
  if (i %% 25 == 0 || i == length(pm_files)) say("PM progress: %d / %d", i, length(pm_files))
  f <- pm_files[i]
  ok <- tryCatch(convert_one(f, "PM"),
                 error = function(e) {
                   log_fail("PM_CONVERT", basename(f), conditionMessage(e))
                   FALSE
                 })
  if (!ok) say("FAILED PM: %s", basename(f))
}

# ----------------------------
# 4B) VALIDATE PARQUET OUTPUTS
# ----------------------------
is_valid_parquet <- function(f) {
  tryCatch({
    arrow::read_parquet(f, as_data_frame = FALSE)
    TRUE
  }, error = function(e) FALSE)
}

era_parqs <- list.files(era_parq_dir, pattern = "\\.parquet$", full.names = TRUE, recursive = TRUE)
pm_parqs  <- list.files(pm_parq_dir,  pattern = "\\.parquet$", full.names = TRUE, recursive = TRUE)

era_ok <- vapply(era_parqs, is_valid_parquet, logical(1))
pm_ok  <- vapply(pm_parqs,  is_valid_parquet, logical(1))

if (any(!era_ok)) file.remove(era_parqs[!era_ok])
if (any(!pm_ok))  file.remove(pm_parqs[!pm_ok])

era_parqs_ok <- era_parqs[era_ok]
pm_parqs_ok  <- pm_parqs[pm_ok]

say("Validated ERA parquet files: %d", length(era_parqs_ok))
say("Validated PM  parquet files: %d", length(pm_parqs_ok))
stopifnot(length(era_parqs_ok) > 0, length(pm_parqs_ok) > 0)

# ----------------------------
# 5) MERGE ERA + PM -> GEE_HOURLY_DS
#    CRITICAL FIX: MERGE ON pid + datetime_hr
# ----------------------------
say("Merging monthly ERA + PM parquet into gee_hourly_ds by pid + datetime_hr ...")

if (exists("con", inherits = TRUE)) {
  try(DBI::dbDisconnect(con, shutdown = FALSE), silent = TRUE)
  rm(con)
}
gc()

con <- DBI::dbConnect(duckdb::duckdb(), db_path_merge)
quiet_db(con, sprintf("PRAGMA temp_directory='%s';", gsub("'", "''", tmp_dir)))
quiet_db(con, "PRAGMA max_temp_directory_size='900GiB';")
quiet_db(con, "PRAGMA preserve_insertion_order=false;")
quiet_db(con, "PRAGMA enable_progress_bar=false;")
quiet_db(con, "PRAGMA threads=4;")
quiet_db(con, "PRAGMA memory_limit='90GB';")

era_months <- list_months(era_parq_dir)
pm_months  <- list_months(pm_parq_dir)
months_all <- sort(unique(c(era_months, pm_months)))
stopifnot(length(months_all) > 0)

say("Months found for merge: %d", length(months_all))
say("First few months: %s", paste(head(months_all, 6), collapse = ", "))

for (m in months_all) {

  out_m_dir  <- file.path(gee_hourly_parq, paste0("yyyymm=", m))
  out_m_file <- file.path(out_m_dir, "part.parquet")
  dir.create(out_m_dir, recursive = TRUE, showWarnings = FALSE)

  if (file.exists(out_m_file) && file.info(out_m_file)$size > 0) {
    say("Skipping existing gee month: %s", m)
    next
  }

  era_glob <- glob_parq(era_parq_dir, m)
  pm_glob  <- glob_parq(pm_parq_dir,  m)

  era_files_m <- Sys.glob(era_glob)
  pm_files_m  <- Sys.glob(pm_glob)

  era_has <- length(era_files_m) > 0
  pm_has  <- length(pm_files_m)  > 0
  if (!era_has && !pm_has) next

  say("Merging month %s | ERA files: %d | PM files: %d", m, length(era_files_m), length(pm_files_m))

  quiet_db(con, "DROP VIEW IF EXISTS era_m;")
  quiet_db(con, "DROP VIEW IF EXISTS pm_m;")
  quiet_db(con, "DROP TABLE IF EXISTS gee_m;")

  if (era_has) {
    quiet_db(con, sprintf(
      "CREATE VIEW era_m AS SELECT * FROM read_parquet('%s', union_by_name=true);",
      gsub("'", "''", era_glob)
    ))
  }
  if (pm_has) {
    quiet_db(con, sprintf(
      "CREATE VIEW pm_m AS SELECT * FROM read_parquet('%s', union_by_name=true);",
      gsub("'", "''", pm_glob)
    ))
  }

  if (era_has && pm_has) {
    quiet_db(con, "
      CREATE TABLE gee_m AS
      SELECT
        COALESCE(e.pid, p.pid) AS gee_pid,
        e.pid AS era_pid,
        p.pid AS pm_pid,
        COALESCE(e.station, p.station) AS station,
        COALESCE(e.datetime_hr, p.datetime_hr) AS datetime_hr,
        e.temperature_2m,
        e.surface_pressure,
        e.dewpoint_temperature_2m,
        p.geoscf_pm25,
        COALESCE(p.lon, e.lon) AS lon,
        COALESCE(p.lat, e.lat) AS lat
      FROM era_m e
      FULL JOIN pm_m p
        ON e.pid = p.pid
       AND e.datetime_hr = p.datetime_hr;
    ")
  } else if (era_has) {
    quiet_db(con, "
      CREATE TABLE gee_m AS
      SELECT
        pid AS gee_pid,
        pid AS era_pid,
        NULL::INTEGER AS pm_pid,
        station,
        datetime_hr,
        temperature_2m,
        surface_pressure,
        dewpoint_temperature_2m,
        NULL::DOUBLE AS geoscf_pm25,
        lon, lat
      FROM era_m;
    ")
  } else {
    quiet_db(con, "
      CREATE TABLE gee_m AS
      SELECT
        pid AS gee_pid,
        NULL::INTEGER AS era_pid,
        pid AS pm_pid,
        station,
        datetime_hr,
        NULL::DOUBLE AS temperature_2m,
        NULL::DOUBLE AS surface_pressure,
        NULL::DOUBLE AS dewpoint_temperature_2m,
        geoscf_pm25,
        lon, lat
      FROM pm_m;
    ")
  }

  diag_merge <- dbGetQuery(con, "
    SELECT
      COUNT(*) AS n_rows,
      COUNT(DISTINCT CAST(gee_pid AS VARCHAR) || '|' || CAST(datetime_hr AS VARCHAR)) AS n_pid_keys,
      COUNT(DISTINCT CAST(station AS VARCHAR) || '|' || CAST(datetime_hr AS VARCHAR)) AS n_station_keys,
      SUM(CASE WHEN era_pid IS NOT NULL THEN 1 ELSE 0 END) AS n_with_era_pid,
      SUM(CASE WHEN pm_pid  IS NOT NULL THEN 1 ELSE 0 END) AS n_with_pm_pid,
      SUM(CASE WHEN gee_pid IS NOT NULL THEN 1 ELSE 0 END) AS n_with_gee_pid,
      SUM(CASE WHEN lon IS NOT NULL AND lat IS NOT NULL THEN 1 ELSE 0 END) AS n_with_coords
    FROM gee_m;
  ")
  print(diag_merge)

  quiet_db(con, sprintf(
    "COPY gee_m TO '%s' (FORMAT PARQUET, COMPRESSION ZSTD);",
    gsub("'", "''", out_m_file)
  ))

  quiet_db(con, "DROP TABLE IF EXISTS gee_m;")
  quiet_db(con, "DROP VIEW IF EXISTS era_m;")
  quiet_db(con, "DROP VIEW IF EXISTS pm_m;")

  say("Wrote gee_hourly month: %s -> %s", m, out_m_file)
}

say("Finished gee_hourly dataset: %s", gee_hourly_parq)

DBI::dbDisconnect(con, shutdown = TRUE)
rm(con)
gc()

# ----------------------------
# 6) JOIN OUTDOOR TO GEE MONTH-BY-MONTH
#    OUTDOOR JOINS ON station + datetime_hr
#    GEE is first collapsed from pid-level to one row per station + datetime_hr
#    IMPORTANT: DO NOT dbWriteTable() the full outdoor object
# ----------------------------
say("Joining outdoor -> gee_hourly_ds month-by-month...")

conj <- DBI::dbConnect(duckdb::duckdb(), db_path_join)
quiet_db(conj, "PRAGMA preserve_insertion_order=false;")
quiet_db(conj, "PRAGMA enable_object_cache=true;")
quiet_db(conj, sprintf("PRAGMA temp_directory='%s';", gsub("'", "''", tmp_dir)))
quiet_db(conj, "PRAGMA threads=4;")
quiet_db(conj, "PRAGMA memory_limit='55GB';")
quiet_db(conj, "PRAGMA enable_progress_bar=false;")

obj_out <- load(outdoor_rdata)
cand_out <- obj_out[vapply(obj_out, function(nm) inherits(get(nm), c("data.frame", "data.table")), logical(1))]
stopifnot(length(cand_out) >= 1)
outdoor <- as.data.table(get(cand_out[1]))

if ("STATION(-)" %in% names(outdoor))     setnames(outdoor, "STATION(-)", "station")
if ("Timestamp(UTC)" %in% names(outdoor)) setnames(outdoor, "Timestamp(UTC)", "datetime_utc")
stopifnot(all(c("station", "datetime_utc") %in% names(outdoor)))

outdoor[, station := as.character(station)]
if (!inherits(outdoor$datetime_utc, "POSIXct")) {
  outdoor[, datetime_utc := as.POSIXct(datetime_utc, tz = "UTC")]
} else {
  attr(outdoor$datetime_utc, "tzone") <- "UTC"
}
outdoor[, datetime_hr := lubridate::floor_date(datetime_utc, "hour")]

say("Outdoor rows loaded: %s", format(nrow(outdoor), big.mark = ","))
say("Outdoor stations: %s", format(uniqueN(outdoor$station), big.mark = ","))

# modest indexing in data.table helps monthly filtering
data.table::setindex(outdoor, datetime_hr)
data.table::setindex(outdoor, station)

start_month  <- as.Date("2018-01-01")
end_month    <- as.Date("2025-10-01")
month_starts <- seq.Date(start_month, end_month, by = "month")
next_month   <- function(d) seq.Date(as.Date(d), by = "month", length.out = 2)[2]

gee_glob_all <- file.path(gee_hourly_parq, "**", "*.parquet")

for (m_start_raw in month_starts) {

  m_start <- as.Date(m_start_raw)
  m_end   <- next_month(m_start)
  yyyymm  <- yyyymm_safe(m_start)

  out_month_dir <- file.path(out_dir, paste0("yyyymm=", yyyymm))
  dir.create(out_month_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_month_dir, "part.parquet")

  if (file.exists(out_file) && file.info(out_file)$size > 0) {
    say("Skipping existing joined month: %s", yyyymm)
    next
  }

  m_start_ts <- strftime(as.POSIXct(m_start, tz = "UTC"), "%Y-%m-%d %H:%M:%S", tz = "UTC")
  m_end_ts   <- strftime(as.POSIXct(m_end,   tz = "UTC"), "%Y-%m-%d %H:%M:%S", tz = "UTC")

  say("Joining month %s ...", yyyymm)

  # ----- FILTER OUTDOOR IN R FIRST, THEN WRITE ONLY THAT MONTH -----
  outdoor_m_r <- outdoor[
    datetime_hr >= as.POSIXct(m_start, tz = "UTC") &
    datetime_hr <  as.POSIXct(m_end,   tz = "UTC") &
    !is.na(station) &
    !is.na(datetime_hr)
  ]

  n_m <- nrow(outdoor_m_r)
  say("  outdoor rows this month: %s", format(n_m, big.mark = ","))

  if (is.na(n_m) || n_m == 0) {
    rm(outdoor_m_r)
    gc()
    next
  }

  quiet_db(conj, "DROP TABLE IF EXISTS outdoor_m;")
  DBI::dbWriteTable(conj, "outdoor_m", outdoor_m_r, overwrite = TRUE)
  rm(outdoor_m_r)
  gc()

  quiet_db(conj, "CREATE INDEX IF NOT EXISTS idx_outdoor_m_key ON outdoor_m(station, datetime_hr);")

  quiet_db(conj, "DROP TABLE IF EXISTS keys_m;")
  quiet_db(conj, "
    CREATE TABLE keys_m AS
    SELECT DISTINCT station, datetime_hr
    FROM outdoor_m;
  ")

  quiet_db(conj, "DROP TABLE IF EXISTS station_typ_loc_m;")
  quiet_db(conj, sprintf("
    CREATE TABLE station_typ_loc_m AS
    SELECT
      CAST(station AS VARCHAR) AS station,
      median(lon) AS typ_lon,
      median(lat) AS typ_lat
    FROM read_parquet('%s')
    WHERE datetime_hr >= TIMESTAMP '%s'
      AND datetime_hr <  TIMESTAMP '%s'
      AND station IS NOT NULL
      AND lon IS NOT NULL
      AND lat IS NOT NULL
    GROUP BY 1;
  ", gsub("'", "''", gee_glob_all), m_start_ts, m_end_ts))

  quiet_db(conj, "DROP TABLE IF EXISTS gee_dedup_m;")
  quiet_db(conj, sprintf("
    CREATE TABLE gee_dedup_m AS
    WITH gee_match AS (
      SELECT
        CAST(g.station AS VARCHAR) AS station,
        g.datetime_hr,
        CAST(g.era_pid AS INTEGER) AS era_pid,
        CAST(g.pm_pid  AS INTEGER) AS pm_pid,
        CAST(g.gee_pid AS INTEGER) AS gee_pid,
        g.temperature_2m,
        g.surface_pressure,
        g.dewpoint_temperature_2m,
        g.geoscf_pm25,
        g.lon,
        g.lat,
        (
          CASE WHEN g.temperature_2m IS NOT NULL THEN 1 ELSE 0 END +
          CASE WHEN g.surface_pressure IS NOT NULL THEN 1 ELSE 0 END +
          CASE WHEN g.dewpoint_temperature_2m IS NOT NULL THEN 1 ELSE 0 END +
          CASE WHEN g.geoscf_pm25 IS NOT NULL THEN 1 ELSE 0 END
        ) AS nn_score
      FROM read_parquet('%s') g
      INNER JOIN keys_m k
        ON CAST(g.station AS VARCHAR) = k.station
       AND g.datetime_hr = k.datetime_hr
      WHERE g.datetime_hr >= TIMESTAMP '%s'
        AND g.datetime_hr <  TIMESTAMP '%s'
    ),
    gee_scored AS (
      SELECT
        gm.*,
        tl.typ_lon,
        tl.typ_lat,
        CASE
          WHEN gm.lon IS NULL OR gm.lat IS NULL OR tl.typ_lon IS NULL OR tl.typ_lat IS NULL
          THEN NULL
          ELSE ((gm.lon - tl.typ_lon)*(gm.lon - tl.typ_lon) + (gm.lat - tl.typ_lat)*(gm.lat - tl.typ_lat))
        END AS dist2
      FROM gee_match gm
      LEFT JOIN station_typ_loc_m tl
        ON gm.station = tl.station
    ),
    ranked AS (
      SELECT
        *,
        row_number() OVER (
          PARTITION BY station, datetime_hr
          ORDER BY
            nn_score DESC,
            CASE WHEN dist2 IS NULL THEN 1 ELSE 0 END ASC,
            dist2 ASC,
            gee_pid ASC
        ) AS rn
      FROM gee_scored
    )
    SELECT
      station, datetime_hr,
      era_pid, pm_pid, gee_pid,
      temperature_2m, surface_pressure, dewpoint_temperature_2m, geoscf_pm25,
      lon, lat
    FROM ranked
    WHERE rn = 1;
  ", gsub("'", "''", gee_glob_all), m_start_ts, m_end_ts))

  diag_gee <- DBI::dbGetQuery(conj, "
    SELECT
      COUNT(*) AS n_rows,
      COUNT(DISTINCT station || '|' || CAST(datetime_hr AS VARCHAR)) AS n_station_keys,
      SUM(CASE WHEN era_pid IS NOT NULL THEN 1 ELSE 0 END) AS n_with_era_pid,
      SUM(CASE WHEN pm_pid  IS NOT NULL THEN 1 ELSE 0 END) AS n_with_pm_pid,
      SUM(CASE WHEN gee_pid IS NOT NULL THEN 1 ELSE 0 END) AS n_with_gee_pid,
      SUM(CASE WHEN lon IS NOT NULL AND lat IS NOT NULL THEN 1 ELSE 0 END) AS n_with_coords
    FROM gee_dedup_m;
  ")
  print(diag_gee)

  quiet_db(conj, "DROP TABLE IF EXISTS joined_m;")
  quiet_db(conj, "
    CREATE TABLE joined_m AS
    SELECT
      o.*,
      g.era_pid,
      g.pm_pid,
      g.gee_pid,
      g.temperature_2m,
      g.surface_pressure,
      g.dewpoint_temperature_2m,
      g.geoscf_pm25,
      g.lon AS gee_longitude,
      g.lat AS gee_latitude
    FROM outdoor_m o
    LEFT JOIN gee_dedup_m g
      ON o.station = g.station
     AND o.datetime_hr = g.datetime_hr;
  ")

  diag_join <- DBI::dbGetQuery(conj, "
    SELECT
      COUNT(*) AS n_rows,
      SUM(CASE WHEN gee_pid IS NOT NULL THEN 1 ELSE 0 END) AS n_matched_gee,
      SUM(CASE WHEN era_pid IS NOT NULL THEN 1 ELSE 0 END) AS n_with_era_pid,
      SUM(CASE WHEN pm_pid  IS NOT NULL THEN 1 ELSE 0 END) AS n_with_pm_pid,
      SUM(CASE WHEN gee_longitude IS NOT NULL AND gee_latitude IS NOT NULL THEN 1 ELSE 0 END) AS n_with_gee_coords
    FROM joined_m;
  ")
  print(diag_join)

  quiet_db(conj, sprintf(
    "COPY joined_m TO '%s' (FORMAT PARQUET, COMPRESSION ZSTD);",
    gsub("'", "''", out_file)
  ))

  quiet_db(conj, "DROP TABLE IF EXISTS joined_m;")
  quiet_db(conj, "DROP TABLE IF EXISTS gee_dedup_m;")
  quiet_db(conj, "DROP TABLE IF EXISTS station_typ_loc_m;")
  quiet_db(conj, "DROP TABLE IF EXISTS keys_m;")
  quiet_db(conj, "DROP TABLE IF EXISTS outdoor_m;")

  say("  wrote joined parquet: %s", out_file)
}

rm(outdoor)
gc()

# ----------------------------
# 7) OPTIONAL MONTHLY RDS
# ----------------------------
say("Building optional monthly RDS files...")

month_dirs <- list.dirs(out_dir, full.names = TRUE, recursive = FALSE)
month_dirs <- month_dirs[grepl("^yyyymm=", basename(month_dirs))]
month_dirs <- month_dirs[order(sub("^yyyymm=", "", basename(month_dirs)))]

for (md in month_dirs) {
  m <- sub("^yyyymm=", "", basename(md))
  out_rds <- file.path(out_store, paste0("outdoor_joined_", m, ".rds"))

  if (file.exists(out_rds) && file.info(out_rds)$size > 0) {
    say("Skipping existing RDS month: %s", m)
    next
  }

  say("Creating RDS for month: %s", m)

  ds_m <- arrow::open_dataset(md, format = "parquet")
  dt_m <- as.data.table(as.data.frame(ds_m))

  say("  rows in RDS month %s: %s", m, format(nrow(dt_m), big.mark = ","))
  saveRDS(dt_m, out_rds, compress = "xz")

  rm(ds_m, dt_m)
  gc()
}

# ----------------------------
# 8) FINAL DIAGNOSTICS
# ----------------------------
say("DONE.")
say("gee_hourly dataset dir: %s", gee_hourly_parq)
say("joined monthly parquet dir: %s", out_dir)
say("monthly RDS dir: %s", out_store)

n_gee_parts  <- length(list.files(gee_hourly_parq, pattern = "\\.parquet$", recursive = TRUE))
n_join_parts <- length(list.files(out_dir, pattern = "\\.parquet$", recursive = TRUE))
n_rds        <- length(list.files(out_store, pattern = "\\.rds$", full.names = TRUE))

say("gee_hourly parquet parts: %d", n_gee_parts)
say("joined parquet parts: %d", n_join_parts)
say("joined monthly RDS files: %d", n_rds)

if (file.exists(fail_log) && file.info(fail_log)$size > 0) {
  say("Some steps failed. See log: %s", fail_log)
  say("First few failure lines:")
  print(readLines(fail_log, n = 10))
} else {
  say("No failures logged.")
}
