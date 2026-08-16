# ==========================================================================
# 05_combine_joined.R
#
# Combine the monthly joined RDS files into the single big joined object
# (outdoor_all_joined.rds).
#
# Extracted VERBATIM from purpleair_degradation2.Rmd (chunk #19 lines 4652-4922).
# Code is unmodified; see REPRODUCING.md for paths and run order.
# ==========================================================================

# ---- source chunk #19 (Rmd lines 4652-4922) ------------------------

# ============================================================
# OPTION B: IF YOU REALLY NEED ONE FILE
# CLEAN + COMBINE MONTHLY JOINED RDS FILES INTO ONE BIG RDS
#
# What this does:
# 1) finds monthly joined RDS files
# 2) drops zero-byte files
# 3) inspects each month and deletes empty/corrupt monthly RDS files
# 4) optionally deletes existing combined RDS
# 5) combines valid monthly RDS files
# 6) prints diagnostics
#
# IMPORTANT:
# - This does NOT rebuild monthly files.
# - It only combines what already exists in outdoor_joined_month_rds.
# - If your monthly files are old/stale and missing PID columns,
#   you must rerun the monthly join pipeline first.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

in_dir  <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_joined_month_rds"
out_rds <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_all_joined.rds"

say <- function(...) cat(sprintf(...), "\n")

stopifnot(dir.exists(in_dir))

# ------------------------------------------------------------
# 0) OPTIONAL CLEANUP
# ------------------------------------------------------------
delete_existing_combined_rds <- TRUE
delete_bad_monthly_rds       <- TRUE

if (delete_existing_combined_rds && file.exists(out_rds)) {
  say("Deleting existing combined RDS: %s", out_rds)
  unlink(out_rds, force = TRUE)
}

# ------------------------------------------------------------
# 1) FIND MONTHLY RDS FILES
# ------------------------------------------------------------
files <- list.files(
  in_dir,
  pattern = "^outdoor_joined_[0-9]{6}\\.rds$",
  full.names = TRUE
)
stopifnot(length(files) > 0)

file_info <- data.table(
  path = files,
  file = basename(files),
  size_bytes = file.info(files)$size
)

file_info[, yyyymm := strcapture(
  "outdoor_joined_([0-9]{6})\\.rds",
  file,
  proto = list(yyyymm = character())
)$yyyymm]

# keep only non-empty files on disk
file_info <- file_info[
  is.finite(size_bytes) & size_bytes > 0 & !is.na(yyyymm)
]

if (nrow(file_info) == 0) {
  stop("No non-empty monthly RDS files found in: ", in_dir)
}

# if duplicate months exist, keep the largest file
setorder(file_info, yyyymm, size_bytes)
file_info <- file_info[, .SD[.N], by = yyyymm]
setorder(file_info, yyyymm)

say("Monthly RDS candidates found: %d", nrow(file_info))
say("Month range: %s -> %s", file_info$yyyymm[1], file_info$yyyymm[nrow(file_info)])

# ------------------------------------------------------------
# 2) INSPECT EACH MONTH AND DROP BAD MONTHS
# ------------------------------------------------------------
say("Inspecting monthly RDS files and removing bad ones if requested...")

inspect_diag <- vector("list", nrow(file_info))
keep_idx <- rep(TRUE, nrow(file_info))

for (i in seq_len(nrow(file_info))) {
  f <- file_info$path[i]
  m <- file_info$yyyymm[i]

  say("Inspecting %s (%d / %d)", basename(f), i, nrow(file_info))

  ok <- TRUE
  dt <- tryCatch(readRDS(f), error = function(e) e)

  if (inherits(dt, "error")) {
    say("  ERROR reading %s: %s", basename(f), conditionMessage(dt))
    ok <- FALSE
    n_rows <- NA_integer_
    n_cols <- NA_integer_
    nm <- character()
  } else {
    setDT(dt)
    n_rows <- nrow(dt)
    n_cols <- ncol(dt)
    nm <- names(dt)

    say("  rows: %s | cols: %d", format(n_rows, big.mark = ","), n_cols)

    if (n_rows == 0L || n_cols == 0L) {
      say("  -> bad monthly file (empty table)")
      ok <- FALSE
    }
  }

  inspect_diag[[i]] <- data.table(
    yyyymm = m,
    file = basename(f),
    n_rows = n_rows,
    n_cols = n_cols,
    ok = ok,
    has_era_pid = "era_pid" %in% nm,
    has_pm_pid  = "pm_pid"  %in% nm,
    has_gee_pid = "gee_pid" %in% nm,
    has_gee_lon = "gee_longitude" %in% nm,
    has_gee_lat = "gee_latitude"  %in% nm
  )

  keep_idx[i] <- ok

  if (!ok && delete_bad_monthly_rds && file.exists(f)) {
    say("  deleting bad monthly RDS: %s", f)
    unlink(f, force = TRUE)
  }

  if (!inherits(dt, "error")) {
    rm(dt)
    gc()
  }
}

inspect_diag <- rbindlist(inspect_diag, use.names = TRUE, fill = TRUE)
say("Inspection summary:")
print(inspect_diag)

file_info <- file_info[keep_idx]

if (nrow(file_info) == 0) {
  stop("No valid monthly RDS files remain after inspection.")
}

say("Valid monthly RDS files remaining: %d", nrow(file_info))
say("Valid month range: %s -> %s", file_info$yyyymm[1], file_info$yyyymm[nrow(file_info)])

# ------------------------------------------------------------
# 3) READ FIRST VALID FILE TO INSPECT SCHEMA
# ------------------------------------------------------------
say("Inspecting first valid file: %s", file_info$file[1])

first_dt <- readRDS(file_info$path[1])
setDT(first_dt)

say("Columns in first valid file: %d", ncol(first_dt))
say("Rows in first valid file: %s", format(nrow(first_dt), big.mark = ","))

pid_cols <- intersect(c("era_pid", "pm_pid", "gee_pid"), names(first_dt))
coord_cols <- intersect(c("gee_longitude", "gee_latitude"), names(first_dt))

say("PID columns present: %s",
    if (length(pid_cols)) paste(pid_cols, collapse = ", ") else "NONE")
say("GEE coord columns present: %s",
    if (length(coord_cols)) paste(coord_cols, collapse = ", ") else "NONE")

rm(first_dt)
gc()

# ------------------------------------------------------------
# 4) COMBINE FILES
# ------------------------------------------------------------
say("Combining valid monthly RDS files...")

big <- NULL
row_diag <- vector("list", nrow(file_info))

for (i in seq_len(nrow(file_info))) {
  f <- file_info$path[i]
  m <- file_info$yyyymm[i]

  say("Reading %s (%d / %d)", basename(f), i, nrow(file_info))

  dt <- readRDS(f)
  setDT(dt)

  row_diag[[i]] <- data.table(
    yyyymm = m,
    n_rows = nrow(dt),
    n_cols = ncol(dt),
    has_era_pid = "era_pid" %in% names(dt),
    has_pm_pid  = "pm_pid"  %in% names(dt),
    has_gee_pid = "gee_pid" %in% names(dt),
    has_gee_lon = "gee_longitude" %in% names(dt),
    has_gee_lat = "gee_latitude"  %in% names(dt)
  )

  say("  rows: %s | cols: %d", format(nrow(dt), big.mark = ","), ncol(dt))

  if (is.null(big)) {
    big <- copy(dt)
  } else {
    big <- rbindlist(list(big, dt), use.names = TRUE, fill = TRUE)
  }

  rm(dt)
  gc()
}

row_diag <- rbindlist(row_diag, use.names = TRUE, fill = TRUE)

# ------------------------------------------------------------
# 5) FINAL DIAGNOSTICS
# ------------------------------------------------------------
say("Combined rows: %s", format(nrow(big), big.mark = ","))
say("Combined cols: %d", ncol(big))

if ("datetime_utc" %in% names(big)) {
  dt_range <- range(as.POSIXct(big$datetime_utc, tz = "UTC"), na.rm = TRUE)
  say("Datetime range: %s -> %s",
      format(dt_range[1], tz = "UTC"),
      format(dt_range[2], tz = "UTC"))
}

if ("station" %in% names(big)) {
  say("Unique stations: %s", format(uniqueN(big$station), big.mark = ","))
}

if ("region" %in% names(big)) {
  reg <- sort(unique(big$region))
  say("Regions: %s", paste(reg, collapse = ", "))
}

if ("gee_pid" %in% names(big)) {
  say("Rows with gee_pid: %s", format(sum(!is.na(big$gee_pid)), big.mark = ","))
}
if ("era_pid" %in% names(big)) {
  say("Rows with era_pid: %s", format(sum(!is.na(big$era_pid)), big.mark = ","))
}
if ("pm_pid" %in% names(big)) {
  say("Rows with pm_pid: %s", format(sum(!is.na(big$pm_pid)), big.mark = ","))
}
if (all(c("gee_longitude", "gee_latitude") %in% names(big))) {
  say(
    "Rows with gee coords: %s",
    format(sum(!is.na(big$gee_longitude) & !is.na(big$gee_latitude)), big.mark = ",")
  )
}

say("Per-month diagnostics:")
print(row_diag)

# ------------------------------------------------------------
# 6) SAVE
# ------------------------------------------------------------
say("Saving combined RDS: %s", out_rds)
saveRDS(big, out_rds, compress = "xz")

say("✅ Combined JOINED RDS written to:")
say("%s", out_rds)
