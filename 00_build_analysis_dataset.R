# ============================================================
# 00_build_analysis_dataset.R
#
# RESTORE outdoor_final.RData TO EXACTLY MATCH THE PAPER'S VERSION
#
# Background: the original outdoor_final.RData (created by the Rmd
# cleaning chunk) was deleted at some point. The add-on scripts rebuilt
# a working substitute from outdoor_all_joined.rds, but that substitute
# differs from the paper's file in two ways that matter if you re-run
# Rmd chunks that load outdoor_final.RData directly:
#   1) it does NOT apply the bad_pm filter (cf1 A or B > 1500), so it
#      has 464,848,642 rows instead of 452,143,148; and
#   2) it lacks the derived columns pm25_cf1_avg, era5_t2m_c,
#      era5_td2m_c, era5_rh_pct, and the gee coords / pid columns that
#      the model-divergence and summary-figure chunks expect.
# (The add-on scripts 01-04 are unaffected: they re-apply the 1500
# filter internally and reproduce the paper's exact sensor counts.)
#
# This script replicates the Rmd cleaning chunk verbatim (same keep
# list, same filter order, same derived columns) and overwrites
# outdoor_final.RData. Built-in verification: row and region counts are
# compared against the numbers recorded in the Rmd comments.
#
# Runtime: ~40-60 min (reads the 16 GB joined rds); needs plenty of RAM.
# ============================================================

suppressPackageStartupMessages(library(data.table))

joined_rds <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_all_joined.rds"
out_rdata  <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_final.RData"

say <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")),
                         paste0(..., collapse = ""), "\n")

a_col   <- "pm25_cf_1_a_hourly(ug/m3)_pm25_cf_1_a"
b_col   <- "pm25_cf_1_b_hourly(ug/m3)_pm25_cf_1_b"
tcol    <- "temperature_hourly(C)_temperature"
rhcol   <- "humidity_hourly(%)_humidity"
era5_tcol  <- "temperature_2m"
era5_tdcol <- "dewpoint_temperature_2m"
geoscf_col <- "geoscf_pm25"

say("Reading: ", joined_rds)
outdoor <- readRDS(joined_rds)
setDT(outdoor)

pa_lon <- intersect(c("LONGITUDE(deg)", "longitude"), names(outdoor))
pa_lat <- intersect(c("LATITUDE(deg)",  "latitude"),  names(outdoor))
gee_lon <- intersect(c("gee_longitude"), names(outdoor))
gee_lat <- intersect(c("gee_latitude"),  names(outdoor))

keep <- unique(intersect(
  c("datetime_utc", "station", "region",
    pa_lon, pa_lat, gee_lon, gee_lat,
    "era_pid", "pm_pid", "gee_pid",
    a_col, b_col, tcol, rhcol,
    era5_tcol, era5_tdcol, geoscf_col),
  names(outdoor)
))
outdoor <- outdoor[, ..keep]
gc()

# --- same filter order as the Rmd chunk ---
outdoor <- outdoor[!is.na(get(a_col)) & !is.na(get(b_col))]
gc()
outdoor <- outdoor[!is.na(get(tcol)) & !is.na(get(rhcol))]
gc()
outdoor[, bad_T  := !is.na(get(tcol))  & (get(tcol) <= -50 | get(tcol) >= 100)]
outdoor[, bad_RH := !is.na(get(rhcol)) & (get(rhcol) > 99)]
outdoor <- outdoor[bad_T == FALSE & bad_RH == FALSE][, c("bad_T","bad_RH") := NULL]
gc()
say("Rows after T/RH filters (expect 464,848,642): ",
    format(nrow(outdoor), big.mark = ","))

outdoor[, bad_pm := get(a_col) > 1500 | get(b_col) > 1500]
outdoor_final <- outdoor[bad_pm == FALSE][, bad_pm := NULL]
rm(outdoor); gc()

# --- verification against the paper's recorded counts ---
n <- nrow(outdoor_final)
say("Final rows (expect 452,143,148): ", format(n, big.mark = ","))
print(table(outdoor_final$region))
say("Expect: Africa 3,134,219 | CONUS 439,768,685 | India 2,798,734 | SouthAmerica 6,441,510")
say("Sensors (expect 36,328): ", format(uniqueN(outdoor_final$station), big.mark = ","))
if (n != 452143148L) warning("ROW COUNT MISMATCH vs paper - investigate before overwriting!")

# --- derived columns, exactly as in the Rmd chunk ---
outdoor_final[, pm25_cf1_avg := (get(a_col) + get(b_col)) / 2]

if (era5_tcol %in% names(outdoor_final)) {
  med_t <- median(outdoor_final[[era5_tcol]], na.rm = TRUE)
  if (is.finite(med_t) && med_t > 150) {
    outdoor_final[, era5_t2m_c := get(era5_tcol) - 273.15]
  } else outdoor_final[, era5_t2m_c := get(era5_tcol)]
}
if (era5_tdcol %in% names(outdoor_final)) {
  med_td <- median(outdoor_final[[era5_tdcol]], na.rm = TRUE)
  if (is.finite(med_td) && med_td > 150) {
    outdoor_final[, era5_td2m_c := get(era5_tdcol) - 273.15]
  } else outdoor_final[, era5_td2m_c := get(era5_tdcol)]
}
if (all(c("era5_t2m_c", "era5_td2m_c") %in% names(outdoor_final))) {
  a <- 17.625; b <- 243.04
  outdoor_final[, era5_rh_pct := {
    T <- era5_t2m_c; Td <- era5_td2m_c
    rh <- 100 * exp((a * Td) / (b + Td) - (a * T) / (b + T))
    pmin(100, pmax(0, rh))
  }]
}

# --- duplicate station-hour sanity check (review flag 5) ---
ndup <- outdoor_final[, .N, by = .(station, datetime_utc)][N > 1L, .N]
say("Duplicate (station, hour) pairs (expect 0): ", format(ndup, big.mark = ","))
if (ndup > 0) warning("Duplicate station-hours present - investigate!")

say("Saving: ", out_rdata)
save(outdoor_final, file = out_rdata)
say("Done. Columns: ", paste(names(outdoor_final), collapse = ", "))
