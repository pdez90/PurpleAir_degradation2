# ============================================================
# 02_ab_asymmetry_tests.R
#
# SHARPER A/B CHANNEL-ASYMMETRY TESTS
#
# You found the simple signed A-B trend null. Two tests below have more
# power to detect channel-specific mechanisms, because they condition on
# the variable that should drive asymmetry (temperature) or isolate
# channel-specific NOISE rather than channel-specific BIAS:
#
#   TEST 1 — Temperature-dependent signed bias.
#     Channel A sits nearest the temperature sensor / electronics. If
#     thermal effects matter, signed bias (A-B) should covary with
#     ambient T *within sensor*, and that coupling should strengthen
#     with operational age. A flat age trend in mean bias can hide this,
#     because positive coupling at high T and negative at low T average
#     to ~zero.
#
#   TEST 2 — Channel-specific noise via lag-1 differencing.
#     sd(diff(A)) vs sd(diff(B)) within station-bin isolates
#     high-frequency channel noise (Allan-variance flavor), independent
#     of shared drift. If one channel's optics degrade faster, its noise
#     ratio drifts from 1 with age.
#
# EITHER WAY THE PAPER WINS: nulls here are evidence that aging acts
# through the SHARED pathway (dust accumulation in the common airflow
# path) rather than channel-specific electronics — worth one sentence in
# the Discussion:
#
#   "Signed channel bias showed no age trend, no strengthening
#    temperature dependence, and no divergence in channel-specific
#    noise (section SX), indicating that aging predominantly affects
#    the shared sampling pathway — consistent with particle
#    accumulation in the common airflow path — rather than
#    channel-specific electronics."
#
# INPUT:  outdoor_final.RData
# OUTPUT: ab_asymmetry/ CSVs + figures
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(fixest)
})

global_rdata <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_final.RData"
outdir <- "/Users/priyanka/Downloads/PurpleAir_EPA/ab_asymmetry"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

bin_w <- 24L
PM_MAX <- 1500
MIN_OBS_BIN_NOISE <- 12L   # need enough hours in a bin for sd(diff)

a_col_default <- "pm25_cf_1_a_hourly(ug/m3)_pm25_cf_1_a"
b_col_default <- "pm25_cf_1_b_hourly(ug/m3)_pm25_cf_1_b"
temp_col <- "temperature_hourly(C)_temperature"

say <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")),
      paste0(..., collapse = ""), "\n")
}


load_outdoor_object <- function(infile) {
  # outdoor_final.RData must already exist and must be the paper's exact
  # analytic dataset (452,143,148 rows / 36,328 sensors). It is NOT rebuilt
  # here: an ad-hoc rebuild risks producing a file that differs from the
  # published dataset (e.g. omitting the cf_1 > 1500 ug/m3 exclusion).
  if (!file.exists(infile)) {
    stop("outdoor_final.RData not found at: ", infile,
         "\n  Run 00_build_analysis_dataset.R first; it reproduces and",
         "\n  verifies the paper's exact analytic dataset.")
  }
  e <- new.env(parent = emptyenv())
  load(infile, envir = e)
  if (exists("outdoor", envir = e)) x <- get("outdoor", envir = e)
  else if (exists("outdoor_final", envir = e)) x <- get("outdoor_final", envir = e)
  else stop("Could not find outdoor or outdoor_final in: ", infile)
  setDT(x)
  x[]
}

say("Loading: ", global_rdata)
outdoor <- load_outdoor_object(global_rdata)

a_col <- if (a_col_default %in% names(outdoor)) a_col_default else "cf1_a"
b_col <- if (b_col_default %in% names(outdoor)) b_col_default else "cf1_b"
has_temp <- temp_col %in% names(outdoor)
if (!has_temp) say("WARNING: temperature column not found; Test 1 will be skipped")

keep_cols <- c("station", "region", "datetime_utc", a_col, b_col,
               if (has_temp) temp_col)
dt <- outdoor[, ..keep_cols]
rm(outdoor); gc()

setnames(dt, c(a_col, b_col), c("cf1_a", "cf1_b"))
if (has_temp) setnames(dt, temp_col, "pa_temp")

dt[, `:=`(
  station = as.character(station),
  ts = as.POSIXct(datetime_utc, tz = "UTC"),
  cf1_a = as.numeric(cf1_a),
  cf1_b = as.numeric(cf1_b)
)]
dt[, datetime_utc := NULL]
dt <- dt[is.finite(cf1_a) & is.finite(cf1_b) &
           cf1_a >= 0 & cf1_b >= 0 & cf1_a <= PM_MAX & cf1_b <= PM_MAX]

setorder(dt, station, ts)
dt[, first_ts := ts[1L], by = station]
dt[, op_hours := as.numeric(difftime(ts, first_ts, units = "hours"))]
dt[, op_bin := floor(op_hours / bin_w) * bin_w]
dt[, op_years := op_bin / (365.25 * 24)]
dt[, bias_cf1 := cf1_a - cf1_b]

# ------------------------------------------------------------
# TEST 1: temperature-dependent signed bias, and its age trend
# ------------------------------------------------------------
if (has_temp) {
  say("Test 1: T-dependence of signed bias within sensor")
  dt[, pa_temp := as.numeric(pa_temp)]

  db <- dt[is.finite(pa_temp),
           .(n_obs = .N,
             mean_bias = mean(bias_cf1),
             mean_temp = mean(pa_temp)),
           by = .(station, region, op_bin, op_years)]

  # (a) Does bias track T at all, within sensor?
  m1 <- feols(mean_bias ~ mean_temp | station,
              data = db, weights = ~n_obs, cluster = ~station)

  # (b) Does the T-coupling strengthen with age? (interaction)
  m2 <- feols(mean_bias ~ mean_temp * op_years | station,
              data = db, weights = ~n_obs, cluster = ~station)

  say("Test 1a  bias ~ T | station:")
  print(coeftable(m1))
  say("Test 1b  bias ~ T x age | station:")
  print(coeftable(m2))

  t1 <- rbind(
    data.table(test = "1a: bias~T within sensor", term = rownames(coeftable(m1)),
               coeftable(m1)),
    data.table(test = "1b: bias~T x age within sensor", term = rownames(coeftable(m2)),
               coeftable(m2))
  )
  fwrite(t1, file.path(outdir, "test1_bias_temperature_coupling.csv"))

  # visual: binned bias vs T by age tertile
  db[, age_group := cut(op_years, c(-Inf, 1, 3, Inf),
                        labels = c("<1 yr", "1-3 yr", ">3 yr"))]
  db[, t_bin := round(mean_temp / 5) * 5]
  viz <- db[, .(bias = weighted.mean(mean_bias, n_obs), n = sum(n_obs)),
            by = .(age_group, t_bin)][n > 1000]
  p1 <- ggplot(viz, aes(t_bin, bias, color = age_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_line() + geom_point(size = 1.5) +
    labs(x = "PurpleAir temperature (°C, 5° bins)",
         y = expression("Mean signed bias A-B (" * mu * "g/" * m^3 * ")"),
         color = "Operational age",
         title = "Signed channel bias vs temperature by operational age") +
    theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(outdir, "bias_vs_temperature_by_age.jpeg"), p1,
         width = 7, height = 4.5, dpi = 300)
}

# ------------------------------------------------------------
# TEST 2: channel-specific noise (lag-1 differencing) vs age
# ------------------------------------------------------------
say("Test 2: channel noise ratio via lag-1 differencing")

# only difference consecutive hours (gaps break the differencing)
dt[, dt_hours := c(NA_real_, diff(op_hours)), by = station]
dt[, `:=`(dA = c(NA_real_, diff(cf1_a)), dB = c(NA_real_, diff(cf1_b))),
   by = station]
dt[dt_hours != 1, `:=`(dA = NA_real_, dB = NA_real_)]

nb <- dt[is.finite(dA) & is.finite(dB),
         .(n_obs = .N,
           sd_dA = sd(dA), sd_dB = sd(dB)),
         by = .(station, region, op_bin, op_years)][n_obs >= MIN_OBS_BIN_NOISE]

nb[, log_noise_ratio := log(pmax(sd_dA, 1e-6) / pmax(sd_dB, 1e-6))]

# Age trend in log noise ratio (0 = symmetric), sensor FE soaks up
# any fixed channel offset; the question is the *trend*.
m3 <- feols(log_noise_ratio ~ op_years | station,
            data = nb, weights = ~n_obs, cluster = ~station)
say("Test 2  log(sd(dA)/sd(dB)) ~ age | station:")
print(coeftable(m3))

t2 <- data.table(test = "2: log noise ratio ~ age within sensor",
                 term = rownames(coeftable(m3)), coeftable(m3))
fwrite(t2, file.path(outdir, "test2_noise_ratio_age_trend.csv"))

viz2 <- nb[, .(ratio = exp(weighted.mean(log_noise_ratio, n_obs)),
               n_sensors = uniqueN(station)),
           by = .(op_year_bin = floor(op_years * 2) / 2)][n_sensors >= 100]
p2 <- ggplot(viz2, aes(op_year_bin, ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_line() + geom_point(size = 1.5) +
  labs(x = "Operational age (years)",
       y = "Geometric mean sd(ΔA)/sd(ΔB)",
       title = "Channel-specific noise ratio by operational age",
       subtitle = "1.0 = symmetric channel noise; bins with ≥100 sensors") +
  theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold"))
ggsave(file.path(outdir, "noise_ratio_by_age.jpeg"), p2,
       width = 7, height = 4.5, dpi = 300)

say("Done. Outputs in: ", outdir)
