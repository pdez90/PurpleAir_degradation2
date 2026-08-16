# ============================================================
# 03_prospective_flag_metric.R
#
# PROSPECTIVELY APPLICABLE FLAGGED-DEGRADATION VARIANTS
#
# Problem: the paper's flag uses the sensor-specific 85th percentile of
# the symmetric percent difference computed over the sensor's FULL
# record — unknowable in real time — while the Discussion recommends
# real-time screening. This script evaluates two variants a network
# operator could actually run, and shows the paper's conclusions are
# insensitive to the choice:
#
#   RETRO  (paper):  pct > sensor-specific 85th pct of full record
#   FIXED  (variant 1): pct > one global constant = the network-wide
#          85th percentile of pctdiff (a number you report once;
#          new sensors can be screened from hour 1)
#   EARLY  (variant 2): pct > sensor-specific 85th pct computed from the
#          first 90 operational days only (a burn-in calibration window,
#          then applicable prospectively for the sensor's whole life)
#
# For each variant: (i) hourly agreement with the retrospective flag
# (percent agreement + Cohen's kappa), (ii) sensor-level flagged-rate
# correlation, (iii) re-estimated global/regional age trends (OR per
# operational year, quasibinomial, sensor-clustered SEs) side by side.
#
# Release this script with the paper (code availability) so operators
# can implement the screening directly.
#
# INPUT:  outdoor_final.RData
# OUTPUT: prospective_flag/ CSVs + Word-ready comparison table
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(sandwich)
  library(lmtest)
  library(flextable)
  library(officer)
})

global_rdata <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_final.RData"
outdir <- "/Users/priyanka/Downloads/PurpleAir_EPA/prospective_flag"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

bin_w <- 24L
PM_MAX <- 1500
FLAG_ABS_CUTOFF <- 5
FLAG_PCTL <- 0.85
EARLY_WINDOW_HOURS <- 90L * 24L   # 90-day burn-in
MIN_EARLY_OBS <- 240L             # >=240 valid hours in burn-in to calibrate

a_col_default <- "pm25_cf_1_a_hourly(ug/m3)_pm25_cf_1_a"
b_col_default <- "pm25_cf_1_b_hourly(ug/m3)_pm25_cf_1_b"

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
dt <- outdoor[, c("station", "region", "datetime_utc", a_col, b_col), with = FALSE]
rm(outdoor); gc()

setnames(dt, c(a_col, b_col), c("cf1_a", "cf1_b"))
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

dt[, `:=`(
  absdiff_cf1 = abs(cf1_a - cf1_b),
  pctdiff_cf1 = fifelse((cf1_a + cf1_b) > 0,
                        2 * abs(cf1_a - cf1_b) / (cf1_a + cf1_b), NA_real_)
)]

# ------------------------------------------------------------
# 1) THREE THRESHOLDS
# ------------------------------------------------------------

# RETRO (paper)
thr_retro <- dt[is.finite(pctdiff_cf1),
                .(thr_retro = quantile(pctdiff_cf1, FLAG_PCTL, na.rm = TRUE, type = 7)),
                by = station]

# FIXED global constant — report this number in the paper
THR_FIXED <- dt[is.finite(pctdiff_cf1),
                quantile(pctdiff_cf1, FLAG_PCTL, na.rm = TRUE, type = 7)]
say(sprintf("Network-wide fixed threshold (85th pct of pctdiff): %.4f", THR_FIXED))

# EARLY sensor-specific baseline (first 90 operational days)
thr_early <- dt[is.finite(pctdiff_cf1) & op_hours <= EARLY_WINDOW_HOURS,
                .(thr_early = quantile(pctdiff_cf1, FLAG_PCTL, na.rm = TRUE, type = 7),
                  n_early = .N),
                by = station][n_early >= MIN_EARLY_OBS]

dt <- thr_retro[dt, on = "station"]
dt <- thr_early[dt, on = "station"]

dt[, flag_retro := as.integer(absdiff_cf1 > FLAG_ABS_CUTOFF & pctdiff_cf1 > thr_retro)]
dt[, flag_fixed := as.integer(absdiff_cf1 > FLAG_ABS_CUTOFF & pctdiff_cf1 > THR_FIXED)]
# Strictly prospective: the threshold is estimated on the first
# EARLY_WINDOW_HOURS of each sensor's record, so those hours are left NA and
# flagging begins only after the burn-in window.
dt[, flag_early := fifelse(
  is.finite(thr_early) & op_hours > EARLY_WINDOW_HOURS,
  as.integer(absdiff_cf1 > FLAG_ABS_CUTOFF & pctdiff_cf1 > thr_early),
  NA_integer_)]

# ------------------------------------------------------------
# 2) AGREEMENT WITH RETROSPECTIVE FLAG
# ------------------------------------------------------------

kappa2 <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]; y <- y[ok]
  po <- mean(x == y)
  pe <- mean(x) * mean(y) + (1 - mean(x)) * (1 - mean(y))
  c(agree = po, kappa = (po - pe) / (1 - pe), n = length(x))
}

agr_fixed <- dt[, as.list(kappa2(flag_retro, flag_fixed))]
agr_early <- dt[!is.na(flag_early), as.list(kappa2(flag_retro, flag_early))]

# sensor-level flagged-rate correlations
sr <- dt[, .(
  rate_retro = mean(flag_retro, na.rm = TRUE),
  rate_fixed = mean(flag_fixed, na.rm = TRUE),
  rate_early = mean(flag_early, na.rm = TRUE),
  n = .N
), by = station][n >= 1000]

agreement <- data.table(
  variant = c("Fixed global threshold", "Early-life baseline threshold"),
  threshold = c(sprintf("pctdiff > %.3f (network 85th pct)", THR_FIXED),
                "pctdiff > sensor 85th pct of first 90 days"),
  pct_agreement = c(agr_fixed$agree, agr_early$agree),
  cohens_kappa = c(agr_fixed$kappa, agr_early$kappa),
  sensor_rate_correlation = c(
    sr[, cor(rate_retro, rate_fixed, use = "complete.obs")],
    sr[, cor(rate_retro, rate_early, use = "complete.obs")]),
  n_hours = c(agr_fixed$n, agr_early$n)
)
fwrite(agreement, file.path(outdir, "flag_variant_agreement.csv"))
print(agreement)

# ------------------------------------------------------------
# 3) RE-ESTIMATED AGE TRENDS UNDER EACH VARIANT
# ------------------------------------------------------------

say("Re-estimating age trends (OR per operational year) per variant")

sens_bin <- dt[, .(
  n_obs = .N,
  n_retro = sum(flag_retro, na.rm = TRUE),
  n_fixed = sum(flag_fixed, na.rm = TRUE),
  n_early = sum(flag_early, na.rm = TRUE),
  n_early_valid = sum(!is.na(flag_early))
), by = .(station, region, op_bin, op_years)]

rm(dt); gc()

or_trend <- function(dd, num_col, den_col, variant, group) {
  dd <- dd[get(den_col) > 0]
  if (nrow(dd) < 20L || uniqueN(dd$station) < 10L) return(data.table())
  fit <- tryCatch(
    glm(cbind(get(num_col), get(den_col) - get(num_col)) ~ op_years,
        data = dd, family = quasibinomial()),
    error = function(e) NULL)
  if (is.null(fit)) return(data.table())
  ct <- tryCatch({
    vc <- sandwich::vcovCL(fit, cluster = dd$station, type = "HC1")
    lmtest::coeftest(fit, vcov. = vc)
  }, error = function(e) NULL)
  if (is.null(ct) || !("op_years" %in% rownames(ct))) return(data.table())
  b <- ct["op_years", "Estimate"]; s <- ct["op_years", "Std. Error"]
  data.table(group = group, variant = variant,
             or = exp(b), ci_lwr = exp(b - 1.96 * s), ci_upr = exp(b + 1.96 * s),
             n_sensors = uniqueN(dd$station))
}

run_group <- function(dd, group) rbindlist(list(
  or_trend(dd, "n_retro", "n_obs", "Retrospective (paper)", group),
  or_trend(dd, "n_fixed", "n_obs", "Fixed global", group),
  or_trend(dd, "n_early", "n_early_valid", "Early-life baseline", group)
))

res <- rbindlist(c(
  list(run_group(sens_bin, "Global")),
  lapply(sort(unique(sens_bin$region)), function(r)
    run_group(sens_bin[region == r], r))
))
fwrite(res, file.path(outdir, "flag_variant_age_trends.csv"))
print(res)

# ------------------------------------------------------------
# 4) WORD TABLE (SI)
# ------------------------------------------------------------

wide <- dcast(res, group ~ variant,
              value.var = c("or", "ci_lwr", "ci_upr"))
fmt <- function(o, l, u) fifelse(is.na(o), "-",
                                 sprintf("%.2f (%.2f, %.2f)", o, l, u))
ft_dt <- wide[, .(
  Group = group,
  `Retrospective (paper)` = fmt(`or_Retrospective (paper)`,
                                `ci_lwr_Retrospective (paper)`,
                                `ci_upr_Retrospective (paper)`),
  `Fixed global` = fmt(`or_Fixed global`, `ci_lwr_Fixed global`,
                       `ci_upr_Fixed global`),
  `Early-life baseline` = fmt(`or_Early-life baseline`,
                              `ci_lwr_Early-life baseline`,
                              `ci_upr_Early-life baseline`)
)]

ft <- autofit(fontsize(flextable(ft_dt), size = 9, part = "all"))
doc <- read_docx()
doc <- body_add_par(doc, paste0(
  "Table SX: Age trends in flagged degradation (OR per operational year, ",
  "95% CI, sensor-clustered SEs) under the retrospective sensor-specific ",
  "threshold used in the main analysis and two prospectively applicable ",
  "variants: a fixed network-wide threshold (pctdiff > ",
  sprintf("%.3f", THR_FIXED), ") and a sensor-specific threshold ",
  "calibrated on the first 90 operational days."), style = "Normal")
doc <- body_add_flextable(doc, ft)
doc <- body_add_par(doc, "")
doc <- body_add_par(doc, paste0(
  "Hourly agreement with the retrospective flag: fixed global ",
  sprintf("%.1f%%", 100 * agreement$pct_agreement[1]),
  " (kappa = ", sprintf("%.2f", agreement$cohens_kappa[1]),
  "); early-life baseline ",
  sprintf("%.1f%%", 100 * agreement$pct_agreement[2]),
  " (kappa = ", sprintf("%.2f", agreement$cohens_kappa[2]), ")."),
  style = "Normal")
print(doc, target = file.path(outdir, "flag_variant_comparison.docx"))
say("Wrote: flag_variant_comparison.docx — done.")
