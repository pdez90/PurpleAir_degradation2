# ============================================================================
# verify_paper_numbers.R
#
# Compares numbers in the analysis OUTPUT files against the values printed in
# the manuscript and SI. Reads existing outputs only — runs no models, so it
# takes seconds. Run it after any script rerun.
#
#   Rscript verify_paper_numbers.R
# ============================================================================
suppressPackageStartupMessages(library(data.table))
ROOT <- "/Users/priyanka/Downloads/PurpleAir_EPA"   # <- edit if your data dir differs

pass <- 0L; fail <- 0L; skip <- 0L
chk <- function(label, actual, expected, tol = 0.005) {
  if (is.null(actual) || length(actual) == 0 || all(is.na(actual))) {
    cat(sprintf("  [SKIP] %-58s (output not found)\n", label)); skip <<- skip + 1L; return(invisible())
  }
  ok <- isTRUE(abs(as.numeric(actual[1]) - expected) <= tol)
  cat(sprintf("  [%s] %-58s paper=%-12s output=%s\n", if (ok) "PASS" else "FAIL",
              label, format(expected), format(round(as.numeric(actual[1]), 4))))
  if (ok) pass <<- pass + 1L else fail <<- fail + 1L
}
rd <- function(...) { p <- file.path(ROOT, ...); if (file.exists(p)) fread(p) else NULL }

cat("\n=== SENSITIVITY / SI ADD-ON ANALYSES ===\n")
x <- rd("age_period_cohort","age_period_cohort_results.csv")
if (!is.null(x)) {
  chk("S7.1 |A-B| age slope, age+period (+0.50)",
      x[group=="Global" & spec=="Age + period" & grepl("Absolute", outcome), beta_age], 0.50, 0.01)
  chk("S7.1 calendar-time effect (~ -0.013)",
      x[group=="Global" & spec=="Age + period" & grepl("Absolute", outcome), beta_period], -0.013, 0.01)
} else chk("S7.1 age slope", NULL, 0.50)

x <- rd("ab_asymmetry","test1_bias_temperature_coupling.csv")
chk("S8 bias~T p-value (0.98)", if (!is.null(x)) x[term=="mean_temp"][1]$`Pr(>|t|)` else NULL, 0.976, 0.02)
x <- rd("ab_asymmetry","test2_noise_ratio_age_trend.csv")
chk("S8 noise-ratio age slope (0.0029)", if (!is.null(x)) x$Estimate else NULL, 0.0029, 0.001)

x <- rd("prospective_flag","flag_variant_age_trends.csv")
if (!is.null(x)) {
  chk("S6.1 Global retrospective OR (1.33)", x[group=="Global" & variant=="Retrospective (paper)", or], 1.327, 0.01)
  chk("S6.1 Global early-life OR (1.29, post-fix)", x[group=="Global" & variant=="Early-life baseline", or], 1.287, 0.01)
  chk("S6.1 India early-life OR (1.01, post-fix)", x[group=="India" & variant=="Early-life baseline", or], 1.015, 0.01)
} else chk("S6.1 ORs", NULL, 1.33)

x <- rd("wildfire_stratification","wildfire_stratified_models.csv")
if (!is.null(x)) {
  chk("S9.1 West fire-season burden (+2.7)",
      x[grepl("Western", stratum) & grepl("Fire", stratum) & outcome=="fut_abs" & term=="burden100", estimate], 2.743, 0.05)
  chk("S9.1 West non-fire burden (+19.0)",
      x[grepl("Western", stratum) & grepl("Non-fire", stratum) & outcome=="fut_abs" & term=="burden100", estimate], 19.02, 0.1)
}
x <- rd("wildfire_stratification","wildfire_interaction_model.csv")
chk("S9 burden x fire-share interaction (-14.9)",
    if (!is.null(x)) x[term=="burden100:fire_share", Estimate] else NULL, -14.94, 0.1)

x <- rd("exposure_window_sensitivity","exposure_window_sensitivity.csv")
if (!is.null(x)) {
  chk("S9.2 7d reported/drop |A-B| (+2.45)",
      x[past_window_days==7 & grepl("Reported", windows) & grepl("drop", exclusion) & outcome=="Future |A-B|", est], 2.453, 0.02)
  chk("S9.2 7d reported/keep |A-B| (+3.45)",
      x[past_window_days==7 & grepl("Reported", windows) & grepl("keep", exclusion) & outcome=="Future |A-B|", est], 3.448, 0.02)
}

x <- rd("pa_epa_cv","pa_epa_calibration_cv.csv")
if (!is.null(x)) {
  chk("S3.2 GAM CV RMSE (6.88)", x[model=="GAM s(PA)" & grepl("CV", eval), rmse], 6.878, 0.02)
  chk("S3.2 tensor CV r (0.810)", x[model=="Age-aware tensor GAM" & grepl("CV", eval), r], 0.8104, 0.005)
}

cat("\n=== MAIN-TEXT: AGE-DEGRADATION MODELS (Fig 2; Tables S3-S5) ===\n")
g <- rd("cf1_degradation_split_inputs","degradation_time_associations_global_cluster_robust.csv")
if (!is.null(g)) {
  gv <- function(oc, col="estimate") g[grepl(oc, outcome)][1][[col]]
  chk("Global |A-B| slope (+0.49)",            gv("Absolute"), 0.4917, 0.01)
  chk("Global |A-B| CI lower (0.44)",          gv("Absolute","ci_lwr"), 0.4429, 0.01)
  chk("Global |A-B| CI upper (0.54)",          gv("Absolute","ci_upr"), 0.5405, 0.01)
  chk("Global signed bias (+0.005, ns)",       gv("Signed"), 0.0048, 0.01)
  chk("Global signed-bias p (0.85)",           gv("Signed","p_value"), 0.845, 0.02)
  chk("Global relative divergence (+0.054)",   gv("Relative"), 0.0542, 0.002)
  chk("Global flagged OR (1.33)",              gv("Flagged"), 1.3268, 0.01)
  chk("Channel A decline (~-0.95)",            gv("Channel A"), -0.9465, 0.01)
  chk("Channel B decline (~-0.95)",            gv("Channel B"), -0.9513, 0.01)
  chk("Analytic sensors (36,328)",             gv("Absolute","n_sensors"), 36328, 0)
}
r <- rd("cf1_degradation_split_inputs","degradation_time_associations_by_region_cluster_robust.csv")
if (!is.null(r)) {
  rv <- function(reg, oc, col="estimate") r[grepl(reg, group) & grepl(oc, outcome)][1][[col]]
  chk("Africa flagged OR (1.55)",  rv("Africa","Flagged"), 1.5508, 0.01)
  chk("CONUS flagged OR (1.35)",   rv("CONUS","Flagged"),  1.3478, 0.01)
  chk("India flagged OR (1.24)",   rv("India","Flagged"),  1.2388, 0.01)
  chk("S.America flagged OR (1.37)", rv("SouthAmerica","Flagged"), 1.3723, 0.01)
  chk("Africa |A-B| (+1.68)",      rv("Africa","Absolute"), 1.6754, 0.02)
  chk("CONUS |A-B| (+0.52)",       rv("CONUS","Absolute"),  0.5154, 0.02)
}
cz <- rd("cf1_degradation_split_inputs","degradation_time_associations_by_climate_zone_cluster_robust.csv")
if (!is.null(cz)) {
  zv <- function(z, oc, col="estimate") cz[grepl(paste0("Climate: ", z, " "), group) & grepl(oc, outcome)][1][[col]]
  chk("Dfb flagged OR (1.25)",  zv("Dfb","Flagged"), 1.252, 0.01)
  chk("Dfc flagged OR (0.98 ns)", zv("Dfc","Flagged"), 0.981, 0.01)
  chk("Cwb |A-B| (+1.93)",      zv("Cwb","Absolute"), 1.929, 0.02)
  chk("Am |A-B| (+1.16)",       zv("Am","Absolute"),  1.157, 0.02)
}

cat("\n=== HIGH-PM TABLE ROW COUNTS (S19-S23) ===\n")
for (tb in c("Table_S19_global_high_pm_associations","Table_S20_regional_joint_28day_models",
             "Table_S21_climate_joint_28day_models","Table_S22_region_x_climate_joint_28day_models",
             "Table_S23_joint_model_signal_summary")) {
  y <- rd("high_pm_future_degradation_combined","word_ready_tables", paste0(tb,".csv"))
  cat(sprintf("  %-46s %s rows\n", tb, if (!is.null(y)) nrow(y) else "MISSING"))
}

cat(sprintf("\n==== PASS: %d   FAIL: %d   SKIP: %d ====\n", pass, fail, skip))
if (fail > 0) cat("Investigate FAIL rows: the output no longer matches the printed value.\n")
if (skip > 0) cat("SKIP = that script's outputs are absent; rerun it if you want it verified.\n")
