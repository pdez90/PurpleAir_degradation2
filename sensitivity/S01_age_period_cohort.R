# ============================================================
# 01_age_vs_calendar_time.R
#
# SEPARATING OPERATIONAL AGE FROM CALENDAR TIME (age-period-cohort)
#
# Motivation (reviewer-proofing): operational age is confounded with
# calendar period and deployment cohort. A 6-year-old sensor observed in
# 2025 was necessarily deployed in 2019 and lived through specific smoke
# seasons, ambient trends, and firmware eras. This script shows the
# age-related degradation signal survives when calendar time is held
# fixed, three complementary ways:
#
#   A) WITHIN-CALENDAR-YEAR age trends: refit the Model 1 / Model 2
#      trend regressions separately inside each calendar year. Within a
#      single year, calendar period is (nearly) constant, so any age
#      gradient is age + cohort, not period.
#   B) JOINT age + period model: op_years and continuous calendar time
#      enter together (identified under the assumption of no cohort
#      effect on the slope).
#   C) COHORT-STRATIFIED age trends: deployment-year cohorts, age trend
#      within cohort, plus a fixest specification with deployment-year
#      fixed effects (identified under the assumption of no period
#      effect).
#
# Age = period - cohort, so no single model separates all three; the
# standard move is triangulation: if the op_years coefficient is stable
# across A, B, and C, the aging interpretation is safe. Report the
# comparison as one SI table + one SI forest-style figure.
#
# Also runs the GEOS-CF divergence outcome (the India concern: ambient/
# model trends masquerading as sensor drift) with month-of-year fixed
# effects to strip seasonality.
#
# INPUT:  outdoor_final.RData  (same object as cf1 degradation chunk)
# OUTPUT: age_period_cohort/ CSVs + figure + Word tables
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(fixest)
  library(sandwich)
  library(lmtest)
  library(flextable)
  library(officer)
})

# ------------------------------------------------------------
# 0) PATHS / SETTINGS (match main pipeline)
# ------------------------------------------------------------

global_rdata <- "/Users/priyanka/Downloads/PurpleAir_EPA/outdoor_final.RData"
outdir <- "/Users/priyanka/Downloads/PurpleAir_EPA/age_period_cohort"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

bin_w <- 24L
min_rows_for_trend <- 20L
min_stations_group <- 10L
MIN_SENSORS_YEAR <- 200L      # only report within-year slopes with support
PM_MAX <- 1500

a_col_default <- "pm25_cf_1_a_hourly(ug/m3)_pm25_cf_1_a"
b_col_default <- "pm25_cf_1_b_hourly(ug/m3)_pm25_cf_1_b"
geos_col <- "geoscf_pm25"

say <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")),
      paste0(..., collapse = ""), "\n")
}

# ------------------------------------------------------------
# 1) LOAD + BUILD HOURLY DT WITH AGE *AND* CALENDAR VARS
# ------------------------------------------------------------


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
has_geos <- geos_col %in% names(outdoor)

keep_cols <- c("station", "region", "datetime_utc", a_col, b_col,
               if (has_geos) geos_col)
dt <- outdoor[, ..keep_cols]
rm(outdoor); gc()

setnames(dt, c(a_col, b_col), c("cf1_a", "cf1_b"))
if (has_geos) setnames(dt, geos_col, "geos")

dt[, `:=`(
  station = as.character(station),
  region = as.character(region),
  ts = as.POSIXct(datetime_utc, tz = "UTC"),
  cf1_a = as.numeric(cf1_a),
  cf1_b = as.numeric(cf1_b)
)]
dt[, datetime_utc := NULL]

dt <- dt[
  !is.na(station) & !is.na(ts) &
    is.finite(cf1_a) & is.finite(cf1_b) &
    cf1_a >= 0 & cf1_b >= 0 & cf1_a <= PM_MAX & cf1_b <= PM_MAX
]

setorder(dt, station, ts)
dt[, first_ts := ts[1L], by = station]
dt[, op_hours := as.numeric(difftime(ts, first_ts, units = "hours"))]
dt[, op_bin := floor(op_hours / bin_w) * bin_w]
dt[, op_years := op_bin / (365.25 * 24)]
dt <- dt[is.finite(op_bin) & op_bin >= 0]

# --- calendar variables ---
dt[, cal_year := year(ts)]
dt[, cal_month := month(ts)]
# continuous calendar time in years since 2018-01-01 (period axis)
dt[, cal_time := as.numeric(difftime(ts, as.POSIXct("2018-01-01", tz = "UTC"),
                                     units = "days")) / 365.25]
# deployment cohort
dt[, deploy_year := year(first_ts)]

# --- degradation metrics (identical to main pipeline) ---
dt[, `:=`(
  absdiff_cf1 = abs(cf1_a - cf1_b),
  bias_cf1 = cf1_a - cf1_b,
  pctdiff_cf1 = fifelse((cf1_a + cf1_b) > 0,
                        2 * abs(cf1_a - cf1_b) / (cf1_a + cf1_b), NA_real_)
)]
thr_dt <- dt[is.finite(pctdiff_cf1),
             .(thr85 = quantile(pctdiff_cf1, 0.85, na.rm = TRUE, type = 7)),
             by = station]
dt <- thr_dt[dt, on = "station"]
dt[, flagged_cf1 := as.integer(absdiff_cf1 > 5 & pctdiff_cf1 > thr85)]
rm(thr_dt); gc()

if (has_geos) {
  dt[, md_cf1 := fifelse(is.finite(geos), (cf1_a + cf1_b) / 2 - geos, NA_real_)]
}

# ------------------------------------------------------------
# 2) STATION x OP-BIN x CALENDAR-YEAR SUMMARIES
#    (bins straddling New Year contribute one row per year — fine,
#     rows are weighted by n_obs)
# ------------------------------------------------------------

say("Aggregating station x op_bin x calendar-year summaries")

sens_bin <- dt[, .(
  n_obs = .N,
  n_flagged = sum(flagged_cf1, na.rm = TRUE),
  mean_abs_s = mean(absdiff_cf1, na.rm = TRUE),
  mean_bias_s = mean(bias_cf1, na.rm = TRUE),
  mean_pct_s = mean(pctdiff_cf1, na.rm = TRUE),
  mean_md_s = if (has_geos) mean(md_cf1, na.rm = TRUE) else NA_real_,
  cal_time = mean(cal_time),
  cal_month = as.integer(round(mean(cal_month)))
), by = .(station, region, deploy_year, cal_year, op_bin, op_years)]

rm(dt); gc()

# ------------------------------------------------------------
# 3) HELPERS: cluster-robust fits (mirror main pipeline)
# ------------------------------------------------------------

cluster_ct <- function(fit, dd) {
  vc <- sandwich::vcovCL(fit, cluster = dd$station, type = "HC1")
  lmtest::coeftest(fit, vcov. = vc)
}

grab <- function(ct, term) {
  if (is.null(ct) || !(term %in% rownames(ct))) return(c(NA_real_, NA_real_, NA_real_))
  p_col <- grep("^Pr", colnames(ct))
  c(ct[term, "Estimate"], ct[term, "Std. Error"],
    if (length(p_col)) ct[term, p_col[1]] else NA_real_)
}

fit_lm_row <- function(dd, yvar, rhs, spec, outcome, group) {
  dd <- dd[is.finite(get(yvar)) & is.finite(op_years) & n_obs > 0]
  if (nrow(dd) < min_rows_for_trend || uniqueN(dd$station) < min_stations_group)
    return(data.table())
  fit <- tryCatch(
    lm(as.formula(paste0(yvar, " ~ ", rhs)), data = dd, weights = n_obs),
    error = function(e) NULL)
  if (is.null(fit)) return(data.table())
  ct <- tryCatch(cluster_ct(fit, dd), error = function(e) NULL)
  age <- grab(ct, "op_years"); per <- grab(ct, "cal_time")
  data.table(
    group = group, spec = spec, outcome = outcome,
    effect_type = "Beta per operational year",
    beta_age = age[1], se_age = age[2], p_age = age[3],
    ci_age_lwr = age[1] - 1.96 * age[2], ci_age_upr = age[1] + 1.96 * age[2],
    beta_period = per[1], se_period = per[2],
    ci_period_lwr = per[1] - 1.96 * per[2], ci_period_upr = per[1] + 1.96 * per[2],
    n_rows = nrow(dd), n_sensors = uniqueN(dd$station)
  )
}

fit_flag_row <- function(dd, rhs, spec, group) {
  dd <- dd[is.finite(op_years) & n_obs > 0 & n_flagged >= 0 & n_flagged <= n_obs]
  if (nrow(dd) < min_rows_for_trend || uniqueN(dd$station) < min_stations_group)
    return(data.table())
  fit <- tryCatch(
    glm(as.formula(paste0("cbind(n_flagged, n_obs - n_flagged) ~ ", rhs)),
        data = dd, family = quasibinomial()),
    error = function(e) NULL)
  if (is.null(fit)) return(data.table())
  ct <- tryCatch(cluster_ct(fit, dd), error = function(e) NULL)
  age <- grab(ct, "op_years"); per <- grab(ct, "cal_time")
  data.table(
    group = group, spec = spec, outcome = "Flagged degradation",
    effect_type = "Odds ratio per operational year",
    beta_age = exp(age[1]), se_age = age[2], p_age = age[3],
    ci_age_lwr = exp(age[1] - 1.96 * age[2]), ci_age_upr = exp(age[1] + 1.96 * age[2]),
    beta_period = exp(per[1]), se_period = per[2],
    ci_period_lwr = exp(per[1] - 1.96 * per[2]), ci_period_upr = exp(per[1] + 1.96 * per[2]),
    n_rows = nrow(dd), n_sensors = uniqueN(dd$station)
  )
}

# ------------------------------------------------------------
# 4) ANALYSIS A: WITHIN-CALENDAR-YEAR AGE TRENDS
# ------------------------------------------------------------

say("Analysis A: within-calendar-year age trends")

res_within <- rbindlist(lapply(sort(unique(sens_bin$cal_year)), function(yy) {
  dd <- sens_bin[cal_year == yy]
  if (uniqueN(dd$station) < MIN_SENSORS_YEAR) return(data.table())
  rbindlist(list(
    fit_lm_row(dd, "mean_abs_s", "op_years",
               paste0("Within ", yy), "Absolute channel divergence |A-B|",
               paste0("Year ", yy)),
    fit_lm_row(dd, "mean_pct_s", "op_years",
               paste0("Within ", yy), "Relative channel divergence",
               paste0("Year ", yy)),
    fit_flag_row(dd, "op_years", paste0("Within ", yy), paste0("Year ", yy))
  ), fill = TRUE)
}), fill = TRUE)

# ------------------------------------------------------------
# 5) ANALYSIS B: JOINT AGE + PERIOD (global and by region)
# ------------------------------------------------------------

say("Analysis B: joint age + continuous calendar time")

groups <- c(list(list(label = "Global", dd = sens_bin)),
            lapply(sort(unique(sens_bin$region)), function(r)
              list(label = r, dd = sens_bin[region == r])))

res_joint <- rbindlist(lapply(groups, function(g) {
  rbindlist(list(
    fit_lm_row(g$dd, "mean_abs_s", "op_years + cal_time",
               "Age + period", "Absolute channel divergence |A-B|", g$label),
    fit_lm_row(g$dd, "mean_pct_s", "op_years + cal_time",
               "Age + period", "Relative channel divergence", g$label),
    fit_flag_row(g$dd, "op_years + cal_time", "Age + period", g$label)
  ), fill = TRUE)
}), fill = TRUE)

# ------------------------------------------------------------
# 6) ANALYSIS C: COHORT — deployment-year FE (fixest) and
#    within-cohort age slopes
# ------------------------------------------------------------

say("Analysis C: deployment-cohort specifications")

# 6a) age slope net of cohort fixed effects (assumes no period effect)
fe_fit <- function(dd, yvar) {
  m <- tryCatch(
    feols(as.formula(paste0("`", yvar, "` ~ op_years | deploy_year")),
          data = dd, weights = ~n_obs, cluster = ~station),
    error = function(e) NULL)
  if (is.null(m)) return(data.table())
  ct <- coeftable(m)
  data.table(
    group = "Global", spec = "Age + cohort FE",
    outcome = yvar, effect_type = "Beta per operational year",
    beta_age = ct["op_years", "Estimate"], se_age = ct["op_years", "Std. Error"],
    p_age = ct["op_years", 4],
    ci_age_lwr = ct["op_years", "Estimate"] - 1.96 * ct["op_years", "Std. Error"],
    ci_age_upr = ct["op_years", "Estimate"] + 1.96 * ct["op_years", "Std. Error"],
    n_rows = m$nobs, n_sensors = uniqueN(dd$station)
  )
}

res_cohort_fe <- rbindlist(list(
  fe_fit(sens_bin, "mean_abs_s"),
  fe_fit(sens_bin, "mean_pct_s")
), fill = TRUE)
res_cohort_fe[outcome == "mean_abs_s", outcome := "Absolute channel divergence |A-B|"]
res_cohort_fe[outcome == "mean_pct_s", outcome := "Relative channel divergence"]

# 6b) age slope within each deployment cohort (descriptive)
res_cohort_within <- rbindlist(lapply(sort(unique(sens_bin$deploy_year)), function(cy) {
  dd <- sens_bin[deploy_year == cy]
  if (uniqueN(dd$station) < MIN_SENSORS_YEAR) return(data.table())
  fit_lm_row(dd, "mean_abs_s", "op_years",
             paste0("Within cohort ", cy),
             "Absolute channel divergence |A-B|", paste0("Cohort ", cy))
}), fill = TRUE)

# ------------------------------------------------------------
# 7) GEOS-CF DIVERGENCE: age + period + month-of-year FE
#    (the India concern — ambient/model secular trends)
# ------------------------------------------------------------

res_md <- data.table()
if (has_geos) {
  say("GEOS-CF divergence: age + period with month FE")
  md_fit <- function(dd, group) {
    dd <- dd[is.finite(mean_md_s)]
    if (nrow(dd) < min_rows_for_trend || uniqueN(dd$station) < min_stations_group)
      return(data.table())
    m <- tryCatch(
      feols(mean_md_s ~ op_years + cal_time | cal_month,
            data = dd, weights = ~n_obs, cluster = ~station),
      error = function(e) NULL)
    if (is.null(m)) return(data.table())
    ct <- coeftable(m)
    data.table(
      group = group, spec = "Age + period + month FE",
      outcome = "PurpleAir - GEOS-CF divergence",
      effect_type = "Beta per operational year",
      beta_age = ct["op_years", "Estimate"], se_age = ct["op_years", "Std. Error"],
      p_age = ct["op_years", 4],
      ci_age_lwr = ct["op_years", "Estimate"] - 1.96 * ct["op_years", "Std. Error"],
      ci_age_upr = ct["op_years", "Estimate"] + 1.96 * ct["op_years", "Std. Error"],
      beta_period = ct["cal_time", "Estimate"], se_period = ct["cal_time", "Std. Error"],
      n_rows = m$nobs, n_sensors = uniqueN(dd$station)
    )
  }
  res_md <- rbindlist(c(
    list(md_fit(sens_bin, "Global")),
    lapply(sort(unique(sens_bin$region)), function(r)
      md_fit(sens_bin[region == r], r))
  ), fill = TRUE)
}

# ------------------------------------------------------------
# 8) SAVE OUTPUTS
# ------------------------------------------------------------

all_res <- rbindlist(list(res_within, res_joint, res_cohort_fe,
                          res_cohort_within, res_md), fill = TRUE)
fwrite(all_res, file.path(outdir, "age_period_cohort_results.csv"))
say("Wrote: age_period_cohort_results.csv (", nrow(all_res), " rows)")

# Forest-style figure: |A-B| age slope across specifications
plot_dt <- all_res[outcome == "Absolute channel divergence |A-B|" &
                     grepl("^Within|Age \\+", spec) & group %in%
                     c("Global", paste0("Year ", 2019:2025),
                       paste0("Cohort ", 2018:2024))]
plot_dt[, label := fifelse(grepl("^Within [0-9]", spec), group,
                           paste0(group, ": ", spec))]
plot_dt[, panel := fifelse(grepl("^Within [0-9]", spec), "Within calendar year",
                    fifelse(grepl("cohort", spec, ignore.case = TRUE),
                            "Within deployment cohort / cohort FE", "Joint age + period"))]

p <- ggplot(plot_dt, aes(x = beta_age, y = label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = ci_age_lwr, xmax = ci_age_upr), height = 0.2) +
  facet_grid(panel ~ ., scales = "free_y", space = "free_y") +
  labs(x = expression("Age slope for |A-B| (" * mu * "g/" * m^3 * " per operational year)"),
       y = NULL,
       title = "Operational-age effect on |A-B| holding calendar time fixed",
       subtitle = "Bars are 95% CIs with sensor-clustered robust SEs") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))

ggsave(file.path(outdir, "age_slope_by_specification.jpeg"), p,
       width = 8, height = 7, dpi = 300)
say("Wrote: age_slope_by_specification.jpeg")

# Word table
ft_dt <- all_res[, .(
  Group = group, Specification = spec, Outcome = outcome,
  `Age effect (95% CI)` = sprintf("%.3f (%.3f, %.3f)", beta_age, ci_age_lwr, ci_age_upr),
  `Period effect (95% CI)` = fifelse(
    is.na(beta_period), "-",
    sprintf("%.3f (%.3f, %.3f)", beta_period, ci_period_lwr, ci_period_upr)),
  `N sensors` = format(n_sensors, big.mark = ",")
)]
ft <- flextable(ft_dt)
ft <- autofit(fontsize(ft, size = 8, part = "all"))
doc <- read_docx()
doc <- body_add_par(doc,
  "Table SX: Operational-age effects on cf_1 channel divergence under specifications separating operational age from calendar time and deployment cohort.",
  style = "Normal")
doc <- body_add_flextable(doc, ft)
print(doc, target = file.path(outdir, "age_period_cohort_table.docx"))
say("Wrote: age_period_cohort_table.docx")

# ------------------------------------------------------------
# Suggested SI text (edit to match your results):
#
# "Because operational age is partially confounded with calendar time
#  and deployment cohort, we re-estimated age trends (i) within each
#  calendar year, (ii) jointly with continuous calendar time, and (iii)
#  within deployment-year cohorts and with cohort fixed effects. Age =
#  period - cohort, so no single model identifies all three effects;
#  instead we verify that the operational-age coefficient is stable
#  across specifications (Table SX; Figure SX). Age-related increases in
#  channel divergence persisted within fixed calendar years and after
#  adjustment for calendar time or cohort, indicating that degradation
#  trajectories are not an artifact of secular trends in ambient
#  concentrations, network composition, or model performance. For
#  PurpleAir-GEOS-CF divergence, adjustment for calendar time and
#  season [reduced/did not materially change] the estimated age trend
#  [notably in India], consistent with [interpretation]."
# ============================================================
