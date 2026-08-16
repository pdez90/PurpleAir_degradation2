# ============================================================
# 07_pa_epa_calibration_cv.R
#
# SENSOR-GROUPED CROSS-VALIDATION OF THE PA-EPA CALIBRATION MODELS
# (review flag 3: the paper's RMSE 17.2 -> 6.7 and age-aware r = 0.83
#  are in-sample; a reviewer will ask for out-of-sample numbers)
#
# 5-fold CV with folds assigned BY SENSOR (whole sensors held out), so
# the metric reflects transfer to unseen sensors — the relevant unit,
# since hours within a sensor are strongly correlated. Replicates the
# paper's model specifications exactly:
#   c1  lm:  EPA ~ PA
#   c2  lm:  EPA ~ PA + PA^2
#   c3  lm:  EPA ~ ns(PA, 4)
#   c4  bam: EPA ~ s(PA, bs="ps")                (fREML, discrete)
#   add bam: EPA ~ s(PA) + s(age/1000)           (age-aware additive)
#   te  bam: EPA ~ te(PA, age/1000)              (age-aware tensor)
# Training rows per fold capped at 2.5M (sampled, seed 123) to mirror
# the paper's MAX_MODEL_ROWS practice; evaluation uses ALL held-out rows.
#
# OUTPUT: pa_epa_cv/ CSV + Word table: in-sample vs cross-validated
#         RMSE and r per model.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(mgcv)
  library(splines)
  library(flextable)
  library(officer)
})

join_rdata <- "/Users/priyanka/Downloads/PurpleAir_EPA/pa_epa_100m_join/pa_epa_joined_within_100m_2018_2025.RData"
outdir <- "/Users/priyanka/Downloads/PurpleAir_EPA/pa_epa_cv"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

K_FOLDS <- 5L
MAX_FIT_ROWS <- 2500000L
say <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")),
      paste0(..., collapse = ""), "\n")
}

say("Loading joined PA-EPA data")
e <- new.env(); load(join_rdata, envir = e)
objs <- ls(e)
joined <- NULL
if ("joined" %in% objs) {
  joined <- get("joined", envir = e)
} else {
  for (o in objs) {
    x <- get(o, envir = e)
    if (is.data.frame(x) && all(c("pa_pm25_cf1", "epa_pm25") %in% names(x))) {
      joined <- x; say("Using object: ", o); break
    }
  }
}
if (is.null(joined)) stop("No object with pa_pm25_cf1 + epa_pm25 found. Objects: ",
                          paste(objs, collapse = ", "))
setDT(joined); rm(e); gc()

joined <- joined[is.finite(pa_pm25_cf1) & is.finite(epa_pm25) & !is.na(station)]
joined[, station := as.character(station)]
if (!("hour" %in% names(joined))) {
  if ("op_hours_since_first_pa" %in% names(joined)) {
    joined[, hour := as.numeric(op_hours_since_first_pa)]
  } else {
    setorder(joined, station, datetime_utc)
    joined[, hour := as.numeric(difftime(datetime_utc, min(datetime_utc), units = "hours")),
           by = station]
  }
}
joined <- joined[is.finite(hour)]
joined[, hour_k := hour / 1000]
say("Rows: ", format(nrow(joined), big.mark = ","),
    " | sensors: ", format(uniqueN(joined$station), big.mark = ","))

# ------------------------------------------------------------
# fold assignment BY SENSOR
# ------------------------------------------------------------
set.seed(123)
stations <- joined[, unique(station)]
folds <- data.table(station = stations,
                    fold = sample(rep_len(1:K_FOLDS, length(stations))))
joined <- folds[joined, on = "station"]

rmse <- function(p, y) sqrt(mean((p - y)^2, na.rm = TRUE))
rr <- function(p, y) { ok <- is.finite(p) & is.finite(y); cor(p[ok], y[ok]) }

fit_predict <- function(train, test, spec) {
  if (nrow(train) > MAX_FIT_ROWS) {
    set.seed(123); train <- train[sample.int(.N, MAX_FIT_ROWS)]
  }
  fit <- switch(spec,
    c1 = lm(epa_pm25 ~ pa_pm25_cf1, data = train),
    c2 = lm(epa_pm25 ~ pa_pm25_cf1 + I(pa_pm25_cf1^2), data = train),
    c3 = lm(epa_pm25 ~ splines::ns(pa_pm25_cf1, df = 4), data = train),
    c4 = mgcv::bam(epa_pm25 ~ s(pa_pm25_cf1, bs = "ps"), data = train,
                   family = gaussian(), method = "fREML", discrete = TRUE),
    add = mgcv::bam(epa_pm25 ~ s(pa_pm25_cf1, bs = "ps") + s(hour_k, bs = "ps"),
                    data = train, family = gaussian(), method = "fREML", discrete = TRUE),
    te = mgcv::bam(epa_pm25 ~ te(pa_pm25_cf1, hour_k, bs = c("ps", "ps")),
                   data = train, family = gaussian(), method = "fREML", discrete = TRUE))
  as.numeric(predict(fit, newdata = test))
}

SPECS <- c(c1 = "Linear", c2 = "Quadratic", c3 = "Natural spline (df=4)",
           c4 = "GAM s(PA)", add = "Age-aware additive GAM",
           te = "Age-aware tensor GAM")

# ------------------------------------------------------------
# in-sample benchmark (paper's procedure) + raw
# ------------------------------------------------------------
say("In-sample benchmark fits")
insample <- rbindlist(lapply(names(SPECS), function(sp) {
  p <- fit_predict(joined, joined, sp)
  data.table(model = SPECS[[sp]], eval = "In-sample (paper)",
             rmse = rmse(p, joined$epa_pm25), r = rr(p, joined$epa_pm25))
}))
raw_row <- data.table(model = "Raw cf_1 (no calibration)", eval = "-",
                      rmse = rmse(joined$pa_pm25_cf1, joined$epa_pm25),
                      r = rr(joined$pa_pm25_cf1, joined$epa_pm25))

# ------------------------------------------------------------
# 5-fold sensor-grouped CV
# ------------------------------------------------------------
cv_preds <- list()
for (k in 1:K_FOLDS) {
  say("Fold ", k, " / ", K_FOLDS)
  train <- joined[fold != k]; test <- joined[fold == k]
  for (sp in names(SPECS)) {
    p <- tryCatch(fit_predict(train, test, sp), error = function(e) rep(NA_real_, nrow(test)))
    cv_preds[[paste(k, sp)]] <- data.table(spec = sp, y = test$epa_pm25, p = p)
  }
}
cvp <- rbindlist(cv_preds)
cv <- cvp[, .(rmse = rmse(p, y), r = rr(p, y)), by = spec]
cv[, model := SPECS[spec]]
cv[, eval := "5-fold CV, sensors held out"]

out <- rbind(raw_row, insample, cv[, .(model, eval, rmse, r)])
fwrite(out, file.path(outdir, "pa_epa_calibration_cv.csv"))
print(out[, .(model, eval, rmse = round(rmse, 2), r = round(r, 3))])

ft_dt <- dcast(out[eval != "-"], model ~ eval, value.var = c("rmse", "r"))
ft_dt_fmt <- ft_dt[, .(
  Model = model,
  `In-sample RMSE (µg/m³)` = sprintf("%.2f", `rmse_In-sample (paper)`),
  `CV RMSE (µg/m³)` = sprintf("%.2f", `rmse_5-fold CV, sensors held out`),
  `In-sample r` = sprintf("%.3f", `r_In-sample (paper)`),
  `CV r` = sprintf("%.3f", `r_5-fold CV, sensors held out`))]
ft <- autofit(fontsize(flextable(ft_dt_fmt), size = 9, part = "all"))
doc <- read_docx()
doc <- body_add_par(doc, paste0(
  "Table SX: In-sample versus sensor-grouped 5-fold cross-validated performance of ",
  "the PA-regulatory calibration models (", format(nrow(joined), big.mark = ","),
  " paired hourly observations; ", format(uniqueN(joined$station), big.mark = ","),
  " sensors; whole sensors held out per fold). Raw cf_1 RMSE = ",
  sprintf("%.1f", raw_row$rmse), " µg/m³, r = ", sprintf("%.2f", raw_row$r), "."),
  style = "Normal")
doc <- body_add_flextable(doc, ft)
print(doc, target = file.path(outdir, "pa_epa_calibration_cv.docx"))
say("Done. Outputs in: ", outdir)
