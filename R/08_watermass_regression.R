# Fig. 5: water-mass regressions with AR1 errors.
#
# (1) nautical area scattering strength S_A = 10 log10(sA) (dB re 1 m^2 nmi^-2;
#     sA = NASC) in each of the six acoustic depth layers against the PSUW,
#     PEW and ENPCW fractions at the closest of the three mooring OMP depths:
#       0-50 m, 50-100 m            -> 75 m
#       150-200 m                   -> 150 m
#       250-300, 350-400, 450-500 m -> 301 m
# (2) the mesopelagic scattering layer (MSL) parameters from the edge-tracing
#     detector (MATLAB seed-and-walk detector, matlab/run_MSL_seedwalk.m; its per-night output MSL_night_seedWalk.csv is part of the data deposit) --
#     nightly top edge depth, core depth and S_A of the integrated NASC below the
#     edge -- against the three fractions at the deepest OMP depth, 301 m.
#     The raw nightly medians are used (not the persistence-filtered or
#     gap-filled columns), so AR1 is not fitted to an already-smoothed series.
#
# Each regression is fit by OLS and by GLS with AR1 residuals
# (nlme::gls, corAR1 indexed by day number, so missing days are spaced
# correctly; ML). Mirroring the GAM pruning rule, a regression is retained
# only if BOTH the OLS p and the AR1 p are < 0.05. Fig. 5 shows every retained
# regression, with the AR1 (GLS) fit and its 95% confidence band.
#
# Outputs (results/fig5_watermass_ar1/):
#   fig5_watermass_ar1_regressions.csv  all 27 regressions + retained flag
#   Fig5_watermass_AR1.pdf              retained regressions only
#   TableS_watermass_AR1.docx / .pdf    Table S-C, all 27 regressions

library(nlme)
source("config.R")
source("00a_load_data.R")
source("00e_word_tables.R")

OUT_DIR <- file.path(RESULTS_DIR, "fig5_watermass_ar1")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

ALPHA <- 0.05

# ---- predictors: water-mass fractions at the three OMP depths ---------------
omp_files <- c(`75` = "x_error_depth75_id1004.csv",
               `150` = "x_error_depth150_id1003.csv",
               `301` = "x_error_depth301_id1002.csv")
water_masses <- c(psuw = "PSUW", pew = "PEW", enpcw = "ENPCW")

load_fractions <- function(omp_depth) {
  d <- read.csv(file.path(OMP_DIR, omp_files[[omp_depth]]))
  out <- data.frame(Date = as.Date(sprintf("%04d-%02d-%02d", d$year, d$month, d$day)))
  for (wm in names(water_masses)) out[[wm]] <- d[[grep(paste0("^", wm, "_"), names(d))[1]]]
  out
}
fractions <- setNames(lapply(names(omp_files), load_fractions), names(omp_files))

# ---- responses ----------------------------------------------------------------
nasc_layers <- data.frame(
  code  = paste0("NASC", c(5, 55, 155, 255, 355, 455)),
  label = c("0–50 m", "50–100 m", "150–200 m", "250–300 m", "350–400 m", "450–500 m"),
  omp   = c("75", "75", "150", "301", "301", "301"),
  stringsAsFactors = FALSE)

daily <- load_cce1_daily()

# MSL nightly series, dated by the night's start. Dates parsed field by field:
# as.Date(format = "%d-%b-%Y") is locale-dependent for "%b".
msl <- read.csv(MSL_CSV, stringsAsFactors = FALSE)
month_num <- c(Jan = 1, Feb = 2, Mar = 3, Apr = 4, May = 5, Jun = 6,
               Jul = 7, Aug = 8, Sep = 9, Oct = 10, Nov = 11, Dec = 12)
parts <- do.call(rbind, strsplit(substr(msl$nightStart, 1, 11), "-"))
msl$Date <- as.Date(sprintf("%04d-%02d-%02d", as.integer(parts[, 3]),
                            month_num[parts[, 2]], as.integer(parts[, 1])))
stopifnot(!anyNA(msl$Date), !any(duplicated(msl$Date)))
msl_params <- data.frame(
  code  = c("medianEdgeM", "medianCoreDepthM", "medianIntegratedNASC"),
  label = c("MSL top edge depth", "MSL core depth", "MSL integrated sA"),
  ylab  = c("Top edge depth (m)", "Core depth (m)", "SA"),
  depth_axis = c(TRUE, TRUE, FALSE),
  to_dB = c(FALSE, FALSE, TRUE),
  stringsAsFactors = FALSE)

# S_A = 10 log10(sA), dB re 1 m^2 nmi^-2; sA <= 0 gives -Inf and is dropped in fit_pair()
to_SA <- function(sA) 10 * log10(sA)
SA_ylab <- expression(italic(S)[A] ~ "(dB re 1 m"^2 ~ "nmi"^-2 * ")")

# ---- fitting --------------------------------------------------------------------
fit_pair <- function(date, y, x) {
  dd <- data.frame(Date = date, y = y, x = x)
  dd <- dd[is.finite(dd$y) & is.finite(dd$x), ]   # drops NASC = 0 days and undetected nights
  dd <- dd[order(dd$Date), ]
  dd$t <- as.numeric(dd$Date)
  ols <- summary(lm(y ~ x, data = dd))
  g <- tryCatch(gls(y ~ x, data = dd, correlation = corAR1(form = ~ t), method = "ML"),
                error = function(e) NULL)
  list(data = dd, gls = g, n = nrow(dd), r2 = ols$r.squared,
       p_ols = ols$coefficients["x", "Pr(>|t|)"],
       slope = if (is.null(g)) NA_real_ else coef(g)[["x"]],
       p_ar1 = if (is.null(g)) NA_real_ else summary(g)$tTable["x", "p-value"],
       phi = if (is.null(g)) NA_real_
             else unname(coef(g$modelStruct$corStruct, unconstrained = FALSE)))
}

fits <- list()
for (i in seq_len(nrow(nasc_layers))) for (wm in names(water_masses)) {
  L <- nasc_layers[i, ]
  dd <- merge(daily[, c("Date", L$code)], fractions[[L$omp]][, c("Date", wm)], by = "Date")
  fits[[length(fits) + 1]] <- c(
    list(group = "Backscatter", response = L$label, ylab = "SA", depth_axis = FALSE,
         water_mass = water_masses[[wm]], omp_depth = L$omp),
    fit_pair(dd$Date, to_SA(dd[[L$code]]), dd[[wm]]))
}
for (i in seq_len(nrow(msl_params))) for (wm in names(water_masses)) {
  P <- msl_params[i, ]
  dd <- merge(msl[, c("Date", P$code)], fractions[["301"]][, c("Date", wm)], by = "Date")
  y <- if (P$to_dB) to_SA(dd[[P$code]]) else dd[[P$code]]
  fits[[length(fits) + 1]] <- c(
    list(group = "MSL", response = P$label, ylab = P$ylab, depth_axis = P$depth_axis,
         water_mass = water_masses[[wm]], omp_depth = "301"),
    fit_pair(dd$Date, y, dd[[wm]]))
}

retained <- vapply(fits, function(f) isTRUE(f$p_ols < ALPHA & f$p_ar1 < ALPHA), logical(1))

# Leave-one-calendar-month-out cross-validation, as for the GAMs
# (09_cross_validation.R): refit the AR1 regression without each month,
# predict that month from the regression line, and compare the squared
# prediction error with that of the training-month mean.
# CV R2 = 1 - SSE(model) / SSE(mean); near or below 0 = no out-of-sample skill.
cv_r2 <- function(dd) {
  month <- format(dd$Date, "%Y-%m")
  sse <- 0; sse0 <- 0
  for (mo in unique(month)) {
    tr <- dd[month != mo, ]; te <- dd[month == mo, ]
    g <- tryCatch(gls(y ~ x, data = tr, correlation = corAR1(form = ~ t), method = "ML"),
                  error = function(e) NULL)
    b <- if (is.null(g)) coef(lm(y ~ x, data = tr)) else coef(g)
    sse <- sse + sum((te$y - (b[1] + b[2] * te$x))^2)
    sse0 <- sse0 + sum((te$y - mean(tr$y))^2)
  }
  1 - sse / sse0
}
for (k in seq_along(fits)) fits[[k]]$cv_r2 <- cv_r2(fits[[k]]$data)

tbl <- do.call(rbind, lapply(fits, function(f) data.frame(
  group = f$group, response = f$response, water_mass = f$water_mass,
  omp_depth_m = f$omp_depth, n = f$n, R2_ols = round(f$r2, 3),
  p_ols = signif(f$p_ols, 3), slope_ar1 = signif(f$slope, 4),
  p_ar1 = signif(f$p_ar1, 3), ar1_phi = round(f$phi, 3), cv_R2 = round(f$cv_r2, 3),
  stringsAsFactors = FALSE)))
tbl$retained <- retained
write.csv(tbl, file.path(OUT_DIR, "fig5_watermass_ar1_regressions.csv"), row.names = FALSE)
print(tbl, row.names = FALSE)

# ---- Fig. 5: retained regressions only -----------------------------------------
fmt_p <- function(p) if (p < 0.0001) "< 0.0001" else paste("=", signif(p, 2))

panel <- function(f, tag) {
  dd <- f$data
  newx <- seq(min(dd$x), max(dd$x), length.out = 100)
  X <- cbind(1, newx)
  fit <- drop(X %*% coef(f$gls))
  se <- sqrt(rowSums((X %*% vcov(f$gls)) * X))
  ylim <- range(dd$y)
  if (f$depth_axis) ylim <- rev(ylim)   # depth increases downward
  plot(dd$x, dd$y, pch = 20, col = "lightblue", ylim = ylim,
       xlab = sprintf("%s fraction at %s m", f$water_mass, f$omp_depth),
       ylab = if (identical(f$ylab, "SA")) SA_ylab else f$ylab,
       main = sprintf("%s × %s", f$response, f$water_mass), cex.main = 0.95)
  polygon(c(newx, rev(newx)), c(fit + 1.96 * se, rev(fit - 1.96 * se)),
          col = adjustcolor("red", 0.15), border = NA)
  lines(newx, fit, lwd = 2)
  legend(if (f$depth_axis) "bottomleft" else "topleft", bty = "n", cex = 0.75,
         legend = c(sprintf("p(OLS) %s, R² = %.2f", fmt_p(f$p_ols), f$r2),
                    sprintf("p(AR1) %s", fmt_p(f$p_ar1))))
  mtext(sprintf("(%s)", tag), side = 3, adj = -0.18, line = 0.9, font = 2)
}

kept <- fits[retained]
ncol <- 2
nrow <- ceiling(length(kept) / ncol)
pdf(file.path(OUT_DIR, "Fig5_watermass_AR1.pdf"), width = 7.5, height = 0.3 + 3 * nrow)
par(mfrow = c(nrow, ncol), mar = c(4, 4.2, 2.5, 1))
for (k in seq_along(kept)) panel(kept[[k]], letters[k])
dev.off()

cat(sprintf("\n%d of %d regressions retained (OLS p and AR1 p both < %.2f)\n",
            sum(retained), length(retained), ALPHA))
cat("Saved:", file.path(OUT_DIR, c("fig5_watermass_ar1_regressions.csv", "Fig5_watermass_AR1.pdf")),
    sep = "\n  ")

# ---- Table S-C: all 27 regressions (Word + PDF) --------------------------------
# Same styling as Tables S-A/S-B (00e_word_tables.R). One block of three rows
# (PSUW, PEW, ENPCW) per response; retained regressions in bold, with their
# Fig. 5 panel letter.
fig5_panel <- character(length(fits))
fig5_panel[retained] <- sprintf("(%s)", letters[seq_len(sum(retained))])

fmt_p_tab <- function(p) if (is.na(p)) "" else if (p < 0.001) "<0.001" else sprintf("%.3f", p)
fmt_slope <- function(x) if (is.na(x)) "" else format(signif(x, 3), scientific = FALSE)

table_rows <- function(idx) {
  f <- fits[idx]
  resp <- vapply(f, `[[`, character(1), "response")
  data.frame(
    group = ifelse(c(TRUE, resp[-1] != resp[-length(resp)]), resp, ""),
    water_mass = vapply(f, `[[`, character(1), "water_mass"),
    omp = vapply(f, `[[`, character(1), "omp_depth"),
    n = vapply(f, function(x) as.character(x$n), character(1)),
    r2 = vapply(f, function(x) sprintf("%.3f", x$r2), character(1)),
    p_ols = vapply(f, function(x) fmt_p_tab(x$p_ols), character(1)),
    slope = vapply(f, function(x) fmt_slope(x$slope), character(1)),
    p_ar1 = vapply(f, function(x) fmt_p_tab(x$p_ar1), character(1)),
    phi = vapply(f, function(x) sprintf("%.2f", x$phi), character(1)),
    cv = vapply(f, function(x) sprintf("%.3f", x$cv_r2), character(1)),
    panel = fig5_panel[idx],
    stringsAsFactors = FALSE)
}

tab_labels <- c(group = "Depth layer / MSL parameter", water_mass = "Water mass",
                omp = "OMP depth (m)", n = "n", r2 = "R² (OLS)", p_ols = "OLS p",
                slope = "AR1 slope", p_ar1 = "AR1 p", phi = "φ (AR1)",
                cv = "CV R²", panel = "Fig. 5 panel")
tab_width <- c(group = 1.9, water_mass = 0.8, omp = 0.8, n = 0.45, r2 = 0.65, p_ols = 0.6,
               slope = 0.75, p_ar1 = 0.6, phi = 0.6, cv = 0.6, panel = 0.65)

table_block <- function(idx, section_label) {
  ft <- booktabs_ft(table_rows(idx), tab_labels, tab_width,
                    left_cols = c("group", "water_mass"), section_label = section_label)
  bold_rows <- which(retained[idx]) + 1   # +1 for the section row
  if (length(bold_rows)) ft <- bold(ft, i = bold_rows, part = "body")
  ft
}

is_msl <- vapply(fits, function(f) f$group == "MSL", logical(1))
tab_note <- paste(
  "Each row is one regression of the response on one water-mass fraction at the stated OMP depth:",
  "the closest of 75, 150 and 301 m for each backscatter layer, and 301 m for the MSL parameters.",
  "R² and OLS p are from ordinary least squares; AR1 slope, AR1 p and φ are from generalized",
  "least squares with first-order autoregressive residuals (nlme; AR1 indexed by day; ML).",
  "Rows in bold have both OLS p and AR1 p < 0.05 and are shown in Fig. 5 (panel given).",
  "CV R² is the out-of-sample R² when each calendar month is predicted by the AR1 regression",
  "refit without it (leave-one-month-out cross-validation), relative to the training-month mean;",
  "values near or below 0 mean no skill in predicting unseen months.",
  "Slopes are per unit water-mass fraction, in dB re 1 m² nmi⁻² (SA = 10 log10 sA) for backscatter and",
  "integrated backscatter, and in m for MSL edge and core depth (negative = shallower). PEW at 75 m is non-zero on only",
  "1 of 390 days, so the PEW regressions for 0–50 and 50–100 m are not meaningful.",
  "MSL, mesopelagic scattering layer; n, days (backscatter) or nights with a detected layer (MSL);",
  "PSUW, Pacific Subarctic Upper Water; PEW, Pacific Equatorial Water; ENPCW, Eastern North",
  "Pacific Central Water.")

doc <- landscape_doc()
doc <- add_caption(doc, "Table S-C.",
  paste("Regressions of acoustic backscatter in six depth layers, and of the mesopelagic scattering",
        "layer's top edge depth, core depth and integrated backscatter, on water-mass fractions,",
        "fit with and without first-order autoregressive (AR1) errors."))
doc <- body_add_flextable(doc, table_block(which(!is_msl), "Backscatter, SA (dB re 1 m² nmi⁻²)"), align = "left")
doc <- add_gap(doc)
doc <- body_add_flextable(doc, table_block(which(is_msl), "Mesopelagic scattering layer (MSL)"),
                          align = "left")
doc <- add_note(doc, tab_note)
save_docx_and_pdf(doc, file.path(OUT_DIR, "TableS_watermass_AR1.docx"))
