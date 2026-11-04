# Out-of-sample check for every final model: leave-one-calendar-month-out
# cross-validation. Each month in turn is removed, the model is refit on the
# remaining months (same family and AR1 procedure as the pipeline, via
# fit_ar1_bam(), which uses gamm() for occurrence models), and the left-out
# month is predicted. Whole months are left out, not single days, because
# neighbouring days are too similar (strong AR1) to count as new data.
#
# Prediction error is the deviance of the left-out observations under the
# trained model. "CV deviance explained" compares it with an intercept-only
# model trained the same way: 1 - CV deviance(model) / CV deviance(intercept).
# Values near 0 mean the model predicts unseen months no better than the
# mean; negative values mean it predicts them worse (overfitting, or
# extrapolation to conditions not seen in the training months).
#
# Besides the final model, every nested model along the forward-selection
# order of the retained terms is cross-validated, to show how many terms
# still improve out-of-sample prediction. Note this is a harsh test for one
# year of data: leaving out a month also removes that month's particular
# combination of conditions, so predictions often extrapolate.
#
# Writes new_models/cross_validation.csv (one row per model x number of
# terms) and prints a one-row-per-model summary.

source("config.R")
source("00a_load_data.R")
source("00b_load_dvm.R")
source("00d_gam_pipeline_functions.R")
source("00f_model_specs.R")


daily <- add_response_columns(merge_dvm_into_daily(load_cce1_daily(), load_dvm_amplitude(validate = FALSE)))
specs <- cetacean_specs()
labels <- c(paste0("NASC", c(5, 55, 155, 255, 355, 455)), paste0("DVM", DVM_DEPTHS), specs$label)

family_for <- function(lab, res) {
  if (grepl("^(NASC|DVM)", lab)) return(gaussian())
  f <- spec_family(specs$family[specs$label == lab])
  if (!is.null(f)) f else switch(res$family, "Tweedie" = tw(), "negative binomial" = nb(), gaussian())
}

#' deviance of observations y at predicted mu. All models of one response are
#' scored with the SAME family object (eval_fam: the final model's family,
#' fitted to all data): negative binomial and Tweedie deviance depend on the
#' dispersion/power parameter, so letting each refit use its own estimate
#' would put the models on different scales and make them incomparable.
pred_deviance <- function(eval_fam, y, mu) {
  if (eval_fam$family == "binomial") mu <- pmin(pmax(mu, 1e-6), 1 - 1e-6)
  sum(eval_fam$dev.resids(y, mu, rep(1, length(y))))
}

cv_one <- function(y_name, tt, d, fam, month, eval_fam) {
  total <- 0
  for (mo in unique(month)) {
    tr <- d[month != mo, ]; te <- d[month == mo, ]
    if (length(tt) == 0) {
      fit <- tryCatch(suppressWarnings(gam(as.formula(paste(y_name, "~ 1")), data = tr, family = fam)),
                      error = function(e) NULL)
      model <- fit
    } else {
      f <- as.formula(paste(y_name, "~", paste(sprintf("s(%s, k=4, bs='cs')", tt), collapse = " + ")))
      fit <- tryCatch(fit_ar1_bam(f, tr, fam), error = function(e) NULL)
      model <- if (is.null(fit)) NULL else fit$model
    }
    if (is.null(model)) return(NA_real_)
    mu <- as.numeric(predict(model, newdata = te, type = "response"))
    total <- total + pred_deviance(eval_fam, te[[y_name]], mu)
  }
  total
}

rows <- list()
for (lab in labels) {
  res <- readRDS(file.path(MODEL_DIR, paste0(lab, "_model.rds")))
  if (is.null(res$final)) {
    rows[[lab]] <- data.frame(response = lab, n_terms = 0, terms = "(no term retained)",
                              insample_dev_expl = NA, cv_deviance = NA, cv_dev_expl = NA)
    next
  }
  fam <- family_for(lab, res)
  final_terms <- res$selected
  ordered <- res$selected_forward[res$selected_forward %in% final_terms]   # forward order
  d <- daily[stats::complete.cases(daily[, c(res$y, final_terms, "Date")]), ]
  d <- d[order(d$Date), ]
  month <- format(d$Date, "%Y-%m")
  t0 <- Sys.time()
  eval_fam <- res$final$model$family
  base <- cv_one(res$y, character(0), d, fam, month, eval_fam)
  for (k in seq_along(ordered)) {
    tt <- ordered[seq_len(k)]
    cvd <- cv_one(res$y, tt, d, fam, month, eval_fam)
    rows[[paste(lab, k)]] <- data.frame(
      response = lab, n_terms = k, terms = paste(tt, collapse = " + "),
      insample_dev_expl = if (k == length(ordered)) round(100 * model_dev_expl(res$final$model), 1) else NA,
      cv_deviance = round(cvd, 2), cv_dev_expl = round(100 * (1 - cvd / base), 1))
  }
  message(sprintf("%-11s done [%.1f min]", lab, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
out <- do.call(rbind, rows); rownames(out) <- NULL
write.csv(out, file.path(MODEL_DIR, "cross_validation.csv"), row.names = FALSE)

# one row per model: final model vs the best nested model
summ <- do.call(rbind, lapply(split(out, out$response), function(x) {
  fin <- x[nrow(x), ]
  best <- if (all(is.na(x$cv_dev_expl))) fin else x[which.max(x$cv_dev_expl), ]
  data.frame(response = fin$response, final_terms = fin$n_terms,
             insample_dev_expl = fin$insample_dev_expl, cv_dev_expl_final = fin$cv_dev_expl,
             best_n_terms = best$n_terms, cv_dev_expl_best = best$cv_dev_expl)
}))
summ <- summ[match(labels, summ$response), ]
print(summ, row.names = FALSE)
