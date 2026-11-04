# Top-detection-day sensitivity check for every final model: refits the same
# final formula (same family, same AR1 procedure via fit_ar1_bam) with the
# 1 and 3 days of highest response value removed, and records whether each
# retained term stays significant and how deviance explained changes. A fit
# whose terms or deviance explained collapse when one to three days are
# removed is resting on those days, not on a general relationship. Used for
# the main-text / supplement / table-only tiering of models and the exploratory flag in
# 04_build_summary_tables.R.
#
# Not applied to the occurrence (presence/absence) models: with a 0/1
# response there are no "top" days to remove.
#
# Also reports detection days within the final model's own rows (matching
# the "Non-zero days" column of Table S-B) and the share of all detections
# falling on the top 3 days. Writes new_models/top_day_sensitivity.csv.

source("config.R")
source("00a_load_data.R")
source("00b_load_dvm.R")
source("00d_gam_pipeline_functions.R")
source("00f_model_specs.R")


daily <- add_response_columns(merge_dvm_into_daily(load_cce1_daily(), load_dvm_amplitude(validate = FALSE)))
specs <- cetacean_specs()

labels <- c(paste0("NASC", c(5, 55, 155, 255, 355, 455)), paste0("DVM", DVM_DEPTHS), specs$label)
spec_of <- function(lab) specs[specs$label == lab, ]

family_for <- function(lab, res) {
  if (grepl("^(NASC|DVM)", lab)) return(gaussian())
  f <- spec_family(spec_of(lab)$family)
  if (!is.null(f)) f else switch(res$family, "Tweedie" = tw(), "negative binomial" = nb(), gaussian())
}

refit_without_top <- function(res, lab, k_drop) {
  f <- formula(res$final$model)
  d <- daily[stats::complete.cases(daily[, unique(c(all.vars(f), "Date")), drop = FALSE]), ]
  if (k_drop >= sum(d[[res$y]] != 0)) return(NULL)   # nothing would be left to model
  d <- d[-order(d[[res$y]], decreasing = TRUE)[seq_len(k_drop)], ]
  fit <- fit_ar1_bam(f, d, family_for(lab, res))
  if (is.null(fit)) return(NULL)
  list(p = summary(fit$model)$s.table[, "p-value"], dev = round(model_dev_expl(fit$model) * 100, 1))
}

rows <- list()
for (lab in labels) {
  res <- readRDS(file.path(MODEL_DIR, paste0(lab, "_model.rds")))
  is_cetacean <- lab %in% specs$label
  part <- if (is_cetacean) spec_of(lab)$part else "single"
  y <- if (!is.null(res$final)) res$final$model$y else daily[[res$y]][!is.na(daily[[res$y]])]
  # detection days and top-3 share, on the original scale of the species
  orig <- if (is_cetacean) daily[[spec_of(lab)$species]] else NULL
  orig_pos <- if (is_cetacean) orig[!is.na(orig) & orig > 0] else NULL
  base <- data.frame(
    response = lab, part = part, family = res$family, n = length(y),
    n_nonzero = if (!is_cetacean) NA else if (part == "occurrence") sum(y == 1)
                else if (part == "when present") length(y) else sum(y > 0),
    pct_nonzero = if (!is_cetacean) NA else if (part == "occurrence") round(100 * mean(y == 1), 1)
                  else if (part == "when present") 100 else round(100 * mean(y > 0), 1),
    top3_share_pct = if (!is_cetacean || part == "occurrence") NA
                     else round(100 * sum(sort(orig_pos, decreasing = TRUE)[1:3]) / sum(orig_pos), 1))
  if (is.null(res$final)) {
    rows[[lab]] <- cbind(base, term = "(none retained)", edf = NA, p = NA, shrinkage_p = NA,
                         dev_expl = NA, p_drop1 = NA, p_drop3 = NA,
                         dev_expl_drop1 = NA, dev_expl_drop3 = NA)
    next
  }
  st <- summary(res$final$model)$s.table
  ss <- summary(res$shrinkage)$s.table
  s1 <- if (part == "occurrence") NULL else refit_without_top(res, lab, 1)
  s3 <- if (part == "occurrence") NULL else refit_without_top(res, lab, 3)
  rows[[lab]] <- cbind(base,
    term = sub("^s\\((.*)\\)$", "\\1", rownames(st)),
    edf = round(st[, "edf"], 2), p = signif(st[, "p-value"], 2),
    shrinkage_p = signif(ss[match(rownames(st), rownames(ss)), "p-value"], 2),
    dev_expl = round(model_dev_expl(res$final$model) * 100, 1),
    p_drop1 = if (is.null(s1)) NA else signif(s1$p, 2),
    p_drop3 = if (is.null(s3)) NA else signif(s3$p, 2),
    dev_expl_drop1 = if (is.null(s1)) NA else s1$dev,
    dev_expl_drop3 = if (is.null(s3)) NA else s3$dev)
}
out <- do.call(rbind, rows); rownames(out) <- NULL
write.csv(out, file.path(MODEL_DIR, "top_day_sensitivity.csv"), row.names = FALSE)
print(out)
