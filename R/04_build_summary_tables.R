# Builds the supplementary GAM summary tables from the final <response>_model.rds
# files, as term-level rows (one row per selected smooth, with edf/F/p from the
# final AR1 model and the matching term's p-value from the shrinkage GAM):
#   Table S-A  NASC backscatter (upper block) + DVM amplitude (lower block)
#   Table S-B  cetacean/odontocete responses (Supporting Information Tables S7, S8)
#
# Table S-B also reports the number (and %) of non-zero detection days per
# response. Models with fewer than EXPLORATORY_MIN_NONZERO non-zero days are
# marked as exploratory (dagger on the taxon name + footnote): their fits can
# hinge on a handful of detection days (see 03_top_day_sensitivity.R), and
# deviance explained overstates predictive skill in that situation. The
# negative-CV and top-day criteria for the dagger in the Supporting Information
# were applied by hand from cross_validation.csv and top_day_sensitivity.csv.

library(mgcv)

source("config.R")
source("00a_load_data.R")
source("00d_gam_pipeline_functions.R")
source("00e_word_tables.R")
source("00f_model_specs.R")


EXPLORATORY_MIN_NONZERO <- 25   # fewer detection days than this -> exploratory

# out-of-sample (leave-one-month-out) deviance explained of each final model,
# from 09_cross_validation.R (run that first); blank if not available
cv_file <- file.path(MODEL_DIR, "cross_validation.csv")
cv_final <- if (file.exists(cv_file)) {
  cv <- read.csv(cv_file)
  cv <- cv[cv$n_terms > 0, ]
  last <- tapply(seq_len(nrow(cv)), cv$response, max)
  setNames(cv$cv_dev_expl[last], names(last))
} else {
  message("cross_validation.csv not found: CV column left blank (run 09_cross_validation.R)")
  numeric(0)
}
cv_for <- function(label) if (label %in% names(cv_final)) cv_final[[label]] else NA_real_
DAGGER <- "†"

# full daily record, used only for the detection-day count of responses with
# no retained term (fitted models count from their own model frame instead)
daily <- add_response_columns(load_cce1_daily())
specs <- cetacean_specs()
part_of <- function(label) if (label %in% specs$label) specs$part[specs$label == label] else "single"

nasc_depth_label <- function(code) {
  # codes are nominal (NASC5, NASC55, ...); the layers span true depths
  # 0-50, 50-100, 150-200, ... m (see 00b_load_dvm.R)
  lo <- as.integer(sub("^(NASC|DVM)", "", code)) - 5
  sprintf("%d–%d m", lo, lo + 50)
}

# family label: with the estimated shape parameter when there is a final
# model, otherwise the same abbreviation without it (null-result rows)
family_label <- function(res, label) {
  if (part_of(label) == "when present" && identical(res$family, "gaussian")) return("Gaussian (logit)")
  if (is.null(res$final)) {
    return(switch(res$family, "negative binomial" = "NB", "Tweedie" = "Tweedie",
                  "gaussian" = "Gaussian", "binomial" = "Binomial", res$family))
  }
  fam <- res$final$model$family
  nm <- fam$family
  if (grepl("^Negative Binomial", nm)) sprintf("NB (θ = %.3f)", fam$getTheta(TRUE))
  else if (grepl("^Tweedie", nm)) sprintf("Tweedie (p = %.3f)", fam$getTheta(TRUE))
  else if (nm == "gaussian") "Gaussian"
  else if (nm == "binomial") "Binomial"
  else nm
}

fmt_p <- function(p) if (is.na(p)) "" else if (p < 0.001) "<0.001" else sprintf("%.3f", p)
fmt_num <- function(x, digits) ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))

#' detection days and n for one response: from the final model's own rows
#' when there is one, otherwise from the full record of that response.
#' Occurrence: days present; when present: every day in the model.
nonzero_info <- function(label, res) {
  y <- if (!is.null(res$final)) res$final$model$y else daily[[res$y]]
  y <- y[!is.na(y)]
  nz <- switch(part_of(label), occurrence = sum(y == 1), "when present" = length(y), sum(y > 0))
  list(n = length(y), nz = nz)
}

#' term-level rows for one response's final model, matched against its
#' shrinkage GAM for the same terms
term_rows <- function(res) {
  m <- res$final$model
  s <- summary(m)
  st <- s$s.table
  labs <- sub("^s\\((.*)\\)$", "\\1", rownames(st))
  shrink_p <- rep(NA_real_, length(labs))
  if (!is.null(res$shrinkage)) {
    ss <- summary(res$shrinkage)$s.table
    idx <- match(rownames(st), rownames(ss))
    shrink_p <- ifelse(is.na(idx), NA_real_, ss[idx, "p-value"])
  }
  data.frame(term = labs, edf = st[, "edf"],
             F = st[, if ("F" %in% colnames(st)) "F" else "Chi.sq"],
             p = st[, "p-value"], shrinkage_p = shrink_p,
             dev_expl = model_dev_expl(m) * 100, adj_r2 = s$r.sq, stringsAsFactors = FALSE)
}

#' one response's block, all columns pre-formatted as text: per-response
#' values only on the first term row
response_block <- function(label, display_name, res, with_nonzero) {
  if (is.null(res)) return(NULL)   # fit never attempted/completed -- a real gap, not a null result
  nzi <- nonzero_info(label, res)
  exploratory <- with_nonzero && nzi$nz < EXPLORATORY_MIN_NONZERO
  name <- if (exploratory) paste0(display_name, DAGGER) else display_name
  nz_txt <- sprintf("%d (%.1f)", nzi$nz, 100 * nzi$nz / nzi$n)

  if (is.null(res$final)) {
    # a real null result (no candidate passed screening, or every selected
    # term was pruned) -- reported explicitly, not dropped
    return(data.frame(group = name, family = family_label(res, label), term = "— (no term retained)",
                      edf = "", F = "", p = "", shrinkage_p = "", n = as.character(nzi$n),
                      nonzero = nz_txt, rho = "", dev_expl = "", cv_dev_expl = "", adj_r2 = "",
                      acf_before = "", acf_after = "", stringsAsFactors = FALSE))
  }
  rows <- term_rows(res)
  k <- nrow(rows)
  first <- function(x) c(x, rep("", k - 1))
  data.frame(
    group = first(name), family = first(family_label(res, label)), term = rows$term,
    edf = fmt_num(rows$edf, 2), F = fmt_num(rows$F, 2),
    p = vapply(rows$p, fmt_p, character(1)),
    shrinkage_p = vapply(rows$shrinkage_p, fmt_p, character(1)),
    n = first(as.character(res$final$n)), nonzero = first(nz_txt),
    rho = first(fmt_num(res$final$rho, 2)),
    dev_expl = first(fmt_num(rows$dev_expl[1], 1)),
    cv_dev_expl = first(fmt_num(cv_for(label), 1)),
    adj_r2 = first(fmt_num(rows$adj_r2[1], 2)),
    acf_before = first(fmt_num(res$final$resid_acf1_before, 3)),
    acf_after = first(fmt_num(res$final$resid_acf1_after, 3)),
    stringsAsFactors = FALSE
  )
}

build_table <- function(labels, display_names, with_nonzero = FALSE) {
  blocks <- Map(function(lab, disp) {
    res <- readRDS(file.path(MODEL_DIR, paste0(lab, "_model.rds")))
    response_block(lab, disp, res, with_nonzero)
  }, labels, display_names)
  dplyr::bind_rows(blocks[!vapply(blocks, is.null, logical(1))])
}

# ---- Table S-A: NASC backscatter (upper) + DVM amplitude (lower) ----------
nasc_codes <- c("NASC5", "NASC55", "NASC155", "NASC255", "NASC355", "NASC455")
dvm_codes  <- c("DVM5", "DVM55", "DVM155", "DVM255", "DVM355", "DVM455")

tblA_nasc <- build_table(nasc_codes, nasc_depth_label(nasc_codes))
tblA_dvm  <- build_table(dvm_codes, nasc_depth_label(dvm_codes))

# ---- Table S-B: cetacean/odontocete responses, grouped by taxon category --
# rows come from 00f_model_specs.R (two parts where defined), grouped by taxon
taxon_group <- list(baleen = c("BmB", "BmD", "Bp20Hz", "Bp40Hz", "Mn"),
                    delphinids = c("Gg", "Lo", "Oo", "Pm", "UO", "Kspp"),
                    beaked = c("Zc", "Mc", "Mg", "Ms"))
group_specs <- function(g) {
  s <- specs[specs$species %in% taxon_group[[g]], ]
  s[order(match(s$species, taxon_group[[g]])), ]
}
build_group <- function(g) {
  s <- group_specs(g)
  build_table(s$label, s$display, with_nonzero = TRUE)
}

tblB_baleen <- build_group("baleen")
tblB_delph  <- build_group("delphinids")
tblB_beaked <- build_group("beaked")

# ---- term-level CSVs (plain data, for anyone who wants to re-derive a
# table in a different tool rather than parse the docx) ----------------------
write.csv(tblA_nasc, file.path(MODEL_DIR, "table_SA_nasc_terms.csv"), row.names = FALSE)
write.csv(tblA_dvm, file.path(MODEL_DIR, "table_SA_dvm_terms.csv"), row.names = FALSE)
write.csv(rbind(cbind(category = "Baleen whales", tblB_baleen),
                cbind(category = "Delphinids and other odontocetes", tblB_delph),
                cbind(category = "Beaked whales", tblB_beaked)),
          file.path(MODEL_DIR, "table_SB_cetacean_terms.csv"), row.names = FALSE)

# ---- Word document (styling from 00e_word_tables.R) -------------------------
header_labels <- c(group = "Depth layer / Taxon", family = "Family", term = "Smooth term",
                   edf = "edf", F = "F", p = "p", shrinkage_p = "Shrinkage p",
                   n = "n", nonzero = "Non-zero days (%)", rho = "ρ (AR1)",
                   dev_expl = "Dev. expl. (%)", cv_dev_expl = "CV dev. expl. (%)",
                   adj_r2 = "Adj. R²",
                   acf_before = "Resid. ACF1 before", acf_after = "Resid. ACF1 after")
col_width <- c(group = 1.6, family = 1.0, term = 1.05, edf = 0.42, F = 0.45, p = 0.5,
               shrinkage_p = 0.62, n = 0.38, nonzero = 0.72, rho = 0.48, dev_expl = 0.55,
               cv_dev_expl = 0.55, adj_r2 = 0.45, acf_before = 0.55, acf_after = 0.55)

make_ft <- function(df, section_label, with_nonzero = FALSE) {
  keys <- names(header_labels)
  if (!with_nonzero) keys <- setdiff(keys, "nonzero")
  booktabs_ft(df, header_labels[keys], col_width[keys],
              left_cols = c("group", "family", "term"), section_label = section_label)
}

footnote_text <- paste(
  "Terms are thin-plate cubic regression splines with shrinkage (k = 4) selected by the",
  "automated screening, concurvity-pruning, forward-AIC and post-selection significance-pruning",
  "procedure (terms with edf < 0.5 or p ≥ 0.05 in the final joint fit, or p ≥ 0.05 in the",
  "all-candidate double-penalty shrinkage model, are dropped and the model refit, iterated to",
  "convergence; see Methods); one row per smooth term retained in the final model. edf, F and p",
  "are for the final AR1 model (bam, discrete fitting; approximate Wald-type tests for penalized",
  "smooths). Shrinkage p is the p-value for the same term in the shrinkage model (< 0.05 for every",
  "retained term by construction). ρ is the lag-1 residual autocorrelation estimated from the",
  "model without a correlation structure and used as the AR1 coefficient (0 = no AR1 term applied).",
  "Dev. expl. and adj. R² refer to the final fitted model. CV dev. expl. is the deviance explained",
  "when each calendar month is predicted by the model refit without it (leave-one-month-out",
  "cross-validation), relative to an intercept-only model; values near or below 0 mean the model",
  "does not predict unseen months better than the mean. Resid. ACF1 is the lag-1 autocorrelation",
  "of deviance residuals before and after applying AR1. Variable abbreviations: T, temperature;",
  "S, salinity; TAnom/SAnom, temperature/salinity anomaly; NO3, nitrate; Chl, chlorophyll a;",
  "DO, dissolved oxygen; FSLE, finite-size Lyapunov exponent; NASC, nautical area scattering",
  "coefficient; DVM, diel vertical migration amplitude (night minus day NASC); numbers give depth (m)."
)
two_part_note <- paste(
  "Two-part models: \"occurrence\" models whether the call type or species was detected on a",
  "day (binomial, all days); \"when present\" models the amount of calling on detection days",
  "only (NB for calls per day; Gaussian on the logit of the proportion of the day with calls for",
  "hours or minutes with calls). Occurrence models were fit with gamm() and an AR1 correlation",
  "(penalized quasi-likelihood), because the bam() AR1 used elsewhere is not reliable for binary",
  "data; as PQL has no true likelihood, their screening and forward selection used term p-values",
  "instead of AIC, their shrinkage model relies on the shrinkage basis alone, and ρ is the",
  "estimated AR1 correlation. With a single year of data, each seasonal species has essentially",
  "one presence season, so occurrence relationships cannot be separated from seasonality",
  "(see Methods)."
)
exploratory_note <- paste0(
  "Non-zero days: days with at least one detection among the n days used in the fit (or, for",
  " responses with no retained term, the full record), with the percentage of n in parentheses;",
  " for when-present models every day in the fit is a detection day. ",
  DAGGER, " Exploratory: fewer than ", EXPLORATORY_MIN_NONZERO, " detection days. With so few",
  " detection days a fit can depend on one or a few days, and deviance explained overstates",
  " predictive skill; these models are reported for completeness and should not be interpreted",
  " as established relationships."
)

doc <- landscape_doc()

doc <- add_caption(doc, "Table S-A.",
  paste("Summary of generalized additive models (GAMs) with AR1 error structure relating mean",
        "acoustic backscatter (NASC; upper block) and diel vertical migration (DVM) amplitude",
        "(night minus day NASC; lower block) in six 50-m depth layers to oceanographic variables.",
        "Backscatter itself was not a candidate predictor in the DVM models."))
doc <- body_add_flextable(doc, make_ft(tblA_nasc, "Mean backscatter (NASC)"), align = "left")
doc <- add_gap(doc)
doc <- body_add_flextable(doc, make_ft(tblA_dvm, "DVM amplitude (night − day NASC)"), align = "left")
doc <- add_note(doc, footnote_text)

doc <- body_add_break(doc)
doc <- add_caption(doc, "Table S-B.",
  paste("Summary of GAMs with AR1 error structure relating cetacean acoustic occurrence and",
        "calling to oceanographic variables, mean backscatter (NASC), and DVM amplitude, which",
        "was a candidate predictor for every response. Species with long presence seasons or",
        "frequent visits are modelled in two parts (occurrence and when present). NB, negative",
        "binomial (θ, dispersion); Tweedie (p, power parameter)."))
doc <- body_add_flextable(doc, make_ft(tblB_baleen, "Baleen whales", TRUE), align = "left")
doc <- add_gap(doc)
doc <- body_add_flextable(doc, make_ft(tblB_delph, "Delphinids and other odontocetes", TRUE), align = "left")
doc <- add_gap(doc)
doc <- body_add_flextable(doc, make_ft(tblB_beaked, "Beaked whales", TRUE), align = "left")
doc <- add_note(doc, footnote_text)
doc <- add_note(doc, two_part_note)
doc <- add_note(doc, exploratory_note)

save_docx_and_pdf(doc, file.path(MODEL_DIR, "GAM_AR1_model_summary_tables.docx"))
