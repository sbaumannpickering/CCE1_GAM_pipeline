# Full batch run of the automated AR1 GAM pipeline across every response
# variable, each fit once, with the post-selection significance pruning from
# 00d_gam_pipeline_functions.R. DVM amplitude is in the candidate set for every
# cetacean/predator response, as an ordinary candidate alongside the
# oceanographic variables and NASC.
#
# DVM is NOT added for the six NASC backscatter-layer responses: DVM at a given
# depth is night-minus-day NASC at that same depth, so using it to predict
# NASC at that depth would be close to circular.
#
# Writes <response>_model.rds and <response>_summary.rds for each response,
# all_model_summary.csv, and dvm_term_summary.csv (which responses retained a
# DVM term in their final model, read off the selected_terms column of
# all_model_summary.csv).

source("config.R")
source("00a_load_data.R")
source("00b_load_dvm.R")
source("00d_gam_pipeline_functions.R")
source("00f_model_specs.R")

OUTPUT_DIR <- MODEL_DIR
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

daily <- load_cce1_daily()
# extract attributes BEFORE merging -- merge() drops them (see 02_dvm_models.R)
depth_ctd <- attr(daily, "depth_ctd")
nasc_depths <- attr(daily, "nasc_sv_depths")

dvm <- load_dvm_amplitude(validate = TRUE)
daily <- merge_dvm_into_daily(daily, dvm)

candidates_physical <- c(
  "Chl40", "DO500", "FSLE", "NO3",   # Chl20 excluded (untrusted data)
  paste0("S", depth_ctd), paste0("SAnom", depth_ctd),
  paste0("T", depth_ctd), paste0("TAnom", depth_ctd)
)
candidates_cetacean <- c(candidates_physical, paste0("NASC", nasc_depths),
                          paste0("DVM", DVM_DEPTHS))

nasc_responses <- paste0("NASC", nasc_depths)
daily <- add_response_columns(daily)
specs <- cetacean_specs()

summarize_result <- function(nm, r, elapsed_min) {
  if (is.null(r) || is.null(r$final)) {
    return(data.frame(response = nm, family = if (is.null(r)) NA else r$family,
                       n_candidates_screened = if (is.null(r)) NA else nrow(r$screen_tbl),
                       n_significant = if (is.null(r)) NA else length(r$sig),
                       n_selected_forward = if (is.null(r)) NA else length(r$selected_forward),
                       n_selected = 0, selected_terms = NA,
                       dropped_nonsig_terms = if (is.null(r)) NA
                         else paste(setdiff(r$selected_forward, r$selected), collapse = "; "),
                       dropped_shrinkage_terms = if (is.null(r)) NA else paste(r$dropped_shrinkage, collapse = "; "),
                       rho = NA, dev_expl = NA, resid_acf1_before = NA, resid_acf1_after = NA,
                       elapsed_min = round(elapsed_min, 1)))
  }
  s <- summary(r$final$model)
  data.frame(response = nm, family = r$family,
             n_candidates_screened = nrow(r$screen_tbl), n_significant = length(r$sig),
             n_selected_forward = length(r$selected_forward),
             n_selected = length(r$selected),
             selected_terms = paste(r$selected, collapse = "; "),
             dropped_nonsig_terms = paste(setdiff(r$selected_forward, r$selected), collapse = "; "),
             dropped_shrinkage_terms = paste(r$dropped_shrinkage, collapse = "; "),
             rho = round(r$final$rho, 3), dev_expl = round(model_dev_expl(r$final$model) * 100, 1),
             resid_acf1_before = r$final$resid_acf1_before,
             resid_acf1_after = r$final$resid_acf1_after,
             elapsed_min = round(elapsed_min, 1))
}

run_and_save <- function(label, y, candidates, data, gaussian_ok = FALSE, family = NULL) {
  model_path <- file.path(OUTPUT_DIR, paste0(label, "_model.rds"))
  summary_path <- file.path(OUTPUT_DIR, paste0(label, "_summary.rds"))

  message(sprintf("\n########## %s ##########", label))
  t0 <- Sys.time()
  res <- tryCatch(
    auto_gam_AR1(y, candidates, data, family = family, gaussian_ok = gaussian_ok, verbose = TRUE),
    error = function(e) { message("TOP-LEVEL ERROR: ", conditionMessage(e)); NULL }
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  message(sprintf("########## %s done in %.1f min ##########", label, elapsed))

  saveRDS(res, model_path)
  summary_row <- summarize_result(label, res, elapsed)
  saveRDS(summary_row, summary_path)

  rm(res); gc(verbose = FALSE)
  summary_row
}

summary_rows <- list()

for (y in nasc_responses) {
  summary_rows[[y]] <- run_and_save(y, y, candidates_physical, daily, gaussian_ok = TRUE)
}

# Cetacean models as defined in 00f_model_specs.R: two-part (occurrence +
# when present) for the baleen whales, sperm whale and goose-beaked whale,
# occurrence only for Risso's dolphin, single models otherwise. Families are
# set per model there; "auto" falls back to pick_family(). Fin whale 20 Hz is
# Gaussian: it is a continuous call index (20 Hz band energy relative to the
# surrounding bands), never zero, so Tweedie's zero mass does not apply.
# The two-part models replace the former calendar-based May-Nov blue whale D
# model: the when-present part uses the actual detection days instead.
for (i in seq_len(nrow(specs))) {
  sp <- specs[i, ]
  summary_rows[[sp$label]] <- run_and_save(sp$label, sp$label, candidates_cetacean, daily,
                                           family = spec_family(sp$family))
}

# ---- consolidated summary table --------------------------------------------
summary_tbl <- dplyr::bind_rows(summary_rows)
write.csv(summary_tbl, file.path(OUTPUT_DIR, "all_model_summary.csv"), row.names = FALSE)
print(summary_tbl)

# ---- which responses retained a DVM term, read off the table above --------
dvm_term_summary <- summary_tbl[grepl("DVM", summary_tbl$selected_terms), c("response", "selected_terms", "dev_expl")]
write.csv(dvm_term_summary, file.path(OUTPUT_DIR, "dvm_term_summary.csv"), row.names = FALSE)
cat("\nResponses retaining a DVM term in their final model:\n")
print(dvm_term_summary, row.names = FALSE)
