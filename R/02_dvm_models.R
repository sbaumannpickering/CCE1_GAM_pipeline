# DVM amplitude ~ water-mass GAMs: quantifies which oceanographic conditions
# predict strong vs. weak diel vertical migration at each depth layer,
# through the same automated AR1 pipeline used for every other response.
# Candidates are the physical water-mass variables only (not other NASC
# layers), since DVM amplitude is itself derived from backscatter.

source("config.R")
source("00a_load_data.R")
source("00b_load_dvm.R")
source("00d_gam_pipeline_functions.R")

OUTPUT_DIR <- MODEL_DIR
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

daily <- load_cce1_daily()
# extract depth_ctd BEFORE merging: merge() drops custom data.frame
# attributes, so attr(daily, "depth_ctd") would be NULL after the merge below,
# and paste0("S", NULL) returns "S" rather than an empty result (the model
# would screen a bare "S"/"T"/etc candidate instead of the 8 depth-expanded ones)
depth_ctd <- attr(daily, "depth_ctd")
candidates_physical <- c(
  "Chl40", "DO500", "FSLE", "NO3",
  paste0("S", depth_ctd), paste0("SAnom", depth_ctd),
  paste0("T", depth_ctd), paste0("TAnom", depth_ctd)
)

dvm <- load_dvm_amplitude(validate = TRUE)
daily <- merge_dvm_into_daily(daily, dvm)

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
             rho = round(r$final$rho, 3), dev_expl = round(s$dev.expl * 100, 1),
             resid_acf1_before = r$final$resid_acf1_before,
             resid_acf1_after = r$final$resid_acf1_after,
             elapsed_min = round(elapsed_min, 1))
}

run_and_save <- function(label, y, candidates, data) {
  model_path <- file.path(OUTPUT_DIR, paste0(label, "_model.rds"))
  summary_path <- file.path(OUTPUT_DIR, paste0(label, "_summary.rds"))
  # always refit (no skip-if-cached) so results reflect the current code and data
  message(sprintf("\n########## %s ##########", label))
  t0 <- Sys.time()
  res <- tryCatch(
    auto_gam_AR1(y, candidates, data, gaussian_ok = TRUE, verbose = TRUE),
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

dvm_responses <- paste0("DVM", DVM_DEPTHS)
summary_rows <- list()
for (y in dvm_responses) {
  summary_rows[[y]] <- run_and_save(y, y, candidates_physical, daily)
}

summary_tbl <- dplyr::bind_rows(summary_rows)
write.csv(summary_tbl, file.path(OUTPUT_DIR, "dvm_model_summary.csv"), row.names = FALSE)
print(summary_tbl)
