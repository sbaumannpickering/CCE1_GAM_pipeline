# Automated GAM pipeline for the CCE1 analysis.
#   - automated univariate screening (loop)
#   - automated, iterative collinearity pruning via mgcv::concurvity()
#     (the GAM-appropriate analogue of VIF)
#   - automated forward-stepwise term addition by AIC
#   - autocorrelation handled inside the fitted model via bam(..., rho=...),
#     at native daily resolution
#   - automatic family selection (nb()/tw() estimate their shape parameter
#   - a double-penalty shrinkage GAM (select=TRUE) fit on the same pruned
#     candidate set, as an independent check on the stepwise pick: any
#     selected term with shrinkage p >= 0.05 is dropped from the final model
#
# Core design choice: bam(..., discrete=TRUE) with an AR1 correlation on
# the residuals, with rho estimated from a pilot fit and AR.start derived
# automatically from actual date gaps in whatever rows survive complete-case
# filtering for a given formula (robust to the handful of scattered NAs in
# Chl/NO3/temp/psal -- see 00a_load_data.R). The full CCE1 record has no
# deployment gaps (386 consecutive days, verified), so AR.start only needs
# to re-trigger when a specific model's predictors drop a day.

library(mgcv)
library(dplyr)

# ---- helpers --------------------------------------------------------------

#' AR.start flags from a (possibly gappy, after NA-dropping) date vector
make_ar_start <- function(dates) {
  dates <- as.Date(dates)
  gap <- c(TRUE, diff(dates) > 1)
  gap
}

#' lag-1 autocorrelation of residuals, respecting AR.start breaks
lag1_acf_safe <- function(resid, ar_start) {
  keep <- which(!ar_start)                 # rows that have a valid "previous" row
  if (length(keep) < 10) return(0)
  r_t  <- resid[keep]
  r_tm1 <- resid[keep - 1L]
  suppressWarnings(cor(r_t, r_tm1, use = "complete.obs"))
}

#' lag-1 ACF of the GLS-WHITENED residuals e_t = r_t - rho*r_{t-1}.
#' This, not the raw residual ACF, is the correct post-fit diagnostic: a
#' GLS/AR1 fit does not and should not remove autocorrelation from the raw
#' residuals (y - fitted) -- that correlation is structural to the assumed
#' error process. What should go to ~0 is the whitened residual sequence.
whitened_acf1 <- function(resid, ar_start, rho) {
  keep <- which(!ar_start)
  if (length(keep) < 10) return(NA_real_)
  e <- resid[keep] - rho * resid[keep - 1L]
  if (length(e) < 3 || stats::sd(e) == 0) return(NA_real_)
  suppressWarnings(cor(e[-1], e[-length(e)]))
}

#' Fit mgcv::bam() from a named list of arguments; returns NULL if the fit fails.
#'
#' The arguments are passed with do.call() rather than forwarded through `...`.
#' bam() evaluates AR.start with substitute(), which fails with
#' "..N used in an incorrect context, no ... to look in" when the call is built
#' from a wrapper's own `...`. do.call() supplies already-evaluated values and
#' avoids this.
safe_bam <- function(args) {
  tryCatch(suppressWarnings(do.call(bam, args)), error = function(e) NULL)
}

#' choose a response family automatically
#'  - integer counts -> negative binomial (nb(); theta estimated internally)
#'  - non-negative continuous (e.g. echolocation-minutes, NASC) -> Tweedie
#'    (tw(); power parameter estimated internally) unless gaussian_ok=TRUE
#'  - gaussian_ok=TRUE forces identity-link Gaussian (used for NASC/sA
#'    response models, matching the family already reported in the MS)
pick_family <- function(y, gaussian_ok = FALSE) {
  if (gaussian_ok) return(gaussian())
  is_integer_like <- all(abs(y - round(y)) < 1e-8, na.rm = TRUE)
  if (is_integer_like) nb() else tw()
}

#' fit one bam() with an AR1 correlation structure. rho is a single ACF-based
#' estimate from a pilot (rho=0) fit's residuals, following the standard
#' mgcv/itsadug workflow, not an AIC grid search. Each extra bam(rho=...) call
#' is another chance to hit an mgcv discrete-AR1 PIRLS non-convergence hang at
#' certain rho values, which R-level timeouts cannot interrupt (setTimeLimit
#' cannot stop a hang inside compiled code). Two bam() calls per model (pilot
#' plus one rho fit) keep that exposure low.
#' AR.start is derived from this formula's own complete-case dates, so a
#' handful of scattered predictor NAs correctly break the AR1 chain instead
#' of silently splicing non-adjacent days together.
fit_ar1_bam <- function(formula, data, family, date_col = "Date", k = 4) {
  if (is_pql(family)) return(fit_ar1_gamm(formula, data, family, date_col))
  vars <- all.vars(formula)
  cc <- stats::complete.cases(data[, unique(c(vars, date_col)), drop = FALSE])
  d <- data[cc, , drop = FALSE]
  d <- d[order(d[[date_col]]), ]
  ar_start <- make_ar_start(d[[date_col]])

  base_args <- list(formula = formula, data = d, family = family,
                     method = "fREML", discrete = TRUE, AR.start = ar_start)

  pilot <- safe_bam(base_args)
  if (is.null(pilot)) return(NULL)
  resid_acf1_before <- lag1_acf_safe(resid(pilot, type = "deviance"), ar_start)
  rho <- max(min(resid_acf1_before, 0.95), 0)
  if (rho < 0.02) rho <- 0                  # no meaningful autocorrelation

  fit <- if (rho == 0) pilot else safe_bam(c(base_args, list(rho = rho)))
  if (is.null(fit)) return(NULL)
  resid_acf1_after <- whitened_acf1(resid(fit, type = "deviance"), ar_start, rho)

  list(model = fit, rho = rho, n = nrow(d),
       resid_acf1_before = round(resid_acf1_before, 3),
       resid_acf1_after  = round(resid_acf1_after, 3))
}

#' Binary (presence/absence) responses are fit differently. bam()'s rho-based
#' AR1 acts on the working model and collapses for binary data with strong
#' persistence (tested on humpback occurrence: rho 0.80 -> edf 0, p reported
#' as 0, deviance explained 0%). These models are instead fit with gamm()
#' (PQL) and a proper AR1 correlation indexed by day number, so missing days
#' are spaced correctly. PQL has no true likelihood, so for these models
#' screening and forward selection use p-values instead of AIC (see
#' screen_univariate() and forward_select()).
is_pql <- function(family) identical(family$family, "binomial")

fit_ar1_gamm <- function(formula, data, family, date_col = "Date") {
  vars <- all.vars(formula)
  cc <- stats::complete.cases(data[, unique(c(vars, date_col)), drop = FALSE])
  d <- data[cc, , drop = FALSE]
  d <- d[order(d[[date_col]]), ]
  d$day_t <- as.numeric(d[[date_col]])
  ar_start <- make_ar_start(d[[date_col]])

  pilot <- tryCatch(suppressWarnings(gam(formula, data = d, family = family, method = "REML")),
                    error = function(e) NULL)
  if (is.null(pilot)) return(NULL)
  resid_acf1_before <- lag1_acf_safe(resid(pilot, type = "deviance"), ar_start)

  g <- tryCatch(suppressWarnings(gamm(formula, data = d, family = family,
                                      correlation = corAR1(form = ~ day_t), verbosePQL = FALSE)),
                error = function(e) NULL)
  if (is.null(g)) return(NULL)
  phi <- unname(coef(g$lme$modelStruct$corStruct, unconstrained = FALSE))
  resid_acf1_after <- lag1_acf_safe(resid(g$lme, type = "normalized"), ar_start)
  # deviance explained of the fitted mean, against the intercept-only binomial
  null_dev <- sum(family$dev.resids(d[[vars[1]]], mean(d[[vars[1]]]), rep(1, nrow(d))))
  fit_dev <- sum(family$dev.resids(d[[vars[1]]], fitted(g$gam), rep(1, nrow(d))))
  g$gam$pql_dev_expl <- 1 - fit_dev / null_dev

  list(model = g$gam, rho = phi, n = nrow(d),
       resid_acf1_before = round(resid_acf1_before, 3),
       resid_acf1_after = round(resid_acf1_after, 3))
}

#' build a y ~ s(x, k=k, bs="cs") formula referring to columns of `data`
.s_formula <- function(y_name, x_name, k = 4) {
  as.formula(sprintf("%s ~ s(%s, k=%d, bs='cs')", y_name, x_name, k))
}

# ---- 1. automated univariate screening ------------------------------------

#' fit y ~ s(x) for every candidate x (AR1-aware), return a tidy ranked table
screen_univariate <- function(y_name, candidates, data, family,
                               date_col = "Date", k = 4, verbose = TRUE) {
  rows <- vector("list", length(candidates))
  for (i in seq_along(candidates)) {
    x <- candidates[i]
    t0 <- Sys.time()
    f <- .s_formula(y_name, x, k)
    fit <- tryCatch(fit_ar1_bam(f, data, family, date_col, k),
                     error = function(e) NULL)
    if (verbose) message(sprintf("  screening %-12s (%d/%d) [%.1fs]", x, i,
                                  length(candidates),
                                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    if (is.null(fit)) { rows[[i]] <- NULL; next }
    s <- summary(fit$model)
    rows[[i]] <- data.frame(
      predictor = x, aic = if (is_pql(family)) NA_real_ else AIC(fit$model), edf = s$edf[1],
      p_value = s$s.table[1, "p-value"], dev_expl = model_dev_expl(fit$model),
      rho = fit$rho, n = fit$n
    )
    rm(fit); gc(verbose = FALSE)   # bam() fit objects are large; force collection
                                    # each iteration -- repeated fits in one long
                                    # R session can stall (likely GC pressure)
  }
  out <- dplyr::bind_rows(rows)
  # ranked by AIC, or by p-value for PQL (binary) models, which have no AIC
  out[order(if (is_pql(family)) out$p_value else out$aic), ]
}

#' deviance explained; PQL (gamm) fits carry their own, computed in fit_ar1_gamm()
model_dev_expl <- function(model) {
  if (!is.null(model$pql_dev_expl)) model$pql_dev_expl else summary(model)$dev.expl
}

# ---- 2. automated, iterative concurvity pruning ----------------------------

#' among `candidates` (already filtered to significant ones), iteratively
#' drop the weaker member (by univariate AIC rank in `screen_tbl`) of any
#' pair whose worst-case concurvity exceeds `threshold`, until none remain.
#' Uses mgcv::concurvity() on one combined model -- the GAM analogue of VIF.
prune_concurvity <- function(y_name, candidates, screen_tbl, data, family,
                              date_col = "Date", k = 4, threshold = 0.6,
                              verbose = TRUE) {
  keep <- candidates
  repeat {
    if (length(keep) < 2) break
    f <- as.formula(paste(
      y_name, "~", paste(sprintf("s(%s, k=%d, bs='cs')", keep, k), collapse = " + ")
    ))
    fit <- tryCatch(fit_ar1_bam(f, data, family, date_col, k),
                     error = function(e) NULL)
    if (is.null(fit)) { keep <- keep[-length(keep)]; next }  # drop last on failure, retry
    cm <- concurvity(fit$model, full = FALSE)$worst
    # strip "s(" / ")" around each smooth's single variable name (a bare s(X)
    # has no comma, so a comma-anchored strip would leave a trailing ")" and no
    # label would match `keep`, so drop_var would be empty and the pruning
    # loop would never shrink `keep`)
    labs <- sub("^s\\(", "", rownames(cm))
    labs <- sub("\\)$", "", labs)
    labs <- sub(",.*$", "", labs)           # defensive: handles multi-arg smooths too
    keep_idx <- labs != "para"               # drop the parametric/intercept row-col
    cm <- cm[keep_idx, keep_idx, drop = FALSE]
    labs <- labs[keep_idx]
    rownames(cm) <- colnames(cm) <- labs
    diag(cm) <- 0
    offending <- which(cm > threshold, arr.ind = TRUE)
    offending <- offending[offending[, 1] < offending[, 2], , drop = FALSE]
    if (nrow(offending) == 0) break
    pair <- labs[offending[1, ]]
    ranks <- match(pair, screen_tbl$predictor)
    drop_var <- pair[which.max(ranks)]      # worse (higher) AIC rank = weaker evidence
    if (verbose) message(sprintf("  dropping %s (concurvity %.2f with %s)",
                                  drop_var, cm[offending[1, 1], offending[1, 2]],
                                  setdiff(pair, drop_var)))
    rm(fit); gc(verbose = FALSE)
    keep <- setdiff(keep, drop_var)
  }
  keep
}

# ---- 3. automated forward-stepwise selection -------------------------------

#' greedy forward addition by AIC; stops when no remaining candidate improves
#' the current model's AIC by at least min_improve (a bare numerical
#' improvement accepts noise-level AIC fluctuations and overfits -- requiring
#' a minimum drop is the same logic as the conventional deltaAIC > 2 rule of
#' thumb for "meaningfully better"), or when max_terms is reached (mirrors
#' the ~5-term ceiling the manuscript's manual stepwise process used in
#' practice)
forward_select <- function(y_name, candidates, data, family,
                            date_col = "Date", k = 4, verbose = TRUE,
                            max_terms = 6, min_improve = 2) {
  # PQL (binary) models have no AIC: add the candidate whose own term has the
  # smallest p-value, and stop when none is below alpha
  use_p <- is_pql(family)
  alpha <- 0.05
  selected <- character(0)
  remaining <- candidates
  cur_aic <- Inf
  log <- list()
  repeat {
    if (length(remaining) == 0 || length(selected) >= max_terms) break
    trial_score <- numeric(length(remaining))
    for (j in seq_along(remaining)) {
      x <- remaining[j]
      t0 <- Sys.time()
      terms <- c(selected, x)
      f <- as.formula(paste(
        y_name, "~", paste(sprintf("s(%s, k=%d, bs='cs')", terms, k), collapse = " + ")
      ))
      fit <- tryCatch(fit_ar1_bam(f, data, family, date_col, k),
                       error = function(e) NULL)
      trial_score[j] <- if (is.null(fit)) Inf
                        else if (use_p) summary(fit$model)$s.table[paste0("s(", x, ")"), "p-value"]
                        else AIC(fit$model)
      if (verbose) message(sprintf("    trying + %-12s [%.1fs]", x,
                                    as.numeric(difftime(Sys.time(), t0, units = "secs"))))
      rm(fit); gc(verbose = FALSE)
    }
    best_i <- which.min(trial_score)
    if (use_p) {
      if (trial_score[best_i] >= alpha) break               # nothing significant to add -> stop
    } else if (trial_score[best_i] > cur_aic - min_improve) break   # not meaningfully better -> stop
    selected <- c(selected, remaining[best_i])
    if (!use_p) cur_aic <- trial_score[best_i]
    log[[length(log) + 1]] <- data.frame(step = length(selected), added = remaining[best_i],
                                          aic = if (use_p) NA_real_ else cur_aic,
                                          p_added = if (use_p) trial_score[best_i] else NA_real_)
    if (verbose) message(sprintf("  step %d: + %s (%s=%.3g)", length(selected), remaining[best_i],
                                  if (use_p) "p" else "AIC", trial_score[best_i]))
    remaining <- setdiff(remaining, remaining[best_i])
  }
  list(selected = selected, log = dplyr::bind_rows(log))
}

# ---- 3b. post-selection significance pruning -------------------------------

#' drop terms from the forward-selected model that turn out to contribute
#' ~nothing once fit jointly with the other selected terms. Forward selection
#' only checks whole-model AIC at each step, so a term that clears the AIC
#' bar in isolation can still end up non-significant (or shrunk to ~0 edf by
#' the bs='cs' shrinkage basis) once it shares variance with a later-added
#' term at a concurvity below the 0.6 pruning threshold -- that threshold
#' screens out strong collinearity but deliberately allows some. One term is
#' dropped per iteration (the weakest by edf among those failing either
#' criterion), and the remaining terms refit, since removing one term can
#' change whether another is still weak. Returns both the pruned term list
#' and its refit; an empty result (all terms pruned) is reported the same
#' way as "no candidate significant" upstream.
prune_final_terms <- function(y_name, selected, data, family, date_col = "Date",
                               k = 4, edf_floor = 0.5, alpha = 0.05, verbose = TRUE) {
  terms <- selected
  repeat {
    if (length(terms) == 0) return(list(selected = terms, final = NULL))
    f <- as.formula(paste(
      y_name, "~", paste(sprintf("s(%s, k=%d, bs='cs')", terms, k), collapse = " + ")
    ))
    fit <- fit_ar1_bam(f, data, family, date_col, k)
    if (is.null(fit)) return(list(selected = terms, final = NULL))
    s <- summary(fit$model)
    edf <- s$s.table[, "edf"]
    pval <- s$s.table[, "p-value"]
    weak <- edf < edf_floor | pval >= alpha
    if (!any(weak)) return(list(selected = terms, final = fit))
    drop_i <- which(weak)[which.min(edf[which(weak)])]
    if (verbose) message(sprintf("  pruning %s (edf=%.2f, p=%.3f) from final model",
                                  terms[drop_i], edf[drop_i], pval[drop_i]))
    rm(fit); gc(verbose = FALSE)
    terms <- terms[-drop_i]
  }
}

# ---- 4. top-level wrapper ---------------------------------------------------

#' full automated pipeline for one response variable
auto_gam_AR1 <- function(y_name, candidates, data, family = NULL,
                          gaussian_ok = FALSE, date_col = "Date", k = 4,
                          screen_alpha = 0.05, concurvity_threshold = 0.6,
                          verbose = TRUE) {
  if (is.null(family)) family <- pick_family(data[[y_name]], gaussian_ok)

  if (verbose) message(sprintf("\n=== %s: univariate screening (%d candidates) ===",
                                y_name, length(candidates)))
  screen_tbl <- screen_univariate(y_name, candidates, data, family, date_col, k, verbose)
  sig <- screen_tbl$predictor[screen_tbl$p_value < screen_alpha]

  if (length(sig) == 0) {
    warning(sprintf("%s: no candidate significant at p<%.2f; returning intercept-only",
                     y_name, screen_alpha))
    return(list(y = y_name, family = family$family, screen_tbl = screen_tbl,
                pruned = character(0), selected_forward = character(0),
                selected = character(0), final = NULL, shrinkage = NULL))
  }

  # cap the candidate set relative to the response's effective (non-zero)
  # sample size before concurvity pruning. Fitting many smooth terms at once
  # against very few informative observations (e.g. Kspp: 14 non-zero days
  # of 386) is statistically unjustified and can trigger an mgcv PIRLS
  # non-convergence hang in the combined-term pruning fit, deterministically,
  # independent of system resources. screen_tbl is
  # already AIC-sorted, so this keeps the strongest candidates.
  # informative = non-zero days for counts/Tweedie; for a Gaussian response
  # (which can be negative) every observed day; for occurrence (binary), the
  # rarer of presence and absence days
  yv <- data[[y_name]]
  n_informative <- if (identical(family$family, "gaussian")) sum(!is.na(yv))
                   else if (is_pql(family)) min(sum(yv == 1, na.rm = TRUE), sum(yv == 0, na.rm = TRUE))
                   else sum(yv > 0, na.rm = TRUE)
  max_candidates <- max(2, min(15, floor(n_informative / 5)))
  sig <- screen_tbl$predictor[screen_tbl$predictor %in% sig]   # re-order by AIC
  if (length(sig) > max_candidates) {
    if (verbose) message(sprintf(
      "  capping to top %d of %d significant candidates (n_informative=%d)",
      max_candidates, length(sig), n_informative))
    sig <- sig[seq_len(max_candidates)]
  }

  if (verbose) message(sprintf("=== %s: concurvity pruning (%d significant) ===",
                                y_name, length(sig)))
  pruned <- prune_concurvity(y_name, sig, screen_tbl, data, family, date_col, k,
                              concurvity_threshold, verbose)

  if (verbose) message(sprintf("=== %s: forward stepwise selection (%d candidates) ===",
                                y_name, length(pruned)))
  fs <- forward_select(y_name, pruned, data, family, date_col, k, verbose)

  # double-penalty shrinkage GAM on the full concurvity-pruned candidate set,
  # with the AR1 rho of the current final model
  fit_shrinkage <- function(rho) {
    if (length(pruned) == 0) return(NULL)
    f <- as.formula(paste(
      y_name, "~", paste(sprintf("s(%s, k=%d, bs='cs')", pruned, k), collapse = " + ")
    ))
    # PQL (binary) models: gamm() has no select=TRUE, so the all-candidate
    # model relies on the cs basis' own shrinkage, with its own AR1 estimate
    if (is_pql(family)) {
      fit <- fit_ar1_gamm(f, data, family, date_col)
      return(if (is.null(fit)) NULL else fit$model)
    }
    vars <- all.vars(f)
    cc <- stats::complete.cases(data[, unique(c(vars, date_col)), drop = FALSE])
    d <- data[cc, , drop = FALSE]; d <- d[order(d[[date_col]]), ]
    ar_start <- make_ar_start(d[[date_col]])
    safe_bam(list(formula = f, data = d, family = family, method = "fREML",
                  discrete = TRUE, select = TRUE, rho = rho, AR.start = ar_start))
  }

  # Post-selection pruning, two criteria: (1) edf < 0.5 or p >= 0.05 in the
  # final joint fit (prune_final_terms), and (2) p >= 0.05 for the same term
  # in the shrinkage GAM -- a stepwise-selected term that the all-candidate
  # double-penalty fit does not support is dropped too. Criterion (2) removes
  # every failing term at once (the shrinkage fit does not depend on which
  # terms were selected, only on rho), then (1) is re-run on what remains;
  # repeated until neither criterion drops anything, since the refit can
  # change rho and hence the shrinkage p-values.
  if (verbose) message(sprintf("=== %s: post-selection significance pruning (%d selected) ===",
                                y_name, length(fs$selected)))
  sel <- fs$selected
  dropped_shrinkage <- character(0)
  repeat {
    pf <- prune_final_terms(y_name, sel, data, family, date_col, k, verbose = verbose)
    final <- pf$final
    shrinkage <- fit_shrinkage(if (!is.null(final)) final$rho else 0)
    if (is.null(shrinkage) || length(pf$selected) == 0) break
    ss <- summary(shrinkage)$s.table
    sp <- ss[match(paste0("s(", pf$selected, ")"), rownames(ss)), "p-value"]
    bad <- pf$selected[!is.na(sp) & sp >= screen_alpha]
    if (length(bad) == 0) break
    if (verbose) message(sprintf("  pruning %s (shrinkage p=%s) from final model",
                                  paste(bad, collapse = ", "),
                                  paste(sprintf("%.3f", sp[pf$selected %in% bad]), collapse = ", ")))
    dropped_shrinkage <- c(dropped_shrinkage, bad)
    sel <- setdiff(pf$selected, bad)
  }

  list(y = y_name, family = family$family, screen_tbl = screen_tbl,
       sig = sig, pruned = pruned, selected_forward = fs$selected,
       selected = pf$selected, dropped_shrinkage = dropped_shrinkage,
       step_log = fs$log, final = final, shrinkage = shrinkage)
}
