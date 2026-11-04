# Response definitions for the cetacean models, shared by 01_run_all_models.R,
# 03_top_day_sensitivity.R, 04_build_summary_tables.R and 05_summary_collages.R.
#
# Species with long presence seasons or frequent visits are modelled in two
# parts, each through the full selection and pruning pipeline:
#   <sp>_occ  occurrence: detected (1) or not (0) on each day, full record;
#             binomial, fit with gamm() + AR1 (see fit_ar1_gamm() in 00d)
#   <sp>_pos  when present: amount of calling on detection days only
#             - call counts (blue whale B, D): negative binomial
#             - hours (fin whale 40 Hz, of 24) or minutes (of 1440) with calls:
#               Gaussian on the logit of the proportion of the day, squeezed
#               away from 0 and 1 (Smithson & Verkuilen 2006) so full days
#               are allowed; Gaussian is the case bam()'s AR1 handles best
# Risso's dolphin is occurrence only (81 detection days in 49 short visits is
# too thin for a when-present model). Fin whale 20 Hz (a continuous call
# index, present every day) is a single Gaussian model. Common dolphin
# (absent on 9 days) and the rare odontocetes (3-22 isolated detection
# days) keep a single model with the automatically chosen family.

two_part <- c(BmB = "count", BmD = "count", Bp40Hz = "hours", Mn = "minutes",
              Pm = "minutes", Zc = "minutes")
occurrence_only <- c("Gg")
single_model <- c(Bp20Hz = "gaussian", UO = "auto", Kspp = "auto", Lo = "auto",
                  Mc = "auto", Mp = "auto", Ms = "auto", Oo = "auto")

species_name <- c(BmB = "Blue whale B calls", BmD = "Blue whale D calls",
                  Bp20Hz = "Fin whale 20 Hz calls", Bp40Hz = "Fin whale 40 Hz calls",
                  Mn = "Humpback whale", Gg = "Risso's dolphin", Kspp = "Kogia spp.",
                  Lo = "Pacific white-sided dolphin", Mc = "Hubbs' beaked whale",
                  Mp = "Ginkgo-toothed beaked whale", Ms = "Stejneger's beaked whale",
                  Oo = "Killer whale", Pm = "Sperm whale", UO = "Common dolphin",
                  Zc = "Goose-beaked whale")
pos_unit <- c(count = "calls/day", hours = "hours with calls", minutes = "minutes with calls")

#' add <sp>_occ and <sp>_pos response columns to the daily data
add_response_columns <- function(daily) {
  for (sp in c(names(two_part), occurrence_only)) {
    daily[[paste0(sp, "_occ")]] <- as.integer(daily[[sp]] > 0)
  }
  for (sp in names(two_part)) {
    y <- daily[[sp]]
    present <- !is.na(y) & y > 0
    pos <- rep(NA_real_, length(y))
    if (two_part[[sp]] == "count") {
      pos[present] <- y[present]
    } else {
      prop <- y[present] / if (two_part[[sp]] == "hours") 24 else 1440
      n <- length(prop)
      prop <- (prop * (n - 1) + 0.5) / n   # squeeze into (0, 1)
      pos[present] <- qlogis(prop)
    }
    daily[[paste0(sp, "_pos")]] <- pos
  }
  daily
}

#' one row per cetacean model: label (= response column), source species,
#' part, family, and display name
cetacean_specs <- function() {
  rows <- list()
  add <- function(label, sp, part, family, display) {
    rows[[length(rows) + 1]] <<- data.frame(label = label, species = sp, part = part,
                                            family = family, display = display,
                                            stringsAsFactors = FALSE)
  }
  for (sp in names(species_name)) {
    if (sp %in% names(two_part)) {
      add(paste0(sp, "_occ"), sp, "occurrence", "binomial", paste0(species_name[[sp]], ": occurrence"))
      add(paste0(sp, "_pos"), sp, "when present",
          if (two_part[[sp]] == "count") "nb" else "gaussian",
          paste0(species_name[[sp]], ": ", pos_unit[[two_part[[sp]]]], " when present"))
    } else if (sp %in% occurrence_only) {
      add(paste0(sp, "_occ"), sp, "occurrence", "binomial", paste0(species_name[[sp]], ": occurrence"))
    } else {
      add(sp, sp, "single", single_model[[sp]], species_name[[sp]])
    }
  }
  do.call(rbind, rows)
}

#' family object for a spec's family key (NULL = choose automatically)
spec_family <- function(key) {
  switch(key, binomial = binomial(), nb = mgcv::nb(), gaussian = gaussian(), auto = NULL)
}

#' how a spec's family is shown in the tables
spec_family_note <- function(spec_row) {
  if (spec_row$part == "when present" && spec_row$family == "gaussian") "Gaussian, logit(prop. of day)"
  else NA_character_
}
