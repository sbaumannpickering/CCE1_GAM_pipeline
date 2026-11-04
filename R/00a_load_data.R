# Unified data loader for the CCE1 GAM pipeline.
#
# Reads combined_diel_daily.mat directly and returns one tidy daily data frame
# covering the full, continuous 386-day record at native (1-day) resolution.
# Autocorrelation is handled later, inside the GAM fit itself (see
# 00d_gam_pipeline_functions.R), not by pre-binning here.
#
# The North-South water-mass index is not a candidate predictor: the .mat file
# contains no such variable (the loader returns ns_available = FALSE).

source("config.R")
library(R.matlab)
library(dplyr)

load_cce1_daily <- function(
  mat_path = COMBINED_DAILY_MAT
) {
  m <- readMat(mat_path)

  depth_ctd <- as.numeric(m$depth.ctd)         # 20 40 75 150 300 500 750 1000
  wbat_depth <- as.numeric(m$WBAT.depth)       # 5 55 155 255 355 455

  matlab2date <- function(x) as.Date(x - 719529, origin = "1970-01-01")

  full_day <- as.numeric(m$full.day)
  stopifnot(all(diff(full_day) == 1))          # continuous daily series, no gaps

  n <- length(full_day)
  dates <- matlab2date(full_day)

  # period index (1-11) from boundaries (MATLAB datenum start/end per period)
  bnd <- m$boundaries
  period <- rep(NA_integer_, n)
  for (p in seq_len(nrow(bnd))) {
    period[full_day >= bnd[p, 1] & full_day <= bnd[p, 2]] <- p
  }

  mat_cols <- function(mat, depths, prefix) {
    mat <- as.matrix(mat)
    colnames(mat) <- paste0(prefix, depths)
    as.data.frame(mat)
  }

  temp     <- mat_cols(m$temp.perDay,    depth_ctd, "T")
  temp_an  <- mat_cols(m$t.anom.perDay,  depth_ctd, "TAnom")
  psal     <- mat_cols(m$psal.perDay,    depth_ctd, "S")
  psal_an  <- mat_cols(m$s.anom.perDay,  depth_ctd, "SAnom")
  spice    <- mat_cols(m$spice.perDay,   depth_ctd, "Spice")

  sv   <- mat_cols(m$Sv.perDay,   wbat_depth, "Sv")
  nasc <- mat_cols(m$NASC.perDay, wbat_depth, "NASC")

  responses <- data.frame(
    BmB    = as.numeric(m$BmB.perDay),
    BmD    = as.numeric(m$BmD.perDay),
    Bp20Hz = as.numeric(m$Bp20Hz.perDay),
    Bp40Hz = as.numeric(m$Bp40Hz.perDay),
    Mn     = as.numeric(m$Mn.perDay),
    Gg     = as.numeric(m$Gg.perDay),
    Kspp   = as.numeric(m$Kspp.perDay),
    Lo     = as.numeric(m$Lo.perDay),
    Mc     = as.numeric(m$Mc.perDay),
    Mp     = as.numeric(m$Mp.perDay),   # likely ginkgo-toothed beaked whale (M. ginkgodens),
                                         # provisionally coded pre-Henderson et al. 2025 -- confirm with SBP
    Ms     = as.numeric(m$Ms.perDay),
    Oo     = as.numeric(m$Oo.perDay),
    Pm     = as.numeric(m$Pm.perDay),
    UO     = as.numeric(m$UO.perDay),   # common dolphin proxy (per ModelSelection_UO.R usage)
    Zc     = as.numeric(m$Zc.perDay)
  )

  daily <- data.frame(
    Date = dates,
    JD = as.integer(format(dates, "%j")),
    Period = period
  )
  daily <- cbind(daily, temp, temp_an, psal, psal_an, spice,
                  DO500 = as.numeric(m$O2.509.perDay),
                  NO3   = as.numeric(m$NO3.40.perDay),
                  # chlorophyll at 20 m (21 m sensor) is excluded: its data
                  # are not trusted, and the field was removed from
                  # combined_diel_daily.mat on 2026-10-11
                  Chl40 = as.numeric(m$Chl.41.perDay),
                  FSLE    = as.numeric(m$FSLE),
                  FSLEstd = as.numeric(m$FSLEstd),
                  sv, nasc, responses)

  attr(daily, "ns_available") <- FALSE
  attr(daily, "response_vars") <- names(responses)
  attr(daily, "depth_predictor_stems") <- c("T", "TAnom", "S", "SAnom", "Spice")
  attr(daily, "depth_ctd") <- depth_ctd
  attr(daily, "scalar_predictors") <- c("DO500", "NO3", "Chl40", "FSLE", "FSLEstd")
  attr(daily, "nasc_sv_depths") <- wbat_depth

  daily
}

# restricted blue-whale-D frame (May-Nov presence window; Methods L285-286),
# built from the same unified loader rather than a separately exported file
load_cce1_bmd_restricted <- function(daily = load_cce1_daily()) {
  dplyr::filter(daily, as.integer(format(Date, "%m")) %in% 5:11)
}

if (sys.nframe() == 0) {
  daily <- load_cce1_daily()
  cat("Loaded", nrow(daily), "rows x", ncol(daily), "cols\n")
  cat("Date range:", format(range(daily$Date)), "\n")
  # validated against Lo_daily.txt row 1 (387 rows = unbinned full-resolution
  # reference; NASC255_daily.txt/Zc_daily.txt are ITS-block-averaged and are
  # NOT valid single-day references)
  cat("Spot check (row 1) against unbinned Lo_daily.txt reference:\n")
  cat("  T20     =", daily$T20[1],     " (expect 16.87368751)\n")
  cat("  T1000   =", daily$T1000[1],   " (expect 3.937226226)\n")
  cat("  S500    =", daily$S500[1],    " (expect 34.25913334)\n")
  cat("  NASC255 =", daily$NASC255[1], " (expect 157.7768188)\n")
}
