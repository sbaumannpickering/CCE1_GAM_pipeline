# Unified data loader for the CCE1 GAM pipeline.
#
# Reads CCE1_WBAT_HARP_daily.mat directly and returns one tidy daily data frame
# covering the full, continuous 386-day record at native (1-day) resolution.
# Autocorrelation is handled later, inside the GAM fit itself (see
# 00d_gam_pipeline_functions.R), not by pre-binning here.

source("config.R")
library(R.matlab)
library(dplyr)

load_cce1_daily <- function(
  mat_path = DAILY_MAT
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
    BmB    = as.numeric(m$BmB.perDay),  # blue whale B call
    BmD    = as.numeric(m$BmD.perDay),  # blue whale D call
    Bp20Hz = as.numeric(m$Bp20Hz.perDay),  # fin whale 20-Hz call index
    Bp40Hz = as.numeric(m$Bp40Hz.perDay),  # find whale 40-Hz call
    Mn     = as.numeric(m$Mn.perDay),  # humpback whale calls
    Gg     = as.numeric(m$Gg.perDay),  # Risso's dolphin echolocation
    Kspp   = as.numeric(m$Kspp.perDay),  # Kogia spp. echolocation
    Lo     = as.numeric(m$Lo.perDay),  # Pacific white-sided dolphin echolocation
    Mc     = as.numeric(m$Mc.perDay),  # Hubb's' beaked whale echolocation
    Mg     = as.numeric(m$Mg.perDay),   # ginkgo-toothed beaked whale echolocation
    Ms     = as.numeric(m$Ms.perDay),  # Stejneger's beaked whale echolocation
    Oo     = as.numeric(m$Oo.perDay),  # Killer whale echolocation
    Pm     = as.numeric(m$Pm.perDay),  # sperm whale echolocation
    UO     = as.numeric(m$UO.perDay),   # common dolphin proxy echolocation
    Zc     = as.numeric(m$Zc.perDay)  # Goose-beaked whale echolocation
  )

  daily <- data.frame(
    Date = dates,
    JD = as.integer(format(dates, "%j")),
    Period = period
  )
  daily <- cbind(daily, temp, temp_an, psal, psal_an, spice,
                  DO500 = as.numeric(m$O2.509.perDay), #dissolved oxygen at 500 m (509 m sensor)
                  NO3   = as.numeric(m$NO3.40.perDay), #nitrate at 40 m (41 m sensor)
                  Chl40 = as.numeric(m$Chl.41.perDay), #chlorophyll a at 40 m (41 m sensor)
                  FSLE    = as.numeric(m$FSLE), #finite-size Lyapunov exponent (FSLE)
                  FSLEstd = as.numeric(m$FSLEstd), #standard deviation of FSLE
                  sv, nasc, responses)

  attr(daily, "response_vars") <- names(responses)
  attr(daily, "depth_predictor_stems") <- c("T", "TAnom", "S", "SAnom", "Spice")
  attr(daily, "depth_ctd") <- depth_ctd
  attr(daily, "scalar_predictors") <- c("DO500", "NO3", "Chl40", "FSLE", "FSLEstd")
  attr(daily, "nasc_sv_depths") <- wbat_depth

  daily
}

if (sys.nframe() == 0) {
  daily <- load_cce1_daily()
  cat("Loaded", nrow(daily), "rows x", ncol(daily), "cols\n")
  cat("Date range:", format(range(daily$Date)), "\n")
}
