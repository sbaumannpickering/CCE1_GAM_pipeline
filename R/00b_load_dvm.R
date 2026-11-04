# DVM (diel vertical migration) amplitude from the raw hourly WBAT export.
#
# Source: CCE1_SB_01_all.mat. NASC.All is hourly backscatter in 12 fine bins
# binned by RANGE FROM THE UPWARD-LOOKING TRANSDUCER (~504 m), not by depth, so
# column order is inverted relative to true depth (true depth ~ 504 - range).
# Averaging each raw column by calendar day reproduces NASC.perDay in
# combined_diel_daily.mat (median relative difference 0, 385/386 days), which
# fixes the mapping:
#
#   NASC bin (code)  raw column  depth bin (m)
#   NASC5               11          0-50
#   NASC55              10          50-100
#   NASC155              8          150-200
#   NASC255              6          250-300
#   NASC355              4          350-400
#   NASC455              2          450-500
#
# Column 1 (range -45 to 5 m) is transducer ringdown at ~504 m and column 12
# (range 505-555 m) is the sea-surface reflection; both are excluded. The
# response codes (5, 55, ..., 455) are internal identifiers; display labels
# use the depth bins above.
#
# DVM is computed at native hourly resolution and aggregated into
# day-mean/night-mean using the sunrise/sunset windows (day/night matrices)
# stored in combined_diel_daily.mat. `validate` reproduces the NASC.day /
# NASC.night split from raw columns 3,5,7,9,11; this checks the window logic
# only, not depth labels.

source("config.R")
library(R.matlab)

WBAT_EXPORT_PATH <- WBAT_EXPORT_MAT

# NASC.All column index -> response code (mapping above)
NASC_ALL_COL_TO_BIN <- c(`11` = "NASC5", `10` = "NASC55", `8` = "NASC155",
                          `6` = "NASC255", `4` = "NASC355", `2` = "NASC455")

#' day/night-window mean of NASC.All for one depth column, for every day
#' index in `windows` (a 385x2 matrix of [start,end] MATLAB datenums)
window_means <- function(hourly_dates, hourly_vals, windows) {
  vapply(seq_len(nrow(windows)), function(i) {
    in_window <- hourly_dates >= windows[i, 1] & hourly_dates <= windows[i, 2]
    if (!any(in_window)) return(NA_real_)
    mean(hourly_vals[in_window], na.rm = TRUE)
  }, numeric(1))
}

load_dvm_amplitude <- function(mat_path = WBAT_EXPORT_PATH, validate = TRUE) {
  m <- readMat(mat_path)
  hourly_dates <- as.numeric(m$dates.All)
  nasc_all <- m$NASC.All
  day_win <- m$day      # 385 x 2: [sunrise, sunset] per day index
  night_win <- m$night  # 385 x 2: [preceding sunset, this day's sunrise]

  # as.Date() on a non-integer numeric keeps the fractional part internally
  # even though format() only ever displays the calendar date -- two dates
  # built this way can look identical when printed yet compare unequal
  # (17087 vs 17087.59). floor() before constructing the Date avoids that;
  # 00a_load_data.R's full_day is already integer so this only matters here,
  # where the anchor is a sunrise TIMESTAMP (fractional day).
  matlab2date <- function(x) as.Date(floor(x) - 719529, origin = "1970-01-01")
  n_days <- nrow(day_win)
  # day/night windows are indexed to the sunrise date; anchor each row to
  # that calendar date via the sunrise timestamp (day_win[,1])
  dates <- matlab2date(day_win[, 1])

  if (validate) {
    # reproduce the OLD 5-bin NASC.day/NASC.night exactly to confirm the
    # day/night window logic (raw columns only; says nothing about depth labels)
    old_cols <- c(3, 5, 7, 9, 11)   # NASC.day/NASC.night's 5 columns, in order
    for (k in seq_along(old_cols)) {
      col <- old_cols[k]
      day_mean <- window_means(hourly_dates, nasc_all[, col], day_win)
      night_mean <- window_means(hourly_dates, nasc_all[, col], night_win)
      day_ok <- isTRUE(all.equal(day_mean, as.numeric(m$NASC.day[, k]), tolerance = 1e-6))
      night_ok <- isTRUE(all.equal(night_mean, as.numeric(m$NASC.night[, k]), tolerance = 1e-6))
      if (!day_ok || !night_ok) {
        stop(sprintf("DVM validation FAILED for NASC.All col %d (day_ok=%s, night_ok=%s)",
                      col, day_ok, night_ok))
      }
    }
    message("DVM validation passed: reproduced existing NASC.day/NASC.night exactly ",
            "from raw columns 3,5,7,9,11 (window logic only).")
  }

  out <- data.frame(Date = dates)
  for (col_str in names(NASC_ALL_COL_TO_BIN)) {
    col <- as.integer(col_str)
    bin <- NASC_ALL_COL_TO_BIN[[col_str]]
    day_mean <- window_means(hourly_dates, nasc_all[, col], day_win)
    night_mean <- window_means(hourly_dates, nasc_all[, col], night_win)
    eps <- 0.01 * stats::median(nasc_all[, col], na.rm = TRUE)   # floor for log-ratio
    out[[paste0(bin, "_day")]] <- day_mean
    out[[paste0(bin, "_night")]] <- night_mean
    out[[paste0("DVM", sub("^NASC", "", bin))]] <- night_mean - day_mean
    out[[paste0("DVMlog", sub("^NASC", "", bin))]] <-
      log((night_mean + eps) / (day_mean + eps))
  }
  out
}

#' left-join DVM columns onto the main daily dataset by Date
merge_dvm_into_daily <- function(daily, dvm = load_dvm_amplitude(validate = FALSE)) {
  merge(daily, dvm[, setdiff(names(dvm), grep("_day$|_night$", names(dvm), value = TRUE))],
        by = "Date", all.x = TRUE)
}

DVM_DEPTHS <- c(5, 55, 155, 255, 355, 455)

if (sys.nframe() == 0) {
  dvm <- load_dvm_amplitude()
  cat("rows:", nrow(dvm), " cols:", ncol(dvm), "\n")
  cat("date range:", format(range(dvm$Date)), "\n")
  print(head(dvm[, c("Date", "DVM5", "DVM255", "DVMlog5", "DVMlog255")]))
  cat("\nsummary of DVM (difference) columns:\n")
  print(summary(dvm[, grep("^DVM[0-9]", names(dvm))]))
}
