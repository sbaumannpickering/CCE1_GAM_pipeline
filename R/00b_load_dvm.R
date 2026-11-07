# DVM (diel vertical migration) amplitude from the hourly WBAT export.
#
# Source: CCE1_WBAT_hourly.mat. NASC.All is hourly backscatter in 12 bins of 50 m.
# The file is ordered SURFACE-FIRST, like CCE1_WBAT_HARP_daily.mat: column 1 is
# the shallowest bin and `depth` gives each bin's top and bottom in m below the
# surface (the transducer is upward-looking at ~504 m). Averaging each column by calendar day reproduces
# NASC.perDay in CCE1_WBAT_HARP_daily.mat (median relative difference ~1e-16),
# which fixes the mapping:
#
#   NASC bin (code)  column  depth bin (m)
#   NASC5               2       0-50
#   NASC55              3       50-100
#   NASC155             5       150-200
#   NASC255             7       250-300
#   NASC355             9       350-400
#   NASC455            11       450-500
#
# Column 1 (-50 to 0 m) is the sea-surface reflection and column 12 (500-550 m)
# is transducer ringdown; both are excluded. The response codes (5, 55, ...,
# 455) are internal identifiers; display labels use the depth bins above.
#
# DVM is computed at native hourly resolution and aggregated into
# day-mean/night-mean using the sunrise/sunset windows (day/night matrices)
# stored in the same file. `validate` reproduces the NASC.day / NASC.night
# split from columns 2,4,6,8,10 (nominal bin tops 5-405 m, shallow to deep);
# this checks the window logic only, not depth labels.

source("config.R")
library(R.matlab)


# NASC.All column index -> response code (mapping above)
NASC_ALL_COL_TO_BIN <- c(`2` = "NASC5", `3` = "NASC55", `5` = "NASC155",
                          `7` = "NASC255", `9` = "NASC355", `11` = "NASC455")

#' day/night-window mean of NASC.All for one depth column, for every day
#' index in `windows` (a 385x2 matrix of [start,end] MATLAB datenums)
window_means <- function(hourly_dates, hourly_vals, windows) {
  vapply(seq_len(nrow(windows)), function(i) {
    in_window <- hourly_dates >= windows[i, 1] & hourly_dates <= windows[i, 2]
    if (!any(in_window)) return(NA_real_)
    mean(hourly_vals[in_window], na.rm = TRUE)
  }, numeric(1))
}

load_dvm_amplitude <- function(mat_path = HOURLY_MAT, validate = TRUE) {
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
    # reproduce the stored 5-bin NASC.day/NASC.night exactly to confirm the
    # day/night window logic (columns only; says nothing about depth labels)
    old_cols <- c(2, 4, 6, 8, 10)   # NASC.day/NASC.night's 5 columns, in order (shallow to deep)
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
            "from columns 2,4,6,8,10 (window logic only).")
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
