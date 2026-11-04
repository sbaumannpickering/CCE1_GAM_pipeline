source("config.R")

# OMP closure-error time series (Supplementary Fig. S2, Table S2): the
# uncertainty behind the water-mass fractions in Fig. 2 (fractions summing to
# more than 100% indicate OMP-analysis uncertainty). Source files are the
# minimum-error OMP output (temperature + salinity only, three water masses:
# PSUW, PEW, ENPCW; climatological T-S end-members from Bograd et al. 2019) at
# the three mooring depths. At each depth, three candidate end-member triads
# (upper/center/lower levels of the water-mass T-S envelopes) were tested and
# the minimum-error triad kept: upper level at 75 m, center level at 150 m,
# lower level at 301 m. "error" is (fPSUW + fPEW + fENPCW) - 1, i.e. how far
# the three fractions miss summing to 1.

PLOT_DIR <- file.path(RESULTS_DIR, "omp_uncertainty")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

depth_files <- c(`75` = "x_error_depth75_id1004.csv",
                  `150` = "x_error_depth150_id1003.csv",
                  `301` = "x_error_depth301_id1002.csv")
depth_endmember_label <- c(`75` = "upper end-members (min. error)",
                            `150` = "center end-members (min. error)",
                            `301` = "lower end-members (min. error)")

load_omp_error <- function(depth) {
  d <- read.csv(file.path(OMP_DIR, depth_files[[depth]]))
  d$Date <- as.Date(sprintf("%04d-%02d-%02d", d$year, d$month, d$day))
  d[, c("Date", "error")]
}

depths <- names(depth_files)
omp <- lapply(depths, load_omp_error)
names(omp) <- depths

summary_tbl <- do.call(rbind, lapply(depths, function(z) {
  e <- omp[[z]]$error
  data.frame(depth_m = z, n = length(e),
             mean_pct = round(mean(e) * 100, 1),
             median_pct = round(median(e) * 100, 1),
             sd_pct = round(sd(e) * 100, 1),
             min_pct = round(min(e) * 100, 1),
             max_pct = round(max(e) * 100, 1))
}))
write.csv(summary_tbl, file.path(PLOT_DIR, "omp_error_summary.csv"), row.names = FALSE)
print(summary_tbl, row.names = FALSE)

ERR_COL <- "#006D77"  # same teal used for NASC/backscatter elsewhere in the figure set

pdf(file.path(PLOT_DIR, "OMP_closure_error_timeseries.pdf"), width = 7, height = 7.5)
par(mfrow = c(length(depths), 1), oma = c(4, 1, 3, 1), mar = c(0.4, 6, 0.4, 1.5))

for (i in seq_along(depths)) {
  z <- depths[i]
  first <- i == 1
  last <- i == length(depths)
  d <- omp[[z]]
  par(mar = c(if (last) 3.5 else 0.4, 6, if (first) 1.8 else 0.4, 1.5))
  plot(d$Date, d$error * 100, type = "n", ylim = c(0, max(d$error * 100, na.rm = TRUE) * 1.08),
       xaxt = if (last) "s" else "n",
       xlab = "", ylab = paste0(z, " m\nerror (%)"), cex.lab = 0.85, cex.axis = 0.75)
  lines(d$Date, d$error * 100, col = ERR_COL, lwd = 1.2)
  mtext(depth_endmember_label[[z]], side = 3, line = -1.1, adj = 0.02, cex = 0.65, font = 3)
  s <- summary_tbl[summary_tbl$depth_m == z, ]
  abline(h = s$mean_pct, lty = 2, col = "grey40", lwd = 0.8)
  box()
}

mtext("OMP closure error by mooring depth (fPSUW + fPEW + fENPCW - 1), minimum-error end-member set",
      outer = TRUE, side = 3, line = 1.2, font = 2, cex = 0.95)
mtext("Date", outer = TRUE, side = 1, line = 2.5, cex = 0.85)

dev.off()
cat("Saved: ", file.path(PLOT_DIR, "OMP_closure_error_timeseries.pdf"), "\n")
cat("Saved: ", file.path(PLOT_DIR, "omp_error_summary.csv"), "\n")
