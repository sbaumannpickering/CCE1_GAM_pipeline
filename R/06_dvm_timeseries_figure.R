# Standalone water-column DVM-amplitude figure (Supplementary Fig. S5). One
# panel per depth layer (shallow -> deep), each showing night-minus-day
# backscatter for every day of the deployment, colored by sign: positive
# (more backscatter at night, consistent with migrators ascending into the
# layer) in blue, negative (more by day) in red. Period shading and styling
# match the model collages (05_summary_collages.R).

source("config.R")
source("00a_load_data.R")
source("00b_load_dvm.R")
source("00c_dvm_plot.R")


daily <- load_cce1_daily()
dvm <- load_dvm_amplitude(validate = TRUE)
daily <- merge_dvm_into_daily(daily, dvm)

period_runs <- local({
  rl <- rle(daily$Period)
  ends <- cumsum(rl$lengths)
  starts <- c(1, utils::head(ends, -1) + 1)
  data.frame(period = rl$values, start = daily$Date[starts], end = daily$Date[ends])
})

shade_periods <- function(show_labels = FALSE) {
  usr <- par("usr")
  for (i in seq_len(nrow(period_runs))) {
    col <- if (period_runs$period[i] %% 2 == 1) "grey88" else "white"
    rect(as.numeric(period_runs$start[i]), usr[3], as.numeric(period_runs$end[i]), usr[4],
         col = col, border = NA)
    if (show_labels) {
      mid <- (as.numeric(period_runs$start[i]) + as.numeric(period_runs$end[i])) / 2
      text(mid, usr[4], period_runs$period[i], xpd = NA, font = 2, cex = 0.7, pos = 3, offset = 0.15)
    }
  }
}

depths <- c(5, 55, 155, 255, 355, 455)
depth_bin_label <- c(`5` = "0-50 m", `55` = "50-100 m", `155` = "150-200 m",
                      `255` = "250-300 m", `355` = "350-400 m", `455` = "450-500 m")
plot_dvm_timeseries <- function(file, log_scale = FALSE, shared_ylim = FALSE) {
  ylim <- if (shared_ylim) dvm_shared_ylim(daily, log_scale, depths) else NULL
  pdf(file, width = 7, height = 9)
  # Identical margins on every panel (date axis and period labels go in the
  # outer margin) so all plot areas are the same height.
  par(mfrow = c(length(depths), 1), oma = c(4, 1, 4, 1), mar = c(0.4, 6, 0.4, 1.5))

  for (i in seq_along(depths)) {
    d <- depths[i]
    first <- i == 1
    last <- i == length(depths)
    dvm_panel(daily$Date, daily[[paste0("DVM", d)]], log_scale = log_scale, ylim = ylim,
              ylab = paste0(depth_bin_label[[as.character(d)]],
                            if (log_scale) "\nNASC night-day\n(symmetric log)" else "\nNASC night-day"),
              shade = function() shade_periods(show_labels = first))
    if (last) axis.Date(1, x = daily$Date, cex.axis = 0.75, xpd = NA)
  }

  mtext(paste0("Diel vertical migration amplitude (NASC night - day) by depth layer",
               if (log_scale) ", symmetric log scale" else ""),
        outer = TRUE, side = 3, line = 2, font = 2, cex = 1.0)
  mtext("Date", outer = TRUE, side = 1, line = 2.5, cex = 0.85)

  dev.off()
  cat("Saved: ", file, "\n")
}

plot_dvm_timeseries(file.path(PLOT_DIR, "DVM_amplitude_timeseries.pdf"))
plot_dvm_timeseries(file.path(PLOT_DIR, "DVM_amplitude_timeseries_log.pdf"), log_scale = TRUE, shared_ylim = TRUE)
