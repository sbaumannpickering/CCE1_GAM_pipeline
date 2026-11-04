# Shared drawing code for DVM-amplitude time series, so the standalone
# figure (06_dvm_timeseries_figure.R) and the model collages
# (05_summary_collages.R) render DVM identically: daily bars colored by sign,
# optional symmetric log y-axis, optional y-limits shared across depths.

DVM_POS_COL <- "#2166AC"   # more backscatter at night (migrators ascending)
DVM_NEG_COL <- "#B2182B"   # more backscatter by day

# Amplitude is signed, so the log variant uses a symmetric log transform,
# sign(x) * log10(1 + |x|), which keeps negative days visible and is ~linear
# near zero. Axis ticks are labelled in the original NASC units.
symlog <- function(x) sign(x) * log10(1 + abs(x))

# Symmetric y-limits (in plotted units) covering every DVM depth layer.
dvm_shared_ylim <- function(daily, log_scale = TRUE, depths = DVM_DEPTHS) {
  tr <- if (log_scale) symlog else identity
  all_v <- unlist(lapply(depths, function(d) tr(daily[[paste0("DVM", d)]])))
  max(abs(all_v), na.rm = TRUE) * c(-1, 1)
}

# One DVM panel. x is raw night-minus-day NASC; ylim is in plotted units
# (NULL = symmetric limits from this series alone). shade is called after the
# empty plot is set up, to draw period bands underneath the bars. tick_step
# labels every Nth power of ten (counting down from the largest), for short
# panels where every decade would crowd.
dvm_panel <- function(dates, x, log_scale = TRUE, ylim = NULL, ylab = "",
                      shade = function() NULL, cex.lab = 0.8, cex.axis = 0.75,
                      cex.axis.y = cex.axis - 0.1, tick_step = 1) {
  v <- if (log_scale) symlog(x) else x
  if (is.null(ylim)) ylim <- max(abs(v), na.rm = TRUE) * c(-1, 1)
  plot(dates, v, type = "n", ylim = ylim, xaxt = "n",
       yaxt = if (log_scale) "n" else "s",
       xlab = "", ylab = ylab, cex.lab = cex.lab, cex.axis = cex.axis)
  shade()
  if (log_scale) {
    e <- 0:ceiling(ylim[2])
    e <- e[symlog(10^e) <= ylim[2]]
    labelled <- (max(e) - e) %% tick_step == 0
    mags <- 10^e[labelled]
    ticks <- c(-rev(mags), 0, mags)
    axis(2, at = symlog(ticks), labels = format(ticks, big.mark = ",", scientific = FALSE, trim = TRUE),
         las = 1, cex.axis = cex.axis.y)
    # skipped decades still get an unlabelled tick (except +/-1, too close to 0)
    minor <- 10^e[!labelled & e >= 1]
    if (length(minor)) axis(2, at = symlog(c(-minor, minor)), labels = FALSE)
  }
  bar_hw <- 0.5   # half-width of each daily bar, in days
  pos <- !is.na(v) & v >= 0
  neg <- !is.na(v) & v < 0
  if (any(pos)) rect(as.numeric(dates[pos]) - bar_hw, 0, as.numeric(dates[pos]) + bar_hw, v[pos],
                     col = DVM_POS_COL, border = NA)
  if (any(neg)) rect(as.numeric(dates[neg]) - bar_hw, 0, as.numeric(dates[neg]) + bar_hw, v[neg],
                     col = DVM_NEG_COL, border = NA)
  abline(h = 0, lwd = 0.8)
  box()
}
