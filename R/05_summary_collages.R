# Per-model summary collages (manuscript Figs. 6 and 7 and Supplementary
# Figs. S6-S8 are assembled from these): top = GAM partial-effect panels
# (dashed 95% CI, significance stars, depth label), bottom = matching
# predictor time series with period shading.
#
# Design notes:
#   - The header reports rho, the AR1 coefficient fitted in
#     00d_gam_pipeline_functions.R: the data are not binned to remove
#     autocorrelation.
#   - Cetacean (non-Gaussian) response panels are weekly-summed bar charts
#     ("cumul./week", as in Fig. 3), because daily detection series are mostly
#     zero with sparse spikes and read poorly as a line. NASC/backscatter
#     (Gaussian) response panels are a continuous daily line on a log scale.

library(mgcv)
library(grid)

source("config.R")
source("00a_load_data.R")
source("00b_load_dvm.R")
source("00c_dvm_plot.R")
source("00d_gam_pipeline_functions.R")
source("00f_model_specs.R")


daily <- load_cce1_daily()
daily <- merge_dvm_into_daily(daily, load_dvm_amplitude(validate = FALSE))
daily <- add_response_columns(daily)
specs <- cetacean_specs()
# DVM series (as response or predictor) are drawn exactly as in
# DVM_amplitude_timeseries_log.pdf: signed bars, symmetric log axis, and one
# y-range shared by all depth layers so collages compare across depths.
DVM_YLIM <- dvm_shared_ylim(daily, log_scale = TRUE)
dvm_collage_panel <- function(v, cex.lab) {
  # two-line label like the other time series panels; the full
  # "DVM amplitude <bin>" label is too long for a 0.6in panel
  ylab <- sub("^DVM amplitude ", "DVM amplitude\n", dvm_bin_label[[v]])
  dvm_panel(daily$Date, daily[[v]], log_scale = TRUE, ylim = DVM_YLIM, ylab = ylab,
            shade = function() shade_periods(show_labels = FALSE, draw_box = FALSE),
            cex.lab = cex.lab, cex.axis = 0.75, cex.axis.y = 0.6, tick_step = 2)
}

responses <- c("NASC5", "NASC55", "NASC155", "NASC255", "NASC355", "NASC455",
               specs$label, paste0("DVM", DVM_DEPTHS))

# NASC<d> response code -> the 50 m depth bin it represents (depth-bin convention)
nasc_bin_label <- c(NASC5 = "sA 0-50 m", NASC55 = "sA 50-100 m", NASC155 = "sA 150-200 m",
                     NASC255 = "sA 250-300 m", NASC355 = "sA 350-400 m", NASC455 = "sA 450-500 m")
dvm_bin_label <- c(DVM5 = "DVM amplitude 0-50 m", DVM55 = "DVM amplitude 50-100 m",
                    DVM155 = "DVM amplitude 150-200 m", DVM255 = "DVM amplitude 250-300 m",
                    DVM355 = "DVM amplitude 350-400 m", DVM455 = "DVM amplitude 450-500 m")
response_label <- function(label) {
  if (label %in% names(nasc_bin_label)) nasc_bin_label[[label]]
  else if (label %in% names(dvm_bin_label)) dvm_bin_label[[label]]
  else specs$display[specs$label == label]
}
spec_part <- function(label) if (label %in% specs$label) specs$part[specs$label == label] else "single"

# variable code -> human-readable stem + unit (prefix-matched)
var_info <- function(v) {
  # NO3 ends in a digit that is part of its NAME, not a depth suffix -- the
  # generic trailing-digit strip below would read it as "depth 3 m" and
  # corrupt stem to "NO" (silently breaking both the depth tag AND the
  # color lookup, since VAR_COLOR is keyed by "NO3"). Handled explicitly,
  # before the generic regex ever sees it. FSLE/FSLEstd have no trailing
  # digit so the generic path already leaves them with depth=NA correctly.
  if (v == "NO3") return(list(depth = NA, stem = "NO3", stem_label = "Nitrate", unit = "umol/L"))

  depth <- regmatches(v, regexpr("[0-9]+$", v))
  depth <- if (length(depth) == 0) NA else depth
  stem <- sub("[0-9]+$", "", v)
  info <- switch(stem,
    "T"      = list(label = "Temperature", unit = "degC"),
    "TAnom"  = list(label = "Temp. anomaly", unit = "degC"),
    "S"      = list(label = "Salinity", unit = "psu"),
    "SAnom"  = list(label = "Sal. anomaly", unit = "psu"),
    "DVM"    = list(label = "DVM amplitude", unit = "m2/nmi2"),
    "NASC"   = list(label = "Backscatter", unit = "m2/nmi2"),
    "Chl"    = list(label = "Chlorophyll a", unit = "ug/L"),
    "DO"     = list(label = "Dissolved oxygen", unit = "ml/L"),   # see DISPLAY_SCALE
    list(label = v, unit = "")
  )
  if (v == "FSLE") info <- list(label = "FSLE", unit = "1/day")
  if (v == "FSLEstd") info <- list(label = "FSLE (SD)", unit = "1/day")
  list(depth = depth, stem = stem, stem_label = info$label, unit = info$unit)
}

# Display-only unit conversions (models are fit in the original units): stem
# -> divisor. DO is umol/L in the data; 1 ml/L O2 = 44.66 umol/L.
DISPLAY_SCALE <- c(DO = 44.66)
display_scale_for <- function(v) {
  stem <- var_info(v)$stem
  if (stem %in% names(DISPLAY_SCALE)) DISPLAY_SCALE[[stem]] else 1
}

xlab_for <- function(v) {
  i <- var_info(v)
  if (i$unit == "") i$stem_label else sprintf("%s (%s)", i$stem_label, i$unit)
}
depth_tag_for <- function(v) {
  i <- var_info(v)
  if (is.na(i$depth)) "" else paste0(i$depth, " m")
}
ts_ylab_for <- function(v) {
  i <- var_info(v)
  tag <- if (is.na(i$depth)) "" else paste0("\n", i$depth, " m")
  paste0(i$stem_label, tag)
}

# one fixed color per variable FAMILY (not per depth -- depth is already
# conveyed by the panel's own label), used consistently across every
# collage so e.g. temperature is always the same color whichever model it
# appears in. Response lines/bars use their own fixed color, kept visually
# distinct from every predictor family below.
CETACEAN_BAR_COLOR <- "#4D4D4D"    # dark grey, cetacean weekly-count bars
VAR_COLOR <- c(T = "#0072B2", TAnom = "#0072B2",      # temperature family: blue
               S = "#C0392B", SAnom = "#C0392B",      # salinity family: reddish
               Chl = "#2E8B57",                        # chlorophyll: green
               NASC = "#006D77",                       # backscatter: dark turquoise
               DO = "#6A3D9A",                          # dissolved oxygen: purple
               NO3 = "#8B4513", FSLE = "#000000", FSLEstd = "#444444")
RESPONSE_LINE_COLOR <- VAR_COLOR[["NASC"]]   # same dark turquoise whether
                                               # backscatter is response or predictor
var_color_for <- function(v) {
  stem <- var_info(v)$stem
  if (stem %in% names(VAR_COLOR)) VAR_COLOR[[stem]] else "black"
}

star_for <- function(p) {
  if (is.na(p)) "" else if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "n.s."
}

# period shading: alternating grey/white bands from the Period column,
# labelled 1-11 at the top of the figure
period_runs <- local({
  rl <- rle(daily$Period)
  ends <- cumsum(rl$lengths)
  starts <- c(1, utils::head(ends, -1) + 1)
  data.frame(period = rl$values, start = daily$Date[starts], end = daily$Date[ends])
})

shade_periods <- function(show_labels = FALSE, draw_box = TRUE) {
  usr <- par("usr")
  for (i in seq_len(nrow(period_runs))) {
    col <- if (period_runs$period[i] %% 2 == 1) "grey88" else "white"
    rect(as.numeric(period_runs$start[i]), usr[3], as.numeric(period_runs$end[i]), usr[4],
         col = col, border = NA)
    if (show_labels) {
      mid <- (as.numeric(period_runs$start[i]) + as.numeric(period_runs$end[i])) / 2
      text(mid, usr[4], period_runs$period[i], xpd = NA, font = 2, cex = 0.65, pos = 3, offset = 0.15)
    }
  }
  if (draw_box) box()
}

# weekly totals for cetacean response bar charts (Fig. 3 convention)
weekly_agg <- function(dates, y) {
  week_idx <- as.integer(floor(as.numeric(dates - min(dates)) / 7))
  val <- tapply(y, week_idx, sum, na.rm = TRUE)
  idx <- as.integer(names(val))
  data.frame(start = min(dates) + idx * 7, end = min(dates) + idx * 7 + 7, value = as.numeric(val))
}

# Sized for two-column assembly: a bit less than half of an 8.5in page width
# (and still comfortably under half of A4's 8.27in). GAM panels use a FIXED
# 2-column grid for every collage, so every partial-effect panel has identical width regardless of
# how many terms a given model has; a 1-term model leaves the second cell of
# its one row blank rather than stretching its single panel to fill the row.
# A dedicated thin period-number row sits between the GAM grid and the time
# series stack so every time series panel (including the response) can use
# the SAME inner margins and come out the same height.
PAGE_WIDTH_IN <- 4.0
GAM_NCOL <- 2
GAM_ROW_IN <- 1.75     # panel height plus room below it for the x-axis title,
                        # given the ~1.4in column width at this page size
TS_ROW_IN  <- 0.6      # height of each time series panel
# The bottom-most time series panel carries a much bigger bottom mar (room
# for the date axis tick labels) than every other row. That extra margin
# takes space from the LAST panel's plot area, so extra row height here keeps
# its plotting box the same size as the others.
TS_LAST_ROW_IN <- TS_ROW_IN + 0.75
PERIOD_ROW_IN <- 0.5    # includes a gap above the period numbers for the
                        # x-axis titles of the last row of effect panels

# depth-sort key for ordering both the GAM-panel grid and the matching time
# series panels: ascending by depth, non-depth variables (e.g. FSLE) last
depth_key <- function(v) {
  d <- var_info(v)$depth
  if (is.na(d)) Inf else as.numeric(d)
}

make_collage <- function(label) {
  model_path <- file.path(MODEL_DIR, paste0(label, "_model.rds"))
  if (!file.exists(model_path)) { message("missing: ", label); return(invisible(NULL)) }
  res <- readRDS(model_path)
  if (is.null(res$final) || is.null(res$final$model)) { message("no final model: ", label); return(invisible(NULL)) }

  m <- res$final$model
  terms <- res$selected                 # original fit order -- matches m$smooth / s$s.table
  n <- length(terms)
  depth_order <- order(vapply(terms, depth_key, numeric(1)))  # display order, shallow -> deep
  nrow_gam <- ceiling(n / GAM_NCOL)
  n_ts <- 1 + n   # response panel + one per selected predictor

  page_height <- 0.9 + nrow_gam * GAM_ROW_IN + PERIOD_ROW_IN +
    (n_ts - 1) * TS_ROW_IN + TS_LAST_ROW_IN

  pdf(file.path(PLOT_DIR, paste0(label, "_collage.pdf")),
      width = PAGE_WIDTH_IN, height = page_height)

  gam_mat <- matrix(seq_len(nrow_gam * GAM_NCOL), nrow = nrow_gam, ncol = GAM_NCOL, byrow = TRUE)
  gam_mat[gam_mat > n] <- 0   # blank unused cell(s) in the last row

  period_row_id <- n + 1
  ts_ids <- (n + 2):(n + 1 + n_ts)
  period_mat <- matrix(period_row_id, nrow = 1, ncol = GAM_NCOL)
  ts_mat <- matrix(rep(ts_ids, each = GAM_NCOL), nrow = n_ts, ncol = GAM_NCOL, byrow = TRUE)

  full_mat <- rbind(gam_mat, period_mat, ts_mat)
  heights <- c(rep(GAM_ROW_IN, nrow_gam), PERIOD_ROW_IN,
               rep(TS_ROW_IN, n_ts - 1), TS_LAST_ROW_IN)
  layout(full_mat, heights = heights)
  par(oma = c(4, 4.5, 2.2, 1))

  # --- GAM partial-effect panels (shown shallow -> deep via depth_order) -----
  s <- summary(m)
  for (j in seq_len(n)) {
    i <- depth_order[j]        # i = index into the model's own term order
    v <- terms[i]
    # bottom margin deep enough that the x-axis title (drawn at line 3) is
    # not covered by the row below
    par(mar = c(4.2, 4.6, 2.2, 1.2))
    sc <- display_scale_for(v)
    # scale = 0: each panel gets its own y-axis (mgcv's default shares one
    # across panels, which lets one steep or poorly constrained term flatten
    # all the others). The y label is mgcv's default, s(variable, edf); what
    # the y values mean (centred effect on the link scale) goes in the
    # figure caption -- see the methods write-up.
    plot(m, select = i, se = TRUE, shade = FALSE, rug = TRUE, scale = 0,
         xlab = xlab_for(v), main = "", cex.lab = 0.85, cex.axis = 0.75,
         xaxt = if (sc != 1) "n" else "s")
    if (sc != 1) {   # model x is in original units; label ticks in display units
      lim <- par("usr")[1:2] / sc
      at <- pretty(lim, n = 3)
      at <- at[at >= lim[1] & at <= lim[2]]
      axis(1, at = at * sc, labels = at, cex.axis = 0.75)
    }
    usr <- par("usr")
    tag <- depth_tag_for(v)
    if (nzchar(tag)) text(usr[1], usr[4], tag, adj = c(-0.1, 1.4), font = 2, cex = 0.8, xpd = NA)
    text(usr[2], usr[4], star_for(s$s.table[i, "p-value"]), adj = c(1.1, 1.4),
         font = 2, cex = 0.85, xpd = NA)
  }

  mtext(sprintf("R-sq.(adj) = %.2f    Dev. expl. = %.0f%%    rho = %.2f",
                s$r.sq, model_dev_expl(m) * 100, res$final$rho),
        outer = TRUE, side = 3, line = 0.5, cex = 0.85, font = 2)

  # --- dedicated period-number row (no data, just the 1-11 labels + shading) --
  # draw_box=FALSE: with bty="n" above AND a box from shade_periods(), this
  # row rendered as a bordered box containing only numbers -- reading as an
  # empty panel. No box here; the row below's box reads as the top border.
  par(mar = c(0, 5, 1.6, 1.5))
  plot(daily$Date, rep(0, nrow(daily)), type = "n", xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", bty = "n")
  shade_periods(show_labels = TRUE, draw_box = FALSE)

  # --- response time series --------------------------------------------------
  is_gaussian <- identical(res$family, "gaussian")
  is_dvm <- grepl("^DVM", label)   # Gaussian but signed -- plain log scale would break it
  yvar <- daily[[res$y]]
  par(mar = c(0.3, 5, 0.3, 1.5))
  part <- spec_part(label)
  if (is_dvm) {
    dvm_collage_panel(res$y, cex.lab = 0.75)
  } else if (part == "when present") {
    # the modelled amount of calling, on the original scale, on detection days
    sp <- specs[specs$label == label, ]
    orig <- daily[[sp$species]]
    present <- !is.na(orig) & orig > 0
    plot(daily$Date, orig, type = "n", xaxt = "n", xlab = "",
         ylab = paste0(pos_unit[[two_part[[sp$species]]]], "\n(when present)"),
         ylim = c(0, max(orig, na.rm = TRUE) * 1.05), cex.lab = 0.75, cex.axis = 0.75)
    shade_periods(show_labels = FALSE)
    segments(as.numeric(daily$Date[present]), 0, as.numeric(daily$Date[present]), orig[present],
             col = CETACEAN_BAR_COLOR, lwd = 0.8)
    box()
  } else if (part == "occurrence") {
    wk <- weekly_agg(daily$Date, yvar)
    plot(daily$Date, rep(0, nrow(daily)), type = "n", xaxt = "n", xlab = "",
         ylab = "Days present\nper week", ylim = c(0, 7.3), cex.lab = 0.75, cex.axis = 0.75)
    shade_periods(show_labels = FALSE)
    rect(as.numeric(wk$start), 0, as.numeric(wk$end), wk$value,
         col = CETACEAN_BAR_COLOR, border = NA)
    box()
  } else if (is_gaussian) {
    plot(daily$Date, yvar, type = "n", xaxt = "n", xlab = "",
         ylab = response_label(label), log = "y", cex.lab = 0.8, cex.axis = 0.75)
    shade_periods(show_labels = FALSE)
    lines(daily$Date, yvar, col = RESPONSE_LINE_COLOR, lwd = 1.2)
  } else {
    wk <- weekly_agg(daily$Date, yvar)
    plot(daily$Date, rep(0, nrow(daily)), type = "n", xaxt = "n", xlab = "",
         ylab = response_label(label), ylim = c(0, max(wk$value, na.rm = TRUE) * 1.05),
         cex.lab = 0.8, cex.axis = 0.75)
    shade_periods(show_labels = FALSE)
    rect(as.numeric(wk$start), 0, as.numeric(wk$end), wk$value,
         col = CETACEAN_BAR_COLOR, border = NA)
    box()
  }

  # --- predictor time series (same shallow -> deep order as the GAM panels) --
  for (j in seq_len(n)) {
    i <- depth_order[j]
    v <- terms[i]
    last <- j == n
    par(mar = c(if (last) 4.2 else 0.3, 5, 0.3, 1.5))
    pv <- daily[[v]]
    if (grepl("^DVM", v)) {
      dvm_collage_panel(v, cex.lab = 0.75)
      if (last) {
        axis.Date(1, x = daily$Date, cex.axis = 0.75)
        title(xlab = "Date", cex.lab = 0.75)
      }
      next
    }
    sc <- display_scale_for(v)
    pv <- pv / sc
    is_nasc_var <- grepl("^NASC", v)
    plot(daily$Date, pv, type = "n", xaxt = if (last) "s" else "n",
         xlab = if (last) "Date" else "", log = if (is_nasc_var) "y" else "",
         yaxt = if (sc != 1) "n" else "s",
         ylab = ts_ylab_for(v), cex.lab = 0.75, cex.axis = 0.75)
    if (sc != 1) {   # converted units: coarser ticks, as on the GAM panel
      lim <- par("usr")[3:4]
      at <- pretty(lim, n = 3)
      axis(2, at = at[at >= lim[1] & at <= lim[2]], cex.axis = 0.75)
    }
    shade_periods(show_labels = FALSE)
    lines(daily$Date, pv, col = var_color_for(v), lwd = 1.1)
  }

  # left sidebar label (grey bar + rotated text), drawn last so it overlays cleanly
  grid::grid.rect(x = unit(0.015, "npc"), y = unit(0.5, "npc"),
                   width = unit(0.03, "npc"), height = unit(0.98, "npc"),
                   gp = gpar(fill = "grey75", col = NA))
  grid::grid.text(response_label(label), x = unit(0.015, "npc"), y = unit(0.5, "npc"),
                   rot = 90, gp = gpar(fontface = "bold", cex = 0.9))

  dev.off()
  message(label, ": collage saved (", n, " terms)")
}

for (r in responses) {
  tryCatch(make_collage(r), error = function(e) message("ERROR on ", r, ": ", conditionMessage(e)))
}

cat("\nCollages written to:", PLOT_DIR, "\n")
