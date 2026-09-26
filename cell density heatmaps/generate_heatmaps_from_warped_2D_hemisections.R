###############################################################################
# Average cell-density heatmap from registered (warped) sections
# -----------------------------------------------------------------------------
# Purpose
#   Takes a set of sections that have already been registered/warped into a
#   common coordinate frame (e.g. one representative individual, or a group
#   template built with ANTsPy or similar), and computes:
#     - a kernel density estimate (KDE) of cell positions for EACH section,
#       evaluated on one shared template window, then
#     - the AVERAGE across sections (raw and normalised versions),
#   plus diagnostic figures and a .rds file with everything needed downstream
#   (e.g. by a companion script that renders multiple groups on one shared
#   colour scale).
#
#   This script processes ONE GROUP at a time. To process several groups
#   (e.g. different ages/conditions/species, or different anatomical
#   subregions), edit the CONFIG block below for each group in turn and
#   re-run the script once per group. See "HOW TO USE" at the bottom.
#
# Why average on a shared template window (rather than each section's own
# outline)
#   - Averaging each section's density over its OWN outline and combining
#     with na.rm gives patchy coverage near the edges and inconsistent edge
#     correction from section to section.
#   - Here every section's density is evaluated on ONE shared window, so
#     coverage is uniform (every interior pixel reflects all N sections) and
#     edge correction is identical across sections.
#
# Why a single global y-flip (if your data needs flipping at all)
#   - Flipping each section's y-axis using that section's OWN max(y) can
#     reintroduce a vertical jitter between sections that were otherwise
#     well registered. This script (if flip_y = TRUE) flips every section by
#     ONE shared reference value instead, preserving your registration.
#
# Required R packages: spatstat
###############################################################################

library(spatstat)

## ============================================================================
## ---- CONFIG: EDIT THESE FOR EACH GROUP, THEN RUN THE WHOLE SCRIPT ----------
## ============================================================================

# ---- Input data ---------------------------------------------------------
# Two folders of matched, registered per-section CSVs (they can be the same
# folder, as in the example below):
#   - outline_dir : folder containing one section-outline CSV per section,
#                   named "<outline_prefix><section_id>.csv"
#   - points_dir  : folder containing one cell-centroid CSV per section,
#                   named "<points_prefix><section_id>.csv"
# Each outline CSV needs the columns given by xcol/ycol below, tracing the
# section boundary as an ordered polygon (open or closed).
# Each points CSV needs the same xcol/ycol columns, one row per detected
# cell centroid, in the SAME registered coordinate frame as the outlines.
outline_dir <- "PATH/TO/YOUR/PROJECT/GroupA/warped_csvs"   # folder with outline_*.csv
points_dir  <- "PATH/TO/YOUR/PROJECT/GroupA/warped_csvs"   # folder with registered_centroids_*.csv
outline_prefix <- "outline_"                # filename prefix for outline CSVs
points_prefix  <- "registered_centroids_"    # filename prefix for centroid CSVs

group_name <- "GroupA"                       # used to build all output filenames
xcol <- "x_um"; ycol <- "y_um"               # column names for coordinates (in microns)

# ---- Orientation -----------------------------------------------------------
flip_y <- TRUE        # TRUE:  flip every section vertically using ONE shared
                      #        y reference (see FIX 1 above).
                      # FALSE: use the coordinates exactly as provided.

# ---- Shared template window --------------------------------------------
# Every section's density is evaluated on ONE template window so that all
# sections contribute uniformly and edge correction is consistent.
#   - Point this at a specific section's outline if you warped everything to
#     one representative individual (recommended when you have one):
template_outline_file <- NULL   # e.g. file.path(outline_dir, "outline_example_id.csv"), or NULL
#   - If NULL, a consensus window is built instead: the region covered by at
#     least `template_frac` of all the section outlines (majority vote).
template_frac <- 0.5            # fraction of sections that must cover a pixel
                                # for it to be included in the consensus window

# ---- KDE bandwidth (smoothing) ----------------------------------------
# Two modes:
#   "absolute" - sigma is a fixed number of MICRONS, the same for every group.
#   "relative" - sigma = sigma_frac * sqrt(template window area), i.e. a fixed
#                FRACTION of the structure's size. Use this when comparing
#                groups of very different physical size (e.g. larval vs adult)
#                so every group gets matched relative smoothing/detail.
sigma_mode <- "absolute"    # "absolute" or "relative"
sigma      <- 4             # microns; used when sigma_mode == "absolute"
sigma_frac <- 0.01           # fraction of sqrt(template area); used when
                             # sigma_mode == "relative"
# Tip for matching a look across modes: run one group in "absolute" mode,
# note the "sqrt(area)" value it prints, then set
#   sigma_frac <- your_chosen_sigma / that_printed_sqrt_area
# and switch every group to "relative" mode with that one sigma_frac.

grid_dimyx    <- 512   # pixel resolution of the density grid (per side)
area_unit_um2 <- 100    # density is reported as cells per this many um^2
area_unit_lab <- "100 um^2"   # matching human-readable label for figure axes

# ---- Output location -----------------------------------------------------
# Recommended convention: one output folder per group, nested under a folder
# named for the sigma used, so it's obvious later which files share a
# bandwidth (only files with the SAME sigma / area_unit_um2 should be
# combined or compared downstream).
out_dir <- file.path("PATH/TO/YOUR/PROJECT/AverageCells",
                      paste0("Sigma_", sigma), group_name)

# ---- Figure appearance ----------------------------------------------------
fig_dpi <- 600; fig_w_in <- 6; fig_h_in <- 5; tiff_compress <- "lzw"
sqrt_palette <- "Plasma"; n_cols <- 255
raw_clip_quantile <- 0.97   # raw-intensity figure: colour white point is set
                            # at this percentile of the map's own values
                            # (this group's own diagnostic figure only - the
                            # cross-group SHARED colour scale, if you want one,
                            # is handled by the companion multi-group script)

# ---- Scale bar (drawn on each figure) --------------------------------------
scalebar_len <- 100; scalebar_label <- "100 µm"
scalebar_col <- "black"; scalebar_halo <- "white"
scalebar_lwd <- 3; scalebar_cex <- 0.9

## ============================================================================
## ---- END CONFIG -------------------------------------------------------------
## ============================================================================

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- Helpers ----------------------------------------------------------------

# Flipped, cleaned, correctly-oriented polygon coords for one outline df
outline_xy <- function(df, yref) {
  x <- df[[xcol]]
  y <- if (flip_y) yref - df[[ycol]] else df[[ycol]]
  n <- length(x)
  if (n > 1 && isTRUE(all.equal(x[1], x[n])) && isTRUE(all.equal(y[1], y[n]))) {
    x <- x[-n]; y <- y[-n]                       # drop repeated closing vertex
  }
  nx <- c(x[-1], x[1]); ny <- c(y[-1], y[1])
  if (0.5 * sum(x*ny - nx*y) < 0) { x <- rev(x); y <- rev(y) }  # force anticlockwise
  list(x = x, y = y)
}
make_win <- function(xy) owin(poly = list(x = xy$x, y = xy$y))

flip_points <- function(df, yref) {
  list(x = df[[xcol]], y = if (flip_y) yref - df[[ycol]] else df[[ycol]])
}

add_scalebar <- function(im, len = scalebar_len, label = scalebar_label,
                         col = scalebar_col, halo = scalebar_halo,
                         lwd = scalebar_lwd, cex = scalebar_cex) {
  R <- as.rectangle(im); w <- diff(R$xrange); h <- diff(R$yrange)
  padx <- 0.07*w; pady <- 0.07*h
  x1 <- R$xrange[2]-padx; x0 <- x1-len; yb <- R$yrange[1]+pady; yt <- yb+0.03*h
  segments(x0, yb, x1, yb, col = halo, lwd = lwd+2, lend = "butt")
  segments(x0, yb, x1, yb, col = col,  lwd = lwd,   lend = "butt")
  xm <- (x0+x1)/2; dx <- 0.004*w; dy <- 0.004*h
  for (ox in c(-1,0,1)) for (oy in c(-1,0,1))
    if (ox!=0||oy!=0) text(xm+ox*dx, yt+oy*dy, label, col=halo, cex=cex, adj=c(0.5,0))
  text(xm, yt, label, col = col, cex = cex, adj = c(0.5,0)); invisible(NULL)
}

sqrt_cmap <- function(zmax, pal = sqrt_palette, n = n_cols) {
  br <- seq(0, sqrt(zmax * (1 + 1e-6)), length.out = n + 1)^2
  colourmap(hcl.colors(n, pal), breaks = br)
}

# linear colourmap with the white point clipped at a chosen value, saturating
# above it (gamma 1, clip at the raw_clip_quantile percentile)
gamma_cmap <- function(clip, truemax, pal = sqrt_palette, n = n_cols) {
  br <- seq(min(0, min(as.matrix(avg_raw), na.rm = TRUE)) - 1e-9, clip, length.out = n + 1)
  br[n + 1] <- max(truemax, clip) * (1 + 1e-6)
  colourmap(hcl.colors(n, pal), breaks = br)
}

## ---- Discover + read matched outline/points pairs ---------------------------

outline_files <- list.files(outline_dir,
                            pattern = paste0("^", outline_prefix, ".*\\.csv$"),
                            full.names = TRUE)
if (!length(outline_files)) stop("No outline files found in outline_dir.")
ids <- sub("\\.csv$", "", sub(paste0("^", outline_prefix), "", basename(outline_files)))

outlines <- list(); points <- list()
for (i in seq_along(ids)) {
  id <- ids[i]
  ppath <- file.path(points_dir, paste0(points_prefix, id, ".csv"))
  if (!file.exists(ppath)) { warning("No points for ", id, " - skipping"); next }
  outlines[[id]] <- read.csv(outline_files[i])
  points[[id]]   <- read.csv(ppath)
}
ids <- names(outlines)
if (!length(ids)) stop("No valid outline/points pairs found.")
message(sprintf("Loaded %d sections for group '%s'.", length(ids), group_name))

## ---- Shared y reference (only matters if flip_y = TRUE) ---------------------

yref <- max(vapply(outlines, function(d) max(d[[ycol]]), numeric(1)))

# Report per-section y-extent spread (a large spread would indicate the
# sections are not well registered vertically)
fish_xy   <- lapply(outlines, outline_xy, yref = yref)
fish_wins <- lapply(fish_xy, make_win)
yext <- vapply(fish_wins, function(w) diff(w$yrange), numeric(1))
message(sprintf("y-extent across sections: min %.1f, max %.1f, sd %.2f (um)",
                min(yext), max(yext), sd(yext)))

## ---- Shared template window --------------------------------------------------

if (!is.null(template_outline_file)) {
  message("Template window: representative outline ", basename(template_outline_file))
  template_win <- make_win(outline_xy(read.csv(template_outline_file), yref))
} else {
  message(sprintf("Template window: consensus (covered by >= %.0f%% of sections)",
                  100 * template_frac))
  fr   <- do.call(boundingbox, fish_wins)
  grid <- as.mask(fr, dimyx = grid_dimyx)
  xc <- grid$xcol; yc <- grid$yrow
  masks <- lapply(fish_wins, function(w) as.mask(w, xy = list(x = xc, y = yc)))
  count <- Reduce(`+`, lapply(masks, function(m) m$m * 1L))
  thr   <- max(1, ceiling(template_frac * length(fish_wins)))
  template_win <- levelset(im(count, xcol = xc, yrow = yc), thr, ">=")
}

## ---- Effective sigma (size-proportional option) -----------------------------

template_area <- area(template_win)
if (sigma_mode == "relative") {
  sigma_eff <- sigma_frac * sqrt(template_area)
  message(sprintf("Template area %.0f um^2  (sqrt %.1f um).  relative sigma = %.3f * %.1f = %.2f um",
                  template_area, sqrt(template_area), sigma_frac, sqrt(template_area), sigma_eff))
} else {
  sigma_eff <- sigma
  message(sprintf("Template area %.0f um^2  (sqrt %.1f um).  absolute sigma = %.2f um",
                  template_area, sqrt(template_area), sigma_eff))
}

## ---- Densities on the shared window -----------------------------------------

samples <- list(); dens_list <- list()
for (id in ids) {
  fp <- flip_points(points[[id]], yref)
  inside <- inside.owin(fp$x, fp$y, template_win)
  if (any(!inside))
    message(sprintf("  %-22s %d/%d points outside template, dropped",
                    id, sum(!inside), length(fp$x)))
  pp <- ppp(fp$x, fp$y, window = template_win)
  samples[[id]]   <- pp
  dens_list[[id]] <- density(pp, sigma = sigma_eff, edge = TRUE, dimyx = grid_dimyx)
}

## ---- Averages across sections (uniform coverage) -----------------------------
avg_raw  <- im.apply(dens_list, mean, na.rm = TRUE, fun.handles.na = TRUE) * area_unit_um2
avg_raw$v[avg_raw$v < 0] <- 0
dens_norm<- lapply(dens_list, function(d) d / integral.im(d))
avg_norm <- im.apply(dens_norm, mean, na.rm = TRUE, fun.handles.na = TRUE) * area_unit_um2
coverage <- im.apply(dens_list, function(z) sum(!is.na(z)), fun.handles.na = TRUE)
dens_disp<- lapply(dens_list, function(d) d * area_unit_um2)

rib_cells <- paste0("cells / ", area_unit_lab)
n_samp    <- length(samples)

## ---- Diagnostic: registered-outline overlay ---------------------------------
# Check this figure: tight agreement between the grey section outlines means
# good registration (the template window mainly sharpens the rim). Loose
# agreement means the template window (red) is doing real work, and that
# interior crispness is limited by registration quality, not by averaging.

png(file.path(out_dir, paste0(group_name, "_outline_overlay.png")),
    width = 1200, height = 1000, res = 150)
# border = NA sets up the plot region at the full extent without drawing the box
plot(do.call(boundingbox, fish_wins), border = NA,
     main = sprintf("%s: %d registered section outlines + template",
                    group_name, n_samp))
for (w in fish_wins) plot(w, add = TRUE, border = rgb(0, 0, 0, 0.35))
# template_win is a MASK in the consensus case; plotting it directly fills black.
# Draw its polygonal boundary so only the outline shows.
plot(as.polygonal(template_win), add = TRUE, border = "red", lwd = 2)
legend("topright", c("section outline", "template window"),
       col = c(rgb(0,0,0,0.5), "red"), lwd = c(1, 2), bty = "n")
dev.off()

## ---- Average figures ----------------------------------------------------

rawvals  <- as.matrix(avg_raw); rawvals <- rawvals[is.finite(rawvals) & rawvals > 0]
truemax  <- max(rawvals)
rawclip  <- as.numeric(quantile(rawvals, raw_clip_quantile))
cm       <- gamma_cmap(rawclip, truemax)
message(sprintf("raw figure: max %.3g, clip(%.0f%%) %.3g %s",
                truemax, 100*raw_clip_quantile, rawclip, rib_cells))

tiff(file.path(out_dir, paste0(group_name, "_average_raw_linear_clip97.tif")),
     width = fig_w_in, height = fig_h_in, units = "in",
     res = fig_dpi, compression = tiff_compress)
plot(avg_raw, col = cm,
     main = sprintf("%s average (n = %d, sigma = %.1f um) [linear, clip %.0f%%]",
                    group_name, n_samp, sigma_eff, 100*raw_clip_quantile),
     riblab = rib_cells, useRaster = TRUE)
add_scalebar(avg_raw)
dev.off()

tiff(file.path(out_dir, paste0(group_name, "_average_normalised.tif")),
     width = fig_w_in, height = fig_h_in, units = "in",
     res = fig_dpi, compression = tiff_compress)
plot(avg_norm,
     main = sprintf("%s normalised (n = %d, sigma = %.1f um)", group_name, n_samp, sigma_eff),
     riblab = "probability density", useRaster = TRUE)
add_scalebar(avg_norm)
dev.off()

## ---- Individual per-section intensity maps -----------------------------------
# One PNG per section (cells overlaid, scale bar) plus a combined panel, all on
# the shared template window - useful for spotting any outlier section.

ind_dir <- file.path(out_dir, paste0(group_name, "_individuals"))
dir.create(ind_dir, showWarnings = FALSE, recursive = TRUE)

for (id in names(dens_disp)) {
  png(file.path(ind_dir, paste0(id, "_intensity.png")),
      width = 1200, height = 1000, res = 150)
  plot(dens_disp[[id]], main = paste0(id, "  (sigma = ", round(sigma_eff,1), " um)"),
       riblab = rib_cells)
  plot(samples[[id]], add = TRUE, pch = 16, cex = 0.3, cols = rgb(1, 1, 1, 0.5))
  add_scalebar(dens_disp[[id]])
  dev.off()
}

png(file.path(out_dir, paste0(group_name, "_all_individuals.png")),
    width = 1600, height = 1200, res = 150)
plot(as.solist(dens_disp), main = paste0(group_name, ": individual intensity maps"))
dev.off()
message("wrote ", n_samp, " individual maps to ", ind_dir)

## ---- Save results (consumed by downstream comparison/figure scripts) --------
# This .rds is the input expected by the companion multi-group script
# (e.g. "per_group_tiffs_shared_ribbon.R"): it needs at least average_raw and
# area_unit_um2, with the SAME area_unit_um2 (and ideally the same sigma)
# across every group you plan to compare on one shared colour scale.

saveRDS(list(patterns = samples, densities = dens_disp,
             average_raw = avg_raw, average_normalised = avg_norm,
             coverage = coverage, area_unit_um2 = area_unit_um2,
             template_win = template_win,
             sigma_used = sigma_eff, sigma_mode = sigma_mode,
             sigma_frac = if (sigma_mode == "relative") sigma_frac else NA),
        file.path(out_dir, paste0(group_name, "_intensity_results.rds")))

message("Done. Outputs in: ", normalizePath(out_dir))

###############################################################################
# HOW TO USE
#
# 1. Prepare, for each group you want to average, a folder of matched CSV
#    pairs in a common (registered/warped) coordinate frame:
#      <outline_prefix><section_id>.csv   - polygon outline of the section
#      <points_prefix><section_id>.csv    - detected cell centroids
#    both with x/y columns named as set in xcol/ycol (microns).
#
# 2. Edit the CONFIG block above for the FIRST group: outline_dir, points_dir,
#    group_name, and out_dir at minimum. Leave sigma / area_unit_um2 / grid
#    settings the same across all groups you intend to compare later.
#
# 3. Run the whole script. Check the console messages and the
#    "<group_name>_outline_overlay.png" diagnostic before trusting the result.
#
# 4. Edit CONFIG again for the NEXT group (only outline_dir / points_dir /
#    group_name / out_dir normally need to change) and re-run. Repeat once
#    per group.
#
# 5. Once you have one "<group_name>_intensity_results.rds" per group, feed
#    them into a downstream comparison/figure script (e.g.
#    "per_group_tiffs_shared_ribbon.R") to render them on one shared colour
#    scale.
#
# WHAT TO CHECK
#
# * The console "y-extent ... sd" line: if sd is more than ~1 pixel of real
#   size, flip_y / the shared y reference is correcting a genuine vertical
#   jitter between sections.
#
# * The "_outline_overlay.png": tight agreement of the grey outlines means
#   good registration. Loose agreement means the template window (red) is
#   doing real work, and flags that interior crispness is limited by
#   registration, not by averaging.
#
# * Points dropped "outside template": a few is fine (sections that bulge
#   past the consensus). Many means template_frac is too strict, or a
#   section is poorly registered - check it in the overlay.
###############################################################################
