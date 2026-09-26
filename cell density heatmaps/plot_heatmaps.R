###############################################################################
# Per-group averaged intensity maps -> individual TIFFs + one shared ribbon TIFF
#
# Reads the PRECOMPUTED average_raw (an already-averaged, already-smoothed
# cell-density map) straight out of each group's "*_intensity_results.rds"
# file. No density estimation is re-run here - this script only visualises
# results that another script already computed (e.g. a companion script
# such as "Averaged_heatmaps_from_warps.R"). Flexible for ANY number of
# groups from 4 to 20 (just edit the `groups` list below).
#
# For each group, this script writes its OWN standalone TIFF (bare image +
# scale bar, black background). It also writes ONE separate "ribbon" TIFF
# (a colour-scale legend) that applies to ALL groups, built from a colour
# map shared across every group's data - so the same colour always means the
# same density, no matter which group's image you're looking at.
#
# Requirements
#   - Every input .rds must contain: average_raw (a spatstat "im" object)
#     and area_unit_um2 (a number). area_unit_um2 must be IDENTICAL across
#     all groups you list below (so "cells per area_unit_um2" means the same
#     thing everywhere) - the script stops with an error otherwise.
#   - Ideally, every group's average_raw should also have been computed with
#     the same smoothing bandwidth (sigma) upstream, so that the same colour
#     truly means the same underlying density pattern across groups. This
#     script does not check sigma itself (it isn't stored under a fixed
#     name), so keep track of that when preparing your inputs.
#
# Required R packages: spatstat
###############################################################################

library(spatstat)

## ============================================================================
## ---- CONFIG: GROUPS (4 to 20) -----------------------------------------------
## ============================================================================
# Each group needs:
#   file = path to that group's "*_intensity_results.rds"
#           (produced upstream, e.g. by "Averaged_heatmaps_from_warps.R")
#   name = short label used to build the output filename for that group
#
# `base_dir` below is just a convenience so the list stays readable - each
# `file` path is built from it with file.path(). You can also give fully
# independent paths per group if your files live in different locations.
#
# EDIT the paths and names below to match your own project. The examples
# show two anatomical regions ("RegionA"/"RegionB") x five conditions
# ("Group1".."Group5"), but any grouping scheme works - just list every
# .rds file you want to appear on the shared colour scale.

base_dir <- "PATH/TO/YOUR/PROJECT/AverageCells/Sigma_4"   # folder holding each group's *_intensity_results.rds
out_dir  <- "PATH/TO/YOUR/PROJECT/AverageCells/PerGroupTiffs"  # where output TIFFs are written

groups <- list(
  # --- Region A ---
  list(file = file.path(base_dir, "Group1", "RegionA", "Group1_RegionA_intensity_results.rds"), name = "Group1_RegionA"),
  list(file = file.path(base_dir, "Group2", "RegionA", "Group2_RegionA_intensity_results.rds"), name = "Group2_RegionA"),
  list(file = file.path(base_dir, "Group3", "RegionA", "Group3_RegionA_intensity_results.rds"), name = "Group3_RegionA"),
  list(file = file.path(base_dir, "Group4", "RegionA", "Group4_RegionA_intensity_results.rds"), name = "Group4_RegionA"),
  list(file = file.path(base_dir, "Group5", "RegionA", "Group5_RegionA_intensity_results.rds"), name = "Group5_RegionA"),
  # --- Region B ---
  list(file = file.path(base_dir, "Group1", "RegionB", "Group1_RegionB_intensity_results.rds"), name = "Group1_RegionB"),
  list(file = file.path(base_dir, "Group2", "RegionB", "Group2_RegionB_intensity_results.rds"), name = "Group2_RegionB"),
  list(file = file.path(base_dir, "Group3", "RegionB", "Group3_RegionB_intensity_results.rds"), name = "Group3_RegionB"),
  list(file = file.path(base_dir, "Group4", "RegionB", "Group4_RegionB_intensity_results.rds"), name = "Group4_RegionB"),
  list(file = file.path(base_dir, "Group5", "RegionB", "Group5_RegionB_intensity_results.rds"), name = "Group5_RegionB")
)

## ---- CONFIG: display / clip / output ----------------------------------------

# Shared colourmap: white point is set at this percentile of the POOLED data
# across all "driving" groups (see clip_exclude below). Raise this to make
# only the very brightest hotspots saturate; lower it to boost contrast in
# dimmer regions at the cost of saturating more of the bright areas.
clip_quantile <- 0.97
clip_exclude  <- character(0)   # group `name`s to exclude when computing the
                                 # shared clip value (their data still gets
                                 # coloured with the resulting scale, they
                                 # just don't influence where the scale's
                                 # white point is set) - leave empty to use
                                 # every group

gamma        <- 1.0             # colour-scale gamma; 1.0 = linear
palette_name <- "Plasma"         # any palette name accepted by hcl.colors()
n_cols_lut   <- 255              # number of discrete colour steps in the map

draw_title <- FALSE             # FALSE = bare image (best for assembling into
                                 #         a multi-panel figure later)
                                 # TRUE  = print the group name as a title on
                                 #         each image TIFF

# Per-image TIFF size (inches) and resolution
img_dpi <- 600; img_w_in <- 4; img_h_in <- 4; tiff_compress <- "lzw"
# Ribbon (colour-scale legend) TIFF size (inches) and resolution
rib_dpi <- 600; rib_w_in <- 1.6; rib_h_in <- 4

# Scale bar (drawn in the bottom outer margin so it never overlaps tissue)
scalebar_len <- 100; scalebar_label <- "100 µm"   # length in the same
                                                        # units as your image
                                                        # coordinates (e.g. um)
scalebar_lwd <- 3; scalebar_cex <- 0.9
## ============================================================================
## ---- END CONFIG --------------------------------------------------------------
## ============================================================================

stopifnot(length(groups) >= 4, length(groups) <= 20)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- Helpers ----------------------------------------------------------------

# Scale bar in the bottom OUTER margin band (below the section), white on black.
add_scalebar_outer <- function(im, len = scalebar_len, label = scalebar_label,
                               col = "white", halo = "black",
                               lwd = scalebar_lwd, cex = scalebar_cex) {
  op <- par(xpd = NA)
  usr <- par("usr")
  x1 <- usr[2] - 0.05*diff(usr[1:2]); x0 <- x1 - len
  yb <- usr[3] - 0.06*diff(usr[3:4])
  segments(x0,yb,x1,yb,col=halo,lwd=lwd+2,lend="butt")
  segments(x0,yb,x1,yb,col=col, lwd=lwd,  lend="butt")
  xm <- (x0+x1)/2; dy <- 0.02*diff(usr[3:4])
  text(xm, yb-dy, label, col=col, cex=cex, adj=c(0.5,1))
  par(op); invisible(NULL)
}

gamma_cmap <- function(clip, truemax, g = gamma, n = n_cols_lut, pal = palette_name) {
  cols <- hcl.colors(n, pal)
  br <- clip * (seq(0, 1, length.out = n + 1) ^ (1 / g))
  br[n + 1] <- max(truemax, clip) * (1 + 1e-6)
  colourmap(cols, breaks = br)
}

pos_vals <- function(im){ v <- as.matrix(im); v[is.finite(v) & v > 0] }

## ---- Load precomputed averages ----------------------------------------------

avgs <- lapply(groups, function(g) {
  if (!file.exists(g$file)) stop("File not found: ", g$file)
  d <- readRDS(g$file)
  if (!all(c("average_raw","area_unit_um2") %in% names(d)))
    stop("'", g$file, "' missing average_raw / area_unit_um2.")
  list(im = d$average_raw, area = d$area_unit_um2)
})
names(avgs) <- vapply(groups, `[[`, "", "name")

if (length(unique(vapply(avgs, `[[`, numeric(1), "area"))) != 1)
  stop("area_unit_um2 differs between files; cannot share one ribbon.")
raw_lab <- sprintf("cells / %g um^2", avgs[[1]]$area)

ims <- lapply(avgs, `[[`, "im")

## ---- Shared colourmap across ALL groups -------------------------------------

driving_names <- setdiff(names(ims), clip_exclude)
pooled  <- unlist(lapply(ims[driving_names], pos_vals))
truemax <- max(unlist(lapply(ims, pos_vals)))
clip    <- as.numeric(quantile(pooled, clip_quantile))
cm      <- gamma_cmap(clip, truemax)
message(sprintf("[shared] clip(%.0f%%) %.3g, truemax %.3g %s",
                100*clip_quantile, clip, truemax, raw_lab))

## ---- Write one TIFF per group (bare image + scale bar, black bg) -------------

for (i in seq_along(groups)) {
  g <- groups[[i]]; im <- ims[[i]]

  # Out-of-tissue pixels -> NA so background renders BLACK, not the colourmap floor
  im_disp <- im
  im_disp$v[!is.finite(im_disp$v) | im_disp$v <= 0] <- NA

  fn <- file.path(out_dir, sprintf("avg_%s.tif", g$name))
  tiff(fn, width = img_w_in, height = img_h_in, units = "in",
       res = img_dpi, compression = tiff_compress, bg = "black")
  op <- par(bg = "black", fg = "white", col.main = "white")
  par(mar = if (draw_title) c(2.5,0.2,2.5,0.2) else c(2.5,0.2,0.2,0.2))
  plot(im_disp, col = cm, ribbon = FALSE, useRaster = TRUE,
       main = if (draw_title) g$name else "",
       box = FALSE, na.col = "black")
  add_scalebar_outer(im_disp)
  par(op); dev.off()
  message("wrote ", fn)
}

## ---- Write the shared ribbon as its own TIFF (black bg) ----------------------

rib_fn <- file.path(out_dir, "SHARED_ribbon.tif")
tiff(rib_fn, width = rib_w_in, height = rib_h_in, units = "in",
     res = rib_dpi, compression = tiff_compress, bg = "black")
op <- par(bg = "black", fg = "white", col.axis = "white",
          col.lab = "white", col.main = "white")
par(mar = c(3, 0.5, 3, 4))
plot(cm, vertical = TRUE,
     main = sprintf("%s\n[clip %.0f%%]", raw_lab, 100*clip_quantile))
par(op); dev.off()
message("wrote ", rib_fn)

message("Done. Outputs in: ", normalizePath(out_dir))

###############################################################################
# NOTES
# * This script reads each group's already-averaged average_raw directly; no
#   density estimation is re-run here. If you change how the upstream average
#   is computed (e.g. a different sigma), re-generate that group's .rds and
#   re-run this script.
# * A SINGLE common sigma/pipeline upstream is what makes the shared ribbon a
#   true "same colour = same density" scale across groups. Unit: cells per
#   <area_unit_um2> um^2 (taken from the .rds files themselves).
# * Background is BLACK: out-of-tissue pixels (<=0 / non-finite) are set to
#   NA and rendered black. If a faint edge-correction halo just outside the
#   tissue still shows a coloured rim, consider masking by the upstream
#   template window instead of by value.
# * The scale bar sits in the bottom outer margin (white on black) so it
#   never obscures the tissue; reduce scalebar_len for small panels.
# * na.col is a current-spatstat plot.im() argument; on older spatstat
#   versions you may need a different approach for colouring NA pixels.
###############################################################################
