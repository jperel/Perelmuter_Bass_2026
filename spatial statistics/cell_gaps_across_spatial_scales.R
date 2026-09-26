###############################################################################
# Size-normalised interior void-size curves (2-6 groups x 1+ levels) + void-AUC
# statistical comparison - UNWARPED, per-animal
# -----------------------------------------------------------------------------
# WHAT THIS SCRIPT DOES
#   For each animal, "cell gap" is measured as the distance from every
#   interior point of a traced section to its nearest detected cell - i.e. how
#   large a cell-free gap is centred there. This script:
#
#   1. Loads each animal's traced section outline + detected cell positions.
#   2. For each animal, builds a gap-distance map over the section interior
#      (excluding a configurable rim) and rescales gap radius per animal in
#      three ways so animals of different size/density are comparable:
#        - by sqrt(section area)         (PRIMARY - pure geometric size)
#        - by each animal's own median nearest-neighbour cell distance
#          (cross-check: are groups differently SPACED, not just differently
#          sized?)
#        - left in absolute microns      (for reference)
#   3. Plots, per group x level, V(r) = the fraction of the section that is
#      farther than r from any cell, averaged across animals with a 95% CI
#      band -> depth_profile-style curves:
#        interior_zone_summary.tif   one example animal's interior zone, per
#                                     group x level (sanity-check figure)
#        void_curves_area.tif        V(r) curves, PRIMARY (area-normalised)
#        void_curves_nnd.tif         V(r) curves, spacing-normalised
#        void_curves_micron.tif      V(r) curves, absolute microns
#        gap_q95_per_fish.csv        each animal's 95th-percentile gap radius,
#                                     in all three normalisations
#        voidsize_plot_data.rds      every animal's raw data, for reuse/redraw
#   4. Reduces each animal's curve (on ONE chosen normalisation, `auc_axis_mode`
#      below) to a single number - AUC = the area under its V(r) curve, i.e.
#      that animal's MEAN gap radius - and compares groups statistically:
#        void_auc_results_<mode>.csv    Mann-Whitney U + Cliff's delta, each
#                                        group vs a reference group, within
#                                        level, Holm-corrected; for BOTH the
#                                        AUC metric and a q95 cross-check
#        void_auc_<mode>_<level>.tif    one box-and-jitter figure PER LEVEL
#                                        (e.g. one file for "rostral", one for
#                                        "caudal"), one box per group
#
# GETTING STARTED (read this before editing CONFIG below)
#
#   Requirements: R (4.x+) and the "spatstat" package. Install once with:
#       install.packages("spatstat")
#
#   Your input data: one CSV pair per animal, per group, per level. Every
#   animal needs BOTH of these, in the same folder, sharing the same <id>:
#     <id>_outline_xy.csv   - the traced section outline, as an ORDERED list
#                              of (x, y) points walking around the boundary
#     <id>_centroids.csv    - one row per detected cell, with its (x, y)
#                              position in the SAME coordinate system
#   Both need x/y columns. Accepted column-name variants are listed in
#   .xcand / .ycand further down (x_microns, x_um, x, X - and the y
#   equivalents); rename your columns to one of those, or add your own name
#   to the list. Coordinates are assumed to be in a consistent physical unit
#   (microns, per the axis labels) - relabel those if you use something else.
#
#   Example outline_xy.csv (an open polygon; you'll have many more rows):
#       x_um,y_um
#       12.4,318.9
#       15.1,322.0
#       19.8,325.6
#       ...
#
#   Example centroids.csv:
#       x_um,y_um
#       140.2,210.5
#       138.9,205.1
#       ...
#
#   Expected folder layout (built automatically from CONFIG below as
#   DATA_DIR/<folder code>/<level>/output/):
#
#       DATA_DIR/
#       |-- Reference_code/
#       |   |-- Rostral/
#       |   |   `-- output/
#       |   |       |-- animal01_outline_xy.csv
#       |   |       |-- animal01_centroids.csv
#       |   |       `-- ...
#       |   `-- Caudal/output/   (same pattern)
#       |-- GroupB_code/
#       |   |-- Rostral/output/  (same pattern)
#       |   `-- Caudal/output/   (same pattern)
#       `-- GroupC_code/
#           |-- Rostral/output/  (same pattern)
#           `-- Caudal/output/   (same pattern)
#
#   If your data isn't organised this way, either reorganise it to match, or
#   edit the `groups` list built in CONFIG below to point `folder` directly
#   at wherever each group/level's CSV pairs live.
#
# Required R packages: spatstat (everything else is base R)
###############################################################################

library(spatstat)

## ============================================================================
## ---- CONFIG: EDIT FOR YOUR DATA, THEN RUN THE WHOLE SCRIPT -----------------
## ============================================================================

# ---- Input data ----------------------------------------------------------------
DATA_DIR <- "PATH/TO/YOUR/PROJECT/MasterData"   # <-- edit: root of the data tree (see the folder layout above)

# Folder code -> group label (2 to 6 groups). The label is just what appears
# in your output files/figures/legends - species, genotype, treatment,
# anything. The FIRST entry is treated as the REFERENCE group for the void-AUC
# statistical comparisons in part 2 ("reference vs every other group").
GCODES <- c(Reference_code = "Reference group",
            GroupB_code    = "Comparison group B",
            GroupC_code    = "Comparison group C")

# Subfolders under each group's data (e.g. anatomical subregions imaged
# separately per animal). One entry is fine if you only have a single level.
LEVELS <- c("Rostral", "Caudal")

outline_suffix <- "_outline_xy.csv"
points_suffix  <- "_centroids.csv"

groups <- list()
for (code in names(GCODES)) for (Lv in LEVELS)
  groups[[length(groups)+1]] <- list(
    folder = file.path(DATA_DIR, code, Lv, "output"),
    group  = unname(GCODES[code]),
    level  = tolower(Lv))

# ---- Output location -------------------------------------------------------------
out_dir <- "PATH/TO/YOUR/PROJECT/voidsize_normalized"   # <-- edit: all outputs (part 1 and part 2) go here

# ---- Display order + colours ------------------------------------------------------
# Explicit left-to-right / panel order. Any group that actually appears in
# `groups` but is missing here is appended at the end automatically.
group_order <- c("Reference group", "Comparison group B", "Comparison group C")
level_order <- c("rostral", "caudal")

group_col <- c("Reference group"     = "#50C878",
              "Comparison group B" = "#FF00FF",
              "Comparison group C" = "#E61C24")

# ---- Part 1: void-curve settings ---------------------------------------------------
interior_depth <- 0     # exclude a rim from the interior: 0 = use the whole
                        # section; raise (in the depth map's 0-1 units) to
                        # trim a border/midline strip before measuring gaps
grid_dimyx     <- 512    # pixel grid resolution for the interior depth map
n_rgrid        <- 60     # number of r-values sampled along each V(r) curve
ci_level       <- 0.95   # confidence level for the curve bands
ref_line       <- c(area = NA, nnd = 2, micron = 15)   # optional dashed
                        # reference line per mode (NA = don't draw one)
fig_dpi <- 600; tiff_compress <- "lzw"
fade_col <- adjustcolor("white", 0.62); boundary_col <- "black"
scalebar_len <- 100; scalebar_label <- "100 µm"

# ---- Part 2: void-AUC statistics + figure settings ----------------------------------
auc_axis_mode   <- "area"   # which normalisation to compute AUC on: "area"
                            # (primary - void-richness), "nnd" (spatial
                            # organisation cross-check), or "micron"
auc_p_adjust    <- "holm"   # p.adjust() method, applied within each level
                            # across the (n groups - 1) comparisons
auc_fig_dpi <- 300
auc_fig_w   <- 1500; auc_fig_h <- 1400   # pixels, one TIFF written per level
auc_show_sig  <- FALSE   # TRUE: mark each non-reference group's box with
                         # * / ** / *** (Holm-adjusted, vs the reference
                         # group, that level) using the AUC metric's result
auc_sig_alpha <- 0.05
auc_group_order <- NULL   # NULL = auto (reference group drawn last, others in
                          # group_order); or supply your own character vector
                          # of group labels to set the box order yourself
## ============================================================================
## ---- END CONFIG ------------------------------------------------------------------
## ============================================================================

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

.in_groups_gr <- unique(vapply(groups, function(g) g$group, ""))
.in_groups_lv <- unique(vapply(groups, function(g) g$level, ""))
group_present <- c(intersect(group_order, .in_groups_gr), setdiff(.in_groups_gr, group_order))
levels_present<- c(intersect(level_order, .in_groups_lv), setdiff(.in_groups_lv, level_order))
reference_group <- group_present[1]                       # first entry in GCODES
other_groups    <- setdiff(group_present, reference_group)
if (length(group_present) < 2 || length(group_present) > 6)
  warning(sprintf("This script is written for 2-6 groups; found %d.", length(group_present)))

# Fail loudly if any group lacks a colour (a silent NA colour is a confusing
# way to find out later that a group name doesn't match).
.missing_col <- setdiff(group_present, names(group_col))
if (length(.missing_col))
  stop("group_col is missing colours for: ", paste(.missing_col, collapse=", "))
group_col <- group_col[group_present]

## ---- Loader: reads one animal's outline + centroid CSVs into a spatstat point pattern

.xcand <- c("x_microns", "x_um", "x", "X")
.ycand <- c("y_microns", "y_um", "y", "Y")
.pick <- function(df, cands, what, file) {
  hit <- cands[cands %in% names(df)]
  if (!length(hit)) stop(sprintf("No %s column in %s (have: %s) - add your column name to .xcand/.ycand above",
                                 what, basename(file), paste(names(df), collapse=", ")))
  df[[hit[1]]]
}
.owin_from_outline <- function(df, file) {
  x <- .pick(df, .xcand, "x", file); y <- .pick(df, .ycand, "y", file)
  n <- length(x)
  if (n>1 && isTRUE(all.equal(x[1],x[n])) && isTRUE(all.equal(y[1],y[n]))) { x<-x[-n]; y<-y[-n] }  # drop repeated closing vertex, if any
  nx <- c(x[-1],x[1]); ny <- c(y[-1],y[1])
  if (0.5*sum(x*ny - nx*y) < 0) { x<-rev(x); y<-rev(y) }   # force consistent (anticlockwise) winding
  owin(poly = list(x = x, y = y))
}
load_group_patterns <- function(folder, outline_suffix, points_suffix, verbose=TRUE) {
  all <- list.files(folder, full.names = TRUE)
  ofiles <- all[endsWith(basename(all), outline_suffix)]
  if (!length(ofiles)) stop("No files matching '*", outline_suffix, "' found in: ", folder,
                            "\n  Check DATA_DIR / GCODES / LEVELS in CONFIG match your actual folders.")
  pats <- list()
  for (of in ofiles) {
    b <- basename(of); id <- substr(b, 1, nchar(b)-nchar(outline_suffix))
    pf <- file.path(folder, paste0(id, points_suffix))
    if (!file.exists(pf)) { warning("No matching centroids file for ", id, " in ", folder, " - skipping"); next }
    W <- .owin_from_outline(read.csv(of), of); cdf <- read.csv(pf)
    px <- .pick(cdf, .xcand, "x", pf); py <- .pick(cdf, .ycand, "y", pf)
    ins <- inside.owin(px, py, W)
    if (any(!ins) && verbose)
      message(sprintf("  %-30s %d/%d centroids outside outline, dropped", id, sum(!ins), length(px)))
    pats[[id]] <- ppp(px[ins], py[ins], window = W)
  }
  if (!length(pats)) stop("No valid outline/centroid pairs found in: ", folder)
  pats
}

sqrt_cmap <- function(zmax, pal="Plasma", n=255)
  colourmap(hcl.colors(n, pal), breaks = seq(0, sqrt(zmax*(1+1e-6)), length.out=n+1)^2)
add_scalebar <- function(obj, len=scalebar_len, label=scalebar_label) {
  R <- as.rectangle(obj); w <- diff(R$xrange); h <- diff(R$yrange)
  x1 <- R$xrange[2]-0.07*w; x0 <- x1-len; yb <- R$yrange[1]+0.07*h
  segments(x0,yb,x1,yb,col="white",lwd=5,lend="butt")
  segments(x0,yb,x1,yb,col="black",lwd=3,lend="butt")
  text((x0+x1)/2, yb+0.03*h, label, col="black", cex=0.9, adj=c(0.5,0))
}

## ---- Collect per animal: interior gap distances, spacing, area ----------------------
# Each animal's interior is defined by ITS OWN outline (per-animal depth map).

collect_group <- function(g) {
  message(sprintf("[%s | %s] loading: %s", g$group, g$level, g$folder))
  pats <- load_group_patterns(g$folder, outline_suffix, points_suffix)
  pats <- pats[vapply(pats, npoints, 0L) >= 2]
  dints <- list(); nnds <- numeric(); areas <- numeric()
  rep_pp <- NULL; rep_depth <- NULL
  for (k in seq_along(pats)) {
    pp <- pats[[k]]; W <- Window(pp)
    depth <- bdist.pixels(as.mask(W, dimyx = grid_dimyx))     # 0 at boundary, rising inward
    depth <- depth / max(as.matrix(depth), na.rm = TRUE)      # normalised to [0, 1]
    dm  <- distmap(pp, dimyx = grid_dimyx)                    # distance to nearest cell, per pixel
    hh  <- harmonise(depth = depth, dm = dm)
    dpv <- as.matrix(hh$depth); dmv <- as.matrix(hh$dm)
    sel <- is.finite(dpv) & is.finite(dmv) & (dpv > interior_depth)
    dints[[k]] <- dmv[sel]                                    # this animal's interior gap distances
    nnds[k]    <- median(nndist(pp))                          # typical cell spacing
    areas[k]   <- area(W)
    if (k == 1) { rep_pp <- pp; rep_depth <- depth }           # keep one example animal for the sanity figure
  }
  list(group=g$group, level=g$level, n=length(pats),
       dints=dints, nnds=nnds, areas=areas, rep_pp=rep_pp, rep_depth=rep_depth)
}

prof <- lapply(groups, collect_group)
get_group <- function(gr, lv) {
  i <- which(vapply(prof, function(p) p$group==gr && p$level==lv, TRUE))
  if (length(i)) prof[[i[1]]] else NULL
}
normed <- function(p, mode) lapply(seq_along(p$dints), function(j) {
  d <- p$dints[[j]]
  switch(mode, area = d / sqrt(p$areas[j]), nnd = d / p$nnds[j], micron = d)
})

## ---- Interior-zone summary (representative animal per group) ------------------------

tiff(file.path(out_dir, "interior_zone_summary.tif"),
     width = 4.6*length(levels_present), height = 4.2*length(group_present),
     units = "in", res = fig_dpi, compression = tiff_compress)
op <- par(mfrow = c(length(group_present), length(levels_present)), mar=c(1,1,2.5,1))
for (gr in group_present) for (lv in levels_present) {
  gg <- get_group(gr, lv)
  if (is.null(gg) || is.null(gg$rep_pp)) { plot.new(); title(paste(gr, lv)); next }
  pp <- gg$rep_pp
  plot(Window(pp), main = sprintf("%s %s (example)", gr, lv), border="grey50")
  plot(pp, add=TRUE, pch=16, cex=0.25, cols=adjustcolor("black",0.6))
  plot(levelset(gg$rep_depth, interior_depth, "<="), add=TRUE, col=fade_col, border=NA)
  contour(gg$rep_depth, levels=interior_depth, add=TRUE, drawlabels=FALSE,
          col=boundary_col, lwd=2)
  add_scalebar(pp)
}
par(op); dev.off()
message("wrote ", file.path(out_dir, "interior_zone_summary.tif"))

## ---- Void curves, one figure per normalisation mode ----------------------------------

band <- function(rgrid, s, col) {
  ok <- is.finite(s$lo) & is.finite(s$hi)
  if (any(ok)) polygon(c(rgrid[ok],rev(rgrid[ok])), c(s$lo[ok],rev(s$hi[ok])),
                       col=adjustcolor(col,0.18), border=NA)
  lines(rgrid, s$mean, col=col, lwd=2.5)
}
curve_summary <- function(vecs, rgrid) {
  mat <- sapply(vecs, function(d) vapply(rgrid, function(r) mean(d > r), 0))
  if (is.null(dim(mat))) mat <- matrix(mat, ncol=1)
  m <- rowMeans(mat, na.rm=TRUE); nn <- rowSums(is.finite(mat))
  se <- apply(mat, 1, function(x){x<-x[is.finite(x)]; if(length(x)<2) NA else sd(x)/sqrt(length(x))})
  tc <- ifelse(nn>1, qt(1-(1-ci_level)/2, df=nn-1), NA)
  list(mean=m, lo=m-tc*se, hi=m+tc*se)
}
mode_xlab <- c(area  = "gap radius, section-normalized",
               nnd   = "gap radius, in units of typical cell spacing  (r / median NND)",
               micron= "gap radius, r  (µm from nearest cell)")

# Build + store the curve data for each mode, then plot from the stored data so
# the exact same summaries can be re-loaded and re-plotted/edited later.
mode_data <- list()
render_mode <- function(mode) {
  pooled <- unlist(lapply(prof, function(p) unlist(normed(p, mode))))
  rmax <- as.numeric(quantile(pooled, 0.995)); rgrid <- seq(0, rmax, length.out=n_rgrid)
  gsum <- list()
  for (p in prof)
    gsum[[paste(p$group, p$level)]] <-
      list(group=p$group, level=p$level, summary=curve_summary(normed(p, mode), rgrid))
  mode_data[[mode]] <<- list(rgrid=rgrid, rmax=rmax, xlab=mode_xlab[[mode]], groups=gsum)

  tiff(file.path(out_dir, sprintf("void_curves_%s.tif", mode)),
       width=5*length(levels_present), height=4.5, units="in",
       res=fig_dpi, compression=tiff_compress)
  op <- par(mfrow=c(1,length(levels_present)), mar=c(4.5,4.5,3,1)); first <- TRUE
  for (lv in levels_present) {
    plot(NA, xlim=c(0,rmax), ylim=c(0,1), main=lv,
         xlab=mode_xlab[[mode]],
         ylab="cell free area fraction")
    if (is.finite(ref_line[[mode]])) abline(v=ref_line[[mode]], col="grey70", lty=3)
    for (gr in group_present) { gg <- get_group(gr, lv)
      if (!is.null(gg)) band(rgrid, gsum[[paste(gr, lv)]]$summary, group_col[[gr]]) }
    if (first) { legend("topright", names(group_col), col=unlist(group_col), lwd=2.5, bty="n"); first <- FALSE }
  }
  par(op); dev.off()
  message("wrote ", file.path(out_dir, sprintf("void_curves_%s.tif", mode)))
}
invisible(lapply(c("area","nnd","micron"), render_mode))

## ---- Per-fish q95 (descriptive) -------------------------------------------------------

q95_tab <- do.call(rbind, lapply(prof, function(p) data.frame(
  group=p$group, level=p$level,
  q95_area   = sapply(normed(p,"area"),   quantile, 0.95),
  q95_nnd    = sapply(normed(p,"nnd"),    quantile, 0.95),
  q95_micron = sapply(normed(p,"micron"), quantile, 0.95))))
write.csv(q95_tab, file.path(out_dir, "gap_q95_per_fish.csv"), row.names = FALSE)

## ---- Save editable plot data -----------------------------------------------------------
# Reload with readRDS() and redraw/edit any curve without recomputing. Contents:
#   prof        per-animal raw data (dints, nnds, areas, representative animal)
#   mode_data   per normalisation: rgrid + per-group mean/lo/hi curve summaries
#   styling     group_col, ref_line, interior_depth, levels/group order
saveRDS(list(prof = prof, mode_data = mode_data, q95 = q95_tab,
             group_col = group_col, ref_line = ref_line,
             interior_depth = interior_depth,
             levels_present = levels_present, group_present = group_present,
             band = band, curve_summary = curve_summary, normed_doc =
               "normed(p, mode): per-animal gap vectors; mode in area/nnd/micron"),
        file.path(out_dir, "voidsize_plot_data.rds"))
message("wrote ", file.path(out_dir, "voidsize_plot_data.rds"))

###############################################################################
# ---- PART 2: void-AUC per animal + statistical comparison + per-level figures
#
# AUC = area under an animal's V(r) curve on the `auc_axis_mode` normalisation
# = that animal's MEAN gap radius on that scale. Uses the SAME r-grid as the
# void_curves_<mode>.tif figure above (mode_data[[auc_axis_mode]]$rgrid), so
# the AUC and the curve figure are always talking about the same axis.
###############################################################################

rgrid_auc <- mode_data[[auc_axis_mode]]$rgrid

animal_auc <- function(dvec, rgrid) {
  V <- vapply(rgrid, function(r) mean(dvec > r), 0)   # fraction of section with gap > r
  dr <- diff(rgrid)
  sum((head(V, -1) + tail(V, -1)) / 2 * dr)            # trapezoidal integral
}

## ---- Per-animal table: group, level, AUC, q95 ----------------------------------------
rows <- list()
for (p in prof) {
  vecs <- normed(p, auc_axis_mode)
  for (j in seq_along(vecs)) {
    d <- vecs[[j]]
    rows[[length(rows) + 1L]] <- data.frame(
      group = p$group, level = p$level,
      auc = animal_auc(d, rgrid_auc),
      q95 = as.numeric(quantile(d, 0.95)),
      stringsAsFactors = FALSE)
  }
}
per_fish_auc <- do.call(rbind, rows)

## ---- Stats: Mann-Whitney U + Cliff's delta, reference vs every other group -----------
# Cliff's delta = P(X>Y) - P(X<Y) for X=reference, Y=other; > 0 means the
# reference group tends to have the LARGER value.
cliffs_delta <- function(x, y) {
  gt <- sum(outer(x, y, ">")); lt <- sum(outer(x, y, "<"))
  (gt - lt) / (length(x) * length(y))
}
delta_mag <- function(d) {                    # Romano et al. (2006) thresholds
  a <- abs(d)
  if (a < 0.147) "negligible" else if (a < 0.33) "small" else
    if (a < 0.474) "medium" else "large"
}

run_one <- function(metric) {
  res <- data.frame()
  for (lv in levels_present) {
    sub <- per_fish_auc[per_fish_auc$level == lv, ]
    xf  <- sub[[metric]][sub$group == reference_group]
    if (length(xf) < 2) { message(sprintf("[%s] reference group has <2 animals, skipped", lv)); next }
    block <- data.frame()
    for (gr in other_groups) {
      yo <- sub[[metric]][sub$group == gr]
      if (length(yo) < 2) {
        message(sprintf("[%s] %s vs %s: <2 animals, skipped", lv, reference_group, gr)); next
      }
      w  <- suppressWarnings(wilcox.test(xf, yo, exact = NULL))
      dl <- cliffs_delta(xf, yo)
      block <- rbind(block, data.frame(
        level = lv, metric = metric,
        comparison = sprintf("%s vs %s", reference_group, gr),
        n_ref = length(xf), n_other = length(yo),
        median_ref = median(xf), median_other = median(yo),
        W = unname(w$statistic), p_value = w$p.value,
        cliffs_delta = dl, magnitude = delta_mag(dl)))
    }
    if (nrow(block))
      block$p_holm <- p.adjust(block$p_value, method = auc_p_adjust)   # within level
    res <- rbind(res, block)
  }
  res
}

res_auc <- run_one("auc")
res_q95 <- run_one("q95")
auc_results <- rbind(res_auc, res_q95)

auc_results_csv <- file.path(out_dir, sprintf("void_auc_results_%s.csv", auc_axis_mode))
write.csv(auc_results, auc_results_csv, row.names = FALSE)
message("wrote ", auc_results_csv)

cat(sprintf("\n=== Void AUC (%s axis): %s vs each group ===\n", auc_axis_mode, reference_group))
print(auc_results[, c("level","metric","comparison","p_value","p_holm",
                      "cliffs_delta","magnitude")], row.names = FALSE)

## ---- One box-and-jitter figure PER LEVEL (e.g. one "rostral", one "caudal") ----------
# One box per group, full colour, in auc_group_order (reference group drawn
# last by default). Optionally marks each non-reference group's box with a
# significance asterisk vs the reference group (auc_show_sig above).

auc_ylab <- sprintf("void AUC  (mean gap radius, %s-normalised)", auc_axis_mode)
if (is.null(auc_group_order)) auc_group_order <- c(other_groups, reference_group)
sec_title <- function(lv) tools::toTitleCase(lv)

redraw_auc_section <- function(per_fish, lv, sig_table = NULL,
                               fill_alpha = 0.6,
                               box_lwd    = 1.4,
                               box_wex    = 0.5,     # half-width of each box
                               grp_gap    = 0.6,     # gap between groups
                               show_points= TRUE,
                               point_cex  = 1.05,
                               show_sig   = auc_show_sig,
                               show_xaxis = TRUE,    # FALSE: no x labels (add your own externally)
                               sig_alpha  = auc_sig_alpha,
                               main       = sec_title(lv),
                               ylab       = auc_ylab,
                               cex_lab    = 0.95,
                               cex_axis   = 1.3,
                               y_tick_by  = NULL,    # NULL = auto (~6 ticks); or set a fixed spacing
                               y_digits   = 3,
                               group_order_ = auc_group_order) {
  sp_ord <- c(intersect(group_order_, unique(per_fish$group)),
              setdiff(unique(per_fish$group), group_order_))
  centres <- (seq_along(sp_ord) - 1) * (2*box_wex + grp_gap) + 1
  names(centres) <- sp_ord
  df <- per_fish[per_fish$level == lv, , drop = FALSE]
  if (!nrow(df)) stop("no data for level: ", lv)

  dr0 <- diff(range(df$auc, na.rm = TRUE))
  yl  <- c(min(df$auc, na.rm = TRUE) - 0.05*dr0, max(df$auc, na.rm = TRUE) + 0.09*dr0)
  xl  <- c(centres[1] - (box_wex + grp_gap/2), centres[length(centres)] + (box_wex + grp_gap/2))
  if (is.null(y_tick_by)) y_tick_by <- signif(dr0 / 6, 1)
  if (!is.finite(y_tick_by) || y_tick_by <= 0) y_tick_by <- 1

  op <- par(mar = c(if (isTRUE(show_xaxis)) 5.5 else 1.5, 5.6, 3, 1),
            mgp = c(3.6, 0.8, 0)); on.exit(par(op))
  plot(NA, xlim = xl, ylim = yl, xaxt = "n", yaxt = "n", xlab = "", ylab = "",
       main = main, cex.main = 1.2)
  title(ylab = ylab, line = 4.0, cex.lab = 1.25)
  yticks <- seq(floor(yl[1]/y_tick_by)*y_tick_by, ceiling(yl[2]/y_tick_by)*y_tick_by, by = y_tick_by)
  axis(2, at = yticks, labels = formatC(yticks, format = "f", digits = y_digits),
       las = 1, cex.axis = cex_axis)

  for (gr in sp_ord) {
    v <- df$auc[df$group == gr]
    if (!length(v)) next
    xc <- centres[[gr]]
    boxplot(v, at = xc, add = TRUE, boxwex = 2*box_wex, axes = FALSE,
            col = adjustcolor(group_col[[gr]], fill_alpha), border = "black",
            medlwd = box_lwd + 0.6, lwd = box_lwd, whisklty = 1,
            staplewex = 0.5, outline = FALSE)
    if (show_points) {
      xj <- xc + runif(length(v), -box_wex*0.45, box_wex*0.45)
      points(xj, v, pch = 21, bg = adjustcolor(group_col[[gr]], 0.9),
             col = "black", cex = point_cex, lwd = 0.5)
    }
  }

  if (isTRUE(show_xaxis)) {
    axis(1, at = centres, labels = FALSE)
    mtext(sp_ord, side = 1, line = 1.2, at = centres, cex = cex_lab)
  }

  # significance symbol above each non-reference group's box (vs reference, this level)
  if (show_sig && !is.null(sig_table)) {
    star <- function(p) if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
    s <- sig_table[sig_table$level == lv, , drop = FALSE]
    for (i in seq_len(nrow(s))) {
      gr <- sub("^.* vs ", "", s$comparison[i]); p <- s$p_holm[i]
      if (is.na(p) || p >= sig_alpha || !(gr %in% sp_ord)) next
      v <- df$auc[df$group == gr]; if (!length(v)) next
      text(centres[[gr]], max(v, na.rm = TRUE) + 0.03*dr0, star(p),
           cex = 1.3, font = 2, adj = c(0.5, 0))
    }
  }
  invisible(NULL)
}

for (lv in levels_present) {
  fn <- file.path(out_dir, sprintf("void_auc_%s_%s.tif", auc_axis_mode, lv))
  tiff(fn, width = auc_fig_w, height = auc_fig_h, res = auc_fig_dpi, compression = "lzw")
  redraw_auc_section(per_fish_auc, lv, sig_table = res_auc)
  dev.off()
  message("wrote ", fn)
}

message("Done. Outputs in: ", normalizePath(out_dir))