###############################################################################
# Area-fraction depth profile (two-group comparison) + per-animal decline
# metrics + group-level statistical test
# -----------------------------------------------------------------------------
# WHAT THIS SCRIPT DOES
#   Compares two groups on how cell density changes
#   from the outer edge of a traced section to its center. For each of one or
#   more "levels" per animal (e.g. different anatomical subregions imaged
#   separately), it:
#
#   1. Loads each animal's traced section outline + detected cell positions.
#   2. Computes an AREA-FRACTION depth profile per animal: relative cell
#      density (observed / expected-if-uniform) as a function of depth, where
#      depth 0 = the section boundary and depth 1 = its center, measured as
#      the fraction of the section's AREA at least that close to a boundary.
#      This is computed analytically via polygon erosion (no pixel grid), and
#      makes profiles from differently sized/shaped sections comparable.
#   3. Averages profiles within each group x level and plots them with 95%
#      confidence bands       -> depth_profile_areafraction.tif
#      and saves every animal's individual curve for reuse
#                              -> depth_profile_data.rds
#   4. Reduces each animal's curve to two summary numbers - slope, and the
#      inner-third/outer-third ("center_edge") density ratio - one row per
#      animal                 -> decline_metrics.csv
#   5. Tests those two numbers, reference group vs the other group, within
#      each level (Wilcoxon rank-sum + effect size + Hodges-Lehmann shift +
#      Holm-corrected p-values)
#                              -> decline_test_results.csv
#
# GETTING STARTED (read this before editing CONFIG below)
#
#   Requirements: R (any recent version, 4.x+) and the "spatstat" package.
#   Install spatstat once, from the R console, with:
#       install.packages("spatstat")
#
#   Your input data: one CSV pair per animal, per group, per level. Every
#   animal needs BOTH of these, in the same folder, sharing the same <id>:
#     <id>_outline_xy.csv   - the traced section outline, as an ORDERED list
#                              of (x, y) points walking around the boundary
#                              (closed or open polygon, either is fine)
#     <id>_centroids.csv    - one row per detected cell, with its (x, y)
#                              position in the SAME coordinate system as the
#                              outline
#   Both files need x/y columns. Accepted column-name variants are listed in
#   .xcand / .ycand further down (x_microns, x_um, x, X - and the y
#   equivalents); rename your columns to one of those, or add your own name
#   to the list. Coordinates can be in any consistent physical unit (microns
#   are assumed by the axis labels / cut_in_um default - relabel those if you
#   use something else).
#
#   Example outline_xy.csv (an open polygon; 4 shown, you'll have more):
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
#       |-- Control_code/
#       |   |-- Rostral/
#       |   |   `-- output/
#       |   |       |-- animal01_outline_xy.csv
#       |   |       |-- animal01_centroids.csv
#       |   |       |-- animal02_outline_xy.csv
#       |   |       |-- animal02_centroids.csv
#       |   |       `-- ...
#       |   `-- Caudal/
#       |       `-- output/  (same pattern)
#       `-- Treatment_code/
#           |-- Rostral/output/  (same pattern)
#           `-- Caudal/output/   (same pattern)
#
#   If your data isn't organised this way, either reorganise it to match, or
#   edit the `groups` list built in the CONFIG section below to point
#   `folder` directly at wherever each group/level's CSV pairs live.
#
# Boundary definition
#   "Boundary" = the nearest point on the traced outline. If your outline
#   traces more than one kind of edge (e.g. an outer rim AND an internal
#   midline, both meaningful as "boundary" for your biology), trace both as
#   part of the same polygon so distance-to-NEAREST-edge captures that.
#
# Decline metrics (per animal, per level)
#   slope        OLS coefficient of relative density on area-fraction depth,
#                fit over the whole 0-1 profile. MORE NEGATIVE = steeper
#                decline from boundary to center.
#   center_edge  mean relative density in the inner third (depth >= 2/3) /
#                mean relative density in the outer third (depth <= 1/3).
#                SMALLER = emptier center relative to the edge.
#   The two are computed from the SAME profile and are typically strongly
#   correlated - they are one finding measured two ways, not independent
#   confirmation. This script reports their correlation so that is visible
#   rather than assumed (see the console "correlation" output when you run
#   it).
#
# Statistical test
#   Wilcoxon rank-sum (Mann-Whitney), reference group vs the other group,
#   within each level, for each metric. Effect size is the rank-biserial
#   correlation (= Cliff's delta), in [-1, 1]; NEGATIVE means the reference
#   group has the LOWER value. Location is reported as the Hodges-Lehmann
#   shift (reference minus other) with its 95% CI. P-values get a Holm
#   correction (see `holm_family` below for the multiplicity scope). This is
#   a descriptive/exploratory test, not confirmatory - treat p-values as
#   indicative, especially with a handful of animals per group.
#
# Required R packages: spatstat (everything else is base R)
###############################################################################

library(spatstat)

## ============================================================================
## ---- CONFIG: EDIT FOR YOUR DATA, THEN RUN THE WHOLE SCRIPT -----------------
## ============================================================================

# ---- Input data --------------------------------------------------------------
DATA_DIR <- "PATH/TO/YOUR/PROJECT/MasterData"   # <-- edit: root of the data tree (see the folder layout above)

# Folder code -> group label. The label is just what appears in your output
# files/figures/legends - it can be a species, genotype, treatment, anything.
# The FIRST entry is treated as the REFERENCE group for every statistical
# comparison below ("reference vs other").
GCODES <- c(Control_code   = "Control",
            Treatment_code = "Treatment")

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

# ---- Output location ----------------------------------------------------------
out_dir <- "PATH/TO/YOUR/PROJECT/Radial Depth Tests"   # <-- edit: where the 4 output files (listed at the top) get written

# ---- Group colours (for the profile figure) ------------------------------------
group_col_user <- c("Control"   = "#0072B2",
                    "Treatment" = "#D55E00")

# ---- Profile-computation settings -----------------------------------------------
cut_in_um <- 5    # erode inward by this many microns (or your coordinate unit)
                  # before sampling cells - use a small positive value to skip
                  # the first, noisiest cell layer right at the traced
                  # boundary; tune by checking the "erosion CDF may be
                  # degenerate" warnings this script prints when it runs.
n_bands   <- 20    # equal-area-fraction bins from boundary (0) to center (1) -
                  # higher = finer profile resolution but noisier per bin
n_cdf     <- 200   # distance samples used to build the analytic area-fraction CDF
ci_level  <- 0.95  # confidence level for the profile-figure bands
fig_dpi <- 600; tiff_compress <- "lzw"

# ---- Decline-metric test settings -----------------------------------------------
metrics <- c("center_edge", "slope")   # tested in this order; must match the
                                        # names returned by decline_metrics()
                                        # below if you add/remove metrics
primary <- "center_edge"               # flagged as "primary" in the console
                                        # report (purely descriptive labelling)

# Holm correction family (multiplicity scope for the p-value correction):
#   "none"   - no correction (report raw p; justify primary/supporting framing)
#   "level"  - within each level, across metrics
#   "metric" - within each metric, across levels
#   "all"    - across every comparison in the table
# "level" treats the metrics as the multiplicity, which is usually the
# defensible reading when levels are separate pre-specified questions and the
# metrics are two summaries of the same underlying profile.
holm_family <- "level"
digits_p <- 4
## ============================================================================
## ---- END CONFIG ----------------------------------------------------------------
## ============================================================================

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
levels_present  <- unique(vapply(groups, function(g) g$level, ""))
group_present   <- unique(vapply(groups, function(g) g$group, ""))
reference_group <- group_present[1]                     # first entry in GCODES
other_groups    <- setdiff(group_present, reference_group)

okabe <- c("#0072B2","#D55E00","#009E73","#CC79A7","#E69F00","#56B4E9","#F0E442","#000000")
group_col <- setNames(as.character(group_col_user[group_present]), group_present)
unfilled <- group_present[is.na(group_col)]
group_col[unfilled] <- setdiff(okabe, group_col[!is.na(group_col)])[seq_along(unfilled)]
if (length(group_present) != 2)
  warning(sprintf("This script is written for exactly 2 groups; found %d.", length(group_present)))

## ---- Loader: reads one animal's outline + centroid CSVs into a spatstat point pattern

.xcand <- c("x_microns","x_um","x","X"); .ycand <- c("y_microns","y_um","y","Y")
.pick <- function(df, cands, what, file){ hit<-cands[cands %in% names(df)]
  if(!length(hit)) stop(sprintf("No %s column in %s (looked for: %s) - add your column name to .xcand/.ycand above",
                                what, basename(file), paste(cands, collapse=", "))); df[[hit[1]]] }
.owin_from_outline <- function(df, file){
  x<-.pick(df,.xcand,"x",file); y<-.pick(df,.ycand,"y",file); n<-length(x)
  if(n>1 && isTRUE(all.equal(x[1],x[n])) && isTRUE(all.equal(y[1],y[n]))){x<-x[-n];y<-y[-n]}  # drop repeated closing vertex, if any
  nx<-c(x[-1],x[1]); ny<-c(y[-1],y[1]); if(0.5*sum(x*ny-nx*y)<0){x<-rev(x);y<-rev(y)}         # force consistent (anticlockwise) winding
  owin(poly=list(x=x,y=y)) }
load_group_patterns <- function(folder, outline_suffix, points_suffix){
  all<-list.files(folder, full.names=TRUE); ofiles<-all[endsWith(basename(all),outline_suffix)]
  if(!length(ofiles)) stop("No files matching '*",outline_suffix,"' found in: ",folder,
                           "\n  Check DATA_DIR / GCODES / LEVELS in CONFIG match your actual folders.")
  pats<-list()
  for(of in ofiles){ b<-basename(of); id<-substr(b,1,nchar(b)-nchar(outline_suffix))
    pf<-file.path(folder,paste0(id,points_suffix)); if(!file.exists(pf)){warning("No matching centroids file for ",id," - skipping");next}
    W<-.owin_from_outline(read.csv(of),of); cdf<-read.csv(pf)
    px<-.pick(cdf,.xcand,"x",pf); py<-.pick(cdf,.ycand,"y",pf); ins<-inside.owin(px,py,W)
    pats[[id]]<-ppp(px[ins],py[ins],window=W) }
  if(!length(pats)) stop("No valid outline/centroid pairs found in: ",folder); pats }

## ---- Per-animal area-fraction-depth profile (robust erosion CDF) ---------------
# Distances are computed via polygon erosion (no rasterisation, exact area in
# your outline's original units). Erosion can occasionally fail or become
# degenerate at large depths on finely-traced/noisy outlines, so a failure is
# treated as a missing sample (NA, interpolated over) rather than assumed to
# mean "fully eroded" - a single early failure could otherwise lock in and
# collapse the whole depth curve to a step function.
profile_one <- function(pp, id = NA_character_){
  W    <- Window(pp)
  W_in <- erosion(W, r = cut_in_um)                         # inward to first real cell layer
  W_in <- simplify.owin(W_in, dmin = 1)                     # thin sub-unit vertices -> stable erosion
  A    <- area(W_in)
  pp_in   <- pp[W_in]
  cell_d  <- bdist.points(pp_in)                            # each cell's distance to interior edge
  if (!length(cell_d) || max(cell_d) <= 0) {
    warning(sprintf("no interior cells for %s after erosion - is cut_in_um too large, or the outline too small?", id))
    return(rep(NA_real_, n_bands)) }

  dmax  <- max(cell_d)
  dgrid <- seq(0, dmax, length.out = n_cdf)
  raw <- vapply(dgrid, function(d) {
    if (d <= 0) return(0)
    Wd <- tryCatch(erosion(W_in, r = d), error = function(e) NULL)
    ad <- if (is.null(Wd)) NA_real_ else suppressWarnings(tryCatch(area(Wd), error=function(e) NA_real_))
    if (!is.finite(ad)) return(NA_real_)                    # erosion FAILURE -> NA (interpolate over), not "fully eroded"
    if (ad <= 0)        return(1)                           # genuinely eroded away -> fully within d of the boundary
    1 - ad / A                                              # fraction of area within d of boundary
  }, numeric(1))

  n_fail <- sum(!is.finite(raw))
  ok <- is.finite(raw)
  frac_within <- if (any(ok)) approx(dgrid[ok], raw[ok], xout = dgrid, rule = 2)$y else rep(NA_real_, length(dgrid))
  frac_within <- cummax(pmin(pmax(frac_within, 0), 1))      # monotone, clamped [0,1]

  # Health check: warn if erosion mostly failed, or the CDF is step-like
  # (jumps near 1 within the first tenth of depth -> cells appear to "pile up"
  # at the center - usually a sign the boundary trace is noisy, or cut_in_um
  # needs to be larger).
  early <- frac_within[max(1, ceiling(0.10 * length(dgrid)))]
  if (n_fail > 0.20 * length(dgrid) || early > 0.90)
    warning(sprintf("erosion CDF may be degenerate for %s (%.0f%% erosion failures, CDF@10%%depth=%.2f) -- try raising cut_in_um, or re-check that outline trace",
                    id, 100 * n_fail / length(dgrid), early))

  cell_af <- approx(dgrid, frac_within, xout = cell_d, rule = 2)$y   # each cell's area-fraction depth
  br <- seq(0, 1, length.out = n_bands + 1)
  counts <- hist(cell_af, breaks = br, plot = FALSE)$counts
  (counts / sum(counts)) / (1 / n_bands)                    # relative density per band (1.0 = uniform)
}

collect <- function(g){
  message(sprintf("[%s | %s] loading: %s", g$group, g$level, g$folder))
  pats <- load_group_patterns(g$folder, outline_suffix, points_suffix)
  pats <- pats[vapply(pats, npoints, 0L) >= 2]              # need >=2 cells to profile an animal
  prof <- sapply(names(pats), function(id) profile_one(pats[[id]], id = id))
  list(group=g$group, level=g$level, prof=prof)
}
prof <- lapply(groups, collect)
get_group <- function(gr,lv){ i<-which(vapply(prof,function(p)p$group==gr&&p$level==lv,TRUE))
  if(length(i)) prof[[i[1]]] else NULL }

## ---- Profile figure -------------------------------------------------------------

xmid <- (seq(0,1,length.out=n_bands+1)[-1] + seq(0,1,length.out=n_bands+1)[-(n_bands+1)]) / 2
curve_summary <- function(mat){
  if(is.null(dim(mat))) mat<-matrix(mat,ncol=1)             # single-animal case
  m<-rowMeans(mat,na.rm=TRUE); nn<-rowSums(is.finite(mat))
  se<-apply(mat,1,function(x){x<-x[is.finite(x)]; if(length(x)<2) NA else sd(x)/sqrt(length(x))})
  tc<-ifelse(nn>1, qt(1-(1-ci_level)/2, df=nn-1), NA); list(mean=m, lo=m-tc*se, hi=m+tc*se) }
band <- function(s,col){ ok<-is.finite(s$lo)&is.finite(s$hi)
  if(any(ok)) polygon(c(xmid[ok],rev(xmid[ok])),c(s$lo[ok],rev(s$hi[ok])),col=adjustcolor(col,0.18),border=NA)
  lines(xmid,s$mean,col=col,lwd=2.5) }

ymax <- max(2, ceiling(max(unlist(lapply(prof, function(p) curve_summary(p$prof)$hi)), na.rm=TRUE)))
tiff(file.path(out_dir, "depth_profile_areafraction.tif"),
     width=5*length(levels_present), height=4.5, units="in", res=fig_dpi, compression=tiff_compress)
op<-par(mfrow=c(1,length(levels_present)), mar=c(4.5,4.5,3,1)); first<-TRUE
for(lv in levels_present){
  plot(NA, xlim=c(0,1), ylim=c(0,ymax), main=lv,
       xlab="fraction of interior (0 = boundary, 1 = center)",
       ylab="relative cell density (obs / uniform)")
  abline(h=1, col="grey70", lty=3)
  # shade + mark the two thirds compared by center_edge (edge = outer 1/3, center = inner 1/3)
  rect(0,   0, 1/3, ymax, col=adjustcolor("#E6550D",0.06), border=NA)
  rect(2/3, 0, 1,   ymax, col=adjustcolor("#3182BD",0.06), border=NA)
  abline(v=c(1/3, 2/3), col="grey55", lty=2)
  mtext(c("edge","center"), side=3, line=-1.2, at=c(1/6,5/6), cex=0.8, col="grey45")
  for(gr in group_present){ gg<-get_group(gr,lv)
    if(!is.null(gg)) band(curve_summary(gg$prof), group_col[[gr]]) }
  if(first){ legend("topright", names(group_col), col=unlist(group_col), lwd=2.5, bty="n"); first<-FALSE }
}
par(op); dev.off()
message("wrote ", file.path(out_dir, "depth_profile_areafraction.tif"))

## ---- Save per-animal profiles (for re-plotting, or feeding into other scripts) --
saveRDS(list(prof=prof, xmid=xmid, n_bands=n_bands, group_col=group_col,
             levels_present=levels_present, group_present=group_present),
        file.path(out_dir, "depth_profile_data.rds"))

## ---- Per-animal decline metrics --------------------------------------------------

inner <- xmid >= 2/3   # inner third of the profile (toward center)
outer <- xmid <= 1/3   # outer third of the profile (toward boundary)

decline_metrics <- function(y){
  # y = relative density per band, aligned to xmid
  ok <- is.finite(y)
  slope <- if (sum(ok) >= 2) unname(coef(lm(y[ok] ~ xmid[ok]))[2]) else NA
  ce    <- mean(y[outer & ok], na.rm = TRUE)
  ce    <- if (is.finite(ce) && ce > 0) mean(y[inner & ok], na.rm = TRUE) / ce else NA
  c(slope = slope, center_edge = ce)
}

rows <- list()
for (p in prof) {
  mat <- p$prof                              # bands x animals
  if (is.null(dim(mat))) mat <- matrix(mat, ncol=1)
  for (j in seq_len(ncol(mat))) {
    m <- decline_metrics(mat[, j])
    rows[[length(rows)+1]] <- data.frame(group=p$group, level=p$level,
                                         animal=j, t(m))
  }
}
tab <- do.call(rbind, rows)

decline_metrics_csv <- file.path(out_dir, "decline_metrics.csv")
write.csv(tab, decline_metrics_csv, row.names = FALSE)
message("wrote ", decline_metrics_csv)

## ---- Group-level test: reference vs other group, within level -------------------
# Wilcoxon rank-sum + rank-biserial effect size + Hodges-Lehmann shift/CI,
# with a Holm correction (scope set by holm_family above).

message("Reference group: ", reference_group, "   |   compared against: ",
        paste(other_groups, collapse = ", "))

grab <- function(mtr, gr, lv) {
  v <- tab[[mtr]][tab$group == gr & tab$level == lv]
  v[is.finite(v)]
}

## descriptives
desc <- data.frame()
for (mtr in metrics) for (lv in levels_present) for (gr in c(reference_group, other_groups)) {
  v <- grab(mtr, gr, lv); if (!length(v)) next
  desc <- rbind(desc, data.frame(
    metric = mtr, level = lv, group = gr, n = length(v),
    mean = round(mean(v), 4), sd = round(stats::sd(v), 4),
    median = round(stats::median(v), 4),
    q25 = round(unname(stats::quantile(v, .25)), 4),
    q75 = round(unname(stats::quantile(v, .75)), 4),
    min = round(min(v), 4), max = round(max(v), 4)))
}

## tests
res <- data.frame()
for (mtr in metrics) for (lv in levels_present) for (gr in other_groups) {
  a <- grab(mtr, reference_group, lv)     # reference
  b <- grab(mtr, gr, lv)                  # other group
  if (length(a) < 2 || length(b) < 2) next   # need >=2 animals per side to test

  ties <- any(duplicated(c(a, b)))
  w    <- suppressWarnings(wilcox.test(a, b, conf.int = TRUE))
  U    <- unname(w$statistic)             # U for the FIRST argument (reference)
  rb   <- 2 * U / (length(a) * length(b)) - 1

  res <- rbind(res, data.frame(
    metric        = mtr,
    role          = if (identical(mtr, primary)) "primary" else "supporting",
    level         = lv,
    comparison    = paste(reference_group, "vs", gr),
    n_ref         = length(a),
    n_other       = length(b),
    ref_mean      = round(mean(a), 4),
    other_mean    = round(mean(b), 4),
    ref_median    = round(stats::median(a), 4),
    other_median  = round(stats::median(b), 4),
    U             = U,
    rank_biserial = round(rb, 3),
    hl_shift      = round(unname(w$estimate), 4),   # reference - other
    hl_lo         = round(w$conf.int[1], 4),
    hl_hi         = round(w$conf.int[2], 4),
    exact         = !ties,
    p_value       = signif(w$p.value, digits_p)))
}
if (!nrow(res))
  stop("No comparison had >= 2 finite values in both groups - not enough animals to test. ",
       "Check decline_metrics.csv (just written) to see what data you actually have.")

## Holm correction
if (identical(holm_family, "none")) {
  res$p_holm <- NA_real_
  fam_note <- "no correction applied (report raw p; see primary/supporting framing)"
} else if (identical(holm_family, "level")) {
  res$p_holm <- ave(res$p_value, res$level, FUN = function(p) p.adjust(p, "holm"))
  fam_note <- "Holm within each level, across metrics"
} else if (identical(holm_family, "metric")) {
  res$p_holm <- ave(res$p_value, res$metric, FUN = function(p) p.adjust(p, "holm"))
  fam_note <- "Holm within each metric, across levels"
} else {
  res$p_holm <- p.adjust(res$p_value, "holm")
  fam_note <- "Holm across all comparisons in the table"
}
res$p_holm <- signif(res$p_holm, digits_p)

## how redundant are the two metrics?
corr <- data.frame()
if (length(metrics) >= 2 && all(c("center_edge","slope") %in% metrics)) {
  for (lv in levels_present) for (gr in c(reference_group, other_groups)) {
    k <- tab$group == gr & tab$level == lv
    x <- tab$center_edge[k]; y <- tab$slope[k]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3) next
    corr <- rbind(corr, data.frame(
      level = lv, group = gr, n = sum(ok),
      spearman_rho = round(unname(stats::cor(x[ok], y[ok], method = "spearman")), 3)))
  }
}

## ---- console report ---------------------------------------------------------------

cat("\n=== Per-group descriptives ===\n")
print(desc, row.names = FALSE)

for (mtr in metrics) {
  rr <- res[res$metric == mtr, ]
  cat("\n=== ", mtr, " (", unique(rr$role), "): Wilcoxon rank-sum vs ", reference_group,
      " ===\n", sep = "")
  print(rr[, c("level","comparison","n_ref","n_other","ref_mean","other_mean",
               "U","rank_biserial","p_value","p_holm")], row.names = FALSE)
  cat("  Hodges-Lehmann shift (reference minus other), 95% CI:\n")
  print(rr[, c("level","hl_shift","hl_lo","hl_hi","exact")], row.names = FALSE)
}

if (nrow(corr)) {
  cat("\n=== center_edge vs slope, within group (Spearman) ===\n")
  print(corr, row.names = FALSE)
  cat("  Strong correlation = the two metrics are largely redundant. Present\n")
  cat("  slope as supporting the same finding, not as a second test of it.\n")
}

cat("\nNotes:\n")
cat("  Correction: ", fam_note, ".\n", sep = "")
if (any(!res$exact)) {
  cat("  Ties present in at least one comparison -> normal approximation with\n",
      "  continuity correction used there (exact = FALSE).\n", sep = "")
}
cat("  rank_biserial < 0 means the reference group has the LOWER value. For\n")
cat("  center_edge that is an emptier center; for slope it is a steeper decline.\n")
cat("  A non-significant p means 'not detected', not 'no difference' - treat\n")
cat("  effect sizes as the primary read-out, especially with few animals/group.\n")

decline_test_csv <- file.path(out_dir, "decline_test_results.csv")
write.csv(res, decline_test_csv, row.names = FALSE)
message("\nWrote: ", normalizePath(decline_test_csv, mustWork = FALSE))
message("Done. Outputs in: ", normalizePath(out_dir))

