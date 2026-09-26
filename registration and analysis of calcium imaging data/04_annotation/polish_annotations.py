#!/usr/bin/env python3
r"""
polish_annotations.py
======================
Stage 4 of 6 (annotation) in the Danionella telencephalon registration /
calcium-imaging pipeline.

Post-processes the manually traced brain-region label volume in two
computational steps, with a required manual correction step done by hand
in ITK-SNAP in between them (see "PIPELINE ORDER" below).

  (1) GAP FILL   -- assign the thin unlabelled "grout" that survives between
                    adjacent hand-drawn regions to the nearest region, without
                    labelling the empty space around the brain.

  (2) A-P SMOOTH -- remove the between-section ripple that arises when regions
                    are traced section-by-section in the coronal plane, by
                    smoothing each label along the anterior-posterior (tracing)
                    axis only, leaving the in-plane (traced) boundaries intact.

Both steps operate purely in array/voxel space; the output keeps the input's
affine/header byte-for-byte and its dtype (uint16). The input file is never
modified.

--------------------------------------------------------------------------------
PIPELINE ORDER (this script is run twice, with a manual step in between)
--------------------------------------------------------------------------------
The label volume was manually traced in ITK-SNAP section-by-section in the
coronal plane, then smoothed once isotropically in ITK-SNAP ("Smooth Labels").
That is the input this script expects (`Annotations_V3b.nii.gz` in the
accompanying data directory). From there:

  Annotations_V3b.nii.gz
      -> this script, --no-smooth (gap fill only)
      -> [MANUAL: hand-correct in ITK-SNAP -- NOT reproduced by this script;
          see "MANUAL CORRECTION" below]
      -> Annotations_V3b_gapfilled2.nii.gz
      -> this script, --no-fill (A-P smoothing only)
      -> Annotations_V3b_final.nii.gz   (finished annotation)

Running both steps in one pass (no flags) reproduces the final result except
for the manual-correction voxels and their downstream propagation through
smoothing -- it is provided for applying this method to a new annotation from
scratch, not for regenerating the shipped `Annotations_V3b_final.nii.gz`
(which already has the manual correction baked in).

--------------------------------------------------------------------------------
MANUAL CORRECTION (between the two steps -- not automated)
--------------------------------------------------------------------------------
`Annotations_V3b_gapfilled2.nii.gz` (the shipped input to step 2) contains a
hand edit made in ITK-SNAP that this script does not and cannot reproduce,
because it is anatomical judgement, not an algorithm:

  - ~36,830 voxels, one contiguous posterior region that had been traced as
    diencephalon / POA / Vi, were reclassified to Vp.
  - ~2,000 voxels of small boundary tweaks at the Vv/Vs and ENd/Vl borders.

If the annotation is ever re-traced from scratch, this correction must be
redone by hand in ITK-SNAP between step 1 and step 2 -- it is not automated
and this script does not attempt it.

--------------------------------------------------------------------------------
METHOD
--------------------------------------------------------------------------------
(1) Gap fill. The unlabelled voxels (label 0) form one connected network that
    threads between every region and leaks out to the true exterior through gaps
    at the pial surface, so a plain connected-component test cannot separate
    "interior gap" from "outside brain". Instead the label-0 mask is eroded by
    ERODE_K voxels: thin inter-region channels (< ~2*ERODE_K wide) pinch off,
    while the outside world -- which is continuous with the brain surface --
    survives. A flood from the array border through the eroded mask, dilated
    back by ERODE_K+1, is the "outside"; everything else is interior gap. Of
    that, a voxel is filled only if it is genuinely *between annotations*:
    >= 2 distinct region labels lie within HAS2_RADIUS voxels, OR (optionally)
    it is fully enclosed by a single region (binary_fill_holes). Filled voxels
    take the label of the nearest region by *physical* distance (anisotropic
    Euclidean distance transform using SPACING_UM).

(2) Anisotropic label smoothing. Same algorithm as ITK-SNAP's "Smooth Labels"
    (per-label Gaussian on the binary indicator, then per-voxel argmax across
    all labels), but with a per-axis standard deviation: SIGMA_AP voxels along
    the anterior-posterior axis and SIGMA_INPLANE (~0) along the other two. A
    large SIGMA_AP averages out the few-voxel jitter between successive traced
    sections while the label's overall A-P extent is essentially unchanged
    (only the anterior/posterior caps pull in, by roughly SIGMA_AP/2). Only
    voxels within SIGMA_AP+3 of a label in SMOOTH_LABELS are rewritten.

--------------------------------------------------------------------------------
UNITS GOTCHA
--------------------------------------------------------------------------------
The template and every file derived from it store voxel spacing in
millimeters in the NIfTI header (`get_xyzt_units() == ('mm', ...)`); the raw
pixdim value 0.0004874 is 0.4874 **micrometers**, not 0.4874 millimeters as
the bare number might suggest. Any code that converts `header.get_zooms()` to
microns on these files must multiply by 1000. This script sidesteps the
ambiguity entirely by taking voxel spacing as an explicit constant
(SPACING_UM below) rather than trusting the header.

--------------------------------------------------------------------------------
DATA (resolved relative to DATA_ROOT, see below)
--------------------------------------------------------------------------------
Reads/writes under DATA_ROOT/04_annotation/:
  Annotations_V3.nii.gz               original manual tracing (provenance anchor, not consumed by this script)
  Annotations_V3b.nii.gz              input to step 1 (gap fill)
  Annotations_V3b_gapfilled2.nii.gz   step 1 output + manual correction; input to step 2
  Annotations_V3b_final.nii.gz        step 2 output; finished annotation
  labels.txt                          ITK-SNAP label description file (for the summary/QC only)

The intensity template these annotations were drawn on
(`antsBTPtemplate0_fixed_cropped_deghosted2.nii.gz`) is not duplicated here:
it is byte-identical (aside from header) to
DATA_ROOT/03_structural_to_template_registration/stage0_geometry_fixed/template_ras.nii.gz,
shipped with stage 3. Pass that path via --template if you want it underlaid
in the optional QC montages.

--------------------------------------------------------------------------------
BEFORE YOU RUN -- check these for your data
--------------------------------------------------------------------------------
  * AP_AXIS      : which array axis (0, 1 or 2) is anterior-posterior, i.e. the
                   axis you stepped through while tracing.
  * SPACING_UM   : voxel size (um) along (axis0, axis1, axis2).
  * SMOOTH_LABELS: integer label ids whose boundaries should be de-rippled.

Usage:
    # input/output default to DATA_ROOT/04_annotation/Annotations_V3b.nii.gz
    # and .../Annotations_V3b_final.nii.gz (both steps, run end to end)
    python polish_annotations.py

    # explicit paths, either under DATA_ROOT or arbitrary
    python polish_annotations.py  labels_in.nii.gz  labels_out.nii.gz
    python polish_annotations.py  in.nii.gz out.nii.gz --no-smooth
    python polish_annotations.py  in.nii.gz out.nii.gz --qc-dir qc/

To reproduce the manual-correction step in the pipeline order above:
    python polish_annotations.py Annotations_V3b.nii.gz filled.nii.gz --no-smooth
    #   ... hand-correct in ITK-SNAP -> corrected.nii.gz ...
    python polish_annotations.py corrected.nii.gz final.nii.gz --no-fill

DATA_ROOT (below) is the data_for_upload/ directory, resolved by default from
this script's location (../../data_for_upload) or overridden by the
PIPELINE_DATA_ROOT environment variable. --labels and --template default to
paths under DATA_ROOT; --template's default points at stage 3's copy of the
template rather than duplicating that ~1.3 GB file here, since it is
byte-identical (aside from header) to the intensity image these annotations
were drawn on.

Requires: numpy, scipy, nibabel (matplotlib only for --qc-dir).
"""
from __future__ import annotations

import argparse
import gc
import os
import sys

import numpy as np
import nibabel as nib
from scipy import ndimage as ndi

DATA_ROOT = os.environ.get("PIPELINE_DATA_ROOT") or os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data_for_upload"))

# Default I/O: running the script with no arguments performs both steps
# end-to-end on the shipped data. To reproduce the pipeline exactly as it was
# run for the manuscript (with the manual ITK-SNAP correction in between),
# pass explicit input/output paths together with --no-fill / --no-smooth
# (see "PIPELINE ORDER" above).
DEFAULT_INPUT = os.path.join(DATA_ROOT, "04_annotation", "Annotations_V3b.nii.gz")
DEFAULT_OUTPUT = os.path.join(DATA_ROOT, "04_annotation", "Annotations_V3b_final.nii.gz")
DEFAULT_LABELS = os.path.join(DATA_ROOT, "04_annotation", "labels.txt")
# Stage 3's copy of the template, referenced rather than duplicated here;
# byte-identical (aside from header) to the intensity image this annotation
# was traced on.
DEFAULT_TEMPLATE = os.path.join(
    DATA_ROOT, "03_structural_to_template_registration",
    "stage0_geometry_fixed", "template_ras.nii.gz")

# ============================================================================
# CONFIGURATION
# ============================================================================
SPACING_UM = (0.4874, 0.4874, 2.0)   # voxel size (um) along array axes (i, j, k)
AP_AXIS = 1                          # array axis that is anterior-posterior

# --- step 1: gap fill -------------------------------------------------------
DO_FILL = True
ERODE_K = 3            # erosion radius (vox); label-0 channels thinner than
                       # ~2*ERODE_K are treated as interior gaps to fill,
                       # wider ones as outside-brain and left as label 0
HAS2_RADIUS = 4        # a label-0 voxel is "between annotations" if >= 2
                       # distinct region labels lie within this radius (vox)
FILL_ENCLOSED = True   # also fill label-0 fully enclosed by one region

# --- step 2: anisotropic A-P smoothing -------------------------------------
DO_SMOOTH = True
SIGMA_AP = 24.0        # Gaussian SD along AP_AXIS (vox)
SIGMA_INPLANE = 0.6    # Gaussian SD along the other two axes (vox); ~0 leaves
                       # the traced plane pixel-exact
SMOOTH_LABELS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 14, 16, 20]  # label ids to
                       # de-ripple; empty list -> every foreground label

BBOX_MARGIN = 40       # padding (vox) around the labelled bounding box for the
                       # cropped working volume
# ============================================================================


def _log(*a):
    print(*a, flush=True)


def _border_mask(shape):
    m = np.zeros(shape, bool)
    m[0, :, :] = m[-1, :, :] = True
    m[:, 0, :] = m[:, -1, :] = True
    m[:, :, 0] = m[:, :, -1] = True
    return m


def _read_label_names(path):
    names = {0: "Clear"}
    if not path:
        return names
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(None, 7)
            try:
                names[int(parts[0])] = parts[7].strip().strip('"')
            except (ValueError, IndexError):
                pass
    return names


# ---------------------------------------------------------------------------
# STEP 1: gap fill
# ---------------------------------------------------------------------------
def fill_gaps(lab: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (filled_labels, fill_mask) for a cropped uint label volume."""
    fg = lab > 0
    bg = ~fg
    st = ndi.generate_binary_structure(3, 1)          # 6-connected

    # -- separate interior gaps from outside-brain by width --------------
    bg_core = ndi.binary_erosion(bg, structure=st, iterations=ERODE_K,
                                 border_value=1)
    cc, _ = ndi.label(bg_core, structure=np.ones((3, 3, 3)))
    border = _border_mask(bg.shape)
    seeds = np.unique(cc[border & bg_core])
    seeds = seeds[seeds > 0]
    if seeds.size == 0:                               # fallback: largest comp
        sizes = np.bincount(cc.ravel())
        sizes[0] = 0
        seeds = np.array([int(sizes.argmax())])
    outside = np.isin(cc, seeds)
    del cc
    gc.collect()
    outside = ndi.binary_dilation(outside, structure=st,
                                  iterations=ERODE_K + 1) & bg
    interior = bg & ~outside
    del outside
    gc.collect()

    # -- keep only label-0 that is between >= 2 regions (or enclosed) ----
    w = 2 * HAS2_RADIUS + 1
    hi = np.where(fg, lab, np.array(np.iinfo(lab.dtype).max, lab.dtype))
    minl = ndi.minimum_filter(hi, size=w)
    maxl = ndi.maximum_filter(lab, size=w)
    hival = np.iinfo(lab.dtype).max
    between = (minl != hival) & (maxl.astype(np.int64) != minl.astype(np.int64))
    del hi, minl, maxl
    gc.collect()
    if FILL_ENCLOSED:
        enc = ndi.binary_fill_holes(fg)
        for k in range(enc.shape[2]):                 # per-slice catches 2D-
            enc[:, :, k] = ndi.binary_fill_holes(enc[:, :, k])   # enclosed holes
        between |= enc & bg
        del enc
        gc.collect()
    fill = interior & between
    del interior, between
    gc.collect()

    # -- assign nearest region by physical distance ---------------------
    out = lab.copy()
    if fill.any():
        idx = np.argwhere(fill)
        lo = np.maximum(idx.min(0) - BBOX_MARGIN, 0)
        hi_ = np.minimum(idx.max(0) + BBOX_MARGIN + 1, np.array(lab.shape))
        sl = tuple(slice(int(a), int(b)) for a, b in zip(lo, hi_))
        bg_sub = bg[sl]
        near_idx = ndi.distance_transform_edt(
            bg_sub, sampling=SPACING_UM,
            return_distances=False, return_indices=True)
        nearest = lab[sl][tuple(near_idx)]
        fsub = fill[sl]
        out[sl][fsub] = nearest[fsub]
        del bg_sub, near_idx, nearest, fsub
        gc.collect()
    return out, fill


# ---------------------------------------------------------------------------
# STEP 2: anisotropic A-P label smoothing
# ---------------------------------------------------------------------------
def ap_smooth(lab: np.ndarray, target_labels) -> tuple[np.ndarray, np.ndarray]:
    """Return (smoothed_labels, changed_mask) for a cropped uint label volume."""
    present = [int(v) for v in np.unique(lab) if v != 0]
    targets = list(target_labels) if target_labels else list(present)

    sigma = [SIGMA_INPLANE, SIGMA_INPLANE, SIGMA_INPLANE]
    sigma[AP_AXIS] = SIGMA_AP

    best_val = np.full(lab.shape, -1.0, np.float32)
    best_lab = np.zeros(lab.shape, lab.dtype)
    for L in [0] + present:                           # background competes too
        ind = (lab == L).astype(np.float32)
        sm = ndi.gaussian_filter(ind, sigma=sigma, mode="nearest")
        del ind
        upd = sm > best_val
        best_val[upd] = sm[upd]
        best_lab[upd] = L
        del sm, upd
        gc.collect()
    del best_val
    gc.collect()

    reach = int(round(SIGMA_AP)) + 3
    scope = ndi.binary_dilation(np.isin(lab, targets),
                                structure=ndi.generate_binary_structure(3, 1),
                                iterations=reach)
    out = lab.copy()
    out[scope] = best_lab[scope]
    changed = out != lab
    del best_lab, scope
    gc.collect()
    return out, changed


# ---------------------------------------------------------------------------
# QC montages (optional)
# ---------------------------------------------------------------------------
def qc_montages(before, after, names, spacing, ap_axis, out_dir, template=None):
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    os.makedirs(out_dir, exist_ok=True)
    mx = int(max(before.max(), after.max()))
    rng = np.random.default_rng(0)
    cols = rng.random((mx + 1, 3))
    cols[0] = 0
    cmap = ListedColormap(cols)

    axes = [a for a in (0, 1, 2)]
    names_ax = {0: "axis0", 1: "axis1", 2: "axis2"}
    names_ax[ap_axis] = "A-P (traced-through)"
    for view in axes:
        idxs = np.linspace(0.2, 0.8, 4) * before.shape[view]
        fig, ax = plt.subplots(len(idxs), 2, figsize=(9, 3.2 * len(idxs)))
        for r, f in enumerate(idxs):
            k = int(f)
            sl = [slice(None)] * 3
            sl[view] = k
            for c, (vol, ttl) in enumerate(((before, "before"), (after, "after"))):
                img = np.squeeze(vol[tuple(sl)]).T
                a = ax[r, c]
                if template is not None:
                    bg = np.squeeze(template[tuple(sl)]).T
                    a.imshow(bg, cmap="gray", origin="lower",
                             vmax=np.percentile(bg, 99.5))
                    o = cmap(np.clip(img, 0, mx))
                    o[..., 3] = np.where(img > 0, 0.55, 0)
                    a.imshow(o, origin="lower")
                else:
                    a.imshow(cmap(np.clip(img, 0, mx)), origin="lower")
                a.set_title(f"{names_ax[view]}  slice {k}  {ttl}", fontsize=8)
                a.axis("off")
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, f"qc_{names_ax[view].split()[0]}.png"),
                    dpi=90)
        plt.close(fig)
    _log(f"  QC montages -> {out_dir}")


# ---------------------------------------------------------------------------
def _summary(before, after, names):
    ids = sorted(set(np.unique(before)) | set(np.unique(after)))
    _log(f"  {'id':>3}  {'name':<28s} {'before':>12s} {'after':>12s} {'delta':>10s}")
    for i in ids:
        b = int((before == i).sum())
        a = int((after == i).sum())
        _log(f"  {int(i):>3}  {names.get(int(i), '?'):<28s} "
             f"{b:>12d} {a:>12d} {a - b:>+10d}")
    lost = set(np.unique(before)) - set(np.unique(after))
    if lost:
        _log(f"  WARNING: label(s) removed entirely: {sorted(lost)}")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", nargs="?", default=DEFAULT_INPUT,
                   help=f"input label NIfTI, never modified (default: {DEFAULT_INPUT})")
    p.add_argument("output", nargs="?", default=DEFAULT_OUTPUT,
                   help=f"output label NIfTI (default: {DEFAULT_OUTPUT})")
    p.add_argument("--no-fill", action="store_true", help="skip step 1 (gap fill)")
    p.add_argument("--no-smooth", action="store_true",
                   help="skip step 2 (A-P smoothing)")
    p.add_argument("--labels", metavar="FILE", default=DEFAULT_LABELS,
                   help="ITK-SNAP label description file, for the summary only "
                        f"(default: {DEFAULT_LABELS})")
    p.add_argument("--qc-dir", metavar="DIR",
                   help="write before/after montages here (needs matplotlib)")
    p.add_argument("--template", metavar="FILE", default=DEFAULT_TEMPLATE,
                   help="intensity image to underlay in the QC montages, only used "
                        f"if --qc-dir is set (default: stage 3's template_ras.nii.gz)")
    args = p.parse_args(argv)

    do_fill = DO_FILL and not args.no_fill
    do_smooth = DO_SMOOTH and not args.no_smooth
    names = _read_label_names(args.labels)

    _log(f"reading {args.input}")
    src = nib.load(args.input)
    full = np.asanyarray(src.dataobj)
    if full.ndim == 4:
        full = full[..., 0]
    if not np.issubdtype(full.dtype, np.integer):
        if not np.allclose(full, np.round(full)):
            sys.exit("ERROR: input is not an integer label volume")
    full = np.ascontiguousarray(full).astype(np.uint16)
    orig = full.copy()
    _log(f"  grid {full.shape}  labels {sorted(int(v) for v in np.unique(full))}")
    _log(f"  AP_AXIS = {AP_AXIS}   SPACING_UM = {SPACING_UM}")

    # crop to a padded bounding box of the labelled voxels
    nz = np.argwhere(full > 0)
    if nz.size == 0:
        sys.exit("ERROR: no non-zero labels in input")
    lo = np.maximum(nz.min(0) - BBOX_MARGIN, 0)
    hi = np.minimum(nz.max(0) + BBOX_MARGIN + 1, np.array(full.shape))
    sl = tuple(slice(int(a), int(b)) for a, b in zip(lo, hi))
    lab = np.ascontiguousarray(full[sl])
    del nz

    if do_fill:
        _log(f"\n[1/2] gap fill  (ERODE_K={ERODE_K}, HAS2_RADIUS={HAS2_RADIUS}, "
             f"FILL_ENCLOSED={FILL_ENCLOSED})")
        lab, fmask = fill_gaps(lab)
        _log(f"  filled {int(fmask.sum()):,} voxels")
        del fmask
        gc.collect()
    else:
        _log("\n[1/2] gap fill  -- skipped")

    if do_smooth:
        _log(f"\n[2/2] anisotropic A-P smoothing  (SIGMA_AP={SIGMA_AP}, "
             f"SIGMA_INPLANE={SIGMA_INPLANE}, "
             f"labels={SMOOTH_LABELS or 'ALL'})")
        lab, cmask = ap_smooth(lab, SMOOTH_LABELS)
        _log(f"  changed {int(cmask.sum()):,} voxels")
        del cmask
        gc.collect()
    else:
        _log("\n[2/2] A-P smoothing  -- skipped")

    full[sl] = lab
    del lab

    _log("\nper-label voxel counts:")
    _summary(orig, full, names)

    out = nib.Nifti1Image(full, src.affine, src.header)
    out.header.set_data_dtype(np.uint16)
    nib.save(out, args.output)
    _log(f"\nwrote {args.output}")

    if args.qc_dir:
        _log("\nrendering QC montages ...")
        tmpl = None
        if args.template:
            if os.path.exists(args.template):
                t = np.asanyarray(nib.load(args.template).dataobj)
                tmpl = t[..., 0] if t.ndim == 4 else t
            else:
                _log(f"  --template not found, skipping underlay: {args.template}")
        qc_montages(orig, full, names, SPACING_UM, AP_AXIS, args.qc_dir, tmpl)

    _log("\ndone.")


if __name__ == "__main__":
    main()
