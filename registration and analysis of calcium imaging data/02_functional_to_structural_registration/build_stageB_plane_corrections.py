#!/usr/bin/env python3
"""
build_stageB_plane_corrections.py
=====================================
Stage 2 of 6 (functional-to-structural registration).

Builds and validates the accepted per-plane correction table for the
residual functional-vs-structural local misalignment that remains after the
Stage 2 volume-to-volume affine (register_ministack_to_structural.py) and
that varies with depth -- most likely a small residual rotation-axis error in
that fit. This is the correction actually consumed downstream, by stage 05's
pull_roi_labels.py (stageB_plane_corrections.csv).

Measurement method and why it changed: earlier versions of this measurement
used small (80x80) local block phase-correlation, which turned out to be
unreliable -- confirmed directly at one depth (z=2550), where it reported a
tiny 6.9px residual while the true offset (confirmed independently two ways)
was approximately 54px, and similar-scale discrepancies showed up at several
other planes too. The root cause is a real, reusable methodological pitfall:
cells form a somewhat repeating blob pattern, and phase correlation on small
blocks with few distinguishing features can confidently lock onto the wrong
periodic repeat (matching one cell to a similar neighboring cell) rather than
the true shift. This is not resolved by more iterations or a better initial
guess -- it is a property of the measurement primitive itself.

Fix: replace the primitive with skimage.feature.match_template -- a single
large (200x200) central template from the functional channel matched via
normalized cross-correlation against the ENTIRE structural frame (FFT-based,
~40ms/call, no small-window ambiguity to get trapped in). This reproduces and
validates the true z=2550 offset exactly.

For each of the 32 included native planes (planes 0-3 are excluded -- see
EXCLUDED_IDX), an iterative closed-loop procedure:
  1. Measure the current residual (match_template peak displacement) with the
     cumulative correction applied so far.
  2. If converged (residual < CONVERGE_PX) or out of iterations, stop.
  3. Otherwise compute an incremental correction via an empirically calibrated
     2x2 Jacobian M (relating slab-origin shifts to measured pixel
     displacement -- a geometric property of the slab-construction/rotation
     chain; the same -MINV@residual formula, with residual measured via
     match_template, converged 54px -> 7.8px in a single step at z=2550),
     scaled by a damping factor, added to the cumulative correction, then
     re-measured.
  4. Track the cumulative correction giving the SMALLEST residual seen across
     all iterations, which guards against any single iteration overshooting.

Reads:
  - data_for_upload/01_suite2p_preprocessing/mean_images/meanImg_z*.tif
  - data_for_upload/03_structural_to_template_registration/stage0_geometry_fixed/
    02F_stack_ras_cropped.nii.gz
  - data_for_upload/02_functional_to_structural_registration/stage1_ministack/
    functional_ministack.nii.gz (for header geometry only, via compute_seed_origin)
  - data_for_upload/02_functional_to_structural_registration/stage2_registration/
    fwd_final_affine.mat

Writes:
  - data_for_upload/02_functional_to_structural_registration/stage2_registration/
    stageB_plane_corrections.csv (depth, extra_x_um, extra_y_um, n_iterations,
    match_score, before_px, after_px)

Run order: after register_ministack_to_structural.py. Its companion validation
script, diag_stageB_block_displacement.py, can be run afterward as a QC check
across all included planes.

Environment: antspy env (ants, nibabel, numpy, scikit-image, tifffile).
"""

import csv
import os
import sys

import ants
import nibabel as nib
import numpy as np
import tifffile
from skimage.feature import match_template

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from register_ministack_to_structural import compute_seed_origin  # noqa: E402

DATA_ROOT = os.environ.get("PIPELINE_DATA_ROOT") or os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data_for_upload"))

STAGE2_DIR_NAME = "02_functional_to_structural_registration"
MEAN_IMG_DIR = os.path.join(DATA_ROOT, "01_suite2p_preprocessing", "mean_images")
STRUCTURAL_PATH = os.path.join(DATA_ROOT, "03_structural_to_template_registration",
                                "stage0_geometry_fixed", "02F_stack_ras_cropped.nii.gz")
STAGE_B_AFFINE = os.path.join(DATA_ROOT, STAGE2_DIR_NAME, "stage2_registration", "fwd_final_affine.mat")
STAGE_B_MINISTACK = os.path.join(DATA_ROOT, STAGE2_DIR_NAME, "stage1_ministack", "functional_ministack.nii.gz")
OUT_CSV = os.path.join(DATA_ROOT, STAGE2_DIR_NAME, "stage2_registration", "stageB_plane_corrections.csv")

XY_SPACING_UM = 1.40625
Z_SPACING_UM = 10.0
DIRECTION = np.array([[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]])
EXCLUDED_IDX = {0, 1, 2, 3}

M = np.array([[0.46374997, 0.0025], [0.19624999, -0.7125]])
MINV = np.linalg.inv(M)

CROP_R = 100  # half-size of the central template crop (200x200)
MAX_ITERS = 8
DAMPING = 0.85
CONVERGE_PX = 1.5
MIN_SCORE = 0.15  # below this, treat the match as unreliable


def unrot(a):
    return np.rot90(a.T, k=-1)


def get_ch1_ch2(seed_origin, i, depth, structural, extra_origin_um=(0.0, 0.0)):
    p = os.path.join(MEAN_IMG_DIR, f"meanImg_z{depth}.tif")
    img_yx = tifffile.imread(p)
    plane_xy = np.rot90(img_yx, 1).T.astype(np.float32)
    z_um = float(seed_origin[2] + i * Z_SPACING_UM)
    origin_um = (-float(seed_origin[0]) + extra_origin_um[0],
                 -float(seed_origin[1]) + extra_origin_um[1], z_um)
    origin = tuple(v / 1000.0 for v in origin_um)
    spacing = tuple(v / 1000.0 for v in (XY_SPACING_UM, XY_SPACING_UM, Z_SPACING_UM))
    slab = ants.from_numpy(plane_xy[:, :, np.newaxis], origin=origin, spacing=spacing, direction=DIRECTION)
    ch2 = ants.apply_transforms(
        fixed=slab, moving=structural, transformlist=[STAGE_B_AFFINE], whichtoinvert=[True],
        interpolator="linear",
    ).numpy()[:, :, 0]
    ch1 = unrot(plane_xy)
    ch2 = unrot(ch2)
    return ch1, ch2


def measure(ch1, ch2, crop_r=CROP_R):
    h, w = ch1.shape
    cy, cx = h // 2, w // 2
    template = ch1[cy - crop_r:cy + crop_r, cx - crop_r:cx + crop_r].astype(np.float32)
    if template.std() < 1e-3:
        return None
    response = match_template(ch2.astype(np.float32), template, pad_input=True)
    peak = np.unravel_index(np.argmax(response), response.shape)
    score = float(response[peak])
    if score < MIN_SCORE:
        return None
    dy, dx = peak[0] - cy, peak[1] - cx
    return np.array([float(dy), float(dx)]), score


def correct_plane(seed_origin, i, depth, structural):
    ch1, ch2 = get_ch1_ch2(seed_origin, i, depth, structural)
    out0 = measure(ch1, ch2)
    if out0 is None:
        print(f"z={depth}: no reliable match at baseline -- no correction possible")
        return 0.0, 0.0, 0, np.nan, np.nan, np.nan

    residual0, score0 = out0
    before_mag = np.hypot(*residual0)
    cumulative = np.array([0.0, 0.0])
    best_cumulative = cumulative.copy()
    best_mag = before_mag
    best_score = score0

    n_iter = 0
    for it in range(MAX_ITERS):
        ch1c, ch2c = get_ch1_ch2(seed_origin, i, depth, structural, extra_origin_um=tuple(cumulative))
        out = measure(ch1c, ch2c)
        if out is None:
            break
        residual, score = out
        mag = np.hypot(*residual)
        n_iter = it + 1
        if mag < best_mag:
            best_mag = mag
            best_cumulative = cumulative.copy()
            best_score = score
        if mag < CONVERGE_PX:
            break
        delta = -MINV @ residual * DAMPING
        cumulative = cumulative + delta

    print(f"z={depth}: before={before_mag:.2f}px -> after={best_mag:.2f}px "
          f"({n_iter} iterations, score={best_score:.3f}, "
          f"correction=({best_cumulative[0]:.1f},{best_cumulative[1]:.1f})um)")
    return best_cumulative[0], best_cumulative[1], n_iter, before_mag, best_mag, best_score


def main():
    mini = nib.load(STAGE_B_MINISTACK)
    mini_spacing = np.array(mini.header.get_zooms()[:3])
    seed_origin = compute_seed_origin(mini.shape, mini_spacing)

    structural = ants.image_read(STRUCTURAL_PATH)

    paths = sorted(
        [p for p in os.listdir(MEAN_IMG_DIR) if p.startswith("meanImg_z") and p.endswith(".tif")]
    )
    assert len(paths) == 36

    rows = []
    for i, fn in enumerate(paths):
        if i in EXCLUDED_IDX:
            continue
        depth = int(fn.replace("meanImg_z", "").replace(".tif", ""))

        extra_x, extra_y, n_iter, before_mag, after_mag, score = correct_plane(seed_origin, i, depth, structural)
        rows.append({"depth": depth, "extra_x_um": extra_x, "extra_y_um": extra_y,
                     "n_iterations": n_iter, "match_score": score,
                     "before_px": before_mag, "after_px": after_mag})

    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    with open(OUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["depth", "extra_x_um", "extra_y_um", "n_iterations",
                                                "match_score", "before_px", "after_px"])
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    valid = [r for r in rows if not np.isnan(r["after_px"])]
    print(f"\nWrote {OUT_CSV}")
    print(f"Mean residual: before={np.mean([r['before_px'] for r in valid]):.2f}px "
          f"-> after={np.mean([r['after_px'] for r in valid]):.2f}px")
    print(f"Max residual after correction: {max(r['after_px'] for r in valid):.2f}px")
    print(f"Planes still >{CONVERGE_PX}px after correction: "
          f"{sum(1 for r in valid if r['after_px'] >= CONVERGE_PX)}/{len(valid)}")


if __name__ == "__main__":
    main()
