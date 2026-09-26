#!/usr/bin/env python3
"""
correct_interplane_drift.py
==============================
Stage 2 of 6 (functional-to-structural registration) -- Step 1a, run before
build_functional_ministack.py.

Suite2p's own registration corrects XY motion WITHIN each of the 36
single-plane movies (across that movie's frames), but each plane was acquired
as an independent imaging session, so there is no correction for XY stage
drift ACROSS planes. Naively stacking the 36 mean images therefore assumes
they are already XY-aligned to each other, which is not true, and this
incoherence degrades the registration to the structural stack if left
uncorrected.

Fix: sequential pairwise 2D Rigid registration between adjacent planes (10um
apart, so content is very similar and registration is reliable), anchored at
the MIDDLE plane and propagated outward in both directions. This halves the
worst-case accumulated-error chain length (~18 steps) compared to sweeping
from one edge (~35 steps).

Applies the same 90-degree CCW rotation used in build_functional_ministack.py
so the corrected planes are already in final orientation, and writes them to
aligned_planes/, along with a per-plane valid-data mask (aligned_planes_valid_mask/,
consumed by build_ministack_mask.py) and a log of the correction found at each
plane. The log is itself useful QC: a physically real stage-drift correction
should be small and vary smoothly plane-to-plane, not erratically.

Reads:
  - data_for_upload/01_suite2p_preprocessing/mean_images/meanImg_z*.tif (36 files)

Writes:
  - data_for_upload/02_functional_to_structural_registration/stage1_ministack/
    aligned_planes/meanImg_z*.tif
  - data_for_upload/02_functional_to_structural_registration/stage1_ministack/
    aligned_planes_valid_mask/meanImg_z*.tif
  - data_for_upload/02_functional_to_structural_registration/stage1_ministack/
    interplane_drift_log.txt

Environment: antspy env (ants, numpy, tifffile).
"""

import glob
import os
import re

import ants
import numpy as np
import tifffile

DATA_ROOT = os.environ.get("PIPELINE_DATA_ROOT") or os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data_for_upload"))

MEAN_IMG_DIR = os.path.join(DATA_ROOT, "01_suite2p_preprocessing", "mean_images")
STAGE1_DIR = os.path.join(DATA_ROOT, "02_functional_to_structural_registration", "stage1_ministack")
OUT_DIR = os.path.join(STAGE1_DIR, "aligned_planes")
MASK_DIR = os.path.join(STAGE1_DIR, "aligned_planes_valid_mask")
LOG_PATH = os.path.join(STAGE1_DIR, "interplane_drift_log.txt")

XY_SPACING_UM = 1.40625


def load_rotated_planes():
    paths = sorted(glob.glob(os.path.join(MEAN_IMG_DIR, "meanImg_z*.tif")))
    assert len(paths) == 36, f"expected 36 planes, found {len(paths)}"
    depths, planes = [], []
    for p in paths:
        m = re.search(r"z(\d+)\.tif$", p)
        depths.append(int(m.group(1)))
        img_yx = tifffile.imread(p)
        rotated_yx = np.rot90(img_yx, k=1)  # 90 deg CCW, matches build_functional_ministack.py
        planes.append(rotated_yx.astype(np.float32))
    return depths, planes, [os.path.basename(p) for p in paths]


def register_pair(fixed_np, moving_np):
    """Returns (dx, dy in px, corrected moving array) via 2D Rigid registration."""
    fixed = ants.from_numpy(fixed_np, spacing=(XY_SPACING_UM, XY_SPACING_UM))
    moving = ants.from_numpy(moving_np, spacing=(XY_SPACING_UM, XY_SPACING_UM))
    reg = ants.registration(fixed=fixed, moving=moving, type_of_transform="Rigid", aff_metric="GC")
    corrected = reg["warpedmovout"].numpy()
    tx = ants.read_transform(reg["fwdtransforms"][0])
    params = tx.parameters  # [r00,r01,r10,r11, tx, ty] for a 2D rigid/affine transform
    return params, corrected


def build_valid_mask(corrected_np, orig_shape):
    """A drift-corrected plane's footprint after Rigid resampling leaves empty (zero)
    borders wherever content shifted away from. The correct way to build the matching
    binary valid-data mask is to resample a full-frame mask of ones through the exact
    same transform used for the intensity plane; simpler and equally correct here:
    threshold on nonzero, since the resampled background is exactly 0 and real
    microscopy signal is never exactly 0 after linear interpolation."""
    return (corrected_np != 0).astype(np.uint8)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(MASK_DIR, exist_ok=True)
    depths, planes, names = load_rotated_planes()
    n = len(planes)
    mid = n // 2
    print(f"n={n} planes, anchor (uncorrected) = idx {mid} (z={depths[mid]})")

    corrected = [None] * n
    corrected[mid] = planes[mid]
    log_lines = [f"{mid}\t{depths[mid]}\t{names[mid]}\tANCHOR (uncorrected)"]

    # propagate downward (mid-1 .. 0), each registered to the already-corrected neighbor above
    for i in range(mid - 1, -1, -1):
        params, corrected[i] = register_pair(fixed_np=corrected[i + 1], moving_np=planes[i])
        log_lines.append(f"{i}\t{depths[i]}\t{names[i]}\tparams={list(np.round(params, 4))}")
        print(f"idx {i} (z={depths[i]}) registered to idx {i+1}: tx,ty(px)~{params[4]/XY_SPACING_UM:.2f},{params[5]/XY_SPACING_UM:.2f}")

    # propagate upward (mid+1 .. n-1), each registered to the already-corrected neighbor below
    for i in range(mid + 1, n):
        params, corrected[i] = register_pair(fixed_np=corrected[i - 1], moving_np=planes[i])
        log_lines.append(f"{i}\t{depths[i]}\t{names[i]}\tparams={list(np.round(params, 4))}")
        print(f"idx {i} (z={depths[i]}) registered to idx {i-1}: tx,ty(px)~{params[4]/XY_SPACING_UM:.2f},{params[5]/XY_SPACING_UM:.2f}")

    for i, name in enumerate(names):
        out_path = os.path.join(OUT_DIR, name)
        tifffile.imwrite(out_path, corrected[i].astype(np.float32))

        mask = np.ones_like(corrected[i], dtype=np.uint8) if i == mid else build_valid_mask(corrected[i], planes[i].shape)
        tifffile.imwrite(os.path.join(MASK_DIR, name), mask)

    with open(LOG_PATH, "w") as f:
        f.write("\n".join(log_lines) + "\n")

    print(f"\nWrote {n} aligned planes to {OUT_DIR}")
    print(f"Wrote drift log to {LOG_PATH}")


if __name__ == "__main__":
    main()
