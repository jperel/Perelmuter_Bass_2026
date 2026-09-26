#!/usr/bin/env python3
"""
Stage 3 (of 6) - structural-to-template registration: QC overlay.

Produces two 2-channel TIFF stacks for visual inspection of the registration in FIJI/
ImageJ, viewing the fit from both directions:

1. check_template_vs_warped_2ch_stack.tif -- template (fixed) + structural warped INTO
   template space (forward chain: affine then warp), on the template's own grid, one
   stack of all template planes.

2. check_structural_vs_warped_template_2ch_stack.tif -- the structural stack in its own
   native space (never resampled) + template warped INTO structural space (inverse
   chain: affine inverted, then the inverse warp), limited to the structural stack's own
   planes. This is the complementary view: it shows inverse-warp quality directly on the
   grid that matters for the downstream ROI pull-back stage (05_roi_pullback), which uses
   this same inverse chain to pull template region labels onto native
   functional/structural space, rather than only checking the forward direction on the
   template's grid.

Reads (all under DATA_ROOT/03_structural_to_template_registration/):
  stage0_geometry_fixed/template_ras.nii.gz
  stage0_geometry_fixed/02F_stack_ras_cropped.nii.gz
  stage1_registration/fwd_0GenericAffine.mat
  stage1_registration/fwd_1Warp.nii.gz
  stage1_registration/fwd_1InverseWarp.nii.gz

Writes (under DATA_ROOT/03_structural_to_template_registration/stage1_registration/qc/):
  check_template_vs_warped_2ch_stack.tif
  check_structural_vs_warped_template_2ch_stack.tif

Run order in this stage: fix_template_geometry.py -> register_structural_to_template.py
-> qc_overlay.py (this script).

Environment: antspy env (ants, nibabel, numpy, scipy, scikit-image, tifffile).
"""

import os

import ants
import numpy as np
import tifffile

DATA_ROOT = os.environ.get("PIPELINE_DATA_ROOT") or os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data_for_upload"))

STAGE_DIR = os.path.join(DATA_ROOT, "03_structural_to_template_registration")
GEOM_DIR = os.path.join(STAGE_DIR, "stage0_geometry_fixed")
REG_DIR = os.path.join(STAGE_DIR, "stage1_registration")
OUT_DIR = os.path.join(REG_DIR, "qc")

TEMPLATE_PATH = os.path.join(GEOM_DIR, "template_ras.nii.gz")
STRUCTURAL_PATH = os.path.join(GEOM_DIR, "02F_stack_ras_cropped.nii.gz")

AFFINE_PATH = os.path.join(REG_DIR, "fwd_0GenericAffine.mat")
FWDWARP_PATH = os.path.join(REG_DIR, "fwd_1Warp.nii.gz")
INVWARP_PATH = os.path.join(REG_DIR, "fwd_1InverseWarp.nii.gz")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    print(f"QC for {REG_DIR}")

    template = ants.image_read(TEMPLATE_PATH)
    structural = ants.image_read(STRUCTURAL_PATH)

    # --- 1. structural forward-warped into TEMPLATE space ---
    print("Warping structural -> template space (forward chain)...")
    struct_in_template = ants.apply_transforms(
        fixed=template, moving=structural,
        transformlist=[FWDWARP_PATH, AFFINE_PATH], whichtoinvert=[False, False],
        interpolator="linear",
    )
    t_np = template.numpy()
    s_in_t_np = struct_in_template.numpy()
    nx, ny, nz = t_np.shape
    print(f"template grid shape {t_np.shape}")

    stack1 = np.zeros((nz, 2, ny, nx), dtype="float32")
    for zi in range(nz):
        stack1[zi, 0] = t_np[:, :, zi].T
        stack1[zi, 1] = s_in_t_np[:, :, zi].T
    out1 = os.path.join(OUT_DIR, "check_template_vs_warped_2ch_stack.tif")
    tifffile.imwrite(out1, stack1, imagej=True, metadata={"axes": "ZCYX"})
    print(f"wrote {out1}  ({nz} planes, ch1=template ch2=structural-warped-to-template)")
    del struct_in_template, s_in_t_np, stack1

    # --- 2. template inverse-warped into STRUCTURAL's own native space ---
    print("Warping template -> structural space (inverse chain)...")
    template_in_structural = ants.apply_transforms(
        fixed=structural, moving=template,
        transformlist=[AFFINE_PATH, INVWARP_PATH], whichtoinvert=[True, False],
        interpolator="linear",
    )
    s_np = structural.numpy()
    t_in_s_np = template_in_structural.numpy()
    nx2, ny2, nz2 = s_np.shape
    print(f"structural grid shape {s_np.shape}")

    stack2 = np.zeros((nz2, 2, ny2, nx2), dtype="float32")
    for zi in range(nz2):
        stack2[zi, 0] = s_np[:, :, zi].T
        stack2[zi, 1] = t_in_s_np[:, :, zi].T
    out2 = os.path.join(OUT_DIR, "check_structural_vs_warped_template_2ch_stack.tif")
    tifffile.imwrite(out2, stack2, imagej=True, metadata={"axes": "ZCYX"})
    print(f"wrote {out2}  ({nz2} planes, ch1=structural(native) ch2=template-warped-to-structural)")

    print(f"\nDone. QC outputs in {OUT_DIR}")


if __name__ == "__main__":
    main()
