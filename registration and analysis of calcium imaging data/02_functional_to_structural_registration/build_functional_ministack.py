#!/usr/bin/env python3
"""
build_functional_ministack.py
================================
Stage 2 of 6 (functional-to-structural registration) -- Step 1.

Stacks the 36 mean images (suite2p's registered mean image per plane -- never
the raw recordings) into one coherent 3D "mini-stack" NIfTI, spacing
(1.40625, 1.40625, 10.0) um, verified directly from OME-XML metadata in
verify_pixel_size.py (all 36 source files agreed exactly).

By default this reads directly from stage 01's mean_images/. Set
USE_ALIGNED_PLANES=1 to instead read from this stage's own
stage1_ministack/aligned_planes/ (produced by correct_interplane_drift.py),
which additionally corrects cross-plane XY stage drift. The functional_ministack.nii.gz
shipped in data_for_upload was built with USE_ALIGNED_PLANES=1, i.e. run
correct_interplane_drift.py first, then this script with that environment
variable set.

Axis order and header convention: NIfTI array axis order (X,Y,Z), identity
direction cosines, origin (0,0,0), xyzt_units="micron" -- i.e. the raw pixdim
values ARE literal microns, with no mm scaling. This matches the convention
used by stage 03's geometry-fixing script; other files elsewhere in the
broader project use a mm-scaled convention instead, so this distinction
matters if cross-referencing NIfTI headers across stages.

The functional camera frame is rotated 90 degrees clockwise relative to the
structural stack's orientation, so each plane is rotated 90 degrees
counter-clockwise (np.rot90, k=1) to match. This also resolves an FOV
aspect-ratio mismatch: the raw functional FOV is 720x540um (landscape) versus
the structural stack's 600x800um (portrait); rotated, it becomes 540x720um,
matching the structural crop's aspect ratio exactly (540/720 == 600/800 ==
0.75).

Suite2p's meanImg arrays are (Y,X) per plane; after the 90-degree CCW rotation
each plane is transposed to (X,Y) before stacking along a new 3rd (Z) axis,
matching this project's NIfTI axis-order convention.

Reads:
  - data_for_upload/01_suite2p_preprocessing/mean_images/meanImg_z*.tif, or
  - data_for_upload/02_functional_to_structural_registration/stage1_ministack/
    aligned_planes/meanImg_z*.tif (if USE_ALIGNED_PLANES=1)

Writes:
  - data_for_upload/02_functional_to_structural_registration/stage1_ministack/
    functional_ministack.nii.gz
  - data_for_upload/02_functional_to_structural_registration/stage1_ministack/
    depths_um.txt

Run order: after correct_interplane_drift.py (if USE_ALIGNED_PLANES=1), before
sanity_check_step1.py and build_ministack_mask.py.

Environment: antspy env (nibabel, numpy, tifffile).
"""

import glob
import os
import re

import nibabel as nib
import numpy as np
import tifffile

DATA_ROOT = os.environ.get("PIPELINE_DATA_ROOT") or os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data_for_upload"))

MEAN_IMG_DIR = os.path.join(DATA_ROOT, "01_suite2p_preprocessing", "mean_images")
STAGE1_DIR = os.path.join(DATA_ROOT, "02_functional_to_structural_registration", "stage1_ministack")
ALIGNED_DIR = os.path.join(STAGE1_DIR, "aligned_planes")
OUT_DIR = STAGE1_DIR
OUT_PATH = os.path.join(OUT_DIR, "functional_ministack.nii.gz")

USE_ALIGNED_PLANES = os.environ.get("USE_ALIGNED_PLANES", "0") == "1"

XY_SPACING_UM = 1.40625  # verified via verify_pixel_size.py OME-XML read
Z_SPACING_UM = 10.0      # verified via verify_pixel_size.py OME-XML read (PhysicalSizeZ)


def build_ras_affine(spacing_xyz):
    sx, sy, sz = (float(v) for v in spacing_xyz)
    return np.diag([sx, sy, sz, 1.0])


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    src_dir = ALIGNED_DIR if USE_ALIGNED_PLANES else MEAN_IMG_DIR
    paths = sorted(glob.glob(os.path.join(src_dir, "meanImg_z*.tif")))
    assert len(paths) == 36, f"expected 36 planes, found {len(paths)} in {src_dir}"
    print(f"Source: {src_dir} (USE_ALIGNED_PLANES={USE_ALIGNED_PLANES})")

    depths = []
    planes_xy = []
    for p in paths:
        m = re.search(r"z(\d+)\.tif$", p)
        depths.append(int(m.group(1)))
        img_yx = tifffile.imread(p)
        if USE_ALIGNED_PLANES:
            # aligned_planes/ already has rotation + inter-plane drift correction applied
            planes_xy.append(img_yx.T)  # (Y,X) -> (X,Y)
        else:
            assert img_yx.shape == (384, 512), f"unexpected shape {img_yx.shape} for {p}"
            rotated_yx = np.rot90(img_yx, k=1)  # 90 deg CCW -- matches structural stack orientation
            planes_xy.append(rotated_yx.T)  # (Y,X) -> (X,Y)

    assert depths == sorted(depths), "planes not in ascending depth order"
    steps = np.diff(depths)
    assert np.all(steps == 10), f"expected uniform 10um steps, got {set(steps)}"
    print(f"Depths: {depths[0]}..{depths[-1]} um, step {steps[0]}um, n={len(depths)}")

    volume_xyz = np.stack(planes_xy, axis=2).astype(np.float32)  # (X,Y,Z) = (512,384,36)
    print(f"Mini-stack array shape (X,Y,Z): {volume_xyz.shape}")

    affine = build_ras_affine((XY_SPACING_UM, XY_SPACING_UM, Z_SPACING_UM))
    out_img = nib.Nifti1Image(volume_xyz, affine)
    out_img.header.set_data_dtype("float32")
    out_img.header.set_xyzt_units(xyz="micron")
    out_img.set_sform(affine, code=1)
    out_img.set_qform(affine, code=1)
    nib.save(out_img, OUT_PATH)
    print(f"Wrote {OUT_PATH}")

    # also save the depth list alongside, so later steps know which mini-stack
    # z-index maps to which native functional-plane filename/z-label
    with open(os.path.join(OUT_DIR, "depths_um.txt"), "w") as f:
        for i, d in enumerate(depths):
            f.write(f"{i}\t{d}\t{os.path.basename(paths[i])}\n")
    print(f"Wrote depth index map to {os.path.join(OUT_DIR, 'depths_um.txt')}")


if __name__ == "__main__":
    main()
