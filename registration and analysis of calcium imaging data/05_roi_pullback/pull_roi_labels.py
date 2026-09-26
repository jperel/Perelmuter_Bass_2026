#!/usr/bin/env python3
"""
pull_roi_labels.py
=====================
Stage 5 of 6 (roi_pullback), step 1.

For each included native functional plane, pulls the template region labels
backward onto that plane's own native pixel grid.

Geometric method ("slab" technique)
------------------------------------
A native functional plane is physically flat (it is one optical section), but it does
NOT map to a single flat z-slice in structural/template space: the stage-02 functional
<->structural affine (``fwd_final_affine.mat``) contains a genuine ~19-degree rotation
mixing the Y and Z axes, so each native plane maps to a *tilted* plane in reference
space, spreading across many structural z-indices from one edge of the frame to the
other. A naive "nearest matching z-slice, flat 2D crop" lookup would therefore be wrong,
worst near the edges of the frame.

The correct approach, used here: embed each native plane as a real, 1-voxel-thick 3D
image at its own true (x, y, z) position in the shared coordinate frame (the same frame
stage 02's registration was seeded in), then pull the labels backward through the full
composed inverse transform chain in a single ``ants.apply_transforms`` call:

    template -> [stage 03 inverse: invert fwd_0GenericAffine.mat, apply fwd_1InverseWarp]
             -> structural -> [stage 02 inverse: invert fwd_final_affine.mat]
             -> native functional pixel grid

using ``interpolator="genericLabel"`` (nearest-neighbor) throughout, since these are
integer region IDs that must never be linearly interpolated.

``seed_origin`` is recomputed via ``compute_seed_origin``, the exact same function
stage 02's own registration script used to seed its coordinate frame (imported directly,
not transcribed), so it is guaranteed to reproduce the same origin the fitted transform
actually operates against.

Excluded planes
----------------
The 4 shallowest native planes (z = 2440-2470) are excluded: they fall in a stage-02
depth range with a documented residual correction that was never fully resolved into an
invertible transform, and they mostly overlap planes already excluded from the suite2p
analysis in stage 06 for insufficient signal.

Reads (relative to DATA_ROOT):
    01_suite2p_preprocessing/mean_images/meanImg_z<depth>.tif   (36 planes)
    05_roi_pullback/Annotations_V3b_final_fixed_header.nii.gz  (output of fix_annotation_header.py)
    02_functional_to_structural_registration/stage2_registration/fwd_final_affine.mat
    02_functional_to_structural_registration/stage1_ministack/functional_ministack.nii.gz
    02_functional_to_structural_registration/stage2_registration/stageB_plane_corrections.csv
    03_structural_to_template_registration/stage1_registration/fwd_0GenericAffine.mat
    03_structural_to_template_registration/stage1_registration/fwd_1InverseWarp.nii.gz

Writes:
    05_roi_pullback/roi_labels/roi_labels_z<depth>.tif   (32 planes, uint16 region-ID images)

Run after ``fix_annotation_header.py`` and after stages 02/03 have produced their
transform files; run before ``assign_cells_to_regions.py``.

Requires the ``antspy`` environment: ants (ANTsPy), nibabel, numpy, tifffile. Also
imports ``compute_seed_origin`` from stage 02's
``register_ministack_to_structural.py`` (sibling stage, resolved by relative path
below).
"""

import csv
import glob
import os
import re
import sys

import ants
import nibabel as nib
import numpy as np
import tifffile

DATA_ROOT = os.environ.get("PIPELINE_DATA_ROOT") or os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data_for_upload"))

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "02_functional_to_structural_registration")))
from register_ministack_to_structural import compute_seed_origin  # noqa: E402

MEAN_IMG_DIR = os.path.join(DATA_ROOT, "01_suite2p_preprocessing", "mean_images")
LABELS_PATH = os.path.join(DATA_ROOT, "05_roi_pullback", "Annotations_V3b_final_fixed_header.nii.gz")

STAGE_B_AFFINE = os.path.join(DATA_ROOT, "02_functional_to_structural_registration",
                               "stage2_registration", "fwd_final_affine.mat")
STAGE_B_MINISTACK = os.path.join(DATA_ROOT, "02_functional_to_structural_registration",
                                  "stage1_ministack", "functional_ministack.nii.gz")

STAGE_A_DIR = os.path.join(DATA_ROOT, "03_structural_to_template_registration", "stage1_registration")
STAGE_A_AFFINE = os.path.join(STAGE_A_DIR, "fwd_0GenericAffine.mat")
STAGE_A_INVWARP = os.path.join(STAGE_A_DIR, "fwd_1InverseWarp.nii.gz")

STAGE_B_CORRECTIONS_CSV = os.path.join(DATA_ROOT, "02_functional_to_structural_registration",
                                        "stage2_registration", "stageB_plane_corrections.csv")

OUT_DIR = os.path.join(DATA_ROOT, "05_roi_pullback", "roi_labels")


def load_stageB_corrections(path):
    """depth -> (extra_x_um, extra_y_um), the validated per-plane stage-02 residual
    correction (see build_stageB_plane_corrections.py) -- (0, 0) for planes where the
    correction was measured to make things worse and was skipped."""
    corrections = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            corrections[int(row["depth"])] = (float(row["extra_x_um"]), float(row["extra_y_um"]))
    return corrections


XY_SPACING_UM = 1.40625
Z_SPACING_UM = 10.0
DIRECTION = np.array([[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]])
EXCLUDED_IDX = {0, 1, 2, 3}  # z=2440, 2450, 2460, 2470


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    mini = nib.load(STAGE_B_MINISTACK)
    mini_spacing = np.array(mini.header.get_zooms()[:3])
    seed_origin = compute_seed_origin(mini.shape, mini_spacing)
    print(f"seed_origin (um): {seed_origin}")

    stageB_corrections = load_stageB_corrections(STAGE_B_CORRECTIONS_CSV)
    n_corrected = sum(1 for v in stageB_corrections.values() if v != (0.0, 0.0))
    print(f"Loaded Stage B corrections: {n_corrected}/{len(stageB_corrections)} planes corrected")

    labels = ants.image_read(LABELS_PATH)

    paths = sorted(glob.glob(os.path.join(MEAN_IMG_DIR, "meanImg_z*.tif")))
    assert len(paths) == 36, f"expected 36 planes, found {len(paths)}"

    transformlist = [STAGE_B_AFFINE, STAGE_A_AFFINE, STAGE_A_INVWARP]
    whichtoinvert = [True, True, False]

    coverage_log = []
    for i, p in enumerate(paths):
        m = re.search(r"z(\d+)\.tif$", p)
        depth = int(m.group(1))

        if i in EXCLUDED_IDX:
            print(f"plane {i} (z={depth}): EXCLUDED (shallow-z range, unresolved patch)")
            continue

        img_yx = tifffile.imread(p)
        rotated_yx = np.rot90(img_yx, k=1)  # matches build_functional_ministack.py exactly
        plane_xy = rotated_yx.T.astype(np.float32)[:, :, np.newaxis]

        z_i_um = float(seed_origin[2] + i * Z_SPACING_UM)
        # ants.image_read()'s RAS->LPS conversion negates the X,Y ORIGIN too, not just
        # the direction matrix (confirmed empirically: a round-trip test wrote
        # origin_um=(30,40,40) via nibabel and ants.image_read() reported back
        # (-0.03,-0.04,0.04), not (0.03,0.04,0.04)). ants.from_numpy() does no such
        # conversion, so X,Y must be negated by hand here to match; Z is untouched
        # (confirmed empirically too).
        extra_x_um, extra_y_um = stageB_corrections.get(depth, (0.0, 0.0))
        plane_origin_um = (-float(seed_origin[0]) + extra_x_um, -float(seed_origin[1]) + extra_y_um, z_i_um)
        plane_origin = tuple(v / 1000.0 for v in plane_origin_um)  # ITK-internal scale
        plane_spacing = tuple(v / 1000.0 for v in (XY_SPACING_UM, XY_SPACING_UM, Z_SPACING_UM))

        plane_slab = ants.from_numpy(
            plane_xy, origin=plane_origin, spacing=plane_spacing, direction=DIRECTION,
        )

        pulled = ants.apply_transforms(
            fixed=plane_slab, moving=labels,
            transformlist=transformlist, whichtoinvert=whichtoinvert,
            interpolator="genericLabel",
        )
        pulled_arr = pulled.numpy()[:, :, 0]  # (X, Y)

        # undo rotation: (X,Y) -> (Y,X) -> rotate 90 CW to match the original mean image orientation
        back_yx = pulled_arr.T
        original_orientation = np.rot90(back_yx, k=-1)

        assert original_orientation.shape == img_yx.shape, (original_orientation.shape, img_yx.shape)

        coverage = (original_orientation > 0).mean()
        coverage_log.append((i, depth, coverage))
        print(f"plane {i} (z={depth}): coverage={coverage:.3f}")

        out_path = os.path.join(OUT_DIR, f"roi_labels_z{depth}.tif")
        tifffile.imwrite(out_path, original_orientation.astype(np.uint16))

    print("\n=== Coverage summary ===")
    for i, depth, cov in coverage_log:
        flag = "  <-- LOW COVERAGE" if cov < 0.5 else ""
        print(f"plane {i} z={depth}: {cov:.3f}{flag}")

    print(f"\nWrote {len(coverage_log)} label images to {OUT_DIR}")


if __name__ == "__main__":
    main()
