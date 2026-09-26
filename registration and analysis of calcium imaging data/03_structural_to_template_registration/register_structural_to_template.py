#!/usr/bin/env python3
"""
Stage 3 (of 6) - structural-to-template registration: main registration.

Registers the RAS-converted structural stack to the population-average reference
template using 19 manually-placed landmark pairs to seed a full 12-DOF affine, followed
by deformable (SyN) refinement.

Method:
  1. Load 19 manually-placed landmark pairs (moving = structural stack, fixed =
     template) from landmarks.csv and convert index coordinates to physical points.
  2. Fit a full 12-DOF affine directly from the paired landmarks
     (ants.fit_transform_to_paired_points, transform_type="affine"). This affine is
     saved on its own (landmark_affine.mat) and used as the initial transform for
     everything downstream.
  3. Resample fixed/moving volumes and masks to 4 um isotropic spacing for the
     intensity-based fitting stages (the landmark affine itself is fit at full
     resolution).
  4. Request a Rigid-only refinement on top of the landmark affine. See the CAVEAT
     below: this stage is confirmed to be a no-op for this dataset.
  5. Run SyNOnly deformable refinement (CC metric) on top of that, with total_sigma
     raised to 2.0 (see below).
  6. Apply the composed forward transform to the full-resolution moving volume and
     write it out for inspection (structural_in_template_space.nii.gz).

CAVEAT - the requested Rigid refinement stage is a confirmed no-op for this data, so
the accepted result is in practice landmark-affine + SyN, not landmark + Rigid + SyN:
the landmark-fit affine carries real scale (roughly 1.11-1.15 on the diagonal, since the
structural stack and template are not the same physical size), which breaks the
near-identity assumption that ANTs' Rigid-family optimizers rely on for their first
gradient step. Every Rigid-family preset tested (Rigid, QuickRigid, DenseRigid,
BOLDRigid) throws an internal ITK exception within 1-2 gradient iterations on this
initialization. `ants.registration()` catches that exception internally and silently
returns the unmodified input transform, with no error surfaced to the caller -- so the
call below appears to succeed and returns a "rigid" transform that is actually just the
landmark affine, unchanged. This has been confirmed and is not fixed, because it does
not change the final result (SyN is initialized from the same landmark affine either
way): it is documented here so the fitted transform chain is not misread as containing a
genuine rigid refinement step.

total_sigma=2.0: total_sigma smooths the cumulative SyN displacement field itself after
each update (as opposed to flow_sigma, which only smooths each per-iteration update
before it is applied). At a lower total_sigma, the fitted warp overfit to sharp,
repeatable high-contrast tissue boundaries via the CC metric, producing an
edge-hugging Jacobian artifact: the fraction of the brain with |Jacobian determinant - 1|
> 0.2 was 15.7%. Raising total_sigma to 2.0 caps how sharp the cumulative field can get
regardless of how many iterations accumulate pull toward the same edge, and dropped that
fraction to 0.84%, with the whole-frame midline curvature check also improved.

Reads (all under DATA_ROOT/03_structural_to_template_registration/):
  stage0_geometry_fixed/template_ras.nii.gz        (fixed)
  stage0_geometry_fixed/02F_stack_ras_cropped.nii.gz  (moving)
  stage0_geometry_fixed/mask_ras_cropped.nii.gz       (moving mask)
  landmarks.csv

Writes (all under DATA_ROOT/03_structural_to_template_registration/stage1_registration/):
  landmark_affine.mat
  fwd_0GenericAffine.mat, fwd_1Warp.nii.gz, fwd_1InverseWarp.nii.gz  (SyN outputs)
  structural_in_template_space.nii.gz  (full-resolution QC volume)

Run order in this stage: fix_template_geometry.py -> register_structural_to_template.py
-> qc_overlay.py.

Environment: antspy env (ants, nibabel, numpy, scipy, scikit-image, tifffile).
"""

import csv
import glob
import os
import sys

import numpy as np
import nibabel as nib
import ants

SEED = int(os.environ.get("ANTS_SEED", sys.argv[1] if len(sys.argv) > 1 else 123))
ants.config.set_ants_deterministic(True, seed_value=SEED)

DATA_ROOT = os.environ.get("PIPELINE_DATA_ROOT") or os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data_for_upload"))

STAGE_DIR = os.path.join(DATA_ROOT, "03_structural_to_template_registration")
GEOM_DIR = os.path.join(STAGE_DIR, "stage0_geometry_fixed")
OUT_DIR = os.path.join(STAGE_DIR, "stage1_registration")
LANDMARKS_CSV = os.path.join(STAGE_DIR, "landmarks.csv")

FIXED_PATH = os.path.join(GEOM_DIR, "template_ras.nii.gz")
MOVING_PATH = os.path.join(GEOM_DIR, "02F_stack_ras_cropped.nii.gz")
MOVING_MASK_PATH = os.path.join(GEOM_DIR, "mask_ras_cropped.nii.gz")

WORKING_SPACING_UM = 4.0

LANDMARK_XFM_PATH = os.path.join(OUT_DIR, "landmark_affine.mat")
OUT_XFM_PREFIX = os.path.join(OUT_DIR, "fwd_")
OUT_WARPED_FULLRES = os.path.join(OUT_DIR, "structural_in_template_space.nii.gz")


def clean_prefix(prefix):
    for f in glob.glob(prefix + "*"):
        os.remove(f)
        print(f"  removed stale {f}")


def save_nifti_ras(path, data_xyz, affine, dtype="float32"):
    img = nib.Nifti1Image(np.asarray(data_xyz, dtype=dtype), affine)
    img.header.set_data_dtype(dtype)
    img.header.set_xyzt_units(xyz="micron")
    img.set_sform(affine, code=1)
    img.set_qform(affine, code=1)
    nib.save(img, path)


def load_landmarks(path, moving_img, fixed_img):
    moving_pts, fixed_pts, names = [], [], []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            mv_idx = (int(row["moving_x"]), int(row["moving_y"]), int(row["moving_z"]))
            fx_idx = (int(row["fixed_x"]), int(row["fixed_y"]), int(row["fixed_z"]))
            oob = False
            for ax, (v, dim) in enumerate(zip(mv_idx, moving_img.shape)):
                if not (0 <= v <= dim - 1):
                    print(f"  SKIPPING landmark {row['name']!r}: moving axis {ax} index {v} "
                          f"out of range [0, {dim - 1}]")
                    oob = True
            for ax, (v, dim) in enumerate(zip(fx_idx, fixed_img.shape)):
                if not (0 <= v <= dim - 1):
                    print(f"  SKIPPING landmark {row['name']!r}: fixed axis {ax} index {v} "
                          f"out of range [0, {dim - 1}]")
                    oob = True
            if oob:
                continue
            moving_pts.append(ants.transform_index_to_physical_point(moving_img, mv_idx))
            fixed_pts.append(ants.transform_index_to_physical_point(fixed_img, fx_idx))
            names.append(row["name"])
    return names, np.array(moving_pts), np.array(fixed_pts)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    print(f"Loading fixed  {FIXED_PATH}")
    fixed = ants.image_read(FIXED_PATH)
    print(f"Loading moving {MOVING_PATH}")
    moving = ants.image_read(MOVING_PATH)
    print(f"Loading moving_mask {MOVING_MASK_PATH}")
    moving_mask = ants.image_read(MOVING_MASK_PATH)
    fixed_mask = ants.threshold_image(fixed, 0, 0, 0, 1)

    print(f"Loading landmarks from {LANDMARKS_CSV}")
    names, moving_pts, fixed_pts = load_landmarks(LANDMARKS_CSV, moving, fixed)
    print(f"  {len(names)} landmark pairs: {names}")
    if len(names) < 4:
        raise SystemExit(f"Need at least 4 landmark pairs for an affine fit, got {len(names)}.")

    print("Fitting Affine transform from landmarks (ants.fit_transform_to_paired_points)...")
    landmark_tx = ants.fit_transform_to_paired_points(
        moving_pts, fixed_pts, transform_type="affine", verbose=True,
    )
    ants.write_transform(landmark_tx, LANDMARK_XFM_PATH)
    print(f"  wrote {LANDMARK_XFM_PATH}")

    residuals = []
    for mv, fx in zip(moving_pts, fixed_pts):
        mapped = np.array(landmark_tx.apply_to_point(tuple(fx)))
        residuals.append(np.linalg.norm(mapped - mv))
    residuals = np.array(residuals) * 1000.0
    print(f"  landmark residuals (um): mean={residuals.mean():.1f}  max={residuals.max():.1f}")

    working_spacing_mm = (WORKING_SPACING_UM / 1000.0,) * 3
    print(f"Resampling to {WORKING_SPACING_UM} um isotropic for the fitting stage...")
    fixed_lo = ants.resample_image(fixed, working_spacing_mm, use_voxels=False, interp_type=0)
    moving_lo = ants.resample_image(moving, working_spacing_mm, use_voxels=False, interp_type=0)
    moving_mask_lo = ants.resample_image(moving_mask, working_spacing_mm, use_voxels=False, interp_type=1)
    fixed_mask_lo = ants.resample_image(fixed_mask, working_spacing_mm, use_voxels=False, interp_type=1)

    print("Stage A: Rigid-only refinement on top of the landmark affine "
          "(confirmed no-op for this data -- see module docstring)...")
    clean_prefix(os.path.join(OUT_DIR, "rigid_"))
    reg_rigid = ants.registration(
        fixed=fixed_lo, moving=moving_lo, type_of_transform="Rigid",
        initial_transform=[LANDMARK_XFM_PATH],
        outprefix=os.path.join(OUT_DIR, "rigid_"),
        mask=fixed_mask_lo, moving_mask=moving_mask_lo, mask_all_stages=True,
        aff_metric="mattes",
        verbose=True,
    )
    print(f"  rigid fwdtransforms: {reg_rigid['fwdtransforms']}")

    print("Stage B: SyNOnly deformable refinement (total_sigma=2.0)...")
    clean_prefix(OUT_XFM_PREFIX)
    reg_syn = ants.registration(
        fixed=fixed_lo, moving=moving_lo, type_of_transform="SyNOnly",
        initial_transform=reg_rigid["fwdtransforms"],
        outprefix=OUT_XFM_PREFIX,
        mask=fixed_mask_lo, moving_mask=moving_mask_lo, mask_all_stages=True,
        syn_metric="CC", syn_sampling=3,
        reg_iterations=(60, 40, 30),
        grad_step=0.25, flow_sigma=1.5, total_sigma=2.0,
        verbose=True,
    )
    print(f"  syn fwdtransforms: {reg_syn['fwdtransforms']}")

    print("Applying composed forward transform to the FULL-resolution moving volume...")
    aligned = ants.apply_transforms(
        fixed, moving, transformlist=reg_syn["fwdtransforms"], interpolator="linear",
    )
    fixed_nib_affine = nib.load(FIXED_PATH).affine
    save_nifti_ras(OUT_WARPED_FULLRES, aligned.numpy(), fixed_nib_affine)
    print(f"Wrote {OUT_WARPED_FULLRES}  shape {aligned.shape}")
    print("Done.")


if __name__ == "__main__":
    main()
