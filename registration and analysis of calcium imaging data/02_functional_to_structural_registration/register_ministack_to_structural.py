#!/usr/bin/env python3
"""
register_ministack_to_structural.py
======================================
Stage 2 of 6 (functional-to-structural registration) -- Step 2.

Registers the 36-plane functional mini-stack
(stage1_ministack/functional_ministack.nii.gz) into the structural stack's
native space (stage 03's stage0_geometry_fixed/02F_stack_ras_cropped.nii.gz)
as ONE 3D volume-to-volume fit (Rigid, then Affine seeded from the Rigid
result) -- not 36 independent per-plane fits, since the mini-stack's own
internal geometry (10um plane spacing, fixed by construction) is already
known and trustworthy.

Before running the optimizer, the mini-stack's origin is shifted (a pure
header edit, no resampling) so its physical bounding box is roughly centered
on the structural stack's, and its mid-plane lands near native structural
z-idx 44 (the best visual shape match found in sanity_check_step1.py) -- so
the search starts from a sensible position instead of blind.

Both images are resampled to a common ~2um isotropic working resolution
before the search: full native-resolution registration on a volume this size
is impractically slow for ANTs' single-threaded optimizer. The fitted affine
lives in continuous physical space, so it is resolution-independent and is
applied to the full-resolution data afterward.

Uses functional_ministack_mask.nii.gz (built by build_ministack_mask.py) as
moving_mask with mask_all_stages=True: after inter-plane drift correction,
each plane has a differently shifted empty (zero) border, and without a mask
the correlation metric's sampling would land mostly in those zero regions.

Masking alone was not sufficient -- with a blind Rigid search, the optimizer
repeatedly converged on a wildly wrong (but numerically valid) rotation,
since the mask's irregular shape combined with the GC metric's sparse
"regular" sampling produces a poorly behaved objective landscape at coarse
pyramid levels. The fit is instead seeded from stage2_reference/reference_affine.mat,
a known-good near-identity result obtained by registering the UNCORRECTED
mini-stack (which converges reliably on its own, confirmed by QC overlay), so
the optimizer starts in the right neighborhood instead of searching blind.

The registration schedule for both the Rigid and Affine stages is
SINGLE-LEVEL, fine-resolution-only (aff_iterations/shrink_factors/
smoothing_sigmas all length-1): the default multi-resolution schedule's
coarsest level (6x shrink, sigma=3 smoothing) blurs away the real punctate
signal while leaving the drift-corrected planes' ragged empty borders as the
dominant coherent shape, which is exactly where the optimizer kept abandoning
the good seed for a trivial background-matches-background local optimum.
Skipping straight to a fine-only local refinement gives stable, sensible
(determinant near 1.0) results, unlike any multi-resolution schedule tried,
and keeps the optimizer in the seed's neighborhood.

The GC (gradient correlation) metric is used rather than mutual information:
both volumes are same-modality 2P GCaMP data from the same fish, unlike the
confocal-vs-2P structural-to-template registration in stage 03, where MI is
needed for cross-modality robustness.

Note for downstream consumers of fwd_final_affine.mat: the fitted affine
includes a genuine rotation of roughly 19 degrees that mixes the Y and Z
axes. A native functional plane therefore maps to a TILTED plane in
structural space, not a flat z-slice -- a naive flat-z-slice lookup is wrong,
especially near frame edges. Stage 05 (roi_pullback/pull_roi_labels.py)
handles this correctly via a "slab" technique: each native plane is embedded
as a true 1-voxel-thick 3D image at its own position and pulled through the
full inverse transform chain in one apply_transforms call.

compute_seed_origin() is imported directly (not transcribed) by
build_stageB_plane_corrections.py in this same folder, and by stage 05's
pull_roi_labels.py via sys.path -- keep its name and signature stable.

Reads:
  - data_for_upload/02_functional_to_structural_registration/stage1_ministack/
    functional_ministack.nii.gz, functional_ministack_mask.nii.gz
  - data_for_upload/03_structural_to_template_registration/stage0_geometry_fixed/
    02F_stack_ras_cropped.nii.gz
  - data_for_upload/02_functional_to_structural_registration/stage2_reference/
    reference_affine.mat (fixed optimizer seed; see above)

Writes (to stage2_registration/):
  - fwd_final_affine.mat (the accepted result; consumed by stage 05)
  - functional_ministack_seeded.nii.gz, functional_ministack_mask_seeded.nii.gz
    (intermediate, header-shifted copies of the inputs)
  - rigid_*, fwd_* ANTs transform/log files from the two registration stages
  - ministack_in_structural_space.nii.gz (full-resolution warped mini-stack,
    for visual QC only -- not required by any downstream step)

Run order: after build_ministack_mask.py, before
build_stageB_plane_corrections.py.

Environment: antspy env (ants, nibabel, numpy).
"""

import os
import shutil

import ants
import nibabel as nib
import numpy as np

DATA_ROOT = os.environ.get("PIPELINE_DATA_ROOT") or os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data_for_upload"))

STAGE2_DIR_NAME = "02_functional_to_structural_registration"
STAGE1_DIR = os.path.join(DATA_ROOT, STAGE2_DIR_NAME, "stage1_ministack")
MINISTACK_PATH = os.path.join(STAGE1_DIR, "functional_ministack.nii.gz")
MASK_PATH = os.path.join(STAGE1_DIR, "functional_ministack_mask.nii.gz")
STRUCTURAL_PATH = os.path.join(DATA_ROOT, "03_structural_to_template_registration",
                                "stage0_geometry_fixed", "02F_stack_ras_cropped.nii.gz")
OUT_DIR = os.path.join(DATA_ROOT, STAGE2_DIR_NAME, "stage2_registration")
SEEDED_PATH = os.path.join(OUT_DIR, "functional_ministack_seeded.nii.gz")
SEEDED_MASK_PATH = os.path.join(OUT_DIR, "functional_ministack_mask_seeded.nii.gz")
OUTPREFIX = os.path.join(OUT_DIR, "fwd_")
REFERENCE_TRANSFORM = os.path.join(DATA_ROOT, STAGE2_DIR_NAME, "stage2_reference", "reference_affine.mat")

SEED_STRUCTURAL_Z_IDX = 44  # best visual shape match from sanity_check_step1.py
WORKING_SPACING_UM = 2.0    # common working resolution for the search stage

SEED = 20260826


def compute_seed_origin(mini_data_shape, mini_spacing):
    mini_extent = np.array(mini_data_shape) * mini_spacing

    struct = nib.load(STRUCTURAL_PATH)
    struct_shape = struct.shape[:3]
    struct_spacing = np.array(struct.header.get_zooms()[:3])
    struct_extent = np.array(struct_shape) * struct_spacing

    # center XY on the structural stack's XY bounding box; center Z on the
    # physical position of native structural z-idx SEED_STRUCTURAL_Z_IDX
    seed_origin = np.zeros(3)
    seed_origin[0] = (struct_extent[0] - mini_extent[0]) / 2.0
    seed_origin[1] = (struct_extent[1] - mini_extent[1]) / 2.0
    z_center_um = SEED_STRUCTURAL_Z_IDX * struct_spacing[2]
    seed_origin[2] = z_center_um - mini_extent[2] / 2.0

    print(f"Mini-stack extent (um): {mini_extent}")
    print(f"Structural extent (um): {struct_extent}")
    print(f"Seed origin (um): {seed_origin}")
    return seed_origin


def save_seeded(data, spacing, seed_origin, out_path, dtype):
    affine = np.diag(list(spacing) + [1.0])
    affine[:3, 3] = seed_origin
    img = nib.Nifti1Image(data, affine)
    img.header.set_data_dtype(dtype)
    img.header.set_xyzt_units(xyz="micron")
    img.set_sform(affine, code=1)
    img.set_qform(affine, code=1)
    nib.save(img, out_path)
    print(f"Wrote {out_path}")


def build_seeded_ministack():
    mini = nib.load(MINISTACK_PATH)
    mini_data = np.asarray(mini.dataobj)
    mini_spacing = np.array(mini.header.get_zooms()[:3])
    seed_origin = compute_seed_origin(mini_data.shape, mini_spacing)
    save_seeded(mini_data, mini_spacing, seed_origin, SEEDED_PATH, "float32")

    mask = nib.load(MASK_PATH)
    mask_data = np.asarray(mask.dataobj)
    assert mask_data.shape == mini_data.shape
    save_seeded(mask_data, mini_spacing, seed_origin, SEEDED_MASK_PATH, "uint8")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    build_seeded_ministack()

    ants.config.set_ants_deterministic(True, seed_value=SEED)

    fixed_full = ants.image_read(STRUCTURAL_PATH)
    moving_full = ants.image_read(SEEDED_PATH)
    moving_mask_full = ants.image_read(SEEDED_MASK_PATH)
    print(f"fixed_full spacing {fixed_full.spacing} shape {fixed_full.shape}")
    print(f"moving_full spacing {moving_full.spacing} shape {moving_full.shape}")

    working_spacing_mm = WORKING_SPACING_UM / 1000.0  # ants/ITK internal units, matches this file's mm-equivalent scale
    target_spacing = (working_spacing_mm, working_spacing_mm, working_spacing_mm)
    fixed = ants.resample_image(fixed_full, target_spacing, use_voxels=False, interp_type=0)
    moving = ants.resample_image(moving_full, target_spacing, use_voxels=False, interp_type=0)
    # nearest-neighbor for the mask (interp_type=1) -- must stay strictly binary
    moving_mask = ants.resample_image(moving_mask_full, target_spacing, use_voxels=False, interp_type=1)
    print(f"fixed (working res) spacing {fixed.spacing} shape {fixed.shape}")
    print(f"moving (working res) spacing {moving.spacing} shape {moving.shape}")
    print(f"moving_mask (working res) nonzero frac: {(moving_mask.numpy() > 0).mean():.4f}")

    for f in os.listdir(OUT_DIR):
        if f.startswith("fwd_") and "[0-9]" not in f:
            fp = os.path.join(OUT_DIR, f)
            if os.path.isfile(fp) and f != os.path.basename(SEEDED_PATH):
                os.remove(fp)

    assert os.path.exists(REFERENCE_TRANSFORM), (
        f"missing {REFERENCE_TRANSFORM} -- this is the fixed optimizer seed obtained "
        "once from an independent (blind Rigid) registration of the uncorrected "
        "mini-stack; see the module docstring."
    )
    single_level_kwargs = dict(aff_iterations=(50,), aff_shrink_factors=(1,), aff_smoothing_sigmas=(0,))

    print(f"\n=== Rigid (single-level, seeded from {REFERENCE_TRANSFORM}) ===")
    rigid_reg = ants.registration(
        fixed=fixed, moving=moving,
        type_of_transform="Rigid",
        aff_metric="GC",
        moving_mask=moving_mask,
        mask_all_stages=True,
        initial_transform=REFERENCE_TRANSFORM,
        outprefix=os.path.join(OUT_DIR, "rigid_"),
        **single_level_kwargs,
    )
    print("Rigid transforms:", rigid_reg["fwdtransforms"])

    print("\n=== Affine (single-level, seeded from Rigid) ===")
    affine_reg = ants.registration(
        fixed=fixed, moving=moving,
        type_of_transform="Affine",
        aff_metric="GC",
        moving_mask=moving_mask,
        mask_all_stages=True,
        initial_transform=rigid_reg["fwdtransforms"][0],
        outprefix=OUTPREFIX,
        **single_level_kwargs,
    )
    print("Affine transforms:", affine_reg["fwdtransforms"])

    shutil.copy(affine_reg["fwdtransforms"][0], os.path.join(OUT_DIR, "fwd_final_affine.mat"))

    # the found affine lives in continuous physical space -- apply it directly to the
    # FULL-resolution seeded mini-stack against the FULL-resolution structural stack,
    # not the working-resolution copies used for the search itself
    print("\n=== Applying found transform to full-resolution data ===")
    warped_full = ants.apply_transforms(
        fixed=fixed_full, moving=moving_full,
        transformlist=[os.path.join(OUT_DIR, "fwd_final_affine.mat")],
        interpolator="linear",
    )
    ants.image_write(warped_full, os.path.join(OUT_DIR, "ministack_in_structural_space.nii.gz"))
    print(f"Wrote warped mini-stack to {os.path.join(OUT_DIR, 'ministack_in_structural_space.nii.gz')}")
    print("Done.")


if __name__ == "__main__":
    main()
