# Stage 2 of 6: functional-to-structural registration

Registers the 36 native functional (2P GCaMP) z-planes, treated as one
coherent 3D mini-stack, into the native space of a separate, higher-SNR
structural stack acquired from the same fish. Runs after Stage 1 (needs its
`mean_images/` output) and its output feeds Stage 5 (ROI pullback).

Internally referred to as "Stage B" in the scripts and comments below (the
functional-stack-to-template registration in Stage 3 is "Stage A" in the same
naming scheme).

## Data layout

Reads from and writes to `data_for_upload/` via `DATA_ROOT`, resolved as:

```python
DATA_ROOT = os.environ.get("PIPELINE_DATA_ROOT") or os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data_for_upload"))
```

Relevant subfolders:

- `data_for_upload/01_suite2p_preprocessing/mean_images/` (input, from Stage 1)
- `data_for_upload/02_functional_to_structural_registration/stage1_ministack/` (this stage)
- `data_for_upload/02_functional_to_structural_registration/stage2_registration/` (this stage)
- `data_for_upload/02_functional_to_structural_registration/stage2_reference/` (this stage; fixed optimizer seed)
- `data_for_upload/03_structural_to_template_registration/stage0_geometry_fixed/` (input, from Stage 3)

## Run order

1. **verify_pixel_size.py** -- Step 0, standalone. Re-derives the functional
   camera's true physical pixel size from OME-XML metadata in the raw `.oir`
   files, rather than trusting a comment. Confirms the 1.40625 um/px (XY) and
   10.0 um (Z step) values used everywhere below. Needs the raw acquisition
   files and a local Bio-Formats JVM/jar (not part of `data_for_upload`;
   configured via `FUNCTIONAL_OIR_DIR`, `BIOFORMATS_JAR_PATH`,
   `BIOFORMATS_JVM_PATH`). Writes nothing to disk.

2. **correct_interplane_drift.py** -- Step 1a. Each of the 36 planes was
   acquired as an independent imaging session, so there is no cross-plane XY
   registration a priori (suite2p only corrects motion within a single
   plane's own movie). Runs sequential pairwise 2D rigid registration between
   adjacent planes (10um apart), anchored at the middle plane and propagated
   outward in both directions, and applies the 90-degree CCW rotation needed
   to match the structural stack's orientation.
   Reads: `01_suite2p_preprocessing/mean_images/meanImg_z*.tif`.
   Writes: `stage1_ministack/aligned_planes/`, `stage1_ministack/aligned_planes_valid_mask/`,
   `stage1_ministack/interplane_drift_log.txt`.

3. **build_functional_ministack.py** -- Step 1. Stacks the 36 mean images
   into one 3D NIfTI (spacing 1.40625 x 1.40625 x 10.0 um; array axis order
   X,Y,Z; identity direction cosines; origin (0,0,0); `xyzt_units="micron"`,
   i.e. pixdim values are literal microns with no mm scaling). Set
   `USE_ALIGNED_PLANES=1` to build from Step 1a's drift-corrected planes
   instead of the raw mean images directly -- the shipped
   `functional_ministack.nii.gz` was built this way (run step 2 first).
   Writes: `stage1_ministack/functional_ministack.nii.gz`, `stage1_ministack/depths_um.txt`.

4. **sanity_check_step1.py** -- visual check before any automated fit: plots
   a functional plane against candidate structural slices at matched physical
   scale, to catch a flip/rotation/scale error by eye.
   Writes: `stage1_ministack/sanity_check.png`.

5. **build_ministack_mask.py** -- stacks the per-plane valid-data masks from
   step 2 into a mask NIfTI with the same geometry as the ministack. Used as
   `moving_mask` in the ANTs registration so drift-correction's ragged empty
   borders don't corrupt the correlation metric.
   Writes: `stage1_ministack/functional_ministack_mask.nii.gz` (intermediate,
   not distributed -- cheaply regenerated from step 2's output).

6. **register_ministack_to_structural.py** -- Step 2. Registers the 36-plane
   ministack into the structural stack's native space as ONE 3D
   volume-to-volume fit (Rigid seeding an Affine), not 36 independent
   per-plane fits. Seeded from a fixed reference transform
   (`stage2_reference/reference_affine.mat`, the near-identity result of
   registering the *uncorrected* ministack, which converges reliably on its
   own) because a blind search on the drift-corrected data repeatedly
   converged to a wrong rotation. Defines `compute_seed_origin`, imported by
   step 7 below and by Stage 5's `pull_roi_labels.py` -- its name and
   signature must not change.
   Writes: `stage2_registration/fwd_final_affine.mat` (the accepted result),
   plus intermediate seeded volumes and ANTs transform files, and
   `stage2_registration/ministack_in_structural_space.nii.gz` (full-resolution
   warped ministack, for visual QC only).

   **Important caveat for anyone consuming `fwd_final_affine.mat` directly**:
   the fitted affine includes a genuine rotation of roughly 19 degrees that
   mixes the Y and Z axes. A native functional plane therefore maps to a
   *tilted* plane in structural space, not a flat z-slice. Stage 5 handles
   this correctly with a "tilted slab" technique (each native plane embedded
   as a true 1-voxel-thick 3D image at its own position, pulled through the
   full inverse transform chain in one `apply_transforms` call) -- a naive
   flat-z-slice lookup will be wrong, especially near frame edges.

7. **build_stageB_plane_corrections.py** -- builds the accepted per-plane 2D
   correction table for the residual functional-vs-structural misalignment
   that varies by depth (likely a small residual rotation-axis error in the
   Stage 2 affine). Measures the residual via `skimage.feature.match_template`
   (a large central template matched by FFT-based normalized
   cross-correlation against the entire frame) -- chosen after an earlier
   small-block phase-correlation approach was found to lock onto the wrong
   repeat of suite2p's repeating cell-blob pattern (confirmed at one depth:
   measured 6.9px, true offset approximately 54px). Imports
   `compute_seed_origin` from step 6 (same-folder import).
   Writes: `stage2_registration/stageB_plane_corrections.csv` -- **this is the
   file actually consumed downstream, by Stage 5's `pull_roi_labels.py`.**

8. **diag_stageB_block_displacement.py** -- validation/QC companion to step
   7: dense measurement of the same residual across all 32 included planes
   via local block-matching (reliable here because both channels are the
   same imaging modality). Writes `stage2_registration/diag_stageB_block_results.npy`
   (diagnostic only, not consumed by any downstream step).

## Environment

`antspy` env: `ants`, `nibabel`, `numpy`, `scipy`, `scikit-image`, `tifffile`
(plus `matplotlib` for the two diagnostic/sanity-check scripts).
`verify_pixel_size.py` instead runs in the same env as Stage 1's
`convert_oir_to_tiff.py` (`jpype` against a local Bio-Formats jar/JVM).

## Notes

- 32 of the 36 native planes are carried through to Stage 5; the shallowest 4
  (z = 2440-2470) are excluded there due to an unresolved shallow-depth
  registration issue in that z-range, not handled by any script in this
  stage.
- The excluded planes and the Stage B correction table's exclusion set
  (`EXCLUDED_IDX = {0, 1, 2, 3}`) are consistent across this stage and Stage 5.
