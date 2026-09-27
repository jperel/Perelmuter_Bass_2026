# 03 - Structural-to-template registration

Stage 3 of 5 in the pipeline. Registers the RAS-converted structural stack (output of
stage 02, `02_functional_to_structural_registration`) to a population-average reference
template. A hand-drawn brain-region annotation volume (`annotations.nii.gz`, shipped as
stage 4's input data; tracing/polishing methodology is in the manuscript's Methods) was
drawn on this same template. Registration here is landmark-based affine followed by SyN
deformable registration. The resulting forward/inverse transform chain is what stage 4
(`04_roi_pullback`) uses to pull template region labels back onto native
functional/structural space.

The original project evaluated several parameter variants of the deformable-registration
step. Only the accepted, final variant is included here; nothing in this codebase or
data should be read as leaving other variants out by omission or oversight.

## Run order

```
fix_template_geometry.py           # builds the fixed-image template with a trusted affine
register_structural_to_template.py # landmark affine + SyN registration
qc_overlay.py                      # visual QC overlays (FIJI/ImageJ)
```

## Environment

`antspy` environment: `ants`, `nibabel`, `numpy`, `scipy`, `scikit-image`, `tifffile`.

## Inputs / outputs

All paths below are relative to `DATA_ROOT` (defaults to `data_for_upload/`
next to the code repository; override with the `PIPELINE_DATA_ROOT`
environment variable). **`data_for_upload/` is a flat folder** -- Zenodo does
not preserve directory structure on upload, so every file lives directly
under `DATA_ROOT`, addressed by filename only.

- `landmarks.csv` -- 19 manually-placed landmark pairs (moving = structural stack, fixed
  = template), by voxel index in each volume.
- `template_ras.nii.gz` -- the fixed/reference template volume, with a trusted,
  byte-identical NIfTI affine (see `fix_template_geometry.py`).
- `02F_stack_ras_cropped.nii.gz` -- the moving structural stack, RAS-converted and
  cropped.
- `mask_ras_cropped.nii.gz` -- brain mask for the structural stack, same grid as above.
- `landmark_affine.mat` -- the 12-DOF affine fit directly from the landmark pairs.
- `fwd_0GenericAffine.mat`, `fwd_1Warp.nii.gz`, `fwd_1InverseWarp.nii.gz` -- the final
  forward-transform chain (affine + SyN warp) and its inverse warp, as produced by
  `ants.registration`.
- `qc_structural_to_template/*.tif` -- QC overlay stacks (see below), produced by
  `qc_overlay.py` (local output, not part of the Zenodo data package).


## QC outputs

`qc_overlay.py` writes two 2-channel TIFF stacks to `qc_structural_to_template/`:

- `check_template_vs_warped_2ch_stack.tif` -- template + structural warped into template
  space (forward chain), on the template's grid.
- `check_structural_vs_warped_template_2ch_stack.tif` -- structural stack in its own
  native space + template warped into structural space (inverse chain), on the
  structural stack's grid. This is the view relevant to stage 4's ROI pull-back, which
  uses the same inverse chain.

