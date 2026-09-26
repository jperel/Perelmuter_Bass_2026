#!/usr/bin/env python3
"""
Stage 3 (of 6) - structural-to-template registration: template geometry fix.

Produces the canonical `template_ras.nii.gz` used as the fixed/reference volume for
every registration in this stage, by taking the (deghosted) population-average template
volume and re-saving it with a known-good NIfTI affine that is reused byte-identically
from a previously-validated reference volume, rather than trusting a freshly re-derived
or re-saved header.

This defensive pattern -- copy a trusted affine wholesale instead of recomputing one --
is used throughout this project because earlier runs hit real bugs from mm/um unit
mismatches and RAS/LPS sign-convention confusion when headers were rebuilt from scratch.
Reusing a known-good affine sidesteps that entire class of bug.

Inputs (template-construction intermediates, upstream of this 6-stage pipeline and not
included in this release -- see README.md):
  - the deghosted population-average template volume (raw voxel data only; its own
    header is not trusted)
  - a reference volume carrying the known-good affine to copy onto it

Because its inputs are not part of the released data, this script is included for
provenance/documentation of the exact geometry-fix logic that produced the shipped
template, not as a step you need to re-run: the output below is provided directly.

Output:
  data_for_upload/03_structural_to_template_registration/stage0_geometry_fixed/template_ras.nii.gz

Run order in this stage: fix_template_geometry.py -> register_structural_to_template.py
-> qc_overlay.py. (The structural-side counterparts of this file,
02F_stack_ras_cropped.nii.gz and mask_ras_cropped.nii.gz, are provided directly for the
same reason -- see README.md.)

Environment: antspy env (ants, nibabel, numpy, scipy, scikit-image, tifffile).
"""

import os

import numpy as np
import nibabel as nib

DATA_ROOT = os.environ.get("PIPELINE_DATA_ROOT") or os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data_for_upload"))

STAGE_DIR = os.path.join(DATA_ROOT, "03_structural_to_template_registration")
OUT = os.path.join(STAGE_DIR, "stage0_geometry_fixed", "template_ras.nii.gz")

# Template-construction intermediates upstream of this pipeline (not included in this
# release). Override with the environment variables below if you are re-deriving the
# template from scratch; otherwise this script does not need to be run.
SRC = os.environ.get("DEGHOSTED_TEMPLATE_PATH", "antsBTPtemplate0_fixed_cropped_deghosted.nii.gz")
REFERENCE_AFFINE_SRC = os.environ.get("REFERENCE_AFFINE_TEMPLATE_PATH", "template_ras_reference.nii.gz")


def main():
    img = nib.load(SRC)
    reference_affine = nib.load(REFERENCE_AFFINE_SRC).affine
    print(f"source shape {img.shape}")
    print(f"reusing exact affine from {REFERENCE_AFFINE_SRC}:\n{reference_affine}")

    data = np.asanyarray(img.dataobj)
    if data.ndim == 4 and data.shape[3] == 1:
        data = data[:, :, :, 0]

    out_img = nib.Nifti1Image(data.astype("float32"), reference_affine)
    out_img.header.set_data_dtype("float32")
    out_img.header.set_xyzt_units(xyz="micron")
    out_img.set_sform(reference_affine, code=1)
    out_img.set_qform(reference_affine, code=1)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    nib.save(out_img, OUT)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
