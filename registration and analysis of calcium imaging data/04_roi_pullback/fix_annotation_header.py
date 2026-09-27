#!/usr/bin/env python3
"""
fix_annotation_header.py
===========================
Stage 4 of 5 (roi_pullback), step 0.

The manually-drawn annotation volume (``annotations.nii.gz``, this stage's own input
data -- see the manuscript's Methods for how it was produced) is
voxel-for-voxel identical in shape to the template grid that stage 03's
structural<->template transform was actually fit against
(``template_ras.nii.gz``), but it was saved by ITK-SNAP with a
different (though internally self-consistent) header/affine convention, including a
nonzero origin. Rather than trust the annotation's own saved header, this script copies
the label array directly onto the known-good affine from the template file the
transform chain actually operates against. This is the same defensive pattern used
elsewhere in this pipeline for the same reason: re-derived or re-saved affines have been
a real source of bugs, so known-good affines are reused verbatim wherever possible
instead of re-derived.

Reads (relative to DATA_ROOT):
    template_ras.nii.gz (source of the known-good affine)
    annotations.nii.gz (label array to re-header)

Writes:
    annotations_fixed_header.nii.gz

Run before ``pull_roi_labels.py``, which consumes this script's output.

Note on data layout: `data_for_upload/` is a flat folder (Zenodo does not
preserve directory structure on upload); see the top-level README.

Requires the ``antspy`` environment: nibabel, numpy.
"""

import os

import nibabel as nib
import numpy as np

DATA_ROOT = os.environ.get("PIPELINE_DATA_ROOT") or os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data_for_upload"))

TEMPLATE_PATH = os.path.join(DATA_ROOT, "template_ras.nii.gz")
ANNOTATIONS_PATH = os.path.join(DATA_ROOT, "annotations.nii.gz")
OUT_PATH = os.path.join(DATA_ROOT, "annotations_fixed_header.nii.gz")


def main():
    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)

    known_good = nib.load(TEMPLATE_PATH)
    known_good_data = np.asarray(known_good.dataobj)

    annotations_img = nib.load(ANNOTATIONS_PATH)
    annotations_data = np.asarray(annotations_img.dataobj)

    # Re-verify voxel-grid identity every run rather than assuming it -- this is the
    # premise the whole re-header operation depends on.
    assert annotations_data.shape == known_good_data.shape, (
        annotations_data.shape, known_good_data.shape
    )

    affine = known_good.affine
    out_img = nib.Nifti1Image(annotations_data.astype(np.uint16), affine)
    out_img.header.set_data_dtype("uint16")
    out_img.header.set_xyzt_units(xyz="micron")
    out_img.set_sform(affine, code=1)
    out_img.set_qform(affine, code=1)
    nib.save(out_img, OUT_PATH)
    print(f"Wrote {OUT_PATH}")

    # Confirm the array is untouched -- only the header changed.
    check = np.asarray(nib.load(OUT_PATH).dataobj)
    assert np.array_equal(check, annotations_data)
    print("Confirmed: label array unchanged, only header/affine replaced.")


if __name__ == "__main__":
    main()
