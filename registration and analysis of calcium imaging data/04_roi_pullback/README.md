# 04 - ROI pullback

Stage 4 of 5. This is the glue stage that composes stages 02-03's registration
chain with a hand-traced region-annotation volume into per-cell brain-region
assignments: it pulls `annotations.nii.gz` (this stage's own input data -- a
finished region-label volume manually traced in ITK-SNAP; tracing/polishing
methodology is described in the manuscript's Methods) backward through the
full registration chain (stage 03's structural<->template transform, then
stage 02's functional<->structural transform) onto each native functional
plane's own native pixel grid, then assigns each suite2p-detected cell (from
stage 01/05's suite2p output) a brain region by majority vote over its pixel
footprint.

## Data layout

All paths below are relative to `DATA_ROOT` (defaults to `data_for_upload/`
next to the code repository; override with `PIPELINE_DATA_ROOT`).
**`data_for_upload/` is a flat folder** -- Zenodo does not preserve directory
structure on upload, so every file lives directly under `DATA_ROOT`, with
three exceptions: `mean_images/`, `roi_labels/`, and `suite2p_output/` are
subfolders created locally by unzipping the three `.zip` archives in the data
package.

## Run order

1. `fix_annotation_header.py`
   - Reads: `template_ras.nii.gz` (known-good affine), `annotations.nii.gz` (label array)
   - Writes: `annotations_fixed_header.nii.gz`
   - Copies the annotation's label array onto the known-good template affine rather
     than trusting ITK-SNAP's own saved header (see "Geometry notes" below).

2. `pull_roi_labels.py`
   - Reads: `mean_images/meanImg_z<depth>.tif` (36 planes),
     `annotations_fixed_header.nii.gz`, `fwd_final_affine.mat`,
     `functional_ministack.nii.gz`, `stageB_plane_corrections.csv`,
     `fwd_0GenericAffine.mat`, `fwd_1InverseWarp.nii.gz`
   - Writes: `roi_labels/roi_labels_z<depth>.tif` (32 planes)
   - Imports `compute_seed_origin` from stage 02's `register_ministack_to_structural.py`
     (sibling stage in the code repository; see "Cross-stage import" below).

3. `assign_cells_to_regions.py`
   - Reads: `roi_labels/roi_labels_z<depth>.tif`, `labels.txt`,
     `suite2p_output/02F_2min_<depth>/suite2p/plane0/{stat.npy,iscell.npy}`
   - Writes: `cell_region_assignments.csv`

4. `export_per_plane_qc_with_suite2p_rois.py` (QC, optional but recommended)
   - Reads: `mean_images/meanImg_z<depth>.tif`, `roi_labels/roi_labels_z<depth>.tif`,
     `labels.txt`, `suite2p_output/.../stat.npy`, `cell_region_assignments.csv`
   - Writes: `qc_roi_pullback/roi_qc_with_cells_z<depth>.png`
   - Mean image + pulled-back region boundaries + suite2p's own `iscell==1` ROIs
     painted on top (green if assigned to a real region, red if landing on "Clear
     Label" or "ventricle"). Use this to sanity-check region-assignment quality for a
     given plane before trusting its numbers, and to visually check the known caveat
     below.

Stages 02 and 03 must have already produced their outputs before running this
stage. `annotations.nii.gz` and `labels.txt` (the ITK-SNAP region ID -> name/color
lookup) are included directly in this stage's own data -- they are not produced
by any script in this release; see the manuscript's Methods for how they were made.

## Environment

`antspy` env: `ants` (ANTsPy), `nibabel`, `numpy`, `scipy`, `tifffile`, plus
`matplotlib` and `pandas` for the QC script.

## Why this can't be a naive flat crop (tilted-slab geometry)

A native functional plane is physically flat (it's a single optical section), but it
does **not** map to a single flat z-slice in structural/template space. The stage-02
functional<->structural affine (`fwd_final_affine.mat`) contains a genuine ~19-degree
rotation mixing the Y and Z axes, so each native plane maps to a *tilted* plane in
reference space, spreading across many structural z-indices from one edge of the frame
to the other. Picking "the nearest matching z-slice" and doing a flat 2D crop would
silently assume every pixel in the native frame sits at the same reference depth, which
is wrong -- pixels near the edges of the frame are the most wrong.

The method used instead (in `pull_roi_labels.py`): embed each native plane as a real,
1-voxel-thick 3D image at its own true (x, y, z) position in the shared coordinate
frame, then pull the labels backward through the full composed inverse transform chain
(template -> stage 03 inverse -> structural -> stage 02 inverse) in a single
`ants.apply_transforms` call, with `interpolator="genericLabel"` (nearest-neighbor,
since these are integer region IDs that must never be linearly interpolated). This
evaluates the true tilted 3D geometry per native pixel rather than assuming one shared
depth for the whole frame.

## Geometry notes

- The `annotations.nii.gz` volume is voxel-for-voxel identical in shape to the template
  grid the stage-03 transform was actually fit against, but ITK-SNAP saved it with a
  different (though internally consistent) header/affine convention, including a
  nonzero origin. `fix_annotation_header.py` copies the label array onto the known-good
  template affine rather than trusting the annotation's own saved header -- the same
  defensive pattern used in stage 03 for the same reason (this project hit real
  header bugs from trusting re-derived or re-saved affines).
- `ants.image_read()`'s RAS->LPS conversion negates the X/Y origin, not just the
  direction matrix; `ants.from_numpy()` does not, so `pull_roi_labels.py` negates X/Y by
  hand when constructing each plane's slab image to match.

## Cross-stage import

`pull_roi_labels.py` imports `compute_seed_origin` from stage 02's
`register_ministack_to_structural.py` (not transcribed) so that the coordinate-frame
origin used here is guaranteed identical to the one stage 02's registration was fit
against. The import path is resolved relative to this script's own location:

```python
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "02_functional_to_structural_registration")))
from register_ministack_to_structural import compute_seed_origin
```

## Included planes

32 of the 36 native functional planes are carried through this stage: z = 2480-2790 in
10 um steps. The 4 shallowest planes (z = 2440, 2450, 2460, 2470) are excluded -- they
fall in a stage-02 depth range with a documented residual correction that was never
fully resolved into an invertible transform, and they mostly overlap planes already
excluded from the suite2p analysis (stage 05) for insufficient signal.

## `cell_region_assignments.csv` schema

| column        | meaning                                                                 |
|---------------|--------------------------------------------------------------------------|
| `z`           | plane depth (um)                                                        |
| `cell_id`     | suite2p's own per-plane cell index; `(z, cell_id)` together are the unique key |
| `iscell`      | suite2p's `iscell.npy` classification (0/1)                             |
| `iscell_prob` | suite2p's classifier probability                                        |
| `n_pixels`    | number of valid pixels in the cell's footprint used for the region vote |
| `region_id`   | majority-vote region ID from the pulled-back label image; 0 = "Clear Label" (not a real region) |
| `region_name` | region name looked up from `labels.txt`                   |
| `purity`      | fraction of the cell's pixels agreeing with `region_id`                 |

Region assignment is by majority vote over each cell's full pixel footprint (not just
its centroid), which is more robust for cells straddling a region boundary; `purity`
quantifies how clean that vote was for a given cell.

## Known caveat: top-middle rim misalignment

There is a confirmed, unresolved local misalignment specifically in the top-middle
boundary/rim region of the atlas. It is worst at the deepest included plane (z = 2790)
and tapers off at shallower depth. Cells physically located in that rim area may show
inflated "Clear Label" assignments, or occasional wrong-neighbor region assignment.
This is a localized issue, not a general dataset-wide problem, and can be checked
visually per-depth via `export_per_plane_qc_with_suite2p_rois.py`'s output for the
affected planes before trusting region counts drawn from that area.
