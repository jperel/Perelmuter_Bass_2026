# 04 — Annotation

Stage 4 of 6 in the Danionella telencephalon registration / calcium-imaging
pipeline. This stage produces the brain-region label volume: a manual
tracing in ITK-SNAP, polished by two computational steps (gap fill, then
anisotropic anterior-posterior smoothing) with a required manual correction
step done by hand in ITK-SNAP in between them.

Input: the template produced by stage 3
(`03_structural_to_template_registration`). Output: the finished annotation
volume consumed by stage 5 (`05_roi_pullback`) to pull region labels back
onto native functional data.

## Pipeline order

```
manual tracing in ITK-SNAP, section-by-section in the coronal plane
        |
        v
Annotations_V3.nii.gz                    (original tracing, provenance anchor)
        |
        |  ITK-SNAP "Smooth Labels", isotropic ~2x2x2 voxel SD (manual, in ITK-SNAP)
        v
Annotations_V3b.nii.gz                   (input to step 1)
        |
        |  STEP 1 -- gap fill                 polish_annotations.py --no-smooth
        v
   (gap-filled, intermediate -- not shipped separately)
        |
        |  MANUAL CORRECTION in ITK-SNAP -- NOT automated, see below
        v
Annotations_V3b_gapfilled2.nii.gz        (step 1 output + manual correction; input to step 2)
        |
        |  STEP 2 -- anisotropic A-P smoothing   polish_annotations.py --no-fill
        v
Annotations_V3b_final.nii.gz             (finished annotation)
```

Every step preserves `Annotations_V3`'s affine/header byte-for-byte, keeps
`uint16`, and does no resampling or reorientation. Input files are never
modified in place.

### Step 1 — gap fill

Hand-painted ROIs don't tile perfectly: a thin unlabelled "grout" network
(label 0) threads between every region and leaks out to the true exterior
through gaps in the pial surface, so a plain connected-component test cannot
tell an interior seam from outside-brain.

- Erode the label-0 mask by 3 voxels (`ERODE_K`) — thin inter-region channels
  pinch off while the outside world, continuous with the brain surface,
  survives.
- Flood from the array border through the eroded mask, dilate back -> that is
  "outside".
- Of the remaining interior label-0, fill a voxel only if it is genuinely
  *between annotations*: >= 2 distinct region labels within 4 voxels
  (`HAS2_RADIUS`), or it is fully enclosed by one region.
- Each filled voxel takes the nearest region by physical (anisotropic)
  distance (`scipy.ndimage.distance_transform_edt` with
  `sampling = SPACING_UM`).

### Manual correction (between the two steps — not automated)

`Annotations_V3b_gapfilled2.nii.gz` contains a hand edit made in ITK-SNAP that
`polish_annotations.py` does not and cannot reproduce, because it is
anatomical judgment, not an algorithm:

- ~36,830 voxels, one contiguous posterior region that had been traced as
  diencephalon / POA / Vi, were reclassified to Vp.
- ~2,000 voxels of small boundary tweaks at the Vv/Vs and ENd/Vl borders.

**If the annotation is ever re-traced from scratch, this correction must be
redone by hand in ITK-SNAP between step 1 and step 2; it is not automated.**
`Annotations_V3b_gapfilled2.nii.gz` is the only shipped file in which this
correction is recorded.

### Step 2 — anisotropic A-P smoothing

Tracing was done per coronal section, so in-plane boundaries are reliable but
ripple along the anterior-posterior axis between sections. This step runs
the same algorithm as ITK-SNAP's own "Smooth Labels" (per-label Gaussian on
the binary indicator, then per-voxel argmax over all labels including
background), but made directional: standard deviation 24 voxels (~11.7 um)
along the A-P axis, ~0.6 voxels (effectively off) in-plane. Only voxels
within `SIGMA_AP + 3` of a targeted label can change, so the coronal (traced)
plane is left visually untouched.

Some labels are deliberately excluded from smoothing: Vp, Vi, ENv, Hab, OB,
and ventricle. Their borders with the smoothed labels may shift slightly as
a result.

Running `polish_annotations.py` with both steps enabled in one pass on
`Annotations_V3b.nii.gz` reproduces `Annotations_V3b_final.nii.gz` except for
the manual-correction voxels above and their downstream propagation through
smoothing. To exactly reproduce the manuscript's final annotation, the two
steps must be run separately with the manual correction in between (see
"How to run" below).

## The label volume

| | |
|---|---|
| Grid | 1093 x 1513 x 231, `uint16` |
| Array axes | axis 0 = R-L (0.4874 um), axis 1 = **A-P** (0.4874 um), axis 2 = S-I / dorsoventral (2.0 um) |
| Intensity reference | the template this annotation was drawn on, shipped with stage 3 as `template_ras.nii.gz` (byte-identical aside from header) |
| Labels | 0 = Clear/background, 1-20 = named regions (`labels.txt`) |

### Region labels (`labels.txt`)

```
 1  Dm-r                          11  Vp
 2  Dl-r                          12  Vi
 3  Dl-c                          13  V - unspecified
 4  Dc                            14  ENd
 5  Dp                            15  ENv
 6  Dm-c                          16  POA
 7  Vl                            17  Hab
 8  Vd                            18  OB
 9  Vv                            19  ventricle
10  Vs                            20  Diencephalon - unspecified
```
(label 0 = Clear / background, not part of the 20)

## Units gotcha

The template and every file derived from it store voxel spacing in
**millimeters** in the NIfTI header (`get_xyzt_units() == ('mm', ...)`); the
raw `pixdim` value `0.0004874` is **0.4874 micrometers**, not 0.4874
millimeters as the bare number might suggest. Any code that converts
`header.get_zooms()` to microns on these files must multiply by 1000.
`polish_annotations.py` sidesteps the ambiguity entirely by taking voxel
spacing as an explicit constant (`SPACING_UM`) rather than trusting the
header.

## Data (`data_for_upload/04_annotation/`)

| File | What it is |
|---|---|
| `Annotations_V3.nii.gz` | original manual tracing, untouched — provenance anchor |
| `Annotations_V3b.nii.gz` | `V3` after isotropic ITK-SNAP label smoothing; input to step 1 |
| `Annotations_V3b_gapfilled2.nii.gz` | after step 1 + the manual correction; input to step 2. The only shipped file recording the manual correction. |
| `Annotations_V3b_final.nii.gz` | finished annotation: gap fill + manual correction + A-P smoothing |
| `labels.txt` | ITK-SNAP label description file (region names, ids, display colors) |

The intensity template itself is not duplicated in this stage's data
directory — it is shipped once, with stage 3
(`data_for_upload/03_structural_to_template_registration/stage0_geometry_fixed/template_ras.nii.gz`),
and referenced from there.

## Requirements

Python 3 with `numpy`, `scipy`, `nibabel`. `matplotlib` is only needed if you
pass `--qc-dir`.

## How to run

The script resolves its data root the same way as the other pipeline stages:
`PIPELINE_DATA_ROOT` environment variable if set, otherwise
`<repo>/data_for_upload` relative to the script's own location. All paths
below are relative to that root's `04_annotation/` (and, for the template,
`03_structural_to_template_registration/stage0_geometry_fixed/`) subdirectory.

Run both steps end-to-end on the shipped default input/output (does **not**
reproduce the manual correction — see caveat above):

```bash
python polish_annotations.py
```

Reproduce the exact pipeline used for the manuscript, with the manual
ITK-SNAP correction step in between:

```bash
python polish_annotations.py Annotations_V3b.nii.gz filled.nii.gz --no-smooth
#   ... hand-correct in ITK-SNAP (Vp reclassification, boundary tweaks) -> corrected.nii.gz ...
python polish_annotations.py corrected.nii.gz Annotations_V3b_final.nii.gz --no-fill
```

Arbitrary input/output paths, and optional QC montages (underlaid on stage
3's template by default):

```bash
python polish_annotations.py in.nii.gz out.nii.gz --qc-dir qc/
```

Run `python polish_annotations.py --help` for the full flag list
(`--no-fill`, `--no-smooth`, `--labels`, `--qc-dir`, `--template`) and the
script's module docstring for the full method description and tunable
parameters (voxel spacing, A-P axis, erosion radius, smoothing sigma, target
labels), all defined in a `CONFIGURATION` block near the top of the file.
