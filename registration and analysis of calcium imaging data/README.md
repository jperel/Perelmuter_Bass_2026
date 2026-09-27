# Danionella telencephalon 2-photon calcium imaging: registration and region-based analysis

Code accompanying the manuscript. Danionella (small teleost fish) telencephalon,
2-photon GCaMP calcium imaging across 36 native single-plane z-recordings
(10 um spacing), registered onto a hand-annotated brain-region atlas, and
analyzed region-by-region. This repository contains the five-stage pipeline
that took raw acquisitions through to the calcium-analysis figures used in
the manuscript (Figure 7, panels G-K).

## Getting the data

This repository is code only. The data package it operates on (~2.1 GB —
template volumes, registration transforms, the annotation volume, and the
suite2p outputs needed to run these scripts) is hosted separately, on
Zenodo:

https://zenodo.org/records/22985217?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjNhOWNlYWI3LWM2NDYtNGViNC1iYmJhLWRiNjViYWFiOTM0NCIsImRhdGEiOnt9LCJyYW5kb20iOiI2YzA0OTYzMTA3NTMwOTllNmQzMmE4ZjBkODFiZTUwZCJ9.2JtHyiZCIcNi9Wu5vlNLtekL-kXsg7OfnXQBjsPeLzIsuQGWEJL5FQ9dy8LEsvmXwvVR3S7AiH1X568h6YX1Dg

**Zenodo does not preserve folder structure on upload — every file in the
deposit is flat**, regardless of how it was organized before uploading. This
data package was prepared as a nested folder tree (one subfolder per
pipeline stage), but Zenodo only keeps the individual files; download them
all into a single flat folder named `data_for_upload`, placed next to this
repository's own root (whatever you name the folder this README is in):

```
<a common parent folder>/
  <this repository>/          (i.e. wherever this README.md lives)
  data_for_upload/            (flat -- every downloaded file goes directly here)
```

Every script finds the data root automatically from this layout (or via the
`PIPELINE_DATA_ROOT` environment variable if you want to place it elsewhere):

```python
DATA_ROOT = os.environ.get("PIPELINE_DATA_ROOT") or os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data_for_upload"))
```

This resolves relative to each script's own location (two directories up:
`<this repo>/<stage>/script.py` -> `<parent>/data_for_upload`), so it works
regardless of what this repository's own folder is named.

Every script in this pipeline builds its data paths as `os.path.join(DATA_ROOT,
"<filename>")` directly -- there are no stage-numbered subfolders in
`data_for_upload/`, to match what Zenodo actually stores. Scripts write their
own outputs back into `DATA_ROOT` (or a small number of locally-created
subfolders, listed below), so deleting a downstream file and re-running the
script that produces it regenerates it in place.

**Zenodo also caps deposits at 100 files.** The three groups of many small,
homogeneous files are shipped as single `.zip` archives instead (20 files
total in the deposit: 17 individual files + these 3 archives). After
downloading everything into `data_for_upload/`, unzip each of these **in
place inside that same flat folder** — each archive's own top-level folder
(e.g. `mean_images/`) lands directly inside `data_for_upload/`, giving you
the one exception to the "everything is flat" rule:

| Archive | Unzips to |
|---|---|
| `mean_images.zip` | `data_for_upload/mean_images/` (36 TIFFs) |
| `roi_labels.zip` | `data_for_upload/roi_labels/` (32 TIFFs) |
| `suite2p_output.zip` | `data_for_upload/suite2p_output/` (32 planes x 6 files) |

e.g. from inside `data_for_upload/`: `unzip mean_images.zip` (or, in
PowerShell, `Expand-Archive mean_images.zip .`). Everything else in the data
package is an individual file living directly in `data_for_upload/`.

Scripts also create a handful of their own local-only subfolders as they
run (not part of the Zenodo deposit) for generated results, e.g.
`data_for_upload/outputs/` (stage 5's figures/traces),
`data_for_upload/qc_structural_to_template/` (stage 3's QC overlays), and
`data_for_upload/qc_roi_pullback/` (stage 4's QC overlays) — each stage's
own README documents its specific generated subfolders.

## Pipeline overview

```
 stage 1                stage 2                  stage 3
 suite2p            functional plane -->      structural stack -->
 preprocessing  -->  structural stack          template (region atlas)
 (raw -> traces)      registration              registration
                      ("Stage B")               ("Stage A")
                            \                        /
                             \                      /
                              v                    v
                    stage 4: ROI pullback + cell-region assignment
                    (composes stages 2+3 onto native functional pixels using
                     annotations.nii.gz -- a hand-traced region-label volume
                     included as this stage's own input data, drawn on stage
                     3's template; tracing/polishing methodology is in the
                     manuscript's Methods -- then assigns each suite2p cell
                     a brain region)
                                      |
                                      v
                    stage 5: region-grouped calcium analysis
                    (Figure 7 G-K: representative-plane maps/traces,
                     region activity, transient rate, enrichment)
```

| # | Folder | What it does | Depends on |
|---|---|---|---|
| 1 | `01_suite2p_preprocessing` | Raw `.oir` -> TIFF -> suite2p (motion correction, ROI detection, trace extraction) -> per-plane mean images | raw recordings (not included) |
| 2 | `02_functional_to_structural_registration` | Registers the 36-plane functional mini-stack into the native space of a separate, higher-SNR structural stack ("Stage B") | stage 1 |
| 3 | `03_structural_to_template_registration` | Landmark + deformable (SyN) registration of the structural stack to a population-average template with hand-drawn region annotations ("Stage A") | independent of stages 1-2 (uses its own structural-stack input, see note below) |
| 4 | `04_roi_pullback` | Pulls region labels backward through stages 2+3's composed transform chain onto native functional pixels using `annotations.nii.gz` (a hand-traced brain-region label volume shipped as this stage's own input data -- see the manuscript's Methods for how it was produced); assigns each suite2p-detected cell a region by majority vote | stages 1, 2, 3 |
| 5 | `05_calcium_region_analysis` | ΔF/F extraction, transient detection, and the figures for manuscript Figure 7 G-K (see this stage's README for panels G-H, which are generalized to every plane rather than reproducing the exact hand-picked panel layout) | stages 1, 4 |

Each stage subfolder has its own `README.md` with exact run order,
inputs/outputs, environment requirements, and known caveats for that stage —
read those before running anything.

## Two things worth knowing before you start

**A registration script requests a "Rigid" refinement step it never actually
performs.** In stage 3's `register_structural_to_template.py`, the landmark-affine
fit is followed by a requested Rigid-refinement stage before SyN. That Rigid
optimizer is confirmed to silently fail to converge on this data — ANTs
catches the internal optimizer exception and returns the unmodified input
transform with no error surfaced to the caller. The accepted registration is
therefore, in practice, landmark-affine + SyN, not landmark + Rigid + SyN,
despite what the script requests. See that stage's README for the full
explanation; this doesn't need fixing to reproduce the manuscript's results
(the shipped transforms already reflect this), but it matters if you adapt
the script.

**A native functional plane maps to a *tilted* plane in registered space, not
a flat z-slice.** Stage 2's fitted transform has a genuine ~19-degree rotation
mixing two axes, so a single optical section spreads across many z-indices of
the structural stack once registered. Stage 4 handles this correctly via a
"slab" technique (each native plane embedded as a true 1-voxel-thick 3D image
at its own position, pulled through the full transform chain in one call) —
naively picking "the nearest z-slice" and doing a flat crop would silently
give the wrong answer, worst at the edges of each frame.

## Known limitations (apply across stages 3-5)

- **A local registration gap** exists at one specific atlas boundary (the
  dorsal/top-middle rim), worst at the deepest included imaging planes
  (z >= 2680, worst at z = 2790), tapering off at shallower depth. This is a
  confirmed, unresolved limitation of the structural-to-template registration
  in that one region — not a general problem across the dataset. Cells near
  that boundary at those depths may show inflated "unassigned" labeling or
  occasional wrong-neighbor assignment; see stage 4's and stage 5's READMEs
  for how to check/flag affected cells.
- **Stage 5's analysis is single-animal and descriptive.** No confidence
  intervals or effect sizes are computed anywhere in that stage, by
  deliberate choice — a census of ~all cells in one animal's chosen depth
  series has no population to formally infer to.

## Environments

Different stages use different environments (each stage's own README gives
the exact package list):

- **suite2p** environment (stage 1): `suite2p` and its dependencies.
- **antspy** environment (stages 2, 3, 4, 5): `ants` (ANTsPy), `nibabel`,
  `numpy`, `scipy`, `scikit-image`, `tifffile`, `pandas`, `matplotlib`.
- Stages 1 and 2's `.oir`-reading scripts additionally need a JVM and the
  Bio-Formats `bioformats_package.jar` (obtain separately from the official
  Bio-Formats distribution; not redistributed here).

## Citation

If you use this code or the associated data, please cite the manuscript:
**[MANUSCRIPT CITATION — fill in once available: authors, title, journal, DOI]**
