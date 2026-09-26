# Danionella telencephalon 2-photon calcium imaging: registration and region-based analysis

Code accompanying the manuscript. Danionella (small teleost fish) telencephalon,
2-photon GCaMP calcium imaging across 36 native single-plane z-recordings
(10 um spacing), registered onto a hand-annotated brain-region atlas, and
analyzed region-by-region. This repository contains the six-stage pipeline
that took raw acquisitions through to the calcium-analysis figures used in
the manuscript (Figure 7, panels G-K).

## Getting the data

This repository is code only. The data package it operates on (~2.3 GB —
template volumes, registration transforms, intermediate annotation volumes,
and the suite2p outputs needed to run these scripts) is hosted separately:

**[DATA REPOSITORY LINK / DOI — fill in once the data package is uploaded, e.g. to Zenodo or Dryad]**

Download and extract it so that `data_for_upload/` sits next to `code/`:

```
<wherever you clone/extract this>/
  code/
  data_for_upload/
```

Every script finds the data root automatically from this layout (or via the
`PIPELINE_DATA_ROOT` environment variable if you want to place it elsewhere):

```python
DATA_ROOT = os.environ.get("PIPELINE_DATA_ROOT") or os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data_for_upload"))
```

Scripts write their own outputs back into the same `data_for_upload/<stage>/`
location they read sibling inputs from, so deleting a downstream file and
re-running the script that produces it regenerates it in place.

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
                         stage 4: manual region annotation
                         (drawn on the template, stage 3's output)
                                      |
                                      v
                    stage 5: ROI pullback + cell-region assignment
                    (composes stages 2+3+4 onto native functional pixels,
                     assigns each suite2p cell a brain region)
                                      |
                                      v
                    stage 6: region-grouped calcium analysis
                    (Figure 7 G-K: representative-plane maps/traces,
                     region activity, transient rate, enrichment)
```

| # | Folder | What it does | Depends on |
|---|---|---|---|
| 1 | `01_suite2p_preprocessing` | Raw `.oir` -> TIFF -> suite2p (motion correction, ROI detection, trace extraction) -> per-plane mean images | raw recordings (not included) |
| 2 | `02_functional_to_structural_registration` | Registers the 36-plane functional mini-stack into the native space of a separate, higher-SNR structural stack ("Stage B") | stage 1 |
| 3 | `03_structural_to_template_registration` | Landmark + deformable (SyN) registration of the structural stack to a population-average template with hand-drawn region annotations ("Stage A") | independent of stages 1-2 (uses its own structural-stack input, see note below) |
| 4 | `04_annotation` | Polishes the manually-traced brain-region label volume (gap fill, anisotropic smoothing, one required manual ITK-SNAP correction step) | drawn on stage 3's template |
| 5 | `05_roi_pullback` | Pulls region labels backward through stages 2+3's composed transform chain onto native functional pixels; assigns each suite2p-detected cell a region by majority vote | stages 1, 2, 3, 4 |
| 6 | `06_calcium_region_analysis` | ΔF/F extraction, transient detection, and the specific figures reproducing manuscript Figure 7 G-K | stages 1, 4, 5 |

Each stage folder has its own `README.md` with exact run order, inputs/outputs,
environment requirements, and known caveats for that stage — read those before
running anything.

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
the structural stack once registered. Stage 5 handles this correctly via a
"slab" technique (each native plane embedded as a true 1-voxel-thick 3D image
at its own position, pulled through the full transform chain in one call) —
naively picking "the nearest z-slice" and doing a flat crop would silently
give the wrong answer, worst at the edges of each frame.

## Known limitations (apply across stages 3-6)

- **A local registration gap** exists at one specific atlas boundary (the
  dorsal/top-middle rim), worst at the deepest included imaging planes
  (z >= 2680, worst at z = 2790), tapering off at shallower depth. This is a
  confirmed, unresolved limitation of the structural-to-template registration
  in that one region — not a general problem across the dataset. Cells near
  that boundary at those depths may show inflated "unassigned" labeling or
  occasional wrong-neighbor assignment; see stage 5's and stage 6's READMEs
  for how to check/flag affected cells.
- **Stage 6's analysis is single-animal and descriptive.** No confidence
  intervals or effect sizes are computed anywhere in that stage, by
  deliberate choice — a census of ~all cells in one animal's chosen depth
  series has no population to formally infer to.

## Environments

Different stages use different environments (each stage's own README gives
the exact package list):

- **suite2p** environment (stage 1): `suite2p` and its dependencies.
- **antspy** environment (stages 2, 3, 4, 5, 6): `ants` (ANTsPy), `nibabel`,
  `numpy`, `scipy`, `scikit-image`, `tifffile`, `pandas`, `matplotlib`.
- Stages 1 and 2's `.oir`-reading scripts additionally need a JVM and the
  Bio-Formats `bioformats_package.jar` (obtain separately from the official
  Bio-Formats distribution; not redistributed here).

## Citation

If you use this code or the associated data, please cite the manuscript:
**[MANUSCRIPT CITATION — fill in once available: authors, title, journal, DOI]**
