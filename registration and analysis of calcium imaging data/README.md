# Dracula fish telencephalon 2-photon calcium imaging: registration and region-based analysis

This repository contains the five-stage pipeline
that took raw acquisitions through to the calcium-analysis figures used in
the manuscript (Figure 7, panels G-K).

## Getting the data

This repository is code only. The data package is hosted separately, on Zenodo:

https://zenodo.org/records/22985217?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjNhOWNlYWI3LWM2NDYtNGViNC1iYmJhLWRiNjViYWFiOTM0NCIsImRhdGEiOnt9LCJyYW5kb20iOiI2YzA0OTYzMTA3NTMwOTllNmQzMmE4ZjBkODFiZTUwZCJ9.2JtHyiZCIcNi9Wu5vlNLtekL-kXsg7OfnXQBjsPeLzIsuQGWEJL5FQ9dy8LEsvmXwvVR3S7AiH1X568h6YX1Dg

**to run entire pipeline download all files into a single flat folder named `data_for_upload`, placed next to this
repository's own root (whatever you name the folder this README is in).


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

## Environments

Different stages use different environments (each stage's own README gives
the exact package list):

- **suite2p** environment (stage 1): `suite2p` and its dependencies.
- **antspy** environment (stages 2, 3, 4, 5): `ants` (ANTsPy), `nibabel`,
  `numpy`, `scipy`, `scikit-image`, `tifffile`, `pandas`, `matplotlib`.
- Stages 1 and 2's `.oir`-reading scripts additionally need a JVM and the
  Bio-Formats `bioformats_package.jar` (obtain separately from the official
  Bio-Formats distribution; not redistributed here).

