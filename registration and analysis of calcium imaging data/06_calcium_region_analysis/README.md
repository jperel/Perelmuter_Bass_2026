# Stage 6 of 6: calcium / region analysis

Final analysis and figure-generation stage. Combines suite2p per-plane
fluorescence traces (stage 1) with cell-to-region assignments (stage 5) to
compute per-neuron dF/F, z-scores, and significant-transient rates, and to
produce the figures for manuscript Figure 7.

**Scope note:** the original project directory this was extracted from
contains roughly 25 analysis/figure scripts. This release includes only the
5 scripts that actually produced Figure 7, panels G-K. Other exploratory
analyses (region synchrony, cohesion, spatial gradients, additional
raincloud/enrichment variants, per-depth summaries, etc.) are not part of
this release.

## Scripts and run order

1. **`run_calcium_region_analysis.py`** -- base pipeline. Must run first;
   every other script in this stage depends on its `traces.npz` /
   `traces_metadata.csv` output. For each of the 32 included planes
   (z = 2480-2790 um, 10 um steps), computes dF/F from raw suite2p
   fluorescence for cells that pass `iscell==1` and are assigned to a
   retained telencephalon region (5,292 cells, 13 regions retained overall).
   Also writes per-plane trace-stack and circle-map figures (not required
   for Figure 7, but part of this script's normal output and left in place).
2. **`transient_rate.py`** -- significant-transient (event) rate per neuron.
   Must run second; `manuscript_panels.py` and `enrichment_compact.py` both
   need its `transient_rate.csv` output.
3. **`z2620_figs.py`**, **`manuscript_panels.py`**, **`enrichment_compact.py`**
   -- can run in any order relative to each other once steps 1-2 are done.

## Figure 7 panel mapping

| Panel | Script | Output file(s) |
|---|---|---|
| G | `z2620_figs.py` | `circlemap_z2620.tif` (`.pdf`) |
| H | `z2620_figs.py` | `z2620_full.tif` / `z2620_rep24.tif` (`.pdf`) |
| I | `manuscript_panels.py` | `fig_mean_dff_scatter_violin_A.tif` |
| J | `manuscript_panels.py` | `raincloud_transients_violin.tif` |
| K | `enrichment_compact.py` | `fig_enrichment_dumbbell_transients.tif` (`.pdf`) |

`manuscript_panels.py` and `enrichment_compact.py` each write several
additional files that are not used in the manuscript (plain dot-scatter
panels, skew-based variants, the mean-ΔF/F+skew dumbbell, the strip
heat-map). These are left in place because they share plotting code with
the panels that are used, not because they are needed for Figure 7. Note
also that panels I and J come from the same function
(`scatter_panel_violin`) but the two output files follow different naming
conventions (`fig_mean_dff_scatter_violin_A.tif` vs.
`raincloud_transients_violin.tif`) -- this inconsistency is carried over
unchanged from the original analysis rather than corrected.

## Inputs / outputs

All paths are resolved relative to a shared data root: the
`PIPELINE_DATA_ROOT` environment variable if set, otherwise `<repo>/data_for_upload`
next to `<repo>/code`. Every script builds its data paths with
`os.path.join(DATA_ROOT, "<stage-folder-name>", ...)`.

Reads:
- `06_calcium_region_analysis/suite2p_output/02F_2min_<z>/suite2p/plane0/`
  -- `F.npy`, `Fneu.npy`, `spks.npy`, `stat.npy`, `iscell.npy`, `ops.npy`,
  32 planes (z = 2480-2790 um). Included in this release.
- `05_roi_pullback/cell_region_assignments.csv`, `05_roi_pullback/roi_labels/roi_labels_z<z>.tif`
- `01_suite2p_preprocessing/mean_images/meanImg_z<z>.tif`
- `04_annotation/labels.txt` (region IDs, colors, names -- ITK-SNAP label
  description format)

Writes (all under `06_calcium_region_analysis/outputs/`, created by
`run_calcium_region_analysis.py` on first run):
- `traces_metadata.csv`, `traces.npz` -- pooled per-neuron dF/F/z-score
  matrices and metadata; the load-bearing output of this stage.
- `trace_plots/`, `circle_maps/` -- per-plane QC figures.
- `transient_rate.csv`, `transient_rate.png`
- The Figure 7 panel files listed above, plus the additional
  `manuscript_panels.py` / `enrichment_compact.py` variants noted above.

## Environment / packages

All 5 scripts run in the same environment (referred to here as `antspy`,
matching the environment name used for registration in stages 2-5):
`numpy`, `scipy`, `pandas`, `matplotlib`, `tifffile`. `z2620_figs.py`
additionally requires `scikit-image` (contour extraction for the circle-map
region boundaries).

The functional recordings in this dataset use a nuclear-localized GCaMP
indicator (see stage 1's README for the tau/diameter parameters this
implies), which is why `run_calcium_region_analysis.py` reads raw `F.npy`
directly with no neuropil subtraction, and why the baseline windows used
throughout this stage are wide relative to what is typical for fast
cytosolic indicators.

## Method summary

**dF/F** (`run_calcium_region_analysis.py`): primary baseline is a centered
rolling 8th-percentile filter, 40 s window, computed on raw F (percentile
filter, reflect-padded edges), following the ~8th-percentile sliding-window
convention used by Dombeck et al. 2010 (Nat Neurosci), CaImAn's
`detrend_df_f` (`quantileMin=8`), and suite2p's own constant-percentile
baseline (`prctile_baseline=8`). The 40 s window is much wider than typical
for fast cytosolic indicators; it is chosen here because the indicator is
slow (nuclear-localized GCaMP6s, tau ~3.8 s, events lasting 5-15 s) relative
to the ~120 s recording length. A cross-check baseline
(multiplicative linear detrend for bleaching, then a static percentile F0)
is also computed. Both are median-filtered (3 frames, ~2.4 s) after
converting to percent.

**Significant-transient detection** (`transient_rate.py`): re-derives dF/F
from raw F using a *symmetric* baseline (linear detrend, then static median
F0) instead of the rolling-percentile baseline above, because the
rolling-percentile baseline clips the negative deflections needed to
calibrate noise. Noise sigma is taken from the noise-dominated lower half of
the dF/F distribution (median minus the 15.87th percentile, i.e. ~1 sigma
below center) rather than from a frame-to-frame difference estimator, which
is deflated by trace autocorrelation by roughly 4x and yields a
false-positive rate around 50% if used naively. Events are detected with
`scipy.signal.find_peaks` at `height >= K*sigma`, `prominence >= 1.5*sigma`,
minimum spacing of 4 frames, minimum width 1 frame; the same detector
applied to the negated trace gives the false-positive floor. K is swept and
the smallest K with a dataset-wide false-positive fraction <= 0.10 is
selected (K = 3 sigma, ~4% false-positive fraction in this dataset).

## Caveats

- **Single animal, descriptive analysis only.** No confidence intervals or
  effect sizes are computed anywhere in this stage. This is deliberate: a
  single census of essentially all cells in one chosen, non-random depth
  series has no population to infer to, so no inferential statistics are
  reported.
- **Deep-plane region-assignment uncertainty.** Cells in the deepest
  included planes (z >= 2680 um, worst at z = 2790 um) sit near a known
  unresolved atlas-registration gap (see stage 5's README) and may show
  inflated "Clear Label" assignment or occasional wrong-neighbor region
  assignment. Treat region assignments in that depth range with more
  caution than shallower planes.
