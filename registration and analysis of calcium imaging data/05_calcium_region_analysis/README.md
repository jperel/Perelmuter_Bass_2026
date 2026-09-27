# Stage 5 of 5: calcium / region analysis

Final analysis and figure-generation stage. Combines suite2p per-plane
fluorescence traces (stage 1) with cell-to-region assignments (stage 4) to
compute per-neuron dF/F, z-scores, and significant-transient rates, and to
produce the figures for manuscript Figure 7.


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
3. **`per_plane_figs.py`**, **`manuscript_panels.py`**, **`enrichment_compact.py`**
   -- can run in any order relative to each other once steps 1-2 are done.

## Figure 7 panel mapping

| Panel | Script | Output file(s) |
|---|---|---|
| G, H | `per_plane_figs.py` | `circlemap_z<depth>.tif`, `traces_z<depth>.tif` (`.pdf`), one pair per plane |
| I | `manuscript_panels.py` | `fig_mean_dff_scatter_violin_A.tif` |
| J | `manuscript_panels.py` | `raincloud_transients_violin.tif` |
| K | `enrichment_compact.py` | `fig_enrichment_dumbbell_transients.tif` (`.pdf`) |

## Inputs / outputs

All paths are resolved relative to a shared data root: the
`PIPELINE_DATA_ROOT` environment variable if set, otherwise `<repo>/data_for_upload`
next to the code repository. **`data_for_upload/` is a flat folder** -- Zenodo
does not preserve directory structure on upload, so every script builds its
data paths with `os.path.join(DATA_ROOT, "<filename>")` directly, with three
exceptions: `mean_images/`, `roi_labels/`, and `suite2p_output/` are
subfolders created locally by unzipping the three `.zip` archives in the data
package.

Reads:
- `suite2p_output/02F_2min_<z>/suite2p/plane0/` -- `F.npy`, `Fneu.npy`,
  `spks.npy`, `stat.npy`, `iscell.npy`, `ops.npy`, 32 planes
  (z = 2480-2790 um). Included in this release.
- `cell_region_assignments.csv`, `roi_labels/roi_labels_z<z>.tif`
- `mean_images/meanImg_z<z>.tif`
- `labels.txt` (region IDs, colors, names -- ITK-SNAP label description
  format)

Writes (all under `outputs/`, a local-only subfolder created by
`run_calcium_region_analysis.py` on first run -- not part of the Zenodo data
package):
- `traces_metadata.csv`, `traces.npz` -- pooled per-neuron dF/F/z-score
  matrices and metadata; the load-bearing output of this stage.
- `trace_plots/`, `circle_maps/` -- per-plane QC figures (lower-resolution
  PNGs, written directly by `run_calcium_region_analysis.py`).
- `per_plane_figs/` -- the publication-quality per-plane TIFF/PDF pairs
  written by `per_plane_figs.py` (`traces_z<depth>.*`, `circlemap_z<depth>.*`,
  32 planes each).
- `transient_rate.csv`, `transient_rate.png`
- The Figure 7 panel files listed above, plus the additional
  `manuscript_panels.py` / `enrichment_compact.py` variants noted above.

