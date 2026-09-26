# Stage 1 of 6: suite2p preprocessing

Converts raw Olympus .oir 2-photon recordings to TIFF, runs suite2p
(motion correction + ROI detection + trace extraction) independently on each
of 36 single-plane recordings, and exports suite2p's registered mean image
per plane for use by the later registration stages.

Each of the 36 recordings is a separate ~2-minute acquisition at a distinct
native z-depth (typically 10 um spacing between planes, spanning z = 2440-2790
in this dataset), not a simultaneous multi-plane volume. Each therefore gets
its own independent suite2p run.

## Scripts and run order

1. **`convert_oir_to_tiff.py`** -- converts each raw `.oir` file to a
   `(T, Y, X)` TIFF stack via Bio-Formats (jpype + `bioformats_package.jar`).
2. **`run_suite2p_batch.py`** -- runs suite2p on each per-plane TIFF
   independently, writing suite2p's standard output tree (`F.npy`, `Fneu.npy`,
   `spks.npy`, `stat.npy`, `iscell.npy`, `ops.npy`) per plane.
3. **`export_mean_images.py`** -- pulls `ops["meanImg"]` (the registered mean
   image) out of each suite2p output folder and writes it as a standalone
   TIFF, named by the plane's native z-depth.

Run the three scripts in this order. Re-running any script regenerates its
outputs in place.

## Inputs / outputs

- `convert_oir_to_tiff.py` reads raw `.oir` recordings from a local directory
  (default `functional_z-planes_oir/` next to the script; override with the
  `RAW_OIR_DIR` environment variable) and writes per-plane TIFF stacks to a
  local directory (default `functional_z-planes_tiff/`; override with
  `FUNCTIONAL_TIFF_DIR`). Raw recordings and the intermediate TIFFs are not
  part of the public data release (2P raw movies are large); supply your own
  copy of the raw `.oir` files to reproduce this step.
- `run_suite2p_batch.py` reads the TIFFs from the same `FUNCTIONAL_TIFF_DIR`
  and writes suite2p's output tree to
  `data_for_upload/06_calcium_region_analysis/suite2p_output/<recording_name>/suite2p/plane0/`.
  This output is shared with, and consumed by, the stage 6
  (calcium/region analysis) scripts.
- `export_mean_images.py` reads `ops.npy` from that same suite2p output tree
  and writes `data_for_upload/01_suite2p_preprocessing/mean_images/meanImg_z<depth>.tif`
  (36 files), which are consumed by stage 2
  (functional-to-structural registration).

All scripts resolve the shared data root the same way (`PIPELINE_DATA_ROOT`
environment variable if set, otherwise `<repo>/data_for_upload` next to
`<repo>/code`), and build every published path under it with
`os.path.join(DATA_ROOT, "<stage-folder-name>", ...)`.

## Environments / packages

- **`convert_oir_to_tiff.py`**: requires `jpype1`, `numpy`, `tifffile`, a JVM
  (JDK/JRE 8+), and `bioformats_package.jar` (the self-contained Bio-Formats
  distribution, obtainable from the official Bio-Formats releases). Place the
  jar at `bioformats_jar/bioformats_package.jar` next to the script, or point
  at it with `BIOFORMATS_JAR_PATH`. By default the script asks jpype to find
  a JVM automatically; override with `BIOFORMATS_JVM_PATH` if you need a
  specific one.
- **`run_suite2p_batch.py`**: requires `suite2p` (`pip install suite2p`).
- **`export_mean_images.py`**: requires `numpy`, `tifffile`.

## Recording / processing parameters

- 36 native z-planes, ~10 um spacing, ~2 minutes per plane (148 frames per
  recording in this dataset).
- `fs` = 1000 / 813.89 Hz (~1.229 Hz) -- measured directly from OME metadata
  frame timestamps, not assumed from the "2min" naming or frame count alone.
- `tau` = 3.8 s -- appropriate for nuclear-localized GCaMP6s in
  zebrafish/Danionella (reported literature range 3.5-4.1 s). This is much
  slower than cytosolic GCaMP6s (~1.25-1.5 s); do not reuse cytosolic
  defaults for nuclear-localized indicators.
- `diameter` = 4 px -- based on a physical pixel size of 1.40625 um/px and an
  assumed ~5-6 um true nuclear diameter for these small teleost neurons. This
  is a starting-point value only: check detected ROI sizes from an initial
  run before trusting it across all 36 planes.

## Data included in this release

`data_for_upload/01_suite2p_preprocessing/mean_images/meanImg_z<depth>.tif`
(36 files) -- the suite2p-registered mean image for each z-plane, saved as
float32 with no rescaling/clipping, so pixel intensities remain quantitatively
meaningful.
