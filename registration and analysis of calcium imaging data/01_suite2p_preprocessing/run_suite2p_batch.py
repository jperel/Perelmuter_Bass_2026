#!/usr/bin/env python3
"""
run_suite2p_batch.py
======================
Stage 1 of 6 (suite2p preprocessing) -- step 2 of 3 in this folder.

Runs suite2p (motion correction + ROI detection + trace extraction)
independently on each of the 36 single-plane functional recordings. Each
recording is a separate 2-minute acquisition at a different native z-value --
not a simultaneous multi-plane volume -- so each gets its own suite2p
run/output folder, keyed by the z-value already present in its filename
(e.g. "02F_2min_2440").

Reads:
    Per-plane TIFF stacks from TIFF_DIR (output of convert_oir_to_tiff.py, the
    previous script in this folder; local/intermediate, not part of the
    public data release).

Writes:
    data_for_upload/06_calcium_region_analysis/suite2p_output/<name>/suite2p/plane0/
    (F.npy, Fneu.npy, spks.npy, stat.npy, iscell.npy, ops.npy), one folder per
    input recording. This is suite2p's own output tree; it is consumed both by
    the calcium/region analysis in stage 6 and by export_mean_images.py (the
    next script in this folder), which reads ops["meanImg"] from it.

Run order in this folder:
    1. convert_oir_to_tiff.py
    2. run_suite2p_batch.py       (this script)
    3. export_mean_images.py

Acquisition/processing parameters (fixed for all 36 planes; do not change
without re-checking the reasoning below against the new data):
    fs       : measured directly from OME metadata frame timestamps (~1.229 Hz,
               i.e. ~813.9 ms between frames) -- not assumed from the
               "2min"/frame-count naming alone.
    tau      : 3.8 s -- nuclear-localized GCaMP6s in zebrafish/Danionella,
               reported literature range 3.5-4.1 s. Cytosolic GCaMP6s decays
               much faster (~1.25-1.5 s); do not reuse that default for
               nuclear-localized indicators.
    diameter : 4 px -- nuclear signal at 1.40625 um/px physical pixel size,
               assuming a ~5-6 um true nuclear diameter for these small
               teleost neurons. This is a starting-point value: inspect the
               detected ROI sizes from an initial run against this assumption
               before trusting it across all planes.

Requires suite2p (pip install suite2p).
"""

import glob
import os

import suite2p

DATA_ROOT = os.environ.get("PIPELINE_DATA_ROOT") or os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data_for_upload"))

HERE = os.path.dirname(os.path.abspath(__file__))

# Local/intermediate input: the per-plane TIFFs produced by convert_oir_to_tiff.py.
# Not part of the public data release; override with FUNCTIONAL_TIFF_DIR if needed.
TIFF_DIR = os.environ.get("FUNCTIONAL_TIFF_DIR", os.path.join(HERE, "functional_z-planes_tiff"))

# suite2p's raw output tree is published as part of stage 6's data.
OUT_BASE = os.path.join(DATA_ROOT, "06_calcium_region_analysis", "suite2p_output")

FS = 1000.0 / 813.89  # frame interval measured from OME metadata (see docstring)
TAU = 3.8
DIAMETER = [4.0, 4.0]


def main():
    os.makedirs(OUT_BASE, exist_ok=True)
    files = sorted(glob.glob(os.path.join(TIFF_DIR, "*.tif")))
    print(f"Found {len(files)} tiff files. fs={FS:.5f} Hz tau={TAU} diameter={DIAMETER}")

    for i, path in enumerate(files):
        name = os.path.splitext(os.path.basename(path))[0]
        save_path0 = os.path.join(OUT_BASE, name)
        print(f"\n=== [{i + 1}/{len(files)}] {name} ===")

        db = {
            "data_path": [TIFF_DIR],
            "file_list": [path],
            "save_path0": save_path0,
            "nplanes": 1,
            "nchannels": 1,
        }
        settings = suite2p.default_settings()
        settings["fs"] = FS
        settings["tau"] = TAU
        settings["diameter"] = DIAMETER

        suite2p.run_s2p(db=db, settings=settings)

    print("\nDone with all files.")


if __name__ == "__main__":
    main()
