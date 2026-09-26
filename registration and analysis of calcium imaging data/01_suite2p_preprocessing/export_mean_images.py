#!/usr/bin/env python3
"""
export_mean_images.py
=======================
Stage 1 of 6 (suite2p preprocessing) -- step 3 of 3 in this folder.

Exports suite2p's registered mean image (ops["meanImg"]) for each of the 36
functional z-planes into a single folder as individual TIFFs, one per plane,
with filenames labeled by the plane's native z-depth (parsed from the
"02F_2min_<z>" suite2p output folder name).

Reads:
    data_for_upload/06_calcium_region_analysis/suite2p_output/<name>/suite2p/plane0/ops.npy
    (written by run_suite2p_batch.py, the previous script in this folder).

Writes:
    data_for_upload/01_suite2p_preprocessing/mean_images/meanImg_z<depth>.tif
    (36 files). Values are saved as-is (float32, no rescaling/clipping) so
    pixel intensities remain quantitatively meaningful for the downstream
    functional-to-structural registration step (stage 2).

Run order in this folder:
    1. convert_oir_to_tiff.py
    2. run_suite2p_batch.py
    3. export_mean_images.py     (this script)

Requires numpy and tifffile.
"""

import os
import re

import numpy as np
import tifffile

DATA_ROOT = os.environ.get("PIPELINE_DATA_ROOT") or os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data_for_upload"))

SUITE2P_OUTPUT_DIR = os.path.join(DATA_ROOT, "06_calcium_region_analysis", "suite2p_output")
OUT_DIR = os.path.join(DATA_ROOT, "01_suite2p_preprocessing", "mean_images")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    names = sorted(os.listdir(SUITE2P_OUTPUT_DIR))

    n_written = 0
    for name in names:
        m = re.search(r"(\d+)$", name)
        if m is None:
            print(f"Skipping {name}: no trailing depth value found")
            continue
        depth = int(m.group(1))

        ops_path = os.path.join(SUITE2P_OUTPUT_DIR, name, "suite2p", "plane0", "ops.npy")
        if not os.path.exists(ops_path):
            print(f"Skipping {name}: no ops.npy found")
            continue

        ops = np.load(ops_path, allow_pickle=True).item()
        mean_img = ops["meanImg"].astype(np.float32)

        out_path = os.path.join(OUT_DIR, f"meanImg_z{depth:04d}.tif")
        tifffile.imwrite(out_path, mean_img)
        n_written += 1
        print(f"{name} (z={depth}) -> {out_path}  shape={mean_img.shape}")

    print(f"\nWrote {n_written} mean image TIFFs to {OUT_DIR}")


if __name__ == "__main__":
    main()
