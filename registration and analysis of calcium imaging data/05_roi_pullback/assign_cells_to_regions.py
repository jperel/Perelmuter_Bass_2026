#!/usr/bin/env python3
"""
assign_cells_to_regions.py
=============================
Stage 5 of 6 (roi_pullback), step 2.

For each included plane's suite2p-detected cells, assigns a brain region by majority
vote over each cell's own pixel footprint against the step-1 pulled-back label image
(more robust than a centroid-only lookup for cells straddling a region boundary), and
records ``purity``: the fraction of the cell's pixels that agree with the assigned
region.

Reads (relative to DATA_ROOT):
    05_roi_pullback/roi_labels/roi_labels_z<depth>.tif        (from pull_roi_labels.py)
    04_annotation/labels.txt                                  (region ID -> name lookup)
    06_calcium_region_analysis/suite2p_output/02F_2min_<depth>/suite2p/plane0/stat.npy
    06_calcium_region_analysis/suite2p_output/02F_2min_<depth>/suite2p/plane0/iscell.npy

Writes:
    05_roi_pullback/cell_region_assignments.csv

Columns: z, cell_id (suite2p's own per-plane index; (z, cell_id) together are the
unique key), iscell, iscell_prob, n_pixels, region_id (0 = "Clear Label", not a real
region), region_name, purity.

Run after ``pull_roi_labels.py``; run before ``export_per_plane_qc_with_suite2p_rois.py``.

Requires: numpy, tifffile.
"""

import csv
import glob
import os
import re

import numpy as np
import tifffile

DATA_ROOT = os.environ.get("PIPELINE_DATA_ROOT") or os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data_for_upload"))

LABEL_DIR = os.path.join(DATA_ROOT, "05_roi_pullback", "roi_labels")
SUITE2P_OUTPUT_DIR = os.path.join(DATA_ROOT, "06_calcium_region_analysis", "suite2p_output")
LABELS_TXT = os.path.join(DATA_ROOT, "04_annotation", "labels.txt")
OUT_CSV = os.path.join(DATA_ROOT, "05_roi_pullback", "cell_region_assignments.csv")


def parse_labels_txt(path):
    names = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(None, 7)
            idx = int(parts[0])
            name = parts[7].strip('"') if len(parts) > 7 else str(idx)
            names[idx] = name
    return names


def main():
    names = parse_labels_txt(LABELS_TXT)

    label_paths = sorted(glob.glob(os.path.join(LABEL_DIR, "roi_labels_z*.tif")))
    assert label_paths, f"no label images found in {LABEL_DIR}"

    rows = []
    for label_path in label_paths:
        m = re.search(r"roi_labels_z(\d+)\.tif$", label_path)
        depth = int(m.group(1))
        label_img = tifffile.imread(label_path)  # (Y, X), matches suite2p's own convention

        plane_dir = os.path.join(SUITE2P_OUTPUT_DIR, f"02F_2min_{depth}", "suite2p", "plane0")
        stat_path = os.path.join(plane_dir, "stat.npy")
        iscell_path = os.path.join(plane_dir, "iscell.npy")
        if not os.path.exists(stat_path):
            print(f"z={depth}: SKIPPED, no stat.npy at {stat_path}")
            continue

        stat = np.load(stat_path, allow_pickle=True)
        iscell = np.load(iscell_path, allow_pickle=True)

        for cell_id, s in enumerate(stat):
            ypix, xpix = s["ypix"], s["xpix"]
            valid = (ypix >= 0) & (ypix < label_img.shape[0]) & (xpix >= 0) & (xpix < label_img.shape[1])
            ypix, xpix = ypix[valid], xpix[valid]
            if len(ypix) == 0:
                region_id, purity = -1, 0.0
            else:
                pixel_labels = label_img[ypix, xpix]
                vals, counts = np.unique(pixel_labels, return_counts=True)
                top = np.argmax(counts)
                region_id = int(vals[top])
                purity = float(counts[top]) / len(pixel_labels)

            rows.append({
                "z": depth,
                "cell_id": cell_id,
                "iscell": int(iscell[cell_id, 0]),
                "iscell_prob": float(iscell[cell_id, 1]),
                "n_pixels": len(ypix),
                "region_id": region_id,
                "region_name": names.get(region_id, "outside_coverage" if region_id == 0 else "unknown"),
                "purity": round(purity, 3),
            })

        n_cells = len(stat)
        n_labeled = sum(1 for r in rows if r["z"] == depth and r["region_id"] > 0)
        print(f"z={depth}: {n_cells} cells, {n_labeled} assigned to a labeled region")

    with open(OUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {len(rows)} cell rows to {OUT_CSV}")


if __name__ == "__main__":
    main()
