#!/usr/bin/env python3
"""
export_per_plane_qc_with_suite2p_rois.py
============================================
Stage 5 of 6 (roi_pullback), QC visualization.

For each included plane, renders the mean image with the pulled-back region
boundaries overlaid, plus suite2p's own detected cell ROIs painted on top:
  - only cells with iscell == 1
  - GREEN if the cell's assigned region (majority vote, from
    cell_region_assignments.csv -- the same assignment assign_cells_to_regions.py
    already computed, reused here for consistency rather than recomputed) is a real
    labeled region
  - RED if it landed in "Clear Label" (region 0, i.e. outside any drawn region) or
    "ventricle" (region 19)

Cell footprints (ypix/xpix) come directly from suite2p's stat.npy, filled at low alpha
so the underlying mean image and region boundaries stay visible underneath.

This is the QC step to use for sanity-checking region-assignment quality before
trusting a given plane's numbers -- in particular for checking the known top-middle
rim misalignment described in this stage's README at the affected depths.

Reads (relative to DATA_ROOT):
    01_suite2p_preprocessing/mean_images/meanImg_z<depth>.tif
    05_roi_pullback/roi_labels/roi_labels_z<depth>.tif
    04_annotation/labels.txt
    06_calcium_region_analysis/suite2p_output/02F_2min_<depth>/suite2p/plane0/stat.npy
    05_roi_pullback/cell_region_assignments.csv

Writes:
    05_roi_pullback/qc/roi_qc_with_cells_z<depth>.png

Run after ``assign_cells_to_regions.py``.

Requires: matplotlib, numpy, pandas, tifffile, scipy.
"""

import glob
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from scipy.ndimage import binary_erosion

DATA_ROOT = os.environ.get("PIPELINE_DATA_ROOT") or os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data_for_upload"))

MEAN_IMG_DIR = os.path.join(DATA_ROOT, "01_suite2p_preprocessing", "mean_images")
LABEL_DIR = os.path.join(DATA_ROOT, "05_roi_pullback", "roi_labels")
LABELS_TXT = os.path.join(DATA_ROOT, "04_annotation", "labels.txt")
SUITE2P_OUTPUT_DIR = os.path.join(DATA_ROOT, "06_calcium_region_analysis", "suite2p_output")
ASSIGNMENTS_CSV = os.path.join(DATA_ROOT, "05_roi_pullback", "cell_region_assignments.csv")
OUT_DIR = os.path.join(DATA_ROOT, "05_roi_pullback", "qc")

RED_REGION_IDS = {0, 19}  # "Clear Label", "ventricle"
CELL_ALPHA = 0.55


def parse_labels_txt(path):
    colors, names = {}, {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(None, 7)
            idx = int(parts[0])
            r, g, b = int(parts[1]), int(parts[2]), int(parts[3])
            name = parts[7].strip('"') if len(parts) > 7 else str(idx)
            colors[idx] = (r / 255.0, g / 255.0, b / 255.0)
            names[idx] = name
    return colors, names


def norm(img, pct=99):
    hi = np.percentile(img[img > 0], pct) if np.any(img > 0) else img.max()
    return np.clip(img / max(hi, 1e-6), 0, 1)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    colors, names = parse_labels_txt(LABELS_TXT)
    assignments = pd.read_csv(ASSIGNMENTS_CSV)

    label_paths = sorted(glob.glob(os.path.join(LABEL_DIR, "roi_labels_z*.tif")))
    for label_path in label_paths:
        m = re.search(r"roi_labels_z(\d+)\.tif$", label_path)
        depth = int(m.group(1))

        mean_img = tifffile.imread(os.path.join(MEAN_IMG_DIR, f"meanImg_z{depth}.tif"))
        label_img = tifffile.imread(label_path)

        gray = norm(mean_img.astype(np.float32))
        rgb = np.stack([gray, gray, gray], axis=-1)

        # region boundaries
        present_labels = []
        for lbl in np.unique(label_img):
            if lbl == 0:
                continue
            mask = label_img == lbl
            edge = mask & ~binary_erosion(mask)
            color = colors.get(int(lbl), (1, 1, 1))
            rgb[edge] = color
            present_labels.append(int(lbl))

        # suite2p cell ROIs
        plane_dir = os.path.join(SUITE2P_OUTPUT_DIR, f"02F_2min_{depth}", "suite2p", "plane0")
        stat_path = os.path.join(plane_dir, "stat.npy")
        if not os.path.exists(stat_path):
            print(f"z={depth}: SKIPPED (no stat.npy)")
            continue
        stat = np.load(stat_path, allow_pickle=True)

        plane_assign = assignments[(assignments["z"] == depth) & (assignments["iscell"] == 1)]
        assign_by_cell = {row.cell_id: row.region_id for row in plane_assign.itertuples()}

        n_green, n_red = 0, 0
        for cell_id, s in enumerate(stat):
            if cell_id not in assign_by_cell:
                continue  # not iscell==1
            region_id = assign_by_cell[cell_id]
            cell_color = np.array([1.0, 0.0, 0.0]) if region_id in RED_REGION_IDS else np.array([0.0, 1.0, 0.0])
            if region_id in RED_REGION_IDS:
                n_red += 1
            else:
                n_green += 1

            ypix, xpix = s["ypix"], s["xpix"]
            valid = (ypix >= 0) & (ypix < rgb.shape[0]) & (xpix >= 0) & (xpix < rgb.shape[1])
            ypix, xpix = ypix[valid], xpix[valid]
            rgb[ypix, xpix] = (1 - CELL_ALPHA) * rgb[ypix, xpix] + CELL_ALPHA * cell_color

        fig, ax = plt.subplots(figsize=(9, 7))
        ax.imshow(rgb, origin="upper")
        ax.set_title(f"z = {depth}   (suite2p ROIs: {n_green} green / {n_red} red)")
        ax.axis("off")

        handles = [plt.Line2D([0], [0], color=colors[i], lw=3, label=names[i])
                   for i in present_labels]
        handles.append(plt.Line2D([0], [0], marker="s", color="none", markerfacecolor="lime",
                                   markersize=10, label="cell: in region"))
        handles.append(plt.Line2D([0], [0], marker="s", color="none", markerfacecolor="red",
                                   markersize=10, label="cell: clear/ventricle"))
        ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8)

        fig.tight_layout()
        out_path = os.path.join(OUT_DIR, f"roi_qc_with_cells_z{depth}.png")
        fig.savefig(out_path, dpi=130, bbox_inches="tight")
        plt.close(fig)
        print(f"z={depth}: wrote {out_path}  ({n_green} green, {n_red} red)")


if __name__ == "__main__":
    main()
