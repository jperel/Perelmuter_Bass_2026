#!/usr/bin/env python3
r"""
run_calcium_region_analysis.py
===============================
Stage 6 of 6 (calcium / region analysis) -- run first in this folder.

This stage is scoped to reproduce Figure 7 panels G-K of the manuscript, not
the full original analysis suite. This script computes per-neuron dF/F
traces and z-scores for every suite2p-classified cell assigned to a retained
telencephalon region, then writes the pooled `traces.npz` /
`traces_metadata.csv` outputs that every other script in this stage
(transient_rate.py, z2620_figs.py, manuscript_panels.py,
enrichment_compact.py) depends on.

For every included native plane (z = 2480..2790 um, 10 um steps, 32 planes):
  1. Load suite2p raw fluorescence F.npy (nuclear-localized GCaMP -> no
     neuropil subtraction).
  2. Keep cells with suite2p iscell==1 whose majority-vote region is not
     Clear Label (0), V-unspecified (13), POA (16), Hab (17), ventricle (19),
     or Diencephalon-unspecified (20), and whose region_id >= 1 (a valid
     footprint). This retains 5,292 cells across 32 planes, 13 regions.
  3. dF/F:
       primary (dff_pct)             = (F - F0)/F0, F0 = rolling 8th
                                        percentile, 40 s centered window.
       cross-check (dff_detrend_pct) = multiplicative linear detrend
                                        (bleaching correction), then a
                                        static 8th-percentile F0.
     Both are expressed in percent, then median-filtered (kernel = 3 frames,
     ~2.4 s at this dataset's frame rate).
  4. z-score, computed on the pooled dataset (every neuron, every plane):
       zscore        = per-neuron (x - mean) / std
       zscore_global = per-neuron mean removed, divided by one pooled SD
                       (keeps relative amplitude differences between
                       neurons, unlike per-neuron zscore)
  5. Per plane (kept for completeness -- not required for Figure 7, but part
     of this script's normal output and not stripped):
       outputs/trace_plots/traces_z<z>.png    stacked % dF/F traces, cells
           sorted region -> centroid, colored along the turbo colormap
       outputs/circle_maps/circlemap_z<z>.png mean image + pulled-back
           region boundaries, one uniform circle per kept cell (radius =
           median area-equivalent ROI radius over all iscell==1 cells)
  6. Dataset-wide outputs (required by every downstream script in this
     stage):
       outputs/traces_metadata.csv   one row per kept neuron, row-aligned
                                     to the .npz matrices
       outputs/traces.npz            time_s + (N, T) matrices: dff_pct,
                                     dff_detrend_pct, zscore, zscore_global

Baseline choice: an ~8th-percentile sliding-window baseline follows Dombeck
et al. 2010 (Nat Neurosci), CaImAn's detrend_df_f (quantileMin=8), and
suite2p's own constant-percentile baseline (prctile_baseline=8). The window
is widened to 40 s here, well beyond what is typical for fast cytosolic
indicators, because the indicator in this dataset is slow (nuclear-localized
GCaMP6s, tau ~3.8 s, events lasting 5-15 s) and each recording is only ~120 s
long, so a short window would ride the transients themselves instead of
tracking the baseline.

Reads (relative to DATA_ROOT):
    06_calcium_region_analysis/suite2p_output/02F_2min_<z>/suite2p/plane0/
        F.npy, stat.npy, iscell.npy, ops.npy
    05_roi_pullback/cell_region_assignments.csv
    05_roi_pullback/roi_labels/roi_labels_z<z>.tif
    01_suite2p_preprocessing/mean_images/meanImg_z<z>.tif
    04_annotation/labels.txt

Writes (relative to DATA_ROOT/06_calcium_region_analysis/outputs/):
    traces_metadata.csv, traces.npz,
    trace_plots/traces_z<z>.png, circle_maps/circlemap_z<z>.png

Run order in this stage:
    1. run_calcium_region_analysis.py   (this script)
    2. transient_rate.py
    3. z2620_figs.py / manuscript_panels.py / enrichment_compact.py (any order)

Requires the antspy environment: numpy, scipy, pandas, matplotlib, tifffile.
"""

import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile
from matplotlib.patches import Circle
from scipy.ndimage import binary_erosion, percentile_filter
from scipy.signal import medfilt

# ----------------------------------------------------------------------------------
# paths
# ----------------------------------------------------------------------------------
DATA_ROOT = os.environ.get("PIPELINE_DATA_ROOT") or os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data_for_upload"))

LABELS_TXT = os.path.join(DATA_ROOT, "04_annotation", "labels.txt")
ASSIGN_CSV = os.path.join(DATA_ROOT, "05_roi_pullback", "cell_region_assignments.csv")
LABEL_DIR = os.path.join(DATA_ROOT, "05_roi_pullback", "roi_labels")
MEAN_IMG_DIR = os.path.join(DATA_ROOT, "01_suite2p_preprocessing", "mean_images")
S2P_DIR = os.path.join(DATA_ROOT, "06_calcium_region_analysis", "suite2p_output")

OUT_DIR = os.path.join(DATA_ROOT, "06_calcium_region_analysis", "outputs")
TRACE_DIR = os.path.join(OUT_DIR, "trace_plots")
CIRCLE_DIR = os.path.join(OUT_DIR, "circle_maps")
for d in (TRACE_DIR, CIRCLE_DIR):
    os.makedirs(d, exist_ok=True)

# ----------------------------------------------------------------------------------
# analysis parameters
# ----------------------------------------------------------------------------------
EXCLUDE_REGION_IDS = {0, 13, 16, 17, 19, 20}   # Clear Label, V-unspecified, POA, Hab, ventricle, Dien-unspec
BASELINE_WIN_S = 40.0                       # rolling-percentile window (seconds)
BASELINE_PCTL = 8.0                         # percentile for F0
MEDFILT_KERNEL = 3                          # frames (odd); ~2.4 s at 1.23 Hz
F0_FLOOR = 1.0                              # guard against divide-by-~0
SCALE_BAR_PCT = 200.0                       # vertical scale bar on the trace plots
BAND_ALPHA = 0.13                           # region background band opacity (trace plots)
MAP_REGION_FILL_ALPHA = 0.12               # region area tint on the circle maps


def region_random_colors(k, seed):
    """suite2p-style maximally-separable random colours for k neurons: evenly-spaced
    hues (guarantees separation) shuffled into a random order, with saturation/value
    jitter kept in a mid range so the colours read on BOTH the light trace-plot bands
    and the dark circle-map background."""
    rng = np.random.default_rng(seed)
    hues = (np.arange(k) / max(k, 1) + rng.random()) % 1.0
    rng.shuffle(hues)
    sat = 0.65 + 0.30 * rng.random(k)
    val = 0.72 + 0.23 * rng.random(k)
    return matplotlib.colors.hsv_to_rgb(np.stack([hues, sat, val], axis=1))


def parse_labels_txt(path):
    colors, names = {}, {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            p = line.split(None, 7)
            idx = int(p[0])
            colors[idx] = (int(p[1]) / 255.0, int(p[2]) / 255.0, int(p[3]) / 255.0)
            names[idx] = p[7].strip('"') if len(p) > 7 else str(idx)
    return colors, names


def norm_img(img, pct=99):
    img = img.astype(np.float32)
    hi = np.percentile(img[img > 0], pct) if np.any(img > 0) else img.max()
    return np.clip(img / max(hi, 1e-6), 0, 1)


def dff_rolling_pctl(F, win_frames):
    """F: (n_cells, n_frames) raw. F0 = centred rolling `BASELINE_PCTL` percentile."""
    F0 = percentile_filter(F, percentile=BASELINE_PCTL, size=(1, win_frames), mode="reflect")
    F0 = np.maximum(F0, F0_FLOOR)
    return (F - F0) / F0 * 100.0


def dff_detrend_static(F, t):
    """Multiplicative linear-detrend (bleaching) per cell, then static percentile F0."""
    out = np.empty_like(F, dtype=np.float64)
    A = np.vstack([t, np.ones_like(t)]).T
    for i, tr in enumerate(F):
        a, b = np.linalg.lstsq(A, tr, rcond=None)[0]
        fit = a * t + b
        if fit.min() > F0_FLOOR:
            tr_dt = tr * fit[0] / fit           # anchor gain = 1 at t=0
        else:
            tr_dt = tr.astype(np.float64)       # pathological fit -> skip detrend
        F0 = max(np.percentile(tr_dt, BASELINE_PCTL), F0_FLOOR)
        out[i] = (tr_dt - F0) / F0 * 100.0
    return out


# ----------------------------------------------------------------------------------
# pass 1: per plane -> dF/F + per-plane figures, accumulate matrices
# ----------------------------------------------------------------------------------
def main():
    colors, names = parse_labels_txt(LABELS_TXT)
    assign = pd.read_csv(ASSIGN_CSV)

    depths = sorted(
        int(re.search(r"roi_labels_z(\d+)\.tif$", f).group(1))
        for f in os.listdir(LABEL_DIR)
        if re.search(r"roi_labels_z(\d+)\.tif$", f)
    )
    print(f"{len(depths)} planes: {depths[0]}..{depths[-1]}")

    # --- median ROI diameter over ALL iscell==1 cells (for the uniform circle) ---
    diam = []
    for z in depths:
        sp = os.path.join(S2P_DIR, f"02F_2min_{z}", "suite2p", "plane0")
        st = np.load(os.path.join(sp, "stat.npy"), allow_pickle=True)
        ic = np.load(os.path.join(sp, "iscell.npy"))[:, 0] == 1
        for s, keep in zip(st, ic):
            if keep:
                diam.append(2.0 * np.sqrt(s["npix"] / np.pi))
    median_diam = float(np.median(diam))
    circle_radius = median_diam / 2.0
    print(f"median area-equivalent ROI diameter (all iscell==1, n={len(diam)}): "
          f"{median_diam:.2f} px  -> circle radius {circle_radius:.2f} px")

    meta_rows = []
    dff_list, dffd_list = [], []          # each entry (n_keep_plane, T)
    plane_cache = []                       # (z, keep_df, stat, fs) for pass-2 figures

    for z in depths:
        sp = os.path.join(S2P_DIR, f"02F_2min_{z}", "suite2p", "plane0")
        F = np.load(os.path.join(sp, "F.npy"))                       # (n_cells, T) raw
        stat = np.load(os.path.join(sp, "stat.npy"), allow_pickle=True)
        ops = np.load(os.path.join(sp, "ops.npy"), allow_pickle=True).item()
        fs = float(ops["fs"])
        n_cells, T = F.shape
        assert len(stat) == n_cells

        pa = assign[(assign.z == z) & (assign.iscell == 1)].copy()
        pa = pa[~pa.region_id.isin(EXCLUDE_REGION_IDS) & (pa.region_id >= 1)]
        pa = pa[pa.cell_id < n_cells]

        # centroids (stat['med'] = [y, x]) and sort: region, then y, then x
        pa["cy"] = [float(stat[c]["med"][0]) for c in pa.cell_id]
        pa["cx"] = [float(stat[c]["med"][1]) for c in pa.cell_id]
        pa = pa.sort_values(["region_id", "cy", "cx"]).reset_index(drop=True)
        n_keep = len(pa)

        idx = pa.cell_id.to_numpy()
        Fk = F[idx].astype(np.float64)
        t = np.arange(T) / fs

        win = int(round(BASELINE_WIN_S * fs))
        win += (win + 1) % 2                                          # force odd
        win = min(win, T if T % 2 else T - 1)

        dff = medfilt(dff_rolling_pctl(Fk, win), kernel_size=(1, MEDFILT_KERNEL))
        dffd = medfilt(dff_detrend_static(Fk, t), kernel_size=(1, MEDFILT_KERNEL))

        dff_list.append(dff)
        dffd_list.append(dffd)

        for order, (_, r) in enumerate(pa.iterrows()):
            meta_rows.append(dict(
                neuron_uid=f"z{z}_c{int(r.cell_id)}",
                z=z, cell_id=int(r.cell_id),
                region_id=int(r.region_id), region_name=r.region_name,
                purity=float(r.purity), n_pixels=int(r.n_pixels),
                centroid_y=r.cy, centroid_x=r.cx,
                plane_plot_order=order, fs_hz=fs,
            ))

        plane_cache.append((z, pa, stat, fs, T))
        print(f"z={z}: {n_cells} cells -> {n_keep} kept (win={win} fr / {win/fs:.1f} s)")

    meta = pd.DataFrame(meta_rows)
    dff_all = np.vstack(dff_list)                                    # (N, T)
    dffd_all = np.vstack(dffd_list)
    N, T = dff_all.shape
    assert len(meta) == N

    # ---- z-scoring across the whole pooled dataset -------------------------------
    mu = dff_all.mean(axis=1, keepdims=True)
    sd = dff_all.std(axis=1, keepdims=True)
    zscore = (dff_all - mu) / np.where(sd > 0, sd, 1.0)
    resid = dff_all - mu
    global_sd = float(resid.std())
    zscore_global = resid / global_sd
    print(f"\npooled: {N} neurons x {T} frames; per-neuron z-score + "
          f"global-scaled (pooled SD = {global_sd:.3f} % dF/F)")

    # ---- colours: suite2p-style random palette, independent within each region ----
    color_rgb = np.zeros((N, 3))
    for (zz, rid), grp in meta.groupby(["z", "region_id"], sort=False):
        idx = grp.index.to_numpy()
        color_rgb[idx] = region_random_colors(len(idx), seed=int(zz) * 100 + int(rid))
    meta["color_hex"] = [matplotlib.colors.to_hex(c) for c in color_rgb]

    time_s = np.arange(T) / float(meta.fs_hz.iloc[0])
    meta.insert(0, "row_index", np.arange(N))
    meta.to_csv(os.path.join(OUT_DIR, "traces_metadata.csv"), index=False)
    np.savez_compressed(
        os.path.join(OUT_DIR, "traces.npz"),
        time_s=time_s,
        neuron_uid=meta.neuron_uid.to_numpy(),
        z=meta.z.to_numpy(),
        region_id=meta.region_id.to_numpy(),
        region_name=meta.region_name.to_numpy(),
        dff_pct=dff_all.astype(np.float32),
        dff_detrend_pct=dffd_all.astype(np.float32),
        zscore=zscore.astype(np.float32),
        zscore_global=zscore_global.astype(np.float32),
    )
    print(f"wrote traces_metadata.csv  +  traces.npz")

    # ------------------------------------------------------------------------------
    # pass 2: per-plane figures
    # ------------------------------------------------------------------------------
    for z, pa, stat, fs, T in plane_cache:
        m = (meta.z == z).to_numpy()
        rows = np.where(m)[0]
        n = len(rows)
        if n == 0:
            continue
        d = dff_all[rows]                                            # (n, T) % dF/F
        cols = meta.color_hex.to_numpy()[rows]
        t = np.arange(T) / fs

        # vertical spacing: ~ one "large transient" (per-trace 96th pct) per lane
        step = float(np.clip(np.median([np.percentile(x, 96) for x in d]), 15.0, 120.0))
        step = round(step / 5.0) * 5.0

        fig_h = float(np.clip(2.6 + 0.10 * n, 4.0, 48.0))
        fig, ax = plt.subplots(figsize=(11.5, fig_h))

        # y-range: cover the full vertical extent the traces occupy (incl. transient
        # peaks/troughs), then draw the top & bottom region bands right out to those
        # edges so no trace peak ever sits on bare white
        y_lo_data = min(i * step + d[i].min() for i in range(n))
        y_hi_data = max(i * step + d[i].max() for i in range(n))
        y_lim_lo = min(-1.2 * step, y_lo_data - 0.35 * step)
        y_lim_hi = max((n - 1) * step + 2.2 * step, y_hi_data + 0.35 * step)

        rid = meta.region_id.to_numpy()[rows]
        rname = meta.region_name.to_numpy()[rows]
        blocks = []                                          # (region_id, y0, y1, y_mid)
        uniq_r = np.unique(rid)                              # rows are region-sorted -> contiguous
        for bi, r in enumerate(uniq_r):
            sel = np.where(rid == r)[0]
            y0 = y_lim_lo if bi == 0 else (sel.min() - 0.5) * step
            y1 = y_lim_hi if bi == len(uniq_r) - 1 else (sel.max() + 0.5) * step
            ax.axhspan(y0, y1, facecolor=colors.get(int(r), "0.5"), alpha=BAND_ALPHA,
                       lw=0, zorder=0)
            if sel.min() != 0:
                ax.axhline((sel.min() - 0.5) * step, color="0.62", lw=1.0, zorder=1)
            blocks.append((int(r), y0, y1, sel.mean() * step))

        for i, (tr, c) in enumerate(zip(d, cols)):
            ax.plot(t, tr + i * step, color=c, lw=0.55, antialiased=True, zorder=2)

        # right-margin region titles; scale bar slotted into the largest gap between them
        label_x = t[-1] * 1.015
        for r, _, _, y_mid in blocks:
            nm = rname[np.where(rid == r)[0][0]]
            lc = tuple(0.62 * v for v in colors.get(int(r), (0, 0, 0)))  # darken for legibility
            ax.text(label_x, y_mid, nm, va="center", ha="left", fontsize=8,
                    fontweight="bold", color=lc, zorder=4)
        mids = np.sort(np.array([b[3] for b in blocks]))
        if len(mids) >= 2:
            k = int(np.argmax(np.diff(mids)))
            gap_c = 0.5 * (mids[k] + mids[k + 1])
        else:
            gap_c = 0.5 * (n - 1) * step
        bx, cap = t[-1] * 1.135, t[-1] * 0.012
        yb0, yb1 = gap_c - SCALE_BAR_PCT / 2, gap_c + SCALE_BAR_PCT / 2
        ax.plot([bx, bx], [yb0, yb1], color="k", lw=3.5, clip_on=False, solid_capstyle="butt")
        ax.plot([bx - cap, bx + cap], [yb0, yb0], color="k", lw=2.5, clip_on=False)
        ax.plot([bx - cap, bx + cap], [yb1, yb1], color="k", lw=2.5, clip_on=False)
        ax.text(bx + cap * 1.8, gap_c, f"{SCALE_BAR_PCT:g}% ΔF/F", rotation=90,
                va="center", ha="left", fontsize=10, clip_on=False)

        ax.set_xlim(0, t[-1])
        ax.set_ylim(y_lim_lo, y_lim_hi)
        ax.set_xlabel("Time (s)", fontsize=11)
        ax.set_yticks([])
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.set_title(f"z = {z} µm    ·    {n} cells    ·    ΔF/F  "
                     f"(rolling p{BASELINE_PCTL:g} / {BASELINE_WIN_S:g}s baseline, medfilt "
                     f"{MEDFILT_KERNEL}, raw F)", fontsize=10)
        fig.text(0.5, 0.004, f"cells sorted by region then centroid · rows offset "
                 f"{step:g}% ΔF/F · band = region colour · neuron colour = "
                 f"random within region (matches circle map)",
                 ha="center", fontsize=8, color="0.35")
        fig.tight_layout(rect=(0, 0.012, 1, 1))
        fig.savefig(os.path.join(TRACE_DIR, f"traces_z{z}.png"), dpi=140, bbox_inches="tight")
        plt.close(fig)

        # ---- circle map -------------------------------------------------------
        mean_img = tifffile.imread(os.path.join(MEAN_IMG_DIR, f"meanImg_z{z}.tif"))
        label_img = tifffile.imread(os.path.join(LABEL_DIR, f"roi_labels_z{z}.tif"))
        gray = norm_img(mean_img)
        rgb = np.stack([gray, gray, gray], axis=-1)
        present = []
        for lbl in np.unique(label_img):
            if lbl == 0:
                continue
            mask = label_img == lbl
            col = np.array(colors.get(int(lbl), (1.0, 1.0, 1.0)))
            rgb[mask] = (1 - MAP_REGION_FILL_ALPHA) * rgb[mask] + MAP_REGION_FILL_ALPHA * col
            edge = mask & ~binary_erosion(mask)
            rgb[edge] = col
            present.append(int(lbl))

        fig, ax = plt.subplots(figsize=(9, 7))
        ax.imshow(rgb, origin="upper", interpolation="nearest")
        cy = meta.centroid_y.to_numpy()[rows]
        cx = meta.centroid_x.to_numpy()[rows]
        for xx, yy, c in zip(cx, cy, cols):
            ax.add_patch(Circle((xx, yy), radius=circle_radius, facecolor=c,
                                edgecolor="white", linewidth=0.5, alpha=0.97))
        ax.set_xlim(0, rgb.shape[1])
        ax.set_ylim(rgb.shape[0], 0)
        ax.axis("off")
        ax.set_title(f"z = {z} µm   ·   {n} cells   ·   "
                     f"circle r = {circle_radius:.1f} px (median ROI size)", fontsize=10)
        handles = [plt.Line2D([0], [0], color=colors[i], lw=3, label=names[i]) for i in present]
        handles.append(plt.Line2D([0], [0], marker="o", color="none", markerfacecolor="0.5",
                                  markersize=8, label="cell (colour = random within region)"))
        ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(CIRCLE_DIR, f"circlemap_z{z}.png"), dpi=140, bbox_inches="tight")
        plt.close(fig)
        print(f"z={z}: traces_z{z}.png + circlemap_z{z}.png  (step={step:.0f}%)")

    print("\nDONE")


if __name__ == "__main__":
    main()
