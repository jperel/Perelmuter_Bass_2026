#!/usr/bin/env python3
r"""
per_plane_figs.py
==================
Stage 5 of 5 (calcium / region analysis). Run after run_calcium_region_analysis.py
(needs traces.npz / traces_metadata.csv); does not depend on transient_rate.py.

For every one of the 32 included native planes (z = 2480-2790 um, 10 um steps),
produces two publication-quality figures showing all of that plane's retained cells:

  traces_z<depth>.tif / .pdf     stacked % dF/F traces for every retained cell in
                                   the plane, sorted by region then position, with
                                   region-colored background bands and a scale bar.
  circlemap_z<depth>.tif / .pdf   mean image with pulled-back region boundaries
                                   (smoothed vector contours) and one dot per
                                   retained cell, colored to match its trace.

This generalizes the single hand-picked representative-plane figure that an earlier
version of this script produced (for z = 2620 only, with a curated 24-cell subset
enlarged/numbered for print layout) to every plane, without that per-plane curation:
here every retained cell is plotted, in each plane's native (unrotated) orientation.

Traces are PCHIP-upsampled and lightly Gaussian-smoothed for display only; the
underlying sample values are unchanged. The dense trace + band layer is rasterized
at 1600 dpi in the PDF (text and axes remain vector) because many thin overlapping
vector lines render broken in PDF viewers. Region boundaries in the circle map are
drawn as smoothed vector curves (skimage find_contours -> Gaussian-smoothed polyline)
so they stay crisp in both the PDF and the TIFF.

Reads (relative to DATA_ROOT):
    outputs/traces_metadata.csv, outputs/traces.npz
    suite2p_output/02F_2min_<z>/suite2p/plane0/stat.npy, iscell.npy
        (all 32 planes, to compute the median area-equivalent ROI radius
        used for the circle map)
    roi_labels/roi_labels_z<z>.tif (all 32 planes)
    mean_images/meanImg_z<z>.tif (all 32 planes)
    labels.txt

Writes (relative to DATA_ROOT/outputs/per_plane_figs/):
    traces_z<depth>.tif/.pdf, circlemap_z<depth>.tif/.pdf  (32 planes each)

Note on data layout: `data_for_upload/` is a flat folder (Zenodo does not
preserve directory structure on upload), except for `mean_images/`,
`roi_labels/`, and `suite2p_output/`, which come from unzipping the three
archives shipped in the data package -- see the top-level README. `outputs/`
is a subfolder run_calcium_region_analysis.py creates locally for generated
results (not part of the Zenodo data package).

Requires the antspy environment: numpy, scipy, pandas, matplotlib, tifffile,
scikit-image.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = ["Arial", "DejaVu Sans"]
matplotlib.rcParams["pdf.fonttype"] = 42
import matplotlib.pyplot as plt
import tifffile
from matplotlib.patches import Circle
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import PchipInterpolator

DATA_ROOT = os.environ.get("PIPELINE_DATA_ROOT") or os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data_for_upload"))

OUT_DIR = os.path.join(DATA_ROOT, "outputs")
FIG_DIR = os.path.join(OUT_DIR, "per_plane_figs")
LABELS_TXT = os.path.join(DATA_ROOT, "labels.txt")
LABEL_DIR = os.path.join(DATA_ROOT, "roi_labels")
MEAN_IMG_DIR = os.path.join(DATA_ROOT, "mean_images")
S2P_DIR = os.path.join(DATA_ROOT, "suite2p_output")
os.makedirs(FIG_DIR, exist_ok=True)

INCLUDED_DEPTHS = list(range(2480, 2800, 10))   # 32 planes
SCALE_BAR_PCT = 200.0
BAND_ALPHA = 0.13
MAP_FILL_ALPHA = 0.12
CM = 1.0 / 2.54
DPI_TIF = 1600
TRACE_W_CM, TRACE_H_CM = 5.3, 11.8
MAP_W_CM = 7.2


def parse_labels_txt(path):
    colors, names = {}, {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            p = line.split(None, 7)
            i = int(p[0])
            colors[i] = (int(p[1]) / 255.0, int(p[2]) / 255.0, int(p[3]) / 255.0)
            names[i] = p[7].strip('"') if len(p) > 7 else str(i)
    return colors, names


def norm_img(img, pct=99):
    img = img.astype(np.float32)
    hi = np.percentile(img[img > 0], pct) if np.any(img > 0) else img.max()
    return np.clip(img / max(hi, 1e-6), 0, 1)


def circle_radius_px():
    d = []
    for zz in INCLUDED_DEPTHS:
        sp = os.path.join(S2P_DIR, f"02F_2min_{zz}", "suite2p", "plane0")
        st = np.load(os.path.join(sp, "stat.npy"), allow_pickle=True)
        ic = np.load(os.path.join(sp, "iscell.npy"))[:, 0] == 1
        d += [2 * np.sqrt(s["npix"] / np.pi) for s, k in zip(st, ic) if k]
    return float(np.median(d)) / 2.0


def _smooth_closed(xy, sig=2.2):
    """Gaussian-smooth a contour (N,2); wrap if it closes on itself, else clamp."""
    closed = np.allclose(xy[0], xy[-1], atol=1.5)
    mode = "wrap" if closed else "nearest"
    x = gaussian_filter1d(xy[:, 1], sig, mode=mode)
    y = gaussian_filter1d(xy[:, 0], sig, mode=mode)
    if closed:
        x = np.append(x, x[0]); y = np.append(y, y[0])
    return x, y


def draw_traces(d, rid, rname, chex, fs, colors, stub):
    """d (m, T) % dF/F rows already in stacked (region-sorted) order; rid/rname/chex
    aligned to d's rows."""
    n, T = d.shape
    t = np.arange(T) / fs
    step = np.median([np.percentile(x, 96) for x in d])
    step = round(float(np.clip(step, 15.0, 400.0)) / 5.0) * 5.0

    # smooth display curve (PCHIP up-sample -> light round; no overshoot, sample values kept)
    tt = np.linspace(t[0], t[-1], (T - 1) * 8 + 1)
    dd = np.array([gaussian_filter1d(PchipInterpolator(t, row)(tt), 2.0, mode="nearest")
                   for row in d])

    fs_lab, fs_ax, lw_tr = 3.2, 3.6, 0.35
    axrect = [0.04, 0.062, 0.80, 0.915]
    lab_x, bar_x = 1.015, 1.07

    fig = plt.figure(figsize=(TRACE_W_CM * CM, TRACE_H_CM * CM))
    ax = fig.add_axes(axrect)

    y_lo = min(i * step + dd[i].min() for i in range(n))
    y_hi = max(i * step + dd[i].max() for i in range(n))
    y_lim_lo = min(-1.2 * step, y_lo - 0.35 * step)
    y_lim_hi = max((n - 1) * step + 1.6 * step, y_hi + 0.35 * step)

    # many thin overlapping vector lines render broken in PDF viewers -> rasterize the
    # trace + band layer (at TIFF dpi) while keeping text/axes vector.
    rast = n > 40

    uniq = list(dict.fromkeys(rid))                     # region order as given (asc)
    blocks = []
    for bi, r in enumerate(uniq):
        sel = np.where(rid == r)[0]
        y0 = y_lim_lo if bi == 0 else (sel.min() - 0.5) * step
        y1 = y_lim_hi if bi == len(uniq) - 1 else (sel.max() + 0.5) * step
        ax.axhspan(y0, y1, facecolor=colors.get(int(r), "0.5"), alpha=BAND_ALPHA, lw=0,
                   zorder=0, rasterized=rast)
        if sel.min() != 0:
            ax.axhline((sel.min() - 0.5) * step, color="0.62", lw=0.6, zorder=1, rasterized=rast)
        blocks.append((int(r), sel.mean() * step))

    for i in range(n):
        ax.plot(tt, dd[i] + i * step, color=chex[i], lw=lw_tr, antialiased=True, zorder=2,
                rasterized=rast, solid_capstyle="round")

    # right-margin region labels + scale bar in the largest gap between them
    for r, ymid in blocks:
        nm = rname[np.where(rid == r)[0][0]]
        lc = tuple(0.6 * v for v in colors.get(int(r), (0, 0, 0)))
        ax.text(t[-1] * lab_x, ymid, nm, va="center", ha="left", fontsize=fs_lab,
                fontweight="bold", color=lc, zorder=4)
    mids = np.sort([b[1] for b in blocks])
    if len(mids) >= 2:
        k = int(np.argmax(np.diff(mids)))
        gap_c = 0.5 * (mids[k] + mids[k + 1])
    else:
        gap_c = 0.5 * (n - 1) * step
    bx = t[-1] * bar_x
    ax.plot([bx, bx], [gap_c - SCALE_BAR_PCT / 2, gap_c + SCALE_BAR_PCT / 2],
            color="k", lw=1.8, clip_on=False, solid_capstyle="butt")
    ax.text(bx + t[-1] * 0.02, gap_c, f"{SCALE_BAR_PCT:g}% ΔF/F", rotation=90,
            va="center", ha="left", fontsize=fs_lab, clip_on=False)

    ax.set_xlim(0.0, t[-1])
    ax.set_ylim(y_lim_lo, y_lim_hi)
    ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.5)
    ax.tick_params(width=0.5, length=1.8, pad=1.5, labelsize=fs_ax)
    ax.set_xlabel("Time (s)", fontsize=fs_ax + 0.5)
    ax.set_xticks([0, 30, 60, 90, 120])

    fig.savefig(os.path.join(FIG_DIR, stub + ".tif"), dpi=DPI_TIF,
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(os.path.join(FIG_DIR, stub + ".pdf"), dpi=DPI_TIF)   # dpi -> rasterized layer only
    plt.close(fig)
    print("  wrote", stub + ".tif / .pdf", f"(n={n}, step {step:g}%, "
          f"traces {'rasterized' if rast else 'vector'} in pdf)")


def draw_circlemap(zm, colors, crad, z):
    from skimage import measure
    mean_img = tifffile.imread(os.path.join(MEAN_IMG_DIR, f"meanImg_z{z}.tif"))
    label_img = tifffile.imread(os.path.join(LABEL_DIR, f"roi_labels_z{z}.tif"))

    g = norm_img(mean_img)
    rgb = np.stack([g, g, g], axis=-1)
    for lbl in np.unique(label_img):
        if lbl == 0:
            continue
        mask = label_img == lbl
        rgb[mask] = (1 - MAP_FILL_ALPHA) * rgb[mask] + MAP_FILL_ALPHA * np.array(
            colors.get(int(lbl), (1.0, 1.0, 1.0)))

    H, W = rgb.shape[:2]
    h_cm = MAP_W_CM * H / W
    fig = plt.figure(figsize=(MAP_W_CM * CM, h_cm * CM))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(rgb, origin="upper", interpolation="nearest", zorder=0)

    # smooth vector region boundaries
    for lbl in np.unique(label_img):
        if lbl == 0:
            continue
        col = colors.get(int(lbl), (1, 1, 1))
        for ct in measure.find_contours((label_img == lbl).astype(float), 0.5):
            xs, ys = _smooth_closed(ct, sig=2.2)
            ax.plot(xs, ys, color=col, lw=0.7, zorder=1, solid_capstyle="round",
                    solid_joinstyle="round")

    # native (unrotated) orientation: centroid_x is the column index, centroid_y the
    # row index, matching mean_img/label_img directly -- same convention used by
    # run_calcium_region_analysis.py's own (lower-resolution) circle_maps/ output.
    cx = zm.centroid_x.to_numpy().astype(float)
    cy = zm.centroid_y.to_numpy().astype(float)
    chex = zm.color_hex.to_numpy()
    for i in range(len(zm)):
        ax.add_patch(Circle((cx[i], cy[i]), radius=crad, facecolor=chex[i],
                            edgecolor="white", linewidth=0.3, alpha=0.97, zorder=2))

    ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.axis("off")
    stub = f"circlemap_z{z}"
    fig.savefig(os.path.join(FIG_DIR, stub + ".tif"), dpi=DPI_TIF,
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(os.path.join(FIG_DIR, stub + ".pdf"))
    plt.close(fig)
    print("  wrote", stub + ".tif / .pdf", f"({MAP_W_CM:.1f} x {h_cm:.1f} cm, n={len(zm)})")


def main():
    colors, names = parse_labels_txt(LABELS_TXT)
    meta = pd.read_csv(os.path.join(OUT_DIR, "traces_metadata.csv"))
    dff = np.load(os.path.join(OUT_DIR, "traces.npz"), allow_pickle=True)["dff_pct"]
    fs = float(meta.fs_hz.iloc[0])
    crad = circle_radius_px()

    for z in INCLUDED_DEPTHS:
        zm = meta[meta.z == z].copy().reset_index(drop=True)
        if len(zm) == 0:
            print(f"z={z}: no retained cells, skipping")
            continue
        print(f"z={z}: {len(zm)} cells")

        # region-sorted like the pipeline: region_id asc, then centroid_y, centroid_x
        full = zm.sort_values(["region_id", "centroid_y", "centroid_x"]).reset_index()
        d_full = dff[full.row_index.to_numpy()].astype(np.float64)
        draw_traces(d_full, full.region_id.to_numpy(), full.region_name.to_numpy(),
                    full.color_hex.to_numpy(), fs, colors, f"traces_z{z}")

        draw_circlemap(zm, colors, crad, z)


if __name__ == "__main__":
    main()
