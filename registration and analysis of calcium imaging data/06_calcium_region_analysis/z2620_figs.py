#!/usr/bin/env python3
r"""
z2620_figs.py
=============
Stage 6 of 6 (calcium / region analysis) -- produces manuscript Figure 7
panels G and H. Run after run_calcium_region_analysis.py (needs traces.npz /
traces_metadata.csv); does not depend on transient_rate.py.

Illustrator-ready figures for plane z = 2620 (199 kept cells, 6 regions plus
one Vs cell):

  circlemap_z2620.tif / .pdf   circle map (mean image + region boundaries +
                                one dot per retained cell), rotated 90
                                degrees counter-clockwise, with only the 24
                                selected cells enlarged/black-edged/numbered
                                1-24 matching z2620_rep24's numbering.
                                Figure 7G.
  z2620_full.tif / .pdf        all 199 traces, exact 5.3 x 11.8 cm, 1600 dpi.
                                The PDF is fully vector -> no pixelation on
                                zoom; use it for print. Figure 7H.
  z2620_rep24.tif / .pdf       24 representative traces (4 per region, 6
                                regions, Vs excluded), numbered 1-24,
                                ~7 x 9 cm. Figure 7H.

Representative-cell selection (per region): standardize peak ΔF/F (95th
percentile) and skew about the region median, then take the 4 cells nearest
4 targets in that (peak_z, skew_z) plane -- two near the median, one
high-peak, one high-skew -- so the selection spans both activity parameters
while staying centered on typical behavior rather than picking tail
outliers.

Circle-map label placement uses force-directed decluttering (labels repel
each other, the cell markers, and the image border, and spring toward a
fixed radius from their own anchor) followed by a swap pass that exchanges
any two label positions whose leader lines cross, to keep the leader lines
readable. Region boundaries are drawn as smoothed vector curves (skimage
find_contours -> Gaussian-smoothed polyline) so they stay crisp in both the
PDF and the TIFF.

Traces are PCHIP-upsampled and lightly Gaussian-smoothed for display only;
the underlying sample values are unchanged. The dense trace + band layer is
rasterized at 1600 dpi in the PDF (text and axes remain vector) because many
thin overlapping vector lines render broken in PDF viewers; figures with
fewer traces stay fully vector.

Reads (relative to DATA_ROOT):
    06_calcium_region_analysis/outputs/traces_metadata.csv, traces.npz
    06_calcium_region_analysis/suite2p_output/02F_2min_<z>/suite2p/plane0/
        stat.npy, iscell.npy (all 32 planes, to recompute the median
        area-equivalent ROI radius used for the circle map)
    05_roi_pullback/roi_labels/roi_labels_z2620.tif
    01_suite2p_preprocessing/mean_images/meanImg_z2620.tif
    04_annotation/labels.txt

Writes (relative to DATA_ROOT/06_calcium_region_analysis/outputs/):
    z2620_full.tif/.pdf, z2620_rep24.tif/.pdf, circlemap_z2620.tif/.pdf

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
import matplotlib.patheffects as pe
import tifffile
import matplotlib.patches as mpatches
from matplotlib.patches import Circle
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import PchipInterpolator
from scipy.stats import skew as sp_skew

DATA_ROOT = os.environ.get("PIPELINE_DATA_ROOT") or os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data_for_upload"))

OUT_DIR = os.path.join(DATA_ROOT, "06_calcium_region_analysis", "outputs")
LABELS_TXT = os.path.join(DATA_ROOT, "04_annotation", "labels.txt")
LABEL_DIR = os.path.join(DATA_ROOT, "05_roi_pullback", "roi_labels")
MEAN_IMG_DIR = os.path.join(DATA_ROOT, "01_suite2p_preprocessing", "mean_images")
S2P_DIR = os.path.join(DATA_ROOT, "06_calcium_region_analysis", "suite2p_output")
os.makedirs(OUT_DIR, exist_ok=True)

Z = 2620
REP_REGION_IDS = [1, 2, 3, 5, 6, 8]          # Dm-r, Dl-r, Dl-c, Dp, Dm-c, Vd  (Vs=10 excluded)
SCALE_BAR_PCT = 200.0
BAND_ALPHA = 0.13
MAP_FILL_ALPHA = 0.12
CM = 1.0 / 2.54
DPI_TIF = 1600


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
    for zz in range(2480, 2800, 10):
        sp = os.path.join(S2P_DIR, f"02F_2min_{zz}", "suite2p", "plane0")
        st = np.load(os.path.join(sp, "stat.npy"), allow_pickle=True)
        ic = np.load(os.path.join(sp, "iscell.npy"))[:, 0] == 1
        d += [2 * np.sqrt(s["npix"] / np.pi) for s, k in zip(st, ic) if k]
    return float(np.median(d)) / 2.0


# ----------------------------------------------------------------------------------
def draw_traces(d, rid, rname, chex, fs, colors, w_cm, h_cm, stub, numbers=None,
                small=False):
    """d (m, T) % ΔF/F rows already in stacked (region-sorted) order; rid/rname/chex
    aligned; optional `numbers` (len m) drawn in a left gutter."""
    n, T = d.shape
    t = np.arange(T) / fs
    if numbers is None:
        step = np.median([np.percentile(x, 96) for x in d])
    else:
        # rep figure: space rows by ~1 typical full trace-swing so peaks stay ~in lane
        # and the vertical scale sits closer to the full plot (less visually "noisy")
        step = 1.35 * np.median([np.percentile(x, 99) - np.percentile(x, 1) for x in d])
    step = round(float(np.clip(step, 15.0, 400.0)) / 5.0) * 5.0

    # smooth display curve (PCHIP up-sample -> light round; no overshoot, sample values kept)
    tt = np.linspace(t[0], t[-1], (T - 1) * 8 + 1)
    dd = np.array([gaussian_filter1d(PchipInterpolator(t, row)(tt), 2.0, mode="nearest")
                   for row in d])

    fs_lab = 3.2 if small else 6.0
    fs_ax = 3.6 if small else 6.5
    lw_tr = 0.35 if small else 0.7
    axrect = [0.04, 0.062, 0.80, 0.915] if numbers is None else [0.055, 0.075, 0.755, 0.90]
    lab_x, bar_x = (1.015, 1.07) if numbers is None else (1.02, 1.13)
    x_left = 0.0 if numbers is None else -0.10 * t[-1]

    fig = plt.figure(figsize=(w_cm * CM, h_cm * CM))
    ax = fig.add_axes(axrect)

    y_lo = min(i * step + dd[i].min() for i in range(n))
    y_hi = max(i * step + dd[i].max() for i in range(n))
    y_lim_lo = min(-1.2 * step, y_lo - 0.35 * step)
    y_lim_hi = max((n - 1) * step + 1.6 * step, y_hi + 0.35 * step)

    # many thin overlapping vector lines render broken in PDF viewers -> rasterise the
    # trace + band layer (at TIFF dpi) while keeping text/axes vector. few-trace figures
    # stay fully vector.
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

    gutter_x = -0.055 * t[-1]
    for i in range(n):
        ax.plot(tt, dd[i] + i * step, color=chex[i], lw=lw_tr, antialiased=True, zorder=2,
                rasterized=rast, solid_capstyle="round")
        if numbers is not None:
            y0i = i * step + dd[i][0]                   # trace value at the left edge
            ax.plot([gutter_x * 0.42, 0.0], [i * step, y0i], color=chex[i], lw=1.0,
                    zorder=4, solid_capstyle="round")
            ax.text(gutter_x, i * step, str(numbers[i]), ha="right", va="center",
                    fontsize=6.2, fontweight="bold", zorder=5,
                    path_effects=[pe.withStroke(linewidth=1.6, foreground="white")])

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

    ax.set_xlim(x_left, t[-1])
    ax.set_ylim(y_lim_lo, y_lim_hi)
    ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.5)
    ax.tick_params(width=0.5, length=1.8, pad=1.5, labelsize=fs_ax)
    ax.set_xlabel("Time (s)", fontsize=fs_ax + 0.5)
    ax.set_xticks([0, 30, 60, 90, 120])

    fig.savefig(os.path.join(OUT_DIR, stub + ".tif"), dpi=DPI_TIF,
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(os.path.join(OUT_DIR, stub + ".pdf"), dpi=DPI_TIF)   # dpi -> rasterised layer only
    plt.close(fig)
    print("wrote", stub + ".tif / .pdf", f"({w_cm} x {h_cm} cm, step {step:g}%, "
          f"traces {'rasterised' if n > 40 else 'vector'} in pdf)")


def _smooth_closed(xy, sig=2.2):
    """gaussian-smooth a contour (N,2); wrap if it closes on itself, else clamp."""
    closed = np.allclose(xy[0], xy[-1], atol=1.5)
    mode = "wrap" if closed else "nearest"
    x = gaussian_filter1d(xy[:, 1], sig, mode=mode)
    y = gaussian_filter1d(xy[:, 0], sig, mode=mode)
    if closed:
        x = np.append(x, x[0]); y = np.append(y, y[0])
    return x, y


def _declutter(anchors, r_lab, sep, bounds, iters=1200):
    """force-directed label placement: repel labels from each other + from every circle
    + from the (W,H) image borders, spring each label toward radius `r_lab` off its own
    anchor. `bounds` = (W, H); labels are kept `pad` inside. anchors: (m,2) px."""
    W, H = bounds
    pad = 0.55 * sep
    rng = np.random.default_rng(0)
    ang = rng.uniform(0, 2 * np.pi, len(anchors))
    pos = anchors + r_lab * np.c_[np.cos(ang), np.sin(ang)]
    for _ in range(iters):
        disp = np.zeros_like(pos)
        for i in range(len(pos)):
            dl = pos[i] - pos                      # label-label
            dd = np.hypot(dl[:, 0], dl[:, 1]) + 1e-6
            near = (dd < sep) & (dd > 0)
            disp[i] += ((dl[near] / dd[near, None]) * (sep - dd[near, None]) * 0.7).sum(0)
            da = pos[i] - anchors                  # label vs all circles
            dda = np.hypot(da[:, 0], da[:, 1]) + 1e-6
            nearc = dda < r_lab
            disp[i] += ((da[nearc] / dda[nearc, None]) * (r_lab - dda[nearc, None]) * 0.5).sum(0)
            v = pos[i] - anchors[i]                # spring to ideal radius off own anchor
            vd = np.hypot(*v) + 1e-6
            disp[i] += (anchors[i] + v / vd * r_lab - pos[i]) * 0.035
            # border repulsion
            disp[i, 0] += max(0.0, pad - pos[i, 0]) - max(0.0, pos[i, 0] - (W - pad))
            disp[i, 1] += max(0.0, pad - pos[i, 1]) - max(0.0, pos[i, 1] - (H - pad))
        pos += np.clip(disp, -sep, sep)
    # hard final pass: force any still-overlapping pair apart
    for _ in range(80):
        for i in range(len(pos)):
            for j in range(i + 1, len(pos)):
                d = pos[i] - pos[j]
                dist = np.hypot(*d) + 1e-6
                if dist < sep:
                    push = (d / dist) * (sep - dist) * 0.5
                    pos[i] += push
                    pos[j] -= push
    return np.clip(pos, pad, [W - pad, H - pad])    # safety clamp inside the image


def draw_circlemap(zm, sel_rows_local, numbers, colors, names, crad, w_cm=7.2):
    from skimage import measure
    mean_img = np.rot90(tifffile.imread(os.path.join(MEAN_IMG_DIR, f"meanImg_z{Z}.tif")), 1)
    label_img = np.rot90(tifffile.imread(os.path.join(LABEL_DIR, f"roi_labels_z{Z}.tif")), 1)
    Worig = tifffile.imread(os.path.join(LABEL_DIR, f"roi_labels_z{Z}.tif")).shape[1]

    g = norm_img(mean_img)
    rgb = np.stack([g, g, g], axis=-1)
    for lbl in np.unique(label_img):
        if lbl == 0:
            continue
        mask = label_img == lbl
        rgb[mask] = (1 - MAP_FILL_ALPHA) * rgb[mask] + MAP_FILL_ALPHA * np.array(
            colors.get(int(lbl), (1.0, 1.0, 1.0)))

    H, W = rgb.shape[:2]
    h_cm = w_cm * H / W
    fig = plt.figure(figsize=(w_cm * CM, h_cm * CM))
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

    # rotate centroids: 90 deg CCW  ->  x' = y ,  y' = (Worig-1) - x
    cx = zm.centroid_y.to_numpy().astype(float)
    cy = (Worig - 1) - zm.centroid_x.to_numpy().astype(float)
    chex = zm.color_hex.to_numpy()
    sel = set(sel_rows_local)
    for i in range(len(zm)):
        s = i in sel
        ax.add_patch(Circle((cx[i], cy[i]), radius=crad * (1.7 if s else 1.0),
                            facecolor=chex[i], edgecolor="black" if s else "white",
                            linewidth=0.7 if s else 0.3, alpha=0.97, zorder=4 if s else 2))

    anchors = np.c_[cx[sel_rows_local], cy[sel_rows_local]]
    labpos = _declutter(anchors, r_lab=crad * 5.5, sep=crad * 8.5, bounds=(W, H))

    # un-cross leader lines: swap any two label positions whose segments intersect
    def _cross(a, b, c, d):
        def ccw(p, q, r):
            return (r[1] - p[1]) * (q[0] - p[0]) > (q[1] - p[1]) * (r[0] - p[0])
        return ccw(a, c, d) != ccw(b, c, d) and ccw(a, b, c) != ccw(a, b, d)
    for _ in range(30):
        done = True
        for i in range(len(labpos)):
            for j in range(i + 1, len(labpos)):
                if _cross(anchors[i], labpos[i], anchors[j], labpos[j]):
                    labpos[[i, j]] = labpos[[j, i]]
                    done = False
        if done:
            break

    for (ax_, ay_), (lx, ly), num in zip(anchors, labpos, numbers):
        conn = mpatches.FancyArrowPatch((ax_, ay_), (lx, ly), arrowstyle="-",
                                        connectionstyle="arc3,rad=0.12", color="white",
                                        lw=0.7, zorder=5, capstyle="round",
                                        path_effects=[pe.withStroke(linewidth=1.7, foreground="0.25")])
        ax.add_patch(conn)
        ax.text(lx, ly, str(num), ha="center", va="center", fontsize=5.2,
                fontweight="bold", color="black", zorder=6,
                path_effects=[pe.withStroke(linewidth=1.9, foreground="white")])
    ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.axis("off")
    fig.savefig(os.path.join(OUT_DIR, "circlemap_z2620.tif"), dpi=DPI_TIF,
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(os.path.join(OUT_DIR, "circlemap_z2620.pdf"))
    plt.close(fig)
    print(f"wrote circlemap_z2620.tif / .pdf ({w_cm:.1f} x {h_cm:.1f} cm)")


def main():
    colors, names = parse_labels_txt(LABELS_TXT)
    meta = pd.read_csv(os.path.join(OUT_DIR, "traces_metadata.csv"))
    dff = np.load(os.path.join(OUT_DIR, "traces.npz"), allow_pickle=True)["dff_pct"]
    fs = float(meta.fs_hz.iloc[0])
    crad = circle_radius_px()

    zm = meta[meta.z == Z].copy().reset_index(drop=True)
    d_all = dff[zm.row_index.to_numpy()].astype(np.float64)      # (199, T) in zm order
    zm["peak"] = np.percentile(d_all, 95, axis=1)
    zm["skew"] = sp_skew(d_all, axis=1)

    # ---- full plot (region-sorted like the pipeline: region_id asc, then cy, cx) ----
    full = zm.sort_values(["region_id", "centroid_y", "centroid_x"]).reset_index()
    d_full = d_all[full["index"].to_numpy()]
    draw_traces(d_full, full.region_id.to_numpy(), full.region_name.to_numpy(),
                full.color_hex.to_numpy(), fs, colors, 5.3, 11.8, "z2620_full", small=True)

    # ---- pick 24 representative (4 / region, 6 regions) ----
    # within each region, standardise peak & skew about the region median, then take the
    # cells nearest 4 targets in that (peak_z, skew_z) plane: near-median, near-median,
    # high-peak, high-skew  -> spans both parameters, 2 of 4 sit at the median.
    TARGETS = [(-0.35, -0.35), (0.15, 0.15), (1.05, -0.10), (-0.10, 1.05)]
    picks = []
    for rid in REP_REGION_IDS:
        sub = zm[zm.region_id == rid].reset_index()      # 'index' = zm-local row
        p, s = sub["peak"].to_numpy(), sub["skew"].to_numpy()
        pz = (p - np.median(p)) / (p.std() + 1e-9)
        sz = (s - np.median(s)) / (s.std() + 1e-9)
        chosen = []
        for tx, ty in TARGETS:
            d2 = (pz - tx) ** 2 + (sz - ty) ** 2
            d2[chosen] = np.inf
            chosen.append(int(np.argmin(d2)))
        g = sub.iloc[chosen].copy()
        g["ord"] = pz[chosen] + sz[chosen]               # order within region: quiet -> active
        picks.append(g)
    pk = pd.concat(picks).sort_values(["region_id", "ord"]).reset_index(drop=True)
    pk["num"] = np.arange(1, len(pk) + 1)
    print("\nselected 24:")
    print(pk[["num", "region_name", "cell_id", "peak", "skew"]].round(2).to_string(index=False))
    print("\nper-region median (all cells):")
    print(zm.groupby("region_name")[["peak", "skew"]].median().round(2).to_string())

    loc = pk["index"].to_numpy()                       # rows into zm / d_all
    d_rep = d_all[loc]
    draw_traces(d_rep, pk.region_id.to_numpy(), pk.region_name.to_numpy(),
                pk.color_hex.to_numpy(), fs, colors, 8.0, 12.5, "z2620_rep24",
                numbers=pk.num.to_numpy())

    draw_circlemap(zm, list(loc), pk.num.to_numpy(), colors, names, crad)


if __name__ == "__main__":
    main()
