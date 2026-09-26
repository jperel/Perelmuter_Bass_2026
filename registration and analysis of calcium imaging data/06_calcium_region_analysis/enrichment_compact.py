#!/usr/bin/env python3
r"""
enrichment_compact.py
======================
Stage 6 of 6 (calcium / region analysis) -- produces manuscript Figure 7
panel K. Run after run_calcium_region_analysis.py and transient_rate.py
(both traces.npz and transient_rate.csv are read).

Top-decile (>= 90th percentile) enrichment dumbbell plot: for each region, a
filled dot = that region's enrichment for one metric, an open dot = its
enrichment for mean ΔF/F, connected by a line; a dashed line at 0.10 marks
chance. Dot area is proportional to the region's total cell count. Regions
are ranked by mean-ΔF/F enrichment. "V - unspecified" is excluded
throughout.

Two variants share all data loading and plotting geometry:
  dumbbell()             mean ΔF/F + skew enrichment.
                          -> fig_enrichment_dumbbell.tif/.pdf
                          (not used in the manuscript)
  dumbbell_transients()   mean ΔF/F + transient/event-rate enrichment.
                          -> fig_enrichment_dumbbell_transients.tif/.pdf
                          (Figure 7K)
A third function, strip(), produces a 2-row heat-map alternative
(fig_enrichment_strip.tif/.pdf), also not used in the manuscript. All three
are kept in one script since dumbbell() and dumbbell_transients() share
nearly all of their plotting code, and strip() shares the same data loading.

Reads (relative to DATA_ROOT):
    06_calcium_region_analysis/outputs/traces.npz
    06_calcium_region_analysis/outputs/transient_rate.csv (from
        transient_rate.py; only dumbbell_transients() needs it)
    04_annotation/labels.txt

Writes (relative to DATA_ROOT/06_calcium_region_analysis/outputs/):
    fig_enrichment_dumbbell.tif/.pdf,
    fig_enrichment_dumbbell_transients.tif/.pdf (Figure 7K),
    fig_enrichment_strip.tif/.pdf

Requires the antspy environment: numpy, scipy, pandas, matplotlib.
"""
import os
import colorsys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = ["Arial", "DejaVu Sans"]
matplotlib.rcParams["pdf.fonttype"] = 42
import matplotlib.pyplot as plt
from scipy.stats import skew as sp_skew

DATA_ROOT = os.environ.get("PIPELINE_DATA_ROOT") or os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data_for_upload"))

OUT_DIR = os.path.join(DATA_ROOT, "06_calcium_region_analysis", "outputs")
LABELS_TXT = os.path.join(DATA_ROOT, "04_annotation", "labels.txt")
os.makedirs(OUT_DIR, exist_ok=True)

DROP_REGIONS = {"V - unspecified"}
FONT_PT, DPI, CM = 5, 1200, 1.0 / 2.54
LW_SPINE, LW_TICK, TICK_LEN, TICK_PAD = 0.5, 0.5, 1.8, 1.5
S_MIN, S_MAX = 3.5, 42.0            # dot area (pt^2) for the smallest / largest region
KEY_N = [100, 500, 2000]           # reference counts in the size key


def parse_labels_txt(path):
    colors = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            p = line.split(None, 7)
            nm = p[7].strip('"') if len(p) > 7 else p[0]
            colors[nm] = (int(p[1]) / 255.0, int(p[2]) / 255.0, int(p[3]) / 255.0)
    return colors


def darken(rgb):
    h, s, v = colorsys.rgb_to_hsv(*rgb)
    return colorsys.hsv_to_rgb(h, min(1.0, s * 1.15 + 0.05), min(v, 0.62))


def load():
    npz = np.load(os.path.join(OUT_DIR, "traces.npz"), allow_pickle=True)
    dff = npz["dff_pct"].astype(np.float64)
    reg = npz["region_name"].astype(str)
    keep = ~np.isin(reg, list(DROP_REGIONS))
    dff, reg = dff[keep], reg[keep]
    d = pd.DataFrame(dict(region=reg, mean_dff=dff.mean(1), skew=sp_skew(dff, axis=1)))
    out = {"n": d.groupby("region").size()}
    for col in ("mean_dff", "skew"):
        thr = np.percentile(d[col], 90)
        out[col] = d.assign(hi=(d[col] >= thr)).groupby("region")["hi"].mean()
    return pd.DataFrame(out)


def load_transients():
    """Same as load(), but mean ΔF/F + transient/event rate (needs transient_rate.csv)
    instead of mean ΔF/F + skew."""
    npz = np.load(os.path.join(OUT_DIR, "traces.npz"), allow_pickle=True)
    dff = npz["dff_pct"].astype(np.float64)
    reg = npz["region_name"].astype(str)
    keep = ~np.isin(reg, list(DROP_REGIONS))
    row_index = np.where(keep)[0]
    dff, reg = dff[keep], reg[keep]
    d = pd.DataFrame(dict(row_index=row_index, region=reg, mean_dff=dff.mean(1)))
    tr = pd.read_csv(os.path.join(OUT_DIR, "transient_rate.csv"))[
        ["row_index", "rate_pos_per_min"]].rename(columns={"rate_pos_per_min": "transients"})
    d = d.merge(tr, on="row_index", how="left")
    out = {"n": d.groupby("region").size()}
    for col in ("mean_dff", "transients"):
        thr = np.percentile(d[col], 90)
        out[col] = d.assign(hi=(d[col] >= thr)).groupby("region")["hi"].mean()
    return pd.DataFrame(out)


def save(fig, stub):
    fig.savefig(os.path.join(OUT_DIR, stub + ".tif"), dpi=DPI, pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(os.path.join(OUT_DIR, stub + ".pdf"))
    plt.close(fig)
    print("wrote", stub + ".tif / .pdf")


def _area(n, n_max):
    return S_MIN + (np.asarray(n) / n_max) * (S_MAX - S_MIN)     # area ∝ n


def dumbbell(E, colors):
    order = E.sort_values("mean_dff", ascending=False)
    regs = list(order.index)
    n_max = float(E["n"].max())
    dcol = {r: darken(colors.get(r, (.4, .4, .4))) for r in regs}

    # PLOT panel = first 5 cm of the figure; KEY sits to the right (crop at 5 cm).
    fig = plt.figure(figsize=(7.3 * CM, 3.5 * CM))
    ax = fig.add_axes([0.095, 0.17, 0.56, 0.80])          # right edge ~ 0.655 -> ~4.8 cm
    kax = fig.add_axes([0.695, 0.05, 0.30, 0.92]); kax.axis("off")

    for i, r in enumerate(regs):
        dff_e, skw_e = order.loc[r, "mean_dff"], order.loc[r, "skew"]
        s = _area(order.loc[r, "n"], n_max)
        c = dcol[r]
        ax.plot([min(dff_e, skw_e), max(dff_e, skw_e)], [i, i], color="0.6", lw=0.9,
                zorder=2, solid_capstyle="round")
        ax.scatter([skw_e], [i], s=s, facecolor=c, edgecolor=c, linewidths=0.8, zorder=3)  # skew
        ax.scatter([dff_e], [i], s=s, facecolor="white", edgecolor=c, linewidths=0.8,
                   zorder=4)                                                            # mean ΔF/F

    ax.axvline(0.10, color="0.55", ls=(0, (2, 2)), lw=0.5, zorder=1)
    ax.set_ylim(len(regs) - 0.5, -0.5)
    ax.set_yticks(range(len(regs)))
    ax.set_yticklabels(regs, fontsize=FONT_PT)
    for tick, r in zip(ax.get_yticklabels(), regs):
        tick.set_color(dcol[r])
        tick.set_fontweight("bold")
    ax.set_xlim(-0.006, 0.27)
    ax.set_xticks([0, 0.1, 0.2])
    ax.set_xticklabels(["0", "0.1", "0.2"], fontsize=FONT_PT)
    ax.set_xlabel("fraction of cells in top 10 %", fontsize=FONT_PT, labelpad=1.5)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_linewidth(LW_SPINE)
    ax.tick_params(width=LW_TICK, length=TICK_LEN, pad=TICK_PAD)

    # ---- key (right of the 5 cm crop line) ---------------------------------
    kax.set_xlim(0, 1); kax.set_ylim(0, 1)
    kax.scatter([0.10], [0.95], s=22, facecolor="0.3", edgecolor="0.3", linewidths=0.7)
    kax.text(0.24, 0.95, "skew", va="center", fontsize=FONT_PT)
    kax.scatter([0.10], [0.85], s=22, facecolor="white", edgecolor="0.3", linewidths=0.7)
    kax.text(0.24, 0.85, "mean ΔF/F", va="center", fontsize=FONT_PT)
    kax.plot([0.02, 0.18], [0.75, 0.75], color="0.55", ls=(0, (2, 2)), lw=0.5)
    kax.text(0.24, 0.75, "chance (0.10)", va="center", fontsize=FONT_PT)

    kax.text(0.0, 0.58, "cells / region", fontsize=FONT_PT, fontweight="bold")
    ys = [0.44, 0.29, 0.12]
    kax.scatter([0.14] * 3, ys, s=_area(KEY_N, n_max), facecolor="0.4", edgecolor="none")
    for y, nn in zip(ys, KEY_N):
        kax.text(0.36, y, f"{nn:,}", va="center", fontsize=FONT_PT)

    save(fig, "fig_enrichment_dumbbell")


def dumbbell_transients(E, colors):
    """Same figure as dumbbell(), with the filled dot = transient/event-rate top-decile
    enrichment instead of skew. Identical geometry/font/dot-size scheme; new filename."""
    order = E.sort_values("mean_dff", ascending=False)
    regs = list(order.index)
    n_max = float(E["n"].max())
    dcol = {r: darken(colors.get(r, (.4, .4, .4))) for r in regs}

    fig = plt.figure(figsize=(7.3 * CM, 3.5 * CM))
    ax = fig.add_axes([0.095, 0.17, 0.56, 0.80])
    kax = fig.add_axes([0.695, 0.05, 0.30, 0.92]); kax.axis("off")

    for i, r in enumerate(regs):
        dff_e, trn_e = order.loc[r, "mean_dff"], order.loc[r, "transients"]
        s = _area(order.loc[r, "n"], n_max)
        c = dcol[r]
        ax.plot([min(dff_e, trn_e), max(dff_e, trn_e)], [i, i], color="0.6", lw=0.9,
                zorder=2, solid_capstyle="round")
        ax.scatter([trn_e], [i], s=s, facecolor=c, edgecolor=c, linewidths=0.8, zorder=3)  # transients
        ax.scatter([dff_e], [i], s=s, facecolor="white", edgecolor=c, linewidths=0.8,
                   zorder=4)                                                            # mean ΔF/F

    ax.axvline(0.10, color="0.55", ls=(0, (2, 2)), lw=0.5, zorder=1)
    ax.set_ylim(len(regs) - 0.5, -0.5)
    ax.set_yticks(range(len(regs)))
    ax.set_yticklabels(regs, fontsize=FONT_PT)
    for tick, r in zip(ax.get_yticklabels(), regs):
        tick.set_color(dcol[r])
        tick.set_fontweight("bold")
    # transient-rate enrichment can exceed the fixed 0.27 x-axis cap used by dumbbell()
    # (observed up to ~0.30 in one region here), so the axis limit is computed from the
    # data instead of hardcoded.
    xmax = max(0.27, float(order[["mean_dff", "transients"]].to_numpy().max()) * 1.10)
    ax.set_xlim(-0.006, xmax)
    xticks = [0, 0.1, 0.2] + ([0.3] if xmax > 0.3 else [])
    ax.set_xticks(xticks)
    ax.set_xticklabels([f"{t:g}" for t in xticks], fontsize=FONT_PT)
    ax.set_xlabel("fraction of cells in top 10 %", fontsize=FONT_PT, labelpad=1.5)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_linewidth(LW_SPINE)
    ax.tick_params(width=LW_TICK, length=TICK_LEN, pad=TICK_PAD)

    # ---- key (right of the 5 cm crop line) ---------------------------------
    kax.set_xlim(0, 1); kax.set_ylim(0, 1)
    kax.scatter([0.10], [0.95], s=22, facecolor="0.3", edgecolor="0.3", linewidths=0.7)
    kax.text(0.24, 0.95, "transient rate", va="center", fontsize=FONT_PT)
    kax.scatter([0.10], [0.85], s=22, facecolor="white", edgecolor="0.3", linewidths=0.7)
    kax.text(0.24, 0.85, "mean ΔF/F", va="center", fontsize=FONT_PT)
    kax.plot([0.02, 0.18], [0.75, 0.75], color="0.55", ls=(0, (2, 2)), lw=0.5)
    kax.text(0.24, 0.75, "chance (0.10)", va="center", fontsize=FONT_PT)

    kax.text(0.0, 0.58, "cells / region", fontsize=FONT_PT, fontweight="bold")
    ys = [0.44, 0.29, 0.12]
    kax.scatter([0.14] * 3, ys, s=_area(KEY_N, n_max), facecolor="0.4", edgecolor="none")
    for y, nn in zip(ys, KEY_N):
        kax.text(0.36, y, f"{nn:,}", va="center", fontsize=FONT_PT)

    save(fig, "fig_enrichment_dumbbell_transients")


def strip(E, colors):
    order = E.sort_values("mean_dff", ascending=False)
    regs = list(order.index)
    M = order[["mean_dff", "skew"]].to_numpy().T
    fig = plt.figure(figsize=(3.85 * CM, 1.25 * CM))
    ax = fig.add_axes([0.135, 0.42, 0.85, 0.45])
    ax.imshow(M, aspect="auto", cmap="RdBu_r", vmin=0.10 - 0.13, vmax=0.10 + 0.13)
    ax.set_xticks(range(len(regs)))
    ax.set_xticklabels(regs, fontsize=FONT_PT, rotation=90, va="top")
    ax.set_yticks([0, 1]); ax.set_yticklabels(["ΔF/F", "skew"], fontsize=FONT_PT)
    ax.tick_params(width=LW_TICK, length=TICK_LEN, pad=1.0)
    for sp in ax.spines.values():
        sp.set_linewidth(LW_SPINE)
    save(fig, "fig_enrichment_strip")


def main():
    colors = parse_labels_txt(LABELS_TXT)
    E = load()
    print(E.round(3).to_string())
    dumbbell(E, colors)
    strip(E, colors)

    Et = load_transients()
    print(Et.round(3).to_string())
    dumbbell_transients(Et, colors)


if __name__ == "__main__":
    main()
