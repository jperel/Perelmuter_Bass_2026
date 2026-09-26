#!/usr/bin/env python3
r"""
manuscript_panels.py
====================
Stage 6 of 6 (calcium / region analysis) -- produces manuscript Figure 7
panels I and J. Run after run_calcium_region_analysis.py and
transient_rate.py (both traces.npz and transient_rate.csv are read).

Tiny standalone TIFF panels sized for a 17.4 x 20 cm manuscript figure
layout. Each panel is drawn at its exact final size (3.85 x 2.70 cm,
Arial 5 pt, 600 dpi, LZW) so text renders at true size once placed into the
layout; no titles or axis labels are added here (added later in layout
software), but tick marks/numbers and the y-axis region names are kept.
Regions are ranked high (top) -> low (bottom) by the plotted statistic;
"V - unspecified" is excluded throughout.

This script writes about seven output files per run (plain scatter panels,
scatter+violin panels, and enrichment-bar panels, across several metrics),
of which only two are used in the manuscript:

  fig_mean_dff_scatter_violin_A.tif   Figure 7I (mean ΔF/F per region: KDE
                                       violin + mean +/- 1 SD)
  raincloud_transients_violin.tif     Figure 7J (transient/event rate per
                                       region: KDE violin + mean +/- 1 SD)

The naming is inconsistent between these two outputs (scatter_violin vs.
raincloud_*_violin); this is carried over unchanged from the original
analysis rather than renamed, since both are produced by the same
scatter_panel_violin() function and the names simply reflect the order the
two panel types were added in. See this stage's README for the full
file-to-panel mapping. The remaining outputs (plain dot-scatter panels,
skew variants, top-decile enrichment bar panels) are not part of Figure 7
but are left in place since they share all plotting code with the two that
are.

Reads (relative to DATA_ROOT):
    06_calcium_region_analysis/outputs/traces.npz
    06_calcium_region_analysis/outputs/transient_rate.csv (from
        transient_rate.py; only the transients-based panels need it)
    04_annotation/labels.txt

Writes (relative to DATA_ROOT/06_calcium_region_analysis/outputs/):
    fig_mean_dff_scatter_A.tif, fig_mean_dff_scatter_violin_A.tif,
    fig_skew_scatter_A.tif, fig_transients_scatter_A.tif,
    raincloud_transients_violin.tif, fig_enriched_mean_dff.tif,
    fig_enriched_skew.tif

Requires the antspy environment: numpy, scipy, pandas, matplotlib.
"""
import os
import colorsys
import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = ["Arial", "DejaVu Sans"]
import matplotlib.pyplot as plt
from scipy.stats import skew as sp_skew, gaussian_kde

DATA_ROOT = os.environ.get("PIPELINE_DATA_ROOT") or os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data_for_upload"))

OUT_DIR = os.path.join(DATA_ROOT, "06_calcium_region_analysis", "outputs")
LABELS_TXT = os.path.join(DATA_ROOT, "04_annotation", "labels.txt")
os.makedirs(OUT_DIR, exist_ok=True)

DROP_REGIONS = {"V - unspecified"}
RNG = np.random.default_rng(42)

# ---- exact panel geometry / type --------------------------------------------
PANEL_W_CM = 3.85
PANEL_H_CM = 2.70
FONT_PT   = 5
DPI       = 600
N_DOTS    = 60          # dots drawn per region (subsampled)
JIT       = 0.22        # +/- y jitter for the dots (in row units)
CM = 1.0 / 2.54
STYLE     = "A"         # scatter style: A = coloured dots + black median/IQR
# region set(s): "" = all 13 (fit at 5 pt in 2.7 cm, ~5.1 pt/row)
VARIANTS = {"": 0}
# axes box as fractions of the figure (room for y-labels left, x-ticks bottom)
LEFT, RIGHT, BOTTOM, TOP = 0.190, 0.975, 0.120, 0.985

LW_SPINE, LW_TICK, TICK_LEN, TICK_PAD = 0.5, 0.5, 1.8, 1.5


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


def new_panel():
    fig = plt.figure(figsize=(PANEL_W_CM * CM, PANEL_H_CM * CM))
    ax = fig.add_axes([LEFT, BOTTOM, RIGHT - LEFT, TOP - BOTTOM])
    return fig, ax


def finish(ax, regs, xlim, xticks):
    ax.set_ylim(len(regs) - 0.5, -0.5)
    ax.set_yticks(range(len(regs)))
    ax.set_yticklabels(regs, fontsize=FONT_PT)
    ax.set_xlim(*xlim)
    ax.set_xticks(xticks)
    ax.set_xticklabels([f"{t:g}" for t in xticks], fontsize=FONT_PT)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(LW_SPINE)
    ax.tick_params(width=LW_TICK, length=TICK_LEN, pad=TICK_PAD)
    ax.set_xlabel(""); ax.set_ylabel(""); ax.set_title("")


def save(fig, name):
    fig.savefig(os.path.join(OUT_DIR, name), dpi=DPI, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print("wrote", name)


def scatter_panel(byreg, order, colors, style, xlim, xticks, name):
    fig, ax = new_panel()
    order = list(order)
    for i, r in enumerate(order):
        v = np.asarray(byreg[r])
        vd = v if len(v) <= N_DOTS else RNG.choice(v, N_DOTS, replace=False)
        yj = i + RNG.uniform(-JIT, JIT, len(vd))
        col = darken(colors.get(r, (.4, .4, .4)))
        dot_c = col if style == "A" else "#4d4d4d"
        mark_c = "black" if style == "A" else col
        ax.scatter(vd, yj, s=1.0, c=[dot_c], alpha=0.55, lw=0, zorder=2, clip_on=True)
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        ax.plot([q1, q3], [i, i], color=mark_c, lw=0.6, solid_capstyle="butt", zorder=3)
        ax.plot(med, i, "o", ms=2.0, color=mark_c, mec="none", zorder=4)
    finish(ax, order, xlim, xticks)
    save(fig, name)


def scatter_panel_violin(byreg, order, colors, xlim, xticks, name, bw=None, vh=0.36):
    """Symmetric (mirrored) violin + mean +/- 1 SD embedded on the row centreline, at the
    exact same panel geometry/font as scatter_panel -- no individual dots. Ranked by
    mean (not median) since mean is the statistic actually drawn."""
    fig, ax = new_panel()
    order = list(order)
    for i, r in enumerate(order):
        v = np.asarray(byreg[r])
        col = darken(colors.get(r, (.4, .4, .4)))
        if len(v) > 2 and np.ptp(v) > 0:
            lo, hi = max(np.percentile(v, 1), xlim[0]), min(np.percentile(v, 99), xlim[1])
            xs = np.linspace(lo, hi, 256)
            try:
                dens = gaussian_kde(v, bw_method=bw)(xs)
                dens = dens / dens.max() * vh
                m = dens > 0.04 * vh
                ax.fill_between(xs[m], i - dens[m], i + dens[m], facecolor=col, alpha=0.8,
                                lw=0.25, edgecolor=col, zorder=2, clip_on=True)
            except Exception:
                pass
        mean, sd = v.mean(), v.std()
        ax.plot([mean - sd, mean + sd], [i, i], color="black", lw=0.6,
                 solid_capstyle="butt", zorder=3)
        ax.plot(mean, i, "o", ms=2.0, color="black", mec="none", zorder=4)
    finish(ax, order, xlim, xticks)
    save(fig, name)


def bar_panel(frac, colors, xlim, xticks, chance, name):
    order = frac.sort_values(ascending=False)
    regs = list(order.index)
    fig, ax = new_panel()
    for i, r in enumerate(regs):
        ax.barh(i, order[r], height=0.68, color=darken(colors.get(r, (.4, .4, .4))),
                lw=0, zorder=2)
    ax.axvline(chance, color="0.5", ls=(0, (2, 2)), lw=0.45, zorder=3)
    finish(ax, regs, xlim, xticks)
    save(fig, name)


def main():
    import pandas as pd
    colors = parse_labels_txt(LABELS_TXT)
    npz = np.load(os.path.join(OUT_DIR, "traces.npz"), allow_pickle=True)
    dff = npz["dff_pct"].astype(np.float64)
    reg = npz["region_name"].astype(str)
    keep = ~np.isin(reg, list(DROP_REGIONS))
    row_index = np.where(keep)[0]
    dff, reg = dff[keep], reg[keep]

    d = pd.DataFrame(dict(row_index=row_index, region=reg, mean_dff=dff.mean(1),
                          skew=sp_skew(dff, axis=1),
                          peak=np.percentile(dff, 95, axis=1)))
    tr_path = os.path.join(OUT_DIR, "transient_rate.csv")
    have_tr = os.path.exists(tr_path)
    if have_tr:
        tr = pd.read_csv(tr_path)[["row_index", "rate_pos_per_min"]].rename(
            columns={"rate_pos_per_min": "transients"})
        d = d.merge(tr, on="row_index", how="left")
    n_by = d.groupby("region").size()

    for vname, min_n in VARIANTS.items():
        sfx = f"_{vname}" if vname else ""
        dv = d[d.region.map(n_by) >= min_n]
        mdf_by = {r: g["mean_dff"].to_numpy() for r, g in dv.groupby("region")}
        skw_by = {r: g["skew"].to_numpy() for r, g in dv.groupby("region")}
        mdf_order = dv.groupby("region")["mean_dff"].median().sort_values(ascending=False).index
        skw_order = dv.groupby("region")["skew"].median().sort_values(ascending=False).index

        # scatter x-axis: mean ΔF/F cropped to 70 (a few % of cells fall off, no pile);
        # skew to the 99.9th pct
        scatter_panel(mdf_by, mdf_order, colors, STYLE, (8, 70), [20, 40, 60],
                      f"fig_mean_dff_scatter_{STYLE}{sfx}.tif")
        mdf_order_mean = dv.groupby("region")["mean_dff"].mean().sort_values(ascending=False).index
        scatter_panel_violin(mdf_by, mdf_order_mean, colors, (8, 70), [20, 40, 60],
                             f"fig_mean_dff_scatter_violin_{STYLE}{sfx}.tif")
        scatter_panel(skw_by, skw_order, colors, STYLE,
                      (-0.6, np.ceil(np.percentile(dv["skew"], 99.9))), [0, 2, 4],
                      f"fig_skew_scatter_{STYLE}{sfx}.tif")
        if have_tr:
            trn_by = {r: g["transients"].to_numpy() for r, g in dv.groupby("region")}
            trn_order = dv.groupby("region")["transients"].median().sort_values(ascending=False).index
            scatter_panel(trn_by, trn_order, colors, STYLE, (-0.15, 3.2), [0, 1, 2, 3],
                          f"fig_transients_scatter_{STYLE}{sfx}.tif")
            trn_order_mean = dv.groupby("region")["transients"].mean().sort_values(ascending=False).index
            # violin fill is KDE-trimmed (density < 4% of peak dropped), so the sparse tail
            # near the true max (~2.99) never renders -- visible violins stop ~2.5, so the
            # dot-scatter panel's 3.2/tick-at-3 range is wasted space here. Cut to just past
            # where the violins actually end.
            scatter_panel_violin(trn_by, trn_order_mean, colors, (-0.15, 2.7), [0, 1, 2],
                                 "raincloud_transients_violin.tif")

        for col, xlim, xticks, tag in [
            ("mean_dff", (0, 0.26), [0, 0.1, 0.2], "mean_dff"),
            ("skew", (0, 0.22), [0, 0.1, 0.2], "skew"),
        ]:
            thr = np.percentile(dv[col], 90)
            frac = dv.assign(hi=(dv[col] >= thr)).groupby("region")["hi"].mean()
            bar_panel(frac, colors, xlim, xticks, 0.10, f"fig_enriched_{tag}{sfx}.tif")
        print(f"[{vname or 'all-13'}] {len(mdf_order)} regions: {list(mdf_order)}")

    print(f"\npanel = {PANEL_W_CM} x {PANEL_H_CM} cm, Arial {FONT_PT} pt, {DPI} dpi")


if __name__ == "__main__":
    main()
