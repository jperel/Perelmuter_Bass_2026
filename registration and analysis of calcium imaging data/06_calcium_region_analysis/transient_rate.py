#!/usr/bin/env python3
r"""
transient_rate.py
=================
Stage 6 of 6 (calcium / region analysis) -- run second in this folder, after
run_calcium_region_analysis.py and before z2620_figs.py / manuscript_panels.py
/ enrichment_compact.py. manuscript_panels.py and enrichment_compact.py both
need this script's transient_rate.csv output.

Significant-transient rate (events / min) per neuron, Dombeck et al.
2007/2010-style event detection.

Two choices this dataset requires, both departures from
run_calcium_region_analysis.py's primary dff_pct baseline:

  1. dff_pct's rolling 8th-percentile baseline clips the negative
     deflections needed to calibrate noise, so ΔF/F is recomputed here from
     raw F with a SYMMETRIC baseline instead: multiplicative linear detrend
     (bleaching correction) -> static median F0 -> 3-frame median filter.
     (A residual-drift subtraction term on top of the linear detrend was
     evaluated and dropped: photobleaching outside the deepest planes
     (z >= 2680) is already handled by the linear detrend, and the extra
     smoothing partially tracked real transients -- tau ~3.8 s events sit on
     a similar timescale -- leaving a negative post-event rebound that
     inflated the measured false-positive rate. Dropping the term left the
     detected event rate essentially unchanged while roughly halving the
     measured false-positive fraction at K = 3 sigma.)
  2. A diff-based noise sigma (e.g. median absolute frame-to-frame
     difference) is deflated by the trace's autocorrelation by roughly 4x,
     which makes a nominal "3 sigma" threshold meaningless (~50% false
     positives in practice). Instead, sigma is taken from the
     noise-dominated lower half of the ΔF/F distribution -- downward
     excursions are assumed to be essentially all noise, following Dombeck
     et al.:
        sigma = median(trace) - 15.87th percentile(trace)   (~1 sigma below center)

Per neuron: positive events = find_peaks(dff, height >= K*sigma,
prominence >= 1.5*sigma, distance >= 4 frames, width >= 1); negative events =
the same applied to -dff, giving the false-positive floor. K is swept over a
fixed grid and the smallest K with a dataset-wide false-positive fraction
<= 0.10 is selected (K = 3 sigma here, ~4% dataset-wide false-positive
fraction).

Caveat carried through the rest of this stage: this is a single-animal,
descriptive analysis (one census of ~all cells in a chosen, non-random depth
series) -- no confidence intervals or effect sizes are computed here or
anywhere else in this stage.

Reads (relative to DATA_ROOT):
    06_calcium_region_analysis/outputs/traces_metadata.csv, traces.npz
        (from run_calcium_region_analysis.py; used for neuron/plane
        bookkeeping and frame rate/duration -- ΔF/F itself is re-derived
        from raw F per plane below, not taken from traces.npz)
    06_calcium_region_analysis/suite2p_output/02F_2min_<z>/suite2p/plane0/F.npy
    04_annotation/labels.txt

Writes (relative to DATA_ROOT/06_calcium_region_analysis/outputs/):
    transient_rate.csv (per neuron, row-aligned to traces_metadata.csv)
    transient_rate.png

Requires the antspy environment: numpy, scipy, pandas, matplotlib.
"""
import os
import colorsys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import medfilt, find_peaks

DATA_ROOT = os.environ.get("PIPELINE_DATA_ROOT") or os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data_for_upload"))

OUT_DIR = os.path.join(DATA_ROOT, "06_calcium_region_analysis", "outputs")
S2P_DIR = os.path.join(DATA_ROOT, "06_calcium_region_analysis", "suite2p_output")
LABELS_TXT = os.path.join(DATA_ROOT, "04_annotation", "labels.txt")
os.makedirs(OUT_DIR, exist_ok=True)

RNG = np.random.default_rng(0)
KS = [2.5, 3.0, 3.5, 4.0]
DIST = 4           # min frames between events (~3.3 s ~ one GCaMP6s decay)
FP_TARGET = 0.10
DROP_REGIONS = {"V - unspecified"}


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


def sym_dff(F, fs):
    n, T = F.shape
    t = np.arange(T, dtype=float)
    A = np.vstack([t, np.ones(T)]).T
    out = np.empty((n, T))
    for i, tr in enumerate(F):
        a, b = np.linalg.lstsq(A, tr, rcond=None)[0]
        fit = a * t + b
        trd = tr * fit[0] / fit if fit.min() > 1.0 else tr.astype(float)
        F0 = max(np.median(trd), 1.0)
        out[i] = (trd - F0) / F0 * 100.0
    return out


def noise_sigma(d):
    """σ from the noise-dominated lower half: median − 15.87th percentile."""
    return np.maximum(np.median(d, axis=1) - np.percentile(d, 15.87, axis=1), 1e-6)


def count(tr, sig, K):
    kw = dict(height=K * sig, prominence=1.5 * sig, distance=DIST, width=1.0)
    p, _ = find_peaks(tr, **kw)
    n, _ = find_peaks(-tr, **kw)
    return len(p), len(n)


def main():
    colors = parse_labels_txt(LABELS_TXT)
    meta = pd.read_csv(os.path.join(OUT_DIR, "traces_metadata.csv"))
    fs = float(meta.fs_hz.iloc[0])
    T = np.load(os.path.join(OUT_DIR, "traces.npz"), allow_pickle=True)["time_s"].shape[0]
    dur_min = T / fs / 60.0
    N = len(meta)
    print(f"{N} neurons, {dur_min:.3f} min/recording")

    cnt = {K: np.zeros((N, 2), int) for K in KS}
    sig_all = np.zeros(N)
    for z, g in meta.groupby("z"):
        F = np.load(os.path.join(S2P_DIR, f"02F_2min_{z}", "suite2p", "plane0", "F.npy"))
        d = medfilt(sym_dff(F[g.cell_id.to_numpy()].astype(float), fs), kernel_size=(1, 3))
        sig = noise_sigma(d)
        for j, gi in enumerate(g.row_index.to_numpy()):
            sig_all[gi] = sig[j]
            for K in KS:
                cnt[K][gi] = count(d[j], sig[j], K)

    print(f"\nnoise σ (neg-side):  median {np.median(sig_all):.1f}%  IQR "
          f"[{np.percentile(sig_all, 25):.1f}, {np.percentile(sig_all, 75):.1f}] ΔF/F")
    print("\n K    pos/min  neg/min  FP frac")
    fp = {}
    for K in KS:
        p, n = cnt[K][:, 0], cnt[K][:, 1]
        fp[K] = n.sum() / max(p.sum() + n.sum(), 1)
        print(f" {K:>3}   {p.mean()/dur_min:6.3f}   {n.mean()/dur_min:6.3f}   {fp[K]:.3f}")
    Ksel = next((K for K in KS if fp[K] <= FP_TARGET), KS[-1])
    print(f"\nselected K = {Ksel}σ  (FP fraction {fp[Ksel]:.2f})")

    p, n = cnt[Ksel][:, 0], cnt[Ksel][:, 1]
    out = meta[["row_index", "neuron_uid", "z", "region_name"]].copy()
    out["noise_sigma_pct"] = sig_all.round(2)
    out["n_pos"] = p
    out["n_neg"] = n
    out["rate_pos_per_min"] = (p / dur_min).round(3)
    out["rate_neg_per_min"] = (n / dur_min).round(3)
    out["has_transient"] = (p >= 1).astype(int)
    out["K_sigma"] = Ksel
    out.to_csv(os.path.join(OUT_DIR, "transient_rate.csv"), index=False)

    # ---- figure: mean rate + FP floor + fraction-with-an-event, by region ----
    d = out[~out.region_name.isin(DROP_REGIONS)]
    order = d.groupby("region_name")["rate_pos_per_min"].mean().sort_values(ascending=False)
    regs = list(order.index)
    nreg = len(regs)

    fig, axes = plt.subplots(1, 2, figsize=(13, 6.5), sharey=True)
    ax = axes[0]
    for i, r in enumerate(regs):
        s = d[d.region_name == r]
        col = darken(colors.get(r, (.4, .4, .4)))
        ax.barh(i, s.rate_pos_per_min.mean(), height=.62, color=col, alpha=.85,
                edgecolor="0.3", lw=.5)
        ax.plot(s.rate_neg_per_min.mean(), i, "x", color="0.15", ms=8, mew=2, zorder=5)
    ax.set_yticks(range(nreg))
    ax.set_yticklabels([f"{r}  (n={int((d.region_name == r).sum())})" for r in regs], fontsize=9)
    ax.set_ylim(nreg - 0.5, -0.5)
    ax.set_xlabel(f"mean significant transients / min   (K = {Ksel}σ)")
    ax.set_title(f"A  ·  transient rate by region   (× = mean noise/FP rate)", loc="left",
                 fontsize=9.5, fontweight="bold")
    ax.grid(axis="x", alpha=.25)

    ax = axes[1]
    for i, r in enumerate(regs):
        s = d[d.region_name == r]
        col = darken(colors.get(r, (.4, .4, .4)))
        ax.barh(i, s.has_transient.mean(), height=.62, color=col, alpha=.85,
                edgecolor="0.3", lw=.5)
        ax.text(s.has_transient.mean() + 0.01, i, f"{s.has_transient.mean():.2f}",
                va="center", fontsize=8, color="0.2")
    ax.set_xlabel(f"fraction of cells with ≥1 significant transient in ~{dur_min:.1f} min")
    ax.set_title("B  ·  share of cells showing any event", loc="left", fontsize=9.5,
                 fontweight="bold")
    ax.grid(axis="x", alpha=.25)
    ax.set_xlim(0, 1)

    fig.suptitle("Significant-transient rate — symmetric-baseline ΔF/F, noise σ from the "
                 f"negative half of the distribution\ndataset-wide false-positive fraction "
                 f"≈ {fp[Ksel]:.0%} at K = {Ksel}σ · one animal · descriptive only", fontsize=9.5)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(os.path.join(OUT_DIR, "transient_rate.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    print("\nby region:  mean pos/min | mean neg/min (FP) | frac cells w/ >=1 event")
    for r in regs:
        s = d[d.region_name == r]
        print(f"  {r:16s} n={len(s):4d}   {s.rate_pos_per_min.mean():5.2f}   "
              f"{s.rate_neg_per_min.mean():5.2f}   {s.has_transient.mean():.2f}")
    print("\nwrote transient_rate.csv + transient_rate.png")


if __name__ == "__main__":
    main()
