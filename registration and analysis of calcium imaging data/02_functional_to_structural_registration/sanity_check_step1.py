#!/usr/bin/env python3
"""
sanity_check_step1.py
========================
Stage 2 of 6 (functional-to-structural registration) -- Step 1 visual check.

Plots a mid-block functional plane next to several candidate structural-stack
slices, all at matched physical (micron) scale via matplotlib's `extent`, so a
flip, rotation, or gross scale error would be visually obvious before any
automated registration search is attempted.

Reads:
  - data_for_upload/02_functional_to_structural_registration/stage1_ministack/
    functional_ministack.nii.gz
  - data_for_upload/03_structural_to_template_registration/stage0_geometry_fixed/
    02F_stack_ras_cropped.nii.gz

Writes:
  - data_for_upload/02_functional_to_structural_registration/stage1_ministack/
    sanity_check.png

Run order: after build_functional_ministack.py, before
register_ministack_to_structural.py.

Environment: antspy env (nibabel, numpy, matplotlib; see this stage's README
for the full package list).
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np

DATA_ROOT = os.environ.get("PIPELINE_DATA_ROOT") or os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data_for_upload"))

MINISTACK_PATH = os.path.join(DATA_ROOT, "02_functional_to_structural_registration",
                               "stage1_ministack", "functional_ministack.nii.gz")
STRUCTURAL_PATH = os.path.join(DATA_ROOT, "03_structural_to_template_registration",
                                "stage0_geometry_fixed", "02F_stack_ras_cropped.nii.gz")
OUT_PNG = os.path.join(DATA_ROOT, "02_functional_to_structural_registration",
                        "stage1_ministack", "sanity_check.png")


def main():
    mini = nib.load(MINISTACK_PATH)
    mini_data = np.asarray(mini.dataobj)  # (X,Y,Z)
    mini_spacing = mini.header.get_zooms()[:3]
    print(f"Mini-stack shape (X,Y,Z): {mini_data.shape}, spacing um: {mini_spacing}")

    struct = nib.load(STRUCTURAL_PATH)
    struct_data = np.asarray(struct.dataobj)  # (X,Y,Z)
    struct_spacing = struct.header.get_zooms()[:3]
    print(f"Structural stack shape (X,Y,Z): {struct_data.shape}, spacing um: {struct_spacing}")

    mid_z = mini_data.shape[2] // 2
    mini_plane = mini_data[:, :, mid_z]  # (X,Y)
    mini_extent_um = (0, mini_data.shape[0] * mini_spacing[0], mini_data.shape[1] * mini_spacing[1], 0)

    struct_extent_um = (0, struct_data.shape[0] * struct_spacing[0], struct_data.shape[1] * struct_spacing[1], 0)

    # candidate structural z-indices to show alongside: full stack midpoint, and
    # top/bottom quartiles, since there is no prior on where the 36-plane block sits
    candidate_z = [
        int(struct_data.shape[2] * 0.15),
        int(struct_data.shape[2] * 0.5),
        int(struct_data.shape[2] * 0.85),
    ]

    fig, axes = plt.subplots(1, 4, figsize=(24, 7))

    ax = axes[0]
    hi = np.percentile(mini_plane[mini_plane > 0], 99) if np.any(mini_plane > 0) else mini_plane.max()
    ax.imshow(np.clip(mini_plane.T / max(hi, 1e-6), 0, 1), cmap="gray", extent=mini_extent_um, origin="upper")
    ax.set_title(f"functional mini-stack\nmid-block plane (idx {mid_z})\n{mini_extent_um[1]:.0f} x {mini_extent_um[2]:.0f} um")
    ax.set_xlabel("X (um)"); ax.set_ylabel("Y (um)")

    for ax, zi in zip(axes[1:], candidate_z):
        sp = struct_data[:, :, zi]
        hi = np.percentile(sp[sp > 0], 99) if np.any(sp > 0) else sp.max()
        ax.imshow(np.clip(sp.T / max(hi, 1e-6), 0, 1), cmap="gray", extent=struct_extent_um, origin="upper")
        ax.set_title(f"structural stack z-idx {zi}\n{struct_extent_um[1]:.0f} x {struct_extent_um[2]:.0f} um")
        ax.set_xlabel("X (um)"); ax.set_ylabel("Y (um)")

    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=130)
    print(f"Wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
