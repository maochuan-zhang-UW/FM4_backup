#!/usr/bin/env python3
"""
plot_focal_mechanisms_16x9.py
=============================
Single-slide 16:9 focal-symbol figure (map + y-slice panels) using simplified
mechanism glyphs (circle + strike tick), similar to the reference layout you
shared.

Default inputs (FM7 repo):
  - 01-scripts/HASH/hashout_Y1_dl_all_1.dat
  - 01-scripts/HASH/hashout_Y2_dl_all_1.dat

Default filter:
  - Quality: A+B+C

Output:
  - 02-data/G_FM/fm_hash_y_slices_16x9_ABC.png

Run (FM_ML env recommended):
  /opt/miniconda3/envs/FM_ML/bin/python plot_focal_mechanisms_16x9.py
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("/tmp") / "matplotlib"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from fault_classification import classify_fault_ptb  # noqa: E402


PROJECT = HERE.parents[1]
HASH_DIR = PROJECT / "01-scripts" / "HASH"
OUT_DIR = PROJECT / "02-data" / "G_FM"


SLICE_CENTERS = np.array([5, 4, 3, 2, 1, 0, -1, -2, -3, -4], dtype=float)
HALF_WIDTH = 0.5
DEPTH_LIM = (0.0, 2.5)
X_LIM = (-2.5, 2.5)
Y_LIM = (-6.8, 6.8)


NORMAL_COLOR = "#0000ff"
REVERSE_COLOR = "#ff0000"
STRIKESLIP_COLOR = "#00aa00"
OBLIQUE_COLOR = "#666666"


CALDERA_RIM_LONLAT = np.array(
    [
        [-130.004785563058, 45.9207755734405],
        [-130.010476202888, 45.9238241104543],
        [-130.018881564079, 45.9351908809594],
        [-130.023946125193, 45.9412238501725],
        [-130.028718653506, 45.9498812001140],
        [-130.030451219380, 45.9511765797916],
        [-130.030679485650, 45.9542732243167],
        [-130.031733279709, 45.9558130656063],
        [-130.031444653500, 45.9586760104296],
        [-130.036188782208, 45.9656647517656],
        [-130.036950110789, 45.9698291665232],
        [-130.039953347000, 45.9750458167927],
        [-130.038595675479, 45.9847117727418],
        [-130.035927416999, 45.9883113986506],
        [-130.018067675296, 45.9933582886740],
        [-130.013629193751, 45.9937552841350],
        [-130.010365710979, 45.9929499241491],
        [-130.008647442296, 45.9924883829037],
        [-130.007262470669, 45.9915471582374],
        [-130.006042022411, 45.9902469280907],
        [-130.005178629490, 45.9897778053610],
        [-130.001868199523, 45.9863506519894],
        [-130.001154359192, 45.9846883853932],
        [-130.000949059432, 45.9827833001814],
        [-129.999393534330, 45.9818434725493],
        [-129.997797388662, 45.9786395525337],
        [-129.995357566829, 45.9760388622191],
        [-129.993956176267, 45.9741441737512],
        [-129.993678708114, 45.9681875631427],
        [-129.993140494256, 45.9667620754035],
        [-129.992087550788, 45.9652218741086],
        [-129.991186410747, 45.9626077204113],
        [-129.989604036931, 45.9601186407320],
        [-129.989238986151, 45.9588108137369],
        [-129.989728217453, 45.9574955894078],
        [-129.985484098670, 45.9494279735802],
        [-129.984788122490, 45.9487188881587],
    ],
    dtype=float,
)


def latlon2xy(lat, lon, *, lat0=45.9547, lon0=-130.0089, rot=-20.0):
    """Rotate lat/lon into local x,y (km) coordinates."""
    xltkm = 111.19
    xlnkm = xltkm * np.cos(np.deg2rad(lat0))
    dlat = (np.asarray(lat) - lat0) * xltkm
    dlon = (np.asarray(lon) - lon0) * xlnkm
    snr = np.sin(np.deg2rad(rot))
    csr = np.cos(np.deg2rad(rot))
    y = csr * dlat + snr * dlon
    x = csr * dlon - snr * dlat
    return x, y


def load_hash_file(path: Path):
    """Parse HASH driver3 output file 1 (summary, one line per event)."""
    rows = []
    with path.open() as f:
        for line in f:
            parts = line.split()
            if len(parts) < 29:
                continue
            try:
                rows.append(
                    {
                        "lat": float(parts[10]),
                        "lon": float(parts[11]),
                        "dep": float(parts[12]),
                        "strike": float(parts[21]),
                        "dip": float(parts[22]),
                        "rake": float(parts[23]),
                        "quality": parts[28],
                    }
                )
            except ValueError:
                continue
    return rows


def mechanism_color_rake3(rake: float) -> str:
    """3-class rake-only scheme (no oblique): Normal / Reverse / Strike-slip."""
    if -135.0 <= rake <= -45.0:
        return NORMAL_COLOR
    if 45.0 <= rake <= 135.0:
        return REVERSE_COLOR
    return STRIKESLIP_COLOR


def edge_colors_for_catalog(strike, dip, rake, scheme: str) -> np.ndarray:
    if scheme == "rake3":
        return np.array([mechanism_color_rake3(rk) for rk in rake])

    # P/T/B plunge scheme (MATLAB-equivalent) with oblique shown separately.
    colors = {"N": NORMAL_COLOR, "R": REVERSE_COLOR, "S": STRIKESLIP_COLOR, "U": OBLIQUE_COLOR}
    return np.array([colors[classify_fault_ptb(s, d, r)] for s, d, r in zip(strike, dip, rake)])


def build_segments(xs, ys, strikes, length_x, length_y):
    theta = np.deg2rad(strikes)
    dx = 0.5 * length_x * np.sin(theta)
    dy = 0.5 * length_y * np.cos(theta)
    return np.stack(
        [
            np.column_stack([xs - dx, ys - dy]),
            np.column_stack([xs + dx, ys + dy]),
        ],
        axis=1,
    )


def add_mechanism_symbols(ax, xs, ys, strikes, edge_colors, marker_size, tick_len_x, tick_len_y, alpha):
    ax.scatter(
        xs,
        ys,
        s=marker_size,
        facecolors="white",
        edgecolors=edge_colors,
        linewidths=0.9,
        alpha=alpha,
        zorder=3,
    )
    segments = build_segments(xs, ys, strikes, tick_len_x, tick_len_y)
    lc = LineCollection(segments, colors="k", linewidths=0.55, alpha=min(1.0, alpha + 0.1), zorder=4)
    ax.add_collection(lc)


def parse_args():
    p = argparse.ArgumentParser(description="Plot HASH focal mechanisms in a 16:9 y-slice layout.")
    p.add_argument(
        "--inputs",
        nargs="+",
        type=Path,
        default=[
            HASH_DIR / "hashout_Y1_dl_all_1.dat",
            HASH_DIR / "hashout_Y2_dl_all_1.dat",
        ],
        help="One or more HASH output files (hashout_*.dat).",
    )
    p.add_argument(
        "--qualities",
        default="ABC",
        help="Quality letters to include (default: ABC).",
    )
    p.add_argument(
        "--scheme",
        choices=["rake3", "ptb4"],
        default="rake3",
        help="Fault-type color scheme: rake3 (Normal/Reverse/Strike-slip) or ptb4 (+Oblique).",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=OUT_DIR / "fm_hash_y_slices_16x9_ABC.png",
        help="Output PNG path.",
    )
    return p.parse_args()


def main():
    args = parse_args()

    OUT_DIR.mkdir(exist_ok=True, parents=True)
    args.out.parent.mkdir(exist_ok=True, parents=True)

    try:
        plt.rcParams["font.family"] = "Arial"
    except Exception:
        pass

    rows = []
    for p in args.inputs:
        if not p.exists():
            print(f"WARNING: missing input: {p}")
            continue
        rows.extend(load_hash_file(p))
    if not rows:
        raise SystemExit("No events loaded — check input paths.")

    catalog = {key: np.array([row[key] for row in rows]) for key in rows[0]}
    q_keep = set(args.qualities.strip().upper())
    mask_q = np.isin(catalog["quality"], sorted(q_keep))
    for k in list(catalog.keys()):
        catalog[k] = catalog[k][mask_q]

    catalog["x"], catalog["y"] = latlon2xy(catalog["lat"], catalog["lon"])
    rim_x, rim_y = latlon2xy(CALDERA_RIM_LONLAT[:, 1], CALDERA_RIM_LONLAT[:, 0])

    edge_colors = edge_colors_for_catalog(catalog["strike"], catalog["dip"], catalog["rake"], scheme=args.scheme)

    fig = plt.figure(figsize=(16, 9), facecolor="white")
    outer = GridSpec(
        1,
        2,
        figure=fig,
        width_ratios=[1.0, 2.45],
        wspace=0.12,
        left=0.05,
        right=0.98,
        top=0.93,
        bottom=0.09,
    )
    map_ax = fig.add_subplot(outer[0, 0])
    right = GridSpecFromSubplotSpec(5, 2, subplot_spec=outer[0, 1], wspace=0.10, hspace=0.18)

    inside_map = (
        (catalog["x"] >= X_LIM[0] - 0.2)
        & (catalog["x"] <= X_LIM[1] + 0.2)
        & (catalog["y"] >= Y_LIM[0])
        & (catalog["y"] <= Y_LIM[1])
    )

    for y0 in SLICE_CENTERS:
        map_ax.axhline(y0, color="0.75", lw=0.8, alpha=0.8, zorder=1)
        map_ax.text(X_LIM[1] + 0.08, y0, f"{y0:+.0f}", color="0.35", fontsize=10, va="center", fontweight="bold")

    add_mechanism_symbols(
        map_ax,
        catalog["x"][inside_map],
        catalog["y"][inside_map],
        catalog["strike"][inside_map],
        edge_colors[inside_map],
        marker_size=18,
        tick_len_x=0.14,
        tick_len_y=0.14,
        alpha=0.65,
    )
    map_ax.plot(rim_x, rim_y, color="0.25", lw=1.2, zorder=5)
    map_ax.set_xlim(X_LIM)
    map_ax.set_ylim(Y_LIM)
    map_ax.set_aspect("equal")
    map_ax.set_xticks(np.arange(-2, 3, 1))
    map_ax.set_yticks(np.arange(-6, 7, 2))
    map_ax.tick_params(labelsize=10)
    map_ax.set_title("a  Focal Mechanisms by Slice", loc="left", fontsize=15, fontweight="bold")
    map_ax.grid(True, lw=0.35, alpha=0.25)

    if args.scheme == "ptb4":
        legend_txt = "Normal = blue   Reverse = red   Strike-slip = green   Oblique = gray"
    else:
        legend_txt = "Normal = blue   Reverse = red   Strike-slip = green"
    map_ax.text(0.02, 0.03, legend_txt, transform=map_ax.transAxes, fontsize=9, va="bottom")

    letters = [chr(ord("b") + i) for i in range(len(SLICE_CENTERS))]
    for idx, y0 in enumerate(SLICE_CENTERS):
        col = 0 if idx < 5 else 1
        row_i = idx if idx < 5 else idx - 5
        ax = fig.add_subplot(right[row_i, col])
        mask = np.abs(catalog["y"] - y0) < HALF_WIDTH
        mask &= (catalog["dep"] >= DEPTH_LIM[0]) & (catalog["dep"] <= DEPTH_LIM[1])
        mask &= (catalog["x"] >= X_LIM[0] - 0.2) & (catalog["x"] <= X_LIM[1] + 0.2)
        add_mechanism_symbols(
            ax,
            catalog["x"][mask],
            catalog["dep"][mask],
            catalog["strike"][mask],
            edge_colors[mask],
            marker_size=14,
            tick_len_x=0.10,
            tick_len_y=0.10,
            alpha=0.62,
        )
        ax.axvline(0, color="0.65", lw=0.8, ls="--", zorder=0)
        ax.set_xlim(X_LIM)
        ax.set_ylim(DEPTH_LIM[1], DEPTH_LIM[0])
        ax.set_aspect("equal", adjustable="box")
        ax.grid(True, lw=0.30, alpha=0.22)
        ax.tick_params(labelsize=8.5, pad=1)
        ax.text(0.01, 0.84, f"{letters[idx]}  y={y0:+.0f} km", transform=ax.transAxes, fontsize=9.5, fontweight="bold")
        ax.set_xticks(np.arange(-2, 3, 1))
        ax.set_yticks([0, 1, 2])
        if col == 1:
            ax.set_yticklabels([])
        if row_i < 4:
            ax.set_xticklabels([])
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)

    fig.suptitle("Axial focal mechanisms by y slice", fontsize=22, y=0.975)
    fig.text(0.50, 0.02, "Single-slide 16:9 layout with all focal-mechanism symbols", ha="center", fontsize=12)
    fig.savefig(args.out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()
