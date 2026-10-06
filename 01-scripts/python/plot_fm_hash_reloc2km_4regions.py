#!/usr/bin/env python3
"""
plot_fm_hash_reloc2km_4regions.py
=================================
Make "large-sphere" focal-mechanism map plots for four spatial regions using
HASH DL-all results rerun with Felix's NLLDD 2 km relocations.

Inputs (default):
  - 01-scripts/HASH/hashout_Y1_dl_all_reloc2km_1.dat
  - 01-scripts/HASH/hashout_Y2_dl_all_reloc2km_1.dat

Outputs (default):
  - 02-data/G_FM/focal_mechanisms_largesphere_group_1.png  (West)
  - 02-data/G_FM/focal_mechanisms_largesphere_group_2.png  (East)
  - 02-data/G_FM/focal_mechanisms_largesphere_group_3.png  (North)
  - 02-data/G_FM/focal_mechanisms_largesphere_group_4.png  (South/ID)

Run with FM_ML env (recommended):
  /opt/miniconda3/envs/FM_ML/bin/python 01-scripts/python/plot_fm_hash_reloc2km_4regions.py
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("/tmp") / "matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from obspy.imaging.beachball import beach


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from fault_classification import classify_fault_ptb  # noqa: E402


PROJECT = HERE.parents[1]
HASH_DIR = PROJECT / "01-scripts" / "HASH"
OUT_DIR = PROJECT / "02-data" / "G_FM"


# Map limits (wide enough to show all OO + 2F stations)
LON_LIM = (-130.07, -129.93)
LAT_LIM = (45.888, 46.025)


# Stations (OO + 2F)
STA_OO = {
    "AS1": (45.93356, -129.99920),
    "AS2": (45.93377, -130.01410),
    "CC1": (45.95468, -130.00890),
    "EC1": (45.94958, -129.97970),
    "EC2": (45.93967, -129.97380),
    "EC3": (45.93607, -129.97850),
    "ID1": (45.92573, -129.97800),
}
STA_2F = {
    "01B": (46.01935, -130.00526),
    "02B": (46.00111, -130.02930),
    "03B": (45.99250, -129.99397),
    "04B": (45.98834, -130.04896),
    "05B": (45.97327, -129.97728),
    "06B": (45.98178, -130.01423),
    "07B": (45.96950, -130.02963),
    "08B": (45.96336, -130.00441),
    "09B": (45.97194, -130.05972),
    "10B": (45.95867, -129.94982),
    "11B": (45.94851, -130.03801),
    "12B": (45.91473, -130.02173),
    "13B": (45.90472, -129.97240),
    "14B": (45.89913, -130.00612),
}


CALDERA_RIM = np.array(
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


FAULT_COLOR = {
    "N": {"AB": [0.15, 0.25, 0.85], "C": [0.50, 0.58, 0.92]},
    "R": {"AB": [0.85, 0.15, 0.15], "C": [0.92, 0.48, 0.48]},
    "S": {"AB": [0.10, 0.70, 0.20], "C": [0.45, 0.80, 0.52]},
    "U": {"AB": [0.40, 0.40, 0.40], "C": [0.62, 0.62, 0.62]},
}


REGIONS = [
    (
        "West",
        lambda lat, lon: (lon < -130.008) and (45.93 <= lat <= 45.96),
        dict(lon=(-130.07, -130.008), lat=(45.93, 45.96)),
    ),
    (
        "East",
        lambda lat, lon: (lon > -130.002) and (45.93 <= lat <= 45.97),
        dict(lon=(-130.002, -129.93), lat=(45.93, 45.97)),
    ),
    ("North", lambda lat, lon: (lat > 45.97), dict(lon=LON_LIM, lat=(45.97, LAT_LIM[1]))),
    ("South/ID", lambda lat, lon: (lat < 45.93), dict(lon=LON_LIM, lat=(LAT_LIM[0], 45.93))),
]


def load_hash_summary(path: Path) -> pd.DataFrame:
    rows = []
    with path.open() as f:
        for line in f:
            parts = line.split()
            if len(parts) < 29:
                continue
            try:
                rows.append(
                    dict(
                        origin_lat=float(parts[10]),
                        origin_lon=float(parts[11]),
                        origin_depth_km=float(parts[12]),
                        strike=float(parts[21]),
                        dip=float(parts[22]),
                        rake=float(parts[23]),
                        quality=str(parts[28]).strip(),
                    )
                )
            except ValueError:
                continue
    return pd.DataFrame(rows)


def parse_args():
    p = argparse.ArgumentParser(description="Plot relocated HASH DL-all focal mechanisms in 4 spatial regions.")
    p.add_argument(
        "--inputs",
        nargs="+",
        type=Path,
        default=[
            HASH_DIR / "hashout_Y1_dl_all_reloc2km_1.dat",
            HASH_DIR / "hashout_Y2_dl_all_reloc2km_1.dat",
        ],
        help="HASH summary outputs (hashout_*_1.dat).",
    )
    p.add_argument("--qualities", default="ABC", help="Quality letters to include (default: ABC).")
    p.add_argument("--bb_width", type=float, default=0.0016, help="Beachball width in degrees (default: 0.0016).")
    p.add_argument("--out_dir", type=Path, default=OUT_DIR, help="Output directory (default: 02-data/G_FM).")
    return p.parse_args()


def region_label(lat: float, lon: float) -> str:
    """
    Assign one of the four regions using ordered precedence.
    Events not matching any region return 'Other'.
    """
    for name, predicate, _ in REGIONS:
        if predicate(lat, lon):
            return name
    return "Other"


def plot_region(df: pd.DataFrame, *, name: str, bb_width: float, out_png: Path) -> None:
    out_png.parent.mkdir(parents=True, exist_ok=True)

    sub = df[df["region"] == name].copy()
    # Draw A/B on top of C
    qrank = {"C": 0, "B": 1, "A": 2}
    sub["qrank"] = sub["quality"].map(qrank).fillna(-1)
    sub = sub.sort_values("qrank")

    # Figure aspect: match lon/lat ratio at this latitude so beach balls stay round
    mean_lat = float(np.mean(LAT_LIM))
    lon_range = abs(LON_LIM[1] - LON_LIM[0])
    lat_range = abs(LAT_LIM[1] - LAT_LIM[0])
    aspect = (lon_range * np.cos(np.radians(mean_lat))) / max(lat_range, 1e-9)
    panel_w = 8.0
    panel_h = panel_w / max(aspect, 1e-9)

    fig, ax = plt.subplots(figsize=(panel_w, panel_h), facecolor="white")
    ax.set_facecolor([0.85, 0.92, 0.97])
    ax.set_aspect(1.0 / np.cos(np.radians(mean_lat)))

    for _, r in sub.iterrows():
        qk = "AB" if r["quality"] in ("A", "B") else "C"
        try:
            bb = beach(
                [r["strike"], r["dip"], r["rake"]],
                xy=(r["origin_lon"], r["origin_lat"]),
                width=bb_width,
                linewidth=0.20 if qk == "AB" else 0.12,
                facecolor=FAULT_COLOR[r["fault_type"]][qk],
                alpha=0.88 if qk == "AB" else 0.62,
            )
            bb.set_transform(ax.transData)
            bb.set_zorder(2 if qk == "AB" else 1)
            ax.add_collection(bb)
        except Exception:
            pass

    # Caldera rim + stations
    ax.plot(CALDERA_RIM[:, 0], CALDERA_RIM[:, 1], "-", color="black", linewidth=1.2, zorder=7)
    ax.plot(
        [v[1] for v in STA_OO.values()],
        [v[0] for v in STA_OO.values()],
        "^",
        color="yellow",
        markeredgecolor="k",
        markeredgewidth=0.4,
        markersize=5,
        zorder=8,
        label="OO",
    )
    ax.plot(
        [v[1] for v in STA_2F.values()],
        [v[0] for v in STA_2F.values()],
        "v",
        color="cyan",
        markeredgecolor="k",
        markeredgewidth=0.4,
        markersize=5,
        zorder=8,
        label="2F",
    )

    # Region box (for reference)
    region_meta = {n: m for n, _, m in REGIONS}.get(name)
    if region_meta:
        x0, x1 = region_meta["lon"]
        y0, y1 = region_meta["lat"]
        ax.add_patch(
            matplotlib.patches.Rectangle(
                (x0, y0),
                x1 - x0,
                y1 - y0,
                linewidth=1.8,
                edgecolor="black",
                facecolor="none",
                linestyle="--",
                zorder=9,
            )
        )

    ax.set_xlim(LON_LIM)
    ax.set_ylim(LAT_LIM)
    ticks = np.arange(-130.06, -129.92, 0.04)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:.2f}" for t in ticks], fontsize=8, rotation=45)
    ax.tick_params(axis="y", labelsize=8)
    ax.grid(True, linewidth=0.35, color="gray", alpha=0.5)
    ax.set_xlabel("Longitude", fontsize=10)
    ax.set_ylabel("Latitude", fontsize=10)

    # Legend: fault types + quality shading
    handles = [
        mpatches.Patch(color=FAULT_COLOR["N"]["AB"], label="Normal (A/B)"),
        mpatches.Patch(color=FAULT_COLOR["R"]["AB"], label="Reverse (A/B)"),
        mpatches.Patch(color=FAULT_COLOR["S"]["AB"], label="Strike-slip (A/B)"),
        mpatches.Patch(color=FAULT_COLOR["U"]["AB"], label="Oblique (A/B)"),
        mpatches.Patch(color=FAULT_COLOR["N"]["C"], label="Normal (C)"),
        mpatches.Patch(color=FAULT_COLOR["R"]["C"], label="Reverse (C)"),
        mpatches.Patch(color=FAULT_COLOR["S"]["C"], label="Strike-slip (C)"),
        mpatches.Patch(color=FAULT_COLOR["U"]["C"], label="Oblique (C)"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=8, framealpha=0.92, ncol=2)

    n_ab = int(sub["quality"].isin(["A", "B"]).sum())
    n_c = int((sub["quality"] == "C").sum())
    fig.suptitle(
        f"Reloc2km HASH DL-all (Y1+Y2) — {name} region\n"
        f"Q=A+B+C | n={len(sub):,}  (A/B={n_ab:,}, C={n_c:,})",
        fontsize=12,
        y=0.98,
    )
    fig.savefig(out_png, dpi=170, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()

    rows = []
    for p in args.inputs:
        if not p.exists():
            print(f"WARNING: missing input: {p}")
            continue
        rows.append(load_hash_summary(p))
    if not rows:
        raise SystemExit("No events loaded — check input paths.")

    df = pd.concat(rows, ignore_index=True)
    q_keep = set(args.qualities.strip().upper())
    df = df[df["quality"].isin(sorted(q_keep))].copy()

    df["fault_type"] = [classify_fault_ptb(s, d, r) for s, d, r in zip(df["strike"], df["dip"], df["rake"])]
    df["region"] = [region_label(lat, lon) for lat, lon in zip(df["origin_lat"], df["origin_lon"])]

    args.out_dir.mkdir(parents=True, exist_ok=True)

    group_names = ["West", "East", "North", "South/ID"]
    for i, name in enumerate(group_names, start=1):
        out = args.out_dir / f"focal_mechanisms_largesphere_group_{i}.png"
        plot_region(df, name=name, bb_width=args.bb_width, out_png=out)
        print(f"Saved {out}")

    n_other = int((df["region"] == "Other").sum())
    if n_other:
        print(f"Note: {n_other:,} events did not match any of the 4 region masks (region='Other').")


if __name__ == "__main__":
    main()
