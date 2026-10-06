#!/usr/bin/env python3
"""
Create a rotating 3D "video" (animated GIF) of Y1 vs Y2 seismicity.

Output:
  - 03-output-graphics/MyPng/seismicity_Y1_Y2_3D.gif

Notes:
  - Uses Pillow (no ffmpeg required).
  - Depth is positive downward (z-axis inverted to show 0 km at top).
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np

DEPTH_MAX_KM = 3.0
N_FRAMES = 72
FPS = 12
MAX_POINTS_Y1 = 12000
MAX_POINTS_Y2 = 25000
RNG_SEED = 7


def _configure_matplotlib() -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(Path("/tmp/matplotlib")))
    import matplotlib

    matplotlib.use("Agg")


def load_mldd_loc(loc_path: Path) -> dict[str, np.ndarray]:
    """
    Expected columns (11):
      YY MM DD HH MIN SS.sss Lat Lon Depth(km) Mag EventID
    """
    data = np.loadtxt(loc_path)
    if data.ndim != 2 or data.shape[1] != 11:
        raise ValueError(f"Unexpected format in {loc_path} (shape={data.shape})")

    return {"lat": data[:, 6], "lon": data[:, 7], "depth_km": data[:, 8]}


def _downsample(mask: np.ndarray, max_points: int, rng: np.random.Generator) -> np.ndarray:
    idx = np.flatnonzero(mask)
    if idx.size <= max_points:
        return mask
    keep = rng.choice(idx, size=max_points, replace=False)
    out = np.zeros_like(mask, dtype=bool)
    out[keep] = True
    return out


def deg_to_km_aspect(
    lon_min: float, lon_max: float, lat_min: float, lat_max: float, depth_range_km: float
) -> tuple[float, float, float]:
    mean_lat_rad = np.deg2rad((lat_min + lat_max) / 2.0)
    km_per_deg_lat = 111.0
    km_per_deg_lon = 111.0 * float(np.cos(mean_lat_rad))
    lon_km = (lon_max - lon_min) * km_per_deg_lon
    lat_km = (lat_max - lat_min) * km_per_deg_lat
    return max(lon_km, 1e-6), max(lat_km, 1e-6), max(depth_range_km, 1e-6)


def main() -> None:
    _configure_matplotlib()
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter

    repo_root = Path(__file__).resolve().parents[2]
    y1_path = repo_root / "02-data/A_all/axial.Y1.mldd.loc.260317"
    y2_path = repo_root / "02-data/A_all/axial.Y2.mldd.loc.260317"

    y1 = load_mldd_loc(y1_path)
    y2 = load_mldd_loc(y2_path)

    y1_mask = (
        np.isfinite(y1["depth_km"])
        & (y1["depth_km"] >= 0.0)
        & (y1["depth_km"] <= DEPTH_MAX_KM)
    )
    y2_mask = (
        np.isfinite(y2["depth_km"])
        & (y2["depth_km"] >= 0.0)
        & (y2["depth_km"] <= DEPTH_MAX_KM)
    )
    if not y1_mask.any():
        y1_mask = np.isfinite(y1["depth_km"])
    if not y2_mask.any():
        y2_mask = np.isfinite(y2["depth_km"])

    rng = np.random.default_rng(RNG_SEED)
    y1_mask = _downsample(y1_mask, MAX_POINTS_Y1, rng)
    y2_mask = _downsample(y2_mask, MAX_POINTS_Y2, rng)

    lon_min = float(min(y1["lon"][y1_mask].min(), y2["lon"][y2_mask].min()))
    lon_max = float(max(y1["lon"][y1_mask].max(), y2["lon"][y2_mask].max()))
    lat_min = float(min(y1["lat"][y1_mask].min(), y2["lat"][y2_mask].min()))
    lat_max = float(max(y1["lat"][y1_mask].max(), y2["lat"][y2_mask].max()))
    pad_lon = (lon_max - lon_min) * 0.05
    pad_lat = (lat_max - lat_min) * 0.05

    station_coords = {}
    stations_legacy: list[str] = []
    stations_main: list[str] = []
    try:
        from pipeline_config import STATION_COORDS, STATIONS_LEGACY, STATIONS_MAIN

        station_coords = dict(STATION_COORDS)
        stations_legacy = list(STATIONS_LEGACY)
        stations_main = list(STATIONS_MAIN)
    except Exception:
        pass

    fig = plt.figure(figsize=(13.0, 6.6), constrained_layout=True)
    ax1 = fig.add_subplot(1, 2, 1, projection="3d")
    ax2 = fig.add_subplot(1, 2, 2, projection="3d", sharex=ax1, sharey=ax1)

    scatter_kwargs = {
        "s": 2,
        "alpha": 0.55,
        "cmap": "viridis_r",
        "vmin": 0.0,
        "vmax": DEPTH_MAX_KM,
        "linewidths": 0,
        "rasterized": True,
    }
    sc1 = ax1.scatter(
        y1["lon"][y1_mask],
        y1["lat"][y1_mask],
        y1["depth_km"][y1_mask],
        c=y1["depth_km"][y1_mask],
        **scatter_kwargs,
    )
    sc2 = ax2.scatter(
        y2["lon"][y2_mask],
        y2["lat"][y2_mask],
        y2["depth_km"][y2_mask],
        c=y2["depth_km"][y2_mask],
        **scatter_kwargs,
    )

    for ax, title, n_plotted in [
        (ax1, "Y1 (2022-09-01–2023-08-31)", int(y1_mask.sum())),
        (ax2, "Y2 (2023-09-01–2024-09-22)", int(y2_mask.sum())),
    ]:
        ax.set_title(f"{title}\n(n={n_plotted:,}, depth≤{DEPTH_MAX_KM:.1f} km)", fontsize=10)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_zlabel("Depth (km)")
        ax.set_xlim(lon_min - pad_lon, lon_max + pad_lon)
        ax.set_ylim(lat_min - pad_lat, lat_max + pad_lat)
        ax.set_zlim(DEPTH_MAX_KM, 0.0)

        if station_coords and (stations_legacy or stations_main):
            legacy_xyz = [
                (station_coords[s][1], station_coords[s][0], 0.0)
                for s in stations_legacy
                if s in station_coords
            ]
            main_xyz = [
                (station_coords[s][1], station_coords[s][0], 0.0)
                for s in stations_main
                if s in station_coords
            ]
            if legacy_xyz:
                xs, ys, zs = zip(*legacy_xyz)
                ax.scatter(
                    xs,
                    ys,
                    zs,
                    marker="o",
                    s=40,
                    facecolors="none",
                    edgecolors="black",
                    linewidths=1.0,
                    zorder=10,
                    label="Stations (legacy)",
                )
            if main_xyz:
                xs, ys, zs = zip(*main_xyz)
                ax.scatter(
                    xs,
                    ys,
                    zs,
                    marker="^",
                    s=55,
                    facecolors="none",
                    edgecolors="black",
                    linewidths=1.0,
                    zorder=10,
                    label="Stations (main)",
                )

        lon_km, lat_km, dep_km = deg_to_km_aspect(
            lon_min - pad_lon, lon_max + pad_lon, lat_min - pad_lat, lat_max + pad_lat, DEPTH_MAX_KM
        )
        ax.set_box_aspect((lon_km, lat_km, dep_km))

    handles, _ = ax1.get_legend_handles_labels()
    if handles:
        ax1.legend(loc="upper left", fontsize=8, frameon=True)

    colorbar = fig.colorbar(sc2, ax=[ax1, ax2], shrink=0.85, pad=0.02)
    colorbar.set_label("Depth (km)")

    fig.suptitle("Axial Seamount 3D seismicity (rotating view)", fontsize=12)

    base_elev = 22
    base_azim = -60

    def update(frame: int) -> None:
        azim = base_azim + (360.0 * frame / N_FRAMES)
        ax1.view_init(elev=base_elev, azim=azim)
        ax2.view_init(elev=base_elev, azim=azim)

    anim = FuncAnimation(fig, update, frames=N_FRAMES, interval=1000 / FPS, blit=False)

    out_gif = repo_root / "03-output-graphics/MyPng/seismicity_Y1_Y2_3D.gif"
    writer = PillowWriter(fps=FPS)
    anim.save(out_gif, writer=writer, dpi=110)

    print(f"Wrote: {out_gif}")


if __name__ == "__main__":
    main()

