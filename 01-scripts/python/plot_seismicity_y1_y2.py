#!/usr/bin/env python3
"""
Plot Axial Seamount Y1 vs Y2 seismicity from MLDD `.loc` files.

Output:
  - 03-output-graphics/MyPng/seismicity_Y1_Y2.png
  - 03-output-graphics/MyPdf/seismicity_Y1_Y2.pdf
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np

DEPTH_MAX_KM = 3.0


def _configure_matplotlib() -> None:
    # Avoid writing to ~/.matplotlib in restricted environments.
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

    return {
        "year": data[:, 0].astype(int),
        "month": data[:, 1].astype(int),
        "day": data[:, 2].astype(int),
        "hour": data[:, 3].astype(int),
        "minute": data[:, 4].astype(int),
        "second": data[:, 5],
        "lat": data[:, 6],
        "lon": data[:, 7],
        "depth_km": data[:, 8],
        "mag": data[:, 9],
        "event_id": data[:, 10].astype(int),
    }


def station_code_to_short(station_code: str) -> str:
    station_code = station_code.strip()
    # Examples in `.pha`:
    # - OOAXAS1  -> AS1
    # - 2FAX07A  -> 07A
    if len(station_code) >= 4 and station_code[2:4] == "AX":
        return station_code[4:]
    if station_code.startswith("AX"):
        return station_code[2:]
    return station_code


def parse_station_p_coverage(pha_path: Path) -> tuple[dict[str, int], int]:
    """
    Count events per station (P picks) from a HypoDD-format `.pha` file.

    Returns:
      - dict: station -> number of events with at least one P pick
      - int: number of event blocks encountered
    """
    counts: dict[str, int] = {}
    stations_in_event: set[str] = set()
    n_events = 0
    has_event = False

    with pha_path.open("rt", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue

            if line.startswith("#"):
                if has_event:
                    for sta in stations_in_event:
                        counts[sta] = counts.get(sta, 0) + 1
                n_events += 1
                has_event = True
                stations_in_event = set()
                continue

            parts = line.split()
            if len(parts) < 4:
                continue
            if parts[-1] != "P":
                continue
            stations_in_event.add(station_code_to_short(parts[0]))

    if has_event:
        for sta in stations_in_event:
            counts[sta] = counts.get(sta, 0) + 1

    return counts, n_events


def main() -> None:
    _configure_matplotlib()
    import matplotlib.pyplot as plt

    repo_root = Path(__file__).resolve().parents[2]
    y1_path = repo_root / "02-data/A_all/axial.Y1.mldd.loc.260317"
    y2_path = repo_root / "02-data/A_all/axial.Y2.mldd.loc.260317"
    y1_pha = repo_root / "02-data/A_all/axial.Y1.mldd.pha.260317"
    y2_pha = repo_root / "02-data/A_all/axial.Y2.mldd.pha.260317"

    y1 = load_mldd_loc(y1_path)
    y2 = load_mldd_loc(y2_path)

    station_coords: dict[str, tuple[float, float, float]] = {}
    stations_legacy: list[str] = []
    stations_main: list[str] = []
    stations_all: list[str] = []
    try:
        from pipeline_config import STATION_COORDS, STATIONS_LEGACY, STATIONS_MAIN, STATIONS_ALL

        station_coords = dict(STATION_COORDS)
        stations_legacy = list(STATIONS_LEGACY)
        stations_main = list(STATIONS_MAIN)
        stations_all = list(STATIONS_ALL)
    except Exception:
        pass

    coverage_y1_counts, _ = parse_station_p_coverage(y1_pha) if y1_pha.exists() else ({}, 0)
    coverage_y2_counts, _ = parse_station_p_coverage(y2_pha) if y2_pha.exists() else ({}, 0)
    if not stations_all and (coverage_y1_counts or coverage_y2_counts):
        stations_all = sorted(set(coverage_y1_counts) | set(coverage_y2_counts))

    # Shallow-depth view for side panels
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

    depth_min = 0.0
    depth_max = DEPTH_MAX_KM
    lon_min = float(min(y1["lon"][y1_mask].min(), y2["lon"][y2_mask].min()))
    lon_max = float(max(y1["lon"][y1_mask].max(), y2["lon"][y2_mask].max()))
    lat_min = float(min(y1["lat"][y1_mask].min(), y2["lat"][y2_mask].min()))
    lat_max = float(max(y1["lat"][y1_mask].max(), y2["lat"][y2_mask].max()))

    pad_lon = (lon_max - lon_min) * 0.05
    pad_lat = (lat_max - lat_min) * 0.05

    fig = plt.figure(figsize=(13.6, 9.0), constrained_layout=True)
    gs = fig.add_gridspec(
        3,
        4,
        height_ratios=[1.0, 0.38, 0.75],
        width_ratios=[1.0, 0.42, 1.0, 0.42],
    )
    ax_map_y1 = fig.add_subplot(gs[0, 0])
    ax_right_y1 = fig.add_subplot(gs[0, 1], sharey=ax_map_y1)
    ax_map_y2 = fig.add_subplot(gs[0, 2], sharex=ax_map_y1, sharey=ax_map_y1)
    ax_right_y2 = fig.add_subplot(gs[0, 3], sharey=ax_map_y1)
    ax_bottom_y1 = fig.add_subplot(gs[1, 0], sharex=ax_map_y1)
    ax_blank_y1 = fig.add_subplot(gs[1, 1])
    ax_blank_y1.axis("off")
    ax_bottom_y2 = fig.add_subplot(gs[1, 2], sharex=ax_map_y1)
    ax_blank_y2 = fig.add_subplot(gs[1, 3])
    ax_blank_y2.axis("off")
    ax_cov = fig.add_subplot(gs[2, 0:2])
    ax_mag = fig.add_subplot(gs[2, 2:4])

    # Maps (lon/lat) colored by depth
    sc_y1 = ax_map_y1.scatter(
        y1["lon"][y1_mask],
        y1["lat"][y1_mask],
        c=y1["depth_km"][y1_mask],
        s=3,
        alpha=0.6,
        cmap="viridis_r",
        vmin=depth_min,
        vmax=depth_max,
        linewidths=0,
    )
    sc_y2 = ax_map_y2.scatter(
        y2["lon"][y2_mask],
        y2["lat"][y2_mask],
        c=y2["depth_km"][y2_mask],
        s=3,
        alpha=0.6,
        cmap="viridis_r",
        vmin=depth_min,
        vmax=depth_max,
        linewidths=0,
    )

    for ax, label, n_plotted in [
        (ax_map_y1, "Y1 (2022-09-01–2023-08-31)", int(y1_mask.sum())),
        (ax_map_y2, "Y2 (2023-09-01–2024-09-22)", int(y2_mask.sum())),
    ]:
        ax.set_title(
            f"{label}\n(n={n_plotted:,}, depth≤{DEPTH_MAX_KM:.1f} km)", fontsize=10
        )
        ax.grid(True, linewidth=0.3, alpha=0.3)
        ax.set_xlim(lon_min - pad_lon, lon_max + pad_lon)
        ax.set_ylim(lat_min - pad_lat, lat_max + pad_lat)
        ax.tick_params(labelbottom=False)

        if station_coords and (stations_legacy or stations_main):
            legacy_xy = [
                (station_coords[s][1], station_coords[s][0])
                for s in stations_legacy
                if s in station_coords
            ]
            main_xy = [
                (station_coords[s][1], station_coords[s][0])
                for s in stations_main
                if s in station_coords
            ]
            if legacy_xy:
                lon_s, lat_s = zip(*legacy_xy)
                ax.scatter(
                    lon_s,
                    lat_s,
                    marker="o",
                    s=40,
                    facecolors="none",
                    edgecolors="black",
                    linewidths=1.0,
                    zorder=6,
                    label="Stations (legacy)",
                )
            if main_xy:
                lon_s, lat_s = zip(*main_xy)
                ax.scatter(
                    lon_s,
                    lat_s,
                    marker="^",
                    s=55,
                    facecolors="none",
                    edgecolors="black",
                    linewidths=1.0,
                    zorder=6,
                    label="Stations (main)",
                )

    ax_map_y1.set_ylabel("Latitude")
    fig.suptitle("Axial Seamount seismicity (MLDD locations)", fontsize=12)

    # Right panels: depth (x) vs latitude (y)
    profile_kwargs = {"s": 2, "alpha": 0.25, "linewidths": 0}
    ax_right_y1.scatter(
        y1["depth_km"][y1_mask], y1["lat"][y1_mask], color="tab:blue", **profile_kwargs
    )
    ax_right_y2.scatter(
        y2["depth_km"][y2_mask],
        y2["lat"][y2_mask],
        color="tab:orange",
        **profile_kwargs,
    )
    for ax in [ax_right_y1, ax_right_y2]:
        ax.set_xlim(depth_min, depth_max)
        ax.set_xlabel("Depth (km)")
        ax.grid(True, linewidth=0.3, alpha=0.3)
        ax.tick_params(labelleft=False)

    # Bottom panels: longitude (x) vs depth (y)
    ax_bottom_y1.scatter(
        y1["lon"][y1_mask], y1["depth_km"][y1_mask], color="tab:blue", **profile_kwargs
    )
    ax_bottom_y2.scatter(
        y2["lon"][y2_mask],
        y2["depth_km"][y2_mask],
        color="tab:orange",
        **profile_kwargs,
    )
    for ax in [ax_bottom_y1, ax_bottom_y2]:
        ax.set_ylim(depth_min, depth_max)
        ax.invert_yaxis()
        ax.set_xlabel("Longitude")
        ax.grid(True, linewidth=0.3, alpha=0.3)
    ax_bottom_y1.set_ylabel("Depth (km)")

    # Shared colorbar for map depth
    colorbar = fig.colorbar(
        sc_y2,
        ax=[ax_map_y1, ax_right_y1, ax_bottom_y1, ax_map_y2, ax_right_y2, ax_bottom_y2],
        shrink=0.98,
        pad=0.01,
    )
    colorbar.set_label("Depth (km)")

    handles, _ = ax_map_y1.get_legend_handles_labels()
    if handles:
        ax_map_y1.legend(loc="lower left", fontsize=8, frameon=True)

    # Station coverage (from `.pha`, P picks)
    if stations_all and (coverage_y1_counts or coverage_y2_counts):
        total_y1 = float(y1["event_id"].size)
        total_y2 = float(y2["event_id"].size)
        cov_y1 = np.array(
            [coverage_y1_counts.get(s, 0) / total_y1 * 100.0 for s in stations_all]
        )
        cov_y2 = np.array(
            [coverage_y2_counts.get(s, 0) / total_y2 * 100.0 for s in stations_all]
        )
        x = np.arange(len(stations_all))
        width = 0.42
        ax_cov.bar(
            x - width / 2, cov_y1, width=width, label="Y1", color="tab:blue", alpha=0.85
        )
        ax_cov.bar(
            x + width / 2,
            cov_y2,
            width=width,
            label="Y2",
            color="tab:orange",
            alpha=0.85,
        )
        ax_cov.set_xticks(x)
        ax_cov.set_xticklabels(stations_all, rotation=90, fontsize=8)
        ax_cov.set_ylabel("Events with P pick (%)")
        ax_cov.set_title("Station coverage (from .pha P picks)", fontsize=10)
        ax_cov.grid(True, axis="y", linewidth=0.3, alpha=0.3)
        ax_cov.set_ylim(0, 100)
        ax_cov.legend(fontsize=8, frameon=True)
    else:
        ax_cov.axis("off")
        ax_cov.text(
            0.5,
            0.5,
            "Station coverage unavailable",
            ha="center",
            va="center",
            fontsize=10,
        )

    # Magnitude distribution (from `.loc` Mag column)
    mag_y1 = y1["mag"]
    mag_y2 = y2["mag"]
    mag_min = float(min(mag_y1.min(), mag_y2.min()))
    mag_max = float(max(mag_y1.max(), mag_y2.max()))
    bins = np.arange(
        np.floor(mag_min * 10) / 10, np.ceil(mag_max * 10) / 10 + 0.1, 0.1
    )
    ax_mag.hist(mag_y1, bins=bins, density=True, alpha=0.55, label="Y1", color="tab:blue")
    ax_mag.hist(
        mag_y2, bins=bins, density=True, alpha=0.55, label="Y2", color="tab:orange"
    )
    ax_mag.set_title("Magnitude distribution (MLDD .loc)", fontsize=10)
    ax_mag.set_xlabel("Magnitude")
    ax_mag.set_ylabel("Density")
    ax_mag.grid(True, linewidth=0.3, alpha=0.3)
    ax_mag.legend(fontsize=8, frameon=True)

    out_png = repo_root / "03-output-graphics/MyPng/seismicity_Y1_Y2.png"
    out_pdf = repo_root / "03-output-graphics/MyPdf/seismicity_Y1_Y2.pdf"
    fig.savefig(out_png, dpi=300)
    fig.savefig(out_pdf)

    print(f"Wrote: {out_png}")
    print(f"Wrote: {out_pdf}")


if __name__ == "__main__":
    main()
