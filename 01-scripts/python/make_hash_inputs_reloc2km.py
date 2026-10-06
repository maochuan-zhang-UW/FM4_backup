#!/usr/bin/env python3
"""
make_hash_inputs_reloc2km.py
============================
Update HASH `phase_*.dat` files to use Felix's NLLDD 2 km relocation locations,
while keeping the *same* DL polarities and amplitudes (i.e., only hypocenter
lat/lon/depth changes).

Creates new files (does not overwrite the originals):
  - 01-scripts/HASH/phase_Y1_dl_all_reloc2km.dat
  - 01-scripts/HASH/hash.input_Y1_dl_all_reloc2km
  - 01-scripts/HASH/phase_Y2_dl_all_reloc2km.dat
  - 01-scripts/HASH/hash.input_Y2_dl_all_reloc2km

Inputs:
  - Existing HASH inputs in 01-scripts/HASH/
  - Existing cluster_idx ↔ event_id maps in 02-data/G_FM/
  - Relocation files downloaded to 02-data/A_all/reloc2km/

Run:
  python3 01-scripts/python/make_hash_inputs_reloc2km.py
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[2]
HASH_DIR = PROJECT / "01-scripts" / "HASH"
DATA_DIR = PROJECT / "02-data"


@dataclass(frozen=True)
class HashHeader:
    year: int
    month: int
    day: int
    hour: int
    minute: int
    second: float
    eh: float
    ez: float
    mag: float
    cluster_idx: int


def deg_to_degmin(decimal_deg: float) -> tuple[int, float]:
    d = abs(float(decimal_deg))
    deg = int(d)
    mins = 60.0 * (d - deg)
    return deg, mins


def parse_hash_phase_header(line: str) -> HashHeader:
    """
    Parse HASH phase.dat fixed-width event header line written by our pipeline.

    Format (no separators, fixed widths):
      year(4) mo(2) da(2) hr(2) mn(2) sec(5)
      lat_deg(2) ns(1) lat_min(5)
      lon_deg(3) ew(1) lon_min(5)
      depth(5) eh(6) ez(6) mag(5) cluster_idx(16) space newline
    """
    if len(line) < 72:
        raise ValueError(f"Header line too short: {len(line)}")

    return HashHeader(
        year=int(line[0:4]),
        month=int(line[4:6]),
        day=int(line[6:8]),
        hour=int(line[8:10]),
        minute=int(line[10:12]),
        second=float(line[12:17]),
        eh=float(line[39:45]),
        ez=float(line[45:51]),
        mag=float(line[51:56]),
        cluster_idx=int(line[56:72]),
    )


def format_hash_phase_header(h: HashHeader, *, lat: float, lon: float, depth_km: float) -> str:
    ilat, mlat = deg_to_degmin(lat)
    ilon, mlon = deg_to_degmin(lon)
    ns = "N" if lat >= 0 else "S"
    ew = "W" if lon < 0 else "E"
    # HASH expects positive degrees + E/W flag; minutes always positive.
    ilon = int(abs(int(ilon)))
    ilat = int(abs(int(ilat)))

    return (
        f"{h.year:4d}{h.month:2d}{h.day:2d}{h.hour:2d}{h.minute:2d}{h.second:5.2f}"
        f"{ilat:2d}{ns}{mlat:5.2f}{ilon:3d}{ew}{mlon:5.2f}"
        f"{depth_km:5.2f}{h.eh:6.2f}{h.ez:6.2f}{h.mag:5.2f}"
        f"{h.cluster_idx:>16} \n"
    )


def load_eid_map(path: Path) -> dict[int, int]:
    out: dict[int, int] = {}
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            out[int(row["cluster_idx"])] = int(row["event_id"])
    return out


def load_reloc_locations(path: Path) -> dict[int, tuple[float, float, float]]:
    """
    Load Felix relocation file into event_id -> (lat, lon, depth_km).

    Both Y1 and Y2 files share the first 4 columns:
      event_id, lat, lon, depth_km, ...
    """
    out: dict[int, tuple[float, float, float]] = {}
    with path.open() as f:
        for line in f:
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                eid = int(parts[0])
                lat = float(parts[1])
                lon = float(parts[2])
                dep = float(parts[3])
            except ValueError:
                continue
            out[eid] = (lat, lon, dep)
    return out


def update_phase_file(
    phase_in: Path,
    phase_out: Path,
    *,
    eid_by_cluster: dict[int, int],
    reloc_by_eid: dict[int, tuple[float, float, float]],
) -> tuple[int, int, int]:
    n_total = n_updated = n_missing = 0
    phase_out.parent.mkdir(parents=True, exist_ok=True)

    with phase_in.open() as fin, phase_out.open("w") as fout:
        for line in fin:
            if line and line[0].isdigit():
                n_total += 1
                try:
                    h = parse_hash_phase_header(line)
                except Exception:
                    fout.write(line)
                    continue

                eid = eid_by_cluster.get(h.cluster_idx)
                if eid is None:
                    n_missing += 1
                    fout.write(line)
                    continue

                loc = reloc_by_eid.get(eid)
                if loc is None:
                    n_missing += 1
                    fout.write(line)
                    continue

                lat, lon, dep = loc
                dep = max(float(dep), 0.01)
                fout.write(format_hash_phase_header(h, lat=lat, lon=lon, depth_km=dep))
                n_updated += 1
            else:
                fout.write(line)

    return n_total, n_updated, n_missing


def write_hash_input_with_new_phase(
    hash_input_in: Path,
    hash_input_out: Path,
    *,
    new_phase_filename: str,
    new_hashout_prefix: str,
) -> None:
    lines = hash_input_in.read_text().splitlines()
    if len(lines) < 9:
        raise ValueError(f"Unexpected HASH input file format: {hash_input_in}")

    # Line 5 (1-based) is the phase file.
    lines[4] = new_phase_filename
    # Lines 6–9 (1-based) are output filenames.
    lines[5] = f"{new_hashout_prefix}_1.dat"
    lines[6] = f"{new_hashout_prefix}_2.dat"
    lines[7] = f"{new_hashout_prefix}_3.dat"
    lines[8] = f"{new_hashout_prefix}_4.dat"

    hash_input_out.write_text("\n".join(lines) + "\n")


def main() -> None:
    reloc_dir = DATA_DIR / "A_all" / "reloc2km"
    if not reloc_dir.exists():
        raise SystemExit(f"Missing reloc dir: {reloc_dir}")

    jobs = [
        {
            "tag": "Y1_dl_all",
            "reloc": reloc_dir / "Y1.nlldd.2km.reloc",
            "eid_map": DATA_DIR / "G_FM" / "hash_Y1_dl_all_eid_map.csv",
            "phase_in": HASH_DIR / "phase_Y1_dl_all.dat",
            "phase_out": HASH_DIR / "phase_Y1_dl_all_reloc2km.dat",
            "hash_in": HASH_DIR / "hash.input_Y1_dl_all",
            "hash_out": HASH_DIR / "hash.input_Y1_dl_all_reloc2km",
            "hashout_prefix": "hashout_Y1_dl_all_reloc2km",
        },
        {
            "tag": "Y2_dl_all",
            "reloc": reloc_dir / "Y2.nlldd.2km.reloc",
            "eid_map": DATA_DIR / "G_FM" / "hash_Y2_dl_all_eid_map.csv",
            "phase_in": HASH_DIR / "phase_Y2_dl_all.dat",
            "phase_out": HASH_DIR / "phase_Y2_dl_all_reloc2km.dat",
            "hash_in": HASH_DIR / "hash.input_Y2_dl_all",
            "hash_out": HASH_DIR / "hash.input_Y2_dl_all_reloc2km",
            "hashout_prefix": "hashout_Y2_dl_all_reloc2km",
        },
    ]

    for job in jobs:
        for key in ("reloc", "eid_map", "phase_in", "hash_in"):
            p = Path(job[key])
            if not p.exists():
                raise SystemExit(f"Missing required input for {job['tag']}: {p}")

        eid_by_cluster = load_eid_map(Path(job["eid_map"]))
        reloc_by_eid = load_reloc_locations(Path(job["reloc"]))

        n_total, n_updated, n_missing = update_phase_file(
            Path(job["phase_in"]),
            Path(job["phase_out"]),
            eid_by_cluster=eid_by_cluster,
            reloc_by_eid=reloc_by_eid,
        )
        print(
            f"{job['tag']}: headers={n_total:,} updated={n_updated:,} missing_loc_or_id={n_missing:,} -> {Path(job['phase_out']).name}"
        )

        write_hash_input_with_new_phase(
            Path(job["hash_in"]),
            Path(job["hash_out"]),
            new_phase_filename=Path(job["phase_out"]).name,
            new_hashout_prefix=str(job["hashout_prefix"]),
        )
        print(f"{job['tag']}: wrote {Path(job['hash_out']).name}")


if __name__ == "__main__":
    main()

