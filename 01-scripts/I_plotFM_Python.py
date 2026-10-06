"""
I_plotFM_Python.py  —  Plot focal mechanisms on Axial Seamount basemap.

Replicates basemap_2015v2.m:
  - Caldera rim (Chadwick)
  - 2015 / 2011 / 1998 eruptive fissures
  - OBS station triangles
  - Geographic aspect ratio: pbaspect([diff(lon)*cos(lat), diff(lat), 1])
  - Only A and B quality focal mechanisms

Usage:
    python I_plotFM_Python.py                          # plots G_3D.mat (SKHASH)
    python I_plotFM_Python.py G_2F_HASH_DL_TMSF001    # plots DL HASH results
"""

import sys, os
import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from obspy.imaging.beachball import beach

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE     = "/Users/mczhang/Documents/GitHub/FM4"
FMDIR    = f"{BASE}/02-data/G_FM"
PNG_DIR  = f"{BASE}/03-output-graphics/MyPng"
PDF_DIR  = f"{BASE}/03-output-graphics/MyPdf"
DATA_DIR = "/Users/mczhang/Documents/GitHub/FM6_RealTime/01-scripts/subcode/data"

stem    = sys.argv[1] if len(sys.argv) > 1 else "G_3D"
matfile = f"{FMDIR}/{stem}.mat"

# ── Map limits (match I_plotFM_SKHASH.m / I_plotFM_HASH.m) ───────────────────
LON_LIM = [-130.04, -129.97]
LAT_LIM = [ 45.908,  46.001]

mean_lat = np.mean(LAT_LIM)
# Geographic aspect: replicate pbaspect([diff(lon)*cos(lat), diff(lat), 1])
# matplotlib set_aspect = dy_display / dx_display per data unit = 1/cos(lat)
GEO_ASPECT = 1.0 / np.cos(np.radians(mean_lat))

# ── Caldera rim (from axial_calderaRim.m) ─────────────────────────────────────
CALDERA_RIM = np.array([
    [-130.004785563058, 45.9207755734405],
    [-130.010476202888, 45.9238241104543],
    [-130.018881564079, 45.9351908809594],
    [-130.023946125193, 45.9412238501725],
    [-130.028718653506, 45.949881200114 ],
    [-130.03045121938,  45.9511765797916],
    [-130.03067948565,  45.9542732243167],
    [-130.031733279709, 45.9558130656063],
    [-130.0314446535,   45.9586760104296],
    [-130.036188782208, 45.9656647517656],
    [-130.036950110789, 45.9698291665232],
    [-130.039953347,    45.9750458167927],
    [-130.038595675479, 45.9847117727418],
    [-130.035927416999, 45.9883113986506],
    [-130.018067675296, 45.993358288674 ],
    [-130.013629193751, 45.993755284135 ],
    [-130.010365710979, 45.9929499241491],
    [-130.008647442296, 45.9924883829037],
    [-130.007262470669, 45.9915471582374],
    [-130.006042022411, 45.9902469280907],
    [-130.00517862949,  45.989777805361 ],
    [-130.001868199523, 45.9863506519894],
    [-130.001154359192, 45.9846883853932],
    [-130.000949059432, 45.9827833001814],
    [-129.99939353433,  45.9818434725493],
    [-129.997797388662, 45.9786395525337],
    [-129.995357566829, 45.9760388622191],
    [-129.993956176267, 45.9741441737512],
    [-129.993678708114, 45.9681875631427],
    [-129.993140494256, 45.9667620754035],
    [-129.992087550788, 45.9652218741086],
    [-129.991186410747, 45.9626077204113],
    [-129.989604036931, 45.960118640732 ],
    [-129.989238986151, 45.9588108137369],
    [-129.989728217453, 45.9574955894078],
    [-129.98548409867,  45.9494279735802],
    [-129.98478812249,  45.9487188881587],
])

# ── Stations (axial_stationsNewOrder.m, all 7 cabled OBS) ─────────────────────
STATIONS = [
    ("AXCC1", -130.0089,  45.95468),
    ("AXEC1", -129.9797,  45.94958),
    ("AXEC2", -129.9738,  45.93967),
    ("AXEC3", -129.9785,  45.93607),
    ("AXAS1", -129.9992,  45.93356),
    ("AXAS2", -130.0141,  45.93377),
    ("AXID1", -129.978,   45.92573),
]

# ── Fissure loader ────────────────────────────────────────────────────────────
def load_fissures_tab(path):
    """Load tab-delimited fissure file: ORIG_FID  LONGITUDE  LATITUDE."""
    segs = {}
    with open(path) as f:
        next(f)                         # skip header
        for line in f:
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            try:
                fid  = int(parts[0])
                lon  = float(parts[1])
                lat  = float(parts[2])
                segs.setdefault(fid, []).append((lon, lat))
            except ValueError:
                continue
    return segs


def load_fissures_space(path):
    """Load space/comma-delimited fissure file: id  lon  lat (or id lon,lat)."""
    segs = {}
    with open(path) as f:
        for line in f:
            line = line.strip().replace(',', ' ')
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                fid = int(float(parts[0]))
                lon = float(parts[1])
                lat = float(parts[2])
                segs.setdefault(fid, []).append((lon, lat))
            except ValueError:
                continue
    return segs


# ── FM data ───────────────────────────────────────────────────────────────────
print(f"Loading: {matfile}")
mat    = sio.loadmat(matfile, simplify_cells=True)
events = mat["event1"]
if not isinstance(events, list):
    events = list(np.asarray(events).ravel())
print(f"  {len(events)} events total.")

def safe_scalar(v, default=None):
    if v is None:
        return default
    a = np.asarray(v).ravel()
    return float(a[0]) if a.size > 0 else default

GOOD_QUAL = {"A", "B"}
FT_COLOR  = {
    "N": [0.12, 0.47, 0.71],    # blue  — Normal
    "R": [0.84, 0.15, 0.16],    # red   — Reverse
    "S": [0.17, 0.63, 0.17],    # green — Strike-slip
    "U": [0.50, 0.50, 0.50],    # grey  — Undefined
}

filtered = []
for ev in events:
    q  = str(ev.get("mechqual", "D")).strip()
    if q not in GOOD_QUAL:
        continue
    lat = safe_scalar(ev.get("lat"))
    lon = safe_scalar(ev.get("lon"))
    if lat is None or lon is None:
        continue
    if not (LAT_LIM[0] <= lat <= LAT_LIM[1] and LON_LIM[0] <= lon <= LON_LIM[1]):
        continue
    avmech = np.asarray(ev.get("avmech", [])).ravel()
    if avmech.size < 3:
        continue
    ft = str(ev.get("faultType", "U")).strip()
    filtered.append({"lat": lat, "lon": lon,
                     "strike": avmech[0], "dip": avmech[1], "rake": avmech[2],
                     "ft": ft, "qual": q})

print(f"  {len(filtered)} A/B quality events to plot.")

# Sort so A plots on top of B
filtered.sort(key=lambda e: {"A": 1, "B": 0}.get(e["qual"], 0))

# ── Figure & axes ─────────────────────────────────────────────────────────────
dlon = LON_LIM[1] - LON_LIM[0]   # 0.07°
dlat = LAT_LIM[1] - LAT_LIM[0]   # 0.093°

# Physical width : height = dlon*cos(lat) : dlat  →  same ratio for figure
fig_h = 9.0
fig_w = fig_h * (dlon * np.cos(np.radians(mean_lat))) / dlat
fig, ax = plt.subplots(figsize=(fig_w, fig_h))

ax.set_xlim(LON_LIM)
ax.set_ylim(LAT_LIM)
ax.set_aspect(GEO_ASPECT)          # = 1/cos(lat): geographic distortion correction

# ── Basemap elements ──────────────────────────────────────────────────────────
# 2011 fissures
try:
    segs2011 = load_fissures_tab(f"{DATA_DIR}/Axial-2011-fissures-points-geo-v2.txt")
    for pts in segs2011.values():
        xs, ys = zip(*pts)
        ax.plot(xs, ys, "-", color="dimgray", linewidth=0.8, zorder=2)
except Exception as e:
    print(f"  Warning: 2011 fissures not loaded ({e})")

# 2015 fissures
try:
    segs2015 = load_fissures_tab(f"{DATA_DIR}/Axial-2015-fissures-points-geo-v2.txt")
    for pts in segs2015.values():
        xs, ys = zip(*pts)
        ax.plot(xs, ys, "-", color="dimgray", linewidth=0.8, zorder=2)
except Exception as e:
    print(f"  Warning: 2015 fissures not loaded ({e})")

# 1998 fissures
try:
    segs1998 = load_fissures_space(f"{DATA_DIR}/Axial-1998-Fissures.txt")
    for pts in segs1998.values():
        xs, ys = zip(*pts)
        ax.plot(xs, ys, "-", color="dimgray", linewidth=0.8, zorder=2)
except Exception as e:
    print(f"  Warning: 1998 fissures not loaded ({e})")

# Caldera rim
ax.plot(CALDERA_RIM[:, 0], CALDERA_RIM[:, 1], "-k", linewidth=2.5, zorder=3)

# Stations
sta_lons = [s[1] for s in STATIONS]
sta_lats = [s[2] for s in STATIONS]
ax.plot(sta_lons, sta_lats, "^k",
        markerfacecolor="black", markersize=9, zorder=10)

# ── Beachballs ────────────────────────────────────────────────────────────────
# Replicate MATLAB radius = 0.0005, scale_event = 1.3  →  width ≈ 0.0013°
# Slightly increased to 0.0015° for Python rasterisation clarity
BB_WIDTH = 0.0015

for ev in filtered:
    color = FT_COLOR.get(ev["ft"], FT_COLOR["U"])
    try:
        b = beach([ev["strike"], ev["dip"], ev["rake"]],
                  xy=(ev["lon"], ev["lat"]),
                  width=BB_WIDTH,
                  linewidth=0.3,
                  facecolor=color,
                  zorder=6)
        ax.add_collection(b)
    except Exception:
        continue

# ── Labels & legend ───────────────────────────────────────────────────────────
ax.set_xlabel("Longitude (°)", fontsize=14)
ax.set_ylabel("Latitude (°)",  fontsize=14)
ax.set_title(f"{len(filtered)} Focal Mechanisms  [{stem}]\n"
             f"Quality A / B  |  N=blue  R=red  S=green", fontsize=13)
ax.grid(True, alpha=0.3, linewidth=0.5)

# x-tick spacing 0.01° (matches MATLAB set(ax,'XTick',-130.03:0.01:-129.97))
xticks = np.arange(np.ceil(LON_LIM[0] / 0.01) * 0.01,
                   LON_LIM[1] + 0.001, 0.01)
ax.set_xticks(xticks)
ax.set_xticklabels([f"{x:.2f}" for x in xticks], fontsize=9, rotation=30)
ax.tick_params(axis="y", labelsize=9)

legend_handles = [
    mpatches.Patch(color=FT_COLOR["N"], label="Normal"),
    mpatches.Patch(color=FT_COLOR["R"], label="Reverse"),
    mpatches.Patch(color=FT_COLOR["S"], label="Strike-slip"),
    mpatches.Patch(color=FT_COLOR["U"], label="Undefined"),
    Line2D([0], [0], marker="^", color="w", markerfacecolor="k",
           markersize=9, label="OBS station"),
    Line2D([0], [0], color="k", linewidth=2.5, label="Caldera rim"),
]
ax.legend(handles=legend_handles, loc="upper right",
          fontsize=9, framealpha=0.85, edgecolor="k")

# ── Save ──────────────────────────────────────────────────────────────────────
os.makedirs(PNG_DIR, exist_ok=True)
os.makedirs(PDF_DIR, exist_ok=True)
out_png = f"{PNG_DIR}/I_FM_{stem}.png"
out_pdf = f"{PDF_DIR}/I_FM_{stem}.pdf"
plt.tight_layout()
plt.savefig(out_png, dpi=150, bbox_inches="tight")
plt.savefig(out_pdf,           bbox_inches="tight")
print(f"Saved: {out_png}")
print(f"Saved: {out_pdf}")
