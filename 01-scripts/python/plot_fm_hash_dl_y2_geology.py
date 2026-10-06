"""
plot_fm_hash_dl_y2_geology.py
==============================
Y2 HASH DL-hq Q=A/B focal mechanisms with geology overlays:
  fissures (1998/2011/2015), lava flows (1998/2011/2015), caldera rim, all stations.

Output: 02-data/G_FM/fm_hash_dl_y2_geology.png  (dpi=300)

Run with FM_ML env:
  /opt/miniconda3/envs/FM_ML/bin/python plot_fm_hash_dl_y2_geology.py
"""
import sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0, '.')
from load_axial_geology import load_fissures, load_lavas

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from obspy.imaging.beachball import beach

from fault_classification import classify_fault_ptb

# ── Paths ──────────────────────────────────────────────────────────────────────
HASH_DIR = Path('/Users/mczhang/Documents/GitHub/FM7/01-scripts/HASH')
OUT_DIR  = Path('/Users/mczhang/Documents/GitHub/FM7/02-data/G_FM')

# ── Map bounds ─────────────────────────────────────────────────────────────────
LON_LIM  = (-130.07, -129.93)
LAT_LIM  = (45.888,  46.025)
BB_WIDTH = 0.0015

# ── Colors ─────────────────────────────────────────────────────────────────────
FAULT_COLOR = {
    'N': [0.15, 0.25, 0.85],
    'R': [0.85, 0.15, 0.15],
    'S': [0.10, 0.70, 0.20],
    'U': [0.40, 0.40, 0.40],
}
LAVA_COLOR   = {'1998': [0.55, 0.10, 0.10],
                '2011': [0.10, 0.20, 0.75],
                '2015': [0.10, 0.55, 0.15]}
LAVA_ALPHA   = 0.35
FISS_COLOR   = {'1998': '#8B4513', '2011': '#4169E1', '2015': '#228B22'}

CALDERA_RIM = np.array([
    [-130.004785563058, 45.9207755734405],[-130.010476202888, 45.9238241104543],
    [-130.018881564079, 45.9351908809594],[-130.023946125193, 45.9412238501725],
    [-130.028718653506, 45.9498812001140],[-130.030451219380, 45.9511765797916],
    [-130.030679485650, 45.9542732243167],[-130.031733279709, 45.9558130656063],
    [-130.031444653500, 45.9586760104296],[-130.036188782208, 45.9656647517656],
    [-130.036950110789, 45.9698291665232],[-130.039953347000, 45.9750458167927],
    [-130.038595675479, 45.9847117727418],[-130.035927416999, 45.9883113986506],
    [-130.018067675296, 45.9933582886740],[-130.013629193751, 45.9937552841350],
    [-130.010365710979, 45.9929499241491],[-130.008647442296, 45.9924883829037],
    [-130.007262470669, 45.9915471582374],[-130.006042022411, 45.9902469280907],
    [-130.005178629490, 45.9897778053610],[-130.001868199523, 45.9863506519894],
    [-130.001154359192, 45.9846883853932],[-130.000949059432, 45.9827833001814],
    [-129.999393534330, 45.9818434725493],[-129.997797388662, 45.9786395525337],
    [-129.995357566829, 45.9760388622191],[-129.993956176267, 45.9741441737512],
    [-129.993678708114, 45.9681875631427],[-129.993140494256, 45.9667620754035],
    [-129.992087550788, 45.9652218741086],[-129.991186410747, 45.9626077204113],
    [-129.989604036931, 45.9601186407320],[-129.989238986151, 45.9588108137369],
    [-129.989728217453, 45.9574955894078],[-129.985484098670, 45.9494279735802],
    [-129.984788122490, 45.9487188881587],
])

STA_OO = {'AS1':(45.93356,-129.99920),'AS2':(45.93377,-130.01410),
           'CC1':(45.95468,-130.00890),'EC1':(45.94958,-129.97970),
           'EC2':(45.93967,-129.97380),'EC3':(45.93607,-129.97850),
           'ID1':(45.92573,-129.97800)}
STA_2F = {'01B':(46.01935,-130.00526),'02B':(46.00111,-130.02930),
           '03B':(45.99250,-129.99397),'04B':(45.98834,-130.04896),
           '05B':(45.97327,-129.97728),'06B':(45.98178,-130.01423),
           '07B':(45.96950,-130.02963),'08B':(45.96336,-130.00441),
           '09B':(45.97194,-130.05972),'10B':(45.95867,-129.94982),
           '11B':(45.94851,-130.03801),'12B':(45.91473,-130.02173),
           '13B':(45.90472,-129.97240),'14B':(45.89913,-130.00612)}

# ── Load geology ───────────────────────────────────────────────────────────────
print('Loading geology...')
fissures = load_fissures()
lavas    = load_lavas()

# ── Load Y2 HASH DL-hq A/B ────────────────────────────────────────────────────
def classify_fault(rake):
    r = rake % 360
    if r > 180: r -= 360
    if -120 <= r <= -60:  return 'N'
    if   60 <= r <=  120: return 'R'
    if abs(r) <= 30 or abs(r) >= 150: return 'S'
    return 'U'

rows = []
with open(HASH_DIR / 'hashout_Y2_dl_hq_1.dat') as fh:
    for line in fh:
        p = line.split()
        if len(p) < 31: continue
        try:
            rows.append(dict(lat=float(p[10]), lon=float(p[11]),
                             strike=float(p[21]), dip=float(p[22]),
                             rake=float(p[23]), quality=p[28]))
        except: continue
df = pd.DataFrame(rows)
df = df[df['quality'].isin(['A','B'])].copy()
df = df[(df['lat'].between(*LAT_LIM)) & (df['lon'].between(*LON_LIM))]
df['fault_type'] = [
    classify_fault_ptb(s, d, r)
    for s, d, r in zip(df['strike'], df['dip'], df['rake'])
]
df['qrank'] = df['quality'].map({'A':4,'B':3})
df = df.sort_values('qrank')
print(f'Y2 HASH DL-hq A/B in bounds: {len(df)}')

# ── Figure ─────────────────────────────────────────────────────────────────────
lon_r = abs(LON_LIM[1]-LON_LIM[0])
lat_r = abs(LAT_LIM[1]-LAT_LIM[0])
pw = 9.0; ph = pw * lat_r / lon_r

fig, ax = plt.subplots(figsize=(pw, ph+1.2), facecolor='white')
ax.set_facecolor([0.88, 0.94, 0.98])
ax.set_aspect('equal')

# ── Draw lava flows (bottom layer) ────────────────────────────────────────────
for yr, polys in lavas.items():
    fc = LAVA_COLOR[yr]
    for poly in polys:
        lons = poly[:,0]; lats = poly[:,1]
        mask = ((lons >= LON_LIM[0]) & (lons <= LON_LIM[1]) &
                (lats >= LAT_LIM[0]) & (lats <= LAT_LIM[1]))
        if mask.sum() < 3: continue
        ax.fill(lons, lats, color=fc, alpha=LAVA_ALPHA,
                edgecolor='none', zorder=1)

# ── Draw fissures ─────────────────────────────────────────────────────────────
for yr, segs in fissures.items():
    fc = FISS_COLOR[yr]
    for seg in segs:
        lons = seg[:,0]; lats = seg[:,1]
        mask = ((lons >= LON_LIM[0]) & (lons <= LON_LIM[1]) &
                (lats >= LAT_LIM[0]) & (lats <= LAT_LIM[1]))
        if mask.sum() < 2: continue
        ax.plot(lons, lats, '-', color=fc, linewidth=0.8, alpha=0.85, zorder=2)

# ── Beach balls ───────────────────────────────────────────────────────────────
for _, r in df.iterrows():
    try:
        bb = beach([r['strike'], r['dip'], r['rake']],
                   xy=(r['lon'], r['lat']), width=BB_WIDTH,
                   linewidth=0.2, facecolor=FAULT_COLOR[r['fault_type']], alpha=0.88)
        bb.set_transform(ax.transData); bb.set_zorder(4)
        ax.add_collection(bb)
    except: pass

# ── Caldera rim ───────────────────────────────────────────────────────────────
ax.plot(CALDERA_RIM[:,0], CALDERA_RIM[:,1], '-k', linewidth=1.8, zorder=6)

# ── Stations ──────────────────────────────────────────────────────────────────
ax.plot([v[1] for v in STA_OO.values()], [v[0] for v in STA_OO.values()],
        '^', color='yellow', markeredgecolor='k', markeredgewidth=0.7,
        markersize=8, zorder=8)
ax.plot([v[1] for v in STA_2F.values()], [v[0] for v in STA_2F.values()],
        'v', color='cyan', markeredgecolor='k', markeredgewidth=0.7,
        markersize=8, zorder=8)

# ── Axes ──────────────────────────────────────────────────────────────────────
ax.set_xlim(LON_LIM); ax.set_ylim(LAT_LIM)
ticks = np.arange(-130.06, -129.92, 0.04)
ax.set_xticks(ticks)
ax.set_xticklabels([f'{t:.2f}' for t in ticks], fontsize=9, rotation=45)
ax.tick_params(axis='y', labelsize=9)
ax.grid(True, lw=0.3, color='gray', alpha=0.4)
ax.set_xlabel('Longitude', fontsize=10)
ax.set_ylabel('Latitude', fontsize=10)
ax.set_title(f'Y2 DL-hq Focal Mechanisms (HASH 1D)  |  Q = A+B only  |  '
             f'n = {len(df)}  |  with fissures & lava flows',
             fontsize=11, pad=6)

# ── Legend ────────────────────────────────────────────────────────────────────
legend_handles = [
    # FM fault types
    Patch(color=FAULT_COLOR['N'], label='Normal FM'),
    Patch(color=FAULT_COLOR['R'], label='Reverse FM'),
    Patch(color=FAULT_COLOR['S'], label='Strike-slip FM'),
    Patch(color=FAULT_COLOR['U'], label='Oblique FM'),
    # Lava flows
    Patch(facecolor=LAVA_COLOR['1998'], alpha=0.7, label='Lava 1998'),
    Patch(facecolor=LAVA_COLOR['2011'], alpha=0.7, label='Lava 2011'),
    Patch(facecolor=LAVA_COLOR['2015'], alpha=0.7, label='Lava 2015'),
    # Fissures
    Line2D([0],[0], color=FISS_COLOR['1998'], lw=1.5, label='Fissures 1998'),
    Line2D([0],[0], color=FISS_COLOR['2011'], lw=1.5, label='Fissures 2011'),
    Line2D([0],[0], color=FISS_COLOR['2015'], lw=1.5, label='Fissures 2015'),
    # Stations
    Line2D([0],[0], marker='^', color='w', markerfacecolor='yellow',
           markeredgecolor='k', markersize=9, label='OO stations'),
    Line2D([0],[0], marker='v', color='w', markerfacecolor='cyan',
           markeredgecolor='k', markersize=9, label='2F stations'),
]
ax.legend(handles=legend_handles, loc='lower right',
          fontsize=7.5, framealpha=0.92, ncol=2)

plt.tight_layout()
out = OUT_DIR / 'fm_hash_dl_y2_geology.png'
fig.savefig(out, dpi=300, bbox_inches='tight')
plt.close(fig)
print(f'Saved → {out}')
