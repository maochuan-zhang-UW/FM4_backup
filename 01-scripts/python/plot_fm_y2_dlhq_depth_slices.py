"""
plot_fm_y2_dlhq_depth_slices.py
================================
Depth-slice beach-ball maps for Y2 DL-hq (SKHASH 1D), Q=A/B only.
Depth bins: 0–0.25, 0.25–0.50, ..., 2.25–2.50 km  (10 panels, 2×5 grid).
All 21 Y2 stations shown in every panel.

Output: 02-data/G_FM/fm_y2_dlhq_depth_slices.png

Run with FM_ML env:
  /opt/miniconda3/envs/FM_ML/bin/python plot_fm_y2_dlhq_depth_slices.py
"""
import warnings; warnings.filterwarnings('ignore')

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from obspy.imaging.beachball import beach

from fault_classification import classify_fault_ptb

# ── Paths ──────────────────────────────────────────────────────────────────────
SKHASH_OUT = Path('/Users/mczhang/Documents/GitHub/SKHASH/SKHASH7/examples/hash3'
                  '/OUT_Y2_dl_hq_sk/out.txt')
OUT_DIR    = Path('/Users/mczhang/Documents/GitHub/FM7/02-data/G_FM')

# ── Map constants ──────────────────────────────────────────────────────────────
LON_LIM  = (-130.07, -129.93)   # wider to fit all stations
LAT_LIM  = (45.888,  46.030)
BB_WIDTH = 0.0014

FAULT_COLOR = {
    'N': [0.15, 0.25, 0.85],
    'R': [0.85, 0.15, 0.15],
    'S': [0.10, 0.70, 0.20],
    'U': [0.40, 0.40, 0.40],
}

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

# All 21 Y2 stations (from axial.Y2.sta)
STATIONS_Y2 = {
    'AS1': (45.93356, -129.99920), 'AS2': (45.93377, -130.01410),
    'CC1': (45.95468, -130.00890), 'EC1': (45.94958, -129.97970),
    'EC2': (45.93967, -129.97380), 'EC3': (45.93607, -129.97850),
    'ID1': (45.92573, -129.97800),
    '01B': (46.01935, -130.00526), '02B': (46.00111, -130.02930),
    '03B': (45.99250, -129.99397), '04B': (45.98834, -130.04896),
    '05B': (45.97327, -129.97728), '06B': (45.98178, -130.01423),
    '07B': (45.96950, -130.02963), '08B': (45.96336, -130.00441),
    '09B': (45.97194, -130.05972), '10B': (45.95867, -129.94982),
    '11B': (45.94851, -130.03801), '12B': (45.91473, -130.02173),
    '13B': (45.90472, -129.97240), '14B': (45.89913, -130.00612),
}

# ── Load & filter ──────────────────────────────────────────────────────────────
def classify_fault(rake):
    r = rake % 360
    if r > 180: r -= 360
    if -120 <= r <= -60:  return 'N'
    if   60 <= r <=  120: return 'R'
    if abs(r) <= 30 or abs(r) >= 150: return 'S'
    return 'U'

df = pd.read_csv(SKHASH_OUT)
df.columns = df.columns.str.strip()
df = df[df['quality'].isin(['A', 'B'])].copy()
df = df[(df['origin_lat'].between(*LAT_LIM)) &
        (df['origin_lon'].between(*LON_LIM))].copy()
df['fault_type'] = [
    classify_fault_ptb(s, d, r)
    for s, d, r in zip(df['strike'], df['dip'], df['rake'])
]
df['qrank']      = df['quality'].map({'A': 4, 'B': 3})
df = df.sort_values('qrank')   # A on top
print(f"Y2 DL-hq Q=A/B in map bounds: {len(df)} events")

# Depth bins: 0–0.25, 0.25–0.50, …, 2.25–2.50
edges = np.arange(0, 2.51, 0.25)
bins  = [(edges[i], edges[i+1]) for i in range(len(edges)-1)]   # 10 bins

# ── Figure layout: 2 rows × 5 cols ────────────────────────────────────────────
lon_range = abs(LON_LIM[1] - LON_LIM[0])
lat_range = abs(LAT_LIM[1] - LAT_LIM[0])
panel_w   = 3.2
panel_h   = panel_w * lat_range / lon_range

n_cols, n_rows = 5, 2
fig, axes = plt.subplots(n_rows, n_cols,
                          figsize=(panel_w * n_cols + 0.6, panel_h * n_rows + 1.2),
                          facecolor='white',
                          gridspec_kw={'hspace': 0.40, 'wspace': 0.15})

sta_lons = [v[1] for v in STATIONS_Y2.values()]
sta_lats = [v[0] for v in STATIONS_Y2.values()]
# OO stations (first 7) get a different marker
oo_lons = sta_lons[:7];  oo_lats = sta_lats[:7]
tf_lons = sta_lons[7:];  tf_lats = sta_lats[7:]

for idx, (d0, d1) in enumerate(bins):
    row, col = divmod(idx, n_cols)
    ax = axes[row, col]

    sub = df[(df['origin_depth_km'] >= d0) & (df['origin_depth_km'] < d1)]

    ax.set_facecolor([0.85, 0.92, 0.97])
    ax.set_aspect('equal')

    # Beach balls
    for _, row_data in sub.iterrows():
        try:
            bb = beach([row_data['strike'], row_data['dip'], row_data['rake']],
                       xy=(row_data['origin_lon'], row_data['origin_lat']),
                       width=BB_WIDTH, linewidth=0.15,
                       facecolor=FAULT_COLOR[row_data['fault_type']], alpha=0.85)
            bb.set_transform(ax.transData)
            bb.set_zorder(2)
            ax.add_collection(bb)
        except Exception:
            pass

    # Caldera rim
    ax.plot(CALDERA_RIM[:, 0], CALDERA_RIM[:, 1], '-',
            color='black', linewidth=1.0, zorder=7)

    # OO stations (triangles, yellow)
    ax.plot(oo_lons, oo_lats, '^', color='yellow', markeredgecolor='k',
            markeredgewidth=0.5, markersize=5, zorder=9, label='OO')
    # TF/2F stations (inverse triangles, cyan)
    ax.plot(tf_lons, tf_lats, 'v', color='cyan', markeredgecolor='k',
            markeredgewidth=0.5, markersize=5, zorder=9, label='2F')

    ax.set_xlim(LON_LIM); ax.set_ylim(LAT_LIM)
    ticks = np.arange(-130.06, -129.92, 0.04)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f'{t:.2f}' for t in ticks], fontsize=5.5, rotation=45)
    ax.tick_params(axis='y', labelsize=5.5)
    ax.grid(True, linewidth=0.25, color='gray', alpha=0.5)

    if col == 0:
        ax.set_ylabel('Latitude', fontsize=6.5)
    if row == n_rows - 1:
        ax.set_xlabel('Longitude', fontsize=6.5)

    ax.set_title(f'{d0:.2f}–{d1:.2f} km  (n={len(sub)})', fontsize=8, pad=3)

# Legend in last panel
legend_handles = [
    mpatches.Patch(color=FAULT_COLOR['N'], label='Normal'),
    mpatches.Patch(color=FAULT_COLOR['R'], label='Reverse'),
    mpatches.Patch(color=FAULT_COLOR['S'], label='Strike-slip'),
    mpatches.Patch(color=FAULT_COLOR['U'], label='Oblique'),
    plt.Line2D([0],[0], marker='^', color='w', markerfacecolor='yellow',
               markeredgecolor='k', markersize=6, label='OO stations'),
    plt.Line2D([0],[0], marker='v', color='w', markerfacecolor='cyan',
               markeredgecolor='k', markersize=6, label='2F stations'),
]
axes[-1, -1].legend(handles=legend_handles, loc='lower right',
                    fontsize=6.5, framealpha=0.9)

fig.suptitle('Y2 DL-hq Focal Mechanisms (SKHASH 1D)  |  Q = A+B only  |  '
             'Depth slices 0–2.5 km  |  ▲ OO  ▽ 2F stations',
             fontsize=11, y=1.01)
plt.subplots_adjust(top=0.95, bottom=0.06, left=0.05, right=0.99)

out = OUT_DIR / 'fm_y2_dlhq_depth_slices.png'
fig.savefig(out, dpi=150, bbox_inches='tight')
plt.close(fig)
print(f'Saved → {out}')
