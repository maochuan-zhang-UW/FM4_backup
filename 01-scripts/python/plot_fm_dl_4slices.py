"""
plot_fm_dl_4slices.py
=====================
DL-hq FMs (SKHASH 1D) Y1+Y2 combined, 4 depth slices.
Depth bins: 0–1, 1–1.5, 1.5–2, 2–end km.
All qualities: A/B = full colour, C/D = medium-faded colour.
All 21 Y2 stations shown. Map extended to 46.025 N.

Output: 02-data/G_FM/fm_dl_4slices.png

Run with FM_ML env:
  /opt/miniconda3/envs/FM_ML/bin/python plot_fm_dl_4slices.py
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
HASH3_DIR = Path('/Users/mczhang/Documents/GitHub/SKHASH/SKHASH7/examples/hash3')
OUT_DIR   = Path('/Users/mczhang/Documents/GitHub/FM7/02-data/G_FM')

# ── Constants ──────────────────────────────────────────────────────────────────
LON_LIM  = (-130.07, -129.93)   # wide enough to show all stations
LAT_LIM  = (45.888,  46.025)    # extended north to 46.025 to show 01B
BB_WIDTH = 0.0015

# A/B: full colour; C/D: medium-faded (more visible than before)
FAULT_COLOR = {
    'N': {'AB': [0.15, 0.25, 0.85], 'CD': [0.55, 0.65, 0.93]},
    'R': {'AB': [0.85, 0.15, 0.15], 'CD': [0.93, 0.50, 0.50]},
    'S': {'AB': [0.10, 0.70, 0.20], 'CD': [0.50, 0.82, 0.55]},
    'U': {'AB': [0.40, 0.40, 0.40], 'CD': [0.65, 0.65, 0.65]},
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

# All 21 Y2 stations (OO + 2F)
STA_OO = {
    'AS1': (45.93356, -129.99920), 'AS2': (45.93377, -130.01410),
    'CC1': (45.95468, -130.00890), 'EC1': (45.94958, -129.97970),
    'EC2': (45.93967, -129.97380), 'EC3': (45.93607, -129.97850),
    'ID1': (45.92573, -129.97800),
}
STA_2F = {
    '01B': (46.01935, -130.00526), '02B': (46.00111, -130.02930),
    '03B': (45.99250, -129.99397), '04B': (45.98834, -130.04896),
    '05B': (45.97327, -129.97728), '06B': (45.98178, -130.01423),
    '07B': (45.96950, -130.02963), '08B': (45.96336, -130.00441),
    '09B': (45.97194, -130.05972), '10B': (45.95867, -129.94982),
    '11B': (45.94851, -130.03801), '12B': (45.91473, -130.02173),
    '13B': (45.90472, -129.97240), '14B': (45.89913, -130.00612),
}

DEPTH_BINS = [(0, 1), (1, 1.5), (1.5, 2), (2, 999)]
BIN_LABELS = ['0–1 km', '1–1.5 km', '1.5–2 km', '2 km – end']

# ── Helpers ────────────────────────────────────────────────────────────────────
def classify_fault(rake):
    r = rake % 360
    if r > 180: r -= 360
    if -120 <= r <= -60:  return 'N'
    if   60 <= r <=  120: return 'R'
    if abs(r) <= 30 or abs(r) >= 150: return 'S'
    return 'U'

def load_dl(tag):
    df = pd.read_csv(HASH3_DIR / f'OUT_{tag}' / 'out.txt')
    df.columns = df.columns.str.strip()
    df['fault_type'] = [
        classify_fault_ptb(s, d, r)
        for s, d, r in zip(df['strike'], df['dip'], df['rake'])
    ]
    df['ab'] = df['quality'].isin(['A', 'B'])
    df = df[(df['origin_lat'].between(*LAT_LIM)) &
            (df['origin_lon'].between(*LON_LIM))].copy()
    df['zorder_key'] = df['ab'].map({False: 0, True: 1})
    return df.sort_values('zorder_key')

# Combine Y1 + Y2
df_y1 = load_dl('Y1_dl_hq_sk')
df_y2 = load_dl('Y2_dl_hq_sk')
df_all = pd.concat([df_y1, df_y2], ignore_index=True).sort_values('zorder_key')

print(f"Y1: {df_y1['ab'].sum()} A/B + {(~df_y1['ab']).sum()} C/D = {len(df_y1)} total")
print(f"Y2: {df_y2['ab'].sum()} A/B + {(~df_y2['ab']).sum()} C/D = {len(df_y2)} total")
print(f"Combined: {df_all['ab'].sum()} A/B + {(~df_all['ab']).sum()} C/D = {len(df_all)} total")

# ── Figure: 1 row × 4 cols ────────────────────────────────────────────────────
lon_range = abs(LON_LIM[1] - LON_LIM[0])
lat_range = abs(LAT_LIM[1] - LAT_LIM[0])
panel_w   = 3.5
panel_h   = panel_w * lat_range / lon_range

fig, axes = plt.subplots(1, 4,
                          figsize=(panel_w * 4 + 0.6, panel_h + 1.0),
                          facecolor='white',
                          gridspec_kw={'wspace': 0.15})

oo_lons = [v[1] for v in STA_OO.values()]
oo_lats = [v[0] for v in STA_OO.values()]
tf_lons = [v[1] for v in STA_2F.values()]
tf_lats = [v[0] for v in STA_2F.values()]

ticks = np.arange(-130.06, -129.92, 0.04)

for col_j, ((d0, d1), bin_lbl) in enumerate(zip(DEPTH_BINS, BIN_LABELS)):
    ax = axes[col_j]
    sub = df_all[(df_all['origin_depth_km'] >= d0) & (df_all['origin_depth_km'] < d1)]

    ax.set_facecolor([0.85, 0.92, 0.97])
    ax.set_aspect('equal')

    for _, r in sub.iterrows():
        qk = 'AB' if r['ab'] else 'CD'
        try:
            bb = beach([r['strike'], r['dip'], r['rake']],
                       xy=(r['origin_lon'], r['origin_lat']),
                       width=BB_WIDTH,
                       linewidth=0.15 if r['ab'] else 0.10,
                       facecolor=FAULT_COLOR[r['fault_type']][qk],
                       alpha=0.85 if r['ab'] else 0.60)
            bb.set_transform(ax.transData)
            bb.set_zorder(2 if r['ab'] else 1)
            ax.add_collection(bb)
        except Exception:
            pass

    ax.plot(CALDERA_RIM[:, 0], CALDERA_RIM[:, 1], '-',
            color='black', linewidth=1.0, zorder=7)

    # OO stations: yellow triangles up
    ax.plot(oo_lons, oo_lats, '^', color='yellow', markeredgecolor='k',
            markeredgewidth=0.5, markersize=5, zorder=9)
    # 2F stations: cyan triangles down
    ax.plot(tf_lons, tf_lats, 'v', color='cyan', markeredgecolor='k',
            markeredgewidth=0.5, markersize=5, zorder=9)

    ax.set_xlim(LON_LIM); ax.set_ylim(LAT_LIM)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f'{t:.2f}' for t in ticks], fontsize=6, rotation=45)
    ax.tick_params(axis='y', labelsize=6)
    ax.grid(True, linewidth=0.25, color='gray', alpha=0.5)
    ax.set_xlabel('Longitude', fontsize=7)

    n_ab = sub['ab'].sum()
    n_cd = (~sub['ab']).sum()
    ax.set_title(f'{bin_lbl}\nA/B={n_ab}  C/D={n_cd}', fontsize=9, pad=4)

axes[0].set_ylabel('Latitude', fontsize=8)

# Legend
legend_handles = [
    mpatches.Patch(color=FAULT_COLOR['N']['AB'], label='Normal (A/B)'),
    mpatches.Patch(color=FAULT_COLOR['R']['AB'], label='Reverse (A/B)'),
    mpatches.Patch(color=FAULT_COLOR['S']['AB'], label='Strike-slip (A/B)'),
    mpatches.Patch(color=FAULT_COLOR['U']['AB'], label='Oblique (A/B)'),
    mpatches.Patch(color=FAULT_COLOR['N']['CD'], label='Normal (C/D)'),
    mpatches.Patch(color=FAULT_COLOR['R']['CD'], label='Reverse (C/D)'),
    mpatches.Patch(color=FAULT_COLOR['S']['CD'], label='Strike-slip (C/D)'),
    mpatches.Patch(color=FAULT_COLOR['U']['CD'], label='Oblique (C/D)'),
    plt.Line2D([0],[0], marker='^', color='w', markerfacecolor='yellow',
               markeredgecolor='k', markersize=7, label='OO stations'),
    plt.Line2D([0],[0], marker='v', color='w', markerfacecolor='cyan',
               markeredgecolor='k', markersize=7, label='2F stations'),
]
axes[-1].legend(handles=legend_handles, loc='lower right',
                fontsize=6.5, framealpha=0.9, ncol=2)

fig.suptitle('DL-hq Focal Mechanisms (SKHASH 1D)  |  Y1 + Y2 Combined  |  '
             '4 Depth Slices  |  Full colour = A/B,  Faded = C/D',
             fontsize=11, y=1.01)
plt.subplots_adjust(top=0.93, bottom=0.08, left=0.05, right=0.99)

out = OUT_DIR / 'fm_dl_4slices.png'
fig.savefig(out, dpi=150, bbox_inches='tight')
plt.close(fig)
print(f'\nSaved → {out}')
