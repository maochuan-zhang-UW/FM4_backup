"""
plot_fm_cc_vs_dlhq.py
=====================
Single figure: 4 panels comparing CC-SVD vs DL-hq focal mechanisms
for Y1 (2022-23) and Y2 (2023-24), all depths combined. Q=A+B only.

Output: 02-data/G_FM/fm_cc_vs_dlhq.png

Run with FM_ML env:
  /opt/miniconda3/envs/FM_ML/bin/python plot_fm_cc_vs_dlhq.py
"""
import warnings
warnings.filterwarnings('ignore')

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from obspy.imaging.beachball import beach

from fault_classification import classify_fault_ptb

# ── Paths ──────────────────────────────────────────────────────────────────────
HASH_DIR = Path('/Users/mczhang/Documents/GitHub/FM7/01-scripts/HASH')
OUT_DIR  = Path('/Users/mczhang/Documents/GitHub/FM7/02-data/G_FM')
OUT_DIR.mkdir(exist_ok=True)

# ── Constants ──────────────────────────────────────────────────────────────────
LON_LIM     = (-130.04, -129.97)
LAT_LIM     = (45.908,  46.001)
BB_WIDTH    = 0.0012
QUAL_FILTER = ['A', 'B']
QUAL_ORDER  = {'A': 4, 'B': 3, 'C': 2, 'D': 1}

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

STA_OO = {
    'AS1': (45.93356, -129.99920), 'AS2': (45.93377, -130.01410),
    'CC1': (45.95468, -130.00890), 'EC1': (45.94958, -129.97970),
    'EC2': (45.93967, -129.97380), 'EC3': (45.93607, -129.97850),
    'ID1': (45.92573, -129.97800),
}

# ── Load and filter HASH output ────────────────────────────────────────────────
def classify_fault(rake):
    r = rake % 360
    if r > 180: r -= 360
    if -120 <= r <= -60:  return 'N'
    if   60 <= r <=  120: return 'R'
    if abs(r) <= 30 or abs(r) >= 150: return 'S'
    return 'U'

def load_hash(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            p = line.split()
            if len(p) < 31:
                continue
            try:
                rows.append(dict(
                    origin_lat      = float(p[10]),
                    origin_lon      = float(p[11]),
                    origin_depth_km = float(p[12]),
                    strike          = float(p[21]),
                    dip             = float(p[22]),
                    rake            = float(p[23]),
                    quality         = p[28],
                ))
            except (ValueError, IndexError):
                continue
    df = pd.DataFrame(rows)
    if len(df):
        df['fault_type'] = [
            classify_fault_ptb(s, d, r)
            for s, d, r in zip(df['strike'], df['dip'], df['rake'])
        ]
        df['qrank']      = df['quality'].map(QUAL_ORDER)
    return df

def prep(df):
    sub = df[df['quality'].isin(QUAL_FILTER)].copy()
    sub = sub[(sub['origin_lat'].between(*LAT_LIM)) &
              (sub['origin_lon'].between(*LON_LIM))]
    return sub.sort_values('qrank')

tags = {
    'Y1 CC-SVD':        'Y1_cc',
    'Y1 DL-hq (≥0.8)': 'Y1_dl_hq',
    'Y2 CC-SVD':        'Y2_cc',
    'Y2 DL-hq (≥0.8)': 'Y2_dl_hq',
}

datasets = {}
for label, tag in tags.items():
    p = HASH_DIR / f'hashout_{tag}_1.dat'
    if p.exists():
        datasets[label] = prep(load_hash(p))
    else:
        datasets[label] = pd.DataFrame()
    print(f'{label}: {len(datasets[label])} Q=A/B events')

# ── Figure ─────────────────────────────────────────────────────────────────────
lon_range = abs(LON_LIM[1] - LON_LIM[0])   # 0.07
lat_range = abs(LAT_LIM[1] - LAT_LIM[0])   # 0.093
panel_h   = 7.0
panel_w   = panel_h * (lon_range / lat_range)  # ~5.27 in

fig, axes = plt.subplots(1, 4,
                          figsize=(panel_w * 4 + 0.8, panel_h),
                          facecolor='white')

for ax, (label, sub) in zip(axes, datasets.items()):
    ax.set_facecolor([0.85, 0.92, 0.97])
    ax.set_aspect('equal')

    for _, row in sub.iterrows():
        try:
            bb = beach([row['strike'], row['dip'], row['rake']],
                       xy=(row['origin_lon'], row['origin_lat']),
                       width=BB_WIDTH, linewidth=0.15,
                       facecolor=FAULT_COLOR[row['fault_type']], alpha=0.85)
            bb.set_transform(ax.transData)
            bb.set_zorder(2)
            ax.add_collection(bb)
        except Exception:
            pass

    ax.plot(CALDERA_RIM[:, 0], CALDERA_RIM[:, 1], '-',
            color='black', linewidth=1.2, zorder=7)
    ax.plot([v[1] for v in STA_OO.values()],
            [v[0] for v in STA_OO.values()],
            '^', color='yellow', markeredgecolor='k',
            markeredgewidth=0.4, markersize=5, zorder=8)

    ax.set_xlim(LON_LIM)
    ax.set_ylim(LAT_LIM)
    ticks = np.arange(-130.04, -129.96, 0.02)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f'{t:.2f}' for t in ticks], fontsize=6.5, rotation=45)
    ax.tick_params(axis='y', labelsize=6.5)
    ax.grid(True, linewidth=0.3, color='gray', alpha=0.5)
    ax.set_xlabel('Longitude', fontsize=8)
    ax.set_title(f'{label}\nn={len(sub)}', fontsize=9, pad=4)

axes[0].set_ylabel('Latitude', fontsize=8)

legend_handles = [
    mpatches.Patch(color=FAULT_COLOR['N'], label='Normal'),
    mpatches.Patch(color=FAULT_COLOR['R'], label='Reverse'),
    mpatches.Patch(color=FAULT_COLOR['S'], label='Strike-slip'),
    mpatches.Patch(color=FAULT_COLOR['U'], label='Oblique'),
]
axes[-1].legend(handles=legend_handles, loc='lower right', fontsize=7, framealpha=0.9)

fig.suptitle('HASH Single-Event FMs  |  CC-SVD vs DL-hq (conf ≥ 0.8)  |  '
             'Q=A+B  |  All depths',
             fontsize=11, y=1.01)
plt.tight_layout()

out = OUT_DIR / 'fm_cc_vs_dlhq.png'
fig.savefig(out, dpi=150, bbox_inches='tight')
plt.close(fig)
print(f'\nSaved → {out}')
