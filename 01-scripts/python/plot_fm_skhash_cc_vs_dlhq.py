"""
plot_fm_skhash_cc_vs_dlhq.py
============================
4-panel beach-ball map: SKHASH CC-SVD vs DL-hq, Y1 and Y2.
Q=A+B only, all depths, 1D velocity model.

Output: 02-data/G_FM/fm_skhash_cc_vs_dlhq.png

Run with FM_ML env:
  /opt/miniconda3/envs/FM_ML/bin/python plot_fm_skhash_cc_vs_dlhq.py
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
SKHASH_DIR = Path('/Users/mczhang/Documents/GitHub/SKHASH/SKHASH7')
HASH3_DIR  = SKHASH_DIR / 'examples' / 'hash3'
OUT_DIR    = Path('/Users/mczhang/Documents/GitHub/FM7/02-data/G_FM')

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

# ── Load SKHASH output ─────────────────────────────────────────────────────────
def classify_fault(rake):
    r = rake % 360
    if r > 180: r -= 360
    if -120 <= r <= -60:  return 'N'
    if   60 <= r <=  120: return 'R'
    if abs(r) <= 30 or abs(r) >= 150: return 'S'
    return 'U'

def load_skhash(path):
    if not path.exists():
        print(f'  WARNING: {path} not found')
        return pd.DataFrame()
    df = pd.read_csv(path)
    df.columns = df.columns.str.strip()
    df['fault_type'] = [
        classify_fault_ptb(s, d, r)
        for s, d, r in zip(df['strike'], df['dip'], df['rake'])
    ]
    df['qrank']      = df['quality'].map(QUAL_ORDER)
    return df

def prep(df):
    if len(df) == 0: return df
    sub = df[df['quality'].isin(QUAL_FILTER)].copy()
    sub = sub[(sub['origin_lat'].between(*LAT_LIM)) &
              (sub['origin_lon'].between(*LON_LIM))]
    return sub.sort_values('qrank')

tags = {
    'Y1 CC-SVD (SKHASH)':        'OUT_Y1_cc_sk',
    'Y1 DL-hq ≥0.8 (SKHASH)':   'OUT_Y1_dl_hq_sk',
    'Y2 CC-SVD (SKHASH)':        'Y2_cc_sk',
    'Y2 DL-hq ≥0.8 (SKHASH)':   'Y2_dl_hq_sk',
}

datasets = {}
for label, tag in tags.items():
    out_path = HASH3_DIR / tag / 'out.txt'
    # fix tag key for Y2 to match actual directory
    if tag.startswith('Y2'):
        out_path = HASH3_DIR / f'OUT_{tag}' / 'out.txt'
    datasets[label] = prep(load_skhash(out_path))
    print(f'{label}: {len(datasets[label])} Q=A/B events')

# ── Figure ─────────────────────────────────────────────────────────────────────
lon_range = abs(LON_LIM[1] - LON_LIM[0])
lat_range = abs(LAT_LIM[1] - LAT_LIM[0])
panel_h   = 7.0
panel_w   = panel_h * (lon_range / lat_range)

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

    ax.set_xlim(LON_LIM); ax.set_ylim(LAT_LIM)
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

fig.suptitle('SKHASH Single-Event FMs  |  CC-SVD vs DL-hq (conf ≥ 0.8)  |  '
             'Q=A+B  |  All depths  |  1D axial velocity',
             fontsize=11, y=1.01)
plt.tight_layout()

out = OUT_DIR / 'fm_skhash_cc_vs_dlhq.png'
fig.savefig(out, dpi=150, bbox_inches='tight')
plt.close(fig)
print(f'\nSaved → {out}')
