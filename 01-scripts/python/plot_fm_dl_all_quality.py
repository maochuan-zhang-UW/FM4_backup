"""
plot_fm_dl_all_quality.py
=========================
2-panel beach-ball map: Y1 and Y2 DL-hq FMs (SKHASH), all qualities.
A/B events: full colour.  C/D events: light/faded colour.

Output: 02-data/G_FM/fm_dl_all_quality.png

Run with FM_ML env:
  /opt/miniconda3/envs/FM_ML/bin/python plot_fm_dl_all_quality.py
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
LON_LIM  = (-130.04, -129.97)
LAT_LIM  = (45.908,  46.001)
BB_WIDTH = 0.0012

# A/B: full colour; C/D: pastel/faded
FAULT_COLOR = {
    'N': {'AB': [0.15, 0.25, 0.85], 'CD': [0.72, 0.78, 0.97]},
    'R': {'AB': [0.85, 0.15, 0.15], 'CD': [0.97, 0.72, 0.72]},
    'S': {'AB': [0.10, 0.70, 0.20], 'CD': [0.72, 0.93, 0.76]},
    'U': {'AB': [0.40, 0.40, 0.40], 'CD': [0.80, 0.80, 0.80]},
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

# ── Helpers ────────────────────────────────────────────────────────────────────
def classify_fault(rake):
    r = rake % 360
    if r > 180: r -= 360
    if -120 <= r <= -60:  return 'N'
    if   60 <= r <=  120: return 'R'
    if abs(r) <= 30 or abs(r) >= 150: return 'S'
    return 'U'

def load_cc(tag):
    path = HASH3_DIR / f'OUT_{tag}' / 'out.txt'
    df = pd.read_csv(path)
    df.columns = df.columns.str.strip()
    df['fault_type'] = [
        classify_fault_ptb(s, d, r)
        for s, d, r in zip(df['strike'], df['dip'], df['rake'])
    ]
    df['ab'] = df['quality'].isin(['A', 'B'])
    # clip to map bounds
    df = df[(df['origin_lat'].between(*LAT_LIM)) &
            (df['origin_lon'].between(*LON_LIM))].copy()
    # draw C/D first (background), then A/B on top
    df['zorder_key'] = df['ab'].map({False: 0, True: 1})
    return df.sort_values('zorder_key')

# ── Load ───────────────────────────────────────────────────────────────────────
data = {
    'Y1 DL-hq (SKHASH)': load_cc('Y1_dl_hq_sk'),
    'Y2 DL-hq (SKHASH)': load_cc('Y2_dl_hq_sk'),
}
for lbl, df in data.items():
    n_ab = df['ab'].sum(); n_cd = (~df['ab']).sum()
    print(f"{lbl}: {n_ab} A/B  +  {n_cd} C/D  =  {len(df)} total")

# ── Figure ─────────────────────────────────────────────────────────────────────
lon_range = abs(LON_LIM[1] - LON_LIM[0])
lat_range = abs(LAT_LIM[1] - LAT_LIM[0])
panel_h   = 8.0
panel_w   = panel_h * (lon_range / lat_range)

fig, axes = plt.subplots(1, 2,
                          figsize=(panel_w * 2 + 0.8, panel_h),
                          facecolor='white')

for ax, (label, df) in zip(axes, data.items()):
    ax.set_facecolor([0.85, 0.92, 0.97])
    ax.set_aspect('equal')

    for _, row in df.iterrows():
        qual_key = 'AB' if row['ab'] else 'CD'
        fc = FAULT_COLOR[row['fault_type']][qual_key]
        alpha = 0.85 if row['ab'] else 0.45
        lw    = 0.15 if row['ab'] else 0.10
        try:
            bb = beach([row['strike'], row['dip'], row['rake']],
                       xy=(row['origin_lon'], row['origin_lat']),
                       width=BB_WIDTH, linewidth=lw,
                       facecolor=fc, alpha=alpha)
            bb.set_transform(ax.transData)
            bb.set_zorder(2 if row['ab'] else 1)
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
    ax.set_xticklabels([f'{t:.2f}' for t in ticks], fontsize=7, rotation=45)
    ax.tick_params(axis='y', labelsize=7)
    ax.grid(True, linewidth=0.3, color='gray', alpha=0.5)
    ax.set_xlabel('Longitude', fontsize=9)

    n_ab = df['ab'].sum(); n_cd = (~df['ab']).sum()
    ax.set_title(f'{label}\nA/B: {n_ab}   C/D: {n_cd}   total: {len(df)}',
                 fontsize=10, pad=5)

axes[0].set_ylabel('Latitude', fontsize=9)

# Legend: fault type (A/B colours) + faded example
legend_handles = [
    mpatches.Patch(color=FAULT_COLOR['N']['AB'], label='Normal (A/B)'),
    mpatches.Patch(color=FAULT_COLOR['R']['AB'], label='Reverse (A/B)'),
    mpatches.Patch(color=FAULT_COLOR['S']['AB'], label='Strike-slip (A/B)'),
    mpatches.Patch(color=FAULT_COLOR['U']['AB'], label='Oblique (A/B)'),
    mpatches.Patch(color=FAULT_COLOR['N']['CD'], label='Normal (C/D)', alpha=0.6),
    mpatches.Patch(color=FAULT_COLOR['R']['CD'], label='Reverse (C/D)', alpha=0.6),
    mpatches.Patch(color=FAULT_COLOR['S']['CD'], label='Strike-slip (C/D)', alpha=0.6),
    mpatches.Patch(color=FAULT_COLOR['U']['CD'], label='Oblique (C/D)', alpha=0.6),
]
axes[-1].legend(handles=legend_handles, loc='lower right', fontsize=7,
                framealpha=0.9, ncol=2)

fig.suptitle('DL-hq Focal Mechanisms (SKHASH 1D)  |  All Qualities  |  '
             'Full colour = A/B,  Faded = C/D',
             fontsize=11, y=1.01)
plt.tight_layout()

out = OUT_DIR / 'fm_dl_all_quality.png'
fig.savefig(out, dpi=150, bbox_inches='tight')
plt.close(fig)
print(f'\nSaved → {out}')
