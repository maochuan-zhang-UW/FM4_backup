"""
plot_depth_slices_compare.py
============================
Depth slice maps every 0.5 km, comparing 3D multi-model vs 1D axial_1D.
A+B+C quality events only.

Output: 02-data/G_FM/depth_slice_compare_X.X-X.Xkm.png  (one per slice)
        02-data/G_FM/depth_slice_compare_all.png          (all slices in one figure)

Run with FM_ML env:
  /opt/miniconda3/envs/FM_ML/bin/python plot_depth_slices_compare.py
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
from matplotlib.patches import Rectangle
from obspy.imaging.beachball import beach

from fault_classification import classify_fault_ptb

# ── Paths ──────────────────────────────────────────────────────────────────────
SKHASH_DIR = Path('/Users/mczhang/Documents/GitHub/SKHASH/SKHASH7')
OUT_DIR    = Path('/Users/mczhang/Documents/GitHub/FM7/02-data/G_FM')
OUT_DIR.mkdir(exist_ok=True)

# ── Load data ──────────────────────────────────────────────────────────────────
def classify_fault(rake):
    r = rake % 360
    if r > 180: r -= 360
    if -120 <= r <= -60:  return 'N'
    if   60 <= r <=  120: return 'R'
    if abs(r) <= 30 or abs(r) >= 150: return 'S'
    return 'U'

def load_run(y1_path, y2_path):
    dfs = []
    for path, yr in [(y1_path, 'Y1'), (y2_path, 'Y2')]:
        if path.exists():
            df = pd.read_csv(path)
            df['year'] = yr
            dfs.append(df)
    df = pd.concat(dfs, ignore_index=True)
    df['fault_type'] = [
        classify_fault_ptb(s, d, r)
        for s, d, r in zip(df['strike'], df['dip'], df['rake'])
    ]
    return df

df3d = load_run(
    SKHASH_DIR / 'examples/hash3/OUT_Y1_ge8/out.txt',
    SKHASH_DIR / 'examples/hash3/OUT_Y2_ge8/out.txt')
df1d = load_run(
    SKHASH_DIR / 'examples/hash3/OUT_Y1_ge8_1D/out.txt',
    SKHASH_DIR / 'examples/hash3/OUT_Y2_ge8_1D/out.txt')

print(f'3D: {len(df3d):,} FMs   1D: {len(df1d):,} FMs')

# ── Constants ──────────────────────────────────────────────────────────────────
LON_LIM    = (-130.04, -129.97)
LAT_LIM    = (45.908,  46.001)
MEAN_LAT   = np.mean(LAT_LIM)
BB_WIDTH   = 0.0012
QUAL_FILTER = ['A', 'B', 'C']
QUAL_ORDER  = {'A': 4, 'B': 3, 'C': 2, 'D': 1}

FAULT_COLOR = {
    'N': [0.15, 0.25, 0.85],
    'R': [0.85, 0.15, 0.15],
    'S': [0.10, 0.70, 0.20],
    'U': [0.40, 0.40, 0.40],
}

CALDERA_RIM = np.array([
    [-130.004785563058, 45.9207755734405],
    [-130.010476202888, 45.9238241104543],
    [-130.018881564079, 45.9351908809594],
    [-130.023946125193, 45.9412238501725],
    [-130.028718653506, 45.9498812001140],
    [-130.030451219380, 45.9511765797916],
    [-130.030679485650, 45.9542732243167],
    [-130.031733279709, 45.9558130656063],
    [-130.031444653500, 45.9586760104296],
    [-130.036188782208, 45.9656647517656],
    [-130.036950110789, 45.9698291665232],
    [-130.039953347000, 45.9750458167927],
    [-130.038595675479, 45.9847117727418],
    [-130.035927416999, 45.9883113986506],
    [-130.018067675296, 45.9933582886740],
    [-130.013629193751, 45.9937552841350],
    [-130.010365710979, 45.9929499241491],
    [-130.008647442296, 45.9924883829037],
    [-130.007262470669, 45.9915471582374],
    [-130.006042022411, 45.9902469280907],
    [-130.005178629490, 45.9897778053610],
    [-130.001868199523, 45.9863506519894],
    [-130.001154359192, 45.9846883853932],
    [-130.000949059432, 45.9827833001814],
    [-129.999393534330, 45.9818434725493],
    [-129.997797388662, 45.9786395525337],
    [-129.995357566829, 45.9760388622191],
    [-129.993956176267, 45.9741441737512],
    [-129.993678708114, 45.9681875631427],
    [-129.993140494256, 45.9667620754035],
    [-129.992087550788, 45.9652218741086],
    [-129.991186410747, 45.9626077204113],
    [-129.989604036931, 45.9601186407320],
    [-129.989238986151, 45.9588108137369],
    [-129.989728217453, 45.9574955894078],
    [-129.985484098670, 45.9494279735802],
    [-129.984788122490, 45.9487188881587],
])

STA_OO = {
    'AS1': (45.93356, -129.99920), 'AS2': (45.93377, -130.01410),
    'CC1': (45.95468, -130.00890), 'EC1': (45.94958, -129.97970),
    'EC2': (45.93967, -129.97380), 'EC3': (45.93607, -129.97850),
    'ID1': (45.92573, -129.97800),
}

REGIONS = {
    'West': dict(lat=(45.930, 45.953), lon=(-130.029, -130.008)),
    'East': dict(lat=(45.930, 45.971), lon=(-130.0015, -129.975)),
    'ID':   dict(lat=(45.921, 45.929), lon=(-130.004, -129.975)),
}
REG_COLORS = {'West': 'navy', 'East': 'darkgreen', 'ID': 'purple'}

# ── Single-panel plotting helper ───────────────────────────────────────────────
def draw_slice(ax, sub, d0, d1, label):
    """Draw beach balls for events in depth slice [d0, d1) onto ax."""
    ax.set_facecolor([0.85, 0.92, 0.97])
    ax.set_aspect(1.0 / np.cos(np.radians(MEAN_LAT)))

    sub = sub.sort_values('qrank')
    for _, row in sub.iterrows():
        color = FAULT_COLOR[row['fault_type']]
        try:
            bb = beach([row['strike'], row['dip'], row['rake']],
                       xy=(row['origin_lon'], row['origin_lat']),
                       width=BB_WIDTH, linewidth=0.2,
                       facecolor=color, alpha=0.85)
            bb.set_transform(ax.transData)
            bb.set_zorder(2)
            ax.add_collection(bb)
        except Exception:
            pass

    # Caldera
    ax.plot(CALDERA_RIM[:, 0], CALDERA_RIM[:, 1], '-',
            color='black', linewidth=1.2, zorder=7)

    # OO stations
    ax.plot([v[1] for v in STA_OO.values()],
            [v[0] for v in STA_OO.values()],
            '^', color='yellow', markeredgecolor='k',
            markeredgewidth=0.4, markersize=5, zorder=8)

    # Regions
    for name, reg in REGIONS.items():
        x0, x1 = reg['lon']; y0, y1 = reg['lat']
        ax.add_patch(Rectangle((x0, y0), x1-x0, y1-y0,
                                linewidth=1.2, edgecolor=REG_COLORS[name],
                                facecolor='none', linestyle='--', zorder=6))

    ax.set_xlim(LON_LIM); ax.set_ylim(LAT_LIM)
    ticks = np.arange(-130.04, -129.96, 0.02)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f'{t:.2f}' for t in ticks], fontsize=6, rotation=45)
    ax.tick_params(axis='y', labelsize=6)
    ax.grid(True, linewidth=0.3, color='gray', alpha=0.5)
    ax.set_title(f'{label}\nn={len(sub)}', fontsize=8, pad=3)


# ── Depth bins ─────────────────────────────────────────────────────────────────
depth_bins = [(d, d + 0.5) for d in np.arange(0.0, 5.0, 0.5)]  # 0–0.5 … 4.5–5.0

# ── Shared filter ──────────────────────────────────────────────────────────────
def prep(df):
    sub = df[df['quality'].isin(QUAL_FILTER)].copy()
    sub = sub[(sub['origin_lat'].between(*LAT_LIM)) &
              (sub['origin_lon'].between(*LON_LIM))]
    sub['qrank'] = sub['quality'].map(QUAL_ORDER)
    return sub

sub3d_all = prep(df3d)
sub1d_all = prep(df1d)

# ── Individual slice figures ───────────────────────────────────────────────────
lon_range = LON_LIM[1] - LON_LIM[0]
lat_range = LAT_LIM[1] - LAT_LIM[0]
aspect = (lon_range * np.cos(np.radians(MEAN_LAT))) / lat_range
panel_w = 5.0

legend_handles = [
    mpatches.Patch(color=FAULT_COLOR['N'], label='Normal'),
    mpatches.Patch(color=FAULT_COLOR['R'], label='Reverse'),
    mpatches.Patch(color=FAULT_COLOR['S'], label='Strike-slip'),
    mpatches.Patch(color=FAULT_COLOR['U'], label='Oblique'),
]

print('\nGenerating individual depth slice figures...')
for d0, d1 in depth_bins:
    s3 = sub3d_all[(sub3d_all['origin_depth_km'] >= d0) &
                   (sub3d_all['origin_depth_km'] <  d1)]
    s1 = sub1d_all[(sub1d_all['origin_depth_km'] >= d0) &
                   (sub1d_all['origin_depth_km'] <  d1)]

    fig, axes = plt.subplots(1, 2,
                              figsize=(panel_w*2 + 0.3, panel_w / aspect),
                              facecolor='white')
    draw_slice(axes[0], s3, d0, d1, f'3D multi-model')
    draw_slice(axes[1], s1, d0, d1, f'1D axial_1D')

    for ax in axes:
        ax.set_xlabel('Longitude', fontsize=7)
    axes[0].set_ylabel('Latitude', fontsize=7)
    axes[1].legend(handles=legend_handles, loc='lower right',
                   fontsize=6, framealpha=0.9)

    fig.suptitle(f'Depth {d0:.1f} – {d1:.1f} km  |  Q=A+B+C',
                 fontsize=10, y=1.01)
    plt.tight_layout()
    tag = f'{d0:.1f}-{d1:.1f}km'.replace('.', 'p')
    out = OUT_DIR / f'depth_slice_compare_{d0:.1f}-{d1:.1f}km.png'
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  {d0:.1f}–{d1:.1f} km  → 3D:{len(s3)}  1D:{len(s1)}  → {out.name}')

# ── Combined overview figure (all slices in one page) ─────────────────────────
print('\nGenerating combined overview figure...')
n_slices = len(depth_bins)
ncols = 2   # 3D | 1D per slice
nrows = n_slices

fig_h = n_slices * (panel_w / aspect) + 0.3 * n_slices
fig, axes = plt.subplots(nrows, ncols,
                          figsize=(panel_w*2, fig_h),
                          facecolor='white')

for row_i, (d0, d1) in enumerate(depth_bins):
    s3 = sub3d_all[(sub3d_all['origin_depth_km'] >= d0) &
                   (sub3d_all['origin_depth_km'] <  d1)]
    s1 = sub1d_all[(sub1d_all['origin_depth_km'] >= d0) &
                   (sub1d_all['origin_depth_km'] <  d1)]
    draw_slice(axes[row_i, 0], s3, d0, d1,
               f'3D  {d0:.1f}–{d1:.1f} km  (n={len(s3)})')
    draw_slice(axes[row_i, 1], s1, d0, d1,
               f'1D  {d0:.1f}–{d1:.1f} km  (n={len(s1)})')
    axes[row_i, 0].set_ylabel(f'{d0:.1f}–{d1:.1f} km', fontsize=7)

# Column headers
axes[0, 0].set_title('3D multi-model\n' + axes[0, 0].get_title(), fontsize=8, pad=3)
axes[0, 1].set_title('1D axial_1D\n'    + axes[0, 1].get_title(), fontsize=8, pad=3)
axes[-1, 1].legend(handles=legend_handles, loc='lower right',
                   fontsize=6, framealpha=0.9)

fig.suptitle('Depth slice comparison  3D vs 1D  |  Q=A+B+C  |  every 0.5 km',
             fontsize=11, y=1.002)
plt.tight_layout(h_pad=0.8)
out = OUT_DIR / 'depth_slice_compare_all.png'
fig.savefig(out, dpi=120, bbox_inches='tight')
plt.close(fig)
print(f'  → {out}')
print('\nDone.')
