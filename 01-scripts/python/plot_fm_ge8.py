"""
plot_fm_ge8.py
==============
Plot focal mechanism beach balls (map + depth cross-section) for the
≥8-station CC-augmented pipeline results.  Reads SKHASH out.txt for
Y1 and/or Y2, combines them, and produces:

  02-data/G_FM/ge8_FM_map_ABC.png
  02-data/G_FM/ge8_FM_map_AB.png
  02-data/G_FM/ge8_FM_depth_section_ABC.png
  02-data/G_FM/ge8_FM_depth_section_AB.png
  02-data/G_FM/ge8_FM_map_ABC_Y1.png   (Y1 only)
  02-data/G_FM/ge8_FM_map_ABC_Y2.png   (Y2 only)

Run with FM_ML env (obspy required for beach balls):
  /opt/miniconda3/envs/FM_ML/bin/python plot_fm_ge8.py
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

# ── Paths ─────────────────────────────────────────────────────────────────────
SKHASH_DIR = Path('/Users/mczhang/Documents/GitHub/SKHASH/SKHASH7')
OUT_Y1     = SKHASH_DIR / 'examples/hash3/OUT_Y1_ge8/out.txt'
OUT_Y2     = SKHASH_DIR / 'examples/hash3/OUT_Y2_ge8/out.txt'
OUT_DIR    = Path('/Users/mczhang/Documents/GitHub/FM7/02-data/G_FM')
OUT_DIR.mkdir(exist_ok=True)

# ── Load results ──────────────────────────────────────────────────────────────
dfs = []
for path, year in [(OUT_Y1, 'Y1'), (OUT_Y2, 'Y2')]:
    if path.exists():
        df = pd.read_csv(path)
        df['year'] = year
        dfs.append(df)
        print(f'{year}: {len(df)} FMs  quality={df.quality.value_counts().to_dict()}')
    else:
        print(f'{year}: {path} not found — skipping')

if not dfs:
    print('No SKHASH output files found. Run SKHASH first.')
    raise SystemExit(1)

df_all = pd.concat(dfs, ignore_index=True)
print(f'\nTotal: {len(df_all)} FMs')
print(df_all['quality'].value_counts().sort_index())

# ── Fault-type classification ──────────────────────────────────────────────────
def classify_fault(rake):
    r = rake % 360
    if r > 180: r -= 360
    if -120 <= r <= -60:  return 'N'
    if   60 <= r <=  120: return 'R'
    if abs(r) <= 30 or abs(r) >= 150: return 'S'
    return 'U'

FAULT_COLOR = {
    'N': [0.15, 0.25, 0.85],   # blue  — normal
    'R': [0.85, 0.15, 0.15],   # red   — reverse
    'S': [0.10, 0.70, 0.20],   # green — strike-slip
    'U': [0.40, 0.40, 0.40],   # grey  — oblique
}

# ── Caldera rim (Chadwick) ─────────────────────────────────────────────────────
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
df_all['fault_type'] = [
    classify_fault_ptb(s, d, r)
    for s, d, r in zip(df_all['strike'], df_all['dip'], df_all['rake'])
]
print('\nFault type (all):')
print(df_all['fault_type'].value_counts())

# ── Map parameters ─────────────────────────────────────────────────────────────
LON_LIM = (-130.04, -129.97)
LAT_LIM = (45.908,  46.001)

REGIONS = {
    'West': dict(lat=(45.930, 45.953), lon=(-130.029, -130.008)),
    'East': dict(lat=(45.930, 45.971), lon=(-130.0015, -129.975)),
    'ID':   dict(lat=(45.921, 45.929), lon=(-130.004, -129.975)),
}

# OO stations (shared across years)
STA_OO = {
    'AS1': (45.93356, -129.99920), 'AS2': (45.93377, -130.01410),
    'CC1': (45.95468, -130.00890), 'EC1': (45.94958, -129.97970),
    'EC2': (45.93967, -129.97380), 'EC3': (45.93607, -129.97850),
    'ID1': (45.92573, -129.97800),
}
# Y1 2F stations (A-suffix)
STA_Y1 = {
    '01A': (46.01933, -130.005375), '02A': (46.000625, -130.022257),
    '03A': (45.992922, -129.992782), '04A': (45.988743, -130.048395),
    '05A': (45.974390, -129.977972), '06A': (45.982938, -130.014128),
    '07A': (45.969067, -130.029852), '08A': (45.964147, -130.004503),
    '09A': (45.970840, -130.062393), '10A': (45.960338, -129.949987),
    '11A': (45.949590, -130.033977), '12A': (45.915367, -130.021773),
    '13A': (45.907817, -129.971555), '14A': (45.899965, -130.007573),
}
# Y2 2F stations (B-suffix)
STA_Y2 = {
    '01B': (46.01935, -130.00526), '02B': (46.00111, -130.02930),
    '03B': (45.99250, -129.99397), '04B': (45.98834, -130.04896),
    '05B': (45.97327, -129.97728), '06B': (45.98178, -130.01423),
    '07B': (45.96950, -130.02963), '08B': (45.96336, -130.00441),
    '09B': (45.97194, -130.05972), '10B': (45.95867, -129.94982),
    '11B': (45.94851, -130.03801), '12B': (45.91473, -130.02173),
    '13B': (45.90472, -129.97240), '14B': (45.89913, -130.00612),
}
STA_COORDS = {**STA_OO, **STA_Y1, **STA_Y2}

BB_WIDTH = 0.0012   # beach ball display width (degrees)
QUAL_ORDER = {'A': 4, 'B': 3, 'C': 2, 'D': 1}

# ── Map plot function ─────────────────────────────────────────────────────────
def plot_fm_map(events_df, quality_filter, title_suffix, out_png, show_stations=True):
    sub = events_df[events_df['quality'].isin(quality_filter)].copy()
    sub = sub[(sub['origin_lat'].between(*LAT_LIM)) &
              (sub['origin_lon'].between(*LON_LIM))]
    n = len(sub)
    print(f'{out_png.name}: {n} events')

    mean_lat = np.mean(LAT_LIM)
    lon_range = LON_LIM[1] - LON_LIM[0]
    lat_range = LAT_LIM[1] - LAT_LIM[0]
    # Figure width/height so 1° lon = cos(lat) * 1° lat in physical space
    aspect = (lon_range * np.cos(np.radians(mean_lat))) / lat_range
    fig, ax = plt.subplots(figsize=(8.0, 8.0 / aspect), facecolor='white')
    ax.set_facecolor([0.85, 0.92, 0.97])
    # Make beach balls round: tell matplotlib that 1 lon-unit = 1/cos(lat) lat-units
    ax.set_aspect(1.0 / np.cos(np.radians(mean_lat)))

    sub['qrank'] = sub['quality'].map(QUAL_ORDER)
    sub = sub.sort_values('qrank')

    for _, row in sub.iterrows():
        color = FAULT_COLOR[row['fault_type']]
        try:
            bb = beach([row['strike'], row['dip'], row['rake']],
                       xy=(row['origin_lon'], row['origin_lat']),
                       width=BB_WIDTH, linewidth=0.3,
                       facecolor=color, alpha=0.85)
            bb.set_transform(ax.transData)
            bb.set_zorder(2)
            ax.add_collection(bb)
        except Exception:
            pass

    # Caldera rim
    ax.plot(CALDERA_RIM[:, 0], CALDERA_RIM[:, 1], '-',
            color='black', linewidth=1.5, zorder=7, label='Caldera rim')

    if show_stations:
        # OO stations
        oo_lats = [v[0] for k, v in STA_OO.items()]
        oo_lons = [v[1] for k, v in STA_OO.items()]
        ax.plot(oo_lons, oo_lats, '^', color='yellow', markeredgecolor='k',
                markeredgewidth=0.5, markersize=6, zorder=8, label='OO sta')
        # Y1 2F stations
        y1_lats = [v[0] for v in STA_Y1.values()]
        y1_lons = [v[1] for v in STA_Y1.values()]
        ax.plot(y1_lons, y1_lats, 's', color='cyan', markeredgecolor='k',
                markeredgewidth=0.5, markersize=5, zorder=8, label='Y1 2F sta')
        # Y2 2F stations
        y2_lats = [v[0] for v in STA_Y2.values()]
        y2_lons = [v[1] for v in STA_Y2.values()]
        ax.plot(y2_lons, y2_lats, 's', color='magenta', markeredgecolor='k',
                markeredgewidth=0.5, markersize=5, zorder=8, label='Y2 2F sta')

    reg_colors = {'West': 'navy', 'East': 'darkgreen', 'ID': 'purple'}
    for name, reg in REGIONS.items():
        x0, x1 = reg['lon']; y0, y1 = reg['lat']
        ax.add_patch(Rectangle((x0, y0), x1-x0, y1-y0,
                                linewidth=1.5, edgecolor=reg_colors[name],
                                facecolor='none', linestyle='--', zorder=6))
        ax.text(x0, y1+0.001, name, fontsize=8, color=reg_colors[name], va='bottom')

    legend_handles = [
        mpatches.Patch(color=FAULT_COLOR['N'], label='Normal'),
        mpatches.Patch(color=FAULT_COLOR['R'], label='Reverse'),
        mpatches.Patch(color=FAULT_COLOR['S'], label='Strike-slip'),
        mpatches.Patch(color=FAULT_COLOR['U'], label='Oblique'),
    ]
    ax.legend(handles=legend_handles, loc='lower right', fontsize=8,
              framealpha=0.9, edgecolor='gray')

    ax.set_xlim(LON_LIM); ax.set_ylim(LAT_LIM)
    ax.set_xlabel('Longitude', fontsize=12); ax.set_ylabel('Latitude', fontsize=12)
    q_str = '+'.join(quality_filter)
    ax.set_title(f'{n} FMs  (Y1+Y2 ≥8-sta, Q={q_str})  {title_suffix}', fontsize=11)
    ticks = np.arange(-130.04, -129.96, 0.01)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f'{t:.2f}' for t in ticks], fontsize=8)
    ax.tick_params(axis='y', labelsize=8)
    ax.grid(True, linewidth=0.4, color='gray', alpha=0.5)
    plt.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  → {out_png}')

# ── Depth cross-section function ──────────────────────────────────────────────
def plot_depth_section(events_df, quality_filter, out_png, title_tag='Y1+Y2'):
    sub = events_df[events_df['quality'].isin(quality_filter)].copy()
    sub = sub[(sub['origin_lat'].between(*LAT_LIM)) &
              (sub['origin_lon'].between(*LON_LIM))]
    n = len(sub)

    KM_PER_DEG = 111.32 * np.cos(np.radians(np.mean(LAT_LIM)))
    sub['x_km'] = (sub['origin_lon'] - LON_LIM[0]) * KM_PER_DEG
    x_range = (LON_LIM[1] - LON_LIM[0]) * KM_PER_DEG
    depth_max = 5.0
    BB_KM = 0.12

    fig, ax = plt.subplots(figsize=(9, 5), facecolor='white')
    ax.set_facecolor([0.85, 0.92, 0.97])

    sub['qrank'] = sub['quality'].map(QUAL_ORDER)
    sub = sub.sort_values('qrank')

    for _, row in sub.iterrows():
        color = FAULT_COLOR[row['fault_type']]
        try:
            bb = beach([row['strike'], row['dip'], row['rake']],
                       xy=(row['x_km'], -row['origin_depth_km']),
                       width=BB_KM, linewidth=0.3,
                       facecolor=color, alpha=0.85)
            bb.set_transform(ax.transData)
            bb.set_zorder(2)
            ax.add_collection(bb)
        except Exception:
            pass

    xtick_lons = np.arange(-130.04, -129.96, 0.01)
    xtick_km   = [(l - LON_LIM[0]) * KM_PER_DEG for l in xtick_lons]
    ax.set_xticks(xtick_km)
    ax.set_xticklabels([f'{l:.2f}' for l in xtick_lons], fontsize=8)
    ax.set_xlim(0, x_range); ax.set_ylim(-depth_max, 0)
    ax.set_yticks(-np.arange(0, depth_max+1, 1))
    ax.set_yticklabels([f'{d:.0f}' for d in np.arange(0, depth_max+1, 1)])
    ax.set_xlabel('Longitude', fontsize=12); ax.set_ylabel('Depth (km)', fontsize=12)
    q_str = '+'.join(quality_filter)
    ax.set_title(f'{n} FMs — E–W depth section ({title_tag} ≥8-sta, Q={q_str})', fontsize=11)
    ax.grid(True, linewidth=0.4, color='gray', alpha=0.5)
    legend_handles = [
        mpatches.Patch(color=FAULT_COLOR['N'], label='Normal'),
        mpatches.Patch(color=FAULT_COLOR['R'], label='Reverse'),
        mpatches.Patch(color=FAULT_COLOR['S'], label='Strike-slip'),
        mpatches.Patch(color=FAULT_COLOR['U'], label='Oblique'),
    ]
    ax.legend(handles=legend_handles, loc='lower right', fontsize=8,
              framealpha=0.9, edgecolor='gray')
    plt.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  → {out_png}')

# ── Generate all plots ─────────────────────────────────────────────────────────
print('\n=== Generating maps ===')
plot_fm_map(df_all, ['A','B','C'], 'N=normal R=reverse S=SS U=oblique',
            OUT_DIR / 'ge8_FM_map_ABC.png')
plot_fm_map(df_all, ['A','B'], 'A+B quality only',
            OUT_DIR / 'ge8_FM_map_AB.png')

print('\n=== Generating depth sections ===')
plot_depth_section(df_all, ['A','B','C'], OUT_DIR / 'ge8_FM_depth_section_ABC.png')
plot_depth_section(df_all, ['A','B'],     OUT_DIR / 'ge8_FM_depth_section_AB.png')

# Per-year maps if both available
for year in ['Y1', 'Y2']:
    df_yr = df_all[df_all['year'] == year]
    if len(df_yr) > 0:
        plot_fm_map(df_yr, ['A','B','C'], year,
                    OUT_DIR / f'ge8_FM_map_ABC_{year}.png')
        plot_depth_section(df_yr, ['A','B','C'],
                           OUT_DIR / f'ge8_FM_depth_section_ABC_{year}.png',
                           title_tag=year)

# ── Summary stats ──────────────────────────────────────────────────────────────
print('\n=== Summary ===')
print(f'Total FMs: {len(df_all)}')
for year in ['Y1', 'Y2', 'combined']:
    df_s = df_all if year == 'combined' else df_all[df_all['year'] == year]
    if len(df_s) == 0:
        continue
    ab = df_s[df_s.quality.isin(['A','B'])]
    print(f'\n{year} ({len(df_s)} total, {len(ab)} A+B):')
    print(f'  Quality: {df_s.quality.value_counts().sort_index().to_dict()}')
    print(f'  Fault:   {df_s.fault_type.value_counts().to_dict()}')
    if len(ab) > 0:
        print(f'  A+B mean strike: {ab.strike.mean():.1f}±{ab.strike.std():.1f}')
        print(f'  A+B mean dip:    {ab.dip.mean():.1f}±{ab.dip.std():.1f}')
        print(f'  A+B mean rake:   {ab.rake.mean():.1f}±{ab.rake.std():.1f}')

print('\nDone.')
