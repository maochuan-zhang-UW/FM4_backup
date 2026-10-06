"""
plot_fm_compare_1D.py
=====================
Compare focal mechanisms from the 3D multi-model run (ge8) vs the new
single axial_1D model run (ge8_1D).

Produces:
  02-data/G_FM/compare_1D_map_ABC.png        — side-by-side maps
  02-data/G_FM/compare_1D_map_AB.png
  02-data/G_FM/compare_1D_kagan_hist.png     — Kagan angle distribution
  02-data/G_FM/compare_1D_quality_bar.png    — quality grade comparison
  02-data/G_FM/compare_1D_faulttype_bar.png  — fault-type composition

Run with FM_ML env:
  /opt/miniconda3/envs/FM_ML/bin/python plot_fm_compare_1D.py
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
OUT_DIR    = Path('/Users/mczhang/Documents/GitHub/FM7/02-data/G_FM')
OUT_DIR.mkdir(exist_ok=True)

RUNS = {
    '3D':  {'Y1': SKHASH_DIR / 'examples/hash3/OUT_Y1_ge8/out.txt',
             'Y2': SKHASH_DIR / 'examples/hash3/OUT_Y2_ge8/out.txt'},
    '1D':  {'Y1': SKHASH_DIR / 'examples/hash3/OUT_Y1_ge8_1D/out.txt',
             'Y2': SKHASH_DIR / 'examples/hash3/OUT_Y2_ge8_1D/out.txt'},
}

# ── Map / display constants ────────────────────────────────────────────────────
LON_LIM = (-130.04, -129.97)
LAT_LIM = (45.908,  46.001)
BB_WIDTH = 0.0012
QUAL_ORDER = {'A': 4, 'B': 3, 'C': 2, 'D': 1}

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
STA_Y1 = {
    '01A': (46.01933, -130.005375), '02A': (46.000625, -130.022257),
    '03A': (45.992922, -129.992782), '04A': (45.988743, -130.048395),
    '05A': (45.974390, -129.977972), '06A': (45.982938, -130.014128),
    '07A': (45.969067, -130.029852), '08A': (45.964147, -130.004503),
    '09A': (45.970840, -130.062393), '10A': (45.960338, -129.949987),
    '11A': (45.949590, -130.033977), '12A': (45.915367, -130.021773),
    '13A': (45.907817, -129.971555), '14A': (45.899965, -130.007573),
}
STA_Y2 = {
    '01B': (46.01935, -130.00526), '02B': (46.00111, -130.02930),
    '03B': (45.99250, -129.99397), '04B': (45.98834, -130.04896),
    '05B': (45.97327, -129.97728), '06B': (45.98178, -130.01423),
    '07B': (45.96950, -130.02963), '08B': (45.96336, -130.00441),
    '09B': (45.97194, -130.05972), '10B': (45.95867, -129.94982),
    '11B': (45.94851, -130.03801), '12B': (45.91473, -130.02173),
    '13B': (45.90472, -129.97240), '14B': (45.89913, -130.00612),
}

REGIONS = {
    'West': dict(lat=(45.930, 45.953), lon=(-130.029, -130.008)),
    'East': dict(lat=(45.930, 45.971), lon=(-130.0015, -129.975)),
    'ID':   dict(lat=(45.921, 45.929), lon=(-130.004, -129.975)),
}
REG_COLORS = {'West': 'navy', 'East': 'darkgreen', 'ID': 'purple'}

# ── Helpers ────────────────────────────────────────────────────────────────────
def classify_fault(rake):
    r = rake % 360
    if r > 180: r -= 360
    if -120 <= r <= -60:  return 'N'
    if   60 <= r <=  120: return 'R'
    if abs(r) <= 30 or abs(r) >= 150: return 'S'
    return 'U'


def load_run(run_paths):
    dfs = []
    for year, path in run_paths.items():
        if path.exists():
            df = pd.read_csv(path)
            df['year'] = year
            dfs.append(df)
        else:
            print(f'  WARNING: {path} not found')
    if not dfs:
        return None
    df = pd.concat(dfs, ignore_index=True)
    df['fault_type'] = [
        classify_fault_ptb(s, d, r)
        for s, d, r in zip(df['strike'], df['dip'], df['rake'])
    ]
    return df


def kagan_angle(s1, d1, r1, s2, d2, r2):
    """Minimum rotation angle (degrees) between two double-couple mechanisms."""
    def sdr_to_matrix(s, d, r):
        s, d, r = np.radians(s), np.radians(d), np.radians(r)
        T = np.array([
            [-np.sin(d)*np.sin(s),  np.sin(d)*np.cos(s),   np.cos(d)],
            [ np.cos(r)*np.cos(s)+np.sin(r)*np.cos(d)*np.sin(s),
              np.cos(r)*np.sin(s)-np.sin(r)*np.cos(d)*np.cos(s),
              np.sin(r)*np.sin(d)],
            [-np.sin(r)*np.cos(s)+np.cos(r)*np.cos(d)*np.sin(s),
             -np.sin(r)*np.sin(s)-np.cos(r)*np.cos(d)*np.cos(s),
              np.cos(r)*np.sin(d)],
        ])
        return T

    M1 = sdr_to_matrix(s1, d1, r1)
    M2 = sdr_to_matrix(s2, d2, r2)
    # 4 equivalent representations (two planes × two polarities)
    angles = []
    for flip_n in [1, -1]:
        for flip_s in [1, -1]:
            M1f = M1.copy()
            M1f[0] *= flip_n
            M1f[1] *= flip_s
            R = M1f @ M2.T
            trace = np.clip((np.trace(R) - 1) / 2, -1, 1)
            angles.append(np.degrees(np.arccos(trace)))
    return min(angles)


# ── Load data ─────────────────────────────────────────────────────────────────
print('Loading results...')
dfs = {}
for label, paths in RUNS.items():
    df = load_run(paths)
    if df is not None:
        print(f'  {label}: {len(df)} FMs  quality={df.quality.value_counts().sort_index().to_dict()}')
        dfs[label] = df
    else:
        print(f'  {label}: no data found — run SKHASH first')

if len(dfs) < 2:
    print('Need both runs to compare. Exiting.')
    raise SystemExit(1)

df3d = dfs['3D']
df1d = dfs['1D']

# ── Match events by event_id for comparison ───────────────────────────────────
merged = pd.merge(
    df3d[['event_id','strike','dip','rake','quality','fault_type',
          'origin_lat','origin_lon','origin_depth_km']],
    df1d[['event_id','strike','dip','rake','quality','fault_type']],
    on='event_id', suffixes=('_3D','_1D')
)
print(f'\nMatched events: {len(merged):,} (of {len(df3d)} 3D, {len(df1d)} 1D)')

# Compute Kagan angles for matched events
print('Computing Kagan angles...')
kagan = np.array([
    kagan_angle(row.strike_3D, row.dip_3D, row.rake_3D,
                row.strike_1D, row.dip_1D, row.rake_1D)
    for _, row in merged.iterrows()
])
merged['kagan'] = kagan
print(f'  Mean Kagan: {kagan.mean():.1f}°  Median: {np.median(kagan):.1f}°')
print(f'  <15°: {(kagan<15).mean()*100:.1f}%  <30°: {(kagan<30).mean()*100:.1f}%  <45°: {(kagan<45).mean()*100:.1f}%')

# ── Side-by-side map function ─────────────────────────────────────────────────
def make_axes_pair(quality_filter):
    mean_lat = np.mean(LAT_LIM)
    lon_range = LON_LIM[1] - LON_LIM[0]
    lat_range = LAT_LIM[1] - LAT_LIM[0]
    aspect = (lon_range * np.cos(np.radians(mean_lat))) / lat_range
    panel_w = 7.0
    fig, axes = plt.subplots(1, 2, figsize=(panel_w*2 + 0.5, panel_w / aspect),
                              facecolor='white')
    for ax, (label, df) in zip(axes, [('3D multi-model', df3d), ('1D axial_1D', df1d)]):
        sub = df[df['quality'].isin(quality_filter)].copy()
        sub = sub[(sub['origin_lat'].between(*LAT_LIM)) &
                  (sub['origin_lon'].between(*LON_LIM))]
        sub['qrank'] = sub['quality'].map(QUAL_ORDER)
        sub = sub.sort_values('qrank')

        ax.set_facecolor([0.85, 0.92, 0.97])
        ax.set_aspect(1.0 / np.cos(np.radians(mean_lat)))

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

        ax.plot(CALDERA_RIM[:, 0], CALDERA_RIM[:, 1], '-',
                color='black', linewidth=1.5, zorder=7)

        oo_lats = [v[0] for v in STA_OO.values()]
        oo_lons = [v[1] for v in STA_OO.values()]
        ax.plot(oo_lons, oo_lats, '^', color='yellow', markeredgecolor='k',
                markeredgewidth=0.5, markersize=6, zorder=8)
        y1_lats = [v[0] for v in STA_Y1.values()]
        y1_lons = [v[1] for v in STA_Y1.values()]
        ax.plot(y1_lons, y1_lats, 's', color='cyan', markeredgecolor='k',
                markeredgewidth=0.5, markersize=4, zorder=8)
        y2_lats = [v[0] for v in STA_Y2.values()]
        y2_lons = [v[1] for v in STA_Y2.values()]
        ax.plot(y2_lons, y2_lats, 's', color='magenta', markeredgecolor='k',
                markeredgewidth=0.5, markersize=4, zorder=8)

        for name, reg in REGIONS.items():
            x0, x1 = reg['lon']; y0, y1 = reg['lat']
            ax.add_patch(Rectangle((x0, y0), x1-x0, y1-y0,
                                    linewidth=1.5, edgecolor=REG_COLORS[name],
                                    facecolor='none', linestyle='--', zorder=6))
            ax.text(x0, y1+0.001, name, fontsize=7, color=REG_COLORS[name], va='bottom')

        q_str = '+'.join(quality_filter)
        ax.set_xlim(LON_LIM); ax.set_ylim(LAT_LIM)
        ticks = np.arange(-130.04, -129.96, 0.01)
        ax.set_xticks(ticks)
        ax.set_xticklabels([f'{t:.2f}' for t in ticks], fontsize=7, rotation=45)
        ax.tick_params(axis='y', labelsize=7)
        ax.set_xlabel('Longitude', fontsize=10)
        ax.set_ylabel('Latitude', fontsize=10)
        ax.grid(True, linewidth=0.4, color='gray', alpha=0.5)
        ax.set_title(f'{label}\n{len(sub)} FMs  Q={q_str}', fontsize=10)

    legend_handles = [
        mpatches.Patch(color=FAULT_COLOR['N'], label='Normal'),
        mpatches.Patch(color=FAULT_COLOR['R'], label='Reverse'),
        mpatches.Patch(color=FAULT_COLOR['S'], label='Strike-slip'),
        mpatches.Patch(color=FAULT_COLOR['U'], label='Oblique'),
    ]
    axes[1].legend(handles=legend_handles, loc='lower right', fontsize=8,
                   framealpha=0.9, edgecolor='gray')
    return fig


# ── Generate maps ─────────────────────────────────────────────────────────────
print('\n=== Side-by-side maps ===')
for qf, tag in [(['A','B','C'], 'ABC'), (['A','B'], 'AB')]:
    fig = make_axes_pair(qf)
    fig.suptitle(f'3D multi-model  vs  1D axial_1D  |  Q={tag}', fontsize=12, y=1.01)
    plt.tight_layout()
    out = OUT_DIR / f'compare_1D_map_{tag}.png'
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  → {out}')

# ── Kagan angle histogram ──────────────────────────────────────────────────────
print('\n=== Kagan angle histogram ===')
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), facecolor='white')

for ax, qf, tag in [(axes[0], ['A','B','C'], 'A+B+C'), (axes[1], ['A','B'], 'A+B only')]:
    sub = merged[merged['quality_3D'].isin(qf) & merged['quality_1D'].isin(qf)]
    kg = sub['kagan'].values
    bins = np.arange(0, 91, 5)
    ax.hist(kg, bins=bins, color='steelblue', edgecolor='white', linewidth=0.5)
    ax.axvline(np.median(kg), color='red', linestyle='--', linewidth=1.5,
               label=f'Median {np.median(kg):.1f}°')
    ax.axvline(np.mean(kg), color='orange', linestyle='--', linewidth=1.5,
               label=f'Mean {np.mean(kg):.1f}°')
    ax.set_xlabel('Kagan angle (°)', fontsize=11)
    ax.set_ylabel('Count', fontsize=11)
    ax.set_title(f'3D vs 1D FM difference  ({tag}, n={len(kg):,})', fontsize=11)
    ax.legend(fontsize=9)
    pcts = [(15, '<15°'), (30, '<30°'), (45, '<45°')]
    txt = '\n'.join([f'{p[1]}: {(kg<p[0]).mean()*100:.0f}%' for p in pcts])
    ax.text(0.97, 0.97, txt, transform=ax.transAxes, fontsize=9,
            va='top', ha='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    ax.set_xlim(0, 90)
    ax.grid(True, linewidth=0.4, alpha=0.5)

plt.tight_layout()
out = OUT_DIR / 'compare_1D_kagan_hist.png'
fig.savefig(out, dpi=150, bbox_inches='tight')
plt.close(fig)
print(f'  → {out}')

# ── Quality and fault-type comparison bars ─────────────────────────────────────
print('\n=== Quality & fault-type comparison ===')
fig, axes = plt.subplots(1, 2, figsize=(12, 5), facecolor='white')

# Quality
ax = axes[0]
quals = ['A', 'B', 'C', 'D']
x = np.arange(len(quals))
w = 0.35
cnt3d = [int((df3d['quality'] == q).sum()) for q in quals]
cnt1d = [int((df1d['quality'] == q).sum()) for q in quals]
ax.bar(x - w/2, cnt3d, w, label='3D multi-model', color='steelblue', edgecolor='k', linewidth=0.5)
ax.bar(x + w/2, cnt1d, w, label='1D axial_1D',    color='salmon',     edgecolor='k', linewidth=0.5)
for xi, v3, v1 in zip(x, cnt3d, cnt1d):
    ax.text(xi - w/2, v3 + 20, str(v3), ha='center', va='bottom', fontsize=8)
    ax.text(xi + w/2, v1 + 20, str(v1), ha='center', va='bottom', fontsize=8)
ax.set_xticks(x); ax.set_xticklabels(quals)
ax.set_xlabel('Quality grade', fontsize=11); ax.set_ylabel('Count', fontsize=11)
ax.set_title('Quality grade distribution', fontsize=11)
ax.legend(fontsize=9); ax.grid(axis='y', linewidth=0.4, alpha=0.5)

# Fault type
ax = axes[1]
ftypes = ['N', 'R', 'S', 'U']
ft_labels = ['Normal', 'Reverse', 'Strike-slip', 'Oblique']
x = np.arange(len(ftypes))
sub3 = df3d[df3d['quality'].isin(['A','B'])]
sub1 = df1d[df1d['quality'].isin(['A','B'])]
cnt3d = [int((sub3['fault_type'] == f).sum()) for f in ftypes]
cnt1d = [int((sub1['fault_type'] == f).sum()) for f in ftypes]
fc = [FAULT_COLOR[f] for f in ftypes]
bars3 = ax.bar(x - w/2, cnt3d, w, label='3D multi-model',
               color=fc, edgecolor='k', linewidth=0.5, alpha=0.9)
bars1 = ax.bar(x + w/2, cnt1d, w, label='1D axial_1D',
               color=fc, edgecolor='k', linewidth=0.5, alpha=0.5, hatch='//')
for xi, v3, v1 in zip(x, cnt3d, cnt1d):
    ax.text(xi - w/2, v3 + 5, str(v3), ha='center', va='bottom', fontsize=8)
    ax.text(xi + w/2, v1 + 5, str(v1), ha='center', va='bottom', fontsize=8)
ax.set_xticks(x); ax.set_xticklabels(ft_labels)
ax.set_xlabel('Fault type', fontsize=11); ax.set_ylabel('Count (A+B only)', fontsize=11)
ax.set_title('Fault type distribution (A+B quality)', fontsize=11)
import matplotlib.patches as mpatches
ax.legend(handles=[
    mpatches.Patch(facecolor='gray', edgecolor='k', label='3D multi-model'),
    mpatches.Patch(facecolor='gray', edgecolor='k', hatch='//', alpha=0.5, label='1D axial_1D'),
], fontsize=9)
ax.grid(axis='y', linewidth=0.4, alpha=0.5)

plt.tight_layout()
out = OUT_DIR / 'compare_1D_quality_faulttype.png'
fig.savefig(out, dpi=150, bbox_inches='tight')
plt.close(fig)
print(f'  → {out}')

# ── Velocity model comparison ─────────────────────────────────────────────────
print('\n=== Velocity model plot ===')
# Load the new 1D model
mod_data = []
for line in open('/Users/mczhang/Downloads/axial_1D.mod.260317').readlines()[1:]:
    p = line.split()
    if len(p) >= 2:
        d, vp = float(p[0]), float(p[1])
        if d >= 0:
            mod_data.append((d, vp))
mod1d = np.array(mod_data)

# Load a couple of 3D reference models for comparison
def load_velmod(path):
    rows = []
    for line in open(path):
        if line.startswith('#'): continue
        p = line.replace(',', ' ').split()
        if len(p) >= 2:
            d, vp = float(p[0]), float(p[1])
            if d <= 6:
                rows.append((d, vp))
    return np.array(rows)

vmod_dir = SKHASH_DIR / 'examples/velocity_models'
ref_models = {
    'W1': load_velmod(vmod_dir / 'velmod_W11.txt'),
    'W2': load_velmod(vmod_dir / 'velmod_W21.txt'),
    'E1': load_velmod(vmod_dir / 'velmod_E11.txt'),
    'E2': load_velmod(vmod_dir / 'velmod_E21.txt'),
}

fig, ax = plt.subplots(figsize=(5, 7), facecolor='white')
colors_ref = ['#aabbcc', '#99aacc', '#bbccaa', '#aaccaa']
for (name, vm), c in zip(ref_models.items(), colors_ref):
    ax.plot(vm[:, 1], -vm[:, 0], '-', color=c, linewidth=1.0, label=f'3D {name}', alpha=0.7)
ax.plot(mod1d[:, 1], -mod1d[:, 0], 'k-', linewidth=2.5, label='1D axial_1D', zorder=5)
ax.set_xlabel('Vp (km/s)', fontsize=12)
ax.set_ylabel('Depth (km)', fontsize=12)
ax.set_ylim(-5, 0)
ax.set_xlim(1.5, 7.5)
ax.set_yticks(-np.arange(0, 6, 1))
ax.set_yticklabels([f'{d:.0f}' for d in np.arange(0, 6, 1)])
ax.legend(fontsize=8, loc='lower right')
ax.set_title('Velocity models\n(new 1D vs subset of 3D models)', fontsize=10)
ax.grid(True, linewidth=0.4, alpha=0.5)
plt.tight_layout()
out = OUT_DIR / 'compare_1D_velmod.png'
fig.savefig(out, dpi=150, bbox_inches='tight')
plt.close(fig)
print(f'  → {out}')

# ── Summary ───────────────────────────────────────────────────────────────────
print('\n=== Summary ===')
for label, df in [('3D multi-model', df3d), ('1D axial_1D', df1d)]:
    ab = df[df['quality'].isin(['A','B'])]
    print(f'\n{label}  ({len(df):,} total, {len(ab):,} A+B):')
    print(f'  Quality: {df.quality.value_counts().sort_index().to_dict()}')
    print(f'  Fault:   {df.fault_type.value_counts().to_dict()}')
    if len(ab) > 0:
        print(f'  A+B mean strike: {ab.strike.mean():.1f}±{ab.strike.std():.1f}')
        print(f'  A+B mean dip:    {ab.dip.mean():.1f}±{ab.dip.std():.1f}')
        print(f'  A+B mean rake:   {ab.rake.mean():.1f}±{ab.rake.std():.1f}')

print('\nDone.')
