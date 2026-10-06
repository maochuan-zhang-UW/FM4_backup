"""
plot_kagan_skhash1d_vs_hash1d.py
=================================
Compute and plot Kagan angle differences between SKHASH 1D and HASH 1D
focal mechanisms for ge8 Y1+Y2 events.  All quality grades included,
then broken out by SKHASH quality.

Output: 02-data/G_FM/kagan_skhash1d_vs_hash1d.png

Run with FM_ML env:
  /opt/miniconda3/envs/FM_ML/bin/python plot_kagan_skhash1d_vs_hash1d.py
"""
import warnings
warnings.filterwarnings('ignore')

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── Paths ──────────────────────────────────────────────────────────────────────
SKHASH_DIR = Path('/Users/mczhang/Documents/GitHub/SKHASH/SKHASH7')
HASH_DIR   = Path('/Users/mczhang/Documents/GitHub/FM7/01-scripts/HASH')
DATA_DIR   = Path('/Users/mczhang/Documents/GitHub/FM7/02-data')
OUT_DIR    = DATA_DIR / 'G_FM'
OUT_DIR.mkdir(exist_ok=True)

# ── Kagan angle ────────────────────────────────────────────────────────────────
def _sdr2ptb(s, d, r):
    """Strike/dip/rake → P, T, B unit vectors (NED coords)."""
    s, d, r = np.radians(s), np.radians(d), np.radians(r)
    n  = np.array([-np.sin(d)*np.sin(s),
                    np.sin(d)*np.cos(s),
                   -np.cos(d)])
    sl = np.array([ np.cos(r)*np.cos(s) + np.sin(r)*np.cos(d)*np.sin(s),
                    np.cos(r)*np.sin(s) - np.sin(r)*np.cos(d)*np.cos(s),
                   -np.sin(r)*np.sin(d)])
    T = (n + sl) / np.sqrt(2)
    P = (n - sl) / np.sqrt(2)
    B = np.cross(T, P)
    return P, T, B

def kagan_angle(s1, d1, r1, s2, d2, r2):
    """Minimum rotation angle (degrees) between two focal mechanisms."""
    P1, T1, B1 = _sdr2ptb(s1, d1, r1)
    P2, T2, B2 = _sdr2ptb(s2, d2, r2)
    R2 = np.column_stack([P2, T2, B2])
    # 4 valid SO(3) sign combinations for R1
    min_ang = np.inf
    for sp, st, sb in [(1,1,1), (-1,-1,1), (1,-1,-1), (-1,1,-1)]:
        R1v = np.column_stack([sp*P1, st*T1, sb*B1])
        cos_a = np.clip((np.trace(R1v.T @ R2) - 1) / 2, -1, 1)
        min_ang = min(min_ang, np.degrees(np.arccos(cos_a)))
    return min_ang

kagan_vec = np.vectorize(kagan_angle)

# ── Load SKHASH 1D ─────────────────────────────────────────────────────────────
# SKHASH runs Y1 and Y2 separately with sequential event_ids restarting at 1.
# HASH combines them: Y1 cluster_idx 1-N_Y1, Y2 cluster_idx N_Y1+1 onward.
# Use eid_map to translate HASH cluster_idx → original catalog event_id,
# then reconstruct the same mapping for SKHASH via phase file ordering.
#
# Simpler approach: SKHASH event_id == sequential position in that year's input.
# Y1: SKHASH event_id k  ↔  HASH cluster_idx k
# Y2: SKHASH event_id k  ↔  HASH cluster_idx (k + N_Y1)

# Count Y1 events (cluster_idx 1..N_Y1) from eid_map
eid_map  = pd.read_csv(DATA_DIR / 'G_FM' / 'hash_ge8_1D_eid_map.csv')

# Determine N_Y1: Y1 phase file event count
n_y1 = sum(1 for line in open(
    SKHASH_DIR / 'examples/hash3/IN/north2_Y1_ge8.txt')
    if line[:4] in ('2022','2023','2024'))
print(f'Y1 events in SKHASH input: {n_y1:,}')

sk_y1 = pd.read_csv(SKHASH_DIR / 'examples/hash3/OUT_Y1_ge8_1D/out.txt')
sk_y1 = sk_y1[['event_id','strike','dip','rake','quality']].copy()
sk_y1['cluster_idx'] = sk_y1['event_id']          # Y1: direct match

sk_y2 = pd.read_csv(SKHASH_DIR / 'examples/hash3/OUT_Y2_ge8_1D/out.txt')
sk_y2 = sk_y2[['event_id','strike','dip','rake','quality']].copy()
sk_y2['cluster_idx'] = sk_y2['event_id'] + n_y1   # Y2: offset by N_Y1

sk1d = pd.concat([sk_y1, sk_y2], ignore_index=True)
print(f'SKHASH 1D: {len(sk1d):,} FMs  (Y1={len(sk_y1):,} Y2={len(sk_y2):,})')

# ── Load HASH 1D ───────────────────────────────────────────────────────────────
def load_hash(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            p = line.split()
            if len(p) < 31:
                continue
            try:
                rows.append(dict(
                    cluster_idx = int(p[0]),
                    strike      = float(p[21]),
                    dip         = float(p[22]),
                    rake        = float(p[23]),
                    quality     = p[28],
                ))
            except (ValueError, IndexError):
                continue
    return pd.DataFrame(rows)

hash1d = load_hash(HASH_DIR / 'hashout_ge8_1D_1.dat')
print(f'HASH 1D:   {len(hash1d):,} FMs')

# ── Merge on cluster_idx ───────────────────────────────────────────────────────
merged = sk1d[['cluster_idx','strike','dip','rake','quality']].merge(
    hash1d[['cluster_idx','strike','dip','rake','quality']].rename(
        columns={'strike':'h_strike','dip':'h_dip','rake':'h_rake','quality':'h_qual'}),
    on='cluster_idx', how='inner')
print(f'Matched:   {len(merged):,} events')

# ── Compute Kagan angles ───────────────────────────────────────────────────────
print('Computing Kagan angles...')
merged['kagan'] = kagan_vec(
    merged['strike'], merged['dip'],   merged['rake'],
    merged['h_strike'], merged['h_dip'], merged['h_rake'])

print(f'\nKagan angle statistics (all matched):')
print(merged['kagan'].describe().round(1))

# ── Plot ───────────────────────────────────────────────────────────────────────
QUAL_COLORS = {'A': '#d62728', 'B': '#ff7f0e', 'C': '#2ca02c', 'D': '#9467bd'}
QUAL_ORDER  = ['A', 'B', 'C', 'D']
bins = np.arange(0, 91, 5)

fig, axes = plt.subplots(1, 2, figsize=(13, 5), facecolor='white')

# Left: overall histogram coloured by SKHASH quality
ax = axes[0]
data_by_qual = [merged[merged['quality'] == q]['kagan'].values for q in QUAL_ORDER]
ax.hist(data_by_qual, bins=bins, stacked=True,
        color=[QUAL_COLORS[q] for q in QUAL_ORDER],
        label=[f'SKHASH {q} (n={len(d):,})' for q, d in zip(QUAL_ORDER, data_by_qual)],
        edgecolor='none', alpha=0.9)
ax.axvline(merged['kagan'].median(), color='k', ls='--', lw=1.2,
           label=f'Median = {merged["kagan"].median():.1f}°')
ax.set_xlabel('Kagan angle (°)', fontsize=11)
ax.set_ylabel('Count', fontsize=11)
ax.set_title(f'SKHASH 1D vs HASH 1D  |  all matched (n={len(merged):,})', fontsize=11)
ax.legend(fontsize=9, framealpha=0.9)
ax.set_xlim(0, 90)
ax.grid(True, alpha=0.3, linewidth=0.5)

# Right: CDF per quality
ax2 = axes[1]
for q in QUAL_ORDER:
    sub = merged[merged['quality'] == q]['kagan']
    if len(sub) == 0:
        continue
    sorted_k = np.sort(sub)
    cdf = np.arange(1, len(sorted_k)+1) / len(sorted_k)
    ax2.plot(sorted_k, cdf * 100, color=QUAL_COLORS[q], lw=1.5,
             label=f'{q} (n={len(sub):,}, med={sub.median():.0f}°)')

ax2.axvline(30, color='grey', ls=':', lw=1)
ax2.set_xlabel('Kagan angle (°)', fontsize=11)
ax2.set_ylabel('Cumulative %', fontsize=11)
ax2.set_title('CDF by SKHASH quality', fontsize=11)
ax2.legend(fontsize=9, framealpha=0.9)
ax2.set_xlim(0, 90)
ax2.set_ylim(0, 100)
ax2.grid(True, alpha=0.3, linewidth=0.5)

fig.suptitle('Kagan angle: SKHASH 1D vs HASH 1D  |  ge8 Y1+Y2', fontsize=12, y=1.01)
plt.tight_layout()
out = OUT_DIR / 'kagan_skhash1d_vs_hash1d.png'
fig.savefig(out, dpi=150, bbox_inches='tight')
plt.close(fig)
print(f'\nSaved → {out}')

# ── Summary table ──────────────────────────────────────────────────────────────
print('\nMedian Kagan angle by SKHASH quality:')
for q in QUAL_ORDER:
    sub = merged[merged['quality'] == q]['kagan']
    if len(sub):
        pct30 = (sub <= 30).mean() * 100
        print(f'  {q}: n={len(sub):5,}  median={sub.median():.1f}°  '
              f'%≤30°={pct30:.1f}%')
