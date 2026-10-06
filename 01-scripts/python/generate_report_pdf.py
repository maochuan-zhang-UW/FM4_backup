"""
generate_report_pdf.py
======================
Generate a multi-page PDF report comparing SKHASH 1D vs HASH 1D focal mechanisms.

Output: 02-data/G_FM/FM_comparison_report.pdf

Run with FM_ML env:
  /opt/miniconda3/envs/FM_ML/bin/python generate_report_pdf.py
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
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import FancyBboxPatch
import matplotlib.gridspec as gridspec

from fault_classification import classify_fault_ptb

# ── Paths ──────────────────────────────────────────────────────────────────────
SKHASH_DIR = Path('/Users/mczhang/Documents/GitHub/SKHASH/SKHASH7')
HASH_DIR   = Path('/Users/mczhang/Documents/GitHub/FM7/01-scripts/HASH')
DATA_DIR   = Path('/Users/mczhang/Documents/GitHub/FM7/02-data')
OUT_DIR    = DATA_DIR / 'G_FM'
OUT_DIR.mkdir(exist_ok=True)
OUT_PDF    = OUT_DIR / 'FM_comparison_report.pdf'

TITLE_COLOR  = '#1a3a5c'
HEADER_COLOR = '#2c6fad'
ROW_ALT      = '#eef4fb'
ROW_WHITE    = '#ffffff'
ACCENT       = '#e8f0f9'

# ── Helpers ────────────────────────────────────────────────────────────────────
def _sdr2ptb(s, d, r):
    s, d, r = np.radians(s), np.radians(d), np.radians(r)
    n  = np.array([-np.sin(d)*np.sin(s), np.sin(d)*np.cos(s), -np.cos(d)])
    sl = np.array([ np.cos(r)*np.cos(s)+np.sin(r)*np.cos(d)*np.sin(s),
                    np.cos(r)*np.sin(s)-np.sin(r)*np.cos(d)*np.cos(s),
                   -np.sin(r)*np.sin(d)])
    T = (n + sl) / np.sqrt(2)
    P = (n - sl) / np.sqrt(2)
    B = np.cross(T, P)
    return P, T, B

def kagan_angle(s1, d1, r1, s2, d2, r2):
    P1, T1, B1 = _sdr2ptb(s1, d1, r1)
    P2, T2, B2 = _sdr2ptb(s2, d2, r2)
    R2 = np.column_stack([P2, T2, B2])
    min_ang = np.inf
    for sp, st, sb in [(1,1,1), (-1,-1,1), (1,-1,-1), (-1,1,-1)]:
        R1v = np.column_stack([sp*P1, st*T1, sb*B1])
        cos_a = np.clip((np.trace(R1v.T @ R2) - 1) / 2, -1, 1)
        min_ang = min(min_ang, np.degrees(np.arccos(cos_a)))
    return min_ang

kagan_vec = np.vectorize(kagan_angle)

def make_table(ax, col_labels, row_labels, data, col_widths=None,
               highlight_rows=None, highlight_cols=None):
    """Draw a styled table on ax."""
    ax.axis('off')
    n_rows = len(row_labels)
    n_cols = len(col_labels)
    if col_widths is None:
        col_widths = [1.0 / n_cols] * n_cols

    total_rows = n_rows + 1  # +1 for header
    row_height = 1.0 / total_rows
    x_starts = np.cumsum([0] + col_widths[:-1])

    # Header
    for j, (label, x, w) in enumerate(zip(col_labels, x_starts, col_widths)):
        y = 1.0 - row_height
        ax.add_patch(FancyBboxPatch((x, y), w, row_height,
                     boxstyle='square,pad=0', linewidth=0.5,
                     edgecolor='white', facecolor=HEADER_COLOR,
                     transform=ax.transAxes, clip_on=False))
        ax.text(x + w/2, y + row_height/2, label,
                ha='center', va='center', fontsize=8.5, fontweight='bold',
                color='white', transform=ax.transAxes)

    # Rows
    for i, (rl, row) in enumerate(zip(row_labels, data)):
        y = 1.0 - (i + 2) * row_height
        bg = ROW_ALT if i % 2 == 0 else ROW_WHITE
        if highlight_rows and i in highlight_rows:
            bg = '#fff3cd'
        for j, (val, x, w) in enumerate(zip([rl] + list(row), [0] + list(x_starts), [col_widths[0]] + col_widths[1:])):
            fc = bg
            if highlight_cols and j in highlight_cols:
                fc = '#d4edda'
            ax.add_patch(FancyBboxPatch((x, y), w, row_height,
                         boxstyle='square,pad=0', linewidth=0.3,
                         edgecolor='#cccccc', facecolor=fc,
                         transform=ax.transAxes, clip_on=False))
            fw = 'bold' if j == 0 else 'normal'
            ax.text(x + w/2, y + row_height/2, str(val),
                    ha='center', va='center', fontsize=8, fontweight=fw,
                    color='#222222', transform=ax.transAxes)

# ── Load data ──────────────────────────────────────────────────────────────────
print('Loading data...')

sk_y1 = pd.read_csv(SKHASH_DIR / 'examples/hash3/OUT_Y1_ge8_1D/out.txt')
sk_y2 = pd.read_csv(SKHASH_DIR / 'examples/hash3/OUT_Y2_ge8_1D/out.txt')
sk_y1['cluster_idx'] = sk_y1['event_id']
sk_y2['cluster_idx'] = sk_y2['event_id'] + 1473
sk1d = pd.concat([sk_y1, sk_y2], ignore_index=True)

sk3d_y1 = pd.read_csv(SKHASH_DIR / 'examples/hash3/OUT_Y1_ge8/out.txt')
sk3d_y2 = pd.read_csv(SKHASH_DIR / 'examples/hash3/OUT_Y2_ge8/out.txt')
sk3d = pd.concat([sk3d_y1, sk3d_y2], ignore_index=True)

hash_rows = []
with open(HASH_DIR / 'hashout_ge8_1D_1.dat') as fh:
    for line in fh:
        p = line.split()
        if len(p) < 31: continue
        try:
            hash_rows.append(dict(cluster_idx=int(p[0]),
                                  strike=float(p[21]), dip=float(p[22]),
                                  rake=float(p[23]), quality=p[28]))
        except: continue
hash1d = pd.DataFrame(hash_rows)

merged = sk1d[['cluster_idx','strike','dip','rake','quality']].merge(
    hash1d[['cluster_idx','strike','dip','rake','quality']].rename(
        columns={'strike':'h_strike','dip':'h_dip','rake':'h_rake','quality':'h_qual'}),
    on='cluster_idx', how='inner')

print('Computing Kagan angles...')
merged['kagan'] = kagan_vec(merged['strike'], merged['dip'], merged['rake'],
                             merged['h_strike'], merged['h_dip'], merged['h_rake'])

# ── Build PDF ──────────────────────────────────────────────────────────────────
print('Building PDF...')

def add_page_header(fig, title, subtitle=''):
    fig.text(0.5, 0.97, title, ha='center', va='top',
             fontsize=15, fontweight='bold', color=TITLE_COLOR)
    if subtitle:
        fig.text(0.5, 0.935, subtitle, ha='center', va='top',
                 fontsize=10, color='#555555', style='italic')

def add_footer(fig, page_num):
    fig.text(0.5, 0.015, f'Axial Seamount Focal Mechanism Analysis  |  Page {page_num}',
             ha='center', va='bottom', fontsize=7.5, color='#888888')
    fig.text(0.97, 0.015, 'FM7 Pipeline', ha='right', va='bottom',
             fontsize=7.5, color='#aaaaaa')

with PdfPages(OUT_PDF) as pdf:

    # ══════════════════════════════════════════════════════════════
    # PAGE 1 — Title Page
    # ══════════════════════════════════════════════════════════════
    fig = plt.figure(figsize=(8.5, 11), facecolor='white')
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')

    # Top bar
    ax.add_patch(FancyBboxPatch((0, 0.78), 1.0, 0.22,
                 boxstyle='square,pad=0', facecolor=TITLE_COLOR,
                 transform=ax.transAxes, clip_on=False))
    ax.text(0.5, 0.915, 'SKHASH 1D vs HASH 1D',
            ha='center', va='center', fontsize=22, fontweight='bold',
            color='white', transform=ax.transAxes)
    ax.text(0.5, 0.845, 'Focal Mechanism Comparison Report',
            ha='center', va='center', fontsize=14,
            color='#b8d4f0', transform=ax.transAxes)
    ax.text(0.5, 0.800, 'Axial Seamount  |  Y1 + Y2  |  ≥8-station events',
            ha='center', va='center', fontsize=11,
            color='#d0e8ff', transform=ax.transAxes)

    # Summary boxes
    boxes = [
        ('SKHASH 3D', '14,597', 'Total FMs', '2,625', 'A+B (18%)'),
        ('SKHASH 1D', '14,685', 'Total FMs', '3,335', 'A+B (23%)'),
        ('HASH 1D',   '15,509', 'Total FMs', '9,075', 'A+B (59%)'),
    ]
    box_colors = ['#1565c0', '#2196f3', '#4caf50']
    for k, (name, total, tlabel, ab, ablabel) in enumerate(boxes):
        x = 0.08 + k * 0.305
        ax.add_patch(FancyBboxPatch((x, 0.57), 0.27, 0.19,
                     boxstyle='round,pad=0.01', linewidth=1.5,
                     edgecolor=box_colors[k], facecolor=ACCENT,
                     transform=ax.transAxes, clip_on=False))
        ax.text(x + 0.135, 0.745, name, ha='center', va='center',
                fontsize=11, fontweight='bold', color=box_colors[k],
                transform=ax.transAxes)
        ax.text(x + 0.135, 0.695, total, ha='center', va='center',
                fontsize=18, fontweight='bold', color=TITLE_COLOR,
                transform=ax.transAxes)
        ax.text(x + 0.135, 0.663, tlabel, ha='center', va='center',
                fontsize=8.5, color='#555555', transform=ax.transAxes)
        ax.plot([x+0.02, x+0.25], [0.648, 0.648],
                linewidth=0.7, color='#cccccc', transform=ax.transAxes)
        ax.text(x + 0.135, 0.618, ab, ha='center', va='center',
                fontsize=14, fontweight='bold', color=box_colors[k],
                transform=ax.transAxes)
        ax.text(x + 0.135, 0.588, ablabel, ha='center', va='center',
                fontsize=8.5, color='#555555', transform=ax.transAxes)

    # Contents
    ax.text(0.08, 0.535, 'Report Contents', fontsize=11, fontweight='bold',
            color=TITLE_COLOR, transform=ax.transAxes)
    ax.add_patch(FancyBboxPatch((0.08, 0.12), 0.84, 0.40,
                 boxstyle='round,pad=0.01', linewidth=0.8,
                 edgecolor='#cccccc', facecolor='#f8fbff',
                 transform=ax.transAxes, clip_on=False))
    contents = [
        ('1.', 'Event Counts & Quality Distribution', 'Breakdown of FM counts and A/B/C/D grades across all three runs'),
        ('2.', 'Quality Criteria Comparison', 'Why HASH produces more A+B events: 2-criterion vs 4-criterion grading'),
        ('3.', 'Kagan Angle Distribution', 'Overall histogram coloured by SKHASH quality grade'),
        ('4.', 'Kagan Angle Statistics', 'Median and %≤30° by SKHASH quality, HASH quality, and combined A+B subsets'),
        ('5.', 'Key Takeaways', 'Summary of findings and practical recommendations'),
        ('D&M', 'Data & Methods (Part 1)', 'Data collection, waveform processing, DL polarity picking'),
        ('D&M', 'Data & Methods (Part 2)', 'Polarity augmentation: CC transfer and hierarchical clustering'),
        ('D&M', 'Data & Methods (Part 3)', 'FM calculation: velocity models, station corrections, Y1/Y2 differences'),
        ('6.', 'FM Analysis (Part 1)', 'Fault-type distribution, depth structure and rake statistics'),
        ('7.', 'FM Analysis (Part 2)', 'Strike rose diagrams, spatial region and Y1/Y2 breakdown'),
        ('8.', 'FM Analysis (Part 3)', 'Geologic interpretation and suggested follow-up analyses'),
    ]
    for i, (num, title, desc) in enumerate(contents):
        y = 0.48 - i * 0.065
        ax.text(0.12, y, num, fontsize=9, fontweight='bold', color=HEADER_COLOR,
                transform=ax.transAxes)
        ax.text(0.165, y, title, fontsize=9, fontweight='bold', color='#222222',
                transform=ax.transAxes)
        ax.text(0.165, y - 0.025, desc, fontsize=8, color='#666666',
                transform=ax.transAxes)

    ax.text(0.5, 0.07, 'Generated by FM7 Pipeline', ha='center', fontsize=8,
            color='#aaaaaa', transform=ax.transAxes)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    # ══════════════════════════════════════════════════════════════
    # PAGE 2 — Event Counts & Quality Distribution
    # ══════════════════════════════════════════════════════════════
    fig = plt.figure(figsize=(8.5, 11), facecolor='white')
    add_page_header(fig, '1. Event Counts & Quality Distribution',
                    'All ge8 events, Y1+Y2 combined')
    add_footer(fig, 2)

    gs = gridspec.GridSpec(2, 1, figure=fig, top=0.88, bottom=0.08,
                           hspace=0.45, left=0.08, right=0.92)

    # Table
    ax_t = fig.add_subplot(gs[0])
    QUAL_ORDER = ['A', 'B', 'C', 'D', 'A+B', 'Total']
    def qcounts(df):
        vc = df['quality'].value_counts()
        total = len(df)
        ab = len(df[df['quality'].isin(['A','B'])])
        return [f"{vc.get(q,0):,} ({vc.get(q,0)/total*100:.0f}%)" for q in ['A','B','C','D']] + \
               [f"{ab:,} ({ab/total*100:.0f}%)", f"{total:,}"]

    col_labels = ['Grade', 'SKHASH 3D', 'SKHASH 1D', 'HASH 1D']
    row_labels = ['A', 'B', 'C', 'D', 'A+B', 'Total']

    sk3d_total = len(sk3d); sk3d_vc = sk3d['quality'].value_counts()
    sk1d_total = len(sk1d); sk1d_vc = sk1d['quality'].value_counts()
    h1d_total  = len(hash1d); h1d_vc  = hash1d['quality'].value_counts()

    table_data = []
    for q in ['A','B','C','D']:
        table_data.append([
            f"{sk3d_vc.get(q,0):,}  ({sk3d_vc.get(q,0)/sk3d_total*100:.0f}%)",
            f"{sk1d_vc.get(q,0):,}  ({sk1d_vc.get(q,0)/sk1d_total*100:.0f}%)",
            f"{h1d_vc.get(q,0):,}  ({h1d_vc.get(q,0)/h1d_total*100:.0f}%)",
        ])
    sk3d_ab = sk3d_vc.get('A',0)+sk3d_vc.get('B',0)
    sk1d_ab = sk1d_vc.get('A',0)+sk1d_vc.get('B',0)
    h1d_ab  = h1d_vc.get('A',0)+h1d_vc.get('B',0)
    table_data.append([
        f"{sk3d_ab:,}  ({sk3d_ab/sk3d_total*100:.0f}%)",
        f"{sk1d_ab:,}  ({sk1d_ab/sk1d_total*100:.0f}%)",
        f"{h1d_ab:,}  ({h1d_ab/h1d_total*100:.0f}%)",
    ])
    table_data.append([f"{sk3d_total:,}", f"{sk1d_total:,}", f"{h1d_total:,}"])

    make_table(ax_t, col_labels, row_labels, table_data,
               col_widths=[0.12, 0.29, 0.29, 0.30],
               highlight_rows=[4], highlight_cols=[])
    ax_t.set_title('Quality Grade Summary', fontsize=10, fontweight='bold',
                   color=TITLE_COLOR, pad=8)

    # Bar chart
    ax_b = fig.add_subplot(gs[1])
    quals = ['A', 'B', 'C', 'D']
    x = np.arange(len(quals))
    w = 0.25
    colors_run = ['#1565c0', '#2196f3', '#4caf50']
    bars1 = ax_b.bar(x - w, [sk3d_vc.get(q,0) for q in quals], w,
                     color=colors_run[0], label='SKHASH 3D', alpha=0.85)
    bars2 = ax_b.bar(x,     [sk1d_vc.get(q,0) for q in quals], w,
                     color=colors_run[1], label='SKHASH 1D', alpha=0.85)
    bars3 = ax_b.bar(x + w, [h1d_vc.get(q,0)  for q in quals], w,
                     color=colors_run[2], label='HASH 1D',   alpha=0.85)
    ax_b.set_xticks(x); ax_b.set_xticklabels(quals, fontsize=11)
    ax_b.set_xlabel('Quality Grade', fontsize=10)
    ax_b.set_ylabel('Number of Events', fontsize=10)
    ax_b.set_title('Quality Grade Distribution', fontsize=10, fontweight='bold',
                   color=TITLE_COLOR)
    ax_b.legend(fontsize=9, framealpha=0.9)
    ax_b.grid(True, axis='y', alpha=0.3, linewidth=0.5)
    ax_b.tick_params(labelsize=9)
    for bars in [bars1, bars2, bars3]:
        for bar in bars:
            h = bar.get_height()
            if h > 200:
                ax_b.text(bar.get_x() + bar.get_width()/2, h + 50,
                          f'{h:,}', ha='center', va='bottom', fontsize=6.5, color='#333333')

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    # ══════════════════════════════════════════════════════════════
    # PAGE 3 — Quality Criteria Comparison
    # ══════════════════════════════════════════════════════════════
    fig = plt.figure(figsize=(8.5, 11), facecolor='white')
    add_page_header(fig, '2. Quality Criteria Comparison',
                    'Why HASH produces more A+B events')
    add_footer(fig, 3)

    # HASH table
    ax1 = fig.add_axes([0.06, 0.60, 0.88, 0.25])
    col_labels_h = ['Grade', 'prob  >', 'var_avg  ≤']
    row_labels_h = ['A', 'B', 'C', 'D']
    data_h = [['0.80', '25°'], ['0.60', '35°'], ['0.50', '45°'], ['—', '—']]
    make_table(ax1, col_labels_h, row_labels_h, data_h,
               col_widths=[0.12, 0.44, 0.44])
    ax1.set_title('HASH 1D Quality Criteria  (2 metrics, both must pass)',
                  fontsize=10, fontweight='bold', color=TITLE_COLOR, pad=8)

    # SKHASH table
    ax2 = fig.add_axes([0.06, 0.27, 0.88, 0.28])
    col_labels_s = ['Grade', 'prob  ≥', 'var_avg  ≤', 'mfrac  ≤', 'stdr  ≥']
    row_labels_s = ['A', 'B', 'C', 'D']
    data_s = [['0.80', '25°', '0.15', '0.50'],
              ['0.60', '35°', '0.20', '0.40'],
              ['0.50', '45°', '0.30', '0.30'],
              ['0.0',  '∞',   '∞',    '0.0']]
    make_table(ax2, col_labels_s, row_labels_s, data_s,
               col_widths=[0.12, 0.22, 0.22, 0.22, 0.22],
               highlight_cols=[3, 4])
    ax2.set_title('SKHASH Quality Criteria  (4 metrics — final grade = worst of all four)',
                  fontsize=10, fontweight='bold', color=TITLE_COLOR, pad=8)

    # Annotations
    fig.text(0.08, 0.24, 'Metric Definitions', fontsize=10, fontweight='bold',
             color=TITLE_COLOR)
    defs = [
        ('prob', 'Probability that the solution is within the acceptable cone (0–1).  Both algorithms use identical thresholds.'),
        ('var_avg', 'Mean of fault-plane and auxiliary-plane RMS uncertainties (degrees).  Both algorithms use identical thresholds.'),
        ('mfrac', 'Weighted fraction of polarity misfits.  SKHASH-only constraint.  Requires ≤15% for A, ≤20% for B.'),
        ('stdr', 'Station Distribution Ratio — measures azimuthal coverage (0–1).  SKHASH-only constraint.\n'
                 '         At Axial Seamount, the compact OBS array limits azimuthal coverage, causing many events to\n'
                 '         fail the stdr threshold and be downgraded from A/B to C/D.  This is the primary driver\n'
                 '         of the ~6,000-event A+B gap between HASH and SKHASH.'),
    ]
    y0 = 0.205
    for metric, desc in defs:
        fig.text(0.08, y0, f'• {metric}:', fontsize=8.5, fontweight='bold', color='#333333')
        fig.text(0.185, y0, desc, fontsize=8, color='#444444', wrap=True)
        y0 -= 0.052

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    # ══════════════════════════════════════════════════════════════
    # PAGE 4 — Kagan Angle Histogram + CDF
    # ══════════════════════════════════════════════════════════════
    fig = plt.figure(figsize=(8.5, 11), facecolor='white')
    add_page_header(fig, '3. Kagan Angle Distribution',
                    'SKHASH 1D vs HASH 1D  |  all 14,685 matched events')
    add_footer(fig, 4)

    QUAL_COLORS = {'A': '#d62728', 'B': '#ff7f0e', 'C': '#2ca02c', 'D': '#9467bd'}
    bins = np.arange(0, 91, 5)

    gs2 = gridspec.GridSpec(2, 1, figure=fig, top=0.87, bottom=0.08,
                            hspace=0.42, left=0.10, right=0.92)

    # Histogram
    ax_h = fig.add_subplot(gs2[0])
    data_by_qual = [merged[merged['quality'] == q]['kagan'].values for q in ['A','B','C','D']]
    ax_h.hist(data_by_qual, bins=bins, stacked=True,
              color=[QUAL_COLORS[q] for q in ['A','B','C','D']],
              label=[f'SKHASH {q} (n={len(d):,})' for q, d in zip(['A','B','C','D'], data_by_qual)],
              edgecolor='none', alpha=0.88)
    med = merged['kagan'].median()
    ax_h.axvline(med, color='black', ls='--', lw=1.4,
                 label=f'Median = {med:.1f}°')
    ax_h.axvline(30, color='gray', ls=':', lw=1.0, label='30° threshold')
    ax_h.set_xlabel('Kagan Angle (°)', fontsize=10)
    ax_h.set_ylabel('Count', fontsize=10)
    ax_h.set_title('Histogram by SKHASH Quality Grade', fontsize=10,
                   fontweight='bold', color=TITLE_COLOR)
    ax_h.legend(fontsize=8.5, framealpha=0.9)
    ax_h.set_xlim(0, 90)
    ax_h.grid(True, alpha=0.3, linewidth=0.5)

    # CDF
    ax_c = fig.add_subplot(gs2[1])
    for q in ['A', 'B', 'C', 'D']:
        sub = merged[merged['quality'] == q]['kagan']
        if len(sub) == 0: continue
        sorted_k = np.sort(sub)
        cdf = np.arange(1, len(sorted_k)+1) / len(sorted_k)
        ax_c.plot(sorted_k, cdf * 100, color=QUAL_COLORS[q], lw=1.8,
                  label=f'{q} (n={len(sub):,}, med={sub.median():.0f}°)')
    ax_c.axvline(30, color='gray', ls=':', lw=1.0, label='30° threshold')
    ax_c.set_xlabel('Kagan Angle (°)', fontsize=10)
    ax_c.set_ylabel('Cumulative %', fontsize=10)
    ax_c.set_title('CDF by SKHASH Quality Grade', fontsize=10,
                   fontweight='bold', color=TITLE_COLOR)
    ax_c.legend(fontsize=8.5, framealpha=0.9)
    ax_c.set_xlim(0, 90); ax_c.set_ylim(0, 100)
    ax_c.grid(True, alpha=0.3, linewidth=0.5)

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    # ══════════════════════════════════════════════════════════════
    # PAGE 5 — Kagan Angle Statistics Tables
    # ══════════════════════════════════════════════════════════════
    fig = plt.figure(figsize=(8.5, 11), facecolor='white')
    add_page_header(fig, '4. Kagan Angle Statistics',
                    'Agreement between SKHASH 1D and HASH 1D')
    add_footer(fig, 5)

    # Overall stats
    ax_ov = fig.add_axes([0.06, 0.74, 0.88, 0.12])
    ov_labels = ['Count', 'Min', '25th pct', 'Median', 'Mean', '75th pct', 'Max', '%  ≤ 30°']
    k = merged['kagan']
    ov_vals = [f'{len(k):,}', f'{k.min():.1f}°', f'{k.quantile(0.25):.1f}°',
               f'{k.median():.1f}°', f'{k.mean():.1f}°', f'{k.quantile(0.75):.1f}°',
               f'{k.max():.1f}°', f'{(k<=30).mean()*100:.1f}%']
    make_table(ax_ov, ['Statistic'] + ov_labels, ['All'],
               [ov_vals], col_widths=[0.10]+[0.1125]*8)
    ax_ov.set_title('Overall Statistics  (all 14,685 matched events)',
                    fontsize=10, fontweight='bold', color=TITLE_COLOR, pad=8)

    # By SKHASH quality
    ax_sk = fig.add_axes([0.06, 0.50, 0.88, 0.19])
    sk_rows = []
    sk_rl   = []
    for q in ['A', 'B', 'C', 'D']:
        sub = merged[merged['quality'] == q]['kagan']
        sk_rl.append(q)
        sk_rows.append([f'{len(sub):,}',
                        f'{sub.median():.1f}°', f'{sub.mean():.1f}°',
                        f'{sub.quantile(0.25):.1f}°', f'{sub.quantile(0.75):.1f}°',
                        f'{(sub<=30).mean()*100:.1f}%'])
    make_table(ax_sk, ['Grade', 'n', 'Median', 'Mean', '25th pct', '75th pct', '% ≤ 30°'],
               sk_rl, sk_rows, col_widths=[0.10, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15])
    ax_sk.set_title('By SKHASH Quality Grade',
                    fontsize=10, fontweight='bold', color=TITLE_COLOR, pad=8)

    # By HASH quality
    ax_hq = fig.add_axes([0.06, 0.28, 0.88, 0.19])
    hq_rows = []
    hq_rl   = []
    for q in ['A', 'B', 'C', 'D']:
        sub = merged[merged['h_qual'] == q]['kagan']
        hq_rl.append(q)
        hq_rows.append([f'{len(sub):,}',
                        f'{sub.median():.1f}°', f'{sub.mean():.1f}°',
                        f'{sub.quantile(0.25):.1f}°', f'{sub.quantile(0.75):.1f}°',
                        f'{(sub<=30).mean()*100:.1f}%'])
    make_table(ax_hq, ['Grade', 'n', 'Median', 'Mean', '25th pct', '75th pct', '% ≤ 30°'],
               hq_rl, hq_rows, col_widths=[0.10, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15])
    ax_hq.set_title('By HASH Quality Grade',
                    fontsize=10, fontweight='bold', color=TITLE_COLOR, pad=8)

    # Combined A+B
    ax_ab = fig.add_axes([0.06, 0.10, 0.88, 0.14])
    subsets = [
        ('SKHASH A+B only', merged[merged['quality'].isin(['A','B'])]['kagan']),
        ('HASH A+B only',   merged[merged['h_qual'].isin(['A','B'])]['kagan']),
        ('Both A+B',        merged[merged['quality'].isin(['A','B']) &
                                   merged['h_qual'].isin(['A','B'])]['kagan']),
    ]
    ab_rows = []
    ab_rl   = []
    for name, sub in subsets:
        ab_rl.append(name)
        ab_rows.append([f'{len(sub):,}',
                        f'{sub.median():.1f}°', f'{sub.mean():.1f}°',
                        f'{(sub<=30).mean()*100:.1f}%'])
    make_table(ax_ab, ['Subset', 'n', 'Median', 'Mean', '% ≤ 30°'],
               ab_rl, ab_rows, col_widths=[0.30, 0.175, 0.175, 0.175, 0.175],
               highlight_rows=[2])
    ax_ab.set_title('High-confidence A+B Subsets',
                    fontsize=10, fontweight='bold', color=TITLE_COLOR, pad=8)

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    # ══════════════════════════════════════════════════════════════
    # PAGE 6 — Data & Methods: Data Collection + Waveform Processing
    # ══════════════════════════════════════════════════════════════
    fig = plt.figure(figsize=(8.5, 11), facecolor='white')
    add_page_header(fig, 'Data & Methods  (Part 1 of 3)',
                    'Data Collection, Waveform Processing & DL Polarity Picking')
    add_footer(fig, 6)

    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')

    def section_title(ax, y, text):
        ax.add_patch(FancyBboxPatch((0.05, y - 0.005), 0.90, 0.030,
                     boxstyle='square,pad=0', facecolor=HEADER_COLOR,
                     transform=ax.transAxes, clip_on=False))
        ax.text(0.07, y + 0.009, text, fontsize=9.5, fontweight='bold',
                color='white', transform=ax.transAxes)

    def bullet(ax, y, label, value, indent=0.07):
        ax.text(indent, y, f'• {label}:', fontsize=8.5, fontweight='bold',
                color='#1a3a5c', transform=ax.transAxes)
        ax.text(indent + 0.22, y, value, fontsize=8.5, color='#333333',
                transform=ax.transAxes)

    def para(ax, y, text, indent=0.07, fs=8.5):
        ax.text(indent, y, text, fontsize=fs, color='#333333',
                transform=ax.transAxes, wrap=True)

    # ── Section 1: Data Overview ──
    section_title(ax, 0.88, '1. Data Overview')
    items1 = [
        ('Study area',    'Axial Seamount, NE Pacific  (~45.9–46.0°N, 129.95–130.06°W)'),
        ('Data period',   'Y1: 2022–2023  |  Y2: 2023–2024  (two separate OBS deployments)'),
        ('Network',       'OO network (Ocean Observatories Initiative) + temporary 2F/TF array'),
        ('Y1 stations',   '22 total — 7 legacy OO (AS1 AS2 CC1 EC1 EC2 EC3 ID1) + 15 TF (01A–15A)'),
        ('Y2 stations',   '21 total — 7 legacy OO + 14 TF (01B–14B); 15A not deployed in Y2'),
        ('Catalog',       'MLDD double-difference relocations  (axial.Y1/Y2.mldd.loc.260317)'),
        ('Event filter',  '≥8 stations with confident polarity picks; ≥15 picks required for FM'),
        ('Y1 FM events',  '1,473 events  |  Y2 FM events: 13,212–14,036 events'),
    ]
    y = 0.845
    for lbl, val in items1:
        bullet(ax, y, lbl, val)
        y -= 0.031

    # ── Section 2: Waveform Processing ──
    section_title(ax, 0.600, '2. Waveform Processing')
    items2 = [
        ('Source data',      'miniSEED files; 2F stations from FM4 datamseed; OO stations from Axial-AutoLocate'),
        ('Sampling rate',    '200 Hz (HHZ / EHZ / ELZ channels)'),
        ('Window',           '−0.32 s to +1.00 s relative to P pick  (264 samples total)'),
        ('Bandpass filter',  '4–50 Hz, 4th-order zero-phase Butterworth (SOS)'),
        ('Noise window',     '−0.70 s to −0.10 s before P pick  → noise amplitude'),
        ('P window',         '−0.05 s to +0.25 s around P pick  → P-wave amplitude'),
        ('S window',         '−0.10 s to +0.60 s around S pick  → S-wave amplitude'),
        ('Amplitude stored', '[noise_amp, S_amp, P_amp]  →  02-data/A_ID/waves_{year}_ge8.h5'),
    ]
    y = 0.560
    for lbl, val in items2:
        bullet(ax, y, lbl, val)
        y -= 0.031

    # ── Section 3: DL Polarity Picking ──
    section_title(ax, 0.310, '3. Deep Learning Polarity Picking')
    items3 = [
        ('Model',           'PolarPicker  (PolarPicker_unified_TMSF_001.keras, Keras/TensorFlow)'),
        ('Input rate',      '100 Hz (downsampled from 200 Hz)'),
        ('Input window',    '64 samples = 0.64 s starting at P pick; peak-normalized'),
        ('Output',          'Two class probabilities: P(Down), P(Up)'),
        ('Polarity',        '+1 (Up) if argmax=Up; −1 (Down) if argmax=Down; 0 if no waveform'),
        ('Confidence',      'max(class_probs); values 0.5–1.0'),
        ('Conf threshold',  '≥0.70 = "high confidence"; ≥0.50 = acceptable for augmentation'),
        ('Raw output',      '02-data/E_Po/dl_polarity_{year}_ge8.h5  per station per event'),
    ]
    y = 0.270
    for lbl, val in items3:
        bullet(ax, y, lbl, val)
        y -= 0.031

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    # ══════════════════════════════════════════════════════════════
    # PAGE 7 — Data & Methods: Polarity Augmentation
    # ══════════════════════════════════════════════════════════════
    fig = plt.figure(figsize=(8.5, 11), facecolor='white')
    add_page_header(fig, 'Data & Methods  (Part 2 of 3)',
                    'Polarity Augmentation — Cross-Correlation Transfer & Clustering')
    add_footer(fig, 7)

    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')

    # ── Section 4: CC Polarity Transfer ──
    section_title(ax, 0.88, '4. Cross-Correlation Polarity Transfer')
    para(ax, 0.845,
         'For target events lacking a polarity at a given station, nearby template events with high-confidence DL\n'
         'picks are used to transfer polarities via waveform cross-correlation.')
    items4 = [
        ('Template conf min', '≥0.70 (DL confidence required to use event as transfer template)'),
        ('CC threshold',      '≥0.80  (normalized cross-correlation coefficient required to transfer)'),
        ('CC window',         '64 samples at 200 Hz = 0.32 s centered on P pick'),
        ('Max lag',           '10 samples = 50 ms'),
        ('Spatial search',    '≤3.0 km from target event'),
        ('Transferred conf',  'CC × template_confidence (lower bound 0.50)'),
        ('Output',            'Augmented polarity added to target event at that station'),
    ]
    y = 0.795
    for lbl, val in items4:
        bullet(ax, y, lbl, val, indent=0.09)
        y -= 0.032

    # ── Section 5: Hierarchical Clustering ──
    section_title(ax, 0.555, '5. Hierarchical Clustering Polarity Augmentation')
    para(ax, 0.520,
         'Waveforms are grouped into clusters by cross-correlation similarity. Within each cluster a majority vote\n'
         'assigns consistent polarities to events lacking reliable DL picks.')
    items5 = [
        ('CC threshold (clustering)', '≥0.70  (UPGMA linkage on pairwise CC distance d = 1 − CC)'),
        ('CC threshold (transfer)',   '≥0.80  (singleton to cluster transfer)'),
        ('Min cluster size',          '3 events (singletons retain original picks)'),
        ('Vote weight',               'DL confidence (≥0.50); winning polarity = weighted majority'),
        ('Spatial search',            '≤3.0 km radius'),
        ('CC window',                 '64 samples at 200 Hz = 0.32 s'),
        ('Max lag',                   '10 samples = 50 ms'),
        ('Output',                    '02-data/E_Po/dl_polarity_{year}_ge8_aug.h5  (augmented)'),
    ]
    y = 0.470
    for lbl, val in items5:
        bullet(ax, y, lbl, val, indent=0.09)
        y -= 0.032

    # ── Section 6: Polarity HDF5 structure ──
    section_title(ax, 0.225, '6. Final Polarity Input to FM Algorithms')
    items6 = [
        ('File',             'dl_polarity_{year}_ge8_aug.h5  (post-augmentation)'),
        ('Datasets',         '/event_ids  (int64),  /polarity/<STA>  (int8),  /confidence/<STA>  (float32)'),
        ('Polarity values',  '+1 = Up,  −1 = Down,  0 = no pick'),
        ('Conf filter',      'Only picks with confidence ≥0.50 used; Y1/Y2 merged by cluster_idx'),
        ('Min picks',        '≥15 polarity picks required per event to enter FM calculation'),
        ('S/P amplitudes',   'From waves_{year}_ge8.h5  →  [noise, S_amp, P_amp] per station'),
        ('Amp noise filter', 'SNR ≥2  (ratmin=2); noise floor floored at 0.01 to avoid zero division'),
    ]
    y = 0.185
    for lbl, val in items6:
        bullet(ax, y, lbl, val, indent=0.09)
        y -= 0.032

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    # ══════════════════════════════════════════════════════════════
    # PAGE 8 — Data & Methods: FM Calculation (SKHASH & HASH)
    # ══════════════════════════════════════════════════════════════
    fig = plt.figure(figsize=(8.5, 11), facecolor='white')
    add_page_header(fig, 'Data & Methods  (Part 3 of 3)',
                    'Focal Mechanism Calculation — SKHASH & HASH Parameters')
    add_footer(fig, 8)

    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')

    # ── Shared parameters table ──
    section_title(ax, 0.88, '7. Shared FM Algorithm Parameters  (SKHASH & HASH identical)')

    ax_tp = fig.add_axes([0.06, 0.62, 0.88, 0.240])
    shared_cols = ['Parameter', 'Value', 'Description']
    shared_rl   = ['npolmin', 'dang', 'nmc', 'maxout', 'ratmin', 'badfrac',
                   'qbadfrac', 'delmax', 'cangle', 'prob_max']
    shared_data = [
        ['15', 'Min polarity picks required per event'],
        ['5°', 'Grid search spacing for strike/dip/rake'],
        ['30', 'Monte Carlo trials (location perturbation)'],
        ['300', 'Max output mechanisms per event'],
        ['2', 'Min S/P amplitude SNR to use ratio'],
        ['0.10', 'Max allowed bad-polarity fraction'],
        ['0.30', 'Amplitude noise threshold (log₁₀ scale)'],
        ['25 km', 'Max source–receiver distance'],
        ['45°', 'Probability cone half-angle'],
        ['0.75', 'Multiple-solution probability threshold'],
    ]
    make_table(ax_tp, shared_cols, shared_rl, shared_data,
               col_widths=[0.20, 0.15, 0.65])
    ax_tp.set_title('', fontsize=1)

    # ── Velocity models ──
    section_title(ax, 0.590, '8. Velocity Models')
    vmod_items = [
        ('1D model',       'velmod_axial1D.txt — single 1D gradient model for all stations'),
        ('Depth range',    '0.0–40.0 km  (60 nodes; Vp = 2.08 km/s at surface → 8.20 km/s at 40 km)'),
        ('Shallow grad.',  '0.0–0.5 km: 2.08 → 3.83 km/s  (rapid increase through shallow crust)'),
        ('Mid crust',      '1.0–3.6 km: 4.93 → 6.00 km/s  (smooth gradient)'),
        ('Deep',           '3.6–6.2 km: 6.00 → 6.82 km/s  |  Half-space: 8.20 km/s at 40 km'),
        ('3D models',      '15 station/region-specific models used in SKHASH 3D run:'),
        ('',               '  OO stations: velmod_AXAS11, AXAS21, AXCC11, AXEC11, AXEC21, AXEC31, AXID11'),
        ('',               '  Regional: velmod_W11, W21, E11, E21, E31, E41, S11, velmod21'),
        ('Vp/Vs ratio',    'Fixed at √3 (Poisson solid; Vp/Vs ≈ 1.732) for both algorithms'),
        ('Station corr.',  'All corrections = 0.000 s  (north3.statcor.txt — no timing offsets applied)'),
    ]
    y = 0.550
    for lbl, val in vmod_items:
        if lbl:
            bullet(ax, y, lbl, val, indent=0.07)
        else:
            para(ax, y, val, indent=0.09)
        y -= 0.031

    # ── Y1 vs Y2 differences ──
    section_title(ax, 0.240, '9. Y1 vs Y2 Deployment Differences')

    ax_y = fig.add_axes([0.06, 0.06, 0.88, 0.155])
    y12_cols = ['Aspect', 'Y1  (2022–2023)', 'Y2  (2023–2024)']
    y12_rl   = ['TF stations', 'Station names', 'Phase file', 'Station file',
                'SKHASH output', 'Event IDs', 'FM events']
    y12_data = [
        ['15  (01A–15A)',         '14  (01B–14B)  — 15A not deployed'],
        ['01A … 15A + OO',        '01B … 14B + OO  (B-suffix year 2)'],
        ['north2_Y1_ge8.txt',     'north2_Y2_ge8.txt'],
        ['scsn.stations_Y1_ge8', 'scsn.stations_Y2_ge8'],
        ['OUT_Y1_ge8_1D/',        'OUT_Y2_ge8_1D/  (event_id restarts at 1)'],
        ['cluster_idx 1–1,473',   'cluster_idx 1,474–15,509  (offset +1,473 for matching)'],
        ['1,473',                 '13,212–14,036  (larger due to 2023 eruption swarm)'],
    ]
    make_table(ax_y, y12_cols, y12_rl, y12_data,
               col_widths=[0.22, 0.39, 0.39])
    ax_y.set_title('', fontsize=1)

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    # ══════════════════════════════════════════════════════════════
    # PAGE 9 — Key Takeaways (renumbered)
    # ══════════════════════════════════════════════════════════════
    fig = plt.figure(figsize=(8.5, 11), facecolor='white')
    add_page_header(fig, '5. Key Takeaways', '')
    add_footer(fig, 9)

    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')

    takeaways = [
        ('1', '#2196f3', 'High-quality solutions agree well.',
         'Events rated A+B by both SKHASH and HASH have a median Kagan angle of only 15.7°,\n'
         '     with 80% within 30°. This confirms both algorithms converge on the same mechanism\n'
         '     for well-constrained events.'),
        ('2', '#4caf50', 'HASH is more permissive — but captures genuine signal.',
         'HASH A+B events (9,075) outnumber SKHASH A+B (3,335) by ~3×. The "extra" HASH A+B\n'
         '     events show a median Kagan angle of 20° vs SKHASH\'s 17°, suggesting they are\n'
         '     slightly less certain but still useful solutions.'),
        ('3', '#ff9800', 'stdr is the dominant SKHASH penalty at Axial Seamount.',
         'SKHASH\'s two additional criteria (mfrac, stdr) downgrade solutions HASH rates as A/B.\n'
         '     The compact OBS array at Axial Seamount limits azimuthal coverage, causing stdr to\n'
         '     fall below the A/B thresholds for many events — even when prob and var_avg pass.'),
        ('4', '#9c27b0', 'SKHASH quality grade is a strong predictor of agreement.',
         'Median Kagan angle rises monotonically A→D: 14.4° → 17.2° → 21.8° → 33.8°.\n'
         '     SKHASH grade D events diverge substantially; use with caution.'),
        ('5', '#f44336', 'Practical recommendation.',
         'Use the Both A+B intersection (2,720 events) for highest-confidence analysis.\n'
         '     HASH A+B alone (8,593 events) extends coverage for spatial/statistical analysis\n'
         '     with acceptable uncertainty (~20° median Kagan). Avoid grade D events for\n'
         '     mechanism-sensitive interpretations.'),
    ]

    y0 = 0.84
    for num, color, title, body in takeaways:
        # Numbered circle
        circle = plt.Circle((0.09, y0 + 0.005), 0.028, color=color,
                             transform=ax.transAxes, clip_on=False)
        ax.add_patch(circle)
        ax.text(0.09, y0 + 0.005, num, ha='center', va='center',
                fontsize=13, fontweight='bold', color='white',
                transform=ax.transAxes)

        ax.add_patch(FancyBboxPatch((0.13, y0 - 0.075), 0.80, 0.095,
                     boxstyle='round,pad=0.01', linewidth=1.0,
                     edgecolor=color, facecolor='#f9f9f9',
                     transform=ax.transAxes, clip_on=False))
        ax.text(0.155, y0 + 0.003, title, fontsize=10, fontweight='bold',
                color=color, transform=ax.transAxes)
        ax.text(0.155, y0 - 0.028, body, fontsize=8.5, color='#333333',
                va='top', transform=ax.transAxes)

        y0 -= 0.142

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    # ══════════════════════════════════════════════════════════════
    # PAGE 10 — Research Context & Brainstorm
    # ══════════════════════════════════════════════════════════════
    fig = plt.figure(figsize=(8.5, 11), facecolor='white')
    add_page_header(fig, 'Research Context & Exploration Plan  (Part 1 of 2)',
                    'Comparison with Zhang et al. 2015–2021 Study & New Opportunities')
    add_footer(fig, 10)

    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')

    # ── Context comparison table ──
    ax.add_patch(FancyBboxPatch((0.04, 0.875), 0.92, 0.024,
                 boxstyle='square,pad=0', facecolor=TITLE_COLOR,
                 transform=ax.transAxes, clip_on=False))
    ax.text(0.06, 0.886, 'STUDY COMPARISON: Zhang et al. (2015–2021)  vs  This Study (2022–2024)',
            fontsize=9, fontweight='bold', color='white', transform=ax.transAxes)

    ax_cmp = fig.add_axes([0.04, 0.655, 0.92, 0.205])
    cmp_cols = ['Aspect', 'Zhang et al. (2015–2021)', 'This Study (2022–2024)']
    cmp_rl = ['Period', 'Network', 'Method', 'FM catalog', 'Type', 'Velocity', 'Key event']
    cmp_data = [
        ['Jan 2015 – Dec 2021  (7 years)', 'Sep 2022 – Sep 2024  (2 years)'],
        ['7 OOI cabled stations', '21–22 stations: 7 OO + 14–15 TF'],
        ['CC polarity + SVD + hierarchical cluster', 'DL (PolarPicker) + CC augment + cluster'],
        ['3,873 composite (23,221 events)', '14,597–15,509 individual events'],
        ['Composite (cluster shares 1 FM)', 'Individual (each event has own FM)'],
        ['1-D + 3-D (Baillard et al. 2019)', '1-D axial + 3-D 15-model station-specific'],
        ['April 2015 eruption (north rift)', '2023 eruption (Y1 Sep 2022–Aug 2023)'],
    ]
    make_table(ax_cmp, cmp_cols, cmp_rl, cmp_data, col_widths=[0.20, 0.40, 0.40])
    ax_cmp.set_title('', fontsize=1)

    # ── What this study adds ──
    ax.add_patch(FancyBboxPatch((0.04, 0.626), 0.92, 0.022,
                 boxstyle='square,pad=0', facecolor='#2c6fad',
                 transform=ax.transAxes, clip_on=False))
    ax.text(0.06, 0.636, 'WHAT THIS STUDY ADDS  (structural advances)',
            fontsize=9, fontweight='bold', color='white', transform=ax.transAxes)

    advances = [
        ('4× more events', '15,509 individual FMs vs 3,873 composite — resolves finer spatial/temporal structure'),
        ('3× more stations', '21 vs 7 stations → better azimuthal coverage, more polarity picks, higher quality FMs'),
        ('Individual mechanisms', 'Each event has its own FM vs one mechanism shared by a cluster of events'),
        ('DL polarity picking', 'PolarPicker model vs CC-based SVD — enables scalable near-real-time processing'),
        ('New eruption cycle', '2023 eruption captured — independent comparison with the 2015 eruption cycle'),
        ('Post-2021 baseline', 'Covers the inflation period leading up to the next eruption'),
    ]
    y_adv = 0.600
    for short, desc in advances:
        ax.text(0.065, y_adv, f'▸  {short}:', fontsize=8.5, fontweight='bold',
                color='#1a3a5c', transform=ax.transAxes)
        ax.text(0.265, y_adv, desc, fontsize=8.2, color='#333333',
                transform=ax.transAxes)
        y_adv -= 0.028

    # ── Key questions ──
    ax.add_patch(FancyBboxPatch((0.04, 0.437), 0.92, 0.022,
                 boxstyle='square,pad=0', facecolor='#7b3294',
                 transform=ax.transAxes, clip_on=False))
    ax.text(0.06, 0.448, 'OPEN SCIENTIFIC QUESTIONS  (motivated by Zhang et al.)',
            fontsize=9, fontweight='bold', color='white', transform=ax.transAxes)

    questions = [
        ('Q1', 'Does the 2023 eruption show the same normal→reverse→normal FM switch as 2015?'),
        ('Q2', 'Does Shallow East again show pre-eruptive reverse faulting before 2023? (eruption precursor signal)'),
        ('Q3', 'How has the caldera matured in the 2022–2024 period? Is the fault system more/less coupled?'),
        ('Q4', 'Did the International District temperature rise continue post-2021? Does seismicity track it?'),
        ('Q5', 'With 21 stations, can we resolve the west/east wall asymmetry (strongly vs weakly locked) more clearly?'),
        ('Q6', 'Does the stress inversion show the same Mogi-like pattern for 2022–2024 inflation?'),
        ('Q7', 'How do individual FMs differ from composite FMs? Does compositing mask spatial heterogeneity?'),
        ('Q8', 'Does the median Kagan angle (48° between 1D/3D in 2015 paper) improve with 21 stations? (it does: 24°)'),
    ]
    y_q = 0.413
    QCOLORS = ['#2166ac','#1a9850','#d73027','#762a83','#b35806','#2166ac','#1a9850','#d73027']
    for (qid, qtxt), qc in zip(questions, QCOLORS):
        ax.text(0.065, y_q, qid + ':', fontsize=8.5, fontweight='bold',
                color=qc, transform=ax.transAxes)
        ax.text(0.105, y_q, qtxt, fontsize=8.2, color='#333333',
                transform=ax.transAxes)
        y_q -= 0.026

    # ── Notable direct answers already known ──
    ax.add_patch(FancyBboxPatch((0.04, 0.050), 0.92, 0.192,
                 boxstyle='round,pad=0.005', linewidth=0.8,
                 edgecolor='#4caf50', facecolor='#f0fff0',
                 transform=ax.transAxes, clip_on=False))
    ax.text(0.065, 0.228, 'ALREADY ANSWERED BY THIS DATASET:', fontsize=8.5,
            fontweight='bold', color='#1a6b1a', transform=ax.transAxes)
    answered = [
        '▸  Q8 resolved: Median Kagan angle drops from 48° (7 stations) to 24° (21 stations) — 3× improvement in 1D/3D agreement.',
        '▸  2023 eruption visible in Y1: Y1 A+B has more strike-slip (25%) vs Y2 (15%), consistent with pre/syn-eruptive stress changes.',
        '▸  East caldera dominates again (63% of events) — consistent with east wall being more active and weakly coupled.',
        '▸  Shallow seismogenic zone (IQR 1.27–1.62 km) confirmed with higher precision due to 3× more stations.',
        '▸  Normal faulting dominates (60%) consistent with post-eruption inflation phase, same as post-2015 in Zhang et al.',
    ]
    y_ans = 0.207
    for line in answered:
        ax.text(0.068, y_ans, line, fontsize=8, color='#1a6b1a', transform=ax.transAxes)
        y_ans -= 0.028

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    # ══════════════════════════════════════════════════════════════
    # PAGE 11 — Research Plan
    # ══════════════════════════════════════════════════════════════
    fig = plt.figure(figsize=(8.5, 11), facecolor='white')
    add_page_header(fig, 'Research Context & Exploration Plan  (Part 2 of 2)',
                    'Proposed Analyses — Priority, Method & Expected Outcome')
    add_footer(fig, 11)

    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')

    ax.add_patch(FancyBboxPatch((0.04, 0.905), 0.92, 0.024,
                 boxstyle='square,pad=0', facecolor=TITLE_COLOR,
                 transform=ax.transAxes, clip_on=False))
    ax.text(0.06, 0.916, 'PROPOSED RESEARCH PLAN  (8 analyses, prioritized)',
            fontsize=9, fontweight='bold', color='white', transform=ax.transAxes)

    plans = [
        ('#1', '#d73027', 'HIGH',
         '2023 Eruption Cycle FM Analysis  (direct sequel to Zhang et al.)',
         ['Map FM type before/during/after 2023 eruption. Test whether normal→reverse→normal',
          'switch repeats. Compare footprint, magnitude, and fault-type fractions with 2015.',
          'Method: Define eruption window from BOTPT pressure records; subset FMs; reproduce',
          'Zhang et al. Fig.4-style maps for all 3 periods. Critical test of cycle reproducibility.']),
        ('#2', '#d73027', 'HIGH',
         'Stress Inversion for 2022–2024  (6 spatial subregions, 3 time periods)',
         ['Apply Michael (1984) stress inversion to A+B FMs in same 6 subregions as Zhang et al.',
          'Recover σ₁/σ₂/σ₃ orientations and stress ratio R per period. Test whether Shallow East',
          'again shows anomalous R and ~90° rotation. Quantify if the Mogi-like signal persists.',
          'This is the most direct scientific extension of the previous paper.']),
        ('#3', '#ff9800', 'MEDIUM',
         'Pre-eruptive FM Precursor Analysis (Shallow East monitoring)',
         ['Zhang et al. found Shallow East shows pre-eruptive reverse faulting increase. Test',
          'whether the same signal appeared before the 2023 eruption using the 21-station dataset.',
          'Plot FM type fraction vs time (weekly bins) for Shallow East region. If confirmed,',
          'this supports Shallow East as a real-time eruption monitoring target.']),
        ('#4', '#ff9800', 'MEDIUM',
         'International District: Post-2021 Seismicity–Temperature Coupling',
         ['Zhang et al. showed seismicity rate and vent temperature co-increased from late 2017',
          'indicating magma recharge. Extend to 2022–2024 using OOI MJ03C temperature data.',
          'Does the temperature continue to rise? Does FM type shift with temperature?',
          'Quantify the seismicity–temperature correlation with the new individual FM catalog.']),
        ('#5', '#ff9800', 'MEDIUM',
         'Individual vs Composite FM Comparison  (methodological advance)',
         ['Compare individual FMs (this study) against composite FMs (Zhang et al.) for the',
          '7 OO-station subset where both methods apply. Compute Kagan angles per event.',
          'Test whether compositing artificially smooths spatial heterogeneity or whether',
          'individual mechanisms are noisier without adding information.']),
        ('#6', '#2196f3', 'EXPLORATORY',
         'P-axis Stress Map  (spatial stress field at Axial Seamount)',
         ['Plot P-axis (compressional axis) orientations spatially for all A+B FMs.',
          'Look for systematic rotation near caldera walls, rift zones, International District.',
          'Use bootstrap sampling to estimate uncertainty. Compare with Mogi-source prediction.']),
        ('#7', '#2196f3', 'EXPLORATORY',
         'West vs East Wall Coupling Contrast with 21 Stations',
         ['Reproduce Zhang et al. fault surface fitting for west and east walls using the denser',
          'FM catalog. Measure median misfit distance and FM alignment angle. More individual FMs',
          'should improve resolution of the curved west-wall geometry and aseismic gaps.']),
        ('#8', '#2196f3', 'EXPLORATORY',
         'DL Polarity Method Validation Against CC-SVD',
         ['For Y1 events that overlap with the 7-station CC-SVD pipeline, compare per-pick',
          'polarity agreement. Quantify DL accuracy vs CC-SVD as a function of SNR and distance.',
          'Demonstrates the DL approach is production-ready for near-real-time FM catalogs.']),
    ]

    PRIORITY_BG = {'HIGH': '#ffeaea', 'MEDIUM': '#fff8e1', 'EXPLORATORY': '#e8f4fd'}
    y0 = 0.887
    for pid, color, priority, title, body_lines in plans:
        h = 0.093
        ax.add_patch(FancyBboxPatch((0.04, y0 - h), 0.92, h - 0.003,
                     boxstyle='round,pad=0.005', linewidth=0.8,
                     edgecolor=color, facecolor=PRIORITY_BG[priority],
                     transform=ax.transAxes, clip_on=False))
        # Badge
        ax.add_patch(FancyBboxPatch((0.045, y0 - 0.022), 0.028, 0.018,
                     boxstyle='round,pad=0.002', facecolor=color,
                     transform=ax.transAxes, clip_on=False))
        ax.text(0.059, y0 - 0.013, pid, ha='center', va='center',
                fontsize=7.5, fontweight='bold', color='white',
                transform=ax.transAxes)
        # Priority tag
        ax.add_patch(FancyBboxPatch((0.077, y0 - 0.022), 0.072, 0.018,
                     boxstyle='round,pad=0.002', facecolor=color, alpha=0.3,
                     transform=ax.transAxes, clip_on=False))
        ax.text(0.113, y0 - 0.013, priority, ha='center', va='center',
                fontsize=6.5, fontweight='bold', color=color,
                transform=ax.transAxes)
        ax.text(0.155, y0 - 0.013, title, fontsize=8.5, fontweight='bold',
                color='#1a1a1a', transform=ax.transAxes)
        for i, line in enumerate(body_lines):
            ax.text(0.055, y0 - 0.030 - i * 0.016, line,
                    fontsize=7.3, color='#444444', transform=ax.transAxes)
        y0 -= h + 0.004

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    # ══════════════════════════════════════════════════════════════
    # PAGE 12 — FM Result Analysis: Overview + Depth
    # ══════════════════════════════════════════════════════════════
    fig = plt.figure(figsize=(8.5, 11), facecolor='white')
    add_page_header(fig, 'FM Result Analysis  (Part 1 of 3)',
                    'Fault-Type Distribution, Depth Structure & Rake Statistics')
    add_footer(fig, 12)

    # ── re-load data for analysis pages ──
    from matplotlib.gridspec import GridSpec

    sk3d_all = pd.concat([
        pd.read_csv(SKHASH_DIR / 'examples/hash3/OUT_Y1_ge8/out.txt'),
        pd.read_csv(SKHASH_DIR / 'examples/hash3/OUT_Y2_ge8/out.txt')
    ], ignore_index=True)

    def classify_fault(rake):
        r = rake % 360
        if r > 180: r -= 360
        if -120 <= r <= -60:  return 'N'
        if   60 <= r <= 120:  return 'R'
        if abs(r) <= 30 or abs(r) >= 150: return 'S'
        return 'U'

    sk3d_all['fault_type'] = [
        classify_fault_ptb(s, d, r)
        for s, d, r in zip(sk3d_all['strike'], sk3d_all['dip'], sk3d_all['rake'])
    ]
    AB3d = sk3d_all[sk3d_all['quality'].isin(['A','B'])].copy()
    AB3d = AB3d[(AB3d['origin_lat'].between(45.90, 46.01)) &
                (AB3d['origin_lon'].between(-130.05, -129.97))].copy()
    AB3d['rake_adj'] = AB3d['rake'].apply(lambda r: r - 360 if r > 180 else r)

    FAULT_COLORS = {'N': '#2166ac', 'R': '#d73027', 'S': '#1a9850', 'U': '#878787'}
    FAULT_LABELS = {'N': 'Normal', 'R': 'Reverse', 'S': 'Strike-slip', 'U': 'Oblique'}

    gs = GridSpec(2, 2, figure=fig, top=0.88, bottom=0.07,
                  hspace=0.45, wspace=0.35, left=0.10, right=0.93)

    # Top-left: fault type pie
    ax1 = fig.add_subplot(gs[0, 0])
    ft_counts = AB3d['fault_type'].value_counts().reindex(['N', 'U', 'S', 'R'], fill_value=0)
    wedges, texts, autotexts = ax1.pie(
        ft_counts.values,
        labels=[FAULT_LABELS[k] for k in ft_counts.index],
        colors=[FAULT_COLORS[k] for k in ft_counts.index],
        autopct='%1.1f%%', startangle=90,
        textprops={'fontsize': 8},
        wedgeprops={'linewidth': 0.5, 'edgecolor': 'white'})
    for at in autotexts:
        at.set_fontsize(7.5)
    ax1.set_title(f'Fault Type  (n={len(AB3d):,})', fontsize=9, fontweight='bold',
                  color=TITLE_COLOR, pad=6)

    # Top-right: depth histogram
    ax2 = fig.add_subplot(gs[0, 1])
    depth_bins = np.arange(0, 4.25, 0.25)
    for ft in ['N', 'U', 'S', 'R']:
        sub = AB3d[AB3d['fault_type'] == ft]['origin_depth_km']
        ax2.hist(sub, bins=depth_bins, orientation='horizontal',
                 color=FAULT_COLORS[ft], alpha=0.75, label=FAULT_LABELS[ft],
                 edgecolor='none', stacked=False)
    ax2.invert_yaxis()
    ax2.set_ylabel('Depth (km)', fontsize=9)
    ax2.set_xlabel('Count', fontsize=9)
    ax2.set_ylim(4.0, 0)
    ax2.axhline(1.0, color='gray', ls=':', lw=0.8)
    ax2.axhline(2.0, color='gray', ls=':', lw=0.8)
    ax2.set_title('Depth Distribution by Fault Type', fontsize=9, fontweight='bold',
                  color=TITLE_COLOR, pad=6)
    ax2.legend(fontsize=7, framealpha=0.9)
    ax2.tick_params(labelsize=8)
    med_d = AB3d['origin_depth_km'].median()
    ax2.axhline(med_d, color='black', ls='--', lw=1.2,
                label=f'Median {med_d:.2f} km')

    # Bottom-left: rake histogram
    ax3 = fig.add_subplot(gs[1, 0])
    rake_bins = np.arange(-180, 181, 10)
    ax3.hist(AB3d['rake_adj'], bins=rake_bins,
             color='#4e79a7', edgecolor='none', alpha=0.85)
    ax3.axvline(-90, color=FAULT_COLORS['N'], lw=1.2, ls='--', label='Pure Normal (−90°)')
    ax3.axvline(0,   color=FAULT_COLORS['S'], lw=1.2, ls='--', label='Strike-slip (0°)')
    ax3.axvline(90,  color=FAULT_COLORS['R'], lw=1.2, ls='--', label='Pure Reverse (+90°)')
    ax3.set_xlabel('Rake (°)', fontsize=9)
    ax3.set_ylabel('Count', fontsize=9)
    ax3.set_title('Rake Distribution', fontsize=9, fontweight='bold', color=TITLE_COLOR, pad=6)
    ax3.legend(fontsize=7, framealpha=0.9)
    ax3.tick_params(labelsize=8)
    ax3.set_xlim(-180, 180)
    ax3.grid(True, alpha=0.3, linewidth=0.5)

    # Bottom-right: depth vs rake scatter
    ax4 = fig.add_subplot(gs[1, 1])
    for ft in ['N', 'U', 'S']:
        sub = AB3d[AB3d['fault_type'] == ft]
        ax4.scatter(sub['rake_adj'], sub['origin_depth_km'],
                    c=FAULT_COLORS[ft], alpha=0.25, s=3, label=FAULT_LABELS[ft])
    ax4.set_xlabel('Rake (°)', fontsize=9)
    ax4.set_ylabel('Depth (km)', fontsize=9)
    ax4.invert_yaxis()
    ax4.set_ylim(3.5, 0)
    ax4.axvline(-90, color=FAULT_COLORS['N'], lw=0.8, ls='--', alpha=0.5)
    ax4.set_title('Rake vs Depth', fontsize=9, fontweight='bold', color=TITLE_COLOR, pad=6)
    ax4.legend(fontsize=7, framealpha=0.9, markerscale=3)
    ax4.tick_params(labelsize=8)
    ax4.grid(True, alpha=0.2, linewidth=0.5)
    ax4.set_xlim(-180, 180)

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    # ══════════════════════════════════════════════════════════════
    # PAGE 11 — FM Analysis: Strike Rose + Spatial Breakdown
    # ══════════════════════════════════════════════════════════════
    fig = plt.figure(figsize=(8.5, 11), facecolor='white')
    add_page_header(fig, 'FM Result Analysis  (Part 2 of 3)',
                    'Strike Orientation & Spatial / Temporal Breakdown')
    add_footer(fig, 13)

    gs2 = GridSpec(2, 2, figure=fig, top=0.88, bottom=0.07,
                   hspace=0.50, wspace=0.38, left=0.10, right=0.93)

    # Top-left: rose diagram of strikes (normal faults only)
    ax_r = fig.add_subplot(gs2[0, 0], projection='polar')
    nf_strikes = AB3d[AB3d['fault_type'] == 'N']['strike'].values
    # Convert to bidirectional (0-180 + mirror)
    both_dirs = np.concatenate([nf_strikes % 360, (nf_strikes + 180) % 360])
    rbins = np.linspace(0, 2 * np.pi, 37)
    counts_r, _ = np.histogram(np.radians(both_dirs), bins=rbins)
    theta_c = (rbins[:-1] + rbins[1:]) / 2
    width = rbins[1] - rbins[0]
    bars_r = ax_r.bar(theta_c, counts_r, width=width, bottom=0,
                      color='#2166ac', alpha=0.75, edgecolor='white', linewidth=0.3)
    ax_r.set_theta_zero_location('N')
    ax_r.set_theta_direction(-1)
    ax_r.set_title('Normal Fault Strikes\n(bidirectional rose)', fontsize=8.5,
                   fontweight='bold', color=TITLE_COLOR, pad=15)
    ax_r.tick_params(labelsize=7)
    ax_r.set_xticklabels(['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'], fontsize=7)
    # Add ridge axis reference
    for az in [10, 190]:  # Juan de Fuca Ridge axis ~N10E
        ax_r.axvline(np.radians(az), color='red', lw=1.2, ls='--', alpha=0.7)
    ax_r.text(np.radians(15), ax_r.get_ylim()[1]*0.85, 'Ridge\naxis', fontsize=6.5,
              color='red', ha='left')

    # Top-right: strike rose for strike-slip
    ax_r2 = fig.add_subplot(gs2[0, 1], projection='polar')
    sf_strikes = AB3d[AB3d['fault_type'] == 'S']['strike'].values
    both_ss = np.concatenate([sf_strikes % 360, (sf_strikes + 180) % 360])
    counts_ss, _ = np.histogram(np.radians(both_ss), bins=rbins)
    ax_r2.bar(theta_c, counts_ss, width=width, bottom=0,
              color='#1a9850', alpha=0.75, edgecolor='white', linewidth=0.3)
    ax_r2.set_theta_zero_location('N')
    ax_r2.set_theta_direction(-1)
    ax_r2.set_title('Strike-slip Fault Strikes\n(bidirectional rose)', fontsize=8.5,
                    fontweight='bold', color=TITLE_COLOR, pad=15)
    ax_r2.tick_params(labelsize=7)
    ax_r2.set_xticklabels(['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'], fontsize=7)

    # Bottom-left: spatial region breakdown bar
    ax_sp = fig.add_subplot(gs2[1, 0])
    regions = {
        'West\ncaldera': AB3d[(AB3d['origin_lon'] < -130.008) & (AB3d['origin_lat'].between(45.93, 45.96))],
        'East\ncaldera': AB3d[(AB3d['origin_lon'] > -130.002) & (AB3d['origin_lat'].between(45.93, 45.97))],
        'North': AB3d[AB3d['origin_lat'] > 45.97],
        'South\n/ID':   AB3d[AB3d['origin_lat'] < 45.93],
    }
    x_pos = np.arange(len(regions))
    bottom = np.zeros(len(regions))
    for ft in ['N', 'U', 'S', 'R']:
        vals = [len(sub[sub['fault_type'] == ft]) / max(len(sub), 1) * 100
                for sub in regions.values()]
        ax_sp.bar(x_pos, vals, bottom=bottom, color=FAULT_COLORS[ft],
                  alpha=0.85, label=FAULT_LABELS[ft])
        bottom += np.array(vals)
    ax_sp.set_xticks(x_pos)
    ax_sp.set_xticklabels(list(regions.keys()), fontsize=8)
    ax_sp.set_ylabel('Percentage (%)', fontsize=9)
    ax_sp.set_title('Fault Type by Spatial Region', fontsize=9, fontweight='bold',
                    color=TITLE_COLOR, pad=6)
    ax_sp.legend(fontsize=7.5, framealpha=0.9, loc='upper right')
    ax_sp.set_ylim(0, 105)
    ax_sp.tick_params(labelsize=8)
    ax_sp.grid(True, axis='y', alpha=0.3, linewidth=0.5)
    for i, (rname, sub) in enumerate(regions.items()):
        ax_sp.text(i, 102, f'n={len(sub)}', ha='center', fontsize=7, color='#333333')

    # Bottom-right: Y1 vs Y2 comparison
    ax_yr = fig.add_subplot(gs2[1, 1])
    sk_y1_a = pd.read_csv(SKHASH_DIR / 'examples/hash3/OUT_Y1_ge8/out.txt')
    sk_y2_a = pd.read_csv(SKHASH_DIR / 'examples/hash3/OUT_Y2_ge8/out.txt')
    for df in [sk_y1_a, sk_y2_a]:
        df['fault_type'] = [
            classify_fault_ptb(s, d, r)
            for s, d, r in zip(df['strike'], df['dip'], df['rake'])
        ]
    y1_ab = sk_y1_a[sk_y1_a['quality'].isin(['A','B'])]
    y2_ab = sk_y2_a[sk_y2_a['quality'].isin(['A','B'])]
    x_yr = np.arange(2)
    bot_yr = np.zeros(2)
    for ft in ['N', 'U', 'S', 'R']:
        vals_yr = [
            y1_ab['fault_type'].eq(ft).sum() / max(len(y1_ab), 1) * 100,
            y2_ab['fault_type'].eq(ft).sum() / max(len(y2_ab), 1) * 100,
        ]
        ax_yr.bar(x_yr, vals_yr, bottom=bot_yr, color=FAULT_COLORS[ft],
                  alpha=0.85, label=FAULT_LABELS[ft])
        bot_yr += np.array(vals_yr)
    ax_yr.set_xticks([0, 1])
    ax_yr.set_xticklabels([
        f'Y1  2022–2023\n(n={len(y1_ab):,})',
        f'Y2  2023–2024\n(n={len(y2_ab):,})'], fontsize=8)
    ax_yr.set_ylabel('Percentage (%)', fontsize=9)
    ax_yr.set_title('Fault Type: Y1 vs Y2', fontsize=9, fontweight='bold',
                    color=TITLE_COLOR, pad=6)
    ax_yr.legend(fontsize=7.5, framealpha=0.9)
    ax_yr.set_ylim(0, 105)
    ax_yr.tick_params(labelsize=8)
    ax_yr.grid(True, axis='y', alpha=0.3, linewidth=0.5)

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    # ══════════════════════════════════════════════════════════════
    # PAGE 12 — FM Analysis: Interpretation & Suggestions
    # ══════════════════════════════════════════════════════════════
    fig = plt.figure(figsize=(8.5, 11), facecolor='white')
    add_page_header(fig, 'FM Result Analysis  (Part 3 of 3)',
                    'Geologic Interpretation & Suggested Follow-up')
    add_footer(fig, 14)

    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')

    def analysis_block(ax, y, num, color, title, lines):
        ax.add_patch(FancyBboxPatch((0.05, y - 0.005), 0.03, 0.028,
                     boxstyle='round,pad=0.005', facecolor=color,
                     transform=ax.transAxes, clip_on=False))
        ax.text(0.065, y + 0.0085, str(num), ha='center', va='center',
                fontsize=10, fontweight='bold', color='white',
                transform=ax.transAxes)
        ax.text(0.10, y + 0.008, title, fontsize=9.5, fontweight='bold',
                color=color, transform=ax.transAxes)
        for i, line in enumerate(lines):
            ax.text(0.10, y - 0.013 - i * 0.019, line, fontsize=8.2,
                    color='#333333', transform=ax.transAxes)

    # ── Interpretation ──
    ax.add_patch(FancyBboxPatch((0.04, 0.875), 0.92, 0.024,
                 boxstyle='square,pad=0', facecolor=TITLE_COLOR,
                 transform=ax.transAxes, clip_on=False))
    ax.text(0.06, 0.886, 'GEOLOGIC INTERPRETATION', fontsize=9.5,
            fontweight='bold', color='white', transform=ax.transAxes)

    analysis_block(ax, 0.840, '1', '#2166ac',
        'Extensional tectonics dominate (~83% with normal/oblique-normal component)',
        ['Normal faults (60%) + oblique-normal (23%) confirm a mid-ocean ridge extensional setting.',
         'Pure normal rake (−90°) is the single dominant mode. Reverse faults are negligible (1%),',
         'consistent with the divergent Juan de Fuca Ridge spreading environment.'])

    analysis_block(ax, 0.745, '2', '#4393c3',
        'Extremely shallow seismogenic zone: IQR 1.27–1.62 km, 91% between 1–2 km',
        ['Seismicity is confined to a thin brittle layer just above the Axial magma chamber.',
         'The near-absence of events >2.5 km suggests the melt lens suppresses deep brittle failure.',
         'This is one of the thinnest seismogenic zones documented at any submarine volcano.'])

    analysis_block(ax, 0.650, '3', '#1a9850',
        'Normal faults strike broadly E-W to ENE-WSW — oblique to the ridge axis',
        ['The ridge axis (Juan de Fuca) trends ~N10°E. Classic ridge-parallel normal faults',
         'would strike N-S. Instead, the median strike is ~76°–82° (E-W), suggesting these',
         'faults accommodate caldera-related collapse or N-S extension within the magma system,',
         'not just plate spreading. West caldera shows N-S strike (174°) more consistent with',
         'the ridge; East caldera (114°) is more oblique, possibly the eastern rift zone.'])

    analysis_block(ax, 0.530, '4', '#1a9850',
        'Strike-slip faults (16%) are predominantly left-lateral on E-W planes',
        ['62% left-lateral / 34% right-lateral on ~E-W striking planes (median strike 96°).',
         'Left-lateral E-W faults could represent: (a) bookshelf faulting in the N-S rift zones,',
         '(b) conjugate shear on caldera ring structures, or (c) oblique interaction between',
         'the spreading direction and the N-S elongated caldera geometry.'])

    analysis_block(ax, 0.430, '5', '#762a83',
        'Y2 (2023–2024) shows stronger normal-fault dominance than Y1 (60% vs 52%)',
        ['Y2 has 10× more events (14,036 vs 1,473 total FMs), consistent with an eruption swarm.',
         'The shift toward more normal faulting in Y2 suggests dike intrusion and rift opening',
         'dominated the 2023 eruption period. Y1\'s higher strike-slip fraction may reflect',
         'pre-eruptive stress reorganization or caldera inflation-related shear.'])

    analysis_block(ax, 0.325, '6', '#b35806',
        'East caldera is the dominant seismic source zone (63% of A+B events)',
        ['1,641 A+B events in the East caldera vs 449 in the West. The eastern rift zone',
         'appears structurally more active, with a distinct E-W to ENE-WSW strike orientation',
         'compared to the more ridge-parallel N-S faulting in the West caldera.'])

    # ── Suggestions ──
    ax.add_patch(FancyBboxPatch((0.04, 0.282), 0.92, 0.024,
                 boxstyle='square,pad=0', facecolor='#7b3294',
                 transform=ax.transAxes, clip_on=False))
    ax.text(0.06, 0.292, 'SUGGESTED FOLLOW-UP ANALYSES', fontsize=9.5,
            fontweight='bold', color='white', transform=ax.transAxes)

    suggestions = [
        ('a', '#7b3294', 'Temporal evolution',
         'Plot FM type and event rate vs time, aligned with the 2023 eruption date to isolate\n'
         '         pre-eruptive, syn-eruptive, and post-eruptive stress states.'),
        ('b', '#7b3294', 'P-axis stress inversion',
         'Run a stress inversion (e.g., Michael 1984) on the A+B FMs to recover the\n'
         '         principal stress axes and the relative stress ratio R = (σ2−σ3)/(σ1−σ3).'),
        ('c', '#7b3294', 'Compare with caldera geometry',
         'Overlay normal fault strikes on the caldera rim map — check whether E-W striking\n'
         '         faults align with the short axis of the caldera or with the 2015/2023 lava flows.'),
        ('d', '#7b3294', 'Depth vs. seismic reflection',
         'Compare the 1–2 km seismicity floor with published seismic reflection profiles\n'
         '         showing the melt lens reflector depth.'),
    ]
    y_s = 0.250
    for let, color, title, body in suggestions:
        ax.text(0.075, y_s, f'{let})', fontsize=8.5, fontweight='bold', color=color,
                transform=ax.transAxes)
        ax.text(0.105, y_s, title + ':', fontsize=8.5, fontweight='bold', color='#222222',
                transform=ax.transAxes)
        ax.text(0.105, y_s - 0.022, body, fontsize=8, color='#444444',
                transform=ax.transAxes)
        y_s -= 0.062

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    # ── Page 15: CC-SVD vs DL polarity coverage ─────────────────────────────
    from collections import defaultdict

    def count_polarities(year):
        out_e = DATA_DIR / 'E_Po_py'
        cc_pol  = defaultdict(int)
        dl_pol  = defaultdict(int)
        dl_hq   = defaultdict(int)
        for path in sorted(out_e.glob(f'E_{year}_*.npy')):
            arr       = np.load(path)
            eids      = arr[:, 0].astype(int)
            svd_final = arr[:, 2].astype(int)
            dl_p      = arr[:, 3].astype(int)
            dl_c      = arr[:, 4]
            for i, eid in enumerate(eids):
                if svd_final[i] != 0: cc_pol[eid] += 1
                if dl_p[i]      != 0:
                    dl_pol[eid] += 1
                    if dl_c[i] >= 0.8: dl_hq[eid] += 1
        return cc_pol, dl_pol, dl_hq

    y1_cc, y1_dl, y1_hq = count_polarities('Y1')
    y2_cc, y2_dl, y2_hq = count_polarities('Y2')

    def ge_thr(d, t): return sum(1 for v in d.values() if v >= t)

    thrs   = [5, 6, 7, 8, 10, 12, 15]
    rows   = []
    for t in thrs:
        rows.append([f'≥{t}',
                     f'{ge_thr(y1_cc,t):,}', f'{ge_thr(y1_dl,t):,}', f'{ge_thr(y1_hq,t):,}',
                     f'{ge_thr(y2_cc,t):,}', f'{ge_thr(y2_dl,t):,}', f'{ge_thr(y2_hq,t):,}'])

    fig = plt.figure(figsize=(14, 10), facecolor='white')
    add_page_header(fig, 'DL vs CC-SVD Polarity Coverage',
                    'Events with ≥N station polarities — Y1 & Y2')
    add_footer(fig, 15)

    # ── Bar charts (top) ──────────────────────────────────────────────────────
    ax1 = fig.add_axes([0.06, 0.56, 0.40, 0.30])
    ax2 = fig.add_axes([0.54, 0.56, 0.40, 0.30])

    bar_w = 0.25
    x = np.arange(len(thrs))
    colors_cc = '#2166ac'; colors_dl = '#d6604d'; colors_hq = '#4dac26'

    for ax, ycc, ydl, yhq, title, ntot in [
        (ax1, y1_cc, y1_dl, y1_hq, 'Y1  (3,105 events)', 3105),
        (ax2, y2_cc, y2_dl, y2_hq, 'Y2  (20,636 events)', 20636),
    ]:
        vals_cc = [ge_thr(ycc, t) for t in thrs]
        vals_dl = [ge_thr(ydl, t) for t in thrs]
        vals_hq = [ge_thr(yhq, t) for t in thrs]
        ax.bar(x - bar_w, vals_cc, bar_w, label='CC-SVD',       color=colors_cc, alpha=0.85)
        ax.bar(x,          vals_dl, bar_w, label='DL (all)',     color=colors_dl, alpha=0.85)
        ax.bar(x + bar_w,  vals_hq, bar_w, label='DL (≥0.8)',   color=colors_hq, alpha=0.85)
        ax.set_xticks(x); ax.set_xticklabels([f'≥{t}' for t in thrs], fontsize=8)
        ax.set_ylabel('# Events', fontsize=8); ax.set_title(title, fontsize=9, fontweight='bold')
        ax.legend(fontsize=7, loc='upper right')
        ax.grid(True, axis='y', alpha=0.3, linewidth=0.5)
        ax.tick_params(labelsize=7)

    # ── Distribution histograms (middle) ─────────────────────────────────────
    ax3 = fig.add_axes([0.06, 0.27, 0.40, 0.22])
    ax4 = fig.add_axes([0.54, 0.27, 0.40, 0.22])

    for ax, ycc, ydl, title in [
        (ax3, y1_cc, y1_dl, 'Y1 — polarity count distribution per event'),
        (ax4, y2_cc, y2_dl, 'Y2 — polarity count distribution per event'),
    ]:
        max_sta = 21
        bins = np.arange(1, max_sta + 2) - 0.5
        cc_vals = [v for v in ycc.values()]
        dl_vals = [v for v in ydl.values()]
        ax.hist(cc_vals, bins=bins, color=colors_cc, alpha=0.7, label='CC-SVD', density=False)
        ax.hist(dl_vals, bins=bins, color=colors_dl, alpha=0.5, label='DL (all)', density=False)
        ax.axvline(8, color='black', lw=1.2, ls='--', label='≥8 threshold')
        ax.set_xlabel('# Stations with polarity', fontsize=8)
        ax.set_ylabel('# Events', fontsize=8)
        ax.set_title(title, fontsize=8.5, fontweight='bold')
        ax.legend(fontsize=7)
        ax.tick_params(labelsize=7)
        ax.grid(True, axis='y', alpha=0.3, linewidth=0.5)

    # ── Summary table (bottom) ────────────────────────────────────────────────
    ax5 = fig.add_axes([0.04, 0.02, 0.92, 0.20])
    ax5.axis('off')
    col_labels = ['Threshold',
                  'Y1 CC-SVD', 'Y1 DL (all)', 'Y1 DL (≥0.8)',
                  'Y2 CC-SVD', 'Y2 DL (all)', 'Y2 DL (≥0.8)']
    tbl = ax5.table(cellText=rows, colLabels=col_labels,
                    cellLoc='center', loc='center')
    tbl.auto_set_font_size(False); tbl.set_fontsize(8.5); tbl.scale(1, 1.8)
    HDR = '#2c3e50'
    for j in range(len(col_labels)):
        tbl[(0, j)].set_facecolor(HDR)
        tbl[(0, j)].set_text_props(color='white', fontweight='bold')
    for i in range(1, len(rows) + 1):
        bg = '#f2f2f2' if i % 2 == 0 else 'white'
        for j in range(len(col_labels)):
            tbl[(i, j)].set_facecolor(bg)
        # Highlight ≥8 row
        if rows[i-1][0] == '≥8':
            for j in range(len(col_labels)):
                tbl[(i, j)].set_facecolor('#fff9c4')
                tbl[(i, j)].set_text_props(fontweight='bold')

    ax5.set_title('Event counts at each polarity threshold  (≥8 highlighted — HASH/SKHASH minimum)',
                  fontsize=8.5, pad=4)

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    # ── Page 16: HASH single-event FM — CC-SVD vs DL Kagan angles ───────────────
    from collections import defaultdict

    QUAL_COLORS = {'A': '#d62728', 'B': '#ff7f0e', 'C': '#2ca02c', 'D': '#9467bd'}
    QUAL_ORDER  = list('ABCD')
    bins_kag    = np.arange(0, 91, 5)

    def load_hash_out(path):
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

    def _sdr2ptb_k(s, d, r):
        s, d, r = np.radians(s), np.radians(d), np.radians(r)
        n  = np.array([-np.sin(d)*np.sin(s),  np.sin(d)*np.cos(s), -np.cos(d)])
        sl = np.array([ np.cos(r)*np.cos(s)+np.sin(r)*np.cos(d)*np.sin(s),
                         np.cos(r)*np.sin(s)-np.sin(r)*np.cos(d)*np.cos(s),
                        -np.sin(r)*np.sin(d)])
        T = (n + sl) / np.sqrt(2); P = (n - sl) / np.sqrt(2); B = np.cross(T, P)
        return P, T, B

    def kagan_angle_k(s1, d1, r1, s2, d2, r2):
        P1, T1, B1 = _sdr2ptb_k(s1, d1, r1)
        P2, T2, B2 = _sdr2ptb_k(s2, d2, r2)
        R2 = np.column_stack([P2, T2, B2])
        min_ang = np.inf
        for sp, st, sb in [(1,1,1), (-1,-1,1), (1,-1,-1), (-1,1,-1)]:
            R1v = np.column_stack([sp*P1, st*T1, sb*B1])
            cos_a = np.clip((np.trace(R1v.T @ R2) - 1) / 2, -1, 1)
            min_ang = min(min_ang, np.degrees(np.arccos(cos_a)))
        return min_ang

    kagan_vec_k = np.vectorize(kagan_angle_k)

    # Load HASH results for all 6 tags
    hash_results = {}
    for tag in ['Y1_cc', 'Y1_dl_all', 'Y1_dl_hq', 'Y2_cc', 'Y2_dl_all', 'Y2_dl_hq']:
        out_path = HASH_DIR / f'hashout_{tag}_1.dat'
        eid_map_path = OUT_DIR / f'hash_{tag}_eid_map.csv'
        if not out_path.exists() or not eid_map_path.exists():
            hash_results[tag] = pd.DataFrame()
            continue
        df = load_hash_out(out_path)
        eid_df = pd.read_csv(eid_map_path)
        df = df.merge(eid_df, on='cluster_idx', how='left')
        hash_results[tag] = df

    def merge_kagan(tag_a, tag_b):
        a = hash_results.get(tag_a, pd.DataFrame())
        b = hash_results.get(tag_b, pd.DataFrame())
        if len(a) == 0 or len(b) == 0:
            return pd.DataFrame()
        m = a[['event_id','strike','dip','rake','quality']].merge(
            b[['event_id','strike','dip','rake','quality']].rename(
                columns={'strike':'h_s','dip':'h_d','rake':'h_r','quality':'h_q'}),
            on='event_id', how='inner')
        if len(m) == 0:
            return m
        m['kagan'] = kagan_vec_k(m['strike'], m['dip'], m['rake'],
                                   m['h_s'], m['h_d'], m['h_r'])
        return m

    # Build comparisons
    comps = {}
    for year in ['Y1', 'Y2']:
        comps[f'{year}_cc_vs_dlall'] = merge_kagan(f'{year}_cc', f'{year}_dl_all')
        comps[f'{year}_cc_vs_dlhq']  = merge_kagan(f'{year}_cc', f'{year}_dl_hq')
        comps[f'{year}_dl_vs_dlhq']  = merge_kagan(f'{year}_dl_all', f'{year}_dl_hq')

    comp_labels = {
        'cc_vs_dlall': 'CC-SVD vs DL-all',
        'cc_vs_dlhq':  'CC-SVD vs DL-hq (≥0.8)',
        'dl_vs_dlhq':  'DL-all vs DL-hq (≥0.8)',
    }

    fig = plt.figure(figsize=(17, 22), facecolor='white')
    add_page_header(fig, 'HASH Single-Event FMs: CC-SVD vs DL Kagan Angles',
                    '1D axial velocity model · ≥8 polarities · Y1 (2022–23) and Y2 (2023–24)')

    # Layout: 2 rows (Y1/Y2) × 3 cols (3 comparisons) × 2 sub-rows (hist + CDF) = 12 axes
    outer = gridspec.GridSpec(2, 1, figure=fig, top=0.91, bottom=0.04,
                               hspace=0.30)
    for row_i, year in enumerate(['Y1', 'Y2']):
        inner = gridspec.GridSpecFromSubplotSpec(2, 3, subplot_spec=outer[row_i],
                                                  hspace=0.55, wspace=0.35)
        for col_j, ck in enumerate(['cc_vs_dlall', 'cc_vs_dlhq', 'dl_vs_dlhq']):
            key = f'{year}_{ck}'
            merged = comps.get(key, pd.DataFrame())
            label  = comp_labels[ck]

            ax_hist = fig.add_subplot(inner[0, col_j])
            ax_cdf  = fig.add_subplot(inner[1, col_j])

            if len(merged) == 0:
                ax_hist.text(0.5, 0.5, 'No data', ha='center', va='center',
                             transform=ax_hist.transAxes)
                ax_cdf.axis('off')
                continue

            # Histogram
            data_by_q = [merged[merged['quality'] == q]['kagan'].values
                         for q in QUAL_ORDER]
            ax_hist.hist(data_by_q, bins=bins_kag, stacked=True,
                         color=[QUAL_COLORS[q] for q in QUAL_ORDER],
                         label=[f'{q}(n={len(d):,})' for q, d in
                                zip(QUAL_ORDER, data_by_q)],
                         edgecolor='none', alpha=0.9)
            med = merged['kagan'].median()
            ax_hist.axvline(med, color='k', ls='--', lw=1.2,
                            label=f'Med={med:.1f}°')
            ax_hist.set_title(f'{year}: {label}\n(n={len(merged):,})',
                               fontsize=8.5)
            ax_hist.set_xlabel('Kagan angle (°)', fontsize=8)
            ax_hist.set_ylabel('Count', fontsize=8)
            ax_hist.set_xlim(0, 90)
            ax_hist.legend(fontsize=6.5, framealpha=0.85, ncol=2)
            ax_hist.grid(True, alpha=0.3, lw=0.5)
            ax_hist.tick_params(labelsize=7)

            # CDF
            for q in QUAL_ORDER:
                sub = merged[merged['quality'] == q]['kagan']
                if len(sub) == 0:
                    continue
                sk = np.sort(sub)
                cdf = np.arange(1, len(sk)+1) / len(sk)
                ax_cdf.plot(sk, cdf * 100, color=QUAL_COLORS[q], lw=1.3,
                            label=f'{q} med={sub.median():.0f}°')
            pct30 = (merged['kagan'] <= 30).mean() * 100
            ax_cdf.axvline(30, color='grey', ls=':', lw=1,
                           label=f'≤30°: {pct30:.0f}%')
            ax_cdf.set_xlabel('Kagan angle (°)', fontsize=8)
            ax_cdf.set_ylabel('Cumulative %', fontsize=8)
            ax_cdf.set_title(f'CDF by CC quality', fontsize=8)
            ax_cdf.legend(fontsize=6.5, framealpha=0.85)
            ax_cdf.set_xlim(0, 90)
            ax_cdf.set_ylim(0, 100)
            ax_cdf.grid(True, alpha=0.3, lw=0.5)
            ax_cdf.tick_params(labelsize=7)

    # Summary table at bottom
    sum_rows = []
    for year in ['Y1', 'Y2']:
        for ck, label in comp_labels.items():
            key = f'{year}_{ck}'
            m = comps.get(key, pd.DataFrame())
            if len(m) == 0:
                continue
            pct30 = (m['kagan'] <= 30).mean() * 100
            qmed = {q: m[m['quality']==q]['kagan'].median()
                    for q in QUAL_ORDER if len(m[m['quality']==q]) > 0}
            qmed_str = '  '.join(f'{q}:{v:.0f}°' for q, v in qmed.items())
            sum_rows.append([year, label, f'{len(m):,}',
                             f'{m.kagan.median():.1f}°',
                             f'{pct30:.1f}%', qmed_str])

    ax_tbl = fig.add_axes([0.03, 0.005, 0.94, 0.11])
    ax_tbl.axis('off')
    col_labels_t = ['Year', 'Comparison', 'N matched', 'Median Kagan', '% ≤ 30°',
                    'Median by quality (A/B/C/D)']
    tbl = ax_tbl.table(cellText=sum_rows, colLabels=col_labels_t,
                        loc='center', cellLoc='center')
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7.5)
    tbl.scale(1, 1.4)
    for j in range(len(col_labels_t)):
        tbl[(0, j)].set_facecolor(HEADER_COLOR)
        tbl[(0, j)].set_text_props(color='white', fontweight='bold')
    for i in range(1, len(sum_rows) + 1):
        bg = ROW_ALT if i % 2 == 0 else ROW_WHITE
        for j in range(len(col_labels_t)):
            tbl[(i, j)].set_facecolor(bg)
    ax_tbl.set_title('Summary: HASH single-event FM Kagan angles by comparison type',
                      fontsize=8.5, pad=3)

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    # ── Page 17: All-depth FM comparison (CC-SVD vs DL-hq, 4 panels) ────────────
    from obspy.imaging.beachball import beach as _beach

    def _classify_fault(rake):
        r = rake % 360
        if r > 180: r -= 360
        if -120 <= r <= -60:  return 'N'
        if   60 <= r <=  120: return 'R'
        if abs(r) <= 30 or abs(r) >= 150: return 'S'
        return 'U'

    def _load_hash_full(path):
        rows = []
        with open(path) as fh:
            for line in fh:
                p = line.split()
                if len(p) < 31: continue
                try:
                    rows.append(dict(
                        origin_lat=float(p[10]), origin_lon=float(p[11]),
                        origin_depth_km=float(p[12]),
                        strike=float(p[21]), dip=float(p[22]), rake=float(p[23]),
                        quality=p[28]))
                except (ValueError, IndexError):
                    continue
        df = pd.DataFrame(rows)
        if len(df):
            df['fault_type'] = [
                classify_fault_ptb(s, d, r)
                for s, d, r in zip(df['strike'], df['dip'], df['rake'])
            ]
            df['qrank'] = df['quality'].map({'A':4,'B':3,'C':2,'D':1})
        return df

    STA_OO_BB = {
        'AS1':(45.93356,-129.99920),'AS2':(45.93377,-130.01410),
        'CC1':(45.95468,-130.00890),'EC1':(45.94958,-129.97970),
        'EC2':(45.93967,-129.97380),'EC3':(45.93607,-129.97850),
        'ID1':(45.92573,-129.97800),
    }
    LON_LIM_F = (-130.04, -129.97); LAT_LIM_F = (45.908, 46.001)
    BB_W_F    = 0.0012
    FC = {'N':[0.15,0.25,0.85],'R':[0.85,0.15,0.15],'S':[0.10,0.70,0.20],'U':[0.40,0.40,0.40]}
    CALDERA_F = np.array([
        [-130.004785563058,45.9207755734405],[-130.010476202888,45.9238241104543],
        [-130.018881564079,45.9351908809594],[-130.023946125193,45.9412238501725],
        [-130.028718653506,45.9498812001140],[-130.030451219380,45.9511765797916],
        [-130.030679485650,45.9542732243167],[-130.031733279709,45.9558130656063],
        [-130.031444653500,45.9586760104296],[-130.036188782208,45.9656647517656],
        [-130.036950110789,45.9698291665232],[-130.039953347000,45.9750458167927],
        [-130.038595675479,45.9847117727418],[-130.035927416999,45.9883113986506],
        [-130.018067675296,45.9933582886740],[-130.013629193751,45.9937552841350],
        [-130.010365710979,45.9929499241491],[-130.008647442296,45.9924883829037],
        [-130.007262470669,45.9915471582374],[-130.006042022411,45.9902469280907],
        [-130.005178629490,45.9897778053610],[-130.001868199523,45.9863506519894],
        [-130.001154359192,45.9846883853932],[-130.000949059432,45.9827833001814],
        [-129.999393534330,45.9818434725493],[-129.997797388662,45.9786395525337],
        [-129.995357566829,45.9760388622191],[-129.993956176267,45.9741441737512],
        [-129.993678708114,45.9681875631427],[-129.993140494256,45.9667620754035],
        [-129.992087550788,45.9652218741086],[-129.991186410747,45.9626077204113],
        [-129.989604036931,45.9601186407320],[-129.989238986151,45.9588108137369],
        [-129.989728217453,45.9574955894078],[-129.985484098670,45.9494279735802],
        [-129.984788122490,45.9487188881587],
    ])

    panel_tags_f   = ['Y1_cc', 'Y1_dl_hq', 'Y2_cc', 'Y2_dl_hq']
    panel_labels_f = ['Y1 CC-SVD', 'Y1 DL-hq (≥0.8)', 'Y2 CC-SVD', 'Y2 DL-hq (≥0.8)']

    lon_r_f = abs(LON_LIM_F[1]-LON_LIM_F[0])
    lat_r_f = abs(LAT_LIM_F[1]-LAT_LIM_F[0])
    pw_f = 3.8; ph_f = pw_f * lat_r_f / lon_r_f
    fig, axes = plt.subplots(1, 4, figsize=(pw_f*4+0.6, ph_f+1.4), facecolor='white')
    add_page_header(fig,
                    'FM Comparison: CC-SVD vs DL-hq  |  All Depths  |  Y1 & Y2',
                    'Q = A+B only  ·  Blue=Normal  Red=Reverse  Green=Strike-slip  Gray=Oblique')

    for ax, tag, lbl in zip(axes, panel_tags_f, panel_labels_f):
        p = HASH_DIR / f'hashout_{tag}_1.dat'
        sub = pd.DataFrame()
        if p.exists():
            df_f = _load_hash_full(p)
            sub = df_f[df_f['quality'].isin(['A','B'])].copy()
            sub = sub[(sub['origin_lat'].between(*LAT_LIM_F)) &
                      (sub['origin_lon'].between(*LON_LIM_F))]
            sub = sub.sort_values('qrank')
        ax.set_facecolor([0.85, 0.92, 0.97])
        ax.set_aspect('equal')
        for _, row in sub.iterrows():
            try:
                bb = _beach([row['strike'], row['dip'], row['rake']],
                            xy=(row['origin_lon'], row['origin_lat']),
                            width=BB_W_F, linewidth=0.15,
                            facecolor=FC[row['fault_type']], alpha=0.85)
                bb.set_transform(ax.transData); bb.set_zorder(2)
                ax.add_collection(bb)
            except Exception:
                pass
        ax.plot(CALDERA_F[:,0], CALDERA_F[:,1], '-', color='black', lw=1.2, zorder=7)
        ax.plot([v[1] for v in STA_OO_BB.values()],
                [v[0] for v in STA_OO_BB.values()],
                '^', color='yellow', markeredgecolor='k', markeredgewidth=0.3,
                markersize=4, zorder=8)
        ax.set_xlim(LON_LIM_F); ax.set_ylim(LAT_LIM_F)
        tks = np.arange(-130.04, -129.96, 0.02)
        ax.set_xticks(tks)
        ax.set_xticklabels([f'{t:.2f}' for t in tks], fontsize=5.5, rotation=45)
        ax.tick_params(axis='y', labelsize=5.5)
        ax.grid(True, lw=0.3, color='gray', alpha=0.4)
        ax.set_xlabel('Longitude', fontsize=7)
        ax.set_title(f'{lbl}\nn={len(sub)}', fontsize=8.5, pad=3)
    axes[0].set_ylabel('Latitude', fontsize=7)
    _leg_h = [mpatches.Patch(color=FC[k], label=v) for k,v in
              [('N','Normal'),('R','Reverse'),('S','Strike-slip'),('U','Oblique')]]
    axes[-1].legend(handles=_leg_h, loc='lower right', fontsize=6.5, framealpha=0.9)
    plt.subplots_adjust(top=0.90, bottom=0.10, left=0.05, right=0.99, wspace=0.15)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

    # ── Pages 18–19: Depth-slice beach-ball FM maps (CC-SVD vs DL-hq, Y1+Y2) ────

    LON_LIM_BB  = (-130.04, -129.97)
    LAT_LIM_BB  = (45.908,  46.001)
    BB_WIDTH_PG = 0.0012
    QUAL_FILTER_BB = ['A', 'B']
    QUAL_ORDER_BB  = {'A': 4, 'B': 3, 'C': 2, 'D': 1}
    FAULT_COLOR_BB = {
        'N': [0.15, 0.25, 0.85],
        'R': [0.85, 0.15, 0.15],
        'S': [0.10, 0.70, 0.20],
        'U': [0.40, 0.40, 0.40],
    }
    CALDERA_RIM_BB = np.array([
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
    STA_OO_BB = {
        'AS1': (45.93356, -129.99920), 'AS2': (45.93377, -130.01410),
        'CC1': (45.95468, -130.00890), 'EC1': (45.94958, -129.97970),
        'EC2': (45.93967, -129.97380), 'EC3': (45.93607, -129.97850),
        'ID1': (45.92573, -129.97800),
    }

    def classify_fault_bb(rake):
        r = rake % 360
        if r > 180: r -= 360
        if -120 <= r <= -60:  return 'N'
        if   60 <= r <=  120: return 'R'
        if abs(r) <= 30 or abs(r) >= 150: return 'S'
        return 'U'

    def load_hash_bb(path):
        rows = []
        with open(path) as fh:
            for line in fh:
                p = line.split()
                if len(p) < 31: continue
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
            df['fault_type'] = df['rake'].apply(classify_fault_bb)
            df['qrank']      = df['quality'].map(QUAL_ORDER_BB)
        return df

    bb_datasets = {}
    for tag in ['Y1_cc', 'Y1_dl_hq', 'Y2_cc', 'Y2_dl_hq']:
        p = HASH_DIR / f'hashout_{tag}_1.dat'
        if p.exists():
            df = load_hash_bb(p)
            df = df[df['quality'].isin(QUAL_FILTER_BB)].copy()
            df = df[(df['origin_lat'].between(*LAT_LIM_BB)) &
                    (df['origin_lon'].between(*LON_LIM_BB))]
            bb_datasets[tag] = df
        else:
            bb_datasets[tag] = pd.DataFrame()

    def draw_bb_slice(ax, sub, title):
        ax.set_facecolor([0.85, 0.92, 0.97])
        ax.set_aspect('equal')
        if len(sub):
            sub = sub.sort_values('qrank')
            for _, row in sub.iterrows():
                try:
                    bb = _beach([row['strike'], row['dip'], row['rake']],
                                xy=(row['origin_lon'], row['origin_lat']),
                                width=BB_WIDTH_PG, linewidth=0.15,
                                facecolor=FAULT_COLOR_BB[row['fault_type']],
                                alpha=0.85)
                    bb.set_transform(ax.transData)
                    bb.set_zorder(2)
                    ax.add_collection(bb)
                except Exception:
                    pass
        ax.plot(CALDERA_RIM_BB[:, 0], CALDERA_RIM_BB[:, 1], '-',
                color='black', linewidth=1.0, zorder=7)
        ax.plot([v[1] for v in STA_OO_BB.values()],
                [v[0] for v in STA_OO_BB.values()],
                '^', color='yellow', markeredgecolor='k',
                markeredgewidth=0.3, markersize=4, zorder=8)
        ax.set_xlim(LON_LIM_BB); ax.set_ylim(LAT_LIM_BB)
        ticks_x = np.arange(-130.04, -129.96, 0.02)
        ax.set_xticks(ticks_x)
        ax.set_xticklabels([f'{t:.2f}' for t in ticks_x], fontsize=5.5, rotation=45)
        ax.tick_params(axis='y', labelsize=5.5)
        ax.grid(True, linewidth=0.3, color='gray', alpha=0.4)
        ax.set_title(title, fontsize=7.5, pad=2)

    legend_handles_bb = [
        mpatches.Patch(color=FAULT_COLOR_BB['N'], label='Normal'),
        mpatches.Patch(color=FAULT_COLOR_BB['R'], label='Reverse'),
        mpatches.Patch(color=FAULT_COLOR_BB['S'], label='Strike-slip'),
        mpatches.Patch(color=FAULT_COLOR_BB['U'], label='Oblique'),
    ]

    depth_bins_bb = [(d, d + 0.25) for d in np.arange(0.5, 2.5, 0.25)]  # 8 bins
    col_tags    = ['Y1_cc', 'Y1_dl_hq', 'Y2_cc', 'Y2_dl_hq']
    col_titles  = ['Y1 CC-SVD', 'Y1 DL-hq (≥0.8)', 'Y2 CC-SVD', 'Y2 DL-hq (≥0.8)']

    # Split 8 depth bins into two pages of 4 rows each
    lon_range_bb = abs(LON_LIM_BB[1] - LON_LIM_BB[0])  # 0.07
    lat_range_bb = abs(LAT_LIM_BB[1] - LAT_LIM_BB[0])  # 0.093
    panel_w_bb = 3.8
    panel_h_bb = panel_w_bb * lat_range_bb / lon_range_bb  # ~5.05 in

    for page_i, slice_group in enumerate([depth_bins_bb[:4], depth_bins_bb[4:]]):
        n_rows = len(slice_group)
        fig_w = panel_w_bb * 4 + 0.6
        fig_h = panel_h_bb * n_rows + 1.2  # extra for header
        fig, axes = plt.subplots(n_rows, 4,
                                  figsize=(fig_w, fig_h),
                                  facecolor='white',
                                  gridspec_kw={'hspace': 0.45, 'wspace': 0.18})
        depth_start = slice_group[0][0]
        depth_end   = slice_group[-1][1]
        add_page_header(fig,
                        f'Depth-Slice Beach-Ball Maps: CC-SVD vs DL-hq  '
                        f'({depth_start:.2f}–{depth_end:.2f} km)',
                        'Q = A+B only  ·  Blue=Normal  Red=Reverse  Green=Strike-slip  '
                        'Gray=Oblique  ·  Y1 2022–23, Y2 2023–24')

        for row_i, (d0, d1) in enumerate(slice_group):
            for col_j, tag in enumerate(col_tags):
                df_full = bb_datasets.get(tag, pd.DataFrame())
                if len(df_full):
                    sub = df_full[(df_full['origin_depth_km'] >= d0) &
                                  (df_full['origin_depth_km'] <  d1)]
                else:
                    sub = pd.DataFrame()
                title = f'{col_titles[col_j]}\n{d0:.2f}–{d1:.2f} km  (n={len(sub)})'
                draw_bb_slice(axes[row_i, col_j], sub, title)
                if col_j == 0:
                    axes[row_i, col_j].set_ylabel('Latitude', fontsize=6)
                if row_i == n_rows - 1:
                    axes[row_i, col_j].set_xlabel('Longitude', fontsize=6)

        axes[-1, -1].legend(handles=legend_handles_bb, loc='lower right',
                             fontsize=6, framealpha=0.9)
        plt.subplots_adjust(top=0.93, bottom=0.04, left=0.05, right=0.99)
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

    # ── Page 20: SKHASH FM comparison (CC-SVD vs DL-hq, Y1+Y2, all depths) ────────

    SKHASH_DIR_PG = Path('/Users/mczhang/Documents/GitHub/SKHASH/SKHASH7')
    HASH3_DIR_PG  = SKHASH_DIR_PG / 'examples' / 'hash3'
    QUAL_ORDER_SK = {'A': 4, 'B': 3, 'C': 2, 'D': 1}

    def load_skhash_pg(path):
        if not path.exists():
            return pd.DataFrame()
        df = pd.read_csv(path)
        df.columns = df.columns.str.strip()
        df['fault_type'] = df['rake'].apply(classify_fault_bb)
        df['qrank']      = df['quality'].map(QUAL_ORDER_SK)
        return df

    sk_tags = {
        'Y1 CC-SVD (SKHASH)':      HASH3_DIR_PG / 'OUT_Y1_cc_sk'    / 'out.txt',
        'Y1 DL-hq ≥0.8 (SKHASH)': HASH3_DIR_PG / 'OUT_Y1_dl_hq_sk' / 'out.txt',
        'Y2 CC-SVD (SKHASH)':      HASH3_DIR_PG / 'OUT_Y2_cc_sk'    / 'out.txt',
        'Y2 DL-hq ≥0.8 (SKHASH)': HASH3_DIR_PG / 'OUT_Y2_dl_hq_sk' / 'out.txt',
    }

    sk_datasets = {}
    for lbl, path in sk_tags.items():
        df = load_skhash_pg(path)
        if len(df):
            sub = df[df['quality'].isin(QUAL_FILTER_BB)].copy()
            sub = sub[(sub['origin_lat'].between(*LAT_LIM_BB)) &
                      (sub['origin_lon'].between(*LON_LIM_BB))]
            sk_datasets[lbl] = sub.sort_values('qrank')
        else:
            sk_datasets[lbl] = pd.DataFrame()

    pw_sk = 3.8; ph_sk = pw_sk * abs(LAT_LIM_BB[1]-LAT_LIM_BB[0]) / abs(LON_LIM_BB[1]-LON_LIM_BB[0])
    fig, axes = plt.subplots(1, 4, figsize=(pw_sk*4+0.6, ph_sk+1.4), facecolor='white')
    add_page_header(fig,
                    'SKHASH FM Comparison: CC-SVD vs DL-hq  |  All Depths  |  Y1 & Y2',
                    'Q = A+B only  ·  Blue=Normal  Red=Reverse  Green=Strike-slip  Gray=Oblique  ·  1D velocity')

    for ax, (lbl, sub) in zip(axes, sk_datasets.items()):
        ax.set_facecolor([0.85, 0.92, 0.97])
        ax.set_aspect('equal')
        for _, row in sub.iterrows():
            try:
                bb = _beach([row['strike'], row['dip'], row['rake']],
                            xy=(row['origin_lon'], row['origin_lat']),
                            width=BB_WIDTH_PG, linewidth=0.15,
                            facecolor=FAULT_COLOR_BB[row['fault_type']], alpha=0.85)
                bb.set_transform(ax.transData); bb.set_zorder(2)
                ax.add_collection(bb)
            except Exception:
                pass
        ax.plot(CALDERA_RIM_BB[:,0], CALDERA_RIM_BB[:,1], '-', color='black', lw=1.2, zorder=7)
        ax.plot([v[1] for v in STA_OO_BB.values()],
                [v[0] for v in STA_OO_BB.values()],
                '^', color='yellow', markeredgecolor='k', markeredgewidth=0.3,
                markersize=4, zorder=8)
        ax.set_xlim(LON_LIM_BB); ax.set_ylim(LAT_LIM_BB)
        tks = np.arange(-130.04, -129.96, 0.02)
        ax.set_xticks(tks)
        ax.set_xticklabels([f'{t:.2f}' for t in tks], fontsize=5.5, rotation=45)
        ax.tick_params(axis='y', labelsize=5.5)
        ax.grid(True, lw=0.3, color='gray', alpha=0.4)
        ax.set_xlabel('Longitude', fontsize=7)
        ax.set_title(f'{lbl}\nn={len(sub)}', fontsize=8.5, pad=3)
    axes[0].set_ylabel('Latitude', fontsize=7)
    _sk_leg = [mpatches.Patch(color=FAULT_COLOR_BB[k], label=v) for k,v in
               [('N','Normal'),('R','Reverse'),('S','Strike-slip'),('U','Oblique')]]
    axes[-1].legend(handles=_sk_leg, loc='lower right', fontsize=6.5, framealpha=0.9)
    plt.subplots_adjust(top=0.90, bottom=0.10, left=0.05, right=0.99, wspace=0.15)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

print(f'\nSaved → {OUT_PDF}')
