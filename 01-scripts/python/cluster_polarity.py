"""
cluster_polarity.py
===================
Hierarchical clustering of seismic waveforms + DL-polarity-based absolute polarity
assignment.  Mirrors the MATLAB B_CC → C_SVD → F_Cl pipeline but in a single script.

Algorithm
---------
1.  For each station: compute normalised CC between all event pairs whose waveforms
    are both valid and within MAX_DIST_KM of each other.
2.  Aggregate per-station CC into a global pairwise similarity matrix
    (average across stations where BOTH events have valid waveforms, weighted by
    the number of shared valid stations).
3.  Convert to distance (d = 1 - CC_avg), run SciPy hierarchical clustering
    (average-linkage UPGMA), cut at CC_CLUSTER (default 0.70) → cluster labels.
4.  For each cluster at each station:
      - Gather all members' DL polarities (from dl_polarity_<YEAR>_ge8.h5)
      - Weighted majority vote (weight = DL confidence)
      - Assign the winning polarity to every cluster member (even those that had
        no DL pick at that station)
      - Members with 0 net vote weight keep their original DL pick (or 0 if none)
5.  Events not reachable by any cluster (singletons or too far from others) keep
    their original DL picks.
6.  Save augmented polarity H5 and cluster assignment arrays.

Outputs
-------
    02-data/E_Po/dl_polarity_{YEAR}_ge8_aug.h5
        /polarity/<STA>   int8   (augmented)
        /confidence/<STA> float32
        /event_ids        int64
    02-data/F_Cl/clusters_{YEAR}_ge8.npz
        event_ids, cluster_labels, n_clusters, cc_threshold

Usage (FM_ML env — obspy not required, only scipy+h5py+numpy):
    /opt/miniconda3/envs/FM_ML/bin/python cluster_polarity.py Y1
    /opt/miniconda3/envs/FM_ML/bin/python cluster_polarity.py Y2
"""
import sys, time, warnings
from pathlib import Path

import numpy as np
import h5py
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial import cKDTree

warnings.filterwarnings('ignore')

# ── Args ──────────────────────────────────────────────────────────────────────
if len(sys.argv) < 2 or sys.argv[1] not in ('Y1', 'Y2'):
    print('Usage: python cluster_polarity.py Y1|Y2')
    sys.exit(1)
YEAR = sys.argv[1]

# ── Config ────────────────────────────────────────────────────────────────────
PROJECT      = Path('/Users/mczhang/Documents/GitHub/FM7')
DATA         = PROJECT / '02-data'

WAVES_H5     = DATA / 'A_ID' / f'waves_{YEAR}_ge8.h5'
POL_H5_IN    = DATA / 'E_Po' / f'dl_polarity_{YEAR}_ge8.h5'
POL_H5_OUT   = DATA / 'E_Po' / f'dl_polarity_{YEAR}_ge8_aug.h5'
CLUSTER_NPZ  = DATA / 'F_Cl' / f'clusters_{YEAR}_ge8.npz'
(DATA / 'F_Cl').mkdir(exist_ok=True)

WIN          = 64       # waveform window length (samples, @ 200 Hz = 0.32 s)
LAG_MAX      = 10       # ±10 samples max lag for CC (~50 ms)
CC_CLUSTER   = 0.70     # CC threshold for clustering (distance = 1 - CC)
CC_TRANSFER  = 0.80     # CC threshold for singleton-to-cluster transfer
MAX_DIST_KM  = 3.0      # spatial radius for CC computation
CONF_MIN     = 0.50     # minimum DL confidence to count in majority vote
MIN_CLUSTER  = 3        # minimum cluster size (singletons kept separately)

LAT_SCALE = 111.0
LON_SCALE = 111.0 * np.cos(np.radians(45.95))

OO_STAS  = ['AS1', 'AS2', 'CC1', 'EC1', 'EC2', 'EC3', 'ID1']
TF_STAS  = [f'{i:02d}{"A" if YEAR=="Y1" else "B"}' for i in range(1, 15)] + ['15A']
STATIONS = OO_STAS + TF_STAS

print(f'=== cluster_polarity.py — {YEAR} ===')
print(f'CC_CLUSTER={CC_CLUSTER}  CC_TRANSFER={CC_TRANSFER}  '
      f'MAX_DIST={MAX_DIST_KM}km  WIN={WIN}samp')

# ── Load waveforms ────────────────────────────────────────────────────────────
print(f'\nLoading waveforms from {WAVES_H5.name}...')
t0 = time.time()
with h5py.File(WAVES_H5, 'r') as hf:
    event_ids = hf['event_ids'][:]
    waves = {s: hf[f'waveforms/{s}'][:] for s in STATIONS if f'waveforms/{s}' in hf}
N = len(event_ids)
print(f'  {N:,} events, {len(waves)} station arrays  ({time.time()-t0:.1f}s)')

# ── Load DL polarities ────────────────────────────────────────────────────────
print(f'Loading polarities from {POL_H5_IN.name}...')
with h5py.File(POL_H5_IN, 'r') as hf:
    pol_orig  = {s: hf[f'polarity/{s}'][:].copy()   for s in STATIONS if f'polarity/{s}'   in hf}
    conf_orig = {s: hf[f'confidence/{s}'][:].copy() for s in STATIONS if f'confidence/{s}' in hf}

# ── Load event locations ──────────────────────────────────────────────────────
LOC_FILE = DATA / 'A_all' / f'axial.{YEAR}.mldd.loc.260317'
print(f'Loading locations from {LOC_FILE.name}...')
loc_map = {}
with open(LOC_FILE) as f:
    for line in f:
        p = line.split()
        if len(p) < 11:
            continue
        try:
            loc_map[int(p[10])] = (float(p[6]), float(p[7]))
        except ValueError:
            continue

ev_lats = np.full(N, np.nan)
ev_lons = np.full(N, np.nan)
for i, eid in enumerate(event_ids):
    if int(eid) in loc_map:
        ev_lats[i], ev_lons[i] = loc_map[int(eid)]

has_loc = ~np.isnan(ev_lats)
print(f'  {has_loc.sum():,}/{N} events have locations')

# ── Helpers ───────────────────────────────────────────────────────────────────
def count_picks(pol_d, N):
    c = np.zeros(N, dtype=np.int32)
    for arr in pol_d.values():
        c += (arr != 0).astype(np.int32)
    return c

# ── CC helper ─────────────────────────────────────────────────────────────────
def batch_cc(target_win, template_mat, lag_max=LAG_MAX):
    """
    Normalised CC peak: target (WIN,) vs template_mat (M, WIN).
    Returns (M,) float32.
    """
    tn = target_win / (np.linalg.norm(target_win) + 1e-12)
    norms = np.linalg.norm(template_mat, axis=1, keepdims=True) + 1e-12
    tm = template_mat / norms
    pad = int(2 ** np.ceil(np.log2(2 * WIN - 1)))
    FT = np.fft.rfft(tn, n=pad)
    FM = np.fft.rfft(tm, n=pad, axis=1)
    xcorr = np.fft.irfft(FT.conj() * FM, n=pad, axis=1)
    pos  = xcorr[:, :lag_max + 1]
    neg  = xcorr[:, -lag_max:]
    return np.max(np.abs(np.concatenate([pos, neg], axis=1)), axis=1).astype(np.float32)

# ── Step 1: Per-station CC → global similarity ────────────────────────────────
# Store global CC as sum and count for averaging: cc_sum[i,j], cc_cnt[i,j]
# Use sparse dict of (i,j) pairs (i<j) to avoid N² memory for 20k events.

print('\n[Step 1] Computing per-station CC matrices...')
# KD-tree on all located events
loc_idx   = np.where(has_loc)[0]
coords_km = np.column_stack([ev_lats[loc_idx] * LAT_SCALE,
                              ev_lons[loc_idx] * LON_SCALE])
tree = cKDTree(coords_km)

# Global CC accumulator: only keep pairs that could qualify (CC >= CC_CLUSTER).
# Using two plain dicts for sum and count — avoids the memory overhead of
# storing every spatial pair (most of which have CC < threshold).
# For large datasets (Y2, 20k events), this is critical: storing all ~75M
# spatial pairs as Python dict entries would require ~15GB of memory.
cc_sum = {}   # (i,j) → float sum of per-station CC values
cc_cnt = {}   # (i,j) → int   count of stations contributing

# Pre-filter threshold: only accumulate pairs whose per-station CC is ever high.
# A pair can only reach CC_CLUSTER average if at least one station exceeds it.
CC_PREFILTER = CC_CLUSTER * 0.5   # admit pairs where any station CC ≥ this

t_step1 = time.time()
for sta in STATIONS:
    if sta not in waves:
        continue
    wmat  = waves[sta]
    valid = ~np.isnan(wmat[:, 0])

    valid_loc = valid & has_loc
    vl_idx = np.where(valid_loc)[0]
    if len(vl_idx) < 2:
        continue

    vl_km = np.column_stack([ev_lats[vl_idx] * LAT_SCALE,
                              ev_lons[vl_idx] * LON_SCALE])
    vl_tree = cKDTree(vl_km)
    vl_wins = wmat[vl_idx, :WIN].copy()

    n_pairs_sta = 0
    n_admitted  = 0
    for li, gi in enumerate(vl_idx):
        q = np.array([ev_lats[gi] * LAT_SCALE, ev_lons[gi] * LON_SCALE])
        near_li = vl_tree.query_ball_point(q, MAX_DIST_KM)
        near_li = [k for k in near_li if k > li]
        if not near_li:
            continue
        near_li = np.array(near_li, dtype=np.int32)
        near_gi = vl_idx[near_li]
        cc_vals = batch_cc(vl_wins[li], vl_wins[near_li])
        n_pairs_sta += len(near_li)

        # Only accumulate pairs where this station's CC is above prefilter
        high_mask = cc_vals >= CC_PREFILTER
        for k in np.where(high_mask)[0]:
            key = (int(gi), int(near_gi[k]))
            v   = float(cc_vals[k])
            if key in cc_sum:
                cc_sum[key] += v
                cc_cnt[key] += 1
            else:
                cc_sum[key] = v
                cc_cnt[key] = 1
            n_admitted += 1

    print(f'  {sta}: {len(vl_idx):>6,} valid  {n_pairs_sta:>8,} spatial pairs  '
          f'{n_admitted:>8,} admitted  ({time.time()-t_step1:.0f}s)')

print(f'  Total admitted pairs: {len(cc_sum):,}  ({time.time()-t_step1:.0f}s)')

# ── Step 2: Build condensed distance matrix for clustering ───────────────────
print('\n[Step 2] Building global CC averages...')

# Average CC per pair, keep only those above CC_CLUSTER
qualifying = {k: cc_sum[k] / cc_cnt[k]
              for k in cc_sum
              if cc_sum[k] / cc_cnt[k] >= CC_CLUSTER}
print(f'  Pairs with avg CC ≥ {CC_CLUSTER}: {len(qualifying):,}')

if len(qualifying) == 0:
    print('  No qualifying pairs — all events remain as singletons.')
    cluster_labels = np.zeros(N, dtype=np.int32)  # all in cluster 0 (singleton)
else:
    # Collect all events that participate in at least one qualifying pair
    involved = sorted(set(i for i,j in qualifying) | set(j for i,j in qualifying))
    M = len(involved)
    print(f'  Events in qualifying pairs: {M:,}')

    # Map to local indices 0..M-1
    g2l = {g: l for l, g in enumerate(involved)}

    # Build M×M condensed distance matrix (upper triangle, row-major)
    # scipy needs condensed form: n*(n-1)/2 elements
    dist_condensed = np.ones(M*(M-1)//2, dtype=np.float32)  # default dist=1 (no CC)
    for (gi, gj), avg_cc in qualifying.items():
        li = g2l[gi]; lj = g2l[gj]
        if li > lj:
            li, lj = lj, li
        idx = li * M - li*(li+1)//2 + (lj - li - 1)
        dist_condensed[idx] = max(0.0, 1.0 - avg_cc)

    # ── Step 3: Hierarchical clustering ──────────────────────────────────────
    print('\n[Step 3] Running hierarchical clustering (average linkage)...')
    t3 = time.time()
    Z = linkage(dist_condensed, method='average')
    # Cut at distance threshold = 1 - CC_CLUSTER
    local_labels = fcluster(Z, t=1.0 - CC_CLUSTER, criterion='distance')
    print(f'  Done in {time.time()-t3:.1f}s  '
          f'{local_labels.max()} clusters from {M} events')

    # Map back to global event indices; events not in `involved` → label 0 (unclustered)
    cluster_labels = np.zeros(N, dtype=np.int32)
    for l, g in enumerate(involved):
        cluster_labels[g] = int(local_labels[l])

# ── Cluster stats ─────────────────────────────────────────────────────────────
uniq, cnts = np.unique(cluster_labels, return_counts=True)
# Exclude label 0 (unclustered) from cluster stats
mask_real = uniq > 0
uniq_real = uniq[mask_real]; cnts_real = cnts[mask_real]
n_unclustered = int((cluster_labels == 0).sum())
big_clusters  = uniq_real[cnts_real >= MIN_CLUSTER]
print(f'\nCluster summary:')
print(f'  Total clusters (label > 0): {len(uniq_real):,}')
print(f'  Clusters ≥ {MIN_CLUSTER} members: {len(big_clusters):,}')
print(f'  Unclustered (label=0)    : {n_unclustered:,}')
# Size distribution (real clusters only)
size_bins = [1, 2, 5, 10, 20, 50, 100, 500, 99999]
for lo, hi in zip(size_bins, size_bins[1:]):
    n = int(((cnts_real >= lo) & (cnts_real < hi)).sum())
    if n > 0:
        print(f'  size [{lo:4d},{hi:4d}): {n:>5,} clusters')

# ── Step 4: Polarity assignment via cluster majority vote ─────────────────────
print('\n[Step 4] Assigning polarities by cluster majority vote...')

aug_pol  = {s: pol_orig[s].copy()  for s in pol_orig}
aug_conf = {s: conf_orig[s].copy() for s in conf_orig}

# For each cluster with ≥ MIN_CLUSTER members, at each station:
#   weighted majority vote of DL polarity (weight = confidence)
n_assigned = 0
n_events_augmented = 0
augmented_flag = np.zeros(N, dtype=bool)

cluster_polarity_store = {}  # (cluster_id, sta) → (pol, vote_strength)

for cl_id in big_clusters:
    members = np.where(cluster_labels == cl_id)[0]
    for sta in STATIONS:
        if sta not in pol_orig:
            continue
        pol_arr  = pol_orig[sta][members]    # raw DL polarity
        conf_arr = conf_orig[sta][members] if sta in conf_orig else np.zeros(len(members))

        # Weighted vote: only members with confident DL pick
        up_weight = 0.0
        dn_weight = 0.0
        for m_idx, (p, c) in enumerate(zip(pol_arr, conf_arr)):
            if p == 0 or np.isnan(c) or c < CONF_MIN:
                continue
            w = float(c)
            if p == 1:
                up_weight += w
            else:
                dn_weight += w

        total_vote = up_weight + dn_weight
        if total_vote == 0:
            continue   # no DL picks at all for this station in this cluster

        cluster_pol  = 1 if up_weight >= dn_weight else -1
        vote_strength = max(up_weight, dn_weight) / total_vote   # 0.5 – 1.0

        cluster_polarity_store[(int(cl_id), sta)] = (cluster_pol, vote_strength)

        # Assign to all members
        for m in members:
            old_pol  = int(aug_pol[sta][m])
            old_conf = float(aug_conf[sta][m]) if not np.isnan(aug_conf[sta][m]) else 0.0
            new_conf = vote_strength * float(conf_arr[np.where(pol_arr != 0)[0][0]]
                                             if (pol_arr != 0).any() else 0.5)
            # Only update if cluster vote is more confident than existing pick
            if old_pol == 0 or new_conf > old_conf:
                aug_pol[sta][m]  = cluster_pol
                aug_conf[sta][m] = vote_strength
                if old_pol == 0:
                    n_assigned += 1
                    augmented_flag[m] = True

print(f'  New picks assigned via clustering: {n_assigned:,}')
print(f'  Events augmented (gained ≥1 pick): {augmented_flag.sum():,}')

# ── Step 5: Master-CC transfer for singletons ─────────────────────────────────
# For singletons: compute average CC across ALL shared valid stations between
# singleton and nearby cluster-member templates. If avg CC >= CC_TRANSFER,
# transfer the template's polarities to the singleton at stations where the
# singleton has no valid waveform (thus no DL pick).
# This is the correct physical approach: waveform similarity across multiple
# stations implies same mechanism → same polarity at unobserved stations.
print('\n[Step 5] Master-CC transfer for singletons / unclustered events...')
singleton_mask = (cluster_labels == 0) | np.isin(
    cluster_labels,
    uniq_real[cnts_real < MIN_CLUSTER]
)
n_singletons = int(singleton_mask.sum())
print(f'  {n_singletons:,} singleton/small-cluster events to augment')

# Build list of template events: cluster members (≥ MIN_CLUSTER) with ≥12 orig picks
tmpl_pol_count = count_picks(aug_pol, N)
tmpl_cand_mask = ~singleton_mask & (tmpl_pol_count >= 12) & has_loc
tmpl_cand_idx  = np.where(tmpl_cand_mask)[0]
print(f'  {len(tmpl_cand_idx):,} template candidates (cluster members, ≥12 picks)')

# Build template spatial tree
tmpl_km  = np.column_stack([ev_lats[tmpl_cand_idx] * LAT_SCALE,
                              ev_lons[tmpl_cand_idx] * LON_SCALE])
tmpl_tree = cKDTree(tmpl_km)

# Stations with enough waveform coverage to use for CC similarity
cc_stations = [s for s in STATIONS if s in waves and
               int((~np.isnan(waves[s][:, 0])).sum()) > N * 0.3]
print(f'  {len(cc_stations)} stations used for master-CC computation')

t5 = time.time()
singleton_transfers = 0

# Targets: singletons that still have < NPOLMIN picks AND have location
NPOLMIN_TGT = 15
tgt_pol_count = count_picks(aug_pol, N)
tgt_mask = (singleton_mask & has_loc & (tgt_pol_count < NPOLMIN_TGT))
tgt_idx  = np.where(tgt_mask)[0]
print(f'  {len(tgt_idx):,} target singletons need polarity boost (< {NPOLMIN_TGT} picks)')

for ti, gi in enumerate(tgt_idx):
    if ti % 500 == 0 and ti > 0:
        print(f'    {ti}/{len(tgt_idx)}  transfers={singleton_transfers}  ({time.time()-t5:.0f}s)')

    q = np.array([ev_lats[gi] * LAT_SCALE, ev_lons[gi] * LON_SCALE])
    near_local = np.array(tmpl_tree.query_ball_point(q, MAX_DIST_KM), dtype=np.int32)
    if len(near_local) == 0:
        continue
    near_global = tmpl_cand_idx[near_local]   # global indices of nearby templates

    # Compute master CC: average across shared valid stations
    cc_sum  = np.zeros(len(near_local), dtype=np.float32)
    cc_cnt  = np.zeros(len(near_local), dtype=np.int32)

    for sta in cc_stations:
        if sta not in waves:
            continue
        wmat = waves[sta]
        valid_s = ~np.isnan(wmat[:, 0])
        if not valid_s[gi]:
            continue   # target has no waveform at this station → can't compute CC
        tgt_win = wmat[gi, :WIN]

        # Which of the nearby templates have valid waveform at this station?
        has_sta = valid_s[near_global]  # (len(near_local),) bool
        if not has_sta.any():
            continue

        pos         = np.where(has_sta)[0]            # positions in near_local
        global_near = near_global[pos]                # global event indices
        tmpl_wins   = wmat[global_near, :WIN].copy()
        cc_vals     = batch_cc(tgt_win, tmpl_wins)    # (len(pos),)

        cc_sum[pos] += cc_vals
        cc_cnt[pos] += 1

    # Average CC for templates that shared ≥ 3 stations with target
    valid_pairs = cc_cnt >= 3
    if not valid_pairs.any():
        continue
    master_cc = np.full(len(near_local), -1.0, dtype=np.float32)
    master_cc[valid_pairs] = cc_sum[valid_pairs] / cc_cnt[valid_pairs]

    best_local = np.argmax(master_cc)
    if master_cc[best_local] < CC_TRANSFER:
        continue
    best_global = near_global[best_local]
    best_cc     = float(master_cc[best_local])

    # Transfer polarities for all stations where target has no valid waveform
    n_transferred = 0
    for sta in STATIONS:
        if sta not in aug_pol:
            continue
        if sta in waves and not np.isnan(waves[sta][gi, 0]):
            continue   # target already has waveform here → DL pick exists, keep it

        # Assign template polarity at this station to the singleton
        tmpl_pol  = int(aug_pol[sta][best_global])
        tmpl_conf = float(aug_conf[sta][best_global]) if not np.isnan(aug_conf[sta][best_global]) else 0.0
        if tmpl_pol == 0 or tmpl_conf < CONF_MIN:
            continue

        new_conf = best_cc * tmpl_conf
        if aug_pol[sta][gi] == 0 or new_conf > float(aug_conf[sta][gi]):
            aug_pol[sta][gi]  = tmpl_pol
            aug_conf[sta][gi] = new_conf
            n_transferred += 1

    singleton_transfers += n_transferred

print(f'  Singleton master-CC transfers: {singleton_transfers:,}  ({time.time()-t5:.0f}s)')

orig_counts = count_picks(pol_orig, N)
aug_counts  = count_picks(aug_pol,  N)

print('\n=== Pick count distribution ===')
print(f'{"Threshold":>10}  {"Original":>10}  {"Augmented":>10}  {"Gain":>8}')
for thr in [8, 10, 12, 15, 18, 20]:
    no = int((orig_counts >= thr).sum())
    na = int((aug_counts  >= thr).sum())
    print(f'{thr:>10}  {no:>10,}  {na:>10,}  {na-no:>+8,}')

# ── Save augmented polarities ─────────────────────────────────────────────────
print(f'\nSaving augmented polarities → {POL_H5_OUT}')
with h5py.File(POL_H5_OUT, 'w') as hf:
    hf.create_dataset('event_ids', data=event_ids)
    for sta in STATIONS:
        if sta in aug_pol:
            hf.create_dataset(f'polarity/{sta}',   data=aug_pol[sta])
            hf.create_dataset(f'confidence/{sta}', data=aug_conf[sta])
print(f'  {POL_H5_OUT.name}  {POL_H5_OUT.stat().st_size/1e6:.1f} MB')

# ── Save cluster assignments ──────────────────────────────────────────────────
print(f'Saving cluster assignments → {CLUSTER_NPZ}')
np.savez(CLUSTER_NPZ,
         event_ids=event_ids,
         cluster_labels=cluster_labels,
         n_clusters=len(big_clusters),
         cc_threshold=CC_CLUSTER)
print(f'  {CLUSTER_NPZ.name}  {CLUSTER_NPZ.stat().st_size/1e3:.0f} KB')

print(f'\ncluster_polarity.py {YEAR} complete.')
