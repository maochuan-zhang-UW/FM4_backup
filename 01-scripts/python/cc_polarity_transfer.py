"""
cc_polarity_transfer.py
=======================
Cross-correlation polarity transfer: augment sparse events (< NPOLMIN picks) by
transferring polarities from high-CC template events (which have a confident pick
at the station in question).

Algorithm (per-station):
  templates = events with confident DL polarity at this station
  targets   = events without a polarity at this station (missing or low-conf)
  For each target, find templates within MAX_DIST_KM:
    compute normalised CC between target waveform and template waveform
    if CC >= CC_MIN: transfer polarity (marked as CC-transferred, conf = cc * template_conf)

Only templates from the same year are used (different 2F hardware; OO stations
share hardware across years but events occur at different times, so same-year only
is conservative and avoids cross-year drift).

Outputs:
  02-data/E_Po/dl_polarity_Y1_ge8_aug.h5
  02-data/E_Po/dl_polarity_Y2_ge8_aug.h5
  02-data/E_Po/cc_pairs.npz   (diagnostics)

Run with FM_ML env:
  /opt/miniconda3/envs/FM_ML/bin/python cc_polarity_transfer.py
"""
import time, warnings
from pathlib import Path

import numpy as np
import h5py
from scipy.spatial import cKDTree

warnings.filterwarnings('ignore')

# ── Config ────────────────────────────────────────────────────────────────────
PROJECT     = Path('/Users/mczhang/Documents/GitHub/FM7')
DATA        = PROJECT / '02-data'

WAVES_Y1    = DATA / 'A_ID' / 'waves_Y1_ge8.h5'
WAVES_Y2    = DATA / 'A_ID' / 'waves_Y2_ge8.h5'
POL_Y1_IN   = DATA / 'E_Po' / 'dl_polarity_Y1_ge8.h5'
POL_Y2_IN   = DATA / 'E_Po' / 'dl_polarity_Y2_ge8.h5'
POL_Y1_OUT  = DATA / 'E_Po' / 'dl_polarity_Y1_ge8_aug.h5'
POL_Y2_OUT  = DATA / 'E_Po' / 'dl_polarity_Y2_ge8_aug.h5'

CC_MIN      = 0.80   # minimum CC to transfer
CONF_MIN    = 0.70   # minimum template DL confidence to use as source
NPOLMIN     = 15     # target minimum total picks
WIN         = 64     # waveform window samples for CC (64 @ 200Hz = 0.32s)
LAG_MAX     = 10     # +/- max lag samples (~50ms at 200Hz)
MAX_DIST_KM = 3.0    # spatial radius for template search (km)

OO_STAS  = ['AS1', 'AS2', 'CC1', 'EC1', 'EC2', 'EC3', 'ID1']
TF_Y1    = [f'{i:02d}A' for i in range(1, 15)] + ['15A']
TF_Y2    = [f'{i:02d}B' for i in range(1, 15)] + ['15A']
STAS_Y1  = OO_STAS + TF_Y1
STAS_Y2  = OO_STAS + TF_Y2

# ── Helpers ───────────────────────────────────────────────────────────────────
def haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0
    la1, lo1, la2, lo2 = map(np.radians, [lat1, lon1, lat2, lon2])
    a = np.sin((la2-la1)/2)**2 + np.cos(la1)*np.cos(la2)*np.sin((lo2-lo1)/2)**2
    return 2*R*np.arcsin(np.sqrt(np.clip(a, 0, 1)))

# 1 degree latitude ≈ 111 km; longitude correction at Axial latitude (~45.95°N)
LAT_SCALE = 111.0          # km per degree lat
LON_SCALE = 111.0 * np.cos(np.radians(45.95))   # km per degree lon


def load_locs(loc_file):
    """Return arrays: event_ids, lats, lons (parallel arrays)."""
    rows = []
    with open(loc_file) as f:
        for line in f:
            p = line.split()
            if len(p) < 11:
                continue
            try:
                rows.append((int(p[10]), float(p[6]), float(p[7])))
            except ValueError:
                continue
    if not rows:
        return np.array([], dtype=np.int64), np.array([]), np.array([])
    arr = np.array(rows)
    return arr[:, 0].astype(np.int64), arr[:, 1], arr[:, 2]


def load_year(waves_h5, pol_h5, stations):
    with h5py.File(waves_h5, 'r') as hf:
        eids = hf['event_ids'][:]
        waves = {s: hf[f'waveforms/{s}'][:] for s in stations if f'waveforms/{s}' in hf}
    with h5py.File(pol_h5, 'r') as hf:
        pol  = {s: hf[f'polarity/{s}'][:].copy()   for s in stations if f'polarity/{s}'   in hf}
        conf = {s: hf[f'confidence/{s}'][:].copy() for s in stations if f'confidence/{s}' in hf}
    return eids, waves, pol, conf


def count_picks(pol_dict, N):
    c = np.zeros(N, dtype=np.int32)
    for arr in pol_dict.values():
        c += (arr != 0).astype(np.int32)
    return c


def compute_cc_1d(target_w, template_mat, lag_max=LAG_MAX, win=WIN):
    """
    Compute normalised CC between one target waveform and a matrix of templates.

    target_w    : (win,) float32
    template_mat: (M, win) float32
    Returns     : (M,) float32 CC values in [0, 1]
    """
    # Normalise
    tn = target_w / (np.linalg.norm(target_w) + 1e-12)
    norms = np.linalg.norm(template_mat, axis=1, keepdims=True) + 1e-12
    tm = template_mat / norms  # (M, win)

    # Batch FFT CC
    pad = int(2 ** np.ceil(np.log2(2 * win - 1)))
    FT = np.fft.rfft(tn, n=pad)           # (freq,)
    FM = np.fft.rfft(tm, n=pad, axis=1)   # (M, freq)
    xcorr = np.fft.irfft(FT.conj() * FM, n=pad, axis=1)  # (M, pad)

    # Consider lags ±lag_max
    pos  = xcorr[:, :lag_max + 1]
    neg  = xcorr[:, -lag_max:]
    peak = np.max(np.abs(np.concatenate([pos, neg], axis=1)), axis=1)
    return peak.astype(np.float32)


def run_cc_transfer(year, eids, waves, pol, conf, loc_eids, lats, lons, stations):
    """
    Run CC transfer for one year.  Returns augmented copies of pol and conf.
    """
    N = len(eids)
    aug_pol  = {s: pol[s].copy()  for s in pol}
    aug_conf = {s: conf[s].copy() for s in conf}

    # Map event_id → index in eids array
    eid_to_idx = {int(eid): i for i, eid in enumerate(eids)}

    # Map event_id → location
    loc_map = {int(loc_eids[k]): (lats[k], lons[k])
               for k in range(len(loc_eids))
               if int(loc_eids[k]) in eid_to_idx}

    # Build coordinate arrays aligned with eids
    ev_lats = np.full(N, np.nan)
    ev_lons = np.full(N, np.nan)
    for eid, (la, lo) in loc_map.items():
        idx = eid_to_idx.get(eid)
        if idx is not None:
            ev_lats[idx] = la
            ev_lons[idx] = lo

    has_loc = ~np.isnan(ev_lats)
    print(f'  {year}: {has_loc.sum():,}/{N} events have locations')

    # Build KD-tree for spatial lookup (only events with locations)
    loc_idx = np.where(has_loc)[0]
    # Convert to km-like coordinates so radius search works in km units
    coords_km = np.column_stack([
        ev_lats[loc_idx] * LAT_SCALE,
        ev_lons[loc_idx] * LON_SCALE,
    ])
    tree = cKDTree(coords_km)

    total_transfers = 0
    t0 = time.time()

    for sta in stations:
        if sta not in waves:
            continue
        wmat = waves[sta]    # (N, 264) float32
        valid_wave = ~np.isnan(wmat[:, 0])

        if sta not in pol:
            continue
        pol_arr  = pol[sta]   # original
        conf_arr = conf[sta]

        # Templates: events with valid waveform AND confident pick at this station
        tmpl_mask = (valid_wave
                     & (pol_arr != 0)
                     & (~np.isnan(conf_arr))
                     & (conf_arr >= CONF_MIN))
        tmpl_idx = np.where(tmpl_mask)[0]
        if len(tmpl_idx) == 0:
            continue

        # Targets: events with valid waveform but NO polarity at this station
        tgt_mask = (valid_wave
                    & ((aug_pol[sta] == 0)
                       | np.isnan(aug_conf.get(sta, conf_arr))))
        tgt_idx = np.where(tgt_mask)[0]
        if len(tgt_idx) == 0:
            continue

        # For efficient lookup: map template indices to their km coordinates
        tmpl_with_loc = tmpl_idx[has_loc[tmpl_idx]]
        if len(tmpl_with_loc) == 0:
            continue
        tmpl_coords_km = np.column_stack([
            ev_lats[tmpl_with_loc] * LAT_SCALE,
            ev_lons[tmpl_with_loc] * LON_SCALE,
        ])
        tmpl_tree = cKDTree(tmpl_coords_km)

        # Waveform windows for templates
        tmpl_wins = wmat[tmpl_with_loc, :WIN].copy()   # (M_tmpl, WIN)

        sta_transfers = 0
        for gi in tgt_idx:
            if not has_loc[gi]:
                continue
            query_pt = np.array([ev_lats[gi] * LAT_SCALE,
                                  ev_lons[gi] * LON_SCALE])
            # Find templates within MAX_DIST_KM
            near_tmpl_local = tmpl_tree.query_ball_point(query_pt, MAX_DIST_KM)
            if not near_tmpl_local:
                continue

            near_tmpl_local = np.array(near_tmpl_local, dtype=np.int32)
            near_tmpl_global = tmpl_with_loc[near_tmpl_local]   # global indices

            # Compute CC
            target_win   = wmat[gi, :WIN]
            template_mat = tmpl_wins[near_tmpl_local]
            cc_vals      = compute_cc_1d(target_win, template_mat)

            # Find best template with CC >= CC_MIN
            best_local = np.argmax(cc_vals)
            if cc_vals[best_local] < CC_MIN:
                continue

            best_global = near_tmpl_global[best_local]
            tmpl_pol    = int(pol_arr[best_global])
            tmpl_conf   = float(conf_arr[best_global])
            transferred_conf = float(cc_vals[best_local]) * tmpl_conf

            aug_pol[sta][gi]  = tmpl_pol
            aug_conf[sta][gi] = transferred_conf
            sta_transfers += 1

        total_transfers += sta_transfers
        elapsed = time.time() - t0
        n_total_conf = int(((aug_pol[sta] != 0) & (aug_conf[sta] >= 0.5)).sum())
        print(f'  {year} {sta}: {len(tmpl_with_loc):>4} templates  '
              f'{len(tgt_idx):>5} targets  '
              f'{sta_transfers:>5} transfers  '
              f'total_conf={n_total_conf:>6,}  '
              f'({elapsed:.0f}s elapsed)')

    print(f'  {year}: {total_transfers:,} total transfers in {time.time()-t0:.0f}s')
    return aug_pol, aug_conf


# ── Main ──────────────────────────────────────────────────────────────────────
print('=== CC Polarity Transfer ===')
print(f'CC_MIN={CC_MIN}  CONF_MIN={CONF_MIN}  NPOLMIN={NPOLMIN}  '
      f'MAX_DIST={MAX_DIST_KM}km  WIN={WIN}samp\n')

# Load locations
LOC_Y1 = DATA / 'A_all' / 'axial.Y1.mldd.loc.260317'
LOC_Y2 = DATA / 'A_all' / 'axial.Y2.mldd.loc.260317'
loc_eids_y1, lats_y1, lons_y1 = load_locs(LOC_Y1)
loc_eids_y2, lats_y2, lons_y2 = load_locs(LOC_Y2)
print(f'Locations: Y1={len(loc_eids_y1):,}  Y2={len(loc_eids_y2):,}')

# Load waveforms and polarities
print('\nLoading Y1...')
eids_y1, waves_y1, pol_y1, conf_y1 = load_year(WAVES_Y1, POL_Y1_IN, STAS_Y1)
print(f'  Y1: {len(eids_y1):,} events, {len(waves_y1)} wave arrays')

print('Loading Y2...')
eids_y2, waves_y2, pol_y2, conf_y2 = load_year(WAVES_Y2, POL_Y2_IN, STAS_Y2)
print(f'  Y2: {len(eids_y2):,} events, {len(waves_y2)} wave arrays')

# Print initial pick distributions
for year, pol_d, N in [('Y1', pol_y1, len(eids_y1)), ('Y2', pol_y2, len(eids_y2))]:
    c = count_picks(pol_d, N)
    print(f'\nInitial {year} picks:')
    for thr in [8, 10, 12, 15, 18]:
        n = int((c >= thr).sum())
        print(f'  ≥{thr:2d}: {n:>6,}  ({n/N*100:.1f}%)')

# CC transfer
print('\n=== Running CC transfer ===')
print('Y1:')
aug_pol_y1, aug_conf_y1 = run_cc_transfer(
    'Y1', eids_y1, waves_y1, pol_y1, conf_y1,
    loc_eids_y1, lats_y1, lons_y1, STAS_Y1)

print('\nY2:')
aug_pol_y2, aug_conf_y2 = run_cc_transfer(
    'Y2', eids_y2, waves_y2, pol_y2, conf_y2,
    loc_eids_y2, lats_y2, lons_y2, STAS_Y2)

# Post-transfer distribution
for year, pol_d, N in [('Y1', aug_pol_y1, len(eids_y1)), ('Y2', aug_pol_y2, len(eids_y2))]:
    c = count_picks(pol_d, N)
    print(f'\nAugmented {year} picks:')
    for thr in [8, 10, 12, 15, 18]:
        n = int((c >= thr).sum())
        print(f'  ≥{thr:2d}: {n:>6,}  ({n/N*100:.1f}%)')

# Save augmented H5 files
def save_aug(out_path, eids, aug_pol, aug_conf, stations):
    print(f'\nSaving → {out_path}')
    with h5py.File(out_path, 'w') as hf:
        hf.create_dataset('event_ids', data=eids)
        for sta in stations:
            if sta in aug_pol:
                hf.create_dataset(f'polarity/{sta}',   data=aug_pol[sta])
                hf.create_dataset(f'confidence/{sta}', data=aug_conf[sta])
    print(f'  {out_path.name}  {out_path.stat().st_size/1e6:.1f} MB')

save_aug(POL_Y1_OUT, eids_y1, aug_pol_y1, aug_conf_y1, STAS_Y1)
save_aug(POL_Y2_OUT, eids_y2, aug_pol_y2, aug_conf_y2, STAS_Y2)

print('\ncc_polarity_transfer.py complete.')
