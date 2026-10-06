"""
cut_waves_Y2_ge18.py
====================
Cut waveforms for Y2 events with ≥18 P-picks from two sources:
  - 2F stations (01B–14B, 15A): FM4/datamseed/{year}/{month}/{HH}.mseed  (obspy)
  - OO stations (AS1–ID1):      Axial-AutoLocate/{year}/{month}/{HH}.mat  (scipy)

Output: 02-data/A_ID/waves_Y2_ge18.h5
  /waveforms/<STA>   shape (N_events, 264)  float32  — NaN where unavailable
  /amplitudes/<STA>  shape (N_events, 3)    float32  — [noise_amp, S_amp, P_amp]
  /event_ids         shape (N_events,)       int64

Usage:
  /opt/miniconda3/bin/python3 cut_waves_Y2_ge18.py
"""

import sys, time, warnings
from pathlib import Path
from datetime import datetime, timedelta
from collections import defaultdict
import multiprocessing as mp

import numpy as np
import scipy.io as sio
import h5py
import obspy
from obspy import UTCDateTime

warnings.filterwarnings('ignore')

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECT      = Path('/Users/mczhang/Documents/GitHub/FM7')
MSEED_ROOT   = Path('/Users/mczhang/Documents/GitHub/FM4/01-scripts/datamseed')
OO_MAT_ROOT  = Path('/Users/mczhang/Documents/GitHub/Axial-AutoLocate')
PHA_Y2       = PROJECT / '02-data/A_all/axial.Y2.mldd.pha.260317'
LOC_Y2       = PROJECT / '02-data/A_all/axial.Y2.mldd.loc.260317'
OUT_H5       = PROJECT / '02-data/A_ID/waves_Y2_ge18.h5'
OUT_H5.parent.mkdir(parents=True, exist_ok=True)

# ── Parameters ────────────────────────────────────────────────────────────────
MIN_STA   = 18
FS        = 200.0
WIN_START = -0.32    # s before P pick
WIN_END   =  1.00    # s after P pick
N_SAMPLES = int((WIN_END - WIN_START) * FS)   # 264
BP_LOW    =  4.0
BP_HIGH   = 50.0
N_WORKERS =  8

# Amplitude windows (relative to P pick in seconds)
NOISE_WIN = (-0.32, -0.10)   # within our cut window
P_WIN     = (-0.05,  0.25)
S_WIN     = (-0.10,  0.60)   # relative to S pick

# ── Station definitions ───────────────────────────────────────────────────────
# 2F stations — in FM4 datamseed .mseed files
STATIONS_2F   = ['01B','02B','03B','04B','05B','06B','07B',
                  '08B','09B','10B','11B','12B','13B','14B','15A']
STA_SEED_2F   = {s: f'AX{s}' for s in STATIONS_2F}  # '01B' -> 'AX01B'
# Z-channel per 2F station — HHZ for some, ELZ for others
TWO_F_Z_CHAN  = {
    '01B':'HHZ','02B':'ELZ','03B':'ELZ','04B':'ELZ','05B':'ELZ',
    '06B':'ELZ','07B':'ELZ','08B':'ELZ','09B':'HHZ','10B':'HHZ',
    '11B':'ELZ','12B':'HHZ','13B':'ELZ','14B':'ELZ','15A':'HHZ',
}

# OO stations — in Axial-AutoLocate .mat files
STATIONS_OO   = ['AS1','AS2','CC1','EC1','EC2','EC3','ID1']
STA_SEED_OO   = {'AS1':'AXAS1','AS2':'AXAS2','CC1':'AXCC1',
                  'EC1':'AXEC1','EC2':'AXEC2','EC3':'AXEC3','ID1':'AXID1'}
# Z-channel per OO station (varies by instrument)
OO_Z_CHAN     = {'AXAS1':'EHZ','AXAS2':'EHZ','AXCC1':'HHZ',
                  'AXEC1':'EHZ','AXEC2':'HHZ','AXEC3':'EHZ','AXID1':'EHZ'}

STATIONS_ALL  = STATIONS_OO + STATIONS_2F   # 22 total

# For .pha parsing: full seed→short mapping (both years combined)
STA_MAP_PHY   = {
    'AXAS1':'AS1','AXAS2':'AS2','AXCC1':'CC1','AXEC1':'EC1','AXEC2':'EC2',
    'AXEC3':'EC3','AXID1':'ID1',
    'AX01B':'01B','AX02B':'02B','AX03B':'03B','AX04B':'04B','AX05B':'05B',
    'AX06B':'06B','AX07B':'07B','AX08B':'08B','AX09B':'09B','AX10B':'10B',
    'AX11B':'11B','AX12B':'12B','AX13B':'13B','AX14B':'14B','AX15A':'15A',
}


# ── Helpers ───────────────────────────────────────────────────────────────────
def matlab2dt(dn):
    """Convert MATLAB datenum (float) to Python datetime."""
    return datetime.fromordinal(int(dn)) + timedelta(days=dn % 1) - timedelta(days=366)


def hour_file_path(root, year, month, day, hour, ext):
    return root / f'{year}' / f'{month:02d}' / \
           f'{year}-{month:02d}-{day:02d}-{hour:02d}-00-00.{ext}'


# ── Step 1: Parse .pha ────────────────────────────────────────────────────────
def parse_pha(pha_file, min_sta):
    counts = {}
    picks  = {}
    curr   = None
    with open(pha_file) as f:
        for line in f:
            s = line.strip()
            if not s: continue
            if s.startswith('#'):
                p = s[1:].split()
                try:
                    curr = int(p[13])
                    counts[curr] = 0
                    picks[curr]  = {'P': {}, 'S': {}}
                except (IndexError, ValueError):
                    curr = None
            elif curr is not None:
                p = s.split()
                if len(p) < 4: continue
                seed  = p[0][2:]
                short = STA_MAP_PHY.get(seed)
                if short is None: continue
                try:
                    t_rel = float(p[1])
                    phase = p[3]
                except (ValueError, IndexError):
                    continue
                if phase == 'P':
                    counts[curr] += 1
                    picks[curr]['P'][short] = t_rel
                elif phase == 'S':
                    picks[curr]['S'][short] = t_rel
    return {eid: picks[eid] for eid, n in counts.items() if n >= min_sta}


# ── Step 2: Parse .loc ────────────────────────────────────────────────────────
def parse_loc(loc_file, event_ids):
    target = set(event_ids)
    locs   = {}
    with open(loc_file) as f:
        for line in f:
            p = line.split()
            if len(p) < 11: continue
            try:
                eid = int(p[10])
                if eid not in target: continue
                sec = float(p[5]); si = int(sec)
                us  = min(int((sec - si) * 1_000_000), 999_999)
                locs[eid] = datetime(int(p[0]),int(p[1]),int(p[2]),
                                     int(p[3]),int(p[4]),si,us)
            except: pass
    return locs


# ── Step 3: Group by hour ─────────────────────────────────────────────────────
def group_by_hour(event_ids, origin_times):
    groups = defaultdict(list)
    for eid in event_ids:
        t = origin_times.get(eid)
        if t is None: continue
        groups[(t.year, t.month, t.day, t.hour)].append(eid)
    return groups


# ── Step 4: Waveform cut helper ───────────────────────────────────────────────
def cut_window(data_1d, t_start_dt, p_abs_dt, s_abs_dt=None):
    """
    Cut and return (wave_264, noise_amp, s_amp, p_amp) or None on failure.
    data_1d : 1D float array at FS Hz
    t_start_dt : datetime of data[0]
    p_abs_dt   : datetime of P arrival
    s_abs_dt   : datetime of S arrival (or None)
    """
    sr = FS
    # Index of P pick in data array
    offset_s = (p_abs_dt - t_start_dt).total_seconds()
    i_p = int(round(offset_s * sr))

    i0 = i_p + int(round(WIN_START * sr))
    i1 = i0 + N_SAMPLES
    if i0 < 0 or i1 > len(data_1d):
        return None

    wave = data_1d[i0:i1].astype(np.float32)
    if len(wave) != N_SAMPLES:
        return None

    # Demean + bandpass
    wave = wave - wave.mean()
    from scipy.signal import butter, sosfilt
    sos  = butter(4, [BP_LOW, BP_HIGH], btype='band', fs=sr, output='sos')
    wave = sosfilt(sos, wave).astype(np.float32)

    # Amplitudes (all indices relative to wave[0] = WIN_START)
    def widx(t_rel):
        return int(round((t_rel - WIN_START) * sr))

    ns0 = widx(NOISE_WIN[0]); ns1 = widx(NOISE_WIN[1])
    pp0 = widx(P_WIN[0]);     pp1 = widx(P_WIN[1])
    ns0 = max(0, ns0); ns1 = min(N_SAMPLES, ns1)
    pp0 = max(0, pp0); pp1 = min(N_SAMPLES, pp1)

    noise_amp = float(np.std(wave[ns0:ns1]))   if ns1 > ns0 else 0.0
    p_amp     = float(np.max(np.abs(wave[pp0:pp1]))) if pp1 > pp0 else 0.0

    if s_abs_dt is not None:
        s_offset  = (s_abs_dt - p_abs_dt).total_seconds()
        ss0 = widx(s_offset + S_WIN[0]); ss1 = widx(s_offset + S_WIN[1])
        ss0 = max(0, ss0); ss1 = min(N_SAMPLES, ss1)
        s_amp = float(np.max(np.abs(wave[ss0:ss1]))) if ss1 > ss0 else p_amp
    else:
        s_amp = p_amp

    return wave, noise_amp, s_amp, p_amp


# ── Step 5: Worker function ───────────────────────────────────────────────────
def process_hour(args):
    (year, month, day, hour), event_ids, origin_times, picks = args
    results = []   # list of (eid, sta, wave, noise_amp, s_amp, p_amp)

    # ── Load 2F mseed ─────────────────────────────────────────────────────────
    mseed_f = hour_file_path(MSEED_ROOT, year, month, day, hour, 'mseed')
    traces_2f = {}   # seed_sta -> Trace (HHZ)
    if mseed_f.exists():
        try:
            st = obspy.read(str(mseed_f))
            for tr in st:
                # Each 2F station uses either HHZ or ELZ — accept both Z channels
                if tr.stats.channel.endswith('Z'):
                    # Keep the first Z trace per station (prefer HHZ over ELZ if both)
                    seed_sta = tr.stats.station
                    if seed_sta not in traces_2f or tr.stats.channel == 'HHZ':
                        traces_2f[seed_sta] = tr
        except Exception:
            pass

    # ── Load OO .mat ──────────────────────────────────────────────────────────
    mat_f = hour_file_path(OO_MAT_ROOT, year, month, day, hour, 'mat')
    traces_oo = {}   # seed_sta -> (data_1d, t_start_dt)
    if mat_f.exists():
        try:
            mat = sio.loadmat(str(mat_f))
            for tr in mat['trace'][0]:
                net  = str(tr['network'][0])
                sta  = str(tr['station'][0])
                chan = str(tr['channel'][0])
                if net != 'OO': continue
                expected_chan = OO_Z_CHAN.get(sta)
                if chan != expected_chan: continue
                data_1d  = tr['data'].flatten().astype(np.float64)
                t_start  = matlab2dt(float(tr['startTime'][0][0]))
                traces_oo[sta] = (data_1d, t_start)
        except Exception:
            pass

    # ── Cut for each event × station ──────────────────────────────────────────
    for eid in event_ids:
        t0        = origin_times.get(eid)
        if t0 is None: continue
        ev_p = picks[eid]['P']
        ev_s = picks[eid]['S']

        # ── 2F stations ───────────────────────────────────────────────────────
        for sta in STATIONS_2F:
            p_rel = ev_p.get(sta)
            if p_rel is None: continue
            seed  = STA_SEED_2F[sta]
            tr    = traces_2f.get(seed)
            if tr is None: continue

            p_abs = UTCDateTime(t0) + p_rel
            s_rel = ev_s.get(sta)
            s_abs = UTCDateTime(t0) + s_rel if s_rel else None

            try:
                t_cut0 = p_abs + WIN_START - 1.0
                t_cut1 = p_abs + WIN_END   + 1.0
                seg    = tr.slice(t_cut0, t_cut1).copy()
                if seg.stats.npts < N_SAMPLES: continue

                seg.detrend('demean')
                seg.filter('bandpass', freqmin=BP_LOW, freqmax=BP_HIGH,
                            corners=4, zerophase=True)

                sr  = seg.stats.sampling_rate
                i_p = int(round((p_abs - seg.stats.starttime) * sr))
                i0  = i_p + int(round(WIN_START * sr))
                i1  = i0 + N_SAMPLES
                if i0 < 0 or i1 > seg.stats.npts: continue
                wave = seg.data[i0:i1].astype(np.float32)
                if len(wave) != N_SAMPLES: continue

                # Amplitudes
                def wi(t_rel):
                    return max(0, min(N_SAMPLES, int(round((t_rel - WIN_START) * sr))))
                noise_amp = float(np.std(wave[wi(NOISE_WIN[0]):wi(NOISE_WIN[1])]))
                p_amp     = float(np.max(np.abs(wave[wi(P_WIN[0]):wi(P_WIN[1])])))
                if s_abs:
                    sd = (s_abs - p_abs).total_seconds() if hasattr(s_abs,'total_seconds') else float(s_abs - UTCDateTime(t0) - p_rel)
                    ss0 = wi(sd + S_WIN[0]); ss1 = wi(sd + S_WIN[1])
                    s_amp = float(np.max(np.abs(wave[ss0:ss1]))) if ss1 > ss0 else p_amp
                else:
                    s_amp = p_amp
                results.append((eid, sta, wave, noise_amp, s_amp, p_amp))
            except Exception:
                continue

        # ── OO stations ───────────────────────────────────────────────────────
        for sta in STATIONS_OO:
            p_rel = ev_p.get(sta)
            if p_rel is None: continue
            seed  = STA_SEED_OO[sta]
            entry = traces_oo.get(seed)
            if entry is None: continue
            data_1d, t_start = entry

            p_abs_dt = t0 + timedelta(seconds=p_rel)
            s_rel    = ev_s.get(sta)
            s_abs_dt = t0 + timedelta(seconds=s_rel) if s_rel else None

            try:
                result = cut_window(data_1d, t_start, p_abs_dt, s_abs_dt)
                if result is None: continue
                wave, noise_amp, s_amp, p_amp = result
                results.append((eid, sta, wave, noise_amp, s_amp, p_amp))
            except Exception:
                continue

    return results


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    t0 = time.time()
    print('='*60)
    print('cut_waves_Y2_ge18.py  —  all 22 stations')
    print('='*60)

    print(f'\n[1] Parsing Y2 .pha  (≥{MIN_STA} stations)...')
    picks = parse_pha(PHA_Y2, MIN_STA)
    print(f'    Events: {len(picks):,}')

    print('[2] Parsing Y2 .loc...')
    origin_times = parse_loc(LOC_Y2, list(picks.keys()))
    print(f'    Matched: {len(origin_times):,}')

    event_ids = sorted(origin_times.keys())
    N = len(event_ids)
    print(f'    {N:,} events × {len(STATIONS_ALL)} stations')

    print('[3] Grouping by hour...')
    hour_groups = group_by_hour(event_ids, origin_times)
    print(f'    Unique hours: {len(hour_groups):,}')

    # Only process hours where at least one file exists
    valid_hours = [
        (k, v) for k, v in hour_groups.items()
        if hour_file_path(MSEED_ROOT, *k, 'mseed').exists() or
           hour_file_path(OO_MAT_ROOT, *k, 'mat').exists()
    ]
    print(f'    Hours with ≥1 file on disk: {len(valid_hours):,}')

    print('[4] Creating HDF5...')
    eid_to_idx = {eid: i for i, eid in enumerate(event_ids)}
    with h5py.File(OUT_H5, 'w') as hf:
        hf.create_dataset('event_ids', data=np.array(event_ids, dtype=np.int64))
        for sta in STATIONS_ALL:
            hf.create_dataset(f'waveforms/{sta}',  shape=(N, N_SAMPLES),
                              dtype=np.float32, fillvalue=np.nan)
            hf.create_dataset(f'amplitudes/{sta}', shape=(N, 3),
                              dtype=np.float32, fillvalue=np.nan)
    print(f'    Stations: {len(STATIONS_ALL)}  (OO: {len(STATIONS_OO)}, 2F: {len(STATIONS_2F)})')

    print(f'\n[5] Cutting waveforms ({N_WORKERS} workers)...')
    args_list    = [(k, v, origin_times, picks) for k, v in valid_hours]
    n_done = 0; n_cut = 0; t_log = time.time()

    with h5py.File(OUT_H5, 'a') as hf, \
         mp.Pool(processes=N_WORKERS) as pool:

        for batch in pool.imap_unordered(process_hour, args_list, chunksize=4):
            for (eid, sta, wave, noise_amp, s_amp, p_amp) in batch:
                idx = eid_to_idx[eid]
                hf[f'waveforms/{sta}'][idx]  = wave
                hf[f'amplitudes/{sta}'][idx] = [noise_amp, s_amp, p_amp]
                n_cut += 1
            n_done += 1
            if time.time() - t_log >= 20:
                pct = n_done / len(args_list) * 100
                elapsed = time.time() - t0
                eta     = elapsed / max(n_done,1) * (len(args_list) - n_done)
                print(f'    {pct:5.1f}%  files={n_done}/{len(args_list)}'
                      f'  cuts={n_cut:,}  {elapsed/60:.1f}min  ETA {eta/60:.1f}min')
                t_log = time.time()

    elapsed = time.time() - t0
    print(f'\n[6] Done in {elapsed/60:.1f} min')
    print(f'    Waveforms cut: {n_cut:,} / {N * len(STATIONS_ALL):,}  '
          f'({n_cut / max(N*len(STATIONS_ALL),1)*100:.1f}%)')
    print(f'    Output: {OUT_H5}  ({OUT_H5.stat().st_size/1e6:.1f} MB)')

    print('\n    Per-station success rate:')
    with h5py.File(OUT_H5, 'r') as hf:
        for sta in STATIONS_ALL:
            w = hf[f'waveforms/{sta}'][:, 0]
            ok = int(np.sum(~np.isnan(w)))
            bar = '█' * int(ok/N*20)
            print(f'      {sta:>3}: {ok:>5}/{N}  {ok/N*100:5.1f}%  {bar}')


if __name__ == '__main__':
    mp.set_start_method('spawn', force=True)
    main()
