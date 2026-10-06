"""
run_skhash_cc_dlhq.py
=====================
Write SKHASH inputs and run SKHASH for 4 configurations:
  Y1_cc, Y1_dl_hq, Y2_cc, Y2_dl_hq

Polarities sourced from E_Po_py npy files (CC-SVD svd_final / DL-hq).
1D axial velocity model only.  Minimum 8 polarities.

Outputs (SKHASH out.txt CSV):
  SKHASH/SKHASH7/examples/hash3/OUT_Y1_cc_sk/out.txt
  SKHASH/SKHASH7/examples/hash3/OUT_Y1_dl_hq_sk/out.txt
  SKHASH/SKHASH7/examples/hash3/OUT_Y2_cc_sk/out.txt
  SKHASH/SKHASH7/examples/hash3/OUT_Y2_dl_hq_sk/out.txt

Run with FM_ML env:
  /opt/miniconda3/envs/FM_ML/bin/python run_skhash_cc_dlhq.py
"""
import warnings; warnings.filterwarnings('ignore')

from pathlib import Path
from datetime import datetime
import subprocess, sys

import numpy as np
import h5py

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT    = Path('/Users/mczhang/Documents/GitHub/FM7')
DATA       = PROJECT / '02-data'
SKHASH_DIR = Path('/Users/mczhang/Documents/GitHub/SKHASH/SKHASH7')
HASH3_DIR  = SKHASH_DIR / 'examples' / 'hash3'
OUT_E      = DATA / 'E_Po_py'

# ── Parameters ─────────────────────────────────────────────────────────────────
NPOLMIN  = 8
DL_CONF  = 0.8
DEP_MIN  = 0.5
EH, EZ   = 0.3, 0.2
MAG      = 1.0
CHAN     = 'HHZ'
NET      = 'OO'
ONSET    = 'I'

STATION_META = {
    'AS1': (45.93356, -129.99920,  60), 'AS2': (45.93377, -130.01410,  44),
    'CC1': (45.95468, -130.00890,  61), 'EC1': (45.94958, -129.97970,  77),
    'EC2': (45.93967, -129.97380,  70), 'EC3': (45.93607, -129.97850,  73),
    'ID1': (45.92573, -129.97800,  62),
    '01A': (46.01933, -130.00538,   0), '02A': (46.00063, -130.02926,  34),
    '03A': (45.99292, -129.99278, 102), '04A': (45.98874, -130.04840,  53),
    '05A': (45.97439, -129.97797,  73), '06A': (45.98294, -130.01413,   9),
    '07A': (45.96907, -130.02985,   8), '08A': (45.96415, -130.00450,  46),
    '09A': (45.97084, -130.06239, 103), '10A': (45.96034, -129.94999,  20),
    '11A': (45.94959, -130.03938, 158), '12A': (45.91537, -130.02178,  73),
    '13A': (45.90782, -129.97156,  24), '14A': (45.89997, -130.00757,   1),
    '15A': (45.91972, -129.93942,   0),
    '01B': (46.01935, -130.00526,   2), '02B': (46.00111, -130.02930,  32),
    '03B': (45.99250, -129.99397, 107), '04B': (45.98834, -130.04896,  55),
    '05B': (45.97327, -129.97728,  74), '06B': (45.98178, -130.01423,  12),
    '07B': (45.96950, -130.02963,  11), '08B': (45.96336, -130.00441,  57),
    '09B': (45.97194, -130.05972, 111), '10B': (45.95867, -129.94982,  23),
    '11B': (45.94851, -130.03801, 160), '12B': (45.91473, -130.02173,  78),
    '13B': (45.90472, -129.97240,  20), '14B': (45.89913, -130.00612,   1),
}

OO_STAS = ['AS1', 'AS2', 'CC1', 'EC1', 'EC2', 'EC3', 'ID1']

def short_to_sk(s): return 'A' + s

def deg_to_degmin(d):
    d = abs(d); deg = int(d); mins = 60.0 * (d - deg)
    return deg, mins

def load_loc_file(path):
    locs = {}
    with open(path) as f:
        for line in f:
            p = line.split()
            if len(p) < 11: continue
            try:
                yr,mo,da = int(p[0]),int(p[1]),int(p[2])
                hh,mn    = int(p[3]),int(p[4])
                ss       = float(p[5])
                lat,lon,depth = float(p[6]),float(p[7]),float(p[8])
                eid = int(p[10])
                sec_i = int(ss); usec = int(round((ss-sec_i)*1e6))
                locs[eid] = (datetime(yr,mo,da,hh,mn,sec_i,usec), lon, lat, depth)
            except (ValueError, IndexError):
                continue
    return locs

def load_amplitudes(year, stations):
    wave_h5 = DATA / 'A_ID' / f'waves_{year}_ge8.h5'
    eid_idx = {}; amp_data = {}
    if not wave_h5.exists():
        return eid_idx, amp_data
    with h5py.File(wave_h5, 'r') as f:
        eids = f['event_ids'][:]
        eid_idx = {int(e): i for i,e in enumerate(eids)}
        for sta in stations:
            if f'amplitudes/{sta}' in f:
                amp_data[sta] = f[f'amplitudes/{sta}'][:]
    return eid_idx, amp_data

def get_amp(eid, sta, eid_idx, amp_data):
    if sta not in amp_data or eid not in eid_idx:
        return 1.0, 1.0, 1.0
    idx = eid_idx[eid]
    amp = amp_data[sta][idx]
    if np.any(np.isnan(amp)): return 1.0, 1.0, 1.0
    return max(float(amp[0]),0.01), max(float(amp[2]),0.01), max(float(amp[1]),0.01)

# ── Load polarities from E_Po_py npy files ─────────────────────────────────────
def load_polarities(year, stations):
    """Returns {eid: {sta: (cc_pol, dl_pol, dl_conf)}}"""
    pol = {}
    for sta in stations:
        p = OUT_E / f'E_{year}_{sta}.npy'
        if not p.exists(): continue
        arr = np.load(p)
        for row in arr:
            eid = int(row[0])
            if eid not in pol: pol[eid] = {}
            pol[eid][sta] = (int(row[2]), int(row[3]), float(row[4]))
    return pol

# ── Write SKHASH inputs for one configuration ──────────────────────────────────
def write_inputs(tag, year, stations, pol_by_eid, eid_idx, amp_data, loc_by_id, mode):
    """
    mode: 'cc'    → svd_final
          'dl_hq' → dl_pol where dl_conf >= DL_CONF
    """
    out_tag  = f'{tag}_sk'
    in_dir   = HASH3_DIR / 'IN'
    out_dir  = HASH3_DIR / f'OUT_{out_tag}'
    out_dir.mkdir(exist_ok=True)

    pha_path = in_dir / f'north2_{out_tag}.txt'
    amp_path = in_dir / f'north3_{out_tag}.txt'
    sta_path = in_dir / f'scsn.stations_{out_tag}.txt'
    cor_path = in_dir / f'north3.statcor_{out_tag}.txt'
    ctl_path = HASH3_DIR / f'control_file_{out_tag}.txt'

    # Station file (sorted as-is — SKHASH doesn't require alpha sort)
    with open(sta_path, 'w') as f:
        for sta in stations:
            if sta not in STATION_META: continue
            lat, lon, elev = STATION_META[sta]
            sn = short_to_sk(sta)
            f.write(f'{sn} HHZ BIG CHUCKAWALLA MTNS             '
                    f'{lat:.5f} {lon:.5f} {elev:5d} 1997/09/19 3000/01/01 OO\n')

    # Statcor file
    with open(cor_path, 'w') as f:
        for sta in stations:
            f.write(f'{short_to_sk(sta)} HHZ OO  0.0000\n')

    # Phase + amplitude files
    n_written = n_noloc = n_npol = 0
    cluster_idx = 0

    with open(pha_path, 'w') as fpha, open(amp_path, 'w') as famp:
        for eid in sorted(pol_by_eid):
            if eid not in loc_by_id:
                n_noloc += 1; continue
            sta_dict = pol_by_eid[eid]
            dt, lon, lat, depth = loc_by_id[eid]
            depth = max(depth, DEP_MIN)

            sta_data = []
            for sta in stations:
                if sta not in sta_dict: continue
                cc_pol, dl_p, dl_c = sta_dict[sta]
                if mode == 'cc':
                    pol_val = cc_pol
                else:  # dl_hq
                    pol_val = dl_p if (not np.isnan(dl_c) and dl_c >= DL_CONF) else 0
                if pol_val == 0: continue
                pol_char = 'U' if pol_val > 0 else 'D'
                sn = short_to_sk(sta)
                noise, p_amp, s_amp = get_amp(eid, sta, eid_idx, amp_data)
                sta_data.append((sn, pol_char, noise, p_amp, s_amp))

            if len(sta_data) < NPOLMIN:
                n_npol += 1; continue

            cluster_idx += 1
            npicks = len(sta_data)

            yr = dt.year; mo = dt.month; da = dt.day
            hr = dt.hour; mn = dt.minute
            sec = dt.second + dt.microsecond / 1e6
            ilat, mlat = deg_to_degmin(lat)
            ilon, mlon = deg_to_degmin(abs(lon))
            ns = 'N' if lat >= 0 else 'S'
            ew = 'W' if lon < 0 else 'E'
            eid_str = f'{cluster_idx:>22}'

            # SKHASH hash3 phase format
            header = (
                f'{yr:4d}{mo:2d}{da:2d}{hr:2d}{mn:2d}{sec:5.2f}'
                f'{ilat:2d}{ns}{mlat:5.2f}{ilon:3d}{ew}{mlon:5.2f}{depth:5.2f}'
                f'     {npicks:3d}{"":41}'
                f'{EH:5.2f}  {EZ:5.2f}   {"":36}'
                f'{MAG:4.2f}{eid_str}\n'
            )
            fpha.write(header)
            for (sn, pol_char, _, _, _) in sta_data:
                fpha.write(f'{sn:4s} {NET:2s}  {CHAN:3s} {ONSET} {pol_char}\n')
            fpha.write(f'{" "*56}{cluster_idx:16d}\n')

            # SKHASH amplitude format (6 values)
            famp.write(f'{cluster_idx}         {npicks}\n')
            for (sn, _, noise, p_amp, s_amp) in sta_data:
                famp.write(
                    f'{sn:4s} {CHAN:3s} {NET:2s}'
                    f'  {noise:5.2f}   {noise:5.2f}'
                    f'  {noise:10.3f} {noise:10.3f}'
                    f' {p_amp:10.3f} {s_amp:10.3f}\n'
                )

            n_written += 1

    # Control file (1D velocity only)
    with open(ctl_path, 'w') as f:
        f.write(f'## SKHASH control — {out_tag}\n\n')
        f.write(f'$input_format\nhash3\n\n')
        f.write(f'$stfile\nexamples/hash3/IN/scsn.stations_{out_tag}.txt\n\n')
        f.write(f'$plfile\nexamples/hash3/IN/scsn.reverse.txt\n\n')
        f.write(f'$corfile\nexamples/hash3/IN/north3.statcor_{out_tag}.txt\n\n')
        f.write(f'$ampfile\nexamples/hash3/IN/north3_{out_tag}.txt\n\n')
        f.write(f'$fpfile\nexamples/hash3/IN/north2_{out_tag}.txt\n\n')
        f.write(f'$outfile1\nexamples/hash3/OUT_{out_tag}/out.txt\n\n')
        f.write(f'$outfile2\nexamples/hash3/OUT_{out_tag}/out2.txt\n\n')
        f.write(f'$npolmin\n{NPOLMIN}\n\n')
        f.write('$dang\n5\n\n$nmc\n30\n\n$maxout\n300\n\n')
        f.write('$ratmin\n2\n\n$badfrac\n0.1\n\n$qbadfrac\n0.3\n\n')
        f.write('$delmax\n25\n\n$cangle\n45\n\n$prob_max\n0.75\n\n')
        f.write('$vmodel_paths\nexamples/velocity_models/velmod_axial1D.txt\n')

    print(f'  {out_tag}: {n_written:,} events -> {ctl_path.name}  '
          f'(no_loc={n_noloc}, <{NPOLMIN}pol={n_npol})')
    return ctl_path, out_dir / 'out.txt', n_written

# ── Main ───────────────────────────────────────────────────────────────────────
configs = [
    ('Y1_cc',     'Y1', OO_STAS + [f'{i:02d}A' for i in range(1,15)], 'cc'),
    ('Y1_dl_hq',  'Y1', OO_STAS + [f'{i:02d}A' for i in range(1,15)], 'dl_hq'),
    ('Y2_cc',     'Y2', OO_STAS + [f'{i:02d}B' for i in range(1,15)], 'cc'),
    ('Y2_dl_hq',  'Y2', OO_STAS + [f'{i:02d}B' for i in range(1,15)], 'dl_hq'),
]

print('Loading polarities and amplitudes...')
pol_Y1 = load_polarities('Y1', OO_STAS + [f'{i:02d}A' for i in range(1,15)])
pol_Y2 = load_polarities('Y2', OO_STAS + [f'{i:02d}B' for i in range(1,15)])
eid_idx_Y1, amp_Y1 = load_amplitudes('Y1', OO_STAS + [f'{i:02d}A' for i in range(1,15)])
eid_idx_Y2, amp_Y2 = load_amplitudes('Y2', OO_STAS + [f'{i:02d}B' for i in range(1,15)])
loc_Y1 = load_loc_file(DATA / 'A_all' / 'axial.Y1.mldd.loc.260317')
loc_Y2 = load_loc_file(DATA / 'A_all' / 'axial.Y2.mldd.loc.260317')
print(f'Y1: {len(pol_Y1):,} events  Y2: {len(pol_Y2):,} events')

pol_map  = {'Y1': pol_Y1,     'Y2': pol_Y2}
eidx_map = {'Y1': eid_idx_Y1, 'Y2': eid_idx_Y2}
amp_map  = {'Y1': amp_Y1,     'Y2': amp_Y2}
loc_map  = {'Y1': loc_Y1,     'Y2': loc_Y2}

print('\nWriting SKHASH inputs...')
run_configs = []
for tag, year, stations, mode in configs:
    ctl, out_txt, n = write_inputs(
        tag, year, stations,
        pol_map[year], eidx_map[year], amp_map[year],
        loc_map[year], mode)
    run_configs.append((tag, ctl, out_txt, n))

print('\nRunning SKHASH...')
skhash_py = SKHASH_DIR / 'SKHASH.py'
for tag, ctl, out_txt, n_events in run_configs:
    print(f'  {tag} ({n_events:,} events)...')
    result = subprocess.run(
        [sys.executable, str(skhash_py), str(ctl)],
        capture_output=True, text=True, cwd=str(SKHASH_DIR))
    if result.returncode != 0:
        print(f'  ERROR: {result.stderr[:300]}')
    else:
        if out_txt.exists():
            n_out = sum(1 for _ in open(out_txt)) - 1  # minus header
            print(f'  -> {out_txt.name}: {n_out:,} FMs')
        else:
            print(f'  WARNING: {out_txt} not found')

print('\nDone. Results in SKHASH/SKHASH7/examples/hash3/OUT_*_sk/')
