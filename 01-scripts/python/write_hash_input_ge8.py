"""
write_hash_input_ge8.py
=======================
Write HASH-format input files for ≥8-station ge8 events (Y1+Y2 combined),
using the 1D axial_1D velocity model only.

Writes to /Users/mczhang/Documents/GitHub/FM7/01-scripts/HASH/:
  phase_ge8_1D.dat
  amp_ge8_1D.dat
  station_ge8.dat
  hash.input_ge8_1D

The event_id → cluster_idx mapping is saved to:
  02-data/G_FM/hash_ge8_1D_eid_map.csv

Run with FM_ML env:
  /opt/miniconda3/envs/FM_ML/bin/python write_hash_input_ge8.py
"""
from pathlib import Path
from datetime import datetime
import numpy as np
import h5py

# ── Config ────────────────────────────────────────────────────────────────────
PROJECT   = Path('/Users/mczhang/Documents/GitHub/FM7')
DATA      = PROJECT / '02-data'
HASH_DIR  = PROJECT / '01-scripts' / 'HASH'
OUT_DIR   = DATA / 'G_FM'
OUT_DIR.mkdir(exist_ok=True)

NPOLMIN  = 15
CONF_MIN = 0.50
DEP_MIN  = 0.5
EH, EZ   = 0.3, 0.2
MAG      = 1.0
CHAN     = 'HHZ'
NET      = 'OO'
ONSET    = 'I'

STATION_META = {
    'AS1': (45.93356, -129.99920,  60),
    'AS2': (45.93377, -130.01410,  44),
    'CC1': (45.95468, -130.00890,  61),
    'EC1': (45.94958, -129.97970,  77),
    'EC2': (45.93967, -129.97380,  70),
    'EC3': (45.93607, -129.97850,  73),
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

def short_to_hash(s):
    return 'A' + s   # e.g. AS1 → AAS1

def deg_to_degmin(decimal_deg):
    d = abs(decimal_deg)
    deg = int(d)
    mins = 60.0 * (d - deg)
    return deg, mins

def load_loc_file(path):
    locs = {}
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 11: continue
            try:
                yr, mo, da = int(parts[0]), int(parts[1]), int(parts[2])
                hh, mn    = int(parts[3]), int(parts[4])
                ss        = float(parts[5])
                lat, lon, depth = float(parts[6]), float(parts[7]), float(parts[8])
                eid = int(parts[10])
                sec_i = int(ss)
                usec  = int(round((ss - sec_i) * 1e6))
                locs[eid] = (datetime(yr, mo, da, hh, mn, sec_i, usec), lon, lat, depth)
            except (ValueError, IndexError):
                continue
    return locs

# ── Load both years ────────────────────────────────────────────────────────────
years = [
    ('Y1', [f'{i:02d}A' for i in range(1, 15)] + ['15A']),
    ('Y2', [f'{i:02d}B' for i in range(1, 15)]),
]
OO_STAS = ['AS1', 'AS2', 'CC1', 'EC1', 'EC2', 'EC3', 'ID1']

all_events = []   # list of (eid, dt, lon, lat, depth, sta_data)

for year, TF_STAS in years:
    STATIONS = OO_STAS + TF_STAS
    loc_by_id = load_loc_file(DATA / 'A_all' / f'axial.{year}.mldd.loc.260317')
    print(f'{year}: {len(loc_by_id):,} locations')

    with h5py.File(DATA / 'E_Po' / f'dl_polarity_{year}_ge8_aug.h5', 'r') as fp:
        event_ids  = fp['event_ids'][:]
        polarity   = {s: fp[f'polarity/{s}'][:]   for s in STATIONS if f'polarity/{s}'   in fp}
        confidence = {s: fp[f'confidence/{s}'][:] for s in STATIONS if f'confidence/{s}' in fp}

    with h5py.File(DATA / 'A_ID' / f'waves_{year}_ge8.h5', 'r') as fw:
        amplitude  = {s: fw[f'amplitudes/{s}'][:] for s in STATIONS if f'amplitudes/{s}' in fw}

    n_written = n_noloc = n_npol = 0
    for i, eid in enumerate(event_ids):
        eid = int(eid)
        if eid not in loc_by_id:
            n_noloc += 1; continue
        dt, lon, lat, depth = loc_by_id[eid]
        depth = max(depth, DEP_MIN)

        sta_data = []
        for sta in STATIONS:
            if sta not in polarity: continue
            pol_val = int(polarity[sta][i])
            if pol_val == 0: continue
            conf = float(confidence[sta][i]) if sta in confidence else 0.0
            if np.isnan(conf) or conf < CONF_MIN: continue
            pol_char = 'U' if pol_val == 1 else 'D'
            sn = short_to_hash(sta)
            if sta in amplitude:
                amp = amplitude[sta][i]
                noise = max(float(amp[0]), 0.01) if not np.isnan(amp[0]) else 1.0
                p_amp = max(float(amp[2]), 0.01) if not np.isnan(amp[2]) else 1.0
                s_amp = max(float(amp[1]), 0.01) if not np.isnan(amp[1]) else p_amp
            else:
                noise, p_amp, s_amp = 1.0, 1.0, 1.0
            sta_data.append((sn, pol_char, noise, p_amp, s_amp))

        if len(sta_data) < NPOLMIN:
            n_npol += 1; continue

        all_events.append((eid, dt, lon, lat, depth, sta_data))
        n_written += 1

    print(f'  written={n_written:,}  no_loc={n_noloc:,}  <{NPOLMIN}pol={n_npol:,}')

print(f'\nTotal events: {len(all_events):,}')

# ── Write station file ─────────────────────────────────────────────────────────
all_stas = OO_STAS
for _, TF in years:
    all_stas = all_stas + TF
all_stas_unique = list(dict.fromkeys(all_stas))  # preserve order, deduplicate

sta_out = HASH_DIR / 'station_ge8.dat'
print(f'Writing {sta_out.name} ({len(all_stas_unique)} stations)...')
with open(sta_out, 'w') as f:
    for sta in all_stas_unique:
        lat, lon, elev = STATION_META[sta]
        sn = short_to_hash(sta)
        f.write(f'{sn:4s} {CHAN:3s}  {lat:.5f} {lon:10.5f} {elev:5d} {NET}\n')

# ── Write combined phase.dat and amp.dat ──────────────────────────────────────
pha_out = HASH_DIR / 'phase_ge8_1D.dat'
amp_out = HASH_DIR / 'amp_ge8_1D.dat'
eid_map = []   # (cluster_idx, event_id)

print(f'Writing {pha_out.name} and {amp_out.name}...')
with open(pha_out, 'w') as fpha, open(amp_out, 'w') as famp:
    for cluster_idx, (eid, dt, lon, lat, depth, sta_data) in enumerate(all_events, start=1):
        npicks = len(sta_data)
        eid_map.append((cluster_idx, eid))

        # Event header — HASH phase format (same as SKHASH north2 / hash3)
        yr = dt.year; mo = dt.month; da = dt.day
        hr = dt.hour; mn = dt.minute
        sec = dt.second + dt.microsecond / 1e6
        ilat, mlat = deg_to_degmin(lat)
        ilon, mlon = deg_to_degmin(abs(lon))
        ns = 'N' if lat >= 0 else 'S'
        ew = 'W' if lon < 0 else 'E'

        # HASH phase.dat fixed-width format:
        # YYYYMMDDHHMMSS.ssLatDegNSLatMin.mmLonDegEWLonMin.mmDepth Mag EH EV ... event_num
        # HASH phase.dat fixed-width: depth EH EV mag event_num (trailing space)
        header = (
            f'{yr:4d}{mo:2d}{da:2d}{hr:2d}{mn:2d}{sec:5.2f}'
            f'{ilat:2d}{ns}{mlat:5.2f}{ilon:3d}{ew}{mlon:5.2f}'
            f'{depth:5.2f}{EH:6.2f}{EZ:6.2f}{MAG:5.2f}'
            f'{cluster_idx:>16} \n'
        )
        fpha.write(header)
        for (sn, pol_char, _, _, _) in sta_data:
            fpha.write(f'{sn:4s} {NET:2s}  {CHAN:3s} {ONSET} {pol_char}\n')
        # Terminator: 70 spaces + event_num + trailing space
        fpha.write(f'{"":70}{cluster_idx} \n')

        # HASH amp.dat Fortran format: (a4,1x,a3,1x,a2,17x,f10.3,1x,f10.3,1x,f10.3,1x,f10.3)
        # columns: P_noise  S_noise  P_amp  S_amp
        famp.write(f'{cluster_idx}     {npicks}\n')
        for (sn, _, noise, p_amp, s_amp) in sta_data:
            famp.write(
                f'{sn:4s} {CHAN:3s} {NET:2s}{"":17}'
                f'{noise:10.3f} {noise:10.3f}'
                f' {p_amp:10.3f} {s_amp:10.3f}\n'
            )

print(f'  {len(all_events):,} events written')

# ── Save event_id mapping ──────────────────────────────────────────────────────
import csv
map_out = OUT_DIR / 'hash_ge8_1D_eid_map.csv'
with open(map_out, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['cluster_idx', 'event_id'])
    w.writerows(eid_map)
print(f'  → {map_out}')

# ── Write HASH input file ──────────────────────────────────────────────────────
inp_out = HASH_DIR / 'hash.input_ge8_1D'
print(f'Writing {inp_out.name}...')
with open(inp_out, 'w') as f:
    f.write(f"""\
station_ge8.dat
reverse.dat
Acor.dat
amp_ge8_1D.dat
phase_ge8_1D.dat
hashout_ge8_1D_1.dat
hashout_ge8_1D_2.dat
hashout_ge8_1D_3.dat
hashout_ge8_1D_4.dat
360
180
5
30
300
2
0.00
2
0
0.30000
25
45
0.750000
1
velmod_axial1D.dat
""")

print('\nDone. Run HASH:')
print(f'  cd {HASH_DIR}')
print(f'  ./hash_driver3 < hash.input_ge8_1D')
