"""
write_skhash_input_Y2_ge18.py
==============================
Generate SKHASH input files for Y2 ≥18-station events.

Format mirrors the MATLAB pipeline (J_Write_run_SKHASH_22OBS.m):
  north2.phase.txt  — HASH3/NCSN fixed-width polarity file (fpfile)
  north3.amp.txt    — HASH3 amplitude file (ampfile)

Inputs:
  02-data/A_all/axial.Y2.mldd.loc.260317  — event locations
  02-data/E_Po/dl_polarity_Y2_ge18.h5     — DL polarity + confidence
  02-data/A_ID/waves_Y2_ge18.h5           — amplitudes [noise, S_amp, P_amp]

Writes (relative to SKHASH/SKHASH7/examples/hash3/):
  IN/north2_Y2_ge18.txt       — polarity file (fpfile)
  IN/north3_Y2_ge18.txt       — amplitude file (ampfile)
  IN/scsn.stations_Y2.txt     — station list (OO + Y2 B-stations)
  IN/north3.statcor_Y2.txt    — zero station corrections
  control_file_Y2.txt         — SKHASH control file
  OUT_Y2/                     — output directory

Run with FM_ML env (needs h5py, numpy):
  /opt/miniconda3/envs/FM_ML/bin/python write_skhash_input_Y2_ge18.py
"""

from pathlib import Path
from datetime import datetime

import numpy as np
import h5py

# ── Paths ──────────────────────────────────────────────────────────────────
PROJECT    = Path('/Users/mczhang/Documents/GitHub/FM7')
DATA       = PROJECT / '02-data'
SKHASH_DIR = Path('/Users/mczhang/Documents/GitHub/SKHASH/SKHASH7')
HASH3_DIR  = SKHASH_DIR / 'examples' / 'hash3'

LOC_FILE   = DATA / 'A_all' / 'axial.Y2.mldd.loc.260317'
POL_H5     = DATA / 'E_Po'  / 'dl_polarity_Y2_ge18.h5'
WAVES_H5   = DATA / 'A_ID'  / 'waves_Y2_ge18.h5'

OUT_PHA    = HASH3_DIR / 'IN' / 'north2_Y2_ge18.txt'
OUT_AMP    = HASH3_DIR / 'IN' / 'north3_Y2_ge18.txt'
OUT_STA    = HASH3_DIR / 'IN' / 'scsn.stations_Y2.txt'
OUT_COR    = HASH3_DIR / 'IN' / 'north3.statcor_Y2.txt'
OUT_CTL    = HASH3_DIR / 'control_file_Y2.txt'
OUT_DIR    = HASH3_DIR / 'OUT_Y2'

# ── Station definitions ────────────────────────────────────────────────────
STATIONS_Y2 = [
    'AS1', 'AS2', 'CC1', 'EC1', 'EC2', 'EC3', 'ID1',
    '01B', '02B', '03B', '04B', '05B', '06B', '07B',
    '08B', '09B', '10B', '11B', '12B', '13B', '14B', '15A',
]

# Coordinates from axial.Y2.sta
STATION_META = {
    'AS1': (45.93356, -129.99920,  57),
    'AS2': (45.93377, -130.01410,  42),
    'CC1': (45.95468, -130.00890,  58),
    'EC1': (45.94958, -129.97970,  74),
    'EC2': (45.93967, -129.97380,  67),
    'EC3': (45.93607, -129.97850,  70),
    'ID1': (45.92573, -129.97800,  59),
    '01B': (46.01935, -130.00526,   2),
    '02B': (46.00111, -130.02930,  32),
    '03B': (45.99250, -129.99397, 107),
    '04B': (45.98834, -130.04896,  55),
    '05B': (45.97327, -129.97728,  74),
    '06B': (45.98178, -130.01423,  12),
    '07B': (45.96950, -130.02963,  11),
    '08B': (45.96336, -130.00441,  57),
    '09B': (45.97194, -130.05972, 111),
    '10B': (45.95867, -129.94982,  23),
    '11B': (45.94851, -130.03801, 160),
    '12B': (45.91473, -130.02173,  78),
    '13B': (45.90472, -129.97240,  20),
    '14B': (45.89913, -130.00612,   1),
    '15A': (45.91972, -129.93942,   0),
}

CHAN    = 'HHZ'
NET    = 'OO'
ONSET  = 'I'
EH     = 0.3   # horizontal uncertainty (km) — same as MATLAB
EZ     = 0.2   # vertical uncertainty (km)
MAG    = 1.0   # placeholder magnitude
DEP_MIN = 0.5  # minimum depth (km)
CONF_MIN = 0.50  # minimum DL confidence to include a pick


def short_to_skhash(s):
    """Short name → 4-char SKHASH station name: AS1→AAS1, 01B→A01B."""
    return 'A' + s


def dt_to_matlab(dt):
    """Python datetime → MATLAB datenum."""
    return (dt.toordinal() + 366
            + dt.hour / 24.0
            + dt.minute / 1440.0
            + (dt.second + dt.microsecond / 1e6) / 86400.0)


def deg_to_degmin(decimal_deg):
    """Decimal degrees → (integer_degrees, decimal_minutes)."""
    d = abs(decimal_deg)
    deg = int(d)
    mins = 60.0 * (d - deg)
    return deg, mins


# ── Load event locations ───────────────────────────────────────────────────
print('Loading event locations...')
loc_by_id = {}  # event_id (int) → (datetime, lon, lat, depth)
with open(LOC_FILE) as f:
    for line in f:
        parts = line.split()
        if len(parts) < 11:
            continue
        yr, mo, da, hh, mm = int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4])
        ss = float(parts[5])
        lat, lon, depth = float(parts[6]), float(parts[7]), float(parts[8])
        eid = int(parts[10])
        sec_i = int(ss)
        usec  = int(round((ss - sec_i) * 1e6))
        dt = datetime(yr, mo, da, hh, mm, sec_i, usec)
        loc_by_id[eid] = (dt, lon, lat, depth)

print(f'  Loaded {len(loc_by_id)} event locations')

# ── Load polarities and amplitudes ─────────────────────────────────────────
print('Loading polarities and amplitudes...')
with h5py.File(POL_H5, 'r') as fp, h5py.File(WAVES_H5, 'r') as fw:
    event_ids = fw['event_ids'][:]
    N = len(event_ids)
    polarity   = {s: fp['polarity'][s][:]   for s in STATIONS_Y2 if s in fp['polarity']}
    confidence = {s: fp['confidence'][s][:] for s in STATIONS_Y2 if s in fp['confidence']}
    amplitude  = {s: fw['amplitudes'][s][:] for s in STATIONS_Y2 if s in fw['amplitudes']}

print(f'  {N} events, {len(polarity)} stations with polarity data')

# ── Write station file ─────────────────────────────────────────────────────
print(f'Writing station file → {OUT_STA}')
with open(OUT_STA, 'w') as f:
    for sta in STATIONS_Y2:
        lat, lon, elev = STATION_META[sta]
        skhash_name = short_to_skhash(sta)
        f.write(f'{skhash_name} HHZ BIG CHUCKAWALLA MTNS             '
                f'{lat:.5f} {lon:.5f} {elev:5d} 1997/09/19 3000/01/01 OO\n')

# ── Write statcor file ─────────────────────────────────────────────────────
print(f'Writing statcor file → {OUT_COR}')
with open(OUT_COR, 'w') as f:
    for sta in STATIONS_Y2:
        skhash_name = short_to_skhash(sta)
        f.write(f'{skhash_name} HHZ OO  0.0000\n')

# ── Write polarity (north2) and amplitude (north3) files ───────────────────
print(f'Writing polarity file → {OUT_PHA}')
print(f'Writing amplitude file → {OUT_AMP}')

n_written  = 0
n_skipped  = 0
total_picks = 0

with open(OUT_PHA, 'w') as fpha, open(OUT_AMP, 'w') as famp:
    for i, eid in enumerate(event_ids):
        eid = int(eid)
        cluster_idx = i + 1  # 1-based

        if eid not in loc_by_id:
            n_skipped += 1
            continue
        dt, lon, lat, depth = loc_by_id[eid]
        depth = max(depth, DEP_MIN)

        # Collect valid station data
        sta_data = []  # [(skhash_name, pol_char, noise_amp, p_amp, s_amp)]
        for sta in STATIONS_Y2:
            if sta not in polarity:
                continue
            pol_val = int(polarity[sta][i])
            if pol_val == 0:
                continue
            conf = float(confidence[sta][i])
            if np.isnan(conf) or conf < CONF_MIN:
                continue

            pol_char = 'U' if pol_val == 1 else 'D'
            skhash_name = short_to_skhash(sta)

            if sta in amplitude:
                amp = amplitude[sta][i]
                noise = max(float(amp[0]), 0.01)
                p_amp = max(float(amp[2]), 0.01)
                s_amp = max(float(amp[1]), 0.01)
                if np.isnan(noise): noise = 1.0
                if np.isnan(p_amp): p_amp = 1.0
                if np.isnan(s_amp): s_amp = p_amp
            else:
                noise, p_amp, s_amp = 1.0, 1.0, 1.0

            sta_data.append((skhash_name, pol_char, noise, p_amp, s_amp))

        if not sta_data:
            n_skipped += 1
            continue

        # ── Polarity file: event header ─────────────────────────────────
        # Format mirrors MATLAB J_Write_run_SKHASH_22OBS.m line 72:
        # '%4i%2i%2i%2i%2i%5.2f%2i %5.2f%3i %5.2f%5.2f     %3i%41s%5.2f  %5.2f   %36s%4.2f%22s\n'
        yr = dt.year; mo = dt.month; da = dt.day
        hr = dt.hour; mn = dt.minute
        sec = dt.second + dt.microsecond / 1e6
        ilat, mlat = deg_to_degmin(lat)
        ilon, mlon = deg_to_degmin(lon)
        npicks = len(sta_data)

        # event_id at position 143-164 (22 chars, right-justified)
        # SKHASH reads line[149:165].strip() → last 16 chars of this 22-char field
        # right-justify cluster_idx in 22 chars: str(1182) = 4 chars → 18 spaces + "1182"
        eid_str = f'{cluster_idx:>22}'

        header = (
            f'{yr:4d}{mo:2d}{da:2d}{hr:2d}{mn:2d}{sec:5.2f}'
            f'{ilat:2d} {mlat:5.2f}{ilon:3d} {mlon:5.2f}{depth:5.2f}'
            f'     {npicks:3d}{"":41}'
            f'{EH:5.2f}  {EZ:5.2f}   {"":36}'
            f'{MAG:4.2f}{eid_str}\n'
        )
        fpha.write(header)

        # Station polarity lines: '%4s %2s  %3s %c %c\n'
        for (skhash_name, pol_char, _, _, _) in sta_data:
            fpha.write(f'{skhash_name:4s} {NET:2s}  {CHAN:3s} {ONSET} {pol_char}\n')

        # Terminator line: 56 spaces + cluster_idx right-justified in 16 chars
        fpha.write(f'{"":56}{cluster_idx:16d}\n')

        # ── Amplitude file ──────────────────────────────────────────────
        # Event header: '{cluster_idx}         {npicks}'
        famp.write(f'{cluster_idx}         {npicks}\n')

        # Station amplitude lines (MATLAB format line 90):
        # '%4s %3s %2s  %5.2f   %5.2f  %10.3f %10.3f %10.3f %10.3f\n'
        # values: sta chan net  Noip Noip Noip Nois Pamp Samp
        # in_sp.py reads: noise_p=col5, noise_s=col6, amp_p=col7, amp_s=col8
        for (skhash_name, _, noise, p_amp, s_amp) in sta_data:
            famp.write(
                f'{skhash_name:4s} {CHAN:3s} {NET:2s}'
                f'  {noise:5.2f}   {noise:5.2f}'
                f'  {noise:10.3f} {noise:10.3f}'
                f' {p_amp:10.3f} {s_amp:10.3f}\n'
            )

        n_written  += 1
        total_picks += npicks

print(f'  Written {n_written} events ({n_skipped} skipped), {total_picks} picks total')
print(f'  Average picks/event: {total_picks/max(n_written,1):.1f}')

# ── Create output directory ────────────────────────────────────────────────
OUT_DIR.mkdir(exist_ok=True)
(OUT_DIR / 'plots').mkdir(exist_ok=True)

# ── Write control file ─────────────────────────────────────────────────────
print(f'Writing control file → {OUT_CTL}')
with open(OUT_CTL, 'w') as f:
    f.write(f"""\
## Control file for Y2 ≥18-station events (Sep 2023–Sep 2024)

$input_format  # format of input files
hash3

$stfile        # station list filepath
examples/hash3/IN/scsn.stations_Y2.txt

$plfile        # station polarity reversal filepath
examples/hash3/IN/scsn.reverse.txt

$corfile       # station correction filepath
examples/hash3/IN/north3.statcor_Y2.txt

$ampfile       # amplitude input filepath
examples/hash3/IN/north3_Y2_ge18.txt

$fpfile        # P-polarity input filepath
examples/hash3/IN/north2_Y2_ge18.txt

$outfile1      # focal mechanisms output filepath
examples/hash3/OUT_Y2/out.txt

$outfile2      # acceptable plane output filepath
examples/hash3/OUT_Y2/out2.txt

$outfile_pol_agree
examples/hash3/OUT_Y2/out3_po.txt

$outfile_sp_agree
examples/hash3/OUT_Y2/out4_sp.txt

$outfile_pol_info
examples/hash3/OUT_Y2/out3_po2.txt

$outfolder_plots
examples/hash3/OUT_Y2/plots/

$vmodel_paths  # velocity models
examples/velocity_models/velmod_W11.txt
examples/velocity_models/velmod_S11.txt
examples/velocity_models/velmod_E21.txt
examples/velocity_models/velmod_E31.txt
examples/velocity_models/velmod_AXAS11.txt
examples/velocity_models/velmod_AXAS21.txt
examples/velocity_models/velmod_AXCC11.txt
examples/velocity_models/velmod_AXEC11.txt
examples/velocity_models/velmod_AXEC21.txt
examples/velocity_models/velmod21.txt
examples/velocity_models/velmod_W21.txt
examples/velocity_models/velmod_E41.txt
examples/velocity_models/velmod_E11.txt
examples/velocity_models/velmod_AXEC31.txt
examples/velocity_models/velmod_AXID11.txt

$npolmin       # minimum number of polarity data
15

$dang          # minimum grid spacing (degrees)
5

$nmc           # number of trials
30

$maxout        # max num of acceptable focal mech. outputs
300

$ratmin        # minimum allowed signal to noise ratio
2

$badfrac       # fraction polarities assumed bad
0.1

$qbadfrac      # assumed noise in amplitude ratios, log10
0.3

$delmax        # maximum allowed source-receiver distance in km
25

$cangle        # angle for computing mechanisms probability
45

$prob_max      # probability threshold for multiples
0.75
""")

print()
print('All SKHASH input files written.')
print()
print('To run SKHASH:')
print(f'  cd {SKHASH_DIR}')
print(f'  python SKHASH.py examples/hash3/control_file_Y2.txt')
print()
print(f'Output: {OUT_DIR}/out.txt')
