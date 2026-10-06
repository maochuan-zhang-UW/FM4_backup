"""
write_skhash_input_ge8.py
=========================
Generate SKHASH input files for Y1 and Y2 ≥8-station events (CC-augmented).

Uses the CC-augmented polarity files:
  02-data/E_Po/dl_polarity_Y1_ge8_aug.h5
  02-data/E_Po/dl_polarity_Y2_ge8_aug.h5

Writes (for each year YEAR in {Y1, Y2}):
  SKHASH/SKHASH7/examples/hash3/IN/north2_YEAR_ge8.txt   — polarity (fpfile)
  SKHASH/SKHASH7/examples/hash3/IN/north3_YEAR_ge8.txt   — amplitude (ampfile)
  SKHASH/SKHASH7/examples/hash3/IN/scsn.stations_YEAR_ge8.txt
  SKHASH/SKHASH7/examples/hash3/IN/north3.statcor_YEAR_ge8.txt
  SKHASH/SKHASH7/examples/hash3/control_file_YEAR_ge8.txt
  SKHASH/SKHASH7/examples/hash3/OUT_YEAR_ge8/

Only events with ≥ NPOLMIN confident polarities after CC augmentation are written.

Run with FM_ML env:
  /opt/miniconda3/envs/FM_ML/bin/python write_skhash_input_ge8.py
"""

from pathlib import Path
from datetime import datetime

import numpy as np
import h5py

# ── Config ────────────────────────────────────────────────────────────────────
PROJECT    = Path('/Users/mczhang/Documents/GitHub/FM7')
DATA       = PROJECT / '02-data'
SKHASH_DIR = Path('/Users/mczhang/Documents/GitHub/SKHASH/SKHASH7')
HASH3_DIR  = SKHASH_DIR / 'examples' / 'hash3'

NPOLMIN   = 15      # minimum polarities to write event to SKHASH input
CONF_MIN  = 0.50    # minimum confidence to include a pick
DEP_MIN   = 0.5     # minimum depth (km)
EH        = 0.3     # horizontal uncertainty (km)
EZ        = 0.2     # vertical uncertainty (km)
MAG       = 1.0     # placeholder magnitude
CHAN      = 'HHZ'
NET       = 'OO'
ONSET     = 'I'

# ── Station metadata ──────────────────────────────────────────────────────────
# (lat, lon, elev_m)
STATION_META = {
    # OO stations (same both years)
    'AS1': (45.93356, -129.99920,  60),
    'AS2': (45.93377, -130.01410,  44),
    'CC1': (45.95468, -130.00890,  61),
    'EC1': (45.94958, -129.97970,  77),
    'EC2': (45.93967, -129.97380,  70),
    'EC3': (45.93607, -129.97850,  73),
    'ID1': (45.92573, -129.97800,  62),
    # Y1 2F (A-suffix)
    '01A': (46.01933, -130.00538,   0),
    '02A': (46.00063, -130.02926,  34),
    '03A': (45.99292, -129.99278, 102),
    '04A': (45.98874, -130.04840,  53),
    '05A': (45.97439, -129.97797,  73),
    '06A': (45.98294, -130.01413,   9),
    '07A': (45.96907, -130.02985,   8),
    '08A': (45.96415, -130.00450,  46),
    '09A': (45.97084, -130.06239, 103),
    '10A': (45.96034, -129.94999,  20),
    '11A': (45.94959, -130.03938, 158),
    '12A': (45.91537, -130.02178,  73),
    '13A': (45.90782, -129.97156,  24),
    '14A': (45.89997, -130.00757,   1),
    '15A': (45.91972, -129.93942,   0),
    # Y2 2F (B-suffix)
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
}

# ── Helpers ───────────────────────────────────────────────────────────────────
def short_to_skhash(s):
    return 'A' + s

def deg_to_degmin(decimal_deg):
    d = abs(decimal_deg)
    deg = int(d)
    mins = 60.0 * (d - deg)
    return deg, mins

def load_loc_file(path):
    """Return dict event_id -> (datetime, lon, lat, depth_km)."""
    locs = {}
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 11:
                continue
            try:
                yr, mo, da = int(parts[0]), int(parts[1]), int(parts[2])
                hh, mn    = int(parts[3]), int(parts[4])
                ss        = float(parts[5])
                lat       = float(parts[6])
                lon       = float(parts[7])
                depth     = float(parts[8])
                eid       = int(parts[10])
                sec_i     = int(ss)
                usec      = int(round((ss - sec_i) * 1e6))
                dt = datetime(yr, mo, da, hh, mn, sec_i, usec)
                locs[eid] = (dt, lon, lat, depth)
            except (ValueError, IndexError):
                continue
    return locs


def write_year(year):
    print(f'\n{"="*60}')
    print(f'Writing {year} ≥8 SKHASH inputs...')
    print(f'{"="*60}')

    # Paths
    tf_suffix = 'A' if year == 'Y1' else 'B'
    TF_STAS = [f'{i:02d}{tf_suffix}' for i in range(1, 15)] + ['15A']
    OO_STAS = ['AS1', 'AS2', 'CC1', 'EC1', 'EC2', 'EC3', 'ID1']
    STATIONS = OO_STAS + TF_STAS

    LOC_FILE  = DATA / 'A_all' / f'axial.{year}.mldd.loc.260317'
    POL_H5    = DATA / 'E_Po'  / f'dl_polarity_{year}_ge8_aug.h5'
    WAVES_H5  = DATA / 'A_ID'  / f'waves_{year}_ge8.h5'
    OUT_TAG   = f'{year}_ge8'
    OUT_PHA   = HASH3_DIR / 'IN' / f'north2_{OUT_TAG}.txt'
    OUT_AMP   = HASH3_DIR / 'IN' / f'north3_{OUT_TAG}.txt'
    OUT_STA   = HASH3_DIR / 'IN' / f'scsn.stations_{OUT_TAG}.txt'
    OUT_COR   = HASH3_DIR / 'IN' / f'north3.statcor_{OUT_TAG}.txt'
    OUT_CTL   = HASH3_DIR / f'control_file_{OUT_TAG}.txt'
    OUT_DIR   = HASH3_DIR / f'OUT_{OUT_TAG}'
    OUT_DIR.mkdir(exist_ok=True)

    # Load locations
    print(f'Loading locations from {LOC_FILE.name}...')
    loc_by_id = load_loc_file(LOC_FILE)
    print(f'  {len(loc_by_id):,} events')

    # Load polarities and amplitudes
    print(f'Loading polarities from {POL_H5.name}...')
    with h5py.File(POL_H5, 'r') as fp:
        event_ids  = fp['event_ids'][:]
        polarity   = {s: fp[f'polarity/{s}'][:]   for s in STATIONS if f'polarity/{s}'   in fp}
        confidence = {s: fp[f'confidence/{s}'][:] for s in STATIONS if f'confidence/{s}' in fp}

    print(f'Loading amplitudes from {WAVES_H5.name}...')
    with h5py.File(WAVES_H5, 'r') as fw:
        amplitude = {s: fw[f'amplitudes/{s}'][:] for s in STATIONS if f'amplitudes/{s}' in fw}

    N = len(event_ids)
    print(f'  {N:,} events, {len(polarity)} pol arrays, {len(amplitude)} amp arrays')

    # Write station file
    print(f'Writing {OUT_STA.name}...')
    with open(OUT_STA, 'w') as f:
        for sta in STATIONS:
            lat, lon, elev = STATION_META[sta]
            sn = short_to_skhash(sta)
            f.write(f'{sn} HHZ BIG CHUCKAWALLA MTNS             '
                    f'{lat:.5f} {lon:.5f} {elev:5d} 1997/09/19 3000/01/01 OO\n')

    # Write statcor file
    print(f'Writing {OUT_COR.name}...')
    with open(OUT_COR, 'w') as f:
        for sta in STATIONS:
            f.write(f'{short_to_skhash(sta)} HHZ OO  0.0000\n')

    # Write polarity and amplitude files
    print(f'Writing {OUT_PHA.name} and {OUT_AMP.name}...')
    n_written = 0
    n_skipped_noloc = 0
    n_skipped_npol  = 0
    total_picks = 0
    cluster_idx = 0

    with open(OUT_PHA, 'w') as fpha, open(OUT_AMP, 'w') as famp:
        for i, eid in enumerate(event_ids):
            eid = int(eid)
            if eid not in loc_by_id:
                n_skipped_noloc += 1
                continue
            dt, lon, lat, depth = loc_by_id[eid]
            depth = max(depth, DEP_MIN)

            # Collect valid picks
            sta_data = []
            for sta in STATIONS:
                if sta not in polarity:
                    continue
                pol_val = int(polarity[sta][i])
                if pol_val == 0:
                    continue
                conf = float(confidence[sta][i]) if sta in confidence else 0.0
                if np.isnan(conf) or conf < CONF_MIN:
                    continue

                pol_char = 'U' if pol_val == 1 else 'D'
                sn = short_to_skhash(sta)

                if sta in amplitude:
                    amp = amplitude[sta][i]
                    noise = max(float(amp[0]), 0.01) if not np.isnan(amp[0]) else 1.0
                    p_amp = max(float(amp[2]), 0.01) if not np.isnan(amp[2]) else 1.0
                    s_amp = max(float(amp[1]), 0.01) if not np.isnan(amp[1]) else p_amp
                else:
                    noise, p_amp, s_amp = 1.0, 1.0, 1.0

                sta_data.append((sn, pol_char, noise, p_amp, s_amp))

            if len(sta_data) < NPOLMIN:
                n_skipped_npol += 1
                continue

            cluster_idx += 1
            npicks = len(sta_data)

            # Polarity file: event header
            yr = dt.year; mo = dt.month; da = dt.day
            hr = dt.hour; mn = dt.minute
            sec = dt.second + dt.microsecond / 1e6
            ilat, mlat = deg_to_degmin(lat)
            ilon, mlon = deg_to_degmin(lon)
            eid_str = f'{cluster_idx:>22}'

            header = (
                f'{yr:4d}{mo:2d}{da:2d}{hr:2d}{mn:2d}{sec:5.2f}'
                f'{ilat:2d} {mlat:5.2f}{ilon:3d} {mlon:5.2f}{depth:5.2f}'
                f'     {npicks:3d}{"":41}'
                f'{EH:5.2f}  {EZ:5.2f}   {"":36}'
                f'{MAG:4.2f}{eid_str}\n'
            )
            fpha.write(header)
            for (sn, pol_char, _, _, _) in sta_data:
                fpha.write(f'{sn:4s} {NET:2s}  {CHAN:3s} {ONSET} {pol_char}\n')
            fpha.write(f'{"":56}{cluster_idx:16d}\n')

            # Amplitude file
            famp.write(f'{cluster_idx}         {npicks}\n')
            for (sn, _, noise, p_amp, s_amp) in sta_data:
                famp.write(
                    f'{sn:4s} {CHAN:3s} {NET:2s}'
                    f'  {noise:5.2f}   {noise:5.2f}'
                    f'  {noise:10.3f} {noise:10.3f}'
                    f' {p_amp:10.3f} {s_amp:10.3f}\n'
                )

            n_written  += 1
            total_picks += npicks

    print(f'  Written   : {n_written:,} events')
    print(f'  Skipped (no loc): {n_skipped_noloc:,}')
    print(f'  Skipped (<{NPOLMIN} pol): {n_skipped_npol:,}')
    print(f'  Total picks: {total_picks:,}  avg {total_picks/max(n_written,1):.1f}/event')

    # Write control file
    # Velocity model paths: one per line, use the '1'-suffixed filenames
    # that match the working control_file_Y2.txt (e.g. velmod_W11.txt)
    print(f'Writing {OUT_CTL.name}...')
    vmod_dir = 'examples/velocity_models'
    vel_files = '\n'.join([
        f'{vmod_dir}/velmod_W11.txt', f'{vmod_dir}/velmod_W21.txt',
        f'{vmod_dir}/velmod_E11.txt', f'{vmod_dir}/velmod_E21.txt',
        f'{vmod_dir}/velmod_E31.txt', f'{vmod_dir}/velmod_E41.txt',
        f'{vmod_dir}/velmod_S11.txt',
        f'{vmod_dir}/velmod_AXAS11.txt', f'{vmod_dir}/velmod_AXAS21.txt',
        f'{vmod_dir}/velmod_AXCC11.txt',
        f'{vmod_dir}/velmod_AXEC11.txt', f'{vmod_dir}/velmod_AXEC21.txt',
        f'{vmod_dir}/velmod_AXEC31.txt',
        f'{vmod_dir}/velmod_AXID11.txt',
        f'{vmod_dir}/velmod21.txt',
    ])
    with open(OUT_CTL, 'w') as f:
        f.write(f"""\
## SKHASH control file for {year} ≥8-station events (CC-augmented)

$input_format
hash3

$stfile
examples/hash3/IN/scsn.stations_{OUT_TAG}.txt

$plfile
examples/hash3/IN/scsn.reverse.txt

$corfile
examples/hash3/IN/north3.statcor_{OUT_TAG}.txt

$ampfile
examples/hash3/IN/north3_{OUT_TAG}.txt

$fpfile
examples/hash3/IN/north2_{OUT_TAG}.txt

$outfile1
examples/hash3/OUT_{OUT_TAG}/out.txt

$outfile2
examples/hash3/OUT_{OUT_TAG}/out2.txt

$npolmin
{NPOLMIN}

$dang
5

$nmc
30

$maxout
300

$ratmin
2

$badfrac
0.1

$qbadfrac
0.3

$delmax
25

$cangle
45

$prob_max
0.75

$vmodel_paths
{vel_files}
""")
    print(f'  Control file → {OUT_CTL}')
    print(f'\nTo run SKHASH for {year}:')
    print(f'  cd /Users/mczhang/Documents/GitHub/SKHASH/SKHASH7')
    print(f'  python SKHASH.py examples/hash3/control_file_{OUT_TAG}.txt')


# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    import sys
    years = sys.argv[1:] if len(sys.argv) > 1 else ['Y1', 'Y2']
    for year in years:
        if year not in ('Y1', 'Y2'):
            print(f'Unknown year: {year}  (use Y1 or Y2)')
            continue
        write_year(year)
    print('\nwrite_skhash_input_ge8.py complete.')
