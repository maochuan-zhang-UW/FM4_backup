"""
pipeline_config.py
==================
Central configuration for the FM7 Y1+Y2 focal mechanism pipeline.
All paths, parameters, and constants live here so every notebook imports
from one place. If you move the project to a new machine, edit only this file.
"""
from pathlib import Path

# ──────────────────────────────────────────────────────────────────────────────
# Root directories
# ──────────────────────────────────────────────────────────────────────────────
PROJECT    = Path('/Users/mczhang/Documents/GitHub/FM7')
DATA       = PROJECT / '02-data'
SCRIPTS    = PROJECT / '01-scripts'
PYTHON_DIR = SCRIPTS / 'python'

# External repos (siblings)
MSEED_ROOT = Path('/Users/mczhang/Documents/GitHub/FM4/01-scripts/datamseed')
FM6_DIR    = Path('/Users/mczhang/Documents/GitHub/FM6_RealTime')
SKHASH_DIR = Path('/Users/mczhang/Documents/GitHub/SKHASH/SKHASH7')
HASH_DIR   = SCRIPTS / 'HASH'

# ──────────────────────────────────────────────────────────────────────────────
# Input files (in 02-data/A_all/)
# ──────────────────────────────────────────────────────────────────────────────
LOC_Y1 = DATA / 'A_all' / 'axial.Y1.mldd.loc.260317'   # 6,727 events
LOC_Y2 = DATA / 'A_all' / 'axial.Y2.mldd.loc.260317'   # 48,244 events
PHA_Y1 = DATA / 'A_all' / 'axial.Y1.mldd.pha.260317'
PHA_Y2 = DATA / 'A_all' / 'axial.Y2.mldd.pha.260317'

# Date ranges (for reference)
Y1_START = '2022-09-01'   # Sep 1  2022
Y1_END   = '2023-08-31'   # Aug 31 2023
Y2_START = '2023-09-01'   # Sep 1  2023
Y2_END   = '2024-09-22'   # Sep 22 2024

# ──────────────────────────────────────────────────────────────────────────────
# Intermediate output files (02-data/) — one entry per notebook output
# ──────────────────────────────────────────────────────────────────────────────

# Notebook 01 — parse mldd catalog
CATALOG_PQ = DATA / 'A_all' / 'catalog_Y1Y2.parquet'
PHASES_PQ  = DATA / 'A_all' / 'phases_Y1Y2.parquet'

# Notebook 02 — cut waveforms
WAVES_H5   = DATA / 'A_ID'  / 'waves_Y1Y2.h5'
# HDF5 groups: /waveforms/<STA>  shape (N,264) float32
#              /amplitudes/<STA> shape (N,3)   float32 [noise, S_amp, P_amp]
#              /event_ids        shape (N,)     int64

# Notebook 03 — cross-correlation
CC_H5      = DATA / 'B_CC'  / 'cc_Y1Y2.h5'
# HDF5 groups: /<STA>  shape (M,3) [ID1, ID2, cc_value]  only pairs > CC_THRESHOLD

# Notebook 04 — SVD clustering
SVD_PQ     = DATA / 'C_SVD' / 'svd_Y1Y2.parquet'
# columns: event_id, AS1_score, AS2_score, ..., 14A_score

# Notebook 05 — select 20 high-SNR events per station
SELECTED_PQ = DATA / 'D_man' / 'selected_20perSta_Y1Y2.parquet'
# columns: station, event_id, snr, svd_score, rank

# Notebook 06 — DL polarity
DL_SELECTED_MAT = DATA / 'D_man' / 'dl_selected_Y1Y2.mat'   # 20 per station → user review
DL_ALL_MAT      = DATA / 'E_Po'  / 'dl_polarity_Y1Y2.mat'   # all 54,971 events

# Notebook 07 — polarity combine
POLARITY_PQ = DATA / 'E_Po' / 'polarity_combined_Y1Y2.parquet'
# columns: event_id, Po_AS1...Po_14A (±1/2/3), NSP_AS1...NSP_14A (3 vals), PoALL

# Notebook 08 — hierarchical clustering
CLUSTERING_PQ  = DATA / 'F_Cl' / 'clustering_Y1Y2.parquet'
REMOVED_IDS_PQ = DATA / 'F_Cl' / 'removed_wave_IDs_Y1Y2.parquet'

# Notebook 09a — write SKHASH input
SKHASH_IN_DAT  = SKHASH_DIR / 'examples/hash3/IN/Axial_22OBSs_Y1Y2.dat'
SKHASH_OUT_TXT = SKHASH_DIR / 'examples/hash3/OUT/out.txt'
SKHASH_SUMMARY = DATA / 'G_FM' / 'skhash_input_summary_Y1Y2.parquet'

# Notebook 09b — write HASH input
HASH_PHASE_DAT  = HASH_DIR / 'phase_Y1Y2.dat'
HASH_AMP_DAT    = HASH_DIR / 'amp_Y1Y2.dat'
HASH_STATION_DAT = HASH_DIR / 'station.dat'       # already exists, reuse
HASH_OUT_DAT    = HASH_DIR / 'hashout_Y1Y2.dat'   # written by hash_driver3

# Notebook 10a — parse SKHASH output
FM_SKHASH_PQ  = DATA / 'G_FM' / 'focal_mechanisms_SKHASH_Y1Y2.parquet'
FM_SKHASH_MAT = DATA / 'G_FM' / 'G_3D_SKHASH_Y1Y2.mat'

# Notebook 10b — parse HASH output
FM_HASH_PQ  = DATA / 'G_FM' / 'focal_mechanisms_HASH_Y1Y2.parquet'
FM_HASH_MAT = DATA / 'G_FM' / 'G_3D_HASH_Y1Y2.mat'

# Notebook 11 — polarity comparison (DL vs pipeline)
POLARITY_COMP_PQ = DATA / 'G_FM' / 'polarity_comparison_Y1Y2.parquet'

# Notebook 12 — FM comparison (SKHASH vs HASH)
FM_COMP_PQ = DATA / 'G_FM' / 'fm_comparison_SKHASH_HASH_Y1Y2.parquet'

# Pipeline state (for recovery after crash)
PROGRESS_JSON = DATA / 'pipeline_progress.json'

# Email attachments staging folder
EMAIL_PLOTS_DIR = DATA / 'email_plots'

# ──────────────────────────────────────────────────────────────────────────────
# Station definitions
# ──────────────────────────────────────────────────────────────────────────────
STATIONS_LEGACY = ['AS1', 'AS2', 'CC1', 'EC1', 'EC2', 'EC3', 'ID1']
STATIONS_MAIN   = ['01A', '02A', '03A', '04A', '05A', '06A', '07A',
                   '08A', '09A', '10A', '11A', '12A', '13A', '14A']
STATIONS_ALL    = STATIONS_LEGACY + STATIONS_MAIN  # 21 total

# (lat, lon, elev_m)
STATION_COORDS = {
    'AS1': (45.93356, -129.99920,  59.7),
    'AS2': (45.93377, -130.01410,  44.3),
    'CC1': (45.95468, -130.00890,  60.7),
    'EC1': (45.94958, -129.97970,  68.5),
    'EC2': (45.93969, -129.97380,  56.3),
    'EC3': (45.93607, -129.97850,  64.8),
    'ID1': (45.92570, -129.97800,  75.0),
    '01A': (46.01933, -130.00538,   0.0),
    '02A': (46.00063, -130.02226,  34.0),
    '03A': (45.99292, -129.99278,   0.0),
    '04A': (45.98874, -130.04840,   0.0),
    '05A': (45.97440, -129.97797,   0.0),
    '06A': (45.98294, -130.01413,   0.0),
    '07A': (45.96907, -130.02985,   0.0),
    '08A': (45.96415, -130.00450,   0.0),
    '09A': (45.97084, -130.06239,   0.0),
    '10A': (45.96034, -129.94999,   0.0),
    '11A': (45.94959, -130.03398,   0.0),
    '12A': (45.91537, -130.02177,   0.0),
    '13A': (45.90782, -129.97156,   0.0),
    '14A': (45.89997, -130.00757,   0.0),
}

# mseed station ID mapping: our short name → mseed network.station
# Discovered by inspecting the first mseed file in notebook 02.
# Legacy stations: network OO, station AXAS1 etc.
# Main array:      network 2F, station AX01A etc.
MSEED_STA_MAP = {
    'AS1': ('OO', 'AXAS1'),  'AS2': ('OO', 'AXAS2'),
    'CC1': ('OO', 'AXCC1'),  'EC1': ('OO', 'AXEC1'),
    'EC2': ('OO', 'AXEC2'),  'EC3': ('OO', 'AXEC3'),
    'ID1': ('OO', 'AXID1'),
    '01A': ('2F', 'AX01A'),  '02A': ('2F', 'AX02A'),
    '03A': ('2F', 'AX03A'),  '04A': ('2F', 'AX04A'),
    '05A': ('2F', 'AX05A'),  '06A': ('2F', 'AX06A'),
    '07A': ('2F', 'AX07A'),  '08A': ('2F', 'AX08A'),
    '09A': ('2F', 'AX09A'),  '10A': ('2F', 'AX10A'),
    '11A': ('2F', 'AX11A'),  '12A': ('2F', 'AX12A'),
    '13A': ('2F', 'AX13A'),  '14A': ('2F', 'AX14A'),
}
HHZ_CHANNEL = 'HHZ'   # vertical component channel used for P polarity

# ──────────────────────────────────────────────────────────────────────────────
# Waveform parameters
# ──────────────────────────────────────────────────────────────────────────────
FS        = 200            # original sample rate (Hz)
WIN_START = -0.32          # seconds before P pick
WIN_END   =  1.00          # seconds after P pick
N_SAMPLES = int((WIN_END - WIN_START) * FS)   # = 264 samples
BANDPASS  = (4.0, 50.0)   # Hz, same as FM4 pipeline

DL_FS = 100      # DL model sample rate (Hz) — resample from 200 to 100
DL_N  = 64       # DL snippet length: 64 samples at 100 Hz = [-0.32, +0.31]s

# Amplitude windows (relative to P pick, seconds) — from D_SP_wave.m
NOISE_WIN = (-0.70, -0.10)
P_WIN     = (-0.05,  0.25)
S_WIN     = (-0.10,  0.60)   # relative to S pick

# ──────────────────────────────────────────────────────────────────────────────
# Cross-correlation parameters (B_CC_all.m)
# ──────────────────────────────────────────────────────────────────────────────
CC_THRESHOLD = 0.70
CC_LAG_AHEAD = 0.10   # seconds before P to start CC window
CC_WIN_SHORT = 0.20   # seconds
CC_WIN_LONG  = 0.40   # seconds

# ──────────────────────────────────────────────────────────────────────────────
# Polarity / clustering parameters (F_Cl_Real.m)
# ──────────────────────────────────────────────────────────────────────────────
MIN_PO_NUM = 5     # minimum polarity picks per event
MIN_SP_NUM = 5     # minimum S/P amplitude ratios per event
A_VALUE    = 0.2   # hierarchical clustering parameter

# ──────────────────────────────────────────────────────────────────────────────
# DL model
# ──────────────────────────────────────────────────────────────────────────────
DL_MODEL   = FM6_DIR / '02-data' / 'PolarPicker_unified_TMSF_001.keras'
PREDICT_PY = FM6_DIR / '01-scripts' / 'subcode' / 'predict_polarity.py'
PYTHON_EXE = '/opt/miniconda3/envs/FM_RT/bin/python'
DL_CONF_MIN = 0.60   # flag picks below this confidence as uncertain

# ──────────────────────────────────────────────────────────────────────────────
# HASH parameters
# ──────────────────────────────────────────────────────────────────────────────
HASH_BINARY      = HASH_DIR / 'hash_driver3'
HASH_INPUT_FILE  = HASH_DIR / 'hash.input'
HASH_VELMOD_DIR  = HASH_DIR   # velocity model .dat files are in HASH_DIR

# ──────────────────────────────────────────────────────────────────────────────
# Email
# ──────────────────────────────────────────────────────────────────────────────
EMAIL_TO       = 'mczhang26@gmail.com'
EMAIL_INTERVAL = 3600   # seconds between hourly updates
