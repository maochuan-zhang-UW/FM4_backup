"""
dl_polarity_ge8.py
==================
Run PolarPicker DL model on waveforms from waves_Y1_ge8.h5 or waves_Y2_ge8.h5.

Usage:
    /opt/miniconda3/envs/FM_RT/bin/python dl_polarity_ge8.py Y1
    /opt/miniconda3/envs/FM_RT/bin/python dl_polarity_ge8.py Y2

Output:
    02-data/E_Po/dl_polarity_Y1_ge8.h5   or   dl_polarity_Y2_ge8.h5
    /polarity/<STA>   int8  +1=Up, -1=Down, 0=no waveform
    /confidence/<STA> float32  0.5-1.0 (NaN where no waveform)
    /event_ids        int64

Environment: FM_RT (keras 3.12.1, tensorflow)
"""
import sys, time, warnings
from pathlib import Path

import numpy as np
import h5py

warnings.filterwarnings('ignore')

# ── Args ──────────────────────────────────────────────────────────────────────
if len(sys.argv) < 2 or sys.argv[1] not in ('Y1', 'Y2'):
    print('Usage: python dl_polarity_ge8.py Y1|Y2')
    sys.exit(1)
YEAR = sys.argv[1]

# ── Paths ─────────────────────────────────────────────────────────────────────
PROJECT = Path('/Users/mczhang/Documents/GitHub/FM7')
DATA    = PROJECT / '02-data'
WAVES_H5  = DATA / 'A_ID' / f'waves_{YEAR}_ge8.h5'
OUT_H5    = DATA / 'E_Po' / f'dl_polarity_{YEAR}_ge8.h5'
DL_MODEL  = Path('/Users/mczhang/Documents/GitHub/FM6_RealTime/02-data/PolarPicker_unified_TMSF_001.keras')
OUT_H5.parent.mkdir(parents=True, exist_ok=True)

# ── Stations ──────────────────────────────────────────────────────────────────
OO_STAS = ['AS1', 'AS2', 'CC1', 'EC1', 'EC2', 'EC3', 'ID1']
if YEAR == 'Y2':
    TF_STAS = [f'{i:02d}B' for i in range(1, 15)] + ['15A']
else:
    TF_STAS = [f'{i:02d}A' for i in range(1, 15)] + ['15A']
STATIONS = OO_STAS + TF_STAS   # 22 total

# ── DL params ─────────────────────────────────────────────────────────────────
DL_FS   = 100     # model input sample rate (Hz)
DL_N    = 64      # model input length (samples)
SRC_FS  = 200     # waveform sample rate in H5
BATCH   = 512
DL_CONF_MIN = 0.7

print(f'=== DL Polarity — {YEAR} ≥8 stations ===')
print(f'Waves : {WAVES_H5}')
print(f'Output: {OUT_H5}')
print(f'Model : {DL_MODEL}')

# ── Load model ────────────────────────────────────────────────────────────────
import keras
print(f'\nkeras {keras.__version__}')
print('Loading model...')
t0 = time.time()
model = keras.models.load_model(str(DL_MODEL))
print(f'  Loaded in {time.time()-t0:.1f}s')
print(f'  Input  : {model.input_shape}')
print(f'  Outputs: {[o.shape for o in model.outputs]}')

# ── Load waveforms ────────────────────────────────────────────────────────────
print('\nLoading waveforms...')
t0 = time.time()
with h5py.File(WAVES_H5, 'r') as hf:
    event_ids = hf['event_ids'][:]
    waves = {}
    for sta in STATIONS:
        key = f'waveforms/{sta}'
        if key in hf:
            waves[sta] = hf[key][:]   # (N, 264) float32
N = len(event_ids)
print(f'  Events  : {N:,}  ({time.time()-t0:.1f}s)')
print(f'  Stations: {len(waves)} waveform arrays found')
for sta in STATIONS:
    if sta in waves:
        ok = int(np.sum(~np.isnan(waves[sta][:, 0])))
        pct = ok / N * 100
        print(f'    {sta}: {ok:>6,}/{N}  ({pct:.0f}%)')

# ── Preprocessing ─────────────────────────────────────────────────────────────
def prepare_input(wave_264):
    """(N,264)@200Hz -> (N,64,1)@100Hz, normalised; also returns valid mask."""
    valid = ~np.isnan(wave_264[:, 0])
    w100  = wave_264[:, ::2]           # downsample 200→100 Hz
    w64   = w100[:, :DL_N].copy()      # first 64 samples
    peak  = np.max(np.abs(w64), axis=1, keepdims=True)
    peak[peak == 0] = 1.0
    w64   = w64 / peak
    w64[~valid] = 0.0
    return w64[:, :, np.newaxis].astype(np.float32), valid

# ── Predictions ───────────────────────────────────────────────────────────────
print(f'\nRunning DL predictions (batch={BATCH})...')
t0 = time.time()
pol_out  = {}
conf_out = {}

for i, sta in enumerate(STATIONS):
    t1 = time.time()
    if sta not in waves:
        pol_out[sta]  = np.zeros(N, dtype=np.int8)
        conf_out[sta] = np.full(N, np.nan, dtype=np.float32)
        print(f'  [{i+1:2d}/{len(STATIONS)}] {sta}: no waveform array — skipped')
        continue

    X, valid = prepare_input(waves[sta])
    _, class_probs = model.predict(X, batch_size=BATCH, verbose=0)
    # class_probs: (N,2)  [P(Down), P(Up)]

    polarity   = np.where(np.argmax(class_probs, axis=1) == 1, 1, -1).astype(np.int8)
    confidence = np.max(class_probs, axis=1).astype(np.float32)
    polarity[~valid]   = 0
    confidence[~valid] = np.nan

    pol_out[sta]  = polarity
    conf_out[sta] = confidence

    n_picks = int(valid.sum())
    n_hi    = int(((polarity != 0) & (confidence >= DL_CONF_MIN)).sum())
    n_up    = int((polarity == 1).sum())
    n_dn    = int((polarity == -1).sum())
    ratio   = n_up / max(n_dn, 1)
    print(f'  [{i+1:2d}/{len(STATIONS)}] {sta}: {n_picks:>6,} picks  '
          f'hi-conf={n_hi:>6,}  Up={n_up:>5,}  Dn={n_dn:>5,}  '
          f'U/D={ratio:.2f}  ({time.time()-t1:.1f}s)')

elapsed = time.time() - t0
total_picks = sum((pol_out[s] != 0).sum() for s in STATIONS)
high_conf   = sum(((pol_out[s] != 0) & (conf_out[s] >= DL_CONF_MIN)).sum() for s in STATIONS)
print(f'\nTotal elapsed: {elapsed:.0f}s')
print(f'Total picks  : {total_picks:,}')
print(f'High conf    : {high_conf:,}  ({high_conf/max(total_picks,1)*100:.1f}%)')

# ── Save ──────────────────────────────────────────────────────────────────────
print(f'\nSaving → {OUT_H5}')
with h5py.File(OUT_H5, 'w') as hf:
    hf.create_dataset('event_ids', data=event_ids)
    for sta in STATIONS:
        hf.create_dataset(f'polarity/{sta}',   data=pol_out[sta])
        hf.create_dataset(f'confidence/{sta}', data=conf_out[sta])
sz = OUT_H5.stat().st_size / 1e6
print(f'Done. {OUT_H5.name}  {sz:.1f} MB')

# ── Per-event pick count distribution ─────────────────────────────────────────
picks_per_event = np.zeros(N, dtype=np.int32)
hi_per_event    = np.zeros(N, dtype=np.int32)
for sta in STATIONS:
    picks_per_event += (pol_out[sta] != 0).astype(np.int32)
    hi_per_event    += ((pol_out[sta] != 0) & (conf_out[sta] >= DL_CONF_MIN)).astype(np.int32)

print('\n=== Per-event polarity count distribution ===')
for threshold in [8, 10, 12, 15, 18, 20]:
    n_ge = int((picks_per_event >= threshold).sum())
    print(f'  ≥{threshold:2d} picks: {n_ge:>6,}  ({n_ge/N*100:.1f}%)')
print()
for threshold in [8, 10, 12, 15, 18, 20]:
    n_ge = int((hi_per_event >= threshold).sum())
    print(f'  ≥{threshold:2d} hi-conf: {n_ge:>6,}  ({n_ge/N*100:.1f}%)')
