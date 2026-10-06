"""
E_polarcap_compare.py
======================
Run PolarCAP (Karim et al. 2023, https://doi.org/10.1016/j.aiig.2022.08.001;
model + weights from https://github.com/srivastavaresearchgroup/SAIPy) on the
FM7 waveform sets and compare its polarity picks, per station, against the
existing PolarPicker DL polarity already computed for the same events.

Input format is identical for both models: 64 samples @ 100 Hz, P pick at
index 32, peak-normalised, shape (N, 64, 1) -- so the exact same
`prepare_input()` used for PolarPicker is reused here unchanged.

Usage:
    /opt/miniconda3/envs/FM_RT/bin/python E_polarcap_compare.py Y1_ge8
    /opt/miniconda3/envs/FM_RT/bin/python E_polarcap_compare.py Y2_ge8
    /opt/miniconda3/envs/FM_RT/bin/python E_polarcap_compare.py Y2_ge18

Output:
    02-data/E_Po/polarcap_<DATASET>.h5        PolarCAP polarity + confidence
    02-data/E_Po_py/PC_<DATASET>_summary.csv  per-station agreement table

Environment: FM_RT (keras 3.12.1, tensorflow)
"""
import sys, time, warnings
from pathlib import Path

import numpy as np
import h5py
import pandas as pd

warnings.filterwarnings('ignore')

if len(sys.argv) < 2 or sys.argv[1] not in ('Y1_ge8', 'Y2_ge8', 'Y2_ge18'):
    print('Usage: python E_polarcap_compare.py Y1_ge8|Y2_ge8|Y2_ge18')
    sys.exit(1)
DATASET = sys.argv[1]

PROJECT   = Path('/Users/mczhang/Documents/GitHub/FM7')
DATA      = PROJECT / '02-data'
SAIPY     = Path('/Users/mczhang/Documents/GitHub/SAIPy')
WAVES_H5  = DATA / 'A_ID' / f'waves_{DATASET}.h5'
DL_H5     = DATA / 'E_Po' / f'dl_polarity_{DATASET}.h5'
PC_H5     = DATA / 'E_Po' / f'polarcap_{DATASET}.h5'
OUT_E     = DATA / 'E_Po_py'; OUT_E.mkdir(exist_ok=True)

DL_N   = 64
SRC_FS = 200
CONF_MIN = 0.8   # same high-confidence threshold used in the CC/DL comparison

print(f'=== PolarCAP vs PolarPicker DL — {DATASET} ===')
print(f'Waves        : {WAVES_H5}')
print(f'PolarPicker  : {DL_H5}')
print(f'PolarCAP wts : {SAIPY / "saipy/saved_models/PolarCAP.h5"}')

# ── Load PolarCAP (compile=False avoids a Keras-3 legacy-optimizer deserialize error) ──
import keras
print(f'\nkeras {keras.__version__}')
t0 = time.time()
model = keras.models.load_model(str(SAIPY / 'saipy/saved_models/PolarCAP.h5'), compile=False)
print(f'Loaded PolarCAP in {time.time()-t0:.1f}s  input={model.input_shape}')

# ── Load waveforms + existing PolarPicker polarity ─────────────────────────
with h5py.File(WAVES_H5, 'r') as hf:
    event_ids = hf['event_ids'][:]
    stations  = sorted(hf['waveforms'].keys())
    waves     = {s: hf[f'waveforms/{s}'][:] for s in stations}
N = len(event_ids)

with h5py.File(DL_H5, 'r') as hf:
    dl_eids = hf['event_ids'][:]
    dl_pol  = {s: hf[f'polarity/{s}'][:]   for s in stations if f'polarity/{s}'   in hf}
    dl_conf = {s: hf[f'confidence/{s}'][:] for s in stations if f'confidence/{s}' in hf}
assert np.array_equal(event_ids, dl_eids), 'event_id ordering mismatch between waves and dl_polarity files'

print(f'Events: {N:,}   Stations: {len(stations)}')

# ── Preprocessing (identical to dl_polarity_ge8.py) ────────────────────────
def prepare_input(wave_264):
    valid = ~np.isnan(wave_264[:, 0])
    w100  = wave_264[:, ::2]
    w64   = w100[:, :DL_N].copy()
    peak  = np.max(np.abs(w64), axis=1, keepdims=True)
    peak[peak == 0] = 1.0
    w64   = w64 / peak
    w64[~valid] = 0.0
    return w64[:, :, np.newaxis].astype(np.float32), valid

# ── Run PolarCAP per station ────────────────────────────────────────────────
pc_pol, pc_conf = {}, {}
t0 = time.time()
for i, sta in enumerate(stations):
    X, valid = prepare_input(waves[sta])
    _, class_probs = model.predict(X, batch_size=512, verbose=0)
    polarity   = np.where(np.argmax(class_probs, axis=1) == 1, 1, -1).astype(np.int8)
    confidence = np.max(class_probs, axis=1).astype(np.float32)
    polarity[~valid]   = 0
    confidence[~valid] = np.nan
    pc_pol[sta]  = polarity
    pc_conf[sta] = confidence
    print(f'  [{i+1:2d}/{len(stations)}] {sta}: {int(valid.sum()):>6,} waveforms  '
          f'Up={int((polarity==1).sum()):>5,}  Dn={int((polarity==-1).sum()):>5,}')
print(f'PolarCAP inference: {time.time()-t0:.1f}s')

with h5py.File(PC_H5, 'w') as hf:
    hf.create_dataset('event_ids', data=event_ids)
    for sta in stations:
        hf.create_dataset(f'polarity/{sta}',   data=pc_pol[sta])
        hf.create_dataset(f'confidence/{sta}', data=pc_conf[sta])
print(f'Saved -> {PC_H5}')

# ── Per-station agreement: PolarCAP vs PolarPicker ─────────────────────────
SUMMARY = []
for sta in stations:
    if sta not in dl_pol:
        continue
    pcp, pcc = pc_pol[sta], pc_conf[sta]
    dlp, dlc = dl_pol[sta], dl_conf[sta]

    both      = (pcp != 0) & (dlp != 0)
    both_hq   = both & (pcc >= CONF_MIN) & (dlc >= CONF_MIN)
    agree     = both    & (pcp == dlp)
    agree_hq  = both_hq & (pcp == dlp)

    n_pc, n_dl = int((pcp != 0).sum()), int((dlp != 0).sum())
    n_both, n_agree       = int(both.sum()),    int(agree.sum())
    n_both_hq, n_agree_hq = int(both_hq.sum()), int(agree_hq.sum())

    SUMMARY.append(dict(
        station=sta, n_pc=n_pc, n_dl=n_dl,
        n_both=n_both, n_agree=n_agree,
        agree_pct=round(n_agree/max(n_both,1)*100, 1),
        n_both_hq=n_both_hq, n_agree_hq=n_agree_hq,
        agree_hq_pct=round(n_agree_hq/max(n_both_hq,1)*100, 1),
    ))

df = pd.DataFrame(SUMMARY)
print()
print(f'=== PolarCAP vs PolarPicker agreement — {DATASET} (conf >= {CONF_MIN} for hq) ===')
print(df.to_string(index=False))
csv_path = OUT_E / f'PC_{DATASET}_summary.csv'
df.to_csv(csv_path, index=False)
print(f'\nSaved -> {csv_path}')

overall_both      = sum(r['n_both']      for r in SUMMARY)
overall_agree     = sum(r['n_agree']     for r in SUMMARY)
overall_both_hq   = sum(r['n_both_hq']   for r in SUMMARY)
overall_agree_hq  = sum(r['n_agree_hq']  for r in SUMMARY)
print(f'\nOverall: {overall_agree}/{overall_both} = {overall_agree/max(overall_both,1)*100:.1f}%  '
      f'(hq: {overall_agree_hq}/{overall_both_hq} = {overall_agree_hq/max(overall_both_hq,1)*100:.1f}%)')
