# FM7 Y1+Y2 Focal Mechanism Pipeline — Master Plan

## Overview
Full Python pipeline to compute focal mechanisms for 54,971 earthquakes at Axial Seamount.
Y1: Sep 1 2022 – Aug 31 2023 (6,727 events) | Y2: Sep 1 2023 – Sep 22 2024 (48,244 events)
Both years processed together. Results from SKHASH and HASH compared at the end.

## Access Requirements
Before running, confirm these are all present:

| Item | Path | Status |
|------|------|--------|
| Y1 location file | 02-data/A_all/axial.Y1.mldd.loc.260317 | ✓ exists |
| Y2 location file | 02-data/A_all/axial.Y2.mldd.loc.260317 | ✓ exists |
| Y1 phase file | 02-data/A_all/axial.Y1.mldd.pha.260317 | ✓ exists |
| Y2 phase file | 02-data/A_all/axial.Y2.mldd.pha.260317 | ✓ exists |
| mseed waveforms | FM4/01-scripts/datamseed/{year}/{month}/ | ✓ exists |
| DL Keras model | FM6_RealTime/02-data/PolarPicker_unified_TMSF_001.keras | ✓ exists |
| SKHASH7 binary | SKHASH/SKHASH7/ | ✓ exists |
| HASH binary | FM7/01-scripts/HASH/hash_driver3 | ✓ exists |
| FM_RT conda env | /opt/miniconda3/envs/FM_RT | needs: pandas, obspy, matplotlib, seaborn, jupyter, papermill |

## Environment Setup
Run once before starting:
```bash
conda activate FM_RT
pip install pandas obspy matplotlib seaborn jupyter papermill h5py
```

## How to Run
```bash
cd /Users/mczhang/Documents/GitHub/FM7/01-scripts/python
conda activate FM_RT
python run_pipeline.py
```

The runner executes all 12 notebooks in sequence using papermill.
- Sends email to mczhang26@gmail.com after every step + every hour
- Emails include plots as attachments
- Saves progress to 02-data/pipeline_progress.json — restart is safe, done steps skipped
- Pauses at Step 06 (DL review) and Step 09a (before running SKHASH externally)

## Two External Commands (pipeline pauses and emails you)
After notebook 09a completes:
```bash
cd /Users/mczhang/Documents/GitHub/SKHASH/SKHASH7
python SKHASH.py examples/hash3/control_file.txt
```
After notebook 09b completes:
```bash
cd /Users/mczhang/Documents/GitHub/FM7/01-scripts/HASH
./hash_driver3 < hash_Y1Y2.input
```
Then set "awaiting_skhash": false and "awaiting_hash": false in pipeline_progress.json

---

## Step-by-Step Notebook Guide

### Notebook 01 — Parse mldd Catalog (~5 min)
**Input:**
- 02-data/A_all/axial.Y1.mldd.loc.260317
- 02-data/A_all/axial.Y2.mldd.pha.260317
(and Y2 equivalents)

**What it does:**
- Reads the space-delimited .loc files: YY MM DD HH MIN SS.sss Lat Lon Depth Mag EventID
- Reads the HypoDD-format .pha files: # header line + station arrival lines
- Combines Y1 + Y2 into unified DataFrames
- Validates: no duplicate IDs, reasonable lat/lon ranges, time order

**Output:**
- 02-data/A_all/catalog_Y1Y2.parquet  (columns: event_id, year, time, lat, lon, depth, mag)
- 02-data/A_all/phases_Y1Y2.parquet  (columns: event_id, DDt_AS1...DDt_14A, DDSt_AS1...DDSt_14A)

**Email after step includes:**
- Plot 1: Map of all events (scatter plot colored by year Y1/Y2)
- Plot 2: Depth histogram
- Plot 3: Monthly event count time series
- Plot 4: Station coverage heatmap (events per station)
- Stats: total events, Y1/Y2 counts, events with ≥7 stations with picks

**What to check:**
- Are event counts correct? (6727 Y1, 48244 Y2)
- Are there any events with unreasonable lat/lon/depth?
- Is station coverage reasonable (all 21 stations have picks)?

---

### Notebook 02 — Cut Waveforms (~3–5 hours)
**Input:**
- 02-data/A_all/catalog_Y1Y2.parquet
- 02-data/A_all/phases_Y1Y2.parquet
- FM4/01-scripts/datamseed/{year}/{month}/YYYY-MM-DD-HH-00-00.mseed

**What it does:**
- For each event × station:
  - Finds the hourly mseed file containing the P pick time
  - Reads it with ObsPy, selects the HHZ channel for the station
  - Cuts window [-0.32, +1.0]s around the P pick at 200 Hz → 264 samples
  - Bandpass filters 4–50 Hz
  - Calculates noise/S/P amplitudes (NSP)
- Saves every 1000 events as a checkpoint to allow restart
- Shows tqdm progress bar per station

**Output:**
- 02-data/A_ID/waves_Y1Y2.h5
  - /waveforms/AS1    shape (N_events, 264) float32
  - /amplitudes/AS1   shape (N_events, 3)   float32 [noise_amp, S_amp, P_amp]
  - (repeated for all 21 stations)
  - /event_ids        shape (N_events,)      int64

**Email after step includes:**
- Plot 1: % waveforms successfully cut per station (bar chart)
- Plot 2: SNR distribution per station (violin plot)
- Plot 3: 5×4 grid of sample waveforms (one per station, random event)
- Stats: total event-station pairs cut, % success, mean SNR per station

**What to check:**
- Any stations with very low success rate (<50%)? Check mseed naming convention.
- Any stations with very low SNR? Might indicate instrument issue.
- Do waveforms look reasonable (P arrival at t=0)?

---

### Notebook 03 — Cross-Correlation (~1–2 hours)
**Input:** 02-data/A_ID/waves_Y1Y2.h5

**What it does:**
- For each station:
  - Loads waveform matrix (N × 264)
  - Extracts CC window: starts at sample int((0.32-0.25)*200)=14, length int(0.2/dt) = 40 samples
  - Computes pairwise CC using batched matrix multiplication (efficient numpy)
  - Keeps only pairs with CC > 0.70 (threshold from B_CC_all.m)
- Handles both short (0.2s) and long (0.4s) windows

**Output:**
- 02-data/B_CC/cc_Y1Y2.h5
  - /AS1   shape (M_AS1, 3) columns [ID1, ID2, cc_value]
  - (repeated for all 21 stations, M varies per station)

**Email after step includes:**
- Plot 1: CC value distribution histogram per station
- Plot 2: Number of high-CC pairs per station (bar chart)
- Plot 3: CC network graph (nodes=events, edges=CC>0.85, top 500 pairs, colored by year)
- Stats: pairs per station, mean CC value

**What to check:**
- Stations with very few pairs may have noisy data
- Stations with too many pairs (>100k) may have too-low threshold

---

### Notebook 04 — SVD Clustering (~30 min)
**Input:**
- 02-data/B_CC/cc_Y1Y2.h5
- 02-data/A_ID/waves_Y1Y2.h5

**What it does:**
- For each station:
  - Builds sparse CC matrix from pairs
  - Computes leading eigenvector via scipy.sparse.linalg.eigsh
  - SVD score = eigenvector component for each event (positive/negative = polarity proxy)
  - Higher |score| = more reliable polarity indicator

**Output:**
- 02-data/C_SVD/svd_Y1Y2.parquet (columns: event_id, AS1_score...14A_score)

**Email after step includes:**
- Plot 1: SVD score distribution per station (violin plot)
- Plot 2: Waveform comparison — top 5 vs bottom 5 SVD score events per station
- Stats: events with |SVD_score| > 0.1 per station

**What to check:**
- Bimodal SVD score distribution is expected (positive/negative polarity clusters)
- If flat/unimodal, the CC data may be insufficient

---

### Notebook 05 — Select 20 High-SNR Events per Station (~5 min)
**Input:**
- 02-data/C_SVD/svd_Y1Y2.parquet
- 02-data/A_ID/waves_Y1Y2.h5 (amplitudes group)

**What it does:**
- For each station: ranks events by combined score = SNR × |SVD_score|
- Selects top 20 events (for DL verification in next step)
- These are the events you will visually inspect

**Output:**
- 02-data/D_man/selected_20perSta_Y1Y2.parquet

**Email after step includes:**
- Plot 1: Waveform grid — 20 rows × 21 columns = 420 waveforms
  (one panel per selected event per station, P pick at center)
- Plot 2: SNR vs SVD score scatter for all events, selected highlighted
- Stats: table of selected event IDs per station

**⚠️ WHAT TO CHECK (important!):**
- Look at the waveform grid — do the P arrivals look clean?
- Are there any clearly wrong picks (artifact waveforms)?
- If yes: note the event IDs, they can be excluded before Step 06

---

### Notebook 06 — DL Polarity (~1–2 hours)
**Input:**
- 02-data/A_ID/waves_Y1Y2.h5
- 02-data/D_man/selected_20perSta_Y1Y2.parquet
- FM6_RealTime/02-data/PolarPicker_unified_TMSF_001.keras

**What it does:**
STEP A (runs first, then PAUSES):
- Resamples 200 Hz → 100 Hz waveforms, takes first 64 samples [-0.32, +0.31]s
- Runs predict_polarity.py on the 20 selected events per station
- Saves DL polarities + confidence scores
- PAUSES pipeline — sends you review email with waveform+prediction plots

STEP B (runs after you clear the review flag):
- Runs DL on ALL 54,971 events
- Uses --cache flag to skip events already predicted

**Output:**
- 02-data/D_man/dl_selected_Y1Y2.mat   (Step A — 20 per station, for review)
- 02-data/E_Po/dl_polarity_Y1Y2.mat    (Step B — all events)

**Email Step A (REVIEW REQUIRED):**
- Plot 1: 420-panel waveform grid with DL label (U↑/D↓) and confidence bar
- Plot 2: Confidence histogram (are most picks >0.6?)
- Stats: pick counts per station, % high confidence
- ⚠️ Subject line: [FM7 ⚠️ REVIEW NEEDED]
- To continue: set "awaiting_review": false in pipeline_progress.json

**Email Step B:**
- Plot 1: Confidence histogram for all 54,971 events
- Plot 2: Pick rate per station (% events with DL pick)
- Stats: total picks, % confident (>0.6), picks per station

**What to check:**
- Step A: Does DL agree with what you would pick manually?
  Any systematic errors for specific stations?
- Step B: Are low-confidence picks (<30%) at specific stations?
  This might mean DL doesn't generalize well to those stations.

---

### Notebook 07 — Polarity Combine (~10 min)
**Input:**
- 02-data/E_Po/dl_polarity_Y1Y2.mat
- 02-data/C_SVD/svd_Y1Y2.parquet

**What it does:**
- Primary: use DL polarity (±1) where available
- Secondary: for events with no DL pick, propagate polarity using SVD sign
- Encode polarity strength: ±1 = DL confident, ±2 = DL uncertain, ±3 = SVD only, 0 = no pick
- Calculate PoALL = total number of stations with any pick per event

**Output:**
- 02-data/E_Po/polarity_combined_Y1Y2.parquet
  columns: event_id, Po_AS1...Po_14A, NSP_AS1...NSP_14A, PoALL

**Email after step includes:**
- Plot 1: PoALL distribution histogram (how many stations per event)
- Plot 2: Pick source breakdown per station (DL confident / DL uncertain / SVD only)
- Stats: events with PoALL ≥ 5 (will enter clustering), total picks by source

**What to check:**
- How many events have PoALL ≥ 5? This is the pool for focal mechanisms.
- Is the DL source dominant? Good. SVD-only should be minority.

---

### Notebook 08 — Hierarchical Clustering (~20 min)
**Input:** 02-data/E_Po/polarity_combined_Y1Y2.parquet

**What it does:**
- Filters: PoALL ≥ 5 (same as F_Cl_Real.m)
- Builds polarity matrix (N × 21) and log10(S/P ratio) matrix
- Converts lat/lon to local km coordinates
- Runs hierarchical clustering with a=0.2
- Removes events with waveform quality issues (|waveform amplitude| = 0)

**Output:**
- 02-data/F_Cl/clustering_Y1Y2.parquet
- 02-data/F_Cl/removed_wave_IDs_Y1Y2.parquet

**Email after step includes:**
- Plot 1: Cluster size distribution (histogram)
- Plot 2: Map of clusters (each cluster a different color)
- Plot 3: Depth profile of clusters
- Stats: total clusters, median cluster size, events removed

**What to check:**
- Are clusters geographically coherent?
- Are there very large clusters (>100 events)? These are excluded from SKHASH/HASH.

---

### Notebook 09a — Write SKHASH Input (~10 min)
**Input:** 02-data/F_Cl/clustering_Y1Y2.parquet

**What it does:**
- Writes SKHASH block format .dat file
- One block per cluster: # separator, cluster header, event sub-lines, station polarity lines
- Format exactly matches Axial_22OBSs.dat (see CLAUDE.md for spec)

**Output:**
- SKHASH/SKHASH7/examples/hash3/IN/Axial_22OBSs_Y1Y2.dat
- 02-data/G_FM/skhash_input_summary_Y1Y2.parquet

**Email — PAUSE (external command needed):**
- Stats: number of event blocks written, stations per block
- ⚠️ Instructions to run SKHASH:
  cd /Users/mczhang/Documents/GitHub/SKHASH/SKHASH7
  python SKHASH.py examples/hash3/control_file.txt
- Set "awaiting_skhash": false when done

---

### Notebook 09b — Write HASH Input (~10 min)
**Input:** 02-data/F_Cl/clustering_Y1Y2.parquet

**What it does:**
- Writes HASH phase.dat format (different from SKHASH format)
- Writes amp.dat with S/P amplitude ratios
- Reuses existing station.dat

**Output:**
- 01-scripts/HASH/phase_Y1Y2.dat
- 01-scripts/HASH/amp_Y1Y2.dat

**Email — PAUSE (external command needed):**
- ⚠️ Instructions to run HASH:
  cd /Users/mczhang/Documents/GitHub/FM7/01-scripts/HASH
  ./hash_driver3 < hash_Y1Y2.input
- Set "awaiting_hash": false when done

---

### Notebook 10a — Parse SKHASH Output (~5 min)
**Input:** SKHASH/SKHASH7/examples/hash3/OUT/out.txt

**What it does:**
- Reads CSV output with readtable-equivalent (pandas)
- Enriches with auxiliary plane calculation (using Aki-Richards formulas)
- Classifies fault type: N (normal), R (reverse), S (strike-slip), U (undefined)
- Matches to catalog by event_id

**Output:**
- 02-data/G_FM/focal_mechanisms_SKHASH_Y1Y2.parquet
- 02-data/G_FM/G_3D_SKHASH_Y1Y2.mat (MATLAB-compatible)

**Email after step includes:**
- Plot 1: Quality distribution pie chart (A/B/C/D)
- Plot 2: Fault type distribution (N/R/S/U)
- Plot 3: Map with colored dots by fault type
- Plot 4: Sample beach balls (10 A-quality events)
- Stats: total FMs, quality breakdown, fault type breakdown

---

### Notebook 10b — Parse HASH Output (~5 min)
**Input:** 01-scripts/HASH/hashout_Y1Y2.dat

Same structure as 10a but for HASH results.

**Output:**
- 02-data/G_FM/focal_mechanisms_HASH_Y1Y2.parquet
- 02-data/G_FM/G_3D_HASH_Y1Y2.mat

---

### Notebook 11 — DL vs Pipeline Polarity Comparison (~10 min)
**Input:**
- 02-data/E_Po/dl_polarity_Y1Y2.mat (DL only)
- 02-data/E_Po/polarity_combined_Y1Y2.parquet (pipeline combined)

**What it does:**
- For each event-station pair where both have picks:
  computes mismatch = sign(DL) ≠ sign(pipeline)
- Breaks down by: station, year (Y1/Y2), legacy vs main array, DL confidence

**Output:**
- 02-data/G_FM/polarity_comparison_Y1Y2.parquet

**Email after step includes:**
- Plot 1: Mismatch rate per station (bar chart)
- Plot 2: Confusion matrix (DL Up/Down vs pipeline Up/Down) — all stations combined
- Plot 3: Mismatch rate vs DL confidence (scatter)
- Plot 4: Y1 vs Y2 mismatch comparison
- Stats: overall mismatch rate, best/worst stations, % confident DL picks that agree

---

### Notebook 12 — SKHASH vs HASH FM Comparison (~10 min)
**Input:**
- 02-data/G_FM/focal_mechanisms_SKHASH_Y1Y2.parquet
- 02-data/G_FM/focal_mechanisms_HASH_Y1Y2.parquet

**What it does:**
- Matches events present in both solutions
- Computes Kagan angle between each pair of focal mechanisms
- Kagan angle < 30° = good agreement
- Fault type agreement analysis

**Output:**
- 02-data/G_FM/fm_comparison_SKHASH_HASH_Y1Y2.parquet

**Email after step includes:**
- Plot 1: Kagan angle distribution histogram
- Plot 2: Strike/dip/rake scatter SKHASH vs HASH
- Plot 3: Beach ball grid — 12 events with largest Kagan angle disagreement
- Plot 4: Map colored by Kagan angle (where do they disagree most?)
- Stats: n events in both, % with Kagan < 30°, median Kagan angle

---

## Intermediate Data Map

```
02-data/
├── A_all/
│   ├── axial.Y1.mldd.loc.260317     ← INPUT
│   ├── axial.Y2.mldd.loc.260317     ← INPUT
│   ├── axial.Y1.mldd.pha.260317     ← INPUT
│   ├── axial.Y2.mldd.pha.260317     ← INPUT
│   ├── catalog_Y1Y2.parquet         ← NB01 output
│   └── phases_Y1Y2.parquet          ← NB01 output
├── A_ID/
│   └── waves_Y1Y2.h5                ← NB02 output
├── B_CC/
│   └── cc_Y1Y2.h5                   ← NB03 output
├── C_SVD/
│   └── svd_Y1Y2.parquet             ← NB04 output
├── D_man/
│   ├── selected_20perSta_Y1Y2.parquet  ← NB05 output
│   └── dl_selected_Y1Y2.mat            ← NB06A output
├── E_Po/
│   ├── dl_polarity_Y1Y2.mat         ← NB06B output
│   └── polarity_combined_Y1Y2.parquet  ← NB07 output
├── F_Cl/
│   ├── clustering_Y1Y2.parquet      ← NB08 output
│   └── removed_wave_IDs_Y1Y2.parquet  ← NB08 output
├── G_FM/
│   ├── skhash_input_summary_Y1Y2.parquet  ← NB09A output
│   ├── focal_mechanisms_SKHASH_Y1Y2.parquet  ← NB10A output
│   ├── G_3D_SKHASH_Y1Y2.mat            ← NB10A output
│   ├── focal_mechanisms_HASH_Y1Y2.parquet   ← NB10B output
│   ├── G_3D_HASH_Y1Y2.mat              ← NB10B output
│   ├── polarity_comparison_Y1Y2.parquet   ← NB11 output
│   └── fm_comparison_SKHASH_HASH_Y1Y2.parquet  ← NB12 output
└── pipeline_progress.json           ← runner state (auto-updated)
```

## Timeline Estimate

| Step | Estimated time |
|------|---------------|
| 01 Parse catalog | 5 min |
| 02 Cut waveforms | 3–5 hours |
| 03 Cross-correlation | 1–2 hours |
| 04 SVD clustering | 30 min |
| 05 Select 20 events | 5 min |
| 06A DL on 20 selected | 10 min + YOUR REVIEW |
| 06B DL on all events | 1–2 hours |
| 07 Polarity combine | 10 min |
| 08 Hierarchical clustering | 20 min |
| 09A Write SKHASH input | 10 min + SKHASH run (external) |
| 09B Write HASH input | 10 min + HASH run (external) |
| 10A Parse SKHASH | 5 min |
| 10B Parse HASH | 5 min |
| 11 Polarity comparison | 10 min |
| 12 FM comparison | 10 min |
| **Total automated** | **~8–11 hours** |

## Things to Modify Later

This plan is intentionally in a markdown file so you can edit it. Common changes:

- **CC threshold**: change `CC_THRESHOLD` in pipeline_config.py (default 0.70)
- **Min polarity picks**: change `MIN_PO_NUM` (default 5)
- **DL confidence cutoff**: change `DL_CONF_MIN` (default 0.60)
- **Clustering parameter**: change `A_VALUE` (default 0.2)
- **Add a station**: add to `STATIONS_ALL` and `STATION_COORDS` and `MSEED_STA_MAP`
- **Run only specific notebooks**: edit the `STEPS` list in run_pipeline.py
