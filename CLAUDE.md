# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

FM4 is a MATLAB-based focal mechanism (FM) determination pipeline for seismic events at **Axial Seamount** (~45.9–46.0°N, 129.95–130.06°W). It processes waveforms from a 22-station ocean bottom seismometer (OBS) network to compute earthquake focal mechanisms (strike, dip, rake) using SKHASH and HASH inversion algorithms.

## Running the Pipeline

All scripts are run directly in MATLAB. Before running any pipeline script, initialize the MATLAB path:

```matlab
run FM_buildpath4.m
```

This adds `Axial-AutoLocate`, `AutomaticFM`, `FM`, `FM4/01-scripts`, and `FM4/02-data` to the path.

Pipeline scripts are run sequentially by letter stage (A → I). Each script stands alone — run the `.m` file directly, e.g.:

```matlab
run 01-scripts/G_FM_SKHASH_22OBS.m
run 01-scripts/H_read_SKHASH.m
run 01-scripts/I_plotFM_SKHASH.m
```

After running G (which writes SKHASH input `.dat` files), SKHASH must be executed externally before running H:

```bash
cd /Users/mczhang/Documents/GitHub/SKHASH/SKHASH
# run SKHASH on the generated input — output goes to examples/hash3/OUT/out.txt
```

## Architecture

### Data Flow

```
Seismic_Catalog.txt + mseed waveforms
    → A: waveform extraction       → 02-data/A_ID/A_Wave2F_15station.mat
    → B: cross-correlation         → 02-data/B_CC/
    → C: SVD clustering            → 02-data/C_SVD/
    → D: manual polarity picking   → 02-data/D_man/
    → E: polarity combination      → 02-data/E_Po/E_Po_15OBS.mat
    → F: hierarchical clustering   → 02-data/F_Cl/Filter_Felix_combined.mat
    → G: write SKHASH input .dat   → SKHASH/examples/hash3/IN/
    → [run SKHASH externally]      → SKHASH/examples/hash3/OUT/out.txt
    → H: parse SKHASH output       → 02-data/G_FM/G_3D.mat (event1, event1D structs)
    → I: beach ball visualization  → 03-output-graphics/
```

### Key Data Structures

The central data structure throughout the pipeline is an array of structs `Po` (from `Po_Clu`), where each element represents one event-station observation with fields:
- `ID`, `on` (origin time), `lat`, `lon`, `depth`
- `Po_<station>`: polarity value (positive=Up, negative=Down, 0=no pick) for each of 21 stations
- `NSP_<station>`: noise/signal/amplitude vector `[Nos, Nop, Sam, Pam]`
- `cluster`, `ALL` (polarity count)

After H, results are stored in `event1` / `event1D` structs with FM parameters:
- `avmech`: `[strike, dip, rake]`
- `mechqual`: `'A'`/`'B'`/`'C'`/`'D'`
- `faultType`: `'N'` (normal), `'R'` (reverse), `'S'` (strike-slip), `'U'` (undefined)
- `polnum`, `polmisfit`, `stdr`, `namp`, `color2`

### Station Names

21 stations in two groups:
- Legacy (7): `AS1 AS2 CC1 EC1 EC2 EC3 ID1`
- Main array (14): `01A 02A 03A 04A 05A 06A 07A 08A 09A 10A 11A 12A 13A 14A`

Station fields in structs use these names directly (e.g., `Po_01A`, `NSP_14A`).

### External Dependencies (sibling repos)

| Repo | Path | Role |
|------|------|------|
| Axial-AutoLocate | `../Axial-AutoLocate` | Hypocenter location, HASH utilities |
| AutomaticFM | `../AutomaticFM` | Legacy FM code |
| SKHASH | `../SKHASH/SKHASH` | Focal mechanism inversion binary |
| FM | `../FM` | Shared MATLAB utilities |

### SKHASH Input/Output Format

G writes a `.dat` file where each event block starts with `#`, followed by a header line (event index, origin time, lon, lat, depth), then sub-event lines, then per-station lines with: station name, polarity (`U`/`D`), noise levels (Nos, Nop), and amplitudes (Pam, Sam).

H reads `out.txt` via `readtable` → `table2struct`, producing the `events` array, then enriches it with auxiliary plane calculation and fault-type classification.

### Figure Output

`savemyfigureFM4.m` (at repo root) is the standard save utility for figures. Pass a figure handle and output name; it writes to `03-output-graphics/MyPng/` and `03-output-graphics/MyPdf/`.

## Important Hardcoded Paths

Scripts contain absolute paths to the user's local machine. When modifying scripts, update:
- `path_FM` — SKHASH installation directory
- `path` — FM4 data directory (`02-data/`)
- `load(...)` / `filename` — specific `.mat` or `out.txt` paths

---

## File Format Reference

This section documents every non-.mat file format so the pipeline can be reconstructed or debugged on any machine. (`.mat` files are gitignored but their contents are described under each stage below.)

---

### Root-Level Files

#### `Seismic_Catalog.txt`
Fixed-width ASCII, 4 columns:
```
     Date                 Latitude     Longitude    Depth
22-01-2015 02:02:41  45.9413      130.0164     0.77
```
Longitude is stored as positive (130.xx) — scripts negate it internally to -130.xx.

#### `hypo71.dat` (root and `02-data/A_all/`)
Binary Hypo71 format — not human-readable. Contains hypocenter solutions.

#### `Po_Clu_noW.mat` / `Po_Clu_ShallowEast.mat`
MATLAB structs — the `Po_Clu` array (same schema as `Po` described in Key Data Structures above). `ShallowEast` is a spatial/depth subset.

---

### `02-data/A_all/` — Raw Catalog and Phase Data

#### `hypoDD.temp.202209.202308.pha`
HypoDD phase format. Event header line starts with `#`:
```
# 2022  9  1  9 35 49.61  45.9560 -129.9900  2.52 -1.39  0.00  0.00  0.00 224400000
OOAXAS1       0.993   0.5536   P
OOAXAS1       1.768   0.6420   S
```
Header columns: `YY MM DD HH MIN SEC Lat Lon Depth dx dy dz nPh EventID`
Phase columns: `NetworkStation  ArrivalTime(s)  Weight  Phase(P/S)`

#### `nll.temp.202209.202308.txt`
NonLinLoc hypocenter output, tab-delimited with header row:
```
ID  YR  MO  DY  HR  MIN  SEC  latitude  longitude  DEPTH  rms  nPh  Gap
224400000  2022  9  1  9  35  49.916  45.9492  -129.9938  0.478  0.0195  12  150.1
```
Depth in km, Gap in degrees.

#### `.mat` files in `02-data/A_all/` (gitignored, generated by A scripts)
- `Kaiwen_eventData_2F.mat` — event struct array (ID, time, lat, lon, depth)
- `ForMaleenV4.mat` (11 MB) — waveform + event data for external collaborator
- `A_Kaiwen_phase2F.mat` — phase pick data

---

### `02-data/A_ID/` — Waveform Matrices (gitignored)

- `A_Wave2F_15station.mat` (105 MB) — main waveform matrix, struct per station with fields `wave`, `dt`, `tstart`, `ID`

---

### `02-data/G_FM/` — Focal Mechanism Results

#### `Matched_event_pairs.csv`
CSV, 13 columns:
```
Po_Clu_ID, Po_Clu_Time, Po_Clu_Lat, Po_Clu_Lon, Po_Clu_Depth, Po_Clu_faultType,
Event1_ID, Event1_Time, Event1_Lat, Event1_Lon, Event1_Depth, Event1_faultType,
TimeDiff_sec
```
Example row:
```
228000017,07-Oct-2022 14:13:32,45.94777,-129.99291,0.39,N,130,07-Oct-2022 14:13:31,45.955,-129.98801,1.66,N,0.430
```

#### `.mat` files in `02-data/G_FM/` (gitignored)
| File | Size | Contents |
|------|------|----------|
| `G_3D.mat` | 776 KB | `event1` struct array with full FM params (strike/dip/rake + quality) — primary output used by H and I |
| `G_FM_SKHASH_realEle.mat` | 1.8 MB | SKHASH solutions, real elevation |
| `G_FM_SKHASHV2.mat` | 1.8 MB | SKHASH V2 solutions |
| `G_2F_HASH.mat` | 101 MB | HASH (1D) solutions |
| `G_2F_HASH_V2.mat` | 103 MB | HASH V2 solutions |
| `G_2F_HASH_V3_singlevel2.mat` | 102 MB | HASH V3 (excluded from git, >100 MB) |
| `G_FM_15station.mat` | 135 KB | 15-station subset solutions |

---

### `01-scripts/HASH/` — HASH Algorithm I/O Files

#### `phase.dat`
HASH polarity input. Per-event block:
```
2022 9 6 12 29 29.50  45.56 -129.94  1.08  0.30  0.20  1.00   1
AAS1 OO HHZ I U
AAS2 OO HHZ I D
```
Event line: `YY MM DD HH MIN SS.ss Lat Lon Depth HorzErr VertErr RMS EventNum`
Station line: `Station Network Channel Interior/Exterior(I/E) Polarity(U/D)`

#### `station.dat`
Station coordinates:
```
A01A HHZ  45.93360 -129.99920  0  OO
A02A HHZ  45.93380 -130.01410  0  OO
```
Columns: `StationCode  Channel  Lat  Lon  Elev(m)  Network`

#### `amp.dat`
S/P amplitude ratios, event blocks:
```
1     12
AAS1 HHZ OO   40.254   40.254   587.615   750.566
AAS2 HHZ OO   47.351   47.351   349.971   506.049
```
Station columns: `Station Channel Network  Samp1  Samp2  Pamp1  Pamp2`

#### `reverse.dat`
Empty for Axial datasets (no known polarity reversals in this network).

#### `velmod*.dat` (velmod1, velmod2, velmod_E1–E4, velmod_W1–W2, velmod_S1, velmod_AXAS1–2, velmod_AXCC1, velmod_AXEC1–3, velmod_AXID1)
1D velocity model, 2 columns, no header:
```
0      2.11493
0.05   2.21779
0.10   2.33263
...
40     8.2
```
Columns: `Depth(km)  Vp(km/s)`

---

### `01-scripts/SKHASH/examples/hash3/` — Active SKHASH Configuration

#### `control_file.txt`
Key-value driver, `$keyword value` format. Key parameters:
```
$input_format   hash3
$stfile         examples/hash3/IN/scsn.stations.txt
$plfile         examples/hash3/IN/scsn.reverse.txt
$corfile        examples/hash3/IN/north3.statcor.txt
$ampfile        examples/hash3/IN/north3.amp.txt
$fpfile         examples/hash3/IN/north2.phase.txt
$outfile1       examples/hash3/OUT/out.txt
$outfile2       examples/hash3/OUT/out2.txt
$npolmin        15        # min polarity picks required
$dang           5         # grid search spacing (degrees)
$nmc            30        # Monte Carlo trials
$maxout         300       # max output mechanisms
$ratmin         2         # min SNR
$badfrac        0.1       # allowed bad polarity fraction
$qbadfrac       0.3       # amplitude noise threshold (log10)
$delmax         25        # max source-receiver distance (km)
$cangle         45        # probability cone angle
$prob_max       0.75      # multiple-solution probability threshold
$vmodel_paths   examples/velocity_models/velmod_W1.txt  ...
```

#### `IN/scsn.stations.txt`
SCSN station format, fixed-width:
```
AAS1 HHZ BIG CHUCKAWALLA MTNS   45.93360 -129.99920   0  1997/09/19  3000/01/01  OO
```
Columns: `Code  Chan  Description  Lat  Lon  Elev  StartDate  EndDate  Network`

#### `IN/north2.phase.txt`
SKHASH phase input. Same block format as `phase.dat` above but with SKHASH-specific padding in the event header line.

#### `IN/north3.amp.txt`
SKHASH amplitude input. 6 amplitude values per station (vs. 4 in HASH):
```
1         12
AAS1 HHZ OO  59.08  59.08  59.079  59.079  53.675  56.286
```
Columns: `Station Chan Net  Samp Samp  Pamp Pamp  Pamp_corr Samp_corr`

#### `IN/north3.statcor.txt`
Station static timing corrections:
```
AAS1 HHZ OO  0.0000
AAS2 HHZ OO  0.0000
```
Columns: `Station  Channel  Network  Correction(s)` — all zeros for Axial network.

#### `IN/scsn.reverse.txt`
Historical polarity reversal database, 146 entries:
```
AQUA 19920101 19921231
BAC  19950624 19960124
```
Columns: `StationCode  StartDate(YYYYMMDD)  EndDate(0=ongoing)`

#### `IN/Axial_22OBSs.dat` / `Axial_22OBSs_Aquality.dat`
SKHASH block format written by script G. One block per event cluster:
```
#
  1   738770.5204803   -129.9960   45.9240   1.08
 224900009   738770.5204803  -129.9960   45.9240   1.081
AAS1  U   40.2539   40.2539   587.6152   750.5663
AAS2  D   47.3510   47.3510   349.9713   506.0489
```
Line 1: `#` separator
Line 2: `ClusterIndex  JulianDay  Lon  Lat  Depth(km)`
Line 3: `EventID  JulianDay  Lon  Lat  Depth(km)`
Station lines: `Station  Polarity(U/D)  Samp  Samp  Pamp  Pamp`

#### `OUT/out.txt`
SKHASH CSV output, read by script H via `readtable`:
```
event_id, strike, dip, rake, quality, fault_plane_uncertainty, aux_plane_uncertainty,
num_p_pol, num_sp_ratios, polarity_misfit, prob_mech, sta_distribution_ratio,
sp_misfit, mult_solution_flag, time, magnitude, origin_lat, origin_lon,
origin_depth_km, horz_uncert_km, vert_uncert_km
```
Example row:
```
100,333.3,88.1,-167.3,B,12.9,17.4,17,12,19.1,100.0,62.3,-51.1,False,2023-04-02 00:03:26,1.0,45.955,-129.994,1.79,0.3,0.2
```
Quality: `A` (best) → `D` (worst). `mult_solution_flag` = True if multiple solutions exist.

#### `examples/velocity_models/*.txt`
Same 2-column depth/velocity format as `velmod*.dat` above but with `.txt` extension. One file per station group (E1–E4, W1–W2, S1, AXAS1–2, AXCC1, AXEC1–3, AXID1).

---

### What the Other Computer Needs to Reconstruct

After `git clone`, the `.mat` files are absent. To re-run the full pipeline from scratch:
1. Place raw `.mseed` waveform files in `01-scripts/datamseed/` (not in git — download separately)
2. Run stages A → F to regenerate all intermediate `.mat` files in `02-data/`
3. Run G to write SKHASH input into `01-scripts/SKHASH/examples/hash3/IN/`
4. Run SKHASH: `python 01-scripts/SKHASH/SKHASH.py 01-scripts/SKHASH/examples/hash3/control_file.txt`
5. Run H to parse `OUT/out.txt` → `G_3D.mat`
6. Run I to produce figures

The pre-existing `IN/*.dat` files and `OUT/out.txt` committed in git allow steps 4–6 to be run immediately without re-processing from raw data.
