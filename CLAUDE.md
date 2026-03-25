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
