"""
load_axial_geology.py
=====================
Helpers to load Axial Seamount geology overlays from existing MATLAB/text files:
  - Fissures 1998, 2011, 2015
  - Lava flows 1998, 2011, 2015
Returns lists of (N,2) arrays [lon, lat] suitable for matplotlib fill/plot.
"""
import re
import gzip
import numpy as np
import pandas as pd
from pathlib import Path

FM_DATA   = Path('/Users/mczhang/Documents/GitHub/FM/02-data/Alldata')
LAVA_DIR  = Path('/Users/mczhang/Documents/GitHub/Axial-AutoLocate/axial')
LAVA_1998 = Path('/Users/mczhang/Documents/GitHub/FM/04-final-paper/axial_lava1998.m')

# ── Fissures ───────────────────────────────────────────────────────────────────
def _load_fissures_tab(path):
    """Load tab-separated fissure file with ORIG_FID, LONGITUDE, LATITUDE."""
    df = pd.read_csv(path, sep='\t')
    segs = []
    for fid, grp in df.groupby('ORIG_FID'):
        segs.append(grp[['LONGITUDE','LATITUDE']].values)
    return segs

def _load_fissures_1998(path):
    """Load 1998 fissure file: col0=segment, col1='lon,lat'."""
    df = pd.read_csv(path, sep=r'\s+', header=None)
    segs = []
    for fid, grp in df.groupby(0):
        pts = []
        for val in grp[1]:
            lon, lat = val.split(',')
            pts.append([float(lon), float(lat)])
        segs.append(np.array(pts))
    return segs

def load_fissures():
    f2015 = _load_fissures_tab(
        FM_DATA/'Fissures2015/JdF:Axial_Clague/Axial-2015-fissures-points-geo-v2.txt')
    f2011 = _load_fissures_tab(
        FM_DATA/'Fissures2011/JdF:Axial_Clague/Axial-2011-fissures-points-geo-v2.txt')
    f1998 = _load_fissures_1998(FM_DATA/'Axial-1998-Fissures.txt')
    return {'1998': f1998, '2011': f2011, '2015': f2015}

# ── Lava flows — parse MATLAB .m scripts ───────────────────────────────────────
def _parse_matlab_xy_blocks(text, var_name):
    """
    Parse blocks like:
        var_name(i).xy = [ ...
        lon,lat
        ...
        ];
    Returns list of (N,2) arrays.
    """
    polys = []
    pattern = re.compile(
        rf'{re.escape(var_name)}\(\d+\)\.xy\s*=\s*\[\s*\.\.\.(.*?)\];',
        re.DOTALL)
    for m in pattern.finditer(text):
        block = m.group(1)
        pts = []
        for line in block.splitlines():
            line = line.strip().rstrip('...')
            if not line: continue
            parts = line.replace(',', ' ').split()
            if len(parts) >= 2:
                try:
                    pts.append([float(parts[0]), float(parts[1])])
                except ValueError:
                    pass
        if pts:
            polys.append(np.array(pts))
    return polys

def _parse_matlab_lonlat_blocks(text, var_name):
    """
    Parse blocks like:
        a=[lon,lat\n...];  flow2015(i).lon=a(:,1); flow2015(i).lat=a(:,2);
    Returns list of (N,2) arrays.
    """
    polys = []
    # Each flow starts with a=[  ...  ];
    blocks = re.split(r'(?=^a\s*=\s*\[)', text, flags=re.MULTILINE)
    for block in blocks:
        if not block.strip().startswith('a'): continue
        # extract the matrix content
        m = re.search(r'a\s*=\s*\[(.*?)\];', block, re.DOTALL)
        if not m: continue
        pts = []
        for line in m.group(1).splitlines():
            line = line.strip()
            if not line: continue
            parts = line.replace(',', ' ').split()
            if len(parts) >= 2:
                try:
                    pts.append([float(parts[0]), float(parts[1])])
                except ValueError:
                    pass
        if pts:
            polys.append(np.array(pts))
    return polys

def load_lavas():
    lava1998_txt = LAVA_1998.read_text()
    lava2011_txt = (LAVA_DIR / 'axial_lava2011.m').read_text()
    lava2015_txt = (LAVA_DIR / 'axial_lava2015.m').read_text()

    l1998 = _parse_matlab_xy_blocks(lava1998_txt, 'lava')
    l2011 = _parse_matlab_xy_blocks(lava2011_txt, 'lava')
    l2015 = _parse_matlab_lonlat_blocks(lava2015_txt, 'flow2015')

    print(f"Lava polys — 1998: {len(l1998)}  2011: {len(l2011)}  2015: {len(l2015)}")
    return {'1998': l1998, '2011': l2011, '2015': l2015}


if __name__ == '__main__':
    fiss  = load_fissures()
    lavas = load_lavas()
    for yr, segs in fiss.items():
        print(f"Fissures {yr}: {len(segs)} segments")
    for yr, polys in lavas.items():
        print(f"Lava {yr}: {len(polys)} polygons")
