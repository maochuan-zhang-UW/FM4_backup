#!/usr/bin/env bash
# run_ge8_pipeline.sh
# ====================
# Full ≥8-station CC-augmented FM pipeline.
# Run from: /Users/mczhang/Documents/GitHub/FM7/01-scripts/python/
#
# Prerequisites (already done):
#   waves_Y1_ge8.h5 — cut by cut_waves_Y1_ge8.py  ✓
#   waves_Y2_ge8.h5 — cut by cut_waves_Y2_ge8.py  ✓
#
# Steps:
#   1. DL polarity Y1 (FM_RT)
#   2. DL polarity Y2 (FM_RT)
#   3. Hierarchical clustering + polarity assignment Y1 (FM_ML)
#   4. Hierarchical clustering + polarity assignment Y2 (FM_ML)
#   5. Write SKHASH inputs (FM_ML)
#   6. Run SKHASH Y1
#   7. Run SKHASH Y2

set -e
cd /Users/mczhang/Documents/GitHub/FM7/01-scripts/python

FM_RT=/opt/miniconda3/envs/FM_RT/bin/python
FM_ML=/opt/miniconda3/envs/FM_ML/bin/python
SKHASH_DIR=/Users/mczhang/Documents/GitHub/SKHASH/SKHASH7

echo "=== Step 1: DL polarity Y1 ==="
$FM_RT dl_polarity_ge8.py Y1 2>&1 | tee /tmp/dl_polarity_Y1_ge8.log
echo ""

echo "=== Step 2: DL polarity Y2 ==="
$FM_RT dl_polarity_ge8.py Y2 2>&1 | tee /tmp/dl_polarity_Y2_ge8.log
echo ""

echo "=== Step 3: Hierarchical clustering + polarity assignment Y1 ==="
$FM_ML cluster_polarity.py Y1 2>&1 | tee /tmp/cluster_polarity_Y1.log
echo ""

echo "=== Step 4: Hierarchical clustering + polarity assignment Y2 ==="
$FM_ML cluster_polarity.py Y2 2>&1 | tee /tmp/cluster_polarity_Y2.log
echo ""

echo "=== Step 5: Write SKHASH inputs (Y1 + Y2) ==="
$FM_ML write_skhash_input_ge8.py 2>&1 | tee /tmp/write_skhash_ge8.log
echo ""

echo "=== Step 6: Run SKHASH Y1 ==="
cd $SKHASH_DIR
$FM_ML SKHASH.py examples/hash3/control_file_Y1_ge8.txt 2>&1 | tee /tmp/skhash_Y1_ge8.log
echo ""

echo "=== Step 7: Run SKHASH Y2 ==="
$FM_ML SKHASH.py examples/hash3/control_file_Y2_ge8.txt 2>&1 | tee /tmp/skhash_Y2_ge8.log
echo ""

echo "=== Pipeline complete ==="
echo "Y1 results: $SKHASH_DIR/examples/hash3/OUT_Y1_ge8/out.txt"
echo "Y2 results: $SKHASH_DIR/examples/hash3/OUT_Y2_ge8/out.txt"
