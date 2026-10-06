"""
X_compare_DL_polarity.py
Compare DL polarity predictions from two models against manual ground truth
in Filter_Felix_combined.mat (FM4 Stage F output), and save TMSF_001
predictions to DL_polarity_TMSF001.mat for use by G_FM_HASH_22OBS_DL.m.

Ground truth : Po_Clu[i].Po_<sta>  (scalar: +1=Up, -1=Down, 0=no pick)
DL input     : Po_Clu[i].W_<sta>   (250-sample 200 Hz waveform)
Stations     : AS1 AS2 CC1 EC1 EC2 EC3 ID1  (7 legacy stations the model supports)

Usage:
    /opt/miniconda3/envs/FM_RT/bin/python X_compare_DL_polarity.py
"""

import numpy as np
import scipy.io as sio
import scipy.signal as ss
import keras
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os

# ── Paths ─────────────────────────────────────────────────────────────────────
MODELS = {
    "TMSF_001": "/Users/mczhang/Documents/GitHub/FM6_RealTime/02-data/PolarPicker_unified_TMSF_001.keras",
    "20260308":  "/Users/mczhang/Documents/GitHub/FM4/02-data/PolarPicker_h5_20260308_154910.keras",
}
DATA_PATH   = "/Users/mczhang/Documents/GitHub/FM4/02-data/F_Cl/Filter_Felix_combined.mat"
DL_MAT_PATH = "/Users/mczhang/Documents/GitHub/FM4/02-data/F_Cl/Filter_Felix_combined_DL.mat"
OUT_PNG     = "/Users/mczhang/Documents/GitHub/FM4/03-output-graphics/MyPng/X_DL_polarity_comparison.png"

STATIONS = ["AS1", "AS2", "CC1", "EC1", "EC2", "EC3", "ID1"]


# ── Helpers ───────────────────────────────────────────────────────────────────
def preprocess_waveform(w) -> np.ndarray | None:
    """
    Convert a 250-sample 200 Hz waveform (−0.25 → +1.0 s) to the
    64-sample 100 Hz window (−0.32 → +0.31 s) expected by the DL model.

    1. Downsample 200 Hz → 100 Hz (factor 2, zero-phase IIR anti-alias)
       → 125 samples, −0.25 → +1.0 s
    2. Pre-pad with first 7 samples (−0.25 → −0.19 s) to extend start to −0.32 s
    3. Slice first 64 samples → −0.32 → +0.31 s
    """
    w = np.asarray(w).reshape(-1).astype(np.float64)
    if w.size != 250 or not np.all(np.isfinite(w)):
        return None
    w_100 = ss.decimate(w, q=2, ftype="iir", zero_phase=True)   # 125 samples
    w_pad = np.concatenate([w_100[:7], w_100])                    # 132 samples, from −0.32 s
    return w_pad[:64].astype(np.float32)                          # 64 samples


def normalize(X: np.ndarray) -> np.ndarray:
    """Normalize each row by its max absolute value (floor at 1)."""
    norms = np.max(np.abs(X), axis=1, keepdims=True)
    norms = np.maximum(norms, 1.0)
    return X / norms


# ── Load data ─────────────────────────────────────────────────────────────────
print("Loading data ...")
mat = sio.loadmat(DATA_PATH, simplify_cells=True)
events = mat["Po_Clu"]
if isinstance(events, np.ndarray):
    events = list(events.ravel())
n_events = len(events)
print(f"  {n_events} event clusters loaded.")

# Extract event IDs (used for the saved DL mat)
event_ids = np.array([
    float(np.asarray(ev.get("ID", i)).ravel()[0])
    for i, ev in enumerate(events)
], dtype=np.float64)

# ── Collect waveforms + ground truth per station ───────────────────────────────
# ev_idx tracks which event each wave belongs to (for saving predictions)
per_sta = {sta: {"waves": [], "gt": [], "ev_idx": []} for sta in STATIONS}

for ev_idx, ev in enumerate(events):
    for sta in STATIONS:
        wkey = f"W_{sta}"
        pkey = f"Po_{sta}"

        if pkey not in ev or wkey not in ev:
            continue

        gt_raw = float(np.asarray(ev[pkey]).ravel()[0])
        if gt_raw == 0:
            continue                            # no pick — skip

        w = preprocess_waveform(ev[wkey])
        if w is None:
            continue                            # bad waveform — skip

        per_sta[sta]["waves"].append(w)
        per_sta[sta]["gt"].append(1 if gt_raw > 0 else -1)
        per_sta[sta]["ev_idx"].append(ev_idx)


# ── Model runner ─────────────────────────────────────────────────────────────
def run_model(model_path, per_sta):
    """Load model, predict per station. Returns (rows, preds_by_sta).
    preds_by_sta[sta] is a 1-D array aligned to per_sta[sta]['waves']."""
    print(f"Loading model: {model_path}")
    model = keras.models.load_model(model_path)
    print("  Model loaded.")

    rows         = []
    preds_by_sta = {}

    for sta in STATIONS:
        waves = per_sta[sta]["waves"]
        gts   = per_sta[sta]["gt"]

        if len(waves) == 0:
            preds_by_sta[sta] = np.array([], dtype=np.int8)
            continue

        X   = np.stack(waves, axis=0)
        Xn  = normalize(X).astype(np.float32)
        Xn3 = Xn[:, :, np.newaxis]

        y_raw = model.predict(Xn3, verbose=0)
        if isinstance(y_raw, (list, tuple)):
            y_prob = np.asarray(y_raw[1], dtype=np.float64)
        else:
            y_prob = np.asarray(y_raw, dtype=np.float64)

        idx    = np.argmax(y_prob, axis=1)
        pred   = np.where(idx == 0, -1, 1)
        conf   = np.max(y_prob, axis=1)          # max class probability
        preds_by_sta[sta] = pred

        gt_arr = np.array(gts)
        acc    = (pred == gt_arr).mean() * 100
        TP = int(np.sum((pred ==  1) & (gt_arr ==  1)))
        TN = int(np.sum((pred == -1) & (gt_arr == -1)))
        FP = int(np.sum((pred ==  1) & (gt_arr == -1)))
        FN = int(np.sum((pred == -1) & (gt_arr ==  1)))

        # Shannon entropy (nats, 2-class)
        eps     = 1e-12
        entropy = -(y_prob * np.log(y_prob + eps)).sum(axis=1)

        # 70 % confidence subset
        mask70   = conf >= 0.70
        n70      = int(mask70.sum())
        acc70    = (pred[mask70] == gt_arr[mask70]).mean() * 100 if n70 > 0 else float("nan")
        pct70    = n70 / len(waves) * 100

        # 80 % confidence subset
        mask80   = conf >= 0.80
        n80      = int(mask80.sum())
        acc80    = (pred[mask80] == gt_arr[mask80]).mean() * 100 if n80 > 0 else float("nan")
        pct80    = n80 / len(waves) * 100

        # entropy ≤ 0.2 subset
        maske    = entropy <= 0.20
        ne       = int(maske.sum())
        acce     = (pred[maske] == gt_arr[maske]).mean() * 100 if ne > 0 else float("nan")
        pcte     = ne / len(waves) * 100

        rows.append((sta, len(waves), acc, TP, TN, FP, FN, gt_arr, pred,
                     n70, pct70, acc70,
                     n80, pct80, acc80,
                     ne,  pcte,  acce))

    return rows, preds_by_sta


def _fmt(v, w=8):
    return f"{v:>{w}.1f}" if not np.isnan(v) else f"{'—':>{w}}"

def print_results(label, rows):
    print(f"\n── {label} ──")
    hdr = (f"{'Station':>7}  {'N':>5}  {'Acc%':>6}"
           f"  {'N≥70%':>6}  {'%kpt':>5}  {'Acc≥70':>7}"
           f"  {'N≥80%':>6}  {'%kpt':>5}  {'Acc≥80':>7}"
           f"  {'NH≤.2':>6}  {'%kpt':>5}  {'AccH≤.2':>8}")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        print(f"{r[0]:>7}  {r[1]:>5}  {r[2]:>6.1f}"
              f"  {r[9]:>6d}  {r[10]:>5.1f}  {_fmt(r[11],7)}"
              f"  {r[12]:>6d}  {r[13]:>5.1f}  {_fmt(r[14],7)}"
              f"  {r[15]:>6d}  {r[16]:>5.1f}  {_fmt(r[17],8)}")
    overall   = np.mean([r[2]  for r in rows])
    ov70      = np.nanmean([r[11] for r in rows])
    ov80      = np.nanmean([r[14] for r in rows])
    ove       = np.nanmean([r[17] for r in rows])
    total_n   = sum(r[1]  for r in rows)
    tot70     = sum(r[9]  for r in rows)
    tot80     = sum(r[12] for r in rows)
    tote      = sum(r[15] for r in rows)
    print("-" * len(hdr))
    print(f"{'Overall':>7}  {total_n:>5}  {overall:>6.1f}"
          f"  {tot70:>6}  {tot70/total_n*100:>5.1f}  {_fmt(ov70,7)}"
          f"  {tot80:>6}  {tot80/total_n*100:>5.1f}  {_fmt(ov80,7)}"
          f"  {tote:>6}  {tote/total_n*100:>5.1f}  {_fmt(ove,8)}")
    return overall


# ── Run both models ───────────────────────────────────────────────────────────
all_rows     = {}
all_preds_by_sta = {}
for label, path in MODELS.items():
    rows, preds_by_sta = run_model(path, per_sta)
    if not rows:
        print(f"No predictions for {label} — check W_<sta> fields.")
        raise SystemExit(1)
    all_rows[label]         = rows
    all_preds_by_sta[label] = preds_by_sta
    print_results(label, rows)


# ── Save TMSF_001 predictions for MATLAB ──────────────────────────────────────
print(f"\nSaving TMSF_001 DL predictions → {DL_MAT_PATH}")
tmsf_preds = all_preds_by_sta["TMSF_001"]

save_dict = {"ids": event_ids}
for sta in STATIONS:
    arr        = np.zeros(n_events, dtype=np.float64)   # 0 = no prediction
    ev_indices = per_sta[sta]["ev_idx"]
    preds      = tmsf_preds[sta]
    for ev_i, p in zip(ev_indices, preds):
        arr[ev_i] = float(p)
    save_dict[f"Po_{sta}"] = arr

sio.savemat(DL_MAT_PATH, save_dict, do_compression=True)
n_pred = sum(int(np.sum(save_dict[f"Po_{sta}"] != 0)) for sta in STATIONS)
print(f"  Saved {n_pred} station-picks across {n_events} events.")


# ── Figure: per-station accuracy bars + pooled confusion matrices ─────────────
n_models = len(MODELS)
fig, axes = plt.subplots(1, n_models + 1, figsize=(6 * (n_models + 1), 5))

labels_list = list(MODELS.keys())
colors      = ["steelblue", "darkorange"]

# Left panel — per-station accuracy comparison
ax0 = axes[0]
x   = np.arange(len(STATIONS))
w   = 0.35
for i, label in enumerate(labels_list):
    rows    = all_rows[label]
    accs    = [next((r[2]  for r in rows if r[0] == s), 0)           for s in STATIONS]
    accs70  = [next((r[11] for r in rows if r[0] == s), float("nan")) for s in STATIONS]
    overall = np.mean([r[2] for r in rows])
    overall70 = np.nanmean([r[11] for r in rows])
    offset  = (i - (n_models - 1) / 2) * w
    bars = ax0.bar(x + offset, accs, width=w, color=colors[i],
                   edgecolor="black", label=f"{label} (all: {overall:.1f}%, ≥70% conf: {overall70:.1f}%)")
    for bar, acc, acc70 in zip(bars, accs, accs70):
        ax0.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                 f"{acc:.0f}", ha="center", va="bottom", fontsize=7)
        if not np.isnan(acc70):
            ax0.text(bar.get_x() + bar.get_width() / 2, bar.get_height() - 8,
                     f"({acc70:.0f})", ha="center", va="bottom", fontsize=6,
                     color="white", fontweight="bold")

ax0.set_xticks(x)
ax0.set_xticklabels(STATIONS)
ax0.set_ylim(0, 108)
ax0.set_ylabel("Accuracy (%)")
ax0.set_title("Per-station accuracy vs ground truth\n(FM4 Filter_Felix_combined)")
ax0.legend(fontsize=8)
ax0.grid(axis="y", alpha=0.3)

# Right panels — pooled confusion matrix per model
for i, label in enumerate(labels_list):
    rows     = all_rows[label]
    all_gt   = np.concatenate([r[7] for r in rows])
    all_pred = np.concatenate([r[8] for r in rows])

    cm = np.array([
        [np.sum((all_pred == -1) & (all_gt == -1)),
         np.sum((all_pred ==  1) & (all_gt == -1))],
        [np.sum((all_pred == -1) & (all_gt ==  1)),
         np.sum((all_pred ==  1) & (all_gt ==  1))],
    ])
    cm_pct = cm.astype(float) / cm.sum(axis=1, keepdims=True) * 100

    ax = axes[i + 1]
    im = ax.imshow(cm_pct, cmap="Blues", vmin=0, vmax=100)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Pred Down", "Pred Up"], fontsize=10)
    ax.set_yticks([0, 1]); ax.set_yticklabels(["GT Down",   "GT Up"],   fontsize=10)
    overall = np.mean([r[2] for r in rows])
    ax.set_title(f"Confusion Matrix — {label}\n(all stations pooled, mean acc {overall:.1f}%)")
    plt.colorbar(im, ax=ax, label="%")

    cell_labels = [["TN", "FP"], ["FN", "TP"]]
    for ii in range(2):
        for jj in range(2):
            color = "white" if cm_pct[ii, jj] > 55 else "black"
            ax.text(jj, ii,
                    f"{cell_labels[ii][jj]}\n{cm_pct[ii,jj]:.1f}%\n(n={cm[ii,jj]})",
                    ha="center", va="center", fontsize=10, color=color)

plt.tight_layout()
os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)
plt.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
print(f"Figure saved: {OUT_PNG}")
