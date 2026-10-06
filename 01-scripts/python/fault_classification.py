"""
fault_classification.py
=======================
Fault-type classification helpers for focal mechanism plotting.

Implements the same P/T/B-axis plunge thresholding used in the MATLAB pipeline
(e.g., `01-scripts/G_FM_HASH_22OBS.m` and `01-scripts/H_read_SKHASH.m`).

Labels:
  'N' = Normal
  'R' = Reverse
  'S' = Strike-slip
  'U' = Oblique/Unknown
"""

from __future__ import annotations

import numpy as np


def _plunge_deg(v: np.ndarray) -> float:
    """Return plunge (0–90°) from the vertical (Down) component."""
    n = float(np.linalg.norm(v))
    if not np.isfinite(n) or n == 0.0:
        return float("nan")
    z = float(abs(v[2]) / n)
    z = float(np.clip(z, 0.0, 1.0))
    return float(np.degrees(np.arcsin(z)))


def _sdr_to_n_s(strike_deg: float, dip_deg: float, rake_deg: float) -> tuple[np.ndarray, np.ndarray]:
    """
    Convert strike/dip/rake (degrees) to unit fault normal and slip vectors.

    Coordinate system: NED (North, East, Down). Strike is clockwise from North.
    Dip follows the right-hand rule (dip direction is strike + 90°).
    Rake follows the common convention: -90 ~ normal, +90 ~ reverse.
    """
    strike = np.radians(float(strike_deg))
    dip = np.radians(float(dip_deg))
    rake = np.radians(float(rake_deg))

    # Unit vectors along strike and down-dip (NED)
    u_strike = np.array([np.cos(strike), np.sin(strike), 0.0], dtype=float)
    u_dip = np.array(
        [-np.sin(strike) * np.cos(dip), np.cos(strike) * np.cos(dip), np.sin(dip)],
        dtype=float,
    )

    # Fault normal (unit) and slip (unit) vectors
    normal = np.cross(u_strike, u_dip)
    slip = np.cos(rake) * u_strike + np.sin(rake) * u_dip

    # Numerical safety (should already be unit vectors)
    normal /= np.linalg.norm(normal)
    slip /= np.linalg.norm(slip)
    return normal, slip


def ptb_plunges_deg(strike_deg: float, dip_deg: float, rake_deg: float) -> tuple[float, float, float]:
    """Return (P_plunge, T_plunge, B_plunge) in degrees."""
    normal, slip = _sdr_to_n_s(strike_deg, dip_deg, rake_deg)
    b_axis = np.cross(normal, slip)
    t_axis = (normal + slip) / np.sqrt(2.0)
    p_axis = (normal - slip) / np.sqrt(2.0)
    return _plunge_deg(p_axis), _plunge_deg(t_axis), _plunge_deg(b_axis)


def classify_fault_ptb(strike_deg: float, dip_deg: float, rake_deg: float) -> str:
    """
    Classify fault type using P/T/B-axis plunge thresholds (MATLAB-equivalent).

    Thresholds match the MATLAB code:
      N: P>=52 and T<=35
      R: P<=35 and T>=52
      S: P<=40 and B>=45 and T<=40
      oblique reverse: P<=20 and 40<=T<=52 -> R
      oblique normal:  40<=P<=52 and T<=20 -> N
      else -> U
    """
    p_l, t_l, b_l = ptb_plunges_deg(strike_deg, dip_deg, rake_deg)
    if not (np.isfinite(p_l) and np.isfinite(t_l) and np.isfinite(b_l)):
        return "U"

    if p_l >= 52 and t_l <= 35:
        return "N"
    if p_l <= 35 and t_l >= 52:
        return "R"
    if p_l <= 40 and b_l >= 45 and t_l <= 40:
        return "S"
    if p_l <= 20 and 40 <= t_l <= 52:
        return "R"
    if 40 <= p_l <= 52 and t_l <= 20:
        return "N"
    return "U"
