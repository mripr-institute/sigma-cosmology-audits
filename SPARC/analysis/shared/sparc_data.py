"""
Shared SPARC catalog data loading utilities for audit custody procedures.

Provides the MRT table parser and morphological classification used by
the SPARC rotation and RAR audit procedures.

Data source: SPARC (Lelli, McGaugh & Schombert 2016).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


# ── Paths ────────────────────────────────────────────────────────────

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR     = PROJECT_ROOT / "data" / "sparc"
TABLE_PATH   = DATA_DIR / "SPARC_Lelli2016c.mrt"
ROTCURVE_DIR = DATA_DIR / "rotcurves"


# ── Morphological type labels ────────────────────────────────────────

HUBBLE_LABELS = {
    0: "S0", 1: "Sa", 2: "Sab", 3: "Sb", 4: "Sbc", 5: "Sc",
    6: "Scd", 7: "Sd", 8: "Sdm", 9: "Sm", 10: "Im", 11: "BCD",
}


def type_class(T):
    """Broad morphological category for colouring.

    Parameters
    ----------
    T : int or None
        Hubble type index from the SPARC catalog.

    Returns
    -------
    str
        One of 'Early Spiral', 'Spiral', 'Late-type', 'Irregular', 'Unknown'.
    """
    if T is None:
        return "Unknown"
    if T <= 3:
        return "Early Spiral"
    if T <= 6:
        return "Spiral"
    if T <= 9:
        return "Late-type"
    return "Irregular"


TYPE_COLORS = {
    "Early Spiral": "#e07cc5",
    "Spiral":       "#4cc9f0",
    "Late-type":    "#7ef08a",
    "Irregular":    "#ff6b6b",
    "Unknown":      "#888888",
}


# ── MRT table parser ─────────────────────────────────────────────────
# Whitespace-split column indices (0-based):
#   0: Galaxy     1: T      2: D      3: e_D     4: f_D
#   5: Inc        6: e_Inc  7: L[3.6] 8: e_L     9: Reff
#  10: SBeff     11: Rdisk 12: SBdisk 13: MHI    14: RHI
#  15: Vflat     16: e_Vflat  17: Q   18+: Ref


def _parse_float(s):
    """Parse a string as float, returning 0.0 on failure."""
    try:
        return float(s)
    except (ValueError, IndexError):
        return 0.0


def _parse_int(s, default=None):
    """Parse a string as int, returning *default* on failure."""
    try:
        return int(s)
    except (ValueError, IndexError):
        return default


def load_sparc_table():
    """Parse SPARC_Lelli2016c.mrt via whitespace splitting.

    Data lines start after the last ``---`` separator in the file.

    Returns
    -------
    list[dict]
        Each dict has keys: name, T, D_Mpc, Reff_kpc, Rdisk_kpc,
        Vflat_kms, e_Vflat, Q.
    """
    if not TABLE_PATH.exists():
        print(f"ERROR: {TABLE_PATH} not found.")
        print("Run:  python3 sparc/download_sparc.py")
        sys.exit(1)

    lines = TABLE_PATH.read_text().splitlines()

    # Data starts after the LAST '---' separator
    last_sep = max(i for i, l in enumerate(lines) if l.startswith("---"))

    galaxies = []
    for line in lines[last_sep + 1:]:
        fields = line.split()
        if len(fields) < 18:
            continue

        galaxies.append({
            "name":      fields[0],
            "T":         _parse_int(fields[1]),
            "D_Mpc":     _parse_float(fields[2]),
            "Reff_kpc":  _parse_float(fields[9]),
            "Rdisk_kpc": _parse_float(fields[11]),
            "Vflat_kms": _parse_float(fields[15]),
            "e_Vflat":   _parse_float(fields[16]),
            "Q":         _parse_int(fields[17], 3),
        })

    return galaxies


def load_rotcurve(name: str) -> dict | None:
    """Load rotation curve for a single galaxy (Vobs, errV only).

    Parameters
    ----------
    name : str
        Galaxy name from the SPARC catalog.

    Returns
    -------
    dict or None
        Dict with r_kpc, Vobs, errV arrays, or None if unavailable.
    """
    path = ROTCURVE_DIR / f"{name}_rotmod.dat"
    if not path.exists():
        return None
    try:
        data = np.loadtxt(path)
        if data.ndim < 2 or data.shape[1] < 3:
            return None
        return {
            "r_kpc": data[:, 0],
            "Vobs":  data[:, 1],
            "errV":  data[:, 2],
        }
    except Exception:
        return None


def load_rotcurve_full(name: str) -> dict | None:
    """Load rotation curve with ALL baryonic velocity components.

    Parameters
    ----------
    name : str
        Galaxy name from the SPARC catalog.

    Returns
    -------
    dict or None
        Dict with arrays: r_kpc, Vobs, errV, Vgas, Vdisk, Vbul,
        or None if file missing / unreadable.
    """
    path = ROTCURVE_DIR / f"{name}_rotmod.dat"
    if not path.exists():
        return None
    try:
        data = np.loadtxt(path)
        if data.ndim < 2 or data.shape[1] < 6:
            return None
        return {
            "r_kpc": data[:, 0],
            "Vobs":  data[:, 1],
            "errV":  data[:, 2],
            "Vgas":  data[:, 3],
            "Vdisk": data[:, 4],
            "Vbul":  data[:, 5],
        }
    except Exception:
        return None
