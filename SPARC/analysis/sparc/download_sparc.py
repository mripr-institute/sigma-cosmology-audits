#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
download_sparc.py
=================

Download the full SPARC catalog and rotation curves.

Source:
  SPARC database (Lelli, McGaugh, Schombert 2016).
  Hosted at Case Western Reserve University.

Downloads to data/sparc/:
  SPARC_Lelli2016c.mrt     Galaxy properties table (175 galaxies)
  rotcurves/               Individual rotation curve data files

Usage:
    python3 sparc/download_sparc.py
"""

from __future__ import annotations

import os
import sys
import zipfile
from io import BytesIO
from pathlib import Path
from urllib.request import urlretrieve, urlopen


# ── Paths ────────────────────────────────────────────────────────────

_ROOT  = Path(__file__).resolve().parents[2]
OUTDIR = _ROOT / "data" / "sparc"
ROTDIR = OUTDIR / "rotcurves"


# ── SPARC URLs ───────────────────────────────────────────────────────

_BASE = "https://astroweb.cwru.edu/SPARC"

TABLE_FILE     = "SPARC_Lelli2016c.mrt"
TABLE_URL      = f"{_BASE}/{TABLE_FILE}"
TABLE_MIN_SIZE = 10_000   # ~30 KB expected

ROTMOD_ZIP_URL     = f"{_BASE}/Rotmod_LTG.zip"
ROTMOD_MIN_FILES   = 100  # expect ~175


# ── Helpers ──────────────────────────────────────────────────────────

def _progress(block_num: int, block_size: int, total_size: int) -> None:
    downloaded = block_num * block_size
    if total_size > 0:
        pct = min(100.0, downloaded / total_size * 100.0)
        print(
            f"\r  {downloaded / 1e3:,.0f} / {total_size / 1e3:,.0f} KB"
            f"  ({pct:.0f}%)",
            end="", flush=True,
        )
    else:
        print(f"\r  {downloaded / 1e3:,.0f} KB", end="", flush=True)


def _valid(path: Path, min_bytes: int) -> bool:
    return path.exists() and path.stat().st_size >= min_bytes


# ── Download: summary table ─────────────────────────────────────────

def download_table() -> bool:
    dest = OUTDIR / TABLE_FILE
    if _valid(dest, TABLE_MIN_SIZE):
        print(f"VALID  {TABLE_FILE}  ({dest.stat().st_size:,} B)")
        return True

    if dest.exists():
        dest.unlink()

    print(f"Downloading {TABLE_FILE}")
    print(f"  {TABLE_URL}")
    try:
        urlretrieve(TABLE_URL, dest, reporthook=_progress)
        print()
    except Exception as e:
        print(f"\n  FAILED: {e}")
        if dest.exists():
            dest.unlink()
        return False

    if _valid(dest, TABLE_MIN_SIZE):
        print(f"  OK  ({dest.stat().st_size:,} B)")
        return True
    return False


# ── Download: rotation curves ────────────────────────────────────────

def download_rotcurves() -> bool:
    ROTDIR.mkdir(parents=True, exist_ok=True)

    existing = list(ROTDIR.glob("*_rotmod.dat"))
    if len(existing) >= ROTMOD_MIN_FILES:
        print(f"VALID  rotcurves/  ({len(existing)} files)")
        return True

    print(f"Downloading rotation curves archive")
    print(f"  {ROTMOD_ZIP_URL}")
    try:
        data = urlopen(ROTMOD_ZIP_URL).read()
        print(f"  {len(data) / 1e6:.1f} MB downloaded")
    except Exception as e:
        print(f"  FAILED: {e}")
        return False

    try:
        with zipfile.ZipFile(BytesIO(data)) as zf:
            members = [m for m in zf.namelist()
                       if m.endswith("_rotmod.dat")]
            for m in members:
                dest = ROTDIR / os.path.basename(m)
                dest.write_bytes(zf.read(m))
            print(f"  Extracted {len(members)} rotation curve files")
            return len(members) >= ROTMOD_MIN_FILES
    except Exception as e:
        print(f"  EXTRACTION FAILED: {e}")
        return False


# ── Main ─────────────────────────────────────────────────────────────

def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    print(f"Output: {OUTDIR.resolve()}\n")

    table_ok = download_table()
    rot_ok   = download_rotcurves()

    # ── Summary ──────────────────────────────────────────────────────
    print("\n" + "=" * 56)
    print("SPARC Download Summary")
    print("=" * 56)

    tp = OUTDIR / TABLE_FILE
    if tp.exists():
        print(f"  [{'OK':>10s}]  {TABLE_FILE}  ({tp.stat().st_size:,} B)")
    else:
        print(f"  [{'MISSING':>10s}]  {TABLE_FILE}")

    n_rot = len(list(ROTDIR.glob("*_rotmod.dat")))
    tag = "OK" if n_rot >= ROTMOD_MIN_FILES else "INCOMPLETE"
    print(f"  [{tag:>10s}]  rotcurves/  ({n_rot} files)")

    if table_ok and rot_ok:
        print("\nAll SPARC data ready.")
    else:
        print("\nSome downloads failed. Re-run this script.")
        sys.exit(1)


if __name__ == "__main__":
    main()
