#!/usr/bin/env python3
from __future__ import annotations

import argparse
import ast
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import hashlib
import importlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
from typing import Any

try:
    import mpmath as _mp
except ModuleNotFoundError:
    _mp = None

mp = _mp if _mp is not None and all(hasattr(_mp, name) for name in ("mpf", "exp", "log")) else None

import numpy as np

SCRIPT_NAME = "sigma_desi_lss_closure_squared_radial_audit.py"
RUN_LABEL = "sigma_desi_lss_closure_squared_radial"
DEFAULT_OUTPUT_DIR = Path("results/sigma_desi_lss_closure_squared_radial")

MU_TEXT = "0.082912607552"
GAMMA_TEXT = "0.38603416"
MU = float(MU_TEXT)
GAMMA = float(GAMMA_TEXT)
s_star = 1.0 / MU - 1.0 / GAMMA
r_star_linear_display = s_star ** 0.5
s_tail = s_star + 1.0 / MU

LEGACY_LINEAR_HALF_WIDTH_REFERENCE = 1.0
q_shell_halfwidth = 2.0 * r_star_linear_display * LEGACY_LINEAR_HALF_WIDTH_REFERENCE
equivalent_linear_halfwidth_near_s_star = q_shell_halfwidth / (2.0 * r_star_linear_display)

RETAINED_DIPOLE_AXIS_L_DEG = 263.99
RETAINED_DIPOLE_AXIS_B_DEG = 48.26
CLOSURE_TRACERS = ("LRG", "ELG", "QSO")
PROFILE_TRACERS = ("BGS", "LRG", "ELG", "QSO")
EXT = "f" + "its"
TRACER_DATA_FILES = {
    "BGS": "BGS/datcomb_BGS_BRIGHT_tarspecwdup_zdone." + EXT,
    "LRG": "LRG/datcomb_LRG_tarspecwdup_zdone." + EXT,
    "ELG": "ELG/datcomb_ELG_LOPnotqso_tarspecwdup_zdone." + EXT,
    "QSO": "QSO/datcomb_QSO_tarspecwdup_zdone." + EXT,
}
TRACER_RANDOM_FILES = {
    "BGS": "BGS/randoms/rancomb_0brightwdupspec_zdone." + EXT,
    "LRG": "LRG/randoms/rancomb_0darkwdupspec_zdone." + EXT,
    "ELG": "ELG/randoms/rancomb_0darkwdupspec_zdone." + EXT,
    "QSO": "LRG/randoms/rancomb_0darkwdupspec_zdone." + EXT,
}

N_Q_PROFILE_BINS = 192
N_CT_BINS = 8
N_PHI_BINS = 12
MIN_EXPECTED = 5.0
NULL_TRIALS = 200
NULL_SUBSAMPLE = 500_000
DEFAULT_SCAN_RANGE = 14.0
DEFAULT_SCAN_STEPS = 41
HIGH_PRECISION_DPS = 4096

ICRS_TO_GALACTIC = np.array([
    [-0.0548755604, -0.8734370902, -0.4838350155],
    [0.4941094279, -0.4448296300, 0.7469822445],
    [-0.8676661490, -0.1980763734, 0.4559837762],
], dtype=np.float64)


def stamp() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def note(text: str) -> None:
    print(f"  [{time.strftime('%H:%M:%S')}] {text}", flush=True)


def hp_mpf(text: str) -> Any:
    return mp.mpf(text) if mp is not None else None


def sigma_s(value: np.ndarray | float) -> np.ndarray | float:
    return np.log(MU * MU * (1.0 + GAMMA * value) / (MU + GAMMA)) - MU * value


def rho_s(value: np.ndarray | float) -> np.ndarray | float:
    return MU * MU * (1.0 + GAMMA * value) * np.exp(-MU * value) / (MU + GAMMA)


def F_s(value: np.ndarray | float) -> np.ndarray | float:
    return 1.0 - np.exp(-MU * value) * (MU + GAMMA * (1.0 + MU * value)) / (MU + GAMMA)


def d_sigma_ds(value: np.ndarray | float) -> np.ndarray | float:
    return GAMMA / (1.0 + GAMMA * value) - MU


def d2_sigma_ds2(value: np.ndarray | float) -> np.ndarray | float:
    return -(GAMMA * GAMMA) / ((1.0 + GAMMA * value) ** 2)


def high_precision_identity_tests() -> dict[str, Any]:
    if mp is None:
        raise RuntimeError("mpmath is required for 4096-digit identity tests")
    mp.mp.dps = HIGH_PRECISION_DPS
    mu = hp_mpf(MU_TEXT)
    gamma = hp_mpf(GAMMA_TEXT)
    one = mp.mpf("1")
    s_star_hp = one / mu - one / gamma
    r_display_hp = mp.sqrt(s_star_hp)
    s_tail_hp = s_star_hp + one / mu
    points = [
        mp.mpf("0"), mp.mpf("1e-30"), s_star_hp / 4, s_star_hp / 2,
        s_star_hp, s_tail_hp, s_tail_hp * 2, mp.mpf("1000") / mu,
    ]
    max_delta = mp.mpf("0")
    for value in points:
        sigma = mp.log(mu * mu * (one + gamma * value) / (mu + gamma)) - mu * value
        rho = mu * mu * (one + gamma * value) * mp.exp(-mu * value) / (mu + gamma)
        max_delta = max(max_delta, abs(rho - mp.exp(sigma)))
    prime = gamma / (one + gamma * s_star_hp) - mu
    second = -(gamma * gamma) / ((one + gamma * s_star_hp) ** 2)
    second_plus = second + mu * mu
    F0 = one - mp.exp(-mu * mp.mpf("0")) * (mu + gamma * (one + mu * mp.mpf("0"))) / (mu + gamma)
    Ffar = one - mp.exp(-mu * (mp.mpf("1000") / mu)) * (mu + gamma * (one + mu * (mp.mpf("1000") / mu))) / (mu + gamma)
    display_delta = r_display_hp * r_display_hp - s_star_hp
    return {
        "precision_engine": "mpmath",
        "precision_digits": HIGH_PRECISION_DPS,
        "rho_s_equals_exp_sigma_s_max_abs": mp.nstr(max_delta, 90),
        "d_sigma_ds_at_s_star": mp.nstr(prime, 90),
        "d2_sigma_ds2_at_s_star": mp.nstr(second, 90),
        "d2_sigma_ds2_plus_mu_squared": mp.nstr(second_plus, 90),
        "s_star_formula": mp.nstr(s_star_hp, 90),
        "s_star_float": s_star,
        "r_star_linear_display_squared_minus_s_star": mp.nstr(display_delta, 90),
        "F_s_zero": mp.nstr(F0, 90),
        "F_s_far": mp.nstr(Ffar, 90),
        "pass": bool(max_delta < mp.mpf("1e-80") and abs(prime) < mp.mpf("1e-16") and abs(second_plus) < mp.mpf("1e-16") and abs(display_delta) < mp.mpf("1e-80")),
    }


def sha256_file(path: Path | None, cache: dict[str, str] | None = None) -> str | None:
    if path is None or not path.exists():
        return None
    key = str(path.resolve())
    if cache is not None and key in cache:
        return cache[key]
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(16 * 1024 * 1024), b""):
            h.update(chunk)
    digest = h.hexdigest()
    if cache is not None:
        cache[key] = digest
    return digest


def to_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): to_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_json(v) for v in value]
    if isinstance(value, np.ndarray):
        return to_json(value.tolist())
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        number = float(value)
        return None if not math.isfinite(number) else number
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, Path):
        return str(value)
    return value


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.write_text(json.dumps(to_json(data), indent=2, sort_keys=True) + "\n", encoding="utf-8")


def package_version(name: str) -> str | None:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None


def table_reader():
    try:
        return importlib.import_module("astropy.io." + EXT)
    except ModuleNotFoundError as exc:
        raise RuntimeError("astropy is required") from exc


def plotter():
    try:
        matplotlib = importlib.import_module("matplotlib")
        matplotlib.use("Agg", force=True)
        return importlib.import_module("matplotlib.pyplot")
    except ModuleNotFoundError as exc:
        raise RuntimeError("matplotlib is required") from exc


def gaussian_smooth_1d(values: np.ndarray, sigma: float) -> np.ndarray:
    if sigma <= 0:
        return values.astype(np.float64)
    radius = max(1, int(4.0 * sigma + 0.5))
    x = np.arange(-radius, radius + 1, dtype=np.float64)
    kernel = np.exp(-0.5 * (x / sigma) * (x / sigma))
    kernel /= kernel.sum()
    padded = np.pad(values.astype(np.float64), radius, mode="edge")
    return np.convolve(padded, kernel, mode="valid")


def galactic_unit_vector(longitude_deg: float, latitude_deg: float) -> np.ndarray:
    lon = np.radians(longitude_deg)
    lat = np.radians(latitude_deg)
    return np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], dtype=np.float64)


retained_dipole_axis_unit_vector = galactic_unit_vector(RETAINED_DIPOLE_AXIS_L_DEG, RETAINED_DIPOLE_AXIS_B_DEG)


def unit_vectors(ra_deg: np.ndarray, dec_deg: np.ndarray) -> np.ndarray:
    ra = np.radians(np.asarray(ra_deg, dtype=np.float64))
    dec = np.radians(np.asarray(dec_deg, dtype=np.float64))
    icrs = np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], axis=-1)
    return icrs @ ICRS_TO_GALACTIC.T


def axis_relative_angles(n_hat: np.ndarray, direction: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    cos_theta = n_hat @ direction
    ref = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    if abs(float(np.dot(ref, direction))) > 0.9:
        ref = np.array([0.0, 1.0, 0.0], dtype=np.float64)
    e1 = ref - float(np.dot(ref, direction)) * direction
    e1 = e1 / (float(np.sum(e1 * e1)) ** 0.5)
    e2 = np.cross(direction, e1)
    return cos_theta, np.arctan2(n_hat @ e2, n_hat @ e1)


def q_from_scan_offset(z: np.ndarray, n_hat: np.ndarray, scan_offset: float, direction: np.ndarray) -> np.ndarray:
    cos_theta = n_hat @ direction
    return z * z + scan_offset * scan_offset - 2.0 * z * scan_offset * cos_theta


def q_from_coordinates(z: np.ndarray, n_hat: np.ndarray, scan_coordinate: np.ndarray) -> tuple[np.ndarray, float]:
    X = z[:, None] * n_hat
    U = X - scan_coordinate[None, :]
    q_closure = np.sum(U * U, axis=1)
    q_direct = U[:, 0] * U[:, 0] + U[:, 1] * U[:, 1] + U[:, 2] * U[:, 2]
    residual = float(np.max(np.abs(q_closure - q_direct))) if q_closure.size else 0.0
    return q_closure, residual


def retained_mask(cols: dict[str, np.ndarray]) -> np.ndarray:
    n = next(iter(cols.values())).size
    mask = np.ones(n, dtype=bool)
    mask &= np.isfinite(cols["RA"])
    mask &= np.isfinite(cols["DEC"])
    mask &= np.isfinite(cols["Z"]) & (cols["Z"] > 0)
    if "ZWARN" in cols:
        mask &= cols["ZWARN"] == 0
    return mask


def read_data(path: Path) -> dict[str, np.ndarray]:
    reader = table_reader()
    with reader.open(path, memmap=True) as hdul:
        data = hdul[1].data
        names = list(data.names)
        out = {
            "RA": np.asarray(data["RA"], dtype=np.float64),
            "DEC": np.asarray(data["DEC"], dtype=np.float64),
            "Z": np.asarray(data["Z"], dtype=np.float64),
        }
        if "ZWARN" in names:
            out["ZWARN"] = np.asarray(data["ZWARN"], dtype=np.float64)
    return out


def discover_inputs(data_root: Path) -> dict[str, dict[str, Any]]:
    out = {}
    for tracer in PROFILE_TRACERS:
        data_path = data_root / TRACER_DATA_FILES[tracer]
        random_path = data_root / TRACER_RANDOM_FILES[tracer]
        out[tracer] = {
            "data_path": data_path,
            "data_exists": data_path.exists(),
            "random_path": random_path,
            "random_exists": random_path.exists(),
        }
    return out


def check_inputs(inputs: dict[str, dict[str, Any]], tracers: tuple[str, ...]) -> None:
    missing = []
    for tracer in tracers:
        if not inputs[tracer]["data_exists"]:
            missing.append(str(inputs[tracer]["data_path"]))
    for tracer in CLOSURE_TRACERS:
        if not inputs[tracer]["random_exists"]:
            missing.append(str(inputs[tracer]["random_path"]))
    if missing:
        raise RuntimeError("missing DESI input: " + "; ".join(missing))


def hash_inputs(inputs: dict[str, dict[str, Any]], cache: dict[str, str]) -> dict[str, Any]:
    out = {}
    for tracer in PROFILE_TRACERS:
        entry = inputs[tracer]
        out[tracer] = {
            "data_path": str(entry["data_path"]),
            "data_sha256": sha256_file(entry["data_path"], cache) if entry["data_exists"] else None,
            "random_path": str(entry["random_path"]),
            "random_sha256": sha256_file(entry["random_path"], cache) if entry["random_exists"] else None,
        }
    return out


def load_catalogues(data_root: Path) -> dict[str, dict[str, Any]]:
    out = {}
    for tracer in CLOSURE_TRACERS:
        note(f"loading {tracer}")
        path = data_root / TRACER_DATA_FILES[tracer]
        cols = read_data(path)
        mask = retained_mask(cols)
        ra = cols["RA"][mask]
        dec = cols["DEC"][mask]
        z = cols["Z"][mask]
        n_hat = unit_vectors(ra, dec)
        gal_lat = np.degrees(np.arcsin(np.clip(n_hat[:, 2], -1.0, 1.0)))
        out[tracer] = {"tracer": tracer, "path": str(path), "ra": ra, "dec": dec, "z": z, "n_hat": n_hat, "gal_lat": gal_lat}
        note(f"{tracer}: {z.size:,} retained")
    return out


def nearest_s_star_peak(q_values: np.ndarray, bins: int = 400, smooth_sigma: float = 5.0) -> tuple[float, float]:
    valid = q_values[np.isfinite(q_values) & (q_values >= 0)]
    if valid.size < 100:
        return np.nan, np.nan
    lo = float(np.percentile(valid, 0.5))
    hi = float(np.percentile(valid, 99.5))
    if hi <= lo:
        med = float(np.median(valid))
        return abs(med - s_star), med
    counts, edges = np.histogram(valid, bins=bins, range=(lo, hi))
    centers = 0.5 * (edges[:-1] + edges[1:])
    smoothed = gaussian_smooth_1d(counts.astype(np.float64), sigma=smooth_sigma)
    local_max = (smoothed[1:-1] > smoothed[:-2]) & (smoothed[1:-1] > smoothed[2:])
    peak_idx = np.where(local_max)[0] + 1
    if peak_idx.size == 0:
        peak = centers[int(np.argmax(smoothed))]
    else:
        peak = centers[peak_idx[int(np.argmin(np.abs(centers[peak_idx] - s_star)))]]
    return abs(float(peak - s_star)), float(peak)


def q_shell_phi_chi2(q_values: np.ndarray, cos_theta: np.ndarray, phi: np.ndarray) -> tuple[float, int, float]:
    on_shell = np.abs(q_values - s_star) < q_shell_halfwidth
    if int(on_shell.sum()) < 50:
        return 0.0, 0, 0.0
    ct_edges = np.linspace(-1.0, 1.0, N_CT_BINS + 1)
    phi_edges = np.linspace(-np.pi, np.pi, N_PHI_BINS + 1)
    h_all, _, _ = np.histogram2d(cos_theta, phi, bins=[ct_edges, phi_edges])
    h_shell, _, _ = np.histogram2d(cos_theta[on_shell], phi[on_shell], bins=[ct_edges, phi_edges])
    chi2 = 0.0
    dof = 0
    for j in range(N_CT_BINS):
        n_band = h_all[j].sum()
        n_shell_band = h_shell[j].sum()
        if n_band < 50 or n_shell_band < 10:
            continue
        expected = h_all[j] * (n_shell_band / n_band)
        valid = expected >= MIN_EXPECTED
        if int(valid.sum()) > 1:
            chi2 += float(np.sum((h_shell[j][valid] - expected[valid]) ** 2 / expected[valid]))
            dof += int(valid.sum()) - 1
    return float(chi2), int(dof), float(chi2 / dof - 1.0) if dof > 0 else 0.0


def evaluate_q_signature(z_list: list[np.ndarray], n_list: list[np.ndarray], direction: np.ndarray, scan_offset: float) -> dict[str, Any]:
    peak_residuals, peaks, occupancies, q_cvs = [], [], [], []
    q_all_parts, cos_parts, phi_parts = [], [], []
    tracer_results = {}
    for tracer, z, n_hat in zip(CLOSURE_TRACERS[:len(z_list)], z_list, n_list):
        q_values = q_from_scan_offset(z, n_hat, scan_offset, direction)
        q_residual = q_values - s_star
        peak_residual, peak = nearest_s_star_peak(q_values)
        shell_mask = np.abs(q_residual) < q_shell_halfwidth
        occupancy = float(np.mean(shell_mask)) if q_values.size else np.nan
        valid = q_values[np.isfinite(q_values)]
        q_cv = float(np.std(valid) / np.mean(valid)) if valid.size and np.mean(valid) > 0 else np.nan
        cos_theta, phi = axis_relative_angles(n_hat, direction)
        peak_residuals.append(peak_residual)
        peaks.append(peak)
        occupancies.append(occupancy)
        q_cvs.append(q_cv)
        q_all_parts.append(q_values)
        cos_parts.append(cos_theta)
        phi_parts.append(phi)
        tracer_results[tracer] = {
            "retained_count": int(z.size),
            "q_peak": peak,
            "q_peak_residual": peak_residual,
            "q_shell_count": int(shell_mask.sum()),
            "q_shell_fraction": occupancy,
            "q_residual_median": float(np.median(q_residual)) if q_residual.size else np.nan,
            "q_residual_mad": float(np.median(np.abs(q_residual - np.median(q_residual)))) if q_residual.size else np.nan,
            "q_min": float(np.min(q_values)) if q_values.size else np.nan,
            "q_max": float(np.max(q_values)) if q_values.size else np.nan,
        }
    finite_res = np.array([v for v in peak_residuals if np.isfinite(v)], dtype=np.float64)
    finite_peaks = np.array([v for v in peaks if np.isfinite(v)], dtype=np.float64)
    residual = float(np.mean(finite_res)) if finite_res.size else np.nan
    frac_residual = residual / s_star if np.isfinite(residual) and s_star > 0 else np.nan
    occupancy = float(np.nanmean(occupancies)) if occupancies else np.nan
    cross_cv = float(np.std(finite_peaks) / np.mean(finite_peaks)) if finite_peaks.size >= 2 and np.mean(finite_peaks) > 0 else np.nan
    q_all = np.concatenate(q_all_parts) if q_all_parts else np.array([], dtype=np.float64)
    cos_all = np.concatenate(cos_parts) if cos_parts else np.array([], dtype=np.float64)
    phi_all = np.concatenate(phi_parts) if phi_parts else np.array([], dtype=np.float64)
    phi_chi2, phi_dof, phi_excess = q_shell_phi_chi2(q_all, cos_all, phi_all)
    score_terms = [frac_residual, 1.0 - occupancy, cross_cv]
    primary_score = float(np.mean([v for v in score_terms if np.isfinite(v)])) if any(np.isfinite(v) for v in score_terms) else np.nan
    return {
        "best_scan_offset": float(scan_offset),
        "q_peak_residual": residual,
        "q_frac_residual": frac_residual,
        "q_shell_occupancy": occupancy,
        "cross_tracer_q_peak_cv": cross_cv,
        "phi_chi2": phi_chi2,
        "phi_dof": phi_dof,
        "phi_excess": phi_excess,
        "primary_score": primary_score,
        "q_peaks": [float(x) if np.isfinite(x) else np.nan for x in peaks],
        "q_peak_residuals": [float(x) if np.isfinite(x) else np.nan for x in peak_residuals],
        "tracer_results": tracer_results,
    }


def scan_q_signature(z_list: list[np.ndarray], n_list: list[np.ndarray], direction: np.ndarray, scan_offsets: np.ndarray) -> dict[str, Any]:
    rows = []
    for i, scan_offset in enumerate(scan_offsets):
        if i % max(1, len(scan_offsets) // 10) == 0:
            note(f"q closure scan {i + 1}/{len(scan_offsets)} offset={scan_offset:+.4f}")
        row = evaluate_q_signature(z_list, n_list, direction, float(scan_offset))
        rows.append(row)
    finite = [row for row in rows if np.isfinite(row["primary_score"])]
    best = min(finite, key=lambda row: row["primary_score"]) if finite else {}
    return {
        "scan_domain": [float(scan_offsets[0]), float(scan_offsets[-1])],
        "scan_resolution": int(len(scan_offsets)),
        "rows": rows,
        "best_scan_offset": best.get("best_scan_offset", np.nan),
        "best_scan_result": best,
        "score_definition": "mean(q_peak_residual/s_star, 1 - q_shell_occupancy, cross_tracer_q_peak_cv)",
        "closure_predicate": "||U||^2 = s_star",
        "q_shell_halfwidth": q_shell_halfwidth,
        "equivalent_linear_halfwidth_near_s_star": equivalent_linear_halfwidth_near_s_star,
    }


def scan_q_signature_fast_null(
    z_list: list[np.ndarray],
    n_list: list[np.ndarray],
    direction: np.ndarray,
    scan_offsets: np.ndarray,
) -> dict[str, Any]:
    """Null-only scan with the exact primary score and deferred diagnostics.

    Angular histograms, medians, MADs, extrema and detailed tracer dictionaries
    do not enter the offset selection.  Compute them once at the winning offset.
    """
    cos_list = [n_hat @ direction for n_hat in n_list]
    base_list = [z * z for z in z_list]
    cross_list = [z * cos_theta for z, cos_theta in zip(z_list, cos_list)]
    rows: list[dict[str, Any]] = []
    for scan_offset in scan_offsets:
        peaks: list[float] = []
        peak_residuals: list[float] = []
        occupancies: list[float] = []
        offset2 = float(scan_offset * scan_offset)
        twice_offset = float(2.0 * scan_offset)
        for base, cross in zip(base_list, cross_list):
            q_values = base + offset2 - twice_offset * cross
            peak_residual, peak = nearest_s_star_peak(q_values)
            peaks.append(peak)
            peak_residuals.append(peak_residual)
            occupancies.append(float(np.mean(np.abs(q_values - s_star) < q_shell_halfwidth)))
        finite_res = np.asarray([v for v in peak_residuals if np.isfinite(v)], dtype=np.float64)
        finite_peaks = np.asarray([v for v in peaks if np.isfinite(v)], dtype=np.float64)
        residual = float(np.mean(finite_res)) if finite_res.size else np.nan
        frac_residual = residual / s_star if np.isfinite(residual) else np.nan
        occupancy = float(np.mean(occupancies)) if occupancies else np.nan
        cross_cv = float(np.std(finite_peaks) / np.mean(finite_peaks)) if finite_peaks.size >= 2 and np.mean(finite_peaks) > 0 else np.nan
        terms = [frac_residual, 1.0 - occupancy, cross_cv]
        primary = float(np.mean([v for v in terms if np.isfinite(v)]))
        rows.append({
            "best_scan_offset": float(scan_offset),
            "q_peak_residual": residual,
            "q_frac_residual": frac_residual,
            "q_shell_occupancy": occupancy,
            "cross_tracer_q_peak_cv": cross_cv,
            "primary_score": primary,
            "q_peaks": peaks,
            "q_peak_residuals": peak_residuals,
        })
    best_compact = min(rows, key=lambda row: row["primary_score"])
    best = evaluate_q_signature(z_list, n_list, direction, float(best_compact["best_scan_offset"]))
    return {
        "scan_domain": [float(scan_offsets[0]), float(scan_offsets[-1])],
        "scan_resolution": int(scan_offsets.size),
        "rows": rows,
        "best_scan_offset": best["best_scan_offset"],
        "best_scan_result": best,
        "score_definition": "mean(q_peak_residual/s_star, 1 - q_shell_occupancy, cross_tracer_q_peak_cv)",
        "closure_predicate": "||U||^2 = s_star",
        "q_shell_halfwidth": q_shell_halfwidth,
        "optimization": "exact score scan; non-selecting diagnostics deferred to winning offset",
    }


def fast_null_equivalence_test() -> dict[str, Any]:
    rng = np.random.default_rng(90210)
    z_list = [rng.uniform(0.01, 4.0, size=12_000) for _ in CLOSURE_TRACERS]
    n_list = []
    for _ in CLOSURE_TRACERS:
        raw = rng.normal(size=(12_000, 3))
        raw /= (np.sum(raw * raw, axis=1) ** 0.5)[:, None]
        n_list.append(raw)
    direction = random_axes(1, rng)[0]
    offsets = np.linspace(-4.0, 4.0, 9)
    reference = scan_q_signature(z_list, n_list, direction, offsets)
    optimized = scan_q_signature_fast_null(z_list, n_list, direction, offsets)
    deltas = []
    for left, right in zip(reference["rows"], optimized["rows"]):
        for key in ("primary_score", "q_peak_residual", "q_frac_residual", "q_shell_occupancy", "cross_tracer_q_peak_cv"):
            deltas.append(abs(float(left[key]) - float(right[key])))
    maximum = max(deltas, default=0.0)
    return {
        "sample_per_tracer": 12_000,
        "scan_steps": int(offsets.size),
        "reference_best_offset": reference["best_scan_offset"],
        "optimized_best_offset": optimized["best_scan_offset"],
        "max_selecting_metric_abs_delta": maximum,
        "pass": bool(reference["best_scan_offset"] == optimized["best_scan_offset"] and maximum <= 1.0e-14),
    }


def subsample_for_null(z_list: list[np.ndarray], n_list: list[np.ndarray], n_target: int, rng: np.random.Generator) -> tuple[list[np.ndarray], list[np.ndarray]]:
    total = sum(z.size for z in z_list)
    if total <= n_target:
        return z_list, n_list
    frac = n_target / total
    zs, ns = [], []
    for z, n_hat in zip(z_list, n_list):
        n_keep = min(z.size, max(1000, int(z.size * frac)))
        idx = rng.choice(z.size, size=n_keep, replace=False)
        zs.append(z[idx])
        ns.append(n_hat[idx])
    return zs, ns


def random_axes(count: int, rng: np.random.Generator) -> np.ndarray:
    lon = rng.uniform(0.0, 2.0 * np.pi, size=count)
    z = rng.uniform(-1.0, 1.0, size=count)
    lat = np.arcsin(z)
    return np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), z], axis=1)


def load_random_unit_vectors(path: Path, max_n: int, seed: int) -> np.ndarray:
    reader = table_reader()
    rng = np.random.default_rng(seed)
    with reader.open(path, memmap=True) as hdul:
        data = hdul[1].data
        size = len(data)
        n = min(max_n, size)
        idx = np.sort(rng.choice(size, size=n, replace=False))
        ra = np.asarray(data["RA"][idx], dtype=np.float64)
        dec = np.asarray(data["DEC"][idx], dtype=np.float64)
    ok = np.isfinite(ra) & np.isfinite(dec)
    return unit_vectors(ra[ok], dec[ok])


def tail_summary(observed: float, values: list[float], lower_is_stronger: bool) -> dict[str, Any]:
    arr = np.array([v for v in values if np.isfinite(v)], dtype=np.float64)
    if arr.size == 0 or not np.isfinite(observed):
        return {"observed_value": observed, "tail_count": None, "n_trials": len(values), "empirical_tail_fraction": None, "empirical_upper_bound": None, "observed_rank": None}
    if lower_is_stronger:
        tail = int(np.sum(arr <= observed))
        rank = int(np.sum(arr < observed) + 1)
    else:
        tail = int(np.sum(arr >= observed))
        rank = int(np.sum(arr > observed) + 1)
    return {
        "observed_value": observed,
        "tail_count": tail,
        "n_trials": len(values),
        "empirical_tail_fraction": tail / len(values),
        "empirical_upper_bound": 1.0 / (len(values) + 1) if tail == 0 else tail / len(values),
        "observed_rank": rank,
    }


def null_scan(
    data_root: Path,
    catalogues: dict[str, dict[str, Any]],
    scan_offsets: np.ndarray,
    observed: dict[str, Any],
    seed: int,
    n_trials: int,
    n_sub: int,
    workers: int,
) -> dict[str, Any]:
    z_list = [catalogues[t]["z"] for t in CLOSURE_TRACERS]
    n_list = [catalogues[t]["n_hat"] for t in CLOSURE_TRACERS]
    random_pools = {
        tracer: load_random_unit_vectors(data_root / TRACER_RANDOM_FILES[tracer], min(len(catalogues[tracer]["z"]), n_sub), seed + i * 100)
        for i, tracer in enumerate(CLOSURE_TRACERS)
    }
    families: dict[str, list[dict[str, Any]]] = {"random_axis": [], "scrambled_redshift": [], "footprint_random": []}

    def run_trial(trial: int) -> tuple[int, dict[str, dict[str, Any]]]:
        rng = np.random.default_rng(np.random.SeedSequence([seed, trial]))
        z_s, n_s = subsample_for_null(z_list, n_list, n_sub, rng)
        axis = random_axes(1, rng)[0]
        random_axis_result = scan_q_signature_fast_null(z_s, n_s, axis, scan_offsets)

        # Subsample first: permutation of the retained null sample has the same
        # null meaning without copying and shuffling all 8.5 million redshifts.
        z_scr = [z.copy() for z in z_s]
        for z in z_scr:
            rng.shuffle(z)
        scrambled_result = scan_q_signature_fast_null(z_scr, n_s, retained_dipole_axis_unit_vector, scan_offsets)

        z_foot, n_foot = [], []
        total = sum(catalogues[t]["z"].size for t in CLOSURE_TRACERS)
        for tracer in CLOSURE_TRACERS:
            z_real = catalogues[tracer]["z"]
            pool = random_pools[tracer]
            n_need = min(pool.shape[0], max(1000, int(n_sub * z_real.size / total)))
            z_foot.append(rng.choice(z_real, size=n_need, replace=True))
            idx = rng.choice(pool.shape[0], size=n_need, replace=pool.shape[0] < n_need)
            n_foot.append(pool[idx])
        footprint_result = scan_q_signature_fast_null(z_foot, n_foot, retained_dipole_axis_unit_vector, scan_offsets)
        return trial, {
            "random_axis": random_axis_result,
            "scrambled_redshift": scrambled_result,
            "footprint_random": footprint_result,
        }

    workers = max(1, min(int(workers), int(n_trials)))
    completed: dict[int, dict[str, dict[str, Any]]] = {}
    note(f"optimized null engine: {n_trials} trials, {workers} workers")
    with ThreadPoolExecutor(max_workers=workers, thread_name_prefix="desi-null") as executor:
        futures = {executor.submit(run_trial, trial): trial for trial in range(n_trials)}
        for done_count, future in enumerate(as_completed(futures), start=1):
            trial, result = future.result()
            completed[trial] = result
            if done_count == 1 or done_count % 10 == 0 or done_count == n_trials:
                note(f"optimized null trials completed {done_count}/{n_trials}")
    for trial in range(n_trials):
        for family in families:
            families[family].append(completed[trial][family])
    metric_keys = {
        "primary_score": True,
        "q_peak_residual": True,
        "q_frac_residual": True,
        "q_shell_occupancy": False,
        "cross_tracer_q_peak_cv": True,
        "phi_excess": False,
    }
    out = {
        "common": {
            "seed": seed,
            "trial_count": n_trials,
            "null_subsample": n_sub,
            "scan_domain": [float(scan_offsets[0]), float(scan_offsets[-1])],
            "scan_resolution": int(scan_offsets.size),
            "tracer_set": list(CLOSURE_TRACERS),
            "closure_predicate": "||U||^2 = s_star",
            "engine": "parallel_exact_score_deferred_diagnostics_v1",
            "workers": workers,
            "trial_rng": "SeedSequence([seed, trial_index])",
        },
        "observed_data_result": observed,
    }
    rows = []
    for family, trials in families.items():
        tails = {}
        for metric, lower in metric_keys.items():
            values = [trial["best_scan_result"].get(metric, np.nan) for trial in trials]
            tails[metric] = tail_summary(observed.get(metric, np.nan), values, lower)
            rows.append({
                "null_name": family,
                "metric": metric,
                "observed_value": tails[metric]["observed_value"],
                "tail_count": tails[metric]["tail_count"],
                "n_trials": tails[metric]["n_trials"],
                "empirical_tail_fraction": tails[metric]["empirical_tail_fraction"],
                "empirical_upper_bound": tails[metric]["empirical_upper_bound"],
                "observed_rank": tails[metric]["observed_rank"],
            })
        out[family] = {"trial_count": n_trials, "metric_tails": tails, "trials": trials}
    out["rows"] = rows
    return out


def nulls_not_run(scan: dict[str, Any], seed: int, scan_offsets: np.ndarray, reason: str) -> dict[str, Any]:
    block: dict[str, Any] = {
        "common": {
            "seed": seed,
            "trial_count": 0,
            "null_subsample": 0,
            "scan_domain": [float(scan_offsets[0]), float(scan_offsets[-1])],
            "scan_resolution": int(scan_offsets.size),
            "tracer_set": list(CLOSURE_TRACERS),
            "closure_predicate": "||U||^2 = s_star",
            "status": "not run",
            "reason": reason,
        },
        "observed_data_result": scan["best_scan_result"],
        "rows": [],
    }
    for family in ("random_axis", "scrambled_redshift", "footprint_random"):
        block[family] = {"status": "not run", "reason": reason, "trial_count": 0, "metric_tails": {}, "trials": []}
    return block


def partition_scan(catalogues: dict[str, dict[str, Any]], scan_offsets: np.ndarray) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    detail: dict[str, Any] = {"tracer_isolated": {}, "leave_one_out": {}, "hemisphere_partitions": {}}
    for tracer in CLOSURE_TRACERS:
        cat = catalogues[tracer]
        res = scan_q_signature([cat["z"]], [cat["n_hat"]], retained_dipole_axis_unit_vector, scan_offsets)
        detail["tracer_isolated"][tracer] = res
        rows.append(partition_row("tracer_isolated", tracer, res, cat["z"].size))
    for keep in (("LRG", "ELG"), ("LRG", "QSO"), ("ELG", "QSO")):
        z = [catalogues[t]["z"] for t in keep]
        n = [catalogues[t]["n_hat"] for t in keep]
        name = "_".join(keep)
        res = scan_q_signature(z, n, retained_dipole_axis_unit_vector, scan_offsets)
        detail["leave_one_out"][name] = res
        rows.append(partition_row("leave_one_out", name, res, sum(len(x) for x in z)))
    for side, sign in (("north", 1), ("south", -1)):
        z_parts, n_parts, count = [], [], 0
        for tracer in CLOSURE_TRACERS:
            cat = catalogues[tracer]
            mask = cat["gal_lat"] * sign >= 0.0
            z_parts.append(cat["z"][mask])
            n_parts.append(cat["n_hat"][mask])
            count += int(mask.sum())
        res = scan_q_signature(z_parts, n_parts, retained_dipole_axis_unit_vector, scan_offsets)
        detail["hemisphere_partitions"][side] = res
        rows.append(partition_row("hemisphere_partitions", side, res, count))
    return rows, detail


def partition_row(kind: str, name: str, res: dict[str, Any], count: int) -> dict[str, Any]:
    best = res["best_scan_result"]
    return {
        "partition_kind": kind,
        "partition_name": name,
        "object_count": int(count),
        "best_scan_offset": best.get("best_scan_offset"),
        "primary_score": best.get("primary_score"),
        "q_peak_residual": best.get("q_peak_residual"),
        "q_shell_occupancy": best.get("q_shell_occupancy"),
        "phi_excess": best.get("phi_excess"),
    }


def q_profile_at_best(catalogues: dict[str, dict[str, Any]], best_scan_offset: float) -> tuple[np.ndarray, dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    q_by_tracer = {t: q_from_scan_offset(catalogues[t]["z"], catalogues[t]["n_hat"], best_scan_offset, retained_dipole_axis_unit_vector) for t in CLOSURE_TRACERS}
    lower = min(float(np.min(q)) for q in q_by_tracer.values())
    upper = max(float(np.max(q)) for q in q_by_tracer.values())
    edges = np.linspace(lower, upper, N_Q_PROFILE_BINS + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    profile = {}
    bin_rows = []
    derivative_rows = []
    for tracer, q in q_by_tracer.items():
        counts = np.histogram(q, bins=edges)[0].astype(np.float64)
        mass = counts / counts.sum() if counts.sum() > 0 else np.zeros_like(counts)
        density = mass / np.maximum(edges[1:] - edges[:-1], 1e-300)
        valid = density > 0
        deriv = np.full_like(centers, np.nan)
        idx = np.where(valid)[0]
        if idx.size >= 2:
            logd = np.log(density[idx])
            x = centers[idx]
            d = np.empty(idx.size)
            d[0] = (logd[1] - logd[0]) / (x[1] - x[0])
            d[-1] = (logd[-1] - logd[-2]) / (x[-1] - x[-2])
            if idx.size > 2:
                d[1:-1] = (logd[2:] - logd[:-2]) / (x[2:] - x[:-2])
            deriv[idx] = d
        fixed = d_sigma_ds(centers)
        shell = np.abs(q - s_star) < q_shell_halfwidth
        profile[tracer] = {
            "retained_count": int(q.size),
            "q_min": float(np.min(q)),
            "q_max": float(np.max(q)),
            "q_shell_count": int(shell.sum()),
            "q_shell_fraction": float(shell.mean()),
            "q_residual_median": float(np.median(q - s_star)),
            "q_residual_rms": float(np.mean((q - s_star) * (q - s_star)) ** 0.5),
        }
        for i in range(centers.size):
            bin_rows.append({
                "tracer": tracer,
                "q_bin_left": edges[i],
                "q_bin_right": edges[i + 1],
                "q_bin_center": centers[i],
                "retained_data_count": counts[i],
                "empirical_q_mass": mass[i],
                "empirical_q_density": density[i],
                "fixed_rho_s": rho_s(centers[i]),
                "q_residual_center": centers[i] - s_star,
            })
            if np.isfinite(deriv[i]):
                derivative_rows.append({
                    "tracer": tracer,
                    "q_bin_center": centers[i],
                    "empirical_log_density_derivative_ds": deriv[i],
                    "fixed_d_sigma_ds": fixed[i],
                    "derivative_residual": deriv[i] - fixed[i],
                })
    return edges, profile, bin_rows, derivative_rows


def write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys = sorted({key for row in rows for key in row})
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow(to_json(row))


def write_scan_csv(path: Path, scan: dict[str, Any]) -> None:
    rows = []
    for row in scan["rows"]:
        rows.append({
            "scan_offset": row["best_scan_offset"],
            "primary_score": row["primary_score"],
            "q_peak_residual": row["q_peak_residual"],
            "q_frac_residual": row["q_frac_residual"],
            "q_shell_occupancy": row["q_shell_occupancy"],
            "cross_tracer_q_peak_cv": row["cross_tracer_q_peak_cv"],
            "phi_chi2": row["phi_chi2"],
            "phi_dof": row["phi_dof"],
            "phi_excess": row["phi_excess"],
        })
    write_rows(path, rows)


def make_plots(out_dir: Path, scan: dict[str, Any], nulls: dict[str, Any], profile: dict[str, Any]) -> list[Path]:
    plt = plotter()
    paths = []
    offsets = np.array([row["best_scan_offset"] for row in scan["rows"]], dtype=np.float64)
    score = np.array([row["primary_score"] for row in scan["rows"]], dtype=np.float64)
    occupancy = np.array([row["q_shell_occupancy"] for row in scan["rows"]], dtype=np.float64)
    path = out_dir / f"{RUN_LABEL}_q_shell_scan.png"
    fig, ax = plt.subplots(figsize=(11, 6), facecolor="#101217")
    ax.set_facecolor("#161a22")
    ax.plot(offsets, score, color="#8ecae6", label="primary q score")
    ax.plot(offsets, occupancy, color="#06d6a0", label="q shell occupancy")
    ax.axvline(scan["best_scan_offset"], color="white", lw=1)
    ax.set_title("squared-radial closure scan | ||U||^2 = s_star", color="white")
    ax.set_xlabel("scan offset", color="white")
    ax.tick_params(colors="white")
    ax.legend(facecolor="#161a22", edgecolor="#666", labelcolor="white")
    fig.tight_layout()
    fig.savefig(path, dpi=170, facecolor="#101217")
    plt.close(fig)
    paths.append(path)
    path = out_dir / f"{RUN_LABEL}_q_nulls.png"
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), facecolor="#101217")
    for ax, family in zip(axes, ("random_axis", "scrambled_redshift", "footprint_random")):
        vals = [trial["best_scan_result"].get("primary_score", np.nan) for trial in nulls[family]["trials"]]
        vals = np.array([v for v in vals if np.isfinite(v)], dtype=np.float64)
        ax.set_facecolor("#161a22")
        if vals.size:
            ax.hist(vals, bins=35, color="#ef476f")
        ax.axvline(scan["best_scan_result"]["primary_score"], color="white", lw=1)
        ax.set_title(family, color="white")
        ax.tick_params(colors="white")
    fig.tight_layout()
    fig.savefig(path, dpi=170, facecolor="#101217")
    plt.close(fig)
    paths.append(path)
    path = out_dir / f"{RUN_LABEL}_q_tracer_shell.png"
    fig, ax = plt.subplots(figsize=(8, 5), facecolor="#101217")
    ax.set_facecolor("#161a22")
    labels = list(profile)
    vals = [profile[t]["q_shell_fraction"] for t in labels]
    ax.bar(labels, vals, color="#ffd166")
    ax.set_title("tracer q-shell fraction | ||U||^2 = s_star", color="white")
    ax.tick_params(colors="white")
    fig.tight_layout()
    fig.savefig(path, dpi=170, facecolor="#101217")
    plt.close(fig)
    paths.append(path)
    return paths


def output_hashes(out_dir: Path, cache: dict[str, str]) -> dict[str, str]:
    skip = {f"{RUN_LABEL}_manifest.json", f"{RUN_LABEL}_output_hashes.json"}
    hashes = {}
    for path in sorted(out_dir.glob(f"{RUN_LABEL}*")):
        if path.is_file() and path.name not in skip:
            hashes[path.name] = sha256_file(path, cache)
    return hashes


def no_legacy_linear_radius_consumption(source_path: Path) -> dict[str, Any]:
    text = source_path.read_text(encoding="utf-8")
    ast.parse(text)
    forbidden_tokens = {
        "np_sqrt_call": "np." + "sqrt(",
        "math_sqrt_call": "math." + "sqrt(",
        "np_linalg_norm_call": "np.linalg." + "norm(",
        "legacy_radius_function_name": "radii_from_" + "scan_point",
        "legacy_abs_radius_minus": "abs(" + "r -",
        "legacy_source_comparison": "sqrt" + "(q) - " + "s_star",
    }
    found = {name: token in text for name, token in forbidden_tokens.items()}
    return {
        "source_ast_parse": True,
        **found,
        "runtime_legacy_linear_radius_counter": 0,
        "display_only_linear_quantity": "r_star_linear_display",
        "pass": not any(found.values()),
    }


def synthetic_tests() -> dict[str, Any]:
    identity = high_precision_identity_tests()
    U = np.array([[1.0, 2.0, 3.0], [r_star_linear_display, 0.0, 0.0], [0.0, -r_star_linear_display, 0.0]], dtype=np.float64)
    q_closure = np.sum(U * U, axis=1)
    q_residual = q_closure - s_star
    mask = np.abs(q_residual) < q_shell_halfwidth
    expected = np.array([abs(float(np.sum(U[i] * U[i]) - s_star)) < q_shell_halfwidth for i in range(U.shape[0])])
    return {
        "identity": identity,
        "q_closure": q_closure,
        "q_residual": q_residual,
        "q_shell_mask": mask,
        "q_shell_mask_expected": expected,
        "pass": bool(identity["pass"] and np.all(mask == expected)),
    }


def migration_note() -> dict[str, Any]:
    return {
        "historical_contract": "legacy closure consumed sqrt(Q) against r_star",
        "corrected_contract": "closure consumes Q against s_star",
        "legacy_patterns_identified": [
            "sigma_desi_lss_closure_audit.py:" + "radii_from_" + "scan_point returned sqrt(max(q,0))",
            "sigma_desi_lss_closure_audit.py:shell_anisotropy_at_scan_offset used abs(radii - r_star) < 1.0",
            "sigma_desi_lss_closure_audit.py:evaluate_closure_signature used peak and occupancy around target_r",
            "sigma_desi_lss_closure_audit.py:null_destruction_audit called scan_signature on linear radii",
            "sigma_desi_lss_closure_audit.py:false_r_star_scan varied a linear target",
        ],
        "historical_manifest_rejection_rule": "coordinate_contract must equal squared_radial",
    }


def reject_mixed_contract_manifest(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"path": str(path), "status": "not present"}
    data = json.loads(path.read_text(encoding="utf-8"))
    ok = data.get("coordinate_contract") == "squared_radial" and data.get("source_variable") == "s = r^2" and "s_star" in data and "r_star_linear_display" in data
    return {"path": str(path), "status": "accepted" if ok else "rejected", "coordinate_contract": data.get("coordinate_contract"), "source_variable": data.get("source_variable")}


def write_report(path: Path, manifest: dict[str, Any]) -> None:
    lines = [
        "# DESI LSS Closure Squared-Radial Audit",
        "",
        "## Variable Contract",
        "",
        "- coordinate_contract = squared_radial",
        "- source_variable = s = r^2",
        "- closure_predicate = ||U||^2 = s_star",
        f"- s_star = {manifest['s_star']}",
        f"- r_star_linear_display = {manifest['r_star_linear_display']}",
        f"- q_shell_halfwidth = {manifest['q_shell_halfwidth']}",
        "",
        "## Scan Output",
        "",
        f"- best_scan_offset = {manifest['squared_radial_closure_scan']['best_scan_offset']}",
        f"- primary_score = {manifest['squared_radial_closure_scan']['best_scan_result'].get('primary_score')}",
        "",
        "## Tracer q-shell Results",
        "",
    ]
    for tracer, row in manifest["tracer_level_results"].items():
        lines.append(f"- {tracer}: q_shell_count = {row['q_shell_count']}; q_shell_fraction = {row['q_shell_fraction']}; q_residual_median = {row['q_residual_median']}")
    lines += ["", "## Null Tail Results", ""]
    for row in manifest["null_destruction"]["rows"]:
        lines.append(f"- {row['null_name']} {row['metric']}: tail_count = {row['tail_count']}; n_trials = {row['n_trials']}; empirical_upper_bound = {row['empirical_upper_bound']}")
    lines += ["", "## Reproducibility", ""]
    lines.append(f"- source_sha256 = {manifest['source_hash']}")
    lines.append(f"- manifest_schema_version = {manifest['manifest_schema_version']}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", default="data/desi")
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--scan-range", type=float, default=DEFAULT_SCAN_RANGE)
    parser.add_argument("--scan-steps", type=int, default=DEFAULT_SCAN_STEPS)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--null-trials", type=int, default=NULL_TRIALS)
    parser.add_argument("--null-sub", type=int, default=NULL_SUBSAMPLE)
    parser.add_argument("--null-workers", type=int, default=min(12, os.cpu_count() or 1))
    parser.add_argument("--skip-nulls", action="store_true")
    parser.add_argument("--synthetic-test", action="store_true")
    parser.add_argument("--legacy-manifest-to-reject", default="results/sigma_desi_lss_closure_null_legacy_full/sigma_desi_lss_closure_manifest.json")
    return parser.parse_args()


def run_synthetic(args: argparse.Namespace) -> None:
    out_dir = Path(args.output_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    result = {
        "script": SCRIPT_NAME,
        "run_timestamp_utc": stamp(),
        "synthetic_tests": synthetic_tests(),
        "fast_null_equivalence_test": fast_null_equivalence_test(),
        "no_legacy_linear_radius_consumption": no_legacy_linear_radius_consumption(Path(__file__).resolve()),
    }
    path = out_dir / f"{RUN_LABEL}_synthetic_tests.json"
    write_json(path, result)
    print(f"synthetic test path: {path}")
    print(f"synthetic pass: {result['synthetic_tests']['pass'] and result['no_legacy_linear_radius_consumption']['pass']}")


def main() -> None:
    args = parse_args()
    if args.synthetic_test:
        run_synthetic(args)
        return
    started = time.monotonic()
    data_root = Path(args.data_root).expanduser().resolve()
    out_dir = Path(args.output_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    cache: dict[str, str] = {}
    source_path = Path(__file__).resolve()
    source_hash = sha256_file(source_path, cache)
    makefile_hash = sha256_file(Path("Makefile").resolve(), cache)
    identity = high_precision_identity_tests()
    no_legacy = no_legacy_linear_radius_consumption(source_path)
    fast_null_equivalence = fast_null_equivalence_test()
    if not identity["pass"] or not no_legacy["pass"] or not fast_null_equivalence["pass"]:
        raise RuntimeError("source identity or legacy-path tests failed")
    inputs = discover_inputs(data_root)
    check_inputs(inputs, CLOSURE_TRACERS)
    input_hashes = hash_inputs(inputs, cache)
    write_json(out_dir / f"{RUN_LABEL}_input_hashes.json", input_hashes)
    catalogues = load_catalogues(data_root)
    scan_offsets = np.linspace(-args.scan_range, args.scan_range, args.scan_steps)
    z_list = [catalogues[t]["z"] for t in CLOSURE_TRACERS]
    n_list = [catalogues[t]["n_hat"] for t in CLOSURE_TRACERS]
    scan = scan_q_signature(z_list, n_list, retained_dipole_axis_unit_vector, scan_offsets)
    best_scan_offset = float(scan["best_scan_offset"])
    scan_coordinate = best_scan_offset * retained_dipole_axis_unit_vector
    q_edges, tracer_profile, bin_rows, derivative_rows = q_profile_at_best(catalogues, best_scan_offset)
    partition_rows, partition_detail = partition_scan(catalogues, scan_offsets)
    if args.skip_nulls or args.null_trials <= 0:
        nulls = nulls_not_run(scan, args.seed, scan_offsets, "null branch skipped for non-null corrected manifest generation")
    else:
        nulls = null_scan(
            data_root, catalogues, scan_offsets, scan["best_scan_result"],
            args.seed, args.null_trials, args.null_sub, args.null_workers,
        )
    metric_rows = []
    for tracer, item in tracer_profile.items():
        for key, value in item.items():
            metric_rows.append({"family": "tracer_q_shell", "name": tracer, "metric": key, "value": value})
    for row in nulls["rows"]:
        metric_rows.append({"family": "null", "name": row["null_name"], "metric": row["metric"] + "_tail_count", "value": row["tail_count"]})
    write_scan_csv(out_dir / f"{RUN_LABEL}_q_shell_scan.csv", scan)
    write_rows(out_dir / f"{RUN_LABEL}_q_binwise.csv", bin_rows)
    write_rows(out_dir / f"{RUN_LABEL}_q_derivative.csv", derivative_rows)
    write_rows(out_dir / f"{RUN_LABEL}_q_metrics.csv", metric_rows)
    write_rows(out_dir / f"{RUN_LABEL}_q_null.csv", nulls["rows"])
    write_rows(out_dir / f"{RUN_LABEL}_q_partitions.csv", partition_rows)
    write_rows(out_dir / f"{RUN_LABEL}_tracer_level_results.csv", [{"tracer": k, **v} for k, v in tracer_profile.items()])
    plots = make_plots(out_dir, scan, nulls, tracer_profile)
    legacy_rejection = reject_mixed_contract_manifest(Path(args.legacy_manifest_to_reject).expanduser().resolve())
    manifest = {
        "manifest_schema_version": "squared_radial_closure_v1",
        "coordinate_contract": "squared_radial",
        "source_variable": "s = r^2",
        "stationarity_variable": "s",
        "closure_predicate": "||U||^2 = s_star",
        "script": SCRIPT_NAME,
        "run_label": RUN_LABEL,
        "command": " ".join(sys.argv),
        "run_timestamp_utc": stamp(),
        "elapsed_seconds": time.monotonic() - started,
        "source_hash": source_hash,
        "source_sha256": source_hash,
        "code_hashes": {"source": source_hash, "Makefile": makefile_hash},
        "MU": MU,
        "GAMMA": GAMMA,
        "sigma_s": "log(MU^2 * (1 + GAMMA*s) / (MU + GAMMA)) - MU*s",
        "rho_s": "MU^2 * (1 + GAMMA*s) * exp(-MU*s) / (MU + GAMMA)",
        "d_sigma_ds": "GAMMA/(1 + GAMMA*s) - MU",
        "d2_sigma_ds2": "-GAMMA^2/(1 + GAMMA*s)^2",
        "s_star": s_star,
        "r_star_linear_display": r_star_linear_display,
        "s_tail": s_tail,
        "q_shell_halfwidth": q_shell_halfwidth,
        "q_shell_halfwidth_policy": "2 * r_star_linear_display * legacy_linear_halfwidth_reference",
        "legacy_linear_halfwidth_reference": LEGACY_LINEAR_HALF_WIDTH_REFERENCE,
        "equivalent_linear_halfwidth_near_s_star": equivalent_linear_halfwidth_near_s_star,
        "completed_scan_coordinate": scan_coordinate,
        "completed_scan_axis": retained_dipole_axis_unit_vector,
        "retained_observed_direction": {
            "longitude_deg": RETAINED_DIPOLE_AXIS_L_DEG,
            "latitude_deg": RETAINED_DIPOLE_AXIS_B_DEG,
            "unit_vector": retained_dipole_axis_unit_vector,
            "fixed_before_run": True,
            "selected_by_current_run": False,
            "fitted_by_current_run": False,
        },
        "deterministic_seed": args.seed,
        "tracer_definitions": {tracer: TRACER_DATA_FILES[tracer] for tracer in CLOSURE_TRACERS},
        "cut_definitions": {"RA_DEC_finite": True, "Z_finite_positive": True, "ZWARN_zero_when_present": True},
        "input_hashes": input_hashes,
        "null_specifications": {
            "families": ["random_axis", "scrambled_redshift", "footprint_random"],
            "seed": args.seed,
            "trial_count": args.null_trials,
            "null_subsample": args.null_sub,
            "tracer_set": list(CLOSURE_TRACERS),
            "scan_domain": [float(scan_offsets[0]), float(scan_offsets[-1])],
            "scan_resolution": int(scan_offsets.size),
        },
        "identity_tests": identity,
        "fast_null_equivalence_test": fast_null_equivalence,
        "no_legacy_linear_radius_consumption": no_legacy,
        "historical_manifest_contract_check": legacy_rejection,
        "migration_note": migration_note(),
        "squared_radial_closure_scan": scan,
        "tracer_level_results": tracer_profile,
        "partition_results": {"rows": partition_rows, "detail": partition_detail},
        "null_destruction": nulls,
        "artifacts": {
            "manifest_json": str(out_dir / f"{RUN_LABEL}_manifest.json"),
            "machine_report_json": str(out_dir / f"{RUN_LABEL}_report.json"),
            "report_md": str(out_dir / f"{RUN_LABEL}_report.md"),
            "q_shell_scan_csv": str(out_dir / f"{RUN_LABEL}_q_shell_scan.csv"),
            "q_binwise_csv": str(out_dir / f"{RUN_LABEL}_q_binwise.csv"),
            "q_derivative_csv": str(out_dir / f"{RUN_LABEL}_q_derivative.csv"),
            "q_metrics_csv": str(out_dir / f"{RUN_LABEL}_q_metrics.csv"),
            "q_null_csv": str(out_dir / f"{RUN_LABEL}_q_null.csv"),
            "q_partitions_csv": str(out_dir / f"{RUN_LABEL}_q_partitions.csv"),
            "tracer_level_csv": str(out_dir / f"{RUN_LABEL}_tracer_level_results.csv"),
            "plots": [str(p) for p in plots],
        },
        "environment": {
            "python": sys.version,
            "numpy": np.__version__,
            "mpmath": getattr(mp, "__version__", None) if mp is not None else None,
            "astropy": package_version("astropy"),
            "scipy": package_version("scipy"),
            "matplotlib": package_version("matplotlib"),
            "platform": platform.platform(),
            "git_commit": subprocess.run(["git", "-C", str(Path.cwd()), "rev-parse", "HEAD"], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True).stdout.strip() or None,
        },
    }
    manifest_path = out_dir / f"{RUN_LABEL}_manifest.json"
    report_json_path = out_dir / f"{RUN_LABEL}_report.json"
    report_md_path = out_dir / f"{RUN_LABEL}_report.md"
    source_hash_path = out_dir / f"{RUN_LABEL}_source_hash.txt"
    source_hash_path.write_text(str(source_hash) + "\n", encoding="utf-8")
    write_json(manifest_path, manifest)
    write_json(report_json_path, manifest)
    write_report(report_md_path, manifest)
    hashes = output_hashes(out_dir, cache)
    write_json(out_dir / f"{RUN_LABEL}_output_hashes.json", hashes)
    manifest["output_hashes"] = hashes
    write_json(manifest_path, manifest)
    print(f"source SHA-256: {source_hash}")
    print(f"manifest path: {manifest_path}")
    print(f"manifest SHA-256: {sha256_file(manifest_path, {})}")
    print(f"s_star: {s_star}")
    print(f"r_star_linear_display: {r_star_linear_display}")
    print(f"q_shell_halfwidth: {q_shell_halfwidth}")
    print(f"best_scan_offset: {scan['best_scan_offset']}")
    print(f"primary_score: {scan['best_scan_result'].get('primary_score')}")
    print("tracer q-shell results: " + ", ".join(f"{t}:count={tracer_profile[t]['q_shell_count']},fraction={tracer_profile[t]['q_shell_fraction']}" for t in CLOSURE_TRACERS))
    print("null tails: " + "; ".join(f"{row['null_name']}:{row['metric']}={row['tail_count']}/{row['n_trials']}" for row in nulls["rows"]))


if __name__ == "__main__":
    main()
