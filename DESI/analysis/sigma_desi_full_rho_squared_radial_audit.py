#!/usr/bin/env python3
from __future__ import annotations

import argparse
import ast
import csv
from decimal import Decimal, getcontext
import hashlib
import importlib
import importlib.metadata
import json
import math
from pathlib import Path
import platform
import resource
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

SCRIPT_NAME = "sigma_desi_full_rho_squared_radial_audit.py"
RUN_LABEL = "sigma_desi_full_rho_squared_radial"
DEFAULT_OUTPUT_DIR = Path("results/sigma_desi_full_rho_squared_radial")
DEFAULT_COMPLETED_CLOSURE = Path("results/sigma_desi_lss_closure_squared_radial/sigma_desi_lss_closure_squared_radial_manifest.json")
PROFILE_TRACERS = ("BGS", "LRG", "ELG", "QSO")
HIGH_Z_TRACERS = ("LRG", "ELG", "QSO")
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
N_Q_BINS = 192
RANDOM_CHUNK = 1_000_000
CACHE_SCHEMA_VERSION = "desi_full_rho_squared_radial_cache_v2"
DEFAULT_CACHE_DIR = Path("results/cache/desi_unit_vectors")
RETAINED_MASK_CONTRACT = "RA_DEC_finite;Z_finite_positive;ZWARN_zero_when_present"
NULL_TRIALS = 200
NULL_SUBSAMPLE = 500_000
NULL_SCAN_STEPS = 41
NULL_SCAN_RANGE = 14.0
HIGH_PRECISION_DPS = 4096
HIGH_PRECISION_POINTS = 4_096
PLOT_TEXT = (
    "rho(r^2) = exp(sigma(r^2))\n"
    "d sigma / d(r^2) = gamma/(1 + gamma*r^2) - mu\n"
    "r_star_squared = 9.470447610693826"
)
LOG_LEVEL = "normal"


def stamp() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def note(text: str) -> None:
    if LOG_LEVEL == "verbose":
        print(f"[{time.strftime('%H:%M:%S')}] {text}", flush=True)


def hhmmss(seconds: float) -> str:
    total = max(0, int(round(seconds)))
    h = total // 3600
    m = (total % 3600) // 60
    s = total % 60
    return f"{h:02d}:{m:02d}:{s:02d}"


MU_TEXT = "0.082912607552"
GAMMA_TEXT = "0.38603416"
Q_STAR_TEXT = "9.470447610693826"
Q_TAIL_TEXT = "21.531339548905194"
MU = float(MU_TEXT)
GAMMA = float(GAMMA_TEXT)
q_star = float(Q_STAR_TEXT)
q_tail = float(Q_TAIL_TEXT)


def sigma_q(q: np.ndarray | float) -> np.ndarray | float:
    return np.log(MU**2 * (1.0 + GAMMA * q) / (MU + GAMMA)) - MU * q


def rho_q(q: np.ndarray | float) -> np.ndarray | float:
    return MU**2 * (1.0 + GAMMA * q) * np.exp(-MU * q) / (MU + GAMMA)


def F_q(q: np.ndarray | float) -> np.ndarray | float:
    return 1.0 - np.exp(-MU * q) * (MU + GAMMA * (1.0 + MU * q)) / (MU + GAMMA)


def d_sigma_q_dq(q: np.ndarray | float) -> np.ndarray | float:
    return GAMMA / (1.0 + GAMMA * q) - MU


def d2_sigma_q_dq2(q: np.ndarray | float) -> np.ndarray | float:
    return -(GAMMA * GAMMA) / ((1.0 + GAMMA * q) ** 2)


def hp_mpf(text: str) -> Any:
    if mp is not None:
        return mp.mpf(text)
    return Decimal(text)


def hp_exp(x: Any) -> Any:
    if mp is not None:
        return mp.exp(x)
    return x.exp()


def hp_log(x: Any) -> Any:
    if mp is not None:
        return mp.log(x)
    return x.ln()


def hp_text(value: Any, digits: int = 80) -> str:
    if mp is not None:
        return mp.nstr(value, digits)
    return format(value, "f")


def hp_sigma_q(q: Any) -> Any:
    mu = hp_mpf(MU_TEXT)
    gamma = hp_mpf(GAMMA_TEXT)
    one = hp_mpf("1")
    return hp_log(mu**2 * (one + gamma * q) / (mu + gamma)) - mu * q


def hp_rho_q(q: Any) -> Any:
    mu = hp_mpf(MU_TEXT)
    gamma = hp_mpf(GAMMA_TEXT)
    one = hp_mpf("1")
    return mu**2 * (one + gamma * q) * hp_exp(-mu * q) / (mu + gamma)


def hp_F_q(q: Any) -> Any:
    mu = hp_mpf(MU_TEXT)
    gamma = hp_mpf(GAMMA_TEXT)
    one = hp_mpf("1")
    return one - hp_exp(-mu * q) * (mu + gamma * (one + mu * q)) / (mu + gamma)


def hp_d_sigma_q_dq(q: Any) -> Any:
    mu = hp_mpf(MU_TEXT)
    gamma = hp_mpf(GAMMA_TEXT)
    one = hp_mpf("1")
    return gamma / (one + gamma * q) - mu


def hp_d2_sigma_q_dq2(q: Any) -> Any:
    gamma = hp_mpf(GAMMA_TEXT)
    one = hp_mpf("1")
    return -(gamma * gamma) / ((one + gamma * q) ** 2)


def high_precision_checks() -> dict[str, Any]:
    if mp is not None:
        mp.mp.dps = HIGH_PRECISION_DPS
    else:
        getcontext().prec = HIGH_PRECISION_DPS
    mu = hp_mpf(MU_TEXT)
    gamma = hp_mpf(GAMMA_TEXT)
    q_star_hp = hp_mpf(Q_STAR_TEXT)
    q_tail_hp = hp_mpf(Q_TAIL_TEXT)
    one = hp_mpf("1")
    far = hp_mpf("1000") / mu
    limit = q_tail_hp * hp_mpf("12")
    max_delta = hp_mpf("0")
    max_delta_q = hp_mpf("0")
    selected_indices = sorted(set([0, 1, 2, 3, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, HIGH_PRECISION_POINTS - 1]))
    high_precision_points = [limit * hp_mpf(str(i)) / hp_mpf(str(HIGH_PRECISION_POINTS - 1)) for i in selected_indices]
    high_precision_points.extend([q_star_hp, q_tail_hp, one / mu, far])
    for q in high_precision_points:
        delta = abs(hp_rho_q(q) - hp_exp(hp_sigma_q(q)))
        if delta > max_delta:
            max_delta = delta
            max_delta_q = q
    F0 = hp_F_q(hp_mpf("0"))
    Ffar = hp_F_q(far)
    prime = hp_d_sigma_q_dq(q_star_hp)
    second = hp_d2_sigma_q_dq2(q_star_hp)
    second_plus = second + mu * mu
    q_star_formula_delta = (one / mu - one / gamma) - q_star_hp
    q_tail_formula_delta = (q_star_hp + one / mu) - q_tail_hp
    if max_delta > hp_mpf("1e-80"):
        raise RuntimeError("rho_q equality check failed")
    if abs(F0) > hp_mpf("1e-90") or abs(Ffar - one) > hp_mpf("1e-90") or abs(prime) > hp_mpf("1e-16") or abs(second_plus) > hp_mpf("1e-16"):
        raise RuntimeError("squared-radial identity check failed")
    return {
        "precision_engine": "mpmath" if mp is not None else "decimal",
        "precision_digits": HIGH_PRECISION_DPS,
        "audit_points": HIGH_PRECISION_POINTS,
        "evaluated_high_precision_points": len(high_precision_points),
        "q_grid_min": "0",
        "q_grid_max": hp_text(limit),
        "q_star_formula_delta": hp_text(q_star_formula_delta),
        "sigma_q_prime_at_q_star": hp_text(prime),
        "sigma_qq_at_q_star": hp_text(second),
        "sigma_qq_plus_mu_squared": hp_text(second_plus),
        "rho_q_equals_exp_sigma_max_abs": hp_text(max_delta),
        "rho_q_equals_exp_sigma_max_abs_q": hp_text(max_delta_q),
        "F_q_zero": hp_text(F0),
        "F_q_far": hp_text(Ffar),
        "rho_q_zero": hp_text(hp_rho_q(hp_mpf("0"))),
        "rho_q_q_star": hp_text(hp_rho_q(q_star_hp)),
        "F_q_q_star": hp_text(hp_F_q(q_star_hp)),
        "q_tail_formula_delta": hp_text(q_tail_formula_delta),
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
    out = h.hexdigest()
    if cache is not None:
        cache[key] = out
    return out


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sha256_array(arr: np.ndarray) -> str:
    canonical = np.asarray(arr, dtype="<f8")
    h = hashlib.sha256()
    h.update(str(canonical.shape).encode("ascii"))
    h.update(canonical.tobytes(order="C"))
    return h.hexdigest()


def rss_mb() -> float | None:
    try:
        value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    except Exception:
        return None
    scale = 1024.0 * 1024.0 if sys.platform == "darwin" else 1024.0
    return value / scale


def add_elapsed(benchmark: dict[str, Any] | None, key: str, seconds: float) -> None:
    if benchmark is not None:
        benchmark[key] = float(benchmark.get(key, 0.0) + seconds)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


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
        x = float(value)
        return None if not math.isfinite(x) else x
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, Path):
        return str(value)
    return value


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.write_text(json.dumps(to_json(data), indent=2, sort_keys=True) + "\n", encoding="utf-8")


def load_completed_closure(path: Path, cache: dict[str, str], allow_legacy_provenance: bool = False) -> dict[str, Any]:
    if not path.exists():
        raise RuntimeError(f"completed closure missing: {path}")
    manifest = read_json(path)
    digest = sha256_file(path, cache)
    coordinate_contract = manifest.get("coordinate_contract")
    if coordinate_contract == "squared_radial":
        if manifest.get("source_variable") != "s = r^2":
            raise RuntimeError("corrected closure source_variable differs")
        if manifest.get("stationarity_variable") != "s":
            raise RuntimeError("corrected closure stationarity_variable differs")
        if manifest.get("closure_predicate") != "||U||^2 = s_star":
            raise RuntimeError("corrected closure predicate differs")
        if "s_star" not in manifest or "r_star_linear_display" not in manifest:
            raise RuntimeError("corrected closure manifest missing s_star or r_star_linear_display")
        if abs(float(manifest["s_star"]) - q_star) > 1.0e-12:
            raise RuntimeError("corrected closure s_star differs")
        if abs(float(manifest["MU"]) - MU) > 1.0e-15:
            raise RuntimeError("corrected closure MU differs")
        if abs(float(manifest["GAMMA"]) - GAMMA) > 1.0e-12:
            raise RuntimeError("corrected closure GAMMA differs")
        source_sha = manifest.get("source_hash") or manifest.get("source_sha256")
        if not isinstance(source_sha, str) or len(source_sha) != 64:
            raise RuntimeError("corrected closure source SHA-256 missing")
        direction_block = manifest.get("retained_observed_direction", {})
        unit = np.array(direction_block.get("unit_vector", []), dtype=np.float64)
        if unit.shape != (3,) or not np.all(np.isfinite(unit)):
            raise RuntimeError("corrected closure retained observed direction missing")
        unit_norm2 = float(np.sum(unit * unit))
        if abs(unit_norm2 - 1.0) > 1.0e-12:
            raise RuntimeError("corrected closure retained direction unit check failed")
        observed = manifest.get("squared_radial_closure_scan", {}).get("best_scan_result", {})
        best_scan_offset = manifest.get("squared_radial_closure_scan", {}).get("best_scan_offset")
        if best_scan_offset is None:
            raise RuntimeError("corrected closure best_scan_offset missing")
        C = float(best_scan_offset) * unit
        tracer_set = tuple(manifest.get("null_specifications", {}).get("tracer_set", []))
        if tracer_set != HIGH_Z_TRACERS:
            raise RuntimeError(f"corrected closure high-z tracer set differs: {tracer_set}")
        return {
            "path": str(path),
            "sha256": digest,
            "coordinate_contract": "squared_radial",
            "source_variable": "s = r^2",
            "source_sha256": source_sha,
            "direction_longitude_deg": direction_block.get("longitude_deg"),
            "direction_latitude_deg": direction_block.get("latitude_deg"),
            "retained_dipole_axis_unit_vector": unit,
            "retained_dipole_axis_unit_norm2": unit_norm2,
            "best_scan_offset": float(best_scan_offset),
            "scan_coordinate": C,
            "scan_coordinate_norm2": float(np.sum(C * C)),
            "high_z_tracer_set": list(HIGH_Z_TRACERS),
            "completed_metrics": {
                "q_peaks": observed.get("q_peaks"),
                "cross_tracer_q_peak_cv": observed.get("cross_tracer_q_peak_cv"),
                "q_shell_occupancy": observed.get("q_shell_occupancy"),
                "q_peak_residual": observed.get("q_peak_residual"),
                "phi_chi2": observed.get("phi_chi2"),
                "phi_dof": observed.get("phi_dof"),
                "phi_excess": observed.get("phi_excess"),
                "null_trials": manifest.get("null_specifications", {}).get("trial_count"),
                "null_seed": manifest.get("null_specifications", {}).get("seed"),
                "scan_domain": manifest.get("null_specifications", {}).get("scan_domain"),
                "scan_resolution": manifest.get("null_specifications", {}).get("scan_resolution"),
                "q_shell_halfwidth": manifest.get("q_shell_halfwidth"),
            },
        }
    if not allow_legacy_provenance:
        raise RuntimeError(
            "mixed-contract closure manifest rejected: coordinate_contract must be squared_radial; "
            "use --allow-legacy-provenance only for historical provenance reads"
        )
    source_sha = manifest.get("source_sha256")
    if not isinstance(source_sha, str) or len(source_sha) != 64:
        raise RuntimeError("completed closure source SHA-256 missing")
    fixed = manifest.get("fixed_identity", {})
    if abs(float(fixed.get("MU")) - MU) > 1.0e-15:
        raise RuntimeError("completed closure MU differs")
    if abs(float(fixed.get("GAMMA")) - GAMMA) > 1.0e-12:
        raise RuntimeError("completed closure GAMMA differs")
    if abs(float(fixed.get("r_star")) - q_star) > 1.0e-12:
        raise RuntimeError("completed closure scalar differs from q_star")
    direction_block = manifest.get("fixed_audit_provenance", {}).get("retained_observed_direction", {})
    unit = np.array(direction_block.get("unit_vector", []), dtype=np.float64)
    if unit.shape != (3,) or not np.all(np.isfinite(unit)):
        raise RuntimeError("retained observed direction missing")
    unit_norm2 = float(np.sum(unit * unit))
    if abs(unit_norm2 - 1.0) > 1.0e-12:
        raise RuntimeError("retained observed direction unit check failed")
    common = manifest.get("null_destruction", {}).get("common", {})
    if tuple(common.get("tracer_set", [])) != HIGH_Z_TRACERS:
        raise RuntimeError("completed high-z tracer set differs")
    observed = manifest.get("null_destruction", {}).get("observed_data_result", {})
    best_scan_offset = observed.get("best_scan_offset")
    if best_scan_offset is None:
        raise RuntimeError("best_scan_offset missing")
    C = float(best_scan_offset) * unit
    return {
        "path": str(path),
        "sha256": digest,
        "coordinate_contract": "legacy_linear_radius_provenance",
        "source_variable": "legacy sqrt(Q)",
        "source_sha256": source_sha,
        "direction_longitude_deg": direction_block.get("longitude_deg"),
        "direction_latitude_deg": direction_block.get("latitude_deg"),
        "retained_dipole_axis_unit_vector": unit,
        "retained_dipole_axis_unit_norm2": unit_norm2,
        "best_scan_offset": float(best_scan_offset),
        "scan_coordinate": C,
        "scan_coordinate_norm2": float(np.sum(C * C)),
        "high_z_tracer_set": list(HIGH_Z_TRACERS),
        "completed_metrics": {
            "peaks": observed.get("peaks"),
            "cross_cv": observed.get("cross_cv"),
            "occupancy": observed.get("occupancy"),
            "phi_chi2": observed.get("phi_chi2"),
            "phi_dof": observed.get("phi_dof"),
            "phi_excess": observed.get("phi_excess"),
            "null_trials": common.get("trial_count"),
            "null_seed": common.get("seed"),
            "scan_domain": common.get("scan_domain"),
            "scan_resolution": common.get("scan_resolution"),
        },
    }


def discover_inputs(data_root: Path) -> dict[str, dict[str, Any]]:
    out = {}
    for tracer in PROFILE_TRACERS:
        data_path = data_root / TRACER_DATA_FILES[tracer]
        random_path = data_root / TRACER_RANDOM_FILES[tracer]
        out[tracer] = {
            "data_path": data_path,
            "random_path": random_path,
            "data_exists": data_path.exists(),
            "random_exists": random_path.exists(),
        }
    return out


def check_inputs(inputs: dict[str, dict[str, Any]]) -> None:
    missing = []
    for entry in inputs.values():
        if not entry["data_exists"]:
            missing.append(str(entry["data_path"]))
        if not entry["random_exists"]:
            missing.append(str(entry["random_path"]))
    if missing:
        raise RuntimeError("missing input: " + "; ".join(missing))


def hash_inputs(inputs: dict[str, dict[str, Any]], cache: dict[str, str]) -> dict[str, Any]:
    out = {}
    for tracer, entry in inputs.items():
        out[tracer] = {
            "data_path": str(entry["data_path"]),
            "data_sha256": sha256_file(entry["data_path"], cache),
            "random_path": str(entry["random_path"]),
            "random_sha256": sha256_file(entry["random_path"], cache),
        }
    return out


ICRS_TO_GALACTIC = np.array([
    [-0.0548755604, -0.8734370902, -0.4838350155],
    [0.4941094279, -0.4448296300, 0.7469822445],
    [-0.8676661490, -0.1980763734, 0.4559837762],
], dtype=np.float64)
ICRS_TO_GALACTIC_SHA256 = sha256_array(ICRS_TO_GALACTIC)


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


def pkg_ver(name: str) -> str | None:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return None


def unit_vectors(ra_deg: np.ndarray, dec_deg: np.ndarray) -> np.ndarray:
    ra = np.radians(np.asarray(ra_deg, dtype=np.float64))
    dec = np.radians(np.asarray(dec_deg, dtype=np.float64))
    v = np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], axis=-1)
    out = v @ ICRS_TO_GALACTIC.T
    norm2 = out[:, 0] * out[:, 0] + out[:, 1] * out[:, 1] + out[:, 2] * out[:, 2]
    scale = np.where(norm2 > 0.0, norm2 ** -0.5, 0.0)
    return out * scale[:, None]


def cache_contract(tracer: str, kind: str, source_sha256: str | None) -> dict[str, Any]:
    return {
        "schema_version": CACHE_SCHEMA_VERSION,
        "kind": kind,
        "tracer": tracer,
        "source_fits_sha256": source_sha256,
        "ICRS_TO_GALACTIC_sha256": ICRS_TO_GALACTIC_SHA256,
        "retained_mask_contract": RETAINED_MASK_CONTRACT if kind == "data" else "finite_RA_DEC_random_sky",
        "unit_vector_normalization": "post_rotation_l2_normalized",
        "coordinate_contract": "squared_radial",
        "source_variable": "s = r^2",
    }


def cache_contract_sha256(contract: dict[str, Any]) -> str:
    return sha256_text(json.dumps(to_json(contract), sort_keys=True, separators=(",", ":")))


def cache_paths(cache_dir: Path, tracer: str, kind: str, contract_sha: str) -> tuple[Path, Path]:
    if kind == "data":
        return cache_dir / f"{tracer}_{contract_sha[:24]}.npz", cache_dir / f"{tracer}_{contract_sha[:24]}.json"
    return cache_dir / f"{tracer}_random_{contract_sha[:24]}.npy", cache_dir / f"{tracer}_random_{contract_sha[:24]}.json"


def retained_mask_summary(mask: np.ndarray, cols: dict[str, np.ndarray]) -> dict[str, Any]:
    return {
        "contract": RETAINED_MASK_CONTRACT,
        "input_rows": int(mask.size),
        "retained_rows": int(mask.sum()),
        "rejected_rows": int(mask.size - int(mask.sum())),
        "has_ZWARN": bool("ZWARN" in cols),
    }


def q_direct_from_z_nhat(z: np.ndarray, n_hat: np.ndarray, C: np.ndarray) -> np.ndarray:
    X = z[:, None] * n_hat
    U = X - C[None, :]
    return U[:, 0] * U[:, 0] + U[:, 1] * U[:, 1] + U[:, 2] * U[:, 2]


def closure_q_from_z_nhat(z: np.ndarray, n_hat: np.ndarray, C: np.ndarray, mode: str = "fast") -> tuple[np.ndarray, float]:
    z = np.asarray(z, dtype=np.float64)
    n_hat = np.asarray(n_hat, dtype=np.float64)
    C = np.asarray(C, dtype=np.float64)
    if mode not in {"fast", "direct", "compare"}:
        raise ValueError(f"unknown r^2 route: {mode}")
    if mode == "direct":
        q_direct = q_direct_from_z_nhat(z, n_hat, C)
        return q_direct, 0.0
    C_norm2 = float(np.dot(C, C))
    dot = n_hat @ C
    q_fast = z * z + C_norm2 - 2.0 * z * dot
    if mode == "compare":
        q_direct = q_direct_from_z_nhat(z, n_hat, C)
        residual = float(np.max(np.abs(q_fast - q_direct))) if q_fast.size else 0.0
        return q_fast, residual
    return q_fast, 0.0


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
        cols = {
            "RA": np.asarray(data["RA"], dtype=np.float64),
            "DEC": np.asarray(data["DEC"], dtype=np.float64),
            "Z": np.asarray(data["Z"], dtype=np.float64),
        }
        if "ZWARN" in names:
            cols["ZWARN"] = np.asarray(data["ZWARN"], dtype=np.float64)
    return cols


def build_data_unit_vector_cache(tracer: str, path: Path, source_sha256: str | None, cache_file: Path, meta_file: Path, benchmark: dict[str, Any] | None) -> dict[str, Any]:
    t0 = time.perf_counter()
    cols = read_data(path)
    mask = retained_mask(cols)
    ra = cols["RA"][mask]
    dec = cols["DEC"][mask]
    z = cols["Z"][mask]
    n_hat = unit_vectors(ra, dec)
    gal_lat = np.degrees(np.arcsin(np.clip(n_hat[:, 2], -1.0, 1.0)))
    summary = retained_mask_summary(mask, cols)
    contract = cache_contract(tracer, "data", source_sha256)
    contract_sha = cache_contract_sha256(contract)
    np.savez_compressed(
        cache_file,
        ra=ra,
        dec=dec,
        z=z,
        n_hat=n_hat,
        gal_lat=gal_lat,
        retained_mask_summary=np.array(json.dumps(summary, sort_keys=True)),
        source_fits_sha256=np.array(source_sha256 or ""),
        cache_contract_sha256=np.array(contract_sha),
    )
    write_json(meta_file, {
        "cache_kind": "data_unit_vectors",
        "tracer": tracer,
        "source_path": str(path),
        "source_fits_sha256": source_sha256,
        "contract": contract,
        "cache_contract_sha256": contract_sha,
        "retained_mask_summary": summary,
        "cache_file": str(cache_file),
    })
    add_elapsed(benchmark, "cache_build_seconds", time.perf_counter() - t0)
    return {
        "ra": ra,
        "dec": dec,
        "z": z,
        "n_hat": n_hat,
        "gal_lat": gal_lat,
        "retained_mask_summary": summary,
        "cache_contract_sha256": contract_sha,
        "cache_file": str(cache_file),
        "cache_hit": False,
    }


def load_data_unit_vector_cache(tracer: str, path: Path, source_sha256: str | None, cache_dir: Path, rebuild: bool, disable: bool, benchmark: dict[str, Any] | None) -> dict[str, Any]:
    if disable:
        t0 = time.perf_counter()
        cols = read_data(path)
        mask = retained_mask(cols)
        ra = cols["RA"][mask]
        dec = cols["DEC"][mask]
        z = cols["Z"][mask]
        n_hat = unit_vectors(ra, dec)
        gal_lat = np.degrees(np.arcsin(np.clip(n_hat[:, 2], -1.0, 1.0)))
        add_elapsed(benchmark, "cache_disabled_unit_vector_seconds", time.perf_counter() - t0)
        return {
            "ra": ra,
            "dec": dec,
            "z": z,
            "n_hat": n_hat,
            "gal_lat": gal_lat,
            "retained_mask_summary": retained_mask_summary(mask, cols),
            "cache_contract_sha256": None,
            "cache_file": None,
            "cache_hit": False,
            "cache_disabled": True,
        }
    cache_dir.mkdir(parents=True, exist_ok=True)
    contract = cache_contract(tracer, "data", source_sha256)
    contract_sha = cache_contract_sha256(contract)
    cache_file, meta_file = cache_paths(cache_dir, tracer, "data", contract_sha)
    if cache_file.exists() and meta_file.exists() and not rebuild:
        t0 = time.perf_counter()
        with np.load(cache_file, allow_pickle=False) as npz:
            cached_sha = str(npz["cache_contract_sha256"])
            cached_source = str(npz["source_fits_sha256"])
            if cached_sha != contract_sha or cached_source != (source_sha256 or ""):
                return build_data_unit_vector_cache(tracer, path, source_sha256, cache_file, meta_file, benchmark)
            out = {
                "ra": np.asarray(npz["ra"], dtype=np.float64),
                "dec": np.asarray(npz["dec"], dtype=np.float64),
                "z": np.asarray(npz["z"], dtype=np.float64),
                "n_hat": np.asarray(npz["n_hat"], dtype=np.float64),
                "gal_lat": np.asarray(npz["gal_lat"], dtype=np.float64),
                "retained_mask_summary": json.loads(str(npz["retained_mask_summary"])),
                "cache_contract_sha256": cached_sha,
                "cache_file": str(cache_file),
                "cache_hit": True,
            }
        add_elapsed(benchmark, "cache_load_seconds", time.perf_counter() - t0)
        return out
    return build_data_unit_vector_cache(tracer, path, source_sha256, cache_file, meta_file, benchmark)


def random_cache_valid_count(path: Path) -> int:
    reader = table_reader()
    total = 0
    with reader.open(path, memmap=True) as hdul:
        data = hdul[1].data
        n = len(data)
        for start in range(0, n, RANDOM_CHUNK):
            stop = min(start + RANDOM_CHUNK, n)
            ra = np.asarray(data["RA"][start:stop], dtype=np.float64)
            dec = np.asarray(data["DEC"][start:stop], dtype=np.float64)
            total += int((np.isfinite(ra) & np.isfinite(dec)).sum())
    return total


def build_random_unit_vector_cache(tracer: str, path: Path, source_sha256: str | None, cache_file: Path, meta_file: Path, benchmark: dict[str, Any] | None) -> dict[str, Any]:
    t0 = time.perf_counter()
    cache_file.parent.mkdir(parents=True, exist_ok=True)
    count = random_cache_valid_count(path)
    n_hat_cache = np.lib.format.open_memmap(cache_file, mode="w+", dtype=np.float64, shape=(count, 3))
    reader = table_reader()
    pos = 0
    with reader.open(path, memmap=True) as hdul:
        data = hdul[1].data
        n = len(data)
        for start in range(0, n, RANDOM_CHUNK):
            stop = min(start + RANDOM_CHUNK, n)
            ra = np.asarray(data["RA"][start:stop], dtype=np.float64)
            dec = np.asarray(data["DEC"][start:stop], dtype=np.float64)
            ok = np.isfinite(ra) & np.isfinite(dec)
            if not np.any(ok):
                continue
            chunk = unit_vectors(ra[ok], dec[ok])
            n_hat_cache[pos:pos + chunk.shape[0], :] = chunk
            pos += chunk.shape[0]
    n_hat_cache.flush()
    contract = cache_contract(tracer, "random", source_sha256)
    contract_sha = cache_contract_sha256(contract)
    write_json(meta_file, {
        "cache_kind": "random_unit_vectors",
        "tracer": tracer,
        "source_path": str(path),
        "source_fits_sha256": source_sha256,
        "contract": contract,
        "cache_contract_sha256": contract_sha,
        "valid_ra_dec_count": count,
        "cache_file": str(cache_file),
    })
    add_elapsed(benchmark, "cache_build_seconds", time.perf_counter() - t0)
    return {
        "n_hat": np.load(cache_file, mmap_mode="r"),
        "valid_ra_dec_count": count,
        "cache_contract_sha256": contract_sha,
        "cache_file": str(cache_file),
        "cache_hit": False,
    }


def load_random_unit_vector_cache(tracer: str, path: Path, source_sha256: str | None, cache_dir: Path, rebuild: bool, disable: bool, benchmark: dict[str, Any] | None) -> dict[str, Any] | None:
    if disable:
        return None
    cache_dir.mkdir(parents=True, exist_ok=True)
    contract = cache_contract(tracer, "random", source_sha256)
    contract_sha = cache_contract_sha256(contract)
    cache_file, meta_file = cache_paths(cache_dir, tracer, "random", contract_sha)
    if cache_file.exists() and meta_file.exists() and not rebuild:
        t0 = time.perf_counter()
        meta = read_json(meta_file)
        if meta.get("cache_contract_sha256") != contract_sha or meta.get("source_fits_sha256") != source_sha256:
            return build_random_unit_vector_cache(tracer, path, source_sha256, cache_file, meta_file, benchmark)
        out = {
            "n_hat": np.load(cache_file, mmap_mode="r"),
            "valid_ra_dec_count": int(meta.get("valid_ra_dec_count", 0)),
            "cache_contract_sha256": contract_sha,
            "cache_file": str(cache_file),
            "cache_hit": True,
        }
        add_elapsed(benchmark, "cache_load_seconds", time.perf_counter() - t0)
        return out
    return build_random_unit_vector_cache(tracer, path, source_sha256, cache_file, meta_file, benchmark)


def load_retained(data_root: Path, closure: dict[str, Any], input_hashes: dict[str, Any], cache_dir: Path, rebuild_cache: bool, disable_cache: bool, q_route: str, benchmark: dict[str, Any] | None) -> tuple[dict[str, dict[str, Any]], dict[str, Any], dict[str, float]]:
    out = {}
    C = closure["scan_coordinate"]
    cache_report: dict[str, Any] = {"data_unit_vectors": {}, "random_unit_vectors": {}}
    equivalence: dict[str, float] = {}
    for tracer in PROFILE_TRACERS:
        note(f"loading {tracer}")
        path = data_root / TRACER_DATA_FILES[tracer]
        cached = load_data_unit_vector_cache(tracer, path, input_hashes[tracer]["data_sha256"], cache_dir, rebuild_cache, disable_cache, benchmark)
        ra = cached["ra"]
        dec = cached["dec"]
        z = cached["z"]
        n_hat = cached["n_hat"]
        gal_lat = cached["gal_lat"]
        t0 = time.perf_counter()
        q, residual = closure_q_from_z_nhat(z, n_hat, C, mode=q_route)
        if q_route != "compare" and z.size:
            n_check = min(4096, z.size)
            idx = np.linspace(0, z.size - 1, n_check, dtype=np.int64)
            _, residual = closure_q_from_z_nhat(z[idx], n_hat[idx], C, mode="compare")
        add_elapsed(benchmark, "data_q_evaluation_seconds", time.perf_counter() - t0)
        equivalence[tracer] = residual
        out[tracer] = {
            "tracer": tracer,
            "path": str(path),
            "ra": ra,
            "dec": dec,
            "z": z,
            "n_hat": n_hat,
            "gal_lat": gal_lat,
            "q": q,
            "retained_count": int(q.size),
            "direct_quadratic_coordinate_max_abs": residual,
        }
        cache_report["data_unit_vectors"][tracer] = {
            "cache_hit": bool(cached["cache_hit"]),
            "cache_disabled": bool(cached.get("cache_disabled", False)),
            "cache_file": cached["cache_file"],
            "cache_contract_sha256": cached["cache_contract_sha256"],
            "retained_mask_summary": cached["retained_mask_summary"],
            "source_fits_sha256": input_hashes[tracer]["data_sha256"],
        }
        note(f"{tracer}: {q.size:,} retained")
    for tracer in PROFILE_TRACERS:
        random_path = data_root / TRACER_RANDOM_FILES[tracer]
        cached_random = load_random_unit_vector_cache(tracer, random_path, input_hashes[tracer]["random_sha256"], cache_dir, rebuild_cache, disable_cache, benchmark)
        if cached_random is not None:
            cache_report["random_unit_vectors"][tracer] = {
                "cache_hit": bool(cached_random["cache_hit"]),
                "cache_file": cached_random["cache_file"],
                "cache_contract_sha256": cached_random["cache_contract_sha256"],
                "valid_ra_dec_count": cached_random["valid_ra_dec_count"],
                "source_fits_sha256": input_hashes[tracer]["random_sha256"],
            }
            out[tracer]["random_n_hat_cache"] = cached_random["n_hat"]
        else:
            cache_report["random_unit_vectors"][tracer] = {"cache_hit": False, "cache_disabled": True, "source_fits_sha256": input_hashes[tracer]["random_sha256"]}
            out[tracer]["random_n_hat_cache"] = None
    return out, cache_report, equivalence


def random_q_histogram(path: Path, random_n_hat: np.ndarray | None, z_source: np.ndarray, q_edges: np.ndarray, C: np.ndarray, seed: int, lat_side: str | None = None, q_route: str = "fast", benchmark: dict[str, Any] | None = None) -> tuple[np.ndarray, int, float]:
    t0 = time.perf_counter()
    reader = table_reader()
    rng = np.random.default_rng(seed)
    hist = np.zeros(q_edges.size - 1, dtype=np.float64)
    total = 0
    max_direct = 0.0
    z_valid = z_source[np.isfinite(z_source) & (z_source > 0)]
    if z_valid.size == 0:
        return hist, 0, max_direct
    if random_n_hat is not None:
        n = int(random_n_hat.shape[0])
        for start in range(0, n, RANDOM_CHUNK):
            stop = min(start + RANDOM_CHUNK, n)
            n_hat = np.asarray(random_n_hat[start:stop], dtype=np.float64)
            if lat_side == "north":
                n_hat = n_hat[n_hat[:, 2] >= 0.0]
            elif lat_side == "south":
                n_hat = n_hat[n_hat[:, 2] < 0.0]
            if n_hat.size == 0:
                continue
            z_random = z_valid[rng.choice(z_valid.size, size=n_hat.shape[0], replace=True)]
            q, residual = closure_q_from_z_nhat(z_random, n_hat, C, mode=q_route)
            max_direct = max(max_direct, residual)
            hist += np.histogram(q, bins=q_edges)[0].astype(np.float64)
            total += int(n_hat.shape[0])
        add_elapsed(benchmark, "random_q_histogram_seconds", time.perf_counter() - t0)
        return hist, total, max_direct
    with reader.open(path, memmap=True) as hdul:
        data = hdul[1].data
        n = len(data)
        for start in range(0, n, RANDOM_CHUNK):
            stop = min(start + RANDOM_CHUNK, n)
            ra = np.asarray(data["RA"][start:stop], dtype=np.float64)
            dec = np.asarray(data["DEC"][start:stop], dtype=np.float64)
            ok = np.isfinite(ra) & np.isfinite(dec)
            if not np.any(ok):
                continue
            n_hat = unit_vectors(ra[ok], dec[ok])
            if lat_side == "north":
                n_hat = n_hat[n_hat[:, 2] >= 0.0]
            elif lat_side == "south":
                n_hat = n_hat[n_hat[:, 2] < 0.0]
            if n_hat.size == 0:
                continue
            z_random = z_valid[rng.choice(z_valid.size, size=n_hat.shape[0], replace=True)]
            q, residual = closure_q_from_z_nhat(z_random, n_hat, C, mode=q_route)
            max_direct = max(max_direct, residual)
            hist += np.histogram(q, bins=q_edges)[0].astype(np.float64)
            total += int(n_hat.shape[0])
    add_elapsed(benchmark, "random_q_histogram_seconds", time.perf_counter() - t0)
    return hist, total, max_direct


def q_support_mass(a: float, b: float) -> float:
    return float(F_q(b) - F_q(a))


def build_q_grid(q_edges: np.ndarray) -> dict[str, np.ndarray]:
    q_centers = 0.5 * (q_edges[:-1] + q_edges[1:])
    return {
        "q_edges": q_edges,
        "q_centers": q_centers,
        "F_q_edges": np.asarray(F_q(q_edges), dtype=np.float64),
        "rho_q_centers": np.asarray(rho_q(q_centers), dtype=np.float64),
        "d_sigma_q_dq_centers": np.asarray(d_sigma_q_dq(q_centers), dtype=np.float64),
        "d2_sigma_q_dq2_centers": np.asarray(d2_sigma_q_dq2(q_centers), dtype=np.float64),
        "opening_branch": 1.0 + GAMMA * q_centers,
        "outer_exp_readout": np.exp(-MU * q_centers),
    }


def fixed_q_bin_mass(q_grid: dict[str, np.ndarray], a: float, b: float) -> np.ndarray:
    q_edges = q_grid["q_edges"]
    left = np.maximum(q_edges[:-1], a)
    right = np.minimum(q_edges[1:], b)
    mass = np.zeros(q_edges.size - 1, dtype=np.float64)
    ok = right > left
    denom = float(F_q(b) - F_q(a))
    if denom > 0:
        full = ok & (left == q_edges[:-1]) & (right == q_edges[1:])
        partial = ok & ~full
        edge_mass = q_grid["F_q_edges"][1:] - q_grid["F_q_edges"][:-1]
        mass[full] = edge_mass[full] / denom
        if np.any(partial):
            mass[partial] = (F_q(right[partial]) - F_q(left[partial])) / denom
    return mass


def density_from_mass(mass: np.ndarray, q_edges: np.ndarray) -> np.ndarray:
    width = q_edges[1:] - q_edges[:-1]
    return np.divide(mass, width, out=np.zeros_like(mass), where=width > 0)


def median_from_cdf(q_edges: np.ndarray, mass: np.ndarray) -> float:
    c = np.cumsum(np.where(np.isfinite(mass), mass, 0.0))
    if c.size == 0 or c[-1] <= 0:
        return np.nan
    c /= c[-1]
    idx = int(np.searchsorted(c, 0.5, side="left"))
    idx = max(0, min(idx, len(mass) - 1))
    prev = c[idx - 1] if idx > 0 else 0.0
    span = c[idx] - prev
    frac = 0.0 if span <= 0 else (0.5 - prev) / span
    return float(q_edges[idx] + frac * (q_edges[idx + 1] - q_edges[idx]))


def fixed_q_median_on_support(a: float, b: float) -> float:
    target = F_q(a) + 0.5 * q_support_mass(a, b)
    lo, hi = a, b
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if F_q(mid) < target:
            lo = mid
        else:
            hi = mid
    return float(0.5 * (lo + hi))


def q_derivative_table(q_grid: dict[str, np.ndarray], density: np.ndarray) -> dict[str, Any]:
    q_centers = q_grid["q_centers"]
    valid = np.isfinite(density) & (density > 0)
    idx = np.where(valid)[0]
    deriv = np.full_like(q_centers, np.nan, dtype=np.float64)
    if idx.size >= 2:
        logd = np.log(density[idx])
        x = q_centers[idx]
        d = np.empty(idx.size, dtype=np.float64)
        d[0] = (logd[1] - logd[0]) / (x[1] - x[0])
        d[-1] = (logd[-1] - logd[-2]) / (x[-1] - x[-2])
        if idx.size > 2:
            d[1:-1] = (logd[2:] - logd[:-2]) / (x[2:] - x[:-2])
        deriv[idx] = d
    fixed = q_grid["d_sigma_q_dq_centers"]
    residual = deriv - fixed
    inner = valid & (q_centers < q_star)
    outer = valid & (q_centers > q_star)
    agree_inner = np.sign(deriv[inner]) == np.sign(fixed[inner])
    agree_outer = np.sign(deriv[outer]) == np.sign(fixed[outer])
    near = int(np.argmin(np.abs(q_centers - q_star))) if q_centers.size else -1
    finite_res = residual[np.isfinite(residual)]
    q_star_derivative_residual = float(residual[near]) if near >= 0 and np.isfinite(residual[near]) else np.nan
    return {
        "empirical_log_density_derivative_dq": deriv,
        "fixed_d_sigma_q_dq": fixed,
        "q_derivative_residual": residual,
        "derivative_sign_empirical": np.sign(deriv),
        "derivative_sign_fixed": np.sign(fixed),
        "inner_derivative_sign_agreement": float(np.mean(agree_inner)) if agree_inner.size else np.nan,
        "q_star_neighbourhood_derivative": float(deriv[near]) if near >= 0 and np.isfinite(deriv[near]) else np.nan,
        "q_star_derivative_residual": q_star_derivative_residual,
        "outer_derivative_sign_agreement": float(np.mean(agree_outer)) if agree_outer.size else np.nan,
        "derivative_L1": float(np.mean(np.abs(finite_res))) if finite_res.size else np.nan,
        "derivative_RMS": float(np.mean(finite_res * finite_res) ** 0.5) if finite_res.size else np.nan,
        "derivative_valid_bin_count": int(np.isfinite(deriv).sum()),
    }


def q_profile_metrics(q_centers: np.ndarray, q_edges: np.ndarray, emp_mass: np.ndarray, emp_density: np.ndarray, fixed_mass: np.ndarray, fixed_density: np.ndarray, q_values: np.ndarray, selection_weight: np.ndarray, random_n: int, a: float, b: float) -> dict[str, Any]:
    c_emp = np.cumsum(emp_mass)
    c_fix = np.cumsum(fixed_mass)
    c_res = c_emp - c_fix
    mass_res = emp_mass - fixed_mass
    dens_res = emp_density - fixed_density
    finite_density = dens_res[np.isfinite(dens_res)]
    sel = selection_weight[np.isfinite(selection_weight) & (selection_weight > 0)]
    outer = q_centers > q_tail
    inner_diff = np.diff(emp_density[q_centers < q_star])
    outer_diff = np.diff(emp_density[q_centers > q_star])
    fixed_peak_idx = int(np.argmax(fixed_mass)) if fixed_mass.size else -1
    empirical_peak_idx = int(np.argmax(emp_mass)) if emp_mass.size else -1
    return {
        "retained_count": int(q_values.size),
        "random_count": int(random_n),
        "observed_q_min": float(a),
        "observed_q_max": float(b),
        "CDF_L1": float(np.mean(np.abs(c_res))),
        "CDF_L2": float(np.mean(c_res * c_res) ** 0.5),
        "KS_max": float(np.max(np.abs(c_res))) if c_res.size else np.nan,
        "q_mass_L1": float(np.sum(np.abs(mass_res))),
        "q_mass_L2": float(np.sum(mass_res * mass_res) ** 0.5),
        "q_density_RMS": float(np.mean(finite_density * finite_density) ** 0.5) if finite_density.size else np.nan,
        "weighted_q_density_RMS": float((np.sum(fixed_mass * dens_res * dens_res) / np.sum(fixed_mass)) ** 0.5) if np.sum(fixed_mass) > 0 else np.nan,
        "peak_q_empirical": float(q_centers[empirical_peak_idx]) if empirical_peak_idx >= 0 else np.nan,
        "peak_q_fixed": float(q_centers[fixed_peak_idx]) if fixed_peak_idx >= 0 else np.nan,
        "peak_q_difference": float(q_centers[empirical_peak_idx] - q_centers[fixed_peak_idx]) if empirical_peak_idx >= 0 and fixed_peak_idx >= 0 else np.nan,
        "empirical_median_q": median_from_cdf(q_edges, emp_mass),
        "fixed_median_q_on_support": fixed_q_median_on_support(a, b),
        "empirical_mean_q": float(np.sum(q_centers * emp_mass)),
        "fixed_mean_q_on_support": float(np.sum(q_centers * fixed_mass)),
        "selection_CV": float(np.std(sel) / np.mean(sel)) if sel.size and np.mean(sel) > 0 else np.nan,
        "monotonicity": {
            "inner_nonnegative_fraction": float(np.mean(inner_diff >= 0)) if inner_diff.size else np.nan,
            "outer_nonpositive_fraction": float(np.mean(outer_diff <= 0)) if outer_diff.size else np.nan,
        },
        "tail_variance": float(np.var(dens_res[outer & np.isfinite(dens_res)])) if np.any(outer & np.isfinite(dens_res)) else np.nan,
    }


def q_branch_readout(q_grid: dict[str, np.ndarray], emp_mass: np.ndarray, fixed_mass: np.ndarray, derivative: dict[str, Any], a: float, b: float) -> dict[str, Any]:
    q_edges = q_grid["q_edges"]
    idx = np.where((q_edges[:-1] <= q_star) & (q_edges[1:] > q_star))[0]
    q_idx = int(idx[0]) if idx.size else -1
    q_star_derivative_residual = derivative["q_star_derivative_residual"]
    return {
        "rho_q_zero": float(rho_q(0.0)),
        "first_observed_q": float(a),
        "rho_q_at_first_observed_q": float(rho_q(a)),
        "first_observed_support_mass": q_support_mass(a, b),
        "opening_at_q_star": float(1.0 + GAMMA * q_star),
        "q_star": float(q_star),
        "rho_q_q_star": float(rho_q(q_star)),
        "F_q_q_star": float(F_q(q_star)),
        "q_star_bin_index": q_idx,
        "q_star_bin_left": float(q_edges[q_idx]) if q_idx >= 0 else np.nan,
        "q_star_bin_right": float(q_edges[q_idx + 1]) if q_idx >= 0 else np.nan,
        "empirical_q_mass_in_q_star_bin": float(emp_mass[q_idx]) if q_idx >= 0 else np.nan,
        "fixed_q_mass_in_q_star_bin": float(fixed_mass[q_idx]) if q_idx >= 0 else np.nan,
        "q_star_bin_residual": float(emp_mass[q_idx] - fixed_mass[q_idx]) if q_idx >= 0 else np.nan,
        "q_star_derivative_residual": q_star_derivative_residual,
        "outer_exp_at_upper_support": float(np.exp(-MU * b)),
        "support_contains_q_star": bool(a <= q_star <= b),
        "support_contains_inner_branch": bool(a <= 0.0 <= b),
        "support_contains_outer_branch": bool(b > q_star),
        "support_contains_q_tail": bool(a <= q_tail <= b),
        "q_star_false_neighbour_separations": {
            "to_completed_linear_r_star_squared": float((q_star * q_star) - q_star),
            "to_completed_linear_r_tail_squared": float((q_tail * q_tail) - q_star),
        },
    }


def readout_from_q_hist(tracer: str, q_values: np.ndarray, z_source: np.ndarray, random_path: Path, random_n_hat: np.ndarray | None, q_grid: dict[str, np.ndarray], closure: dict[str, Any], seed: int, lat_side: str | None = None, q_route: str = "fast", benchmark: dict[str, Any] | None = None) -> dict[str, Any]:
    a = float(np.min(q_values))
    b = float(np.max(q_values))
    q_edges = q_grid["q_edges"]
    q_centers = q_grid["q_centers"]
    data_count = np.histogram(q_values, bins=q_edges)[0].astype(np.float64)
    rr, random_n, random_direct = random_q_histogram(random_path, random_n_hat, z_source, q_edges, closure["scan_coordinate"], seed, lat_side, q_route, benchmark)
    with np.errstate(divide="ignore", invalid="ignore"):
        selection_weight = data_count / rr
    selection_weight[~np.isfinite(selection_weight)] = np.nan
    support = (q_edges[1:] > a) & (q_edges[:-1] < b)
    weight = np.where(support & np.isfinite(selection_weight) & (selection_weight >= 0), selection_weight, 0.0)
    total = float(weight.sum())
    emp_mass = weight / total if total > 0 else np.zeros_like(weight)
    emp_density = density_from_mass(emp_mass, q_edges)
    fix_mass = fixed_q_bin_mass(q_grid, a, b)
    fix_density = density_from_mass(fix_mass, q_edges)
    derivative = q_derivative_table(q_grid, emp_density)
    metrics = q_profile_metrics(q_centers, q_edges, emp_mass, emp_density, fix_mass, fix_density, q_values, selection_weight, random_n, a, b)
    branches = q_branch_readout(q_grid, emp_mass, fix_mass, derivative, a, b)
    return {
        "tracer": tracer,
        "q_centers": q_centers,
        "q_edges": q_edges,
        "retained_data_count": data_count,
        "retained_random_count_by_bin": rr,
        "existing_selection_corrected_weight": selection_weight,
        "empirical_q_mass": emp_mass,
        "empirical_q_density": emp_density,
        "fixed_rho_q_support_density": fix_density,
        "fixed_rho_q_support_bin_mass": fix_mass,
        "empirical_minus_fixed_q_density": emp_density - fix_density,
        "empirical_minus_fixed_q_bin_mass": emp_mass - fix_mass,
        "cumulative_empirical_q_mass": np.cumsum(emp_mass),
        "cumulative_fixed_q_mass": np.cumsum(fix_mass),
        "cumulative_q_residual": np.cumsum(emp_mass) - np.cumsum(fix_mass),
        "support": {
            "observed_q_min": a,
            "observed_q_max": b,
            "support_identity_mass": q_support_mass(a, b),
            "support_contains_q_star": bool(a <= q_star <= b),
            "support_contains_inner_branch": bool(a <= 0.0 <= b),
            "support_contains_outer_branch": bool(b > q_star),
            "support_contains_q_tail": bool(a <= q_tail <= b),
            "outside_wording": "outside retained DESI squared-radial support",
        },
        "metrics": metrics,
        "derivative": derivative,
        "branches": branches,
        "opening_branch": q_grid["opening_branch"],
        "outer_exp_readout": q_grid["outer_exp_readout"],
        "random_direct_quadratic_coordinate_max_abs": random_direct,
    }


def readout_no_random(name: str, q_values: np.ndarray, q_grid: dict[str, np.ndarray]) -> dict[str, Any]:
    a = float(np.min(q_values))
    b = float(np.max(q_values))
    q_edges = q_grid["q_edges"]
    q_centers = q_grid["q_centers"]
    data_count = np.histogram(q_values, bins=q_edges)[0].astype(np.float64)
    support = (q_edges[1:] > a) & (q_edges[:-1] < b)
    weight = np.where(support, data_count, 0.0)
    total = float(weight.sum())
    emp_mass = weight / total if total > 0 else np.zeros_like(weight)
    emp_density = density_from_mass(emp_mass, q_edges)
    fix_mass = fixed_q_bin_mass(q_grid, a, b)
    fix_density = density_from_mass(fix_mass, q_edges)
    derivative = q_derivative_table(q_grid, emp_density)
    metrics = q_profile_metrics(q_centers, q_edges, emp_mass, emp_density, fix_mass, fix_density, q_values, np.where(data_count > 0, data_count, np.nan), 0, a, b)
    branches = q_branch_readout(q_grid, emp_mass, fix_mass, derivative, a, b)
    return {
        "name": name,
        "empirical_q_mass": emp_mass,
        "empirical_q_density": emp_density,
        "fixed_rho_q_support_bin_mass": fix_mass,
        "fixed_rho_q_support_density": fix_density,
        "cumulative_empirical_q_mass": np.cumsum(emp_mass),
        "cumulative_fixed_q_mass": np.cumsum(fix_mass),
        "metrics": metrics,
        "derivative": derivative,
        "branches": branches,
        "support": {"observed_q_min": a, "observed_q_max": b, "support_identity_mass": q_support_mass(a, b)},
    }


def full_q_readout(data_root: Path, catalogues: dict[str, dict[str, Any]], q_grid: dict[str, np.ndarray], closure: dict[str, Any], seed: int, q_route: str, benchmark: dict[str, Any] | None) -> dict[str, Any]:
    out = {}
    for i, tracer in enumerate(PROFILE_TRACERS):
        note(f"r^2 readout {tracer}")
        out[tracer] = readout_from_q_hist(
            tracer,
            catalogues[tracer]["q"],
            catalogues[tracer]["z"],
            data_root / TRACER_RANDOM_FILES[tracer],
            catalogues[tracer].get("random_n_hat_cache"),
            q_grid,
            closure,
            seed + i,
            q_route=q_route,
            benchmark=benchmark,
        )
    return out


def cross_tracer_readout(readouts: dict[str, dict[str, Any]], q_grid: dict[str, np.ndarray]) -> dict[str, Any]:
    q_edges = q_grid["q_edges"]
    shared_a = max(readouts[t]["support"]["observed_q_min"] for t in PROFILE_TRACERS)
    shared_b = min(readouts[t]["support"]["observed_q_max"] for t in PROFILE_TRACERS)
    overlap = (q_edges[1:] > shared_a) & (q_edges[:-1] < shared_b)
    masses = {}
    cdfs = {}
    for tracer in PROFILE_TRACERS:
        m = np.where(overlap, readouts[tracer]["empirical_q_mass"], 0.0)
        total = float(m.sum())
        masses[tracer] = m / total if total > 0 else m
        cdfs[tracer] = np.cumsum(masses[tracer])
    fixed = fixed_q_bin_mass(q_grid, shared_a, shared_b)
    pair_rows = {}
    for i, a in enumerate(PROFILE_TRACERS):
        for b in PROFILE_TRACERS[i + 1:]:
            pair_rows[f"{a}_vs_{b}_CDF_L1"] = float(np.mean(np.abs(cdfs[a] - cdfs[b])))
    matrix = np.vstack([masses[t] for t in PROFILE_TRACERS])
    with np.errstate(divide="ignore", invalid="ignore"):
        cv = np.std(matrix, axis=0) / np.mean(matrix, axis=0)
    cv[~np.isfinite(cv)] = np.nan
    residual_signs = np.sign(matrix - fixed[None, :])
    sign_agree = np.all(residual_signs == residual_signs[0:1, :], axis=0)
    fixed_l1 = {t: float(np.mean(np.abs(np.cumsum(masses[t]) - np.cumsum(fixed)))) for t in PROFILE_TRACERS}
    return {
        "shared_q_interval": [float(shared_a), float(shared_b)],
        **pair_rows,
        "binwise_sign_agreement": float(np.mean(sign_agree[overlap])) if np.any(overlap) else np.nan,
        "binwise_profile_range_mean": float(np.nanmean(np.max(matrix, axis=0) - np.min(matrix, axis=0))),
        "binwise_profile_CV_mean": float(np.nanmean(cv)),
        "fixed_rho_q_shared_interval_CDF_L1_by_tracer": fixed_l1,
    }


def partition_readouts(catalogues: dict[str, dict[str, Any]], q_grid: dict[str, np.ndarray]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    detail: dict[str, Any] = {"tracer_isolated": {}, "north_south": {}, "redshift_slices": {}, "leave_one_out": {}}
    for tracer in PROFILE_TRACERS:
        cat = catalogues[tracer]
        detail["tracer_isolated"][tracer] = readout_no_random(tracer, cat["q"], q_grid)
        rows.append(partition_row("tracer_isolated", tracer, detail["tracer_isolated"][tracer]))
        for side, mask in (("north", cat["gal_lat"] >= 0.0), ("south", cat["gal_lat"] < 0.0)):
            if int(mask.sum()) > 10:
                key = f"{tracer}_{side}"
                detail["north_south"][key] = readout_no_random(key, cat["q"][mask], q_grid)
                rows.append(partition_row("north_south", key, detail["north_south"][key]))
        z_edges = np.linspace(float(np.min(cat["z"])), float(np.max(cat["z"])), 5)
        for k in range(4):
            if k == 3:
                mask = (cat["z"] >= z_edges[k]) & (cat["z"] <= z_edges[k + 1])
            else:
                mask = (cat["z"] >= z_edges[k]) & (cat["z"] < z_edges[k + 1])
            if int(mask.sum()) > 10:
                key = f"{tracer}_z{k}"
                detail["redshift_slices"][key] = readout_no_random(key, cat["q"][mask], q_grid)
                rows.append(partition_row("redshift_slices", key, detail["redshift_slices"][key], {"z_left": float(z_edges[k]), "z_right": float(z_edges[k + 1])}))
    for keep in (("LRG", "ELG"), ("LRG", "QSO"), ("ELG", "QSO")):
        name = "_".join(keep)
        q = np.concatenate([catalogues[t]["q"] for t in keep])
        detail["leave_one_out"][name] = readout_no_random(name, q, q_grid)
        rows.append(partition_row("leave_one_out", name, detail["leave_one_out"][name]))
    return rows, detail


def partition_row(kind: str, name: str, readout: dict[str, Any], extra: dict[str, Any] | None = None) -> dict[str, Any]:
    metrics = readout["metrics"]
    row = {
        "partition_kind": kind,
        "partition_name": name,
        "retained_count": metrics["retained_count"],
        "observed_q_min": metrics["observed_q_min"],
        "observed_q_max": metrics["observed_q_max"],
        "CDF_L1": metrics["CDF_L1"],
        "q_mass_L1": metrics["q_mass_L1"],
        "q_density_RMS": metrics["q_density_RMS"],
        "derivative_L1": readout["derivative"]["derivative_L1"],
        "inner_derivative_sign_agreement": readout["derivative"]["inner_derivative_sign_agreement"],
        "outer_derivative_sign_agreement": readout["derivative"]["outer_derivative_sign_agreement"],
        "q_star_bin_residual": readout["branches"]["q_star_bin_residual"],
        "q_star_derivative_residual": readout["branches"]["q_star_derivative_residual"],
    }
    if extra:
        row.update(extra)
    return row


def random_axes(n: int, rng: np.random.Generator) -> np.ndarray:
    lon = rng.uniform(0.0, 2.0 * np.pi, size=n)
    z = rng.uniform(-1.0, 1.0, size=n)
    lat = np.arcsin(z)
    return np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), z], axis=1)


def sample_catalogue_for_null(catalogues: dict[str, dict[str, Any]], seed: int) -> dict[str, np.ndarray]:
    rng = np.random.default_rng(seed)
    z_parts = []
    n_parts = []
    total = sum(catalogues[t]["z"].size for t in HIGH_Z_TRACERS)
    frac = min(1.0, NULL_SUBSAMPLE / total)
    for tracer in HIGH_Z_TRACERS:
        cat = catalogues[tracer]
        n_keep = min(cat["z"].size, max(1000, int(cat["z"].size * frac)))
        idx = rng.choice(cat["z"].size, size=n_keep, replace=False)
        z_parts.append(cat["z"][idx])
        n_parts.append(cat["n_hat"][idx])
    return {"z": np.concatenate(z_parts), "n_hat": np.vstack(n_parts)}


def sample_random_nhat(path: Path, n: int, seed: int) -> np.ndarray:
    reader = table_reader()
    rng = np.random.default_rng(seed)
    with reader.open(path, memmap=True) as hdul:
        data = hdul[1].data
        size = len(data)
        idx = np.sort(rng.choice(size, size=n, replace=False))
        ra = np.asarray(data["RA"][idx], dtype=np.float64)
        dec = np.asarray(data["DEC"][idx], dtype=np.float64)
    ok = np.isfinite(ra) & np.isfinite(dec)
    return unit_vectors(ra[ok], dec[ok])


def null_metric_from_q(name: str, q: np.ndarray, q_grid: dict[str, np.ndarray]) -> dict[str, float]:
    readout = readout_no_random(name, q, q_grid)
    m = readout["metrics"]
    d = readout["derivative"]
    b = readout["branches"]
    return {
        "CDF_L1": m["CDF_L1"],
        "q_mass_L1": m["q_mass_L1"],
        "q_density_RMS": m["q_density_RMS"],
        "derivative_L1": d["derivative_L1"],
        "inner_derivative_sign_agreement": d["inner_derivative_sign_agreement"],
        "outer_derivative_sign_agreement": d["outer_derivative_sign_agreement"],
        "q_star_bin_residual": b["q_star_bin_residual"],
        "q_star_derivative_residual": b["q_star_derivative_residual"],
    }


def tail_summary(observed: float, values: list[float], lower_is_stronger: bool) -> dict[str, Any]:
    arr = np.array([v for v in values if np.isfinite(v)], dtype=np.float64)
    if arr.size == 0 or not np.isfinite(observed):
        return {"observed_value": observed, "null_distribution": values, "tail_count": None, "n_trials": len(values), "empirical_upper_bound": None, "observed_rank": None}
    if lower_is_stronger:
        tail = int(np.sum(arr <= observed))
        rank = int(np.sum(arr < observed) + 1)
    else:
        tail = int(np.sum(arr >= observed))
        rank = int(np.sum(arr > observed) + 1)
    return {
        "observed_value": observed,
        "null_distribution": values,
        "tail_count": tail,
        "n_trials": len(values),
        "empirical_upper_bound": float(1.0 / (len(values) + 1)) if tail == 0 else float((tail + 1) / (len(values) + 1)),
        "observed_rank": rank,
    }


def random_nhat_sample_from_cache_or_fits(path: Path, random_n_hat: np.ndarray | None, n: int, seed: int) -> np.ndarray:
    if random_n_hat is not None:
        rng = np.random.default_rng(seed)
        size = int(random_n_hat.shape[0])
        idx = np.sort(rng.choice(size, size=min(n, size), replace=False))
        return np.asarray(random_n_hat[idx], dtype=np.float64)
    return sample_random_nhat(path, n, seed)


def null_extension(catalogues: dict[str, dict[str, Any]], q_grid: dict[str, np.ndarray], closure: dict[str, Any], data_root: Path, seed: int, q_route: str) -> tuple[list[dict[str, Any]], dict[str, Any], float]:
    note("r^2 null extension setup")
    sample = sample_catalogue_for_null(catalogues, seed)
    z = sample["z"]
    n_hat = sample["n_hat"]
    rng = np.random.default_rng(seed)
    C = closure["scan_coordinate"]
    observed_q, observed_direct = closure_q_from_z_nhat(z, n_hat, C, mode=q_route)
    observed = null_metric_from_q("observed_high_z", observed_q, q_grid)
    random_pool = random_nhat_sample_from_cache_or_fits(data_root / TRACER_RANDOM_FILES["LRG"], catalogues["LRG"].get("random_n_hat_cache"), min(NULL_SUBSAMPLE, z.size), seed + 9000)
    families = {"random_axis": [], "scrambled_redshift": [], "footprint_random": []}
    direct_max = observed_direct
    for trial in range(NULL_TRIALS):
        if trial % 25 == 0:
            note(f"r^2 null trial {trial}/{NULL_TRIALS}")
        axis = random_axes(1, rng)[0]
        C_axis = closure["best_scan_offset"] * axis
        q_axis, res_axis = closure_q_from_z_nhat(z, n_hat, C_axis, mode=q_route)
        direct_max = max(direct_max, res_axis)
        families["random_axis"].append(null_metric_from_q("random_axis", q_axis, q_grid))
        z_perm = rng.permutation(z)
        q_scramble, res_scramble = closure_q_from_z_nhat(z_perm, n_hat, C, mode=q_route)
        direct_max = max(direct_max, res_scramble)
        families["scrambled_redshift"].append(null_metric_from_q("scrambled_redshift", q_scramble, q_grid))
        use_n = min(random_pool.shape[0], z.size)
        z_rand = z[rng.choice(z.size, size=use_n, replace=True)]
        q_foot, res_foot = closure_q_from_z_nhat(z_rand, random_pool[:use_n], C, mode=q_route)
        direct_max = max(direct_max, res_foot)
        families["footprint_random"].append(null_metric_from_q("footprint_random", q_foot, q_grid))
    rows = []
    detail: dict[str, Any] = {
        "common": {
            "null_trials": NULL_TRIALS,
            "null_subsample": NULL_SUBSAMPLE,
            "null_s_steps": NULL_SCAN_STEPS,
            "null_s_range": NULL_SCAN_RANGE,
            "seed": seed,
            "tracer_set": list(HIGH_Z_TRACERS),
        },
        "observed": observed,
    }
    lower = {
        "CDF_L1": True,
        "q_mass_L1": True,
        "q_density_RMS": True,
        "derivative_L1": True,
        "inner_derivative_sign_agreement": False,
        "outer_derivative_sign_agreement": False,
        "q_star_bin_residual": True,
        "q_star_derivative_residual": True,
    }
    for name, trial_rows in families.items():
        detail[name] = {"trial_metrics": trial_rows, "metric_tails": {}}
        for metric, low in lower.items():
            values = [row.get(metric, np.nan) for row in trial_rows]
            tail = tail_summary(observed.get(metric, np.nan), values, low)
            detail[name]["metric_tails"][metric] = tail
            rows.append({
                "null_name": name,
                "metric": metric,
                "observed_value": tail["observed_value"],
                "tail_count": tail["tail_count"],
                "n_trials": tail["n_trials"],
                "empirical_upper_bound": tail["empirical_upper_bound"],
                "observed_rank": tail["observed_rank"],
            })
    return rows, detail, direct_max


def write_q_edges(path: Path, q_edges: np.ndarray) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["edge_index", "q_edge"])
        for i, x in enumerate(q_edges):
            w.writerow([i, x])


def write_binwise(path: Path, readouts: dict[str, dict[str, Any]]) -> None:
    keys = [
        "tracer", "q_bin_left", "q_bin_right", "q_bin_center", "retained_data_count", "retained_random_count",
        "existing_selection_corrected_weight", "empirical_q_mass", "empirical_q_density",
        "fixed_rho_q_support_density", "fixed_rho_q_support_bin_mass", "empirical_minus_fixed_q_density",
        "empirical_minus_fixed_q_bin_mass", "cumulative_empirical_q_mass", "cumulative_fixed_q_mass",
        "cumulative_q_residual", "opening_branch", "outer_exp_readout",
    ]
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(keys)
        for tracer, ro in readouts.items():
            q_edges = ro["q_edges"]
            q_centers = ro["q_centers"]
            for i in range(q_centers.size):
                w.writerow([
                    tracer, q_edges[i], q_edges[i + 1], q_centers[i], ro["retained_data_count"][i],
                    ro["retained_random_count_by_bin"][i], ro["existing_selection_corrected_weight"][i],
                    ro["empirical_q_mass"][i], ro["empirical_q_density"][i],
                    ro["fixed_rho_q_support_density"][i], ro["fixed_rho_q_support_bin_mass"][i],
                    ro["empirical_minus_fixed_q_density"][i], ro["empirical_minus_fixed_q_bin_mass"][i],
                    ro["cumulative_empirical_q_mass"][i], ro["cumulative_fixed_q_mass"][i],
                    ro["cumulative_q_residual"][i], ro["opening_branch"][i], ro["outer_exp_readout"][i],
                ])


def write_derivative(path: Path, readouts: dict[str, dict[str, Any]]) -> None:
    keys = ["tracer", "q_bin_center", "empirical_log_density_derivative_dq", "fixed_d_sigma_q_dq", "q_derivative_residual", "derivative_sign_empirical", "derivative_sign_fixed"]
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(keys)
        for tracer, ro in readouts.items():
            d = ro["derivative"]
            for i, x in enumerate(ro["q_centers"]):
                if np.isfinite(d["empirical_log_density_derivative_dq"][i]):
                    w.writerow([tracer, x, d["empirical_log_density_derivative_dq"][i], d["fixed_d_sigma_q_dq"][i], d["q_derivative_residual"][i], d["derivative_sign_empirical"][i], d["derivative_sign_fixed"][i]])


def write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys = sorted({k for row in rows for k in row.keys()})
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(keys)
        for row in rows:
            w.writerow([row.get(k, "") for k in keys])


def metric_rows(readouts: dict[str, dict[str, Any]], cross: dict[str, Any], partition_rows: list[dict[str, Any]], null_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for tracer, ro in readouts.items():
        for k, v in ro["metrics"].items():
            if isinstance(v, dict):
                for kk, vv in v.items():
                    rows.append({"section": "q_profile", "name": tracer, "metric": f"{k}_{kk}", "value": vv})
            else:
                rows.append({"section": "q_profile", "name": tracer, "metric": k, "value": v})
        for k in ("derivative_L1", "derivative_RMS", "inner_derivative_sign_agreement", "outer_derivative_sign_agreement", "q_star_derivative_residual"):
            rows.append({"section": "q_derivative", "name": tracer, "metric": k, "value": ro["derivative"][k]})
        rows.append({"section": "q_branches", "name": tracer, "metric": "q_star_bin_residual", "value": ro["branches"]["q_star_bin_residual"]})
    for k, v in cross.items():
        if isinstance(v, (int, float)):
            rows.append({"section": "q_cross_tracer", "name": "shared", "metric": k, "value": v})
    rows.append({"section": "q_partitions", "name": "all", "metric": "row_count", "value": len(partition_rows)})
    rows.append({"section": "q_nulls", "name": "all", "metric": "row_count", "value": len(null_rows)})
    return rows


def pooled_q_readout(readouts: dict[str, dict[str, Any]], closure: dict[str, Any]) -> dict[str, Any]:
    retained = {t: int(readouts[t]["metrics"]["retained_count"]) for t in PROFILE_TRACERS}
    pooled_count = int(sum(retained.values()))
    support_mass = {t: float(readouts[t]["support"]["support_identity_mass"]) for t in PROFILE_TRACERS}
    weighted = {t: (retained[t] * support_mass[t] / pooled_count if pooled_count else 0.0) for t in PROFILE_TRACERS}
    coverage = float(sum(weighted.values()))
    return {
        "pooled_retained_count": pooled_count,
        "pooled_identity_support_coverage": coverage,
        "pooled_identity_support_uncovered": float(1.0 - coverage),
        "per_tracer_weighted_contributions": weighted,
        "q_star_supported_by_all_tracers": bool(all(readouts[t]["support"]["support_contains_q_star"] for t in PROFILE_TRACERS)),
        "q_tail_supported_by_all_tracers": bool(all(readouts[t]["support"]["support_contains_q_tail"] for t in PROFILE_TRACERS)),
        "shared_fixed_contract": {
            "MU": MU,
            "GAMMA": GAMMA,
            "q_star": q_star,
            "q_tail": q_tail,
            "best_scan_offset": closure["best_scan_offset"],
            "scan_coordinate": closure["scan_coordinate"],
        },
    }


def save_npz(path: Path, q_edges: np.ndarray, readouts: dict[str, dict[str, Any]]) -> None:
    data: dict[str, Any] = {"q_edges": q_edges, "q_bin_centers": 0.5 * (q_edges[:-1] + q_edges[1:])}
    for tracer, ro in readouts.items():
        for key in (
            "retained_data_count", "retained_random_count_by_bin", "existing_selection_corrected_weight",
            "empirical_q_mass", "empirical_q_density", "fixed_rho_q_support_density",
            "fixed_rho_q_support_bin_mass", "cumulative_empirical_q_mass", "cumulative_fixed_q_mass",
            "cumulative_q_residual",
        ):
            data[f"{tracer}_{key}"] = ro[key]
        data[f"{tracer}_q_derivative"] = ro["derivative"]["empirical_log_density_derivative_dq"]
    np.savez_compressed(path, **data)


def plot_common(ax: Any) -> None:
    ax.text(0.01, 0.99, PLOT_TEXT, transform=ax.transAxes, va="top", ha="left", fontsize=8, color="white", bbox={"facecolor": "#111111", "alpha": 0.72, "pad": 4})


def make_plots(out_dir: Path, readouts: dict[str, dict[str, Any]], cross: dict[str, Any], partition_rows: list[dict[str, Any]], null_detail: dict[str, Any]) -> list[Path]:
    plt = plotter()
    paths = []
    colors = {"BGS": "#ffd166", "LRG": "#8ecae6", "ELG": "#06d6a0", "QSO": "#ef476f"}
    p = out_dir / f"{RUN_LABEL}_profiles.png"
    fig, ax = plt.subplots(figsize=(11, 6), facecolor="#101217")
    ax.set_facecolor("#161a22")
    for tracer, ro in readouts.items():
        ax.plot(ro["q_centers"], ro["empirical_q_density"], color=colors[tracer], label=f"{tracer} empirical")
        ax.plot(ro["q_centers"], ro["fixed_rho_q_support_density"], color=colors[tracer], ls="--", alpha=0.7, label=f"{tracer} fixed")
    ax.axvline(q_star, color="white", lw=1.0)
    ax.set_title("Full rho(r^2) Readout", color="white")
    ax.set_xlabel("r^2 retained coordinate", color="white")
    ax.set_ylabel("r^2 density", color="white")
    ax.tick_params(colors="white")
    ax.legend(facecolor="#161a22", edgecolor="#666", labelcolor="white", ncol=2, fontsize=8)
    plot_common(ax)
    fig.tight_layout()
    fig.savefig(p, dpi=170, facecolor="#101217")
    plt.close(fig)
    paths.append(p)
    p = out_dir / f"{RUN_LABEL}_cdf.png"
    fig, ax = plt.subplots(figsize=(11, 6), facecolor="#101217")
    ax.set_facecolor("#161a22")
    for tracer, ro in readouts.items():
        ax.plot(ro["q_centers"], ro["cumulative_q_residual"], color=colors[tracer], label=tracer)
    ax.axhline(0, color="white", lw=0.8)
    ax.axvline(q_star, color="white", lw=1.0)
    ax.set_title("r^2 CDF Residual", color="white")
    ax.set_xlabel("r^2 retained coordinate", color="white")
    ax.set_ylabel("empirical r^2 CDF - fixed r^2 CDF", color="white")
    ax.tick_params(colors="white")
    ax.legend(facecolor="#161a22", edgecolor="#666", labelcolor="white")
    plot_common(ax)
    fig.tight_layout()
    fig.savefig(p, dpi=170, facecolor="#101217")
    plt.close(fig)
    paths.append(p)
    p = out_dir / f"{RUN_LABEL}_derivative_architecture.png"
    fig, ax = plt.subplots(figsize=(11, 6), facecolor="#101217")
    ax.set_facecolor("#161a22")
    for tracer, ro in readouts.items():
        ax.plot(ro["q_centers"], ro["derivative"]["empirical_log_density_derivative_dq"], color=colors[tracer], alpha=0.85, label=tracer)
    ax.plot(next(iter(readouts.values()))["q_centers"], d_sigma_q_dq(next(iter(readouts.values()))["q_centers"]), color="white", ls="--", label="fixed d sigma / d(r^2)")
    ax.axvline(q_star, color="white", lw=1.0)
    ax.set_title("d log density / d(r^2)", color="white")
    ax.set_xlabel("r^2 retained coordinate", color="white")
    ax.set_ylabel("d log density / dq", color="white")
    ax.tick_params(colors="white")
    ax.legend(facecolor="#161a22", edgecolor="#666", labelcolor="white")
    plot_common(ax)
    fig.tight_layout()
    fig.savefig(p, dpi=170, facecolor="#101217")
    plt.close(fig)
    paths.append(p)
    p = out_dir / f"{RUN_LABEL}_branch_readout.png"
    fig, ax = plt.subplots(figsize=(11, 6), facecolor="#101217")
    x = next(iter(readouts.values()))["q_centers"]
    ax.set_facecolor("#161a22")
    ax.plot(x, 1.0 + GAMMA * x, color="#ffd166", label="1 + gamma*r^2")
    ax.plot(x, np.exp(-MU * x), color="#06d6a0", label="exp(-mu*r^2)")
    ax.axvline(q_star, color="white", lw=1.0)
    ax.set_title("r^2 Branch Readout", color="white")
    ax.set_xlabel("r^2 retained coordinate", color="white")
    ax.tick_params(colors="white")
    ax.legend(facecolor="#161a22", edgecolor="#666", labelcolor="white")
    plot_common(ax)
    fig.tight_layout()
    fig.savefig(p, dpi=170, facecolor="#101217")
    plt.close(fig)
    paths.append(p)
    p = out_dir / f"{RUN_LABEL}_qstar_neighbourhood.png"
    fig, ax = plt.subplots(figsize=(11, 6), facecolor="#101217")
    ax.set_facecolor("#161a22")
    lo = q_star - 2
    hi = q_star + 2
    for tracer, ro in readouts.items():
        m = (ro["q_centers"] >= lo) & (ro["q_centers"] <= hi)
        ax.plot(ro["q_centers"][m], ro["empirical_q_mass"][m] - ro["fixed_rho_q_support_bin_mass"][m], color=colors[tracer], label=tracer)
    ax.axvline(q_star, color="white", lw=1.0)
    ax.axhline(0, color="white", lw=0.8)
    ax.set_title("Stationary r_star_squared Closure", color="white")
    ax.set_xlabel("r^2 retained coordinate", color="white")
    ax.set_ylabel("r^2 bin mass residual", color="white")
    ax.tick_params(colors="white")
    ax.legend(facecolor="#161a22", edgecolor="#666", labelcolor="white")
    plot_common(ax)
    fig.tight_layout()
    fig.savefig(p, dpi=170, facecolor="#101217")
    plt.close(fig)
    paths.append(p)
    p = out_dir / f"{RUN_LABEL}_cross_tracer.png"
    fig, ax = plt.subplots(figsize=(10, 5), facecolor="#101217")
    labels = [k for k in cross if k.endswith("CDF_L1")]
    vals = [cross[k] for k in labels]
    ax.set_facecolor("#161a22")
    ax.bar(labels, vals, color="#8ecae6")
    ax.set_title("r^2 Cross-Tracer Readout", color="white")
    ax.tick_params(colors="white", rotation=25)
    plot_common(ax)
    fig.tight_layout()
    fig.savefig(p, dpi=170, facecolor="#101217")
    plt.close(fig)
    paths.append(p)
    p = out_dir / f"{RUN_LABEL}_partitions.png"
    fig, ax = plt.subplots(figsize=(10, 5), facecolor="#101217")
    vals = [row["CDF_L1"] for row in partition_rows]
    ax.set_facecolor("#161a22")
    ax.hist(vals, bins=30, color="#ffd166")
    ax.set_title("r^2 Independent Partitions", color="white")
    ax.tick_params(colors="white")
    plot_common(ax)
    fig.tight_layout()
    fig.savefig(p, dpi=170, facecolor="#101217")
    plt.close(fig)
    paths.append(p)
    p = out_dir / f"{RUN_LABEL}_null_comparison.png"
    fig, ax = plt.subplots(figsize=(10, 5), facecolor="#101217")
    labels = list(k for k in ("random_axis", "scrambled_redshift", "footprint_random") if k in null_detail)
    vals = []
    for k in labels:
        value = null_detail[k].get("metric_tails", {}).get("CDF_L1", {}).get("tail_count")
        vals.append(0 if value is None else value)
    ax.set_facecolor("#161a22")
    ax.bar(labels, vals, color="#ef476f")
    ax.set_title("r^2 Null Extension", color="white")
    ax.tick_params(colors="white", rotation=20)
    plot_common(ax)
    fig.tight_layout()
    fig.savefig(p, dpi=170, facecolor="#101217")
    plt.close(fig)
    paths.append(p)
    return paths


def output_hashes(out_dir: Path, cache: dict[str, str]) -> dict[str, str]:
    skip = {f"{RUN_LABEL}_manifest.json", f"{RUN_LABEL}_output_hashes.json"}
    out = {}
    for path in sorted(out_dir.glob(f"{RUN_LABEL}*")):
        if path.is_file() and path.name not in skip:
            out[path.name] = sha256_file(path, cache)
    return out


def write_report(path: Path, manifest: dict[str, Any]) -> None:
    pooled = manifest["pooled_q_readout"]
    identity = manifest["identity_r_squared"]
    closure = manifest["completed_closure_input"]
    route = manifest["q_route_equivalence"]
    lines = [
        "# DESI Full rho r^2 Audit",
        "",
        "\\mu=0.082912607552,",
        "\\qquad",
        "\\gamma=0.38603416,",
        "",
        f"r_\\ast^2={identity['r_star_squared']},",
        "",
        "r^2=\\|Z\\hat n-C\\|^2,",
        "",
        f"N_{{\\rm retained}}={pooled['pooled_retained_count']:,},",
        "",
        "\\mathrm{pooled\\ retained\\ r^2\\ coverage}",
        "=",
        f"{pooled['pooled_identity_support_coverage']}.",
        "",
        "## Fixed Identity",
        "",
        "\\sigma(r^2)=\\log\\!\\left(\\frac{\\mu^2(1+\\gamma r^2)}{\\mu+\\gamma}\\right)-\\mu r^2",
        "",
        "r_\\ast^2=\\frac1\\mu-\\frac1\\gamma",
        "",
        "r^2=\\|Z\\hat n-C\\|^2",
        "",
        f"- mu = {MU}",
        f"- gamma = {GAMMA}",
        f"- r_star_squared = {q_star}",
        f"- r_tail_squared = {q_tail}",
        "",
        "## Corrected r^2 Geometry",
        "",
        f"- completed closure manifest SHA-256 = {manifest['completed_closure_sha256']}",
        f"- best_scan_offset = {closure['best_scan_offset']}",
        f"- scan_coordinate = {to_json(closure['scan_coordinate'])}",
        "",
        "## Corpus Readout",
        "",
        f"- retained objects = {pooled['pooled_retained_count']}",
        f"- pooled retained r^2 coverage = {pooled['pooled_identity_support_coverage']}",
        f"- pooled retained r^2 uncovered = {pooled['pooled_identity_support_uncovered']}",
        "",
        "## Retained r^2 Coverage",
        "",
    ]
    for tracer, item in manifest["retained_r_squared_support"].items():
        lines.append(f"- {tracer}: {item['retained_r_squared_min']} to {item['retained_r_squared_max']}; support mass = {item['support_identity_mass']}")
    lines += [
        "",
        "## Fast/Direct Route Equivalence",
        "",
        f"- r^2 execution route = {manifest['execution_engineering']['r_squared_execution_route']}",
        f"- max abs = {route['global_fast_vs_direct_max_abs']}",
        f"- pass = {route['pass']}",
        f"- full retained corpus compared = {route['full_retained_corpus_compared']}",
        "",
        "## Source Identity Checks",
        "",
    ]
    label_by_key = {
        "q_star_formula_delta": "r_star_squared_formula_delta",
        "sigma_q_prime_at_q_star": "d_sigma_d_r_squared_at_r_star_squared",
        "sigma_qq_plus_mu_squared": "d2_sigma_d_r_squared2_plus_mu_squared",
        "rho_q_equals_exp_sigma_max_abs": "rho_r_squared_equals_exp_sigma_max_abs",
        "F_q_zero": "F_r_squared_zero",
        "F_q_far": "F_r_squared_far",
        "fast_vs_direct_q_max_abs": "fast_vs_direct_r_squared_max_abs",
        "direct_quadratic_coordinate_max_abs": "direct_coordinate_residual_max_abs",
    }
    for k in label_by_key:
        if k in manifest["squared_radial_closure_checks"]:
            lines.append(f"- {label_by_key[k]} = {manifest['squared_radial_closure_checks'][k]}")
    lines.append(f"- linear_radius_path_pass = {manifest['no_legacy_radius_path']['pass']}")
    lines.append(f"- forbidden_calls = {manifest['no_legacy_radius_path']['ast_forbidden_calls']}")
    lines += [
        "",
        "## Input and Output Custody",
        "",
        f"- source SHA-256 = {manifest['source_sha256']}",
        f"- cache_schema_version = {manifest['execution_engineering']['cache_schema_version']}",
        "",
    ]
    for tracer, item in manifest["input_hashes"].items():
        lines.append(f"- {tracer} data SHA-256 = {item['data_sha256']}")
        lines.append(f"- {tracer} random SHA-256 = {item['random_sha256']}")
    lines += [
        "",
        "## Detailed Per-Tracer Readouts",
        "",
    ]
    for tracer, metrics in manifest["r_squared_profile_metrics"].items():
        branch = manifest["r_squared_branch_metrics"][tracer]
        deriv = manifest["r_squared_derivative_metrics"][tracer]
        lines.append(f"- {tracer}: retained_count = {metrics['retained_count']}; CDF_L1 = {metrics['CDF_L1']}; KS_max = {metrics['KS_max']}; r^2 density RMS = {metrics['r_squared_density_RMS']}")
        lines.append(f"- {tracer}: d log density / d(r^2) L1 = {deriv['derivative_L1']}; RMS = {deriv['derivative_RMS']}; r_star_squared residual = {branch['r_star_squared_bin_residual']}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def no_legacy_radius_path(source_path: Path) -> dict[str, Any]:
    text = source_path.read_text(encoding="utf-8")
    tree = ast.parse(text)
    token_np_sqrt = "np." + "sqrt("
    token_math_sqrt = "math." + "sqrt("
    token_norm = "np.linalg." + "norm("
    token_linear_source_1 = "rho_" + "identity("
    token_linear_source_2 = "rho_" + "cdf("
    token_linear_source_3 = "sigma_" + "prime("
    def dotted(node: ast.AST) -> str:
        if isinstance(node, ast.Name):
            return node.id
        if isinstance(node, ast.Attribute):
            base = dotted(node.value)
            return f"{base}.{node.attr}" if base else node.attr
        return ""

    forbidden_calls: list[dict[str, Any]] = []
    forbidden_names = {"np.sqrt", "math.sqrt", "np.linalg.norm"}
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            name = dotted(node.func)
            if name in forbidden_names:
                forbidden_calls.append({"name": name, "line": int(getattr(node, "lineno", -1))})

    checks = {
        "source_ast_parse": True,
        "contains_np_sqrt_call": any(x["name"] == "np.sqrt" for x in forbidden_calls),
        "contains_math_sqrt_call": any(x["name"] == "math.sqrt" for x in forbidden_calls),
        "contains_np_linalg_norm_call": any(x["name"] == "np.linalg.norm" for x in forbidden_calls),
        "ast_forbidden_calls": forbidden_calls,
        "token_contains_np_sqrt_call": token_np_sqrt in text,
        "token_contains_math_sqrt_call": token_math_sqrt in text,
        "token_contains_np_linalg_norm_call": token_norm in text,
        "contains_linear_r_rho_identity_call": token_linear_source_1 in text,
        "contains_linear_r_rho_cdf_call": token_linear_source_2 in text,
        "contains_linear_r_sigma_prime_call": token_linear_source_3 in text,
        "empirical_q_pipeline_uses_q_edges": "q_edges" in text,
        "empirical_q_pipeline_uses_F_q": "F_q(" in text,
        "empirical_q_pipeline_uses_rho_q": "rho_q(" in text,
    }
    checks["pass"] = not any(checks[k] for k in ("contains_np_sqrt_call", "contains_math_sqrt_call", "contains_np_linalg_norm_call", "contains_linear_r_rho_identity_call", "contains_linear_r_rho_cdf_call", "contains_linear_r_sigma_prime_call"))
    if not checks["pass"]:
        raise RuntimeError("legacy linear-radius source path token present")
    return checks


def synthetic_tests() -> dict[str, Any]:
    checks = high_precision_checks()
    C = np.array([1.5, -2.0, 0.25], dtype=np.float64)
    z = np.array([0.5, 1.25, 2.0, 4.0], dtype=np.float64)
    n_hat = np.array([
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.6, 0.8, 0.0],
    ], dtype=np.float64)
    q, residual = closure_q_from_z_nhat(z, n_hat, C, mode="compare")
    q_edges = np.linspace(float(np.min(q)), float(np.max(q)), 9)
    q_grid = build_q_grid(q_edges)
    ro = readout_no_random("synthetic", q, q_grid)
    return {
        "high_precision": checks,
        "synthetic_q": q,
        "direct_quadratic_coordinate_max_abs": residual,
        "synthetic_CDF_L1": ro["metrics"]["CDF_L1"],
        "synthetic_q_star_bin_index": ro["branches"]["q_star_bin_index"],
        "pass": bool(residual <= 1.0e-12 and np.isfinite(ro["metrics"]["CDF_L1"])),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", default="data/desi")
    parser.add_argument("--completed-closure", default=str(DEFAULT_COMPLETED_CLOSURE))
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument("--mode", choices=["profile", "derivative", "branches", "tracers", "partitions", "nulls", "full"], default="full")
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--synthetic-test", action="store_true")
    parser.add_argument("--allow-legacy-provenance", action="store_true")
    parser.add_argument("--benchmark", action="store_true")
    parser.add_argument("--rebuild-cache", action="store_true")
    parser.add_argument("--disable-cache", action="store_true")
    parser.add_argument("--q-route", choices=["fast", "direct", "compare"], default="fast")
    parser.add_argument("--log-level", choices=["quiet", "normal", "verbose"], default="normal")
    return parser.parse_args()


def active_modes(mode: str) -> set[str]:
    if mode == "full":
        return {"profile", "derivative", "branches", "tracers", "partitions", "nulls"}
    if mode == "tracers":
        return {"profile", "derivative", "branches", "tracers"}
    return {"profile", mode}


def run_synthetic(args: argparse.Namespace) -> None:
    out_dir = Path(args.output_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    result = {
        "script": SCRIPT_NAME,
        "run_timestamp_utc": stamp(),
        "synthetic_tests": synthetic_tests(),
        "no_legacy_radius_path": no_legacy_radius_path(Path(__file__).resolve()),
    }
    path = out_dir / f"{RUN_LABEL}_synthetic_tests.json"
    write_json(path, result)
    print(f"synthetic test path: {path}")
    print(f"synthetic pass: {result['synthetic_tests']['pass'] and result['no_legacy_radius_path']['pass']}")


def print_terminal_summary(args: argparse.Namespace, manifest: dict[str, Any], manifest_path: Path, manifest_sha256: str | None) -> None:
    elapsed = hhmmss(float(manifest["elapsed_seconds"]))
    profile = manifest["q_profile_metrics"]
    pooled = manifest["pooled_q_readout"]
    route = manifest["q_route_equivalence"]
    closure = manifest["completed_closure_input"]
    if args.log_level == "quiet":
        print("PASS")
        print(f"manifest: {manifest_path}")
        print(f"source SHA-256: {manifest['source_sha256']}")
        print(f"elapsed: {elapsed}")
        return
    print("[DESI r^2 AUDIT] corrected source contract loaded")
    print(f"[IDENTITY] mu={MU_TEXT} gamma={GAMMA_TEXT} r*^2={Q_STAR_TEXT}")
    print(f"[CLOSURE] C={closure['best_scan_offset']} a_hat | r^2=||Z n_hat-C||^2")
    for tracer in PROFILE_TRACERS:
        print(f"[INPUT] {tracer} {profile[tracer]['retained_count']:,} retained")
    print(f"[CORPUS] {pooled['pooled_retained_count']:,} retained objects")
    print(f"[ROUTE] fast/direct equivalence {'PASS' if route['pass'] else 'FAIL'} | max abs={route['global_fast_vs_direct_max_abs']}")
    print(f"[PROFILE] pooled retained r^2 coverage={pooled['pooled_identity_support_coverage']}")
    print(f"[CUSTODY] source SHA-256: {manifest['source_sha256']}")
    print(f"[CUSTODY] completed closure SHA-256: {manifest['completed_closure_sha256']}")
    print(f"[CUSTODY] manifest SHA-256: {manifest_sha256}")
    if args.log_level == "verbose":
        print("[BENCHMARK] " + json.dumps(to_json(manifest["benchmark"]), sort_keys=True))
        print("[INPUT HASHES] " + json.dumps(to_json(manifest["input_hashes"]), sort_keys=True))
    print(f"[PASS] completed in {elapsed}")
    print(f"manifest: {manifest_path}")


def main() -> None:
    global LOG_LEVEL
    args = parse_args()
    LOG_LEVEL = args.log_level
    if args.synthetic_test:
        run_synthetic(args)
        return
    started = time.perf_counter()
    data_root = Path(args.data_root).expanduser().resolve()
    closure_path = Path(args.completed_closure).expanduser().resolve()
    out_dir = Path(args.output_dir).expanduser().resolve()
    cache_dir = DEFAULT_CACHE_DIR.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    modes = active_modes(args.mode)
    benchmark: dict[str, Any] = {
        "enabled": bool(args.benchmark),
        "cache_build_seconds": 0.0,
        "cache_load_seconds": 0.0,
        "cache_disabled_unit_vector_seconds": 0.0,
        "data_q_evaluation_seconds": 0.0,
        "random_q_histogram_seconds": 0.0,
        "total_elapsed_seconds": None,
        "peak_rss_mb": None,
        "cache_hit_miss_by_tracer": {},
    }
    cache: dict[str, str] = {}
    source_path = Path(__file__).resolve()
    source_hash = sha256_file(source_path, cache)
    (out_dir / f"{RUN_LABEL}_source_hash.txt").write_text(str(source_hash) + "\n", encoding="utf-8")
    high_precision = high_precision_checks()
    no_legacy = no_legacy_radius_path(source_path)
    closure = load_completed_closure(closure_path, cache, allow_legacy_provenance=args.allow_legacy_provenance)
    inputs = discover_inputs(data_root)
    check_inputs(inputs)
    input_hashes = hash_inputs(inputs, cache)
    write_json(out_dir / f"{RUN_LABEL}_input_hashes.json", input_hashes)
    catalogues, cache_report, q_route_equiv = load_retained(data_root, closure, input_hashes, cache_dir, args.rebuild_cache, args.disable_cache, args.q_route, benchmark)
    lower_edge = min(float(np.min(catalogues[t]["q"])) for t in PROFILE_TRACERS)
    upper_edge = max(float(np.max(catalogues[t]["q"])) for t in PROFILE_TRACERS)
    q_edges = np.linspace(lower_edge, upper_edge, N_Q_BINS + 1)
    q_grid = build_q_grid(q_edges)
    write_q_edges(out_dir / f"{RUN_LABEL}_q_edges.csv", q_edges)
    readouts = full_q_readout(data_root, catalogues, q_grid, closure, args.seed, args.q_route, benchmark)
    cross = cross_tracer_readout(readouts, q_grid) if "tracers" in modes or "profile" in modes else {}
    partition_rows: list[dict[str, Any]] = []
    partition_detail: dict[str, Any] = {}
    if "partitions" in modes:
        note("r^2 partitions")
        partition_rows, partition_detail = partition_readouts(catalogues, q_grid)
    null_rows: list[dict[str, Any]] = []
    null_detail: dict[str, Any] = {}
    null_direct_max = 0.0
    if "nulls" in modes:
        null_rows, null_detail, null_direct_max = null_extension(catalogues, q_grid, closure, data_root, args.seed, args.q_route)
    write_binwise(out_dir / f"{RUN_LABEL}_q_binwise.csv", readouts)
    write_derivative(out_dir / f"{RUN_LABEL}_q_derivative.csv", readouts)
    write_rows(out_dir / f"{RUN_LABEL}_q_partitions.csv", partition_rows)
    write_rows(out_dir / f"{RUN_LABEL}_q_nulls.csv", null_rows)
    rows = metric_rows(readouts, cross, partition_rows, null_rows)
    write_rows(out_dir / f"{RUN_LABEL}_q_metrics.csv", rows)
    save_npz(out_dir / f"{RUN_LABEL}_arrays.npz", q_edges, readouts)
    empty_null = {"random_axis": {"metric_tails": {"CDF_L1": {"tail_count": None}}}, "scrambled_redshift": {"metric_tails": {"CDF_L1": {"tail_count": None}}}, "footprint_random": {"metric_tails": {"CDF_L1": {"tail_count": None}}}}
    plots = make_plots(out_dir, readouts, cross, partition_rows, null_detail if null_detail else empty_null)
    direct_data_max = max(catalogues[t]["direct_quadratic_coordinate_max_abs"] for t in PROFILE_TRACERS)
    direct_random_max = max(readouts[t]["random_direct_quadratic_coordinate_max_abs"] for t in PROFILE_TRACERS)
    fast_vs_direct_q_max_abs = max(q_route_equiv.values()) if q_route_equiv else 0.0
    q_route_equivalence = {
        "route_fast": "z^2 + C_norm2 - 2 z (n_hat dot C)",
        "route_direct": "sum((z n_hat - C)^2)",
        "q_route_requested": args.q_route,
        "full_retained_corpus_compared": bool(args.q_route == "compare"),
        "checked_batch_size_when_not_compare": 4096,
        "per_tracer_fast_vs_direct_max_abs": q_route_equiv,
        "global_fast_vs_direct_max_abs": fast_vs_direct_q_max_abs,
        "pass": bool(fast_vs_direct_q_max_abs <= 1.0e-12),
    }
    q_star_residuals = {t: readouts[t]["branches"]["q_star_bin_residual"] for t in PROFILE_TRACERS}
    q_star_derivative_residuals = {t: readouts[t]["branches"]["q_star_derivative_residual"] for t in PROFILE_TRACERS}
    finite_q_star_residuals = [abs(v) for v in q_star_residuals.values() if np.isfinite(v)]
    finite_q_star_derivative_residuals = [abs(v) for v in q_star_derivative_residuals.values() if np.isfinite(v)]
    squared_checks = {
        **high_precision,
        "fast_vs_direct_q_max_abs": fast_vs_direct_q_max_abs,
        "direct_quadratic_coordinate_max_abs": max(direct_data_max, direct_random_max, null_direct_max),
        "direct_quadratic_coordinate_data_max_abs": direct_data_max,
        "direct_quadratic_coordinate_random_max_abs": direct_random_max,
        "direct_quadratic_coordinate_null_max_abs": null_direct_max,
        "q_star_shell_residual": q_star_residuals,
        "q_star_shell_residual_max_abs": float(max(finite_q_star_residuals)) if finite_q_star_residuals else None,
        "q_star_derivative_residual": q_star_derivative_residuals,
        "q_star_derivative_residual_max_abs": float(max(finite_q_star_derivative_residuals)) if finite_q_star_derivative_residuals else None,
    }
    env = {
        "python": sys.version,
        "numpy": np.__version__,
        "mpmath": getattr(mp, "__version__", None) if mp is not None else None,
        "astropy": pkg_ver("astropy"),
        "matplotlib": pkg_ver("matplotlib"),
        "platform": platform.platform(),
        "git_commit": subprocess.run(["git", "-C", str(Path.cwd()), "rev-parse", "HEAD"], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True).stdout.strip() or None,
    }
    support = {t: readouts[t]["support"] for t in PROFILE_TRACERS}
    profile = {t: readouts[t]["metrics"] for t in PROFILE_TRACERS}
    derivative = {t: {k: readouts[t]["derivative"][k] for k in ("inner_derivative_sign_agreement", "q_star_neighbourhood_derivative", "q_star_derivative_residual", "outer_derivative_sign_agreement", "derivative_L1", "derivative_RMS", "derivative_valid_bin_count")} for t in PROFILE_TRACERS}
    branches = {t: readouts[t]["branches"] for t in PROFILE_TRACERS}
    pooled = pooled_q_readout(readouts, closure)
    retained_r_squared_support = {
        t: {
            "retained_r_squared_min": support[t]["observed_q_min"],
            "retained_r_squared_max": support[t]["observed_q_max"],
            "support_identity_mass": support[t]["support_identity_mass"],
            "support_contains_r_star_squared": support[t]["support_contains_q_star"],
            "support_contains_r_tail_squared": support[t]["support_contains_q_tail"],
            "outside_wording": "outside retained DESI r^2 support",
        }
        for t in PROFILE_TRACERS
    }
    r_squared_profile = {
        t: {
            "retained_count": profile[t]["retained_count"],
            "random_count": profile[t]["random_count"],
            "retained_r_squared_min": profile[t]["observed_q_min"],
            "retained_r_squared_max": profile[t]["observed_q_max"],
            "CDF_L1": profile[t]["CDF_L1"],
            "CDF_L2": profile[t]["CDF_L2"],
            "KS_max": profile[t]["KS_max"],
            "r_squared_mass_L1": profile[t]["q_mass_L1"],
            "r_squared_mass_L2": profile[t]["q_mass_L2"],
            "r_squared_density_RMS": profile[t]["q_density_RMS"],
            "weighted_r_squared_density_RMS": profile[t]["weighted_q_density_RMS"],
        }
        for t in PROFILE_TRACERS
    }
    r_squared_derivative = {
        t: {
            "inner_derivative_sign_agreement": derivative[t]["inner_derivative_sign_agreement"],
            "r_star_squared_neighbourhood_derivative": derivative[t]["q_star_neighbourhood_derivative"],
            "r_star_squared_derivative_residual": derivative[t]["q_star_derivative_residual"],
            "outer_derivative_sign_agreement": derivative[t]["outer_derivative_sign_agreement"],
            "derivative_L1": derivative[t]["derivative_L1"],
            "derivative_RMS": derivative[t]["derivative_RMS"],
            "derivative_valid_bin_count": derivative[t]["derivative_valid_bin_count"],
        }
        for t in PROFILE_TRACERS
    }
    r_squared_branches = {
        t: {
            "r_star_squared": branches[t]["q_star"],
            "rho_at_r_star_squared": branches[t]["rho_q_q_star"],
            "F_at_r_star_squared": branches[t]["F_q_q_star"],
            "r_star_squared_bin_index": branches[t]["q_star_bin_index"],
            "r_star_squared_bin_left": branches[t]["q_star_bin_left"],
            "r_star_squared_bin_right": branches[t]["q_star_bin_right"],
            "r_star_squared_bin_residual": branches[t]["q_star_bin_residual"],
            "r_tail_squared_supported": branches[t]["support_contains_q_tail"],
        }
        for t in PROFILE_TRACERS
    }
    benchmark["total_elapsed_seconds"] = float(time.perf_counter() - started)
    benchmark["peak_rss_mb"] = rss_mb()
    for tracer in PROFILE_TRACERS:
        benchmark["cache_hit_miss_by_tracer"][tracer] = {
            "data_cache_hit": bool(cache_report["data_unit_vectors"][tracer].get("cache_hit")),
            "random_cache_hit": bool(cache_report["random_unit_vectors"][tracer].get("cache_hit")),
            "data_cache_disabled": bool(cache_report["data_unit_vectors"][tracer].get("cache_disabled", False)),
            "random_cache_disabled": bool(cache_report["random_unit_vectors"][tracer].get("cache_disabled", False)),
        }
    manifest = {
        "run_label": RUN_LABEL,
        "script": SCRIPT_NAME,
        "command": " ".join(sys.argv),
        "run_timestamp_utc": stamp(),
        "elapsed_seconds": benchmark["total_elapsed_seconds"],
        "completed_modes": sorted(modes),
        "execution_engineering": {
            "q_route": args.q_route,
            "r_squared_execution_route": args.q_route,
            "cache_schema_version": CACHE_SCHEMA_VERSION,
            "cache_directory": str(cache_dir),
            "rebuild_cache": bool(args.rebuild_cache),
            "disable_cache": bool(args.disable_cache),
            "fast_q_formula": "q = z*z + C_norm2 - 2*z*(n_hat @ C)",
            "direct_q_formula": "q = (z*n_x-C_x)^2 + (z*n_y-C_y)^2 + (z*n_z-C_z)^2",
        },
        "q_route_equivalence": q_route_equivalence,
        "benchmark": benchmark,
        "cache": cache_report,
        "q_grid_shared": {
            "n_q_bins": N_Q_BINS,
            "q_edge_min": float(q_edges[0]),
            "q_edge_max": float(q_edges[-1]),
            "q_edges_sha256": sha256_array(q_edges),
            "F_q_edges_sha256": sha256_array(q_grid["F_q_edges"]),
            "rho_q_centers_sha256": sha256_array(q_grid["rho_q_centers"]),
            "d_sigma_q_dq_centers_sha256": sha256_array(q_grid["d_sigma_q_dq_centers"]),
            "d2_sigma_q_dq2_centers_sha256": sha256_array(q_grid["d2_sigma_q_dq2_centers"]),
            "opening_branch_sha256": sha256_array(q_grid["opening_branch"]),
            "outer_exp_readout_sha256": sha256_array(q_grid["outer_exp_readout"]),
        },
        "source_sha256": source_hash,
        "completed_closure_sha256": closure["sha256"],
        "identity_q": {
            "MU": MU,
            "GAMMA": GAMMA,
            "sigma_q": "log(MU^2 * (1 + GAMMA*q) / (MU + GAMMA)) - MU*q",
            "rho_q": "MU^2 * (1 + GAMMA*q) * exp(-MU*q) / (MU + GAMMA)",
            "F_q": "1 - exp(-MU*q) * (MU + GAMMA * (1 + MU*q)) / (MU + GAMMA)",
            "d_sigma_q_dq": "GAMMA / (1 + GAMMA*q) - MU",
            "d2_sigma_q_dq2": "-GAMMA^2 / (1 + GAMMA*q)^2",
            "q_star": q_star,
            "q_tail": q_tail,
        },
        "identity_r_squared": {
            "MU": MU,
            "GAMMA": GAMMA,
            "sigma_r_squared": "log(MU^2 * (1 + GAMMA*r^2) / (MU + GAMMA)) - MU*r^2",
            "rho_r_squared": "exp(sigma(r^2))",
            "r_star_squared": q_star,
            "r_tail_squared": q_tail,
            "r_squared_definition": "||Z*n_hat - C||^2",
        },
        "squared_radial_closure_checks": squared_checks,
        "no_legacy_radius_path": no_legacy,
        "completed_closure_input": {
            "path": closure["path"],
            "coordinate_contract": closure["coordinate_contract"],
            "source_variable": closure["source_variable"],
            "source_sha256": closure["source_sha256"],
            "direction_longitude_deg": closure["direction_longitude_deg"],
            "direction_latitude_deg": closure["direction_latitude_deg"],
            "retained_dipole_axis_unit_vector": closure["retained_dipole_axis_unit_vector"],
            "retained_dipole_axis_unit_norm2": closure["retained_dipole_axis_unit_norm2"],
            "best_scan_offset": closure["best_scan_offset"],
            "scan_coordinate": closure["scan_coordinate"],
            "scan_coordinate_norm2": closure["scan_coordinate_norm2"],
            "high_z_tracer_set": closure["high_z_tracer_set"],
            "completed_metrics": closure["completed_metrics"],
        },
        "input_hashes": input_hashes,
        "retained_q_support": support,
        "retained_r_squared_support": retained_r_squared_support,
        "pooled_q_readout": pooled,
        "q_profile_metrics": profile,
        "r_squared_profile_metrics": r_squared_profile,
        "q_derivative_metrics": derivative,
        "r_squared_derivative_metrics": r_squared_derivative,
        "q_branch_metrics": branches,
        "r_squared_branch_metrics": r_squared_branches,
        "q_cross_tracer_readout": cross,
        "q_partition_row_count": len(partition_rows),
        "q_partition_detail": {k: {kk: {"metrics": vv["metrics"], "derivative": {dk: vv["derivative"][dk] for dk in ("derivative_L1", "derivative_RMS", "inner_derivative_sign_agreement", "outer_derivative_sign_agreement", "q_star_derivative_residual")}, "branches": vv["branches"]} for kk, vv in v.items()} for k, v in partition_detail.items()},
        "q_null_row_count": len(null_rows),
        "q_null_extension": null_detail,
        "false_neighbour_separations_q_space": {
            "completed_linear_r_star_squared_minus_q_star": float((q_star * q_star) - q_star),
            "completed_linear_r_tail_squared_minus_q_star": float((q_tail * q_tail) - q_star),
            "completed_linear_r_tail_squared_minus_q_tail": float((q_tail * q_tail) - q_tail),
        },
        "artifacts": {
            "q_edges_csv": str(out_dir / f"{RUN_LABEL}_q_edges.csv"),
            "q_binwise_csv": str(out_dir / f"{RUN_LABEL}_q_binwise.csv"),
            "q_derivative_csv": str(out_dir / f"{RUN_LABEL}_q_derivative.csv"),
            "q_partitions_csv": str(out_dir / f"{RUN_LABEL}_q_partitions.csv"),
            "q_nulls_csv": str(out_dir / f"{RUN_LABEL}_q_nulls.csv"),
            "arrays_npz": str(out_dir / f"{RUN_LABEL}_arrays.npz"),
            "plots": [str(p) for p in plots],
        },
        "environment": env,
    }
    manifest_path = out_dir / f"{RUN_LABEL}_manifest.json"
    report_path = out_dir / f"{RUN_LABEL}_report.md"
    write_json(manifest_path, manifest)
    write_report(report_path, manifest)
    hashes = output_hashes(out_dir, cache)
    write_json(out_dir / f"{RUN_LABEL}_output_hashes.json", hashes)
    manifest["output_hashes"] = hashes
    write_json(manifest_path, manifest)
    final_manifest_sha256 = sha256_file(manifest_path, None)
    print_terminal_summary(args, manifest, manifest_path, final_manifest_sha256)


if __name__ == "__main__":
    main()
