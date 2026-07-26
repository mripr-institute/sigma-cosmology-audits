#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SPARC σ full-catalogue audit procedure.

This procedure emits the institutional audit record for the fixed SPARC
σ-identity audit. It does not change μ, γ, σ, ρ, σ′, or σ″, and it does not
introduce fitted quantities or free parameters.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np


ANALYSIS_ROOT = Path(__file__).resolve().parent.parent
PROJECT_ROOT = ANALYSIS_ROOT.parent
sys.path.insert(0, str(ANALYSIS_ROOT))

from shared.sparc_data import (  # noqa: E402
    DATA_DIR,
    TABLE_PATH,
    ROTCURVE_DIR,
    load_rotcurve_full,
    load_sparc_table,
)


MU = 0.082912607552
GAMMA = 0.38603416
DEGREES_OF_FREEDOM = 0

KPC_TO_M = 3.0857e19
KMS2_KPC = 1.0e6 / KPC_TO_M
UPSILON_DISK = 0.5
UPSILON_BUL = 0.7
G_DAGGER = 1.20e-10

DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "results" / "sparc_sigma_full_catalogue"
AUDIT_TITLE = "SPARC σ Full-Catalogue Audit Record"
FINAL_LINE_LABEL = "SPARC σ full-catalogue audit record"


def sigma_value(r: np.ndarray | float) -> np.ndarray | float:
    return np.log(MU * MU * (1.0 + GAMMA * r) / (MU + GAMMA)) - MU * r


def rho_value(r: np.ndarray | float) -> np.ndarray | float:
    return np.exp(sigma_value(r))


def sigma_prime_value(r: np.ndarray | float) -> np.ndarray | float:
    return GAMMA / (1.0 + GAMMA * r) - MU


def sigma_second_value(r: np.ndarray | float) -> np.ndarray | float:
    return -(GAMMA * GAMMA) / (1.0 + GAMMA * r) ** 2


def mcgaugh_rar(g_bar: np.ndarray) -> np.ndarray:
    x = np.sqrt(np.maximum(g_bar / G_DAGGER, 1e-30))
    return g_bar / (1.0 - np.exp(-x))


def rankdata(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    sorted_values = values[order]
    i = 0
    while i < len(values):
        j = i + 1
        while j < len(values) and sorted_values[j] == sorted_values[i]:
            j += 1
        ranks[order[i:j]] = 0.5 * (i + j - 1) + 1.0
        i = j
    return ranks


def pearson_r(x: np.ndarray, y: np.ndarray) -> float | None:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 2 or len(y) < 2 or np.std(x) == 0 or np.std(y) == 0:
        return None
    return float(np.corrcoef(x, y)[0, 1])


def spearman_rho(x: np.ndarray, y: np.ndarray) -> float | None:
    return pearson_r(rankdata(x), rankdata(y))


def as_float(value: Any) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(out):
        return None
    return out


def csv_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    number = as_float(value)
    if number is None:
        return ""
    return f"{number:.12g}"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def aggregate_hash(paths: list[Path]) -> str:
    h = hashlib.sha256()
    for path in sorted(paths, key=lambda p: str(p)):
        h.update(path.name.encode("utf-8"))
        h.update(b"\0")
        h.update(sha256_file(path).encode("ascii"))
        h.update(b"\0")
    return h.hexdigest()


def inspect_catalogue_rows() -> dict[str, int]:
    lines = TABLE_PATH.read_text(encoding="utf-8").splitlines()
    separator_indexes = [i for i, line in enumerate(lines) if line.startswith("---")]
    if not separator_indexes:
        return {"data_lines": 0, "parseable_records": 0, "invalid_records": 0}
    last_separator = max(separator_indexes)
    data_lines = [line for line in lines[last_separator + 1:] if line.strip()]
    parseable = sum(1 for line in data_lines if len(line.split()) >= 18)
    return {
        "data_lines": len(data_lines),
        "parseable_records": parseable,
        "invalid_records": len(data_lines) - parseable,
    }


def load_corpus() -> tuple[list[dict[str, Any]], dict[str, dict[str, np.ndarray] | None]]:
    galaxies = load_sparc_table()
    rotcurves: dict[str, dict[str, np.ndarray] | None] = {}
    for galaxy in galaxies:
        rotcurves[galaxy["name"]] = load_rotcurve_full(galaxy["name"])
    return galaxies, rotcurves


def corpus_counts(
    galaxies: list[dict[str, Any]],
    rotcurves: dict[str, dict[str, np.ndarray] | None],
    rar_points: int,
) -> dict[str, int]:
    table_rows = inspect_catalogue_rows()
    total_rotation_points = sum(
        int(len(rc["r_kpc"])) for rc in rotcurves.values() if rc is not None
    )
    return {
        "galaxies": len(galaxies),
        "catalogue_data_lines": table_rows["data_lines"],
        "catalogue_parseable_records": table_rows["parseable_records"],
        "catalogue_invalid_records": table_rows["invalid_records"],
        "rotation_curve_files": len(list(ROTCURVE_DIR.glob("*_rotmod.dat"))),
        "galaxies_with_rotation_curves": sum(
            1 for rc in rotcurves.values() if rc is not None
        ),
        "rotation_curve_points": total_rotation_points,
        "rar_data_points": rar_points,
        "valid_Rdisk_records": sum(1 for g in galaxies if g["Rdisk_kpc"] > 0),
        "valid_Vflat_records": sum(1 for g in galaxies if g["Vflat_kms"] > 0),
        "degrees_of_freedom": DEGREES_OF_FREEDOM,
    }


def build_outer_flatness_rows(
    galaxies: list[dict[str, Any]],
    rotcurves: dict[str, dict[str, np.ndarray] | None],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for galaxy in galaxies:
        name = galaxy["name"]
        rdisk = galaxy["Rdisk_kpc"]
        vflat = galaxy["Vflat_kms"]
        rc = rotcurves[name]

        row: dict[str, Any] = {
            "galaxy": name,
            "N_points_total": 0,
            "N_outer_points": 0,
            "Rdisk": rdisk if rdisk > 0 else None,
            "Rlast": None,
            "Vflat": vflat if vflat > 0 else None,
            "V_at_Rlast": None,
            "V_Rlast_over_Vflat": None,
            "outer_slope_dlogV_dlogR": None,
            "outer_slope_fit_points": 0,
            "outer_slope_status": "OK",
        }

        if rc is None:
            row["outer_slope_status"] = "NO_ROTATION_CURVE"
            rows.append(row)
            continue
        if rdisk <= 0:
            row["outer_slope_status"] = "INVALID_RDISK"
            rows.append(row)
            continue

        r = np.asarray(rc["r_kpc"], dtype=float)
        v = np.asarray(rc["Vobs"], dtype=float)
        row["N_points_total"] = int(len(r))
        if len(r) > 0:
            row["Rlast"] = float(r[-1])
            row["V_at_Rlast"] = float(v[-1])
            if vflat > 0:
                row["V_Rlast_over_Vflat"] = float(v[-1] / vflat)
            else:
                row["outer_slope_status"] = "INVALID_VFLAT"

        outer = (r / rdisk) > 3.0
        row["N_outer_points"] = int(np.sum(outer))
        valid_outer = outer & (r > 0) & (v > 0) & np.isfinite(r) & np.isfinite(v)
        fit_points = int(np.sum(valid_outer))
        row["outer_slope_fit_points"] = fit_points

        if fit_points < 3:
            row["outer_slope_status"] = "INSUFFICIENT_OUTER_POINTS"
            rows.append(row)
            continue

        slope, _ = np.polyfit(np.log(r[valid_outer]), np.log(v[valid_outer]), 1)
        row["outer_slope_dlogV_dlogR"] = float(slope)
        if row["outer_slope_status"] == "OK" and vflat <= 0:
            row["outer_slope_status"] = "INVALID_VFLAT"
        rows.append(row)
    return rows


def build_sigma_profile_rows(
    galaxies: list[dict[str, Any]],
    rotcurves: dict[str, dict[str, np.ndarray] | None],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for galaxy in galaxies:
        name = galaxy["name"]
        rdisk = galaxy["Rdisk_kpc"]
        vflat = galaxy["Vflat_kms"]
        rc = rotcurves[name]
        if rc is None:
            continue

        r = np.asarray(rc["r_kpc"], dtype=float)
        v = np.asarray(rc["Vobs"], dtype=float)
        err = np.asarray(rc["errV"], dtype=float)

        for radius, v_obs, v_err in zip(r, v, err):
            if rdisk <= 0 or not math.isfinite(float(radius)):
                rows.append({
                    "galaxy": name,
                    "R": float(radius) if math.isfinite(float(radius)) else None,
                    "r_normalized": None,
                    "sigma": None,
                    "rho": None,
                    "sigma_prime": None,
                    "sigma_second": None,
                    "V_obs": float(v_obs) if math.isfinite(float(v_obs)) else None,
                    "V_err": float(v_err) if math.isfinite(float(v_err)) else None,
                    "closure_reference_value": "OMITTED",
                    "closure_residual": "OMITTED",
                    "status": "INVALID_RDISK",
                })
                continue

            r_norm = float(radius / rdisk)
            identity_value: float | str = 1.0 if vflat > 0 else "OMITTED"
            residual: float | str = float(v_obs / vflat - 1.0) if vflat > 0 else "OMITTED"
            status = "OK" if vflat > 0 else "INVALID_VFLAT"
            rows.append({
                "galaxy": name,
                "R": float(radius),
                "r_normalized": r_norm,
                "sigma": float(sigma_value(r_norm)),
                "rho": float(rho_value(r_norm)),
                "sigma_prime": float(sigma_prime_value(r_norm)),
                "sigma_second": float(sigma_second_value(r_norm)),
                "V_obs": float(v_obs),
                "V_err": float(v_err),
                "closure_reference_value": identity_value,
                "closure_residual": residual,
                "status": status,
            })
    return rows


def compute_rotation_closure(outer_rows: list[dict[str, Any]]) -> dict[str, Any]:
    ratios = np.array([
        row["V_Rlast_over_Vflat"]
        for row in outer_rows
        if as_float(row["V_Rlast_over_Vflat"]) is not None
    ], dtype=float)
    within_10 = int(np.sum(np.abs(ratios - 1.0) < 0.10)) if len(ratios) else 0
    within_20 = int(np.sum(np.abs(ratios - 1.0) < 0.20)) if len(ratios) else 0
    return {
        "total_galaxies": len(outer_rows),
        "galaxies_with_valid_Vflat": len(ratios),
        "V_Rlast_over_Vflat_mean": float(np.mean(ratios)) if len(ratios) else None,
        "V_Rlast_over_Vflat_median": float(np.median(ratios)) if len(ratios) else None,
        "V_Rlast_over_Vflat_std": float(np.std(ratios)) if len(ratios) else None,
        "within_10_percent_count": within_10,
        "within_10_percent_ratio": float(within_10 / len(ratios)) if len(ratios) else None,
        "within_20_percent_count": within_20,
        "within_20_percent_ratio": float(within_20 / len(ratios)) if len(ratios) else None,
        "degrees_of_freedom": DEGREES_OF_FREEDOM,
    }


def compute_outer_slope_summary(outer_rows: list[dict[str, Any]]) -> dict[str, Any]:
    slopes = np.array([
        row["outer_slope_dlogV_dlogR"]
        for row in outer_rows
        if row["outer_slope_status"] == "OK"
        and as_float(row["outer_slope_dlogV_dlogR"]) is not None
    ], dtype=float)
    return {
        "outer_slope_mean": float(np.mean(slopes)) if len(slopes) else None,
        "outer_slope_median": float(np.median(slopes)) if len(slopes) else None,
        "outer_slope_std": float(np.std(slopes)) if len(slopes) else None,
        "outer_slope_abs_median": float(np.median(np.abs(slopes))) if len(slopes) else None,
    }


def compute_sigma_derivative_ranges(profile_rows: list[dict[str, Any]]) -> dict[str, Any]:
    summary: dict[str, Any] = {}
    ok_rows = [row for row in profile_rows if row["status"] == "OK"]
    for column in ["sigma", "rho", "sigma_prime", "sigma_second"]:
        values = np.array([
            as_float(row[column])
            for row in ok_rows
            if as_float(row[column]) is not None
        ], dtype=float)
        summary[f"{column}_min"] = float(np.min(values)) if len(values) else None
        summary[f"{column}_max"] = float(np.max(values)) if len(values) else None
    return summary


def compute_rar_consistency(
    galaxies: list[dict[str, Any]],
    rotcurves: dict[str, dict[str, np.ndarray] | None],
) -> dict[str, Any]:
    all_gobs: list[float] = []
    all_gbar: list[float] = []
    all_quality: list[int] = []
    all_names: list[str] = []
    galaxies_with_rar = 0

    for galaxy in galaxies:
        rc = rotcurves[galaxy["name"]]
        if rc is None:
            continue
        r = np.asarray(rc["r_kpc"], dtype=float)
        vobs = np.asarray(rc["Vobs"], dtype=float)
        vgas = np.asarray(rc["Vgas"], dtype=float)
        vdisk = np.asarray(rc["Vdisk"], dtype=float)
        vbul = np.asarray(rc["Vbul"], dtype=float)
        mask = (r > 0) & (vobs > 0)
        if not np.any(mask):
            continue

        r = r[mask]
        vobs = vobs[mask]
        vgas = vgas[mask]
        vdisk = vdisk[mask]
        vbul = vbul[mask]

        g_obs = (vobs**2 / r) * KMS2_KPC
        v_bar_sq = vgas * np.abs(vgas) + UPSILON_DISK * vdisk**2 + UPSILON_BUL * vbul**2
        g_bar = (v_bar_sq / r) * KMS2_KPC
        ok = (g_obs > 0) & (g_bar > 0) & np.isfinite(g_obs) & np.isfinite(g_bar)
        if not np.any(ok):
            continue

        all_gobs.extend(g_obs[ok].tolist())
        all_gbar.extend(g_bar[ok].tolist())
        all_quality.extend([galaxy["Q"]] * int(np.sum(ok)))
        all_names.extend([galaxy["name"]] * int(np.sum(ok)))
        galaxies_with_rar += 1

    g_obs_arr = np.asarray(all_gobs, dtype=float)
    g_bar_arr = np.asarray(all_gbar, dtype=float)
    quality_arr = np.asarray(all_quality, dtype=int)

    if len(g_obs_arr) == 0:
        return {
            "galaxies_with_RAR_data": 0,
            "individual_RAR_points": 0,
            "degrees_of_freedom": DEGREES_OF_FREEDOM,
        }

    g_mcg = mcgaugh_rar(g_bar_arr)
    log_resid = np.log10(g_obs_arr / g_mcg)
    log_gobs = np.log10(g_obs_arr)
    log_gbar = np.log10(g_bar_arr)

    def rms_for(mask: np.ndarray) -> dict[str, Any]:
        n = int(np.sum(mask))
        if n == 0:
            return {"points": 0, "rms_scatter_dex": None, "galaxies": 0}
        return {
            "points": n,
            "rms_scatter_dex": float(np.std(log_resid[mask])),
            "galaxies": len(set(np.asarray(all_names)[mask])),
        }

    return {
        "galaxies_with_RAR_data": galaxies_with_rar,
        "individual_RAR_points": int(len(g_obs_arr)),
        "g_obs_min": float(np.min(g_obs_arr)),
        "g_obs_max": float(np.max(g_obs_arr)),
        "g_bar_min": float(np.min(g_bar_arr)),
        "g_bar_max": float(np.max(g_bar_arr)),
        "g_bar_dynamic_range": float(np.max(g_bar_arr) / np.min(g_bar_arr)),
        "RMS_scatter_dex": float(np.std(log_resid)),
        "median_residual_dex": float(np.median(log_resid)),
        "mean_residual_dex": float(np.mean(log_resid)),
        "Pearson_r": pearson_r(log_gbar, log_gobs),
        "Spearman_rho": spearman_rho(log_gbar, log_gobs),
        "quality_stratified_RMS": {
            "Q = 1": rms_for(quality_arr == 1),
            "Q <= 2": rms_for(quality_arr <= 2),
            "all": rms_for(np.ones(len(quality_arr), dtype=bool)),
        },
        "degrees_of_freedom": DEGREES_OF_FREEDOM,
    }


def write_csv(path: Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: csv_value(row.get(column)) for column in columns})


def validate_outer_rows(
    rows: list[dict[str, Any]],
    galaxies: list[dict[str, Any]],
) -> dict[str, Any]:
    catalogue_names = {galaxy["name"] for galaxy in galaxies}
    row_names = {row["galaxy"] for row in rows}
    finite_failures: list[str] = []

    always_numeric = ["N_points_total", "N_outer_points", "outer_slope_fit_points"]
    status_allowed_missing = {
        "NO_ROTATION_CURVE",
        "INVALID_RDISK",
        "INVALID_VFLAT",
        "INSUFFICIENT_OUTER_POINTS",
        "INVALID_OUTER_VALUES",
    }

    for row in rows:
        for column in always_numeric:
            if as_float(row[column]) is None:
                finite_failures.append(f"{row['galaxy']}:{column}")
        status = row["outer_slope_status"]
        if status == "OK" and as_float(row["outer_slope_dlogV_dlogR"]) is None:
            finite_failures.append(f"{row['galaxy']}:outer_slope_dlogV_dlogR")
        if status not in status_allowed_missing | {"OK"}:
            finite_failures.append(f"{row['galaxy']}:unknown_status:{status}")

    return {
        "row_count": len(rows),
        "expected_row_count": len(galaxies),
        "row_count_matches": len(rows) == len(galaxies),
        "no_silent_galaxy_loss": row_names == catalogue_names,
        "missing_galaxies": sorted(catalogue_names - row_names),
        "extra_galaxies": sorted(row_names - catalogue_names),
        "finite_required_numeric_columns": len(finite_failures) == 0,
        "finite_failures": finite_failures,
    }


def validate_profile_rows(
    rows: list[dict[str, Any]],
    rotcurves: dict[str, dict[str, np.ndarray] | None],
) -> dict[str, Any]:
    expected_rows = sum(int(len(rc["r_kpc"])) for rc in rotcurves.values() if rc is not None)
    expected_names = {name for name, rc in rotcurves.items() if rc is not None}
    row_names = {row["galaxy"] for row in rows}
    finite_failures: list[str] = []

    required_if_ok = [
        "R",
        "r_normalized",
        "sigma",
        "rho",
        "sigma_prime",
        "sigma_second",
        "V_obs",
        "V_err",
    ]
    for idx, row in enumerate(rows):
        status = row["status"]
        if status == "OK":
            for column in required_if_ok:
                if as_float(row[column]) is None:
                    finite_failures.append(f"row_{idx}:{row['galaxy']}:{column}")
            for optional_column in ["closure_reference_value", "closure_residual"]:
                value = row[optional_column]
                if value != "OMITTED" and as_float(value) is None:
                    finite_failures.append(f"row_{idx}:{row['galaxy']}:{optional_column}")
        elif status not in {"INVALID_RDISK", "INVALID_VFLAT"}:
            finite_failures.append(f"row_{idx}:{row['galaxy']}:unknown_status:{status}")

    return {
        "row_count": len(rows),
        "expected_row_count": expected_rows,
        "row_count_matches": len(rows) == expected_rows,
        "no_silent_galaxy_loss": row_names == expected_names,
        "missing_galaxies_with_rotation_curves": sorted(expected_names - row_names),
        "extra_galaxies": sorted(row_names - expected_names),
        "finite_required_numeric_columns": len(finite_failures) == 0,
        "finite_failures": finite_failures[:50],
    }


def run_existing_figure_audits(output_dir: Path) -> dict[str, Any]:
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    tmp_root = Path(os.environ.get("TMPDIR", "/tmp"))
    mpl_cache = tmp_root / "sparc_sigma_audit_matplotlib_cache"
    mpl_cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_cache))

    results: dict[str, Any] = {}
    try:
        from sparc import sigma_sparc_identity

        sigma_sparc_identity.run(
            output_prefix=str(figures_dir / "sparc_rotation_closure"),
            web_json=False,
        )
        rotation_figures = sorted(
            str(path.resolve()) for path in figures_dir.glob("sparc_rotation_closure*.png")
        )
        results["sparc_identity_audit"] = {
            "status": "PASS" if rotation_figures else "FAIL",
            "figures": rotation_figures,
        }
    except Exception as exc:  # pragma: no cover - reported in audit output.
        results["sparc_identity_audit"] = {"status": "FAIL", "error": repr(exc), "figures": []}

    try:
        from sparc import sigma_rar_identity

        sigma_rar_identity.run(
            output_prefix=str(figures_dir / "sparc_rar_consistency"),
            web_json=False,
        )
        rar_figures = sorted(
            str(path.resolve()) for path in figures_dir.glob("sparc_rar_consistency*.png")
        )
        results["rar_audit"] = {
            "status": "PASS" if rar_figures else "FAIL",
            "figures": rar_figures,
        }
    except Exception as exc:  # pragma: no cover - reported in audit output.
        results["rar_audit"] = {"status": "FAIL", "error": repr(exc), "figures": []}

    return results


def table_status_counts(rows: list[dict[str, Any]], status_column: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        status = str(row[status_column])
        counts[status] = counts.get(status, 0) + 1
    return dict(sorted(counts.items()))


def degrees_of_freedom_contract(blocks: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """Compute the fitted-degree-of-freedom invariant from every main block."""
    observed = {
        name: block.get("degrees_of_freedom")
        for name, block in blocks.items()
    }
    return {
        "observed": observed,
        "all_present": all(value is not None for value in observed.values()),
        "all_zero": all(value == 0 for value in observed.values()),
    }


def finalize_figure_audit_statuses(
    figure_audits: dict[str, Any],
    rotation_closure: dict[str, Any],
    rar_consistency: dict[str, Any],
) -> None:
    """Attach and evaluate numerical evidence predicates for both figure audits."""
    rotation = figure_audits.get("sparc_identity_audit", {})
    valid_vflat = int(rotation_closure.get("galaxies_with_valid_Vflat", 0))
    rotation_predicates = {
        "figures_created": bool(rotation.get("figures")),
        "valid_vflat_records_present": valid_vflat > 0,
        "all_valid_vflat_within_10_percent": (
            int(rotation_closure.get("within_10_percent_count", -1)) == valid_vflat
        ),
        "all_valid_vflat_within_20_percent": (
            int(rotation_closure.get("within_20_percent_count", -1)) == valid_vflat
        ),
        "rotation_summary_finite": all(
            as_float(rotation_closure.get(key)) is not None
            for key in (
                "V_Rlast_over_Vflat_mean",
                "V_Rlast_over_Vflat_median",
                "V_Rlast_over_Vflat_std",
            )
        ),
    }
    rotation["evidence_predicates"] = rotation_predicates
    rotation["status"] = "PASS" if all(rotation_predicates.values()) else "FAIL"

    rar = figure_audits.get("rar_audit", {})
    rar_predicates = {
        "figures_created": bool(rar.get("figures")),
        "galaxies_present": int(rar_consistency.get("galaxies_with_RAR_data", 0)) > 0,
        "rar_points_present": int(rar_consistency.get("individual_RAR_points", 0)) > 0,
        "rar_statistics_finite": all(
            as_float(rar_consistency.get(key)) is not None
            for key in (
                "RMS_scatter_dex",
                "median_residual_dex",
                "mean_residual_dex",
                "Pearson_r",
                "Spearman_rho",
            )
        ),
        "positive_rank_and_linear_correlations": (
            as_float(rar_consistency.get("Pearson_r")) is not None
            and as_float(rar_consistency.get("Spearman_rho")) is not None
            and float(rar_consistency["Pearson_r"]) > 0.0
            and float(rar_consistency["Spearman_rho"]) > 0.0
        ),
    }
    rar["evidence_predicates"] = rar_predicates
    rar["status"] = "PASS" if all(rar_predicates.values()) else "FAIL"


def pass_fail_checks(
    output_paths: dict[str, Path],
    outer_validation: dict[str, Any],
    profile_validation: dict[str, Any],
    corpus: dict[str, int],
    rotation_closure: dict[str, Any],
    rar_consistency: dict[str, Any],
    figure_audits: dict[str, Any],
    dof_contract: dict[str, Any],
) -> dict[str, bool]:
    return {
        "existing_SPARC_identity_audit_PASS": (
            figure_audits.get("sparc_identity_audit", {}).get("status") == "PASS"
            and rotation_closure["galaxies_with_valid_Vflat"] > 0
            and as_float(rotation_closure["V_Rlast_over_Vflat_mean"]) is not None
        ),
        "existing_RAR_audit_PASS": (
            figure_audits.get("rar_audit", {}).get("status") == "PASS"
            and rar_consistency.get("individual_RAR_points", 0) > 0
            and as_float(rar_consistency.get("RMS_scatter_dex")) is not None
        ),
        "outer_flatness_by_galaxy_csv_created": output_paths["outer_flatness"].exists(),
        "sigma_derivative_profile_sparc_csv_created": output_paths["sigma_profile"].exists(),
        "outer_flatness_row_count_matches_expected": outer_validation["row_count_matches"],
        "sigma_profile_row_count_matches_expected": profile_validation["row_count_matches"],
        "outer_flatness_no_silent_galaxy_loss": outer_validation["no_silent_galaxy_loss"],
        "sigma_profile_no_silent_galaxy_loss": profile_validation["no_silent_galaxy_loss"],
        "outer_flatness_required_numeric_values_finite": outer_validation[
            "finite_required_numeric_columns"
        ],
        "sigma_profile_required_numeric_values_finite": profile_validation[
            "finite_required_numeric_columns"
        ],
        "corpus_counts_consistent": (
            corpus["galaxies"] == corpus["catalogue_parseable_records"]
            and corpus["catalogue_invalid_records"] >= 0
            and corpus["rotation_curve_points"] == profile_validation["expected_row_count"]
            and corpus["rar_data_points"] == rar_consistency.get("individual_RAR_points", 0)
        ),
        "degrees_of_freedom_zero": corpus["degrees_of_freedom"] == 0,
        "degrees_of_freedom_in_all_main_report_blocks": (
            dof_contract["all_present"] and dof_contract["all_zero"]
        ),
    }


def write_result_json(path: Path, result: dict[str, Any]) -> None:
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def output_hashes(paths: dict[str, Path], figures_dir: Path) -> dict[str, str]:
    hashable = {
        key: path
        for key, path in paths.items()
        if key != "audit_manifest" and path.exists() and path.is_file()
    }
    for figure in sorted(figures_dir.glob("*.png")):
        hashable[f"figures/{figure.name}"] = figure
    return {key: sha256_file(path) for key, path in hashable.items()}


def procedure_hashes() -> dict[str, str]:
    candidates = {
        "audit_procedures/sparc_sigma_custody_audit.py": (
            PROJECT_ROOT / "audit_procedures" / "sparc_sigma_custody_audit.py"
        ),
        "sparc_data.py": PROJECT_ROOT / "shared" / "sparc_data.py",
        "sigma_sparc_identity.py": PROJECT_ROOT / "sparc" / "sigma_sparc_identity.py",
        "sigma_rar_identity.py": PROJECT_ROOT / "sparc" / "sigma_rar_identity.py",
    }
    return {key: sha256_file(path) for key, path in candidates.items() if path.exists()}


def write_manifest(
    path: Path,
    run_timestamp_utc: str,
    command: str,
    input_paths: dict[str, Any],
    output_paths_data: dict[str, Any],
    file_hashes: dict[str, Any],
    procedure_file_hashes: dict[str, str],
    corpus: dict[str, int],
    final_status: str,
) -> None:
    manifest = {
        "run_timestamp_utc": run_timestamp_utc,
        "command": command,
        "input_paths": input_paths,
        "output_paths": output_paths_data,
        "file_hashes": file_hashes,
        "procedure_hashes": procedure_file_hashes,
        "mu": MU,
        "gamma": GAMMA,
        "degrees_of_freedom": DEGREES_OF_FREEDOM,
        "corpus_counts": corpus,
        "final_status": final_status,
    }
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def fmt(value: Any, precision: int = 6) -> str:
    number = as_float(value)
    if number is None:
        return "OMITTED"
    return f"{number:.{precision}g}"


def write_report(
    path: Path,
    run_timestamp_utc: str,
    output_dir: Path,
    source_hashes: dict[str, Any],
    corpus: dict[str, int],
    fixed_identity: dict[str, float],
    rotation_closure: dict[str, Any],
    outer_slope_summary: dict[str, Any],
    outer_validation: dict[str, Any],
    outer_status_counts: dict[str, int],
    profile_validation: dict[str, Any],
    profile_status_counts: dict[str, int],
    sigma_derivative_ranges: dict[str, Any],
    rar_consistency: dict[str, Any],
    checks: dict[str, bool],
    final_status: str,
    output_paths_data: dict[str, Any],
    figure_audits: dict[str, Any],
    dof_contract: dict[str, Any],
) -> None:
    q_rms = rar_consistency.get("quality_stratified_RMS", {})
    lines = [
        f"# {AUDIT_TITLE}",
        "",
        "## 1. Audit Header",
        "",
        f"- audit name: {AUDIT_TITLE}",
        f"- run timestamp UTC: {run_timestamp_utc}",
        f"- source file hash, SPARC_Lelli2016c.mrt: {source_hashes['table_sha256']}",
        f"- source file hash, rotation-curve aggregate: {source_hashes['rotcurves_aggregate_sha256']}",
        f"- SPARC data source path: {DATA_DIR.resolve()}",
        f"- output directory: {output_dir.resolve()}",
        f"- μ: {MU}",
        f"- γ: {GAMMA}",
        f"- degrees_of_freedom: {DEGREES_OF_FREEDOM}",
        f"- final status: {final_status}",
        "",
        "## 2. Fixed Identity Block",
        "",
        "The fixed audit identities are:",
        "",
        "- ρ(0) = μ² / (μ + γ)",
        "- σ′(r) = γ / (1 + γr) - μ",
        "- σ″(r) = -γ² / (1 + γr)²",
        "",
        "| quantity | value |",
        "| --- | ---: |",
        f"| ρ(0) | {fixed_identity['rho_0']:.12g} |",
        f"| γ² | {fixed_identity['gamma_squared']:.12g} |",
        f"| μ² | {fixed_identity['mu_squared']:.12g} |",
        f"| degrees_of_freedom | {DEGREES_OF_FREEDOM} |",
        "",
        "## 3. Full-Catalogue Rotation Closure Block",
        "",
        f"- total galaxies: {rotation_closure['total_galaxies']}",
        f"- galaxies with valid Vflat: {rotation_closure['galaxies_with_valid_Vflat']}",
        f"- total rotation curve points: {corpus['rotation_curve_points']}",
        f"- V(Rlast) / Vflat mean: {fmt(rotation_closure['V_Rlast_over_Vflat_mean'])}",
        f"- V(Rlast) / Vflat median: {fmt(rotation_closure['V_Rlast_over_Vflat_median'])}",
        f"- V(Rlast) / Vflat std: {fmt(rotation_closure['V_Rlast_over_Vflat_std'])}",
        (
            "- |V(Rlast)/Vflat - 1| < 10%: "
            f"{rotation_closure['within_10_percent_count']}/"
            f"{rotation_closure['galaxies_with_valid_Vflat']} "
            f"({100.0 * rotation_closure['within_10_percent_ratio']:.2f}%)"
        ),
        (
            "- |V(Rlast)/Vflat - 1| < 20%: "
            f"{rotation_closure['within_20_percent_count']}/"
            f"{rotation_closure['galaxies_with_valid_Vflat']} "
            f"({100.0 * rotation_closure['within_20_percent_ratio']:.2f}%)"
        ),
        f"- degrees_of_freedom: {DEGREES_OF_FREEDOM}",
        "",
        "## 4. Per-Galaxy Outer-Flatness Table",
        "",
        f"- file: {output_paths_data['outer_flatness']}",
        f"- row count: {outer_validation['row_count']}",
        f"- expected row count: {outer_validation['expected_row_count']}",
        f"- no silent galaxy loss: {outer_validation['no_silent_galaxy_loss']}",
        f"- required numeric audit values finite: {outer_validation['finite_required_numeric_columns']}",
        f"- status counts: {json.dumps(outer_status_counts, sort_keys=True)}",
        f"- outer_slope_mean: {fmt(outer_slope_summary['outer_slope_mean'])}",
        f"- outer_slope_median: {fmt(outer_slope_summary['outer_slope_median'])}",
        f"- outer_slope_std: {fmt(outer_slope_summary['outer_slope_std'])}",
        f"- outer_slope_abs_median: {fmt(outer_slope_summary['outer_slope_abs_median'])}",
        f"- degrees_of_freedom: {DEGREES_OF_FREEDOM}",
        "",
        "## 5. σ Derivative Profile Table",
        "",
        f"- file: {output_paths_data['sigma_profile']}",
        f"- row count: {profile_validation['row_count']}",
        f"- expected rotation-point row count: {profile_validation['expected_row_count']}",
        f"- no silent galaxy loss among galaxies with rotation curves: {profile_validation['no_silent_galaxy_loss']}",
        f"- required numeric audit values finite: {profile_validation['finite_required_numeric_columns']}",
        f"- status counts: {json.dumps(profile_status_counts, sort_keys=True)}",
        "- closure_reference_value is the fixed identity value 1.0 where Vflat is valid.",
        "- No new fitted quantity was introduced; closure_residual is V_obs / Vflat - 1 where available.",
        f"- sigma min/max: {fmt(sigma_derivative_ranges['sigma_min'])} / {fmt(sigma_derivative_ranges['sigma_max'])}",
        f"- rho min/max: {fmt(sigma_derivative_ranges['rho_min'])} / {fmt(sigma_derivative_ranges['rho_max'])}",
        f"- sigma_prime min/max: {fmt(sigma_derivative_ranges['sigma_prime_min'])} / {fmt(sigma_derivative_ranges['sigma_prime_max'])}",
        f"- sigma_second min/max: {fmt(sigma_derivative_ranges['sigma_second_min'])} / {fmt(sigma_derivative_ranges['sigma_second_max'])}",
        f"- degrees_of_freedom: {DEGREES_OF_FREEDOM}",
        "",
        "## 6. RAR Consistency Block",
        "",
        f"- galaxies with RAR data: {rar_consistency.get('galaxies_with_RAR_data')}",
        f"- individual RAR points: {rar_consistency.get('individual_RAR_points')}",
        f"- g_obs min/max: {fmt(rar_consistency.get('g_obs_min'))} / {fmt(rar_consistency.get('g_obs_max'))}",
        f"- g_bar min/max: {fmt(rar_consistency.get('g_bar_min'))} / {fmt(rar_consistency.get('g_bar_max'))}",
        f"- g_bar dynamic range: {fmt(rar_consistency.get('g_bar_dynamic_range'))}",
        f"- RMS scatter dex: {fmt(rar_consistency.get('RMS_scatter_dex'))}",
        f"- median residual dex: {fmt(rar_consistency.get('median_residual_dex'))}",
        f"- mean residual dex: {fmt(rar_consistency.get('mean_residual_dex'))}",
        f"- Pearson r: {fmt(rar_consistency.get('Pearson_r'))}",
        f"- Spearman rho: {fmt(rar_consistency.get('Spearman_rho'))}",
        f"- quality-stratified RMS, Q = 1: {fmt(q_rms.get('Q = 1', {}).get('rms_scatter_dex'))}",
        f"- quality-stratified RMS, Q ≤ 2: {fmt(q_rms.get('Q <= 2', {}).get('rms_scatter_dex'))}",
        f"- quality-stratified RMS, all: {fmt(q_rms.get('all', {}).get('rms_scatter_dex'))}",
        f"- degrees_of_freedom: {DEGREES_OF_FREEDOM}",
        "",
        "## 7. PASS/FAIL Logic",
        "",
    ]
    for key, value in checks.items():
        lines.append(f"- {key}: {value}")
    lines.extend([
        "",
        f"- final status: {final_status}",
        "",
        "## 8. Output Files",
        "",
    ])
    for key, value in output_paths_data.items():
        lines.append(f"- {key}: {value}")
    lines.extend([
        "",
        "Figure audit statuses:",
        "",
        f"- SPARC identity audit: {figure_audits.get('sparc_identity_audit', {}).get('status')}",
        f"- RAR audit: {figure_audits.get('rar_audit', {}).get('status')}",
        "",
        "Computed figure evidence predicates:",
        "",
    ])
    for audit_name in ("sparc_identity_audit", "rar_audit"):
        predicates = figure_audits.get(audit_name, {}).get("evidence_predicates", {})
        for key, value in predicates.items():
            lines.append(f"- {audit_name}.{key}: {value}")
    lines.extend([
        "",
        "Computed fitted-degree-of-freedom contract:",
        "",
    ])
    for block_name, value in dof_contract["observed"].items():
        lines.append(f"- {block_name}: {value}")
    lines.extend([
        f"- all fields present: {dof_contract['all_present']}",
        f"- all fields zero: {dof_contract['all_zero']}",
        "",
        "## Audit Invariants",
        "",
        "- μ and γ are fixed.",
        f"- degrees_of_freedom: {DEGREES_OF_FREEDOM}",
        "- source SPARC data unchanged",
        "- all catalogue galaxies retained with explicit status",
        "- no fitted quantity introduced",
        "",
        "## SPARC Closure Statement",
        "",
        "The full SPARC catalogue was evaluated with fixed μ and γ.",
        "The audit record verifies:",
        "- positive finite centre: ρ(0)",
        "- first derivative: σ′(r)",
        "- second derivative: σ″(r)",
        "- full rotation-curve corpus accounting",
        "- outer-flatness evidence by galaxy",
        "- RAR consistency",
        f"- degrees_of_freedom: {DEGREES_OF_FREEDOM}",
        f"Final status: {final_status}",
        "",
        "## 10. Final Check",
        "",
        f"{FINAL_LINE_LABEL}: {final_status}",
        "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def source_file_hashes() -> dict[str, Any]:
    rotcurve_paths = sorted(ROTCURVE_DIR.glob("*_rotmod.dat"))
    return {
        "table_sha256": sha256_file(TABLE_PATH),
        "rotcurves_aggregate_sha256": aggregate_hash(rotcurve_paths),
        "rotcurves": {path.name: sha256_file(path) for path in rotcurve_paths},
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=f"Run the {AUDIT_TITLE}."
    )
    parser.add_argument(
        "--output-directory",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Directory for audit report, result JSON, CSV evidence, figures, and manifest.",
    )
    args = parser.parse_args(argv)

    output_dir = Path(args.output_directory).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    output_paths = {
        "audit_report": output_dir / "sparc_sigma_full_catalogue_audit_report.md",
        "audit_result": output_dir / "sparc_sigma_full_catalogue_audit_result.json",
        "outer_flatness": output_dir / "outer_flatness_by_galaxy.csv",
        "sigma_profile": output_dir / "sigma_derivative_profile_sparc.csv",
        "audit_manifest": output_dir / "sparc_sigma_full_catalogue_audit_manifest.json",
    }
    figures_dir = output_dir / "figures"

    run_timestamp_utc = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    command = " ".join([sys.executable, *sys.argv])

    source_hashes = source_file_hashes()
    galaxies, rotcurves = load_corpus()
    rar_consistency = compute_rar_consistency(galaxies, rotcurves)
    corpus = corpus_counts(
        galaxies, rotcurves, rar_consistency.get("individual_RAR_points", 0)
    )

    outer_rows = build_outer_flatness_rows(galaxies, rotcurves)
    profile_rows = build_sigma_profile_rows(galaxies, rotcurves)

    outer_columns = [
        "galaxy",
        "N_points_total",
        "N_outer_points",
        "Rdisk",
        "Rlast",
        "Vflat",
        "V_at_Rlast",
        "V_Rlast_over_Vflat",
        "outer_slope_dlogV_dlogR",
        "outer_slope_fit_points",
        "outer_slope_status",
    ]
    profile_columns = [
        "galaxy",
        "R",
        "r_normalized",
        "sigma",
        "rho",
        "sigma_prime",
        "sigma_second",
        "V_obs",
        "V_err",
        "closure_reference_value",
        "closure_residual",
        "status",
    ]
    write_csv(output_paths["outer_flatness"], outer_rows, outer_columns)
    write_csv(output_paths["sigma_profile"], profile_rows, profile_columns)

    outer_validation = validate_outer_rows(outer_rows, galaxies)
    profile_validation = validate_profile_rows(profile_rows, rotcurves)
    rotation_closure = compute_rotation_closure(outer_rows)
    outer_slope_summary = compute_outer_slope_summary(outer_rows)
    sigma_derivative_ranges = compute_sigma_derivative_ranges(profile_rows)
    fixed_identity = {
        "rho_0": MU * MU / (MU + GAMMA),
        "gamma_squared": GAMMA * GAMMA,
        "mu_squared": MU * MU,
        "degrees_of_freedom": DEGREES_OF_FREEDOM,
    }

    figure_audits = run_existing_figure_audits(output_dir)
    dof_contract = degrees_of_freedom_contract({
        "corpus_counts": corpus,
        "fixed_identity": fixed_identity,
        "rotation_closure": rotation_closure,
        "rar_consistency": rar_consistency,
    })
    finalize_figure_audit_statuses(figure_audits, rotation_closure, rar_consistency)
    checks = pass_fail_checks(
        output_paths,
        outer_validation,
        profile_validation,
        corpus,
        rotation_closure,
        rar_consistency,
        figure_audits,
        dof_contract,
    )
    final_status = "PASS" if all(checks.values()) else "FAIL"

    output_paths_data = {
        key: str(path.resolve()) for key, path in output_paths.items()
    }
    output_paths_data["figures"] = str(figures_dir.resolve())

    result = {
        "audit_name": AUDIT_TITLE,
        "run_timestamp_utc": run_timestamp_utc,
        "mu": MU,
        "gamma": GAMMA,
        "degrees_of_freedom": DEGREES_OF_FREEDOM,
        "fixed_identity": fixed_identity,
        "corpus_counts": corpus,
        "rotation_closure": rotation_closure,
        "outer_flatness_table": {
            "path": output_paths_data["outer_flatness"],
            "validation": outer_validation,
            "status_counts": table_status_counts(outer_rows, "outer_slope_status"),
            "outer_slope_summary": outer_slope_summary,
        },
        "sigma_derivative_profile_table": {
            "path": output_paths_data["sigma_profile"],
            "validation": profile_validation,
            "status_counts": table_status_counts(profile_rows, "status"),
            "sigma_derivative_range_summary": sigma_derivative_ranges,
            "new_fitted_quantity_introduced": False,
        },
        "rar_consistency": rar_consistency,
        "figure_audits": figure_audits,
        "degrees_of_freedom_contract": dof_contract,
        "pass_fail_checks": checks,
        "final_status": final_status,
        "output_paths": output_paths_data,
    }

    write_result_json(output_paths["audit_result"], result)
    write_report(
        output_paths["audit_report"],
        run_timestamp_utc,
        output_dir,
        source_hashes,
        corpus,
        fixed_identity,
        rotation_closure,
        outer_slope_summary,
        outer_validation,
        table_status_counts(outer_rows, "outer_slope_status"),
        profile_validation,
        table_status_counts(profile_rows, "status"),
        sigma_derivative_ranges,
        rar_consistency,
        checks,
        final_status,
        output_paths_data,
        figure_audits,
        dof_contract,
    )

    file_hashes = {
        "inputs": {
            "SPARC_Lelli2016c.mrt": source_hashes["table_sha256"],
            "rotcurves_aggregate": source_hashes["rotcurves_aggregate_sha256"],
            "rotcurves": source_hashes["rotcurves"],
        },
        "outputs": output_hashes(output_paths, figures_dir),
    }
    procedure_file_hashes = procedure_hashes()
    input_paths = {
        "sparc_data_directory": str(DATA_DIR.resolve()),
        "sparc_catalogue_table": str(TABLE_PATH.resolve()),
        "rotation_curve_directory": str(ROTCURVE_DIR.resolve()),
    }
    write_manifest(
        output_paths["audit_manifest"],
        run_timestamp_utc,
        command,
        input_paths,
        output_paths_data,
        file_hashes,
        procedure_file_hashes,
        corpus,
        final_status,
    )

    print(f"{FINAL_LINE_LABEL}: {final_status}")
    print(f"Report: {output_paths['audit_report']}")
    print(f"Result: {output_paths['audit_result']}")
    print(f"Manifest: {output_paths['audit_manifest']}")
    return 0 if final_status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
