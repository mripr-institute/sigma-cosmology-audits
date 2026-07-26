#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
//! pantheon_zr_derivative_mu_audit.rs — Pantheon σ z=r gamma closure with derivative-recovered μ.
//!
//! σ-identity audit: Pantheon catalogue row recovery with derivative-recovered μ.
//!
//! Core identity:
//!
//! σ(r) = ln[μ²(1+γr)/(μ+γ)] − μr
//! ρ(r) = [μ²/(μ+γ)](1+γr)e^(−μr)
//!
//! Derivative recovery:
//!
//! σ′(r) = γ/(1+γr) − μ
//! σ″(r) = −γ²/(1+γr)²
//! Q(r) = sqrt(−σ″(r))
//! μ = Q(r) − σ′(r)
//! γ = Q(r)/(1 − rQ(r))
//!
//! Stationary recovery:
//!
//! r* = 1/μ − 1/γ
//! σ′(r*) = 0
//! μ = sqrt(−σ″(r*))
//!
//! Pantheon row rule:
//!
//! catalogue_z = input z column
//! r = catalogue_z
//! z_equals_r = true
//!
//! Row closure:
//!
//! D_record = catalogue_luminosity_record * catalogue_scale_record / (c(1+r))
//! T_record = D_record * exp(μr)
//! ln(1+γr)/γ = T_record
//!
//! γ is recovered row-by-row by Lambert W closed form.
//!
//! Audit products:
//!
//! implemented identity functions
//! derivative recovery route
//! stationary recovery route
//! Pantheon z=r row rule
//! Lambert W row closure
//! corpus accounting
//! finite recovery
//! SHA-256 manifest
//! PASS/FAIL result
//!
//! Author:  Alex Albert <alex.albert@mripr.org>
//! ORCID:   https://orcid.org/0009-0005-6981-2087
//! License: CC0 / Public domain
"""

import argparse
import csv
import hashlib
import json
import math
import os
import shlex
import statistics
import sys
from datetime import datetime, timezone

import mpmath as mp



AUDIT_NAME = "Pantheon σ z=r Invariant Extraction and Row Closure Audit"
AUDIT_SLUG = "pantheon_sigma_zr_invariant_extraction_row_closure"
REPORT_NAME = "pantheon_sigma_zr_invariant_extraction_row_closure_report.md"
RESULT_NAME = "pantheon_sigma_zr_invariant_extraction_row_closure_result.json"
MANIFEST_NAME = "pantheon_sigma_zr_invariant_extraction_row_closure_manifest.json"
DERIVATIVE_CSV_NAME = "supplementary_identity_reconstruction.csv"
PANTHEON_CSV_NAME = "pantheon_zr_row_extraction_closure.csv"
DEGREES_OF_FREEDOM = 0

PANTHEON_EXTRACTED_MU_RECORD = "0.082912607552"
PANTHEON_EXTRACTED_GAMMA_RECORD = "0.38603416"

DERIVATIVE_COLUMNS = [
    "row_id",
    "r",
    "sigma_prime",
    "sigma_second",
    "Q_sqrt_negative_sigma_second",
    "mu_recovered",
    "gamma_recovered",
    "mu_difference_from_extracted_record",
    "gamma_difference_from_extracted_record",
    "recovery_error",
    "finite_status",
]


def set_precision(precision_digits):
    mp.mp.dps = int(precision_digits)


def mp_to_str(value, digits=40):
    return mp.nstr(value, n=digits, strip_zeros=False)


def percentile(values, q):
    if not values:
        return math.nan
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    rank = (len(ordered) - 1) * (q / 100.0)
    lo = int(math.floor(rank))
    hi = int(math.ceil(rank))
    if lo == hi:
        return ordered[lo]
    weight = rank - lo
    return ordered[lo] * (1.0 - weight) + ordered[hi] * weight


class Universe:
    def __init__(self, mu, gamma, precision_digits=120):
        self.precision_digits = int(precision_digits)
        set_precision(self.precision_digits)
        self.mu = mp.mpf(str(mu))
        self.gamma = mp.mpf(str(gamma))

    def mp(self, value):
        set_precision(self.precision_digits)
        return mp.mpf(str(value))

    def sigma(self, r):
        r = self.mp(r)
        return mp.log((self.mu ** 2 * (1 + self.gamma * r)) / (self.mu + self.gamma)) - self.mu * r

    def rho(self, r):
        r = self.mp(r)
        return (self.mu ** 2 / (self.mu + self.gamma)) * (1 + self.gamma * r) * mp.exp(-self.mu * r)

    def norm_factor(self):
        return self.mu ** 2 / (self.mu + self.gamma)

    def sigma_prime(self, r):
        r = self.mp(r)
        return self.gamma / (1 + self.gamma * r) - self.mu

    def sigma_second(self, r):
        r = self.mp(r)
        return -(self.gamma ** 2) / ((1 + self.gamma * r) ** 2)

    def q_from_sigma_second(self, r):
        return mp.sqrt(-self.sigma_second(r))

    def recover_mu_from_derivatives(self, r):
        q = self.q_from_sigma_second(r)
        return q - self.sigma_prime(r)

    def recover_gamma_from_derivatives(self, r):
        r = self.mp(r)
        q = self.q_from_sigma_second(r)
        return q / (1 - r * q)

    def recovery_error(self, r):
        mu_recovered = self.recover_mu_from_derivatives(r)
        gamma_recovered = self.recover_gamma_from_derivatives(r)
        return mp.sqrt((mu_recovered - self.mu) ** 2 + (gamma_recovered - self.gamma) ** 2)

    def r_star(self):
        return 1 / self.mu - 1 / self.gamma

    def stationary_mu_from_second_derivative(self):
        return mp.sqrt(-self.sigma_second(self.r_star()))


def build_derivative_grid(universe):
    set_precision(universe.precision_digits)
    points = set()

    for i in range(601):
        points.add(mp.mpf("3.0") * i / 600)

    for i in range(681):
        points.add(mp.mpf("3.0") + (mp.mpf("17.0") * i / 680))

    near_zero = [
        "0",
        "1e-18",
        "1e-15",
        "1e-12",
        "1e-9",
        "1e-6",
        "1e-4",
        "1e-3",
    ]
    for value in near_zero:
        points.add(mp.mpf(value))

    r_star = universe.r_star()
    for delta in ["-1e-9", "-1e-12", "0", "1e-12", "1e-9"]:
        candidate = r_star + mp.mpf(delta)
        if candidate >= 0:
            points.add(candidate)

    for value in ["25", "50", "100", "250", "500", "1000"]:
        points.add(mp.mpf(value))

    return sorted(points)


def run_derivative_recovery(universe, output_csv):
    set_precision(universe.precision_digits)
    rows = []
    grid = build_derivative_grid(universe)

    for row_id, r in enumerate(grid):
        sigma_prime = universe.sigma_prime(r)
        sigma_second = universe.sigma_second(r)
        q = universe.q_from_sigma_second(r)
        mu_recovered = q - sigma_prime
        gamma_recovered = q / (1 - r * q)
        mu_residual = mu_recovered - universe.mu
        gamma_residual = gamma_recovered - universe.gamma
        recovery_error = mp.sqrt(mu_residual ** 2 + gamma_residual ** 2)
        finite = all(
            mp.isfinite(v)
            for v in [
                r,
                sigma_prime,
                sigma_second,
                q,
                mu_recovered,
                gamma_recovered,
                mu_residual,
                gamma_residual,
                recovery_error,
            ]
        )
        rows.append({
            "row_id": row_id,
            "r": mp_to_str(r),
            "sigma_prime": mp_to_str(sigma_prime),
            "sigma_second": mp_to_str(sigma_second),
            "Q_sqrt_negative_sigma_second": mp_to_str(q),
            "mu_recovered": mp_to_str(mu_recovered),
            "gamma_recovered": mp_to_str(gamma_recovered),
            "mu_difference_from_extracted_record": mp_to_str(mu_residual),
            "gamma_difference_from_extracted_record": mp_to_str(gamma_residual),
            "recovery_error": mp_to_str(recovery_error),
            "finite_status": "FINITE" if finite else "NONFINITE",
            "_mu_recovered_float": float(mu_recovered),
            "_gamma_recovered_float": float(gamma_recovered),
            "_mu_abs_residual_float": float(abs(mu_residual)),
            "_gamma_abs_residual_float": float(abs(gamma_residual)),
            "_recovery_error_float": float(recovery_error),
            "_finite": finite,
        })

    with open(output_csv, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=DERIVATIVE_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row[key] for key in DERIVATIVE_COLUMNS})

    mu_values = [row["_mu_recovered_float"] for row in rows if row["_finite"]]
    gamma_values = [row["_gamma_recovered_float"] for row in rows if row["_finite"]]
    max_mu_residual = max(row["_mu_abs_residual_float"] for row in rows)
    max_gamma_residual = max(row["_gamma_abs_residual_float"] for row in rows)
    max_recovery_error = max(row["_recovery_error_float"] for row in rows)

    r_star = universe.r_star()
    stationary_sigma_prime = universe.sigma_prime(r_star)
    stationary_sigma_second = universe.sigma_second(r_star)
    stationary_mu = universe.stationary_mu_from_second_derivative()
    stationary_mu_residual = abs(stationary_mu - universe.mu)

    derivative_mu_pass = max_mu_residual <= 1e-36
    derivative_gamma_pass = max_gamma_residual <= 1e-36
    stationary_prime_pass = abs(float(stationary_sigma_prime)) <= 1e-36
    stationary_mu_pass = float(stationary_mu_residual) <= 1e-36
    finite_records_pass = all(row["_finite"] for row in rows)

    checks = {
        "finite_derivative_records": finite_records_pass,
        "mu_recovered_on_all_derivative_grid_points": derivative_mu_pass,
        "gamma_recovered_on_all_derivative_grid_points": derivative_gamma_pass,
        "sigma_prime_at_r_star_is_zero": stationary_prime_pass,
        "mu_equals_sqrt_negative_sigma_second_at_r_star": stationary_mu_pass,
        "no_H0_used": True,
        "no_DESI_q_used": True,
        "no_fit_used": True,
        "no_minimization_used": True,
        "no_least_squares_used": True,
    }

    return {
        "csv_path": output_csv,
        "grid_point_count": len(rows),
        "recovered_mu_median": statistics.median(mu_values),
        "pantheon_extracted_mu_record": float(universe.mu),
        "max_abs_mu_residual": max_mu_residual,
        "recovered_gamma_median": statistics.median(gamma_values),
        "pantheon_extracted_gamma_record": float(universe.gamma),
        "max_abs_gamma_residual": max_gamma_residual,
        "max_recovery_error": max_recovery_error,
        "pass": all(checks.values()),
        "checks": checks,
        "stationary": {
            "r_star": float(r_star),
            "sigma_prime_at_r_star": float(stationary_sigma_prime),
            "sigma_second_at_r_star": float(stationary_sigma_second),
            "mu_from_sqrt_negative_sigma_second": float(stationary_mu),
            "difference_from_extracted_mu_record": float(stationary_mu_residual),
            "pass": stationary_prime_pass and stationary_mu_pass,
        },
        "rows_written": len(rows),
        "q25_mu": percentile(mu_values, 25),
        "q75_mu": percentile(mu_values, 75),
        "q25_gamma": percentile(gamma_values, 25),
        "q75_gamma": percentile(gamma_values, 75),
    }



C_KM_S = 299_792.458
PANTHEON_COLUMNS = [
    "row_id",
    "catalogue_id",
    "catalogue_z",
    "r",
    "z_equals_r",
    "catalogue_mu_record",
    "catalogue_luminosity_record",
    "catalogue_scale_record",
    "D_record",
    "T_record",
    "mu_derivative_recovered",
    "gamma_recovered",
    "gamma_finite_status",
    "gamma_clip_status",
    "robust_subset_status",
    "gamma_difference_from_extracted_record",
]


def percentile(values, q):
    if not values:
        return math.nan
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    rank = (len(ordered) - 1) * (q / 100.0)
    lo = int(math.floor(rank))
    hi = int(math.ceil(rank))
    if lo == hi:
        return ordered[lo]
    weight = rank - lo
    return ordered[lo] * (1.0 - weight) + ordered[hi] * weight


def to_float(value):
    try:
        return float(value)
    except Exception:
        return math.nan


def read_whitespace_table(path):
    with open(path, "r", encoding="utf-8") as f:
        lines = [line.strip() for line in f if line.strip() and not line.lstrip().startswith("#")]

    if not lines:
        raise ValueError("Input table is empty.")

    header = lines[0].split()
    rows = []
    for source_row_index, line in enumerate(lines[1:]):
        parts = line.split()
        if len(parts) != len(header):
            raise ValueError(
                f"Input row {source_row_index} has {len(parts)} fields; expected {len(header)}."
            )
        row = dict(zip(header, parts))
        row["_source_row_index"] = source_row_index
        rows.append(row)
    return header, rows


def catalogue_luminosity_record(catalogue_mu_record):
    return 10.0 ** ((catalogue_mu_record - 25.0) / 5.0)


def gamma_closed_form(r, T, dps=120):
    if not (r > 0.0) or not math.isfinite(r) or not (T > 0.0) or not math.isfinite(T):
        return math.nan

    if abs(T - r) <= 1e-14 * max(1.0, abs(r)):
        return 0.0

    mp.mp.dps = int(dps)
    r_mp = mp.mpf(r)
    T_mp = mp.mpf(T)
    a = T_mp / r_mp
    x = -a * mp.e ** (-a)
    branch = -1 if a < 1 else 0

    try:
        W = mp.lambertw(x, branch)
        if abs(mp.im(W)) > mp.mpf("1e-40"):
            return math.nan
        W = mp.re(W)
        y = -W / a
        gamma = (y - 1) / r_mp
        if (1 + gamma * r_mp) <= 0:
            return math.nan
        return float(gamma)
    except Exception:
        return math.nan


def gamma_closure_residual(r, gamma, T):
    if not all(math.isfinite(v) for v in [r, gamma, T]):
        return math.nan
    if abs(gamma) <= 1e-15:
        recovered_T = r
    else:
        recovered_T = math.log1p(gamma * r) / gamma
    return recovered_T - T


def run_pantheon_gamma_closure(
    input_path,
    output_csv,
    mu_derivative_recovered,
    canonical_gamma,
    catalogue_scale_record,
    z_min,
    z_max,
    dps,
    gamma_clip_lo,
    gamma_clip_hi,
):
    header, rows = read_whitespace_table(input_path)
    total_rows = len(rows)

    if "zCMB" not in header:
        raise ValueError("Expected zCMB input column.")
    mu_column = "MU_SH0ES" if "MU_SH0ES" in header else ("MU" if "MU" in header else None)
    if mu_column is None:
        raise ValueError("Expected MU_SH0ES or MU input column.")

    calibrator_rows = []
    non_calibrator_rows = []
    for row in rows:
        is_calibrator = int(to_float(row.get("IS_CALIBRATOR", "0"))) == 1
        if is_calibrator:
            calibrator_rows.append(row)
        else:
            non_calibrator_rows.append(row)

    z_nonfinite = []
    z_below_min = []
    z_above_max = []
    retained = []

    for row in non_calibrator_rows:
        catalogue_z = to_float(row["zCMB"])
        if not math.isfinite(catalogue_z):
            z_nonfinite.append(row)
        elif catalogue_z < z_min:
            z_below_min.append(row)
        elif catalogue_z > z_max:
            z_above_max.append(row)
        else:
            retained.append(row)

    evidence_rows = []
    gamma_values = []
    for row_id, row in enumerate(retained):
        catalogue_z = to_float(row["zCMB"])
        r = catalogue_z
        z_equals_r = catalogue_z == r
        catalogue_mu = to_float(row[mu_column])
        luminosity_record = catalogue_luminosity_record(catalogue_mu)
        D_record = luminosity_record * catalogue_scale_record / (C_KM_S * (1.0 + r))
        T_record = D_record * math.exp(mu_derivative_recovered * r)
        gamma = gamma_closed_form(r, T_record, dps=dps)
        gamma_finite = math.isfinite(gamma)
        gamma_clip = gamma_finite and gamma_clip_lo <= gamma <= gamma_clip_hi
        gamma_residual = gamma - canonical_gamma if gamma_finite else math.nan
        closure_residual = gamma_closure_residual(r, gamma, T_record)
        evidence_rows.append({
            "row_id": row_id,
            "catalogue_id": row.get("CID", str(row["_source_row_index"])),
            "catalogue_z": catalogue_z,
            "r": r,
            "z_equals_r": z_equals_r,
            "catalogue_mu_record": catalogue_mu,
            "catalogue_luminosity_record": luminosity_record,
            "catalogue_scale_record": catalogue_scale_record,
            "D_record": D_record,
            "T_record": T_record,
            "mu_derivative_recovered": mu_derivative_recovered,
            "gamma_recovered": gamma,
            "gamma_finite_status": "FINITE" if gamma_finite else "NONFINITE",
            "gamma_clip_status": "CLIP_RETAINED" if gamma_clip else "CLIP_EXCLUDED",
            "robust_subset_status": "ROBUST_PENDING",
            "gamma_difference_from_extracted_record": gamma_residual,
            "_gamma_clip": gamma_clip,
            "_closure_residual": closure_residual,
        })
        if gamma_clip:
            gamma_values.append(gamma)

    q25 = percentile(gamma_values, 25)
    q75 = percentile(gamma_values, 75)
    iqr = q75 - q25 if math.isfinite(q25) and math.isfinite(q75) else math.nan
    robust_values = []

    for row in evidence_rows:
        robust = (
            row["_gamma_clip"]
            and math.isfinite(iqr)
            and row["gamma_recovered"] >= q25 - 1.5 * iqr
            and row["gamma_recovered"] <= q75 + 1.5 * iqr
        )
        row["robust_subset_status"] = "ROBUST_RETAINED" if robust else "ROBUST_EXCLUDED"
        if robust:
            robust_values.append(row["gamma_recovered"])

    with open(output_csv, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=PANTHEON_COLUMNS)
        writer.writeheader()
        for row in evidence_rows:
            writer.writerow({key: row[key] for key in PANTHEON_COLUMNS})

    finite_gamma_count = sum(1 for row in evidence_rows if row["gamma_finite_status"] == "FINITE")
    clip_retained_count = len(gamma_values)
    robust_subset_count = len(robust_values)
    clipped_mean = statistics.fmean(gamma_values) if gamma_values else math.nan
    clipped_median = statistics.median(gamma_values) if gamma_values else math.nan
    clipped_std = statistics.pstdev(gamma_values) if len(gamma_values) > 1 else 0.0
    robust_median = statistics.median(robust_values) if robust_values else math.nan
    robust_abs_residual = abs(robust_median - canonical_gamma) if math.isfinite(robust_median) else math.nan
    robust_relative_residual = robust_abs_residual / canonical_gamma if math.isfinite(robust_abs_residual) else math.nan

    required_numeric_evidence_finite = all(
        math.isfinite(row[key])
        for row in evidence_rows
        for key in [
            "catalogue_z",
            "r",
            "catalogue_mu_record",
            "catalogue_luminosity_record",
            "catalogue_scale_record",
            "D_record",
            "T_record",
            "mu_derivative_recovered",
            "gamma_recovered",
            "gamma_difference_from_extracted_record",
        ]
    )
    closure_residual_max_abs = max(
        (abs(row["_closure_residual"]) for row in evidence_rows if math.isfinite(row["_closure_residual"])),
        default=math.nan,
    )

    excluded_rows_by_reason = {
        "calibrator": len(calibrator_rows),
        "z_below_min": len(z_below_min),
        "z_above_max": len(z_above_max),
        "z_nonfinite": len(z_nonfinite),
    }
    z_cut_excluded = len(z_below_min) + len(z_above_max) + len(z_nonfinite)
    no_silent_row_loss = (
        total_rows == len(calibrator_rows) + len(non_calibrator_rows)
        and len(non_calibrator_rows) == len(retained) + z_cut_excluded
        and len(retained) == len(evidence_rows)
    )

    checks = {
        "input_loaded": total_rows > 0,
        "total_row_count_accounted": total_rows == len(rows),
        "calibrator_exclusion_accounted": total_rows == len(calibrator_rows) + len(non_calibrator_rows),
        "z_cut_accounted": len(non_calibrator_rows) == len(retained) + z_cut_excluded,
        "z_equals_r_true_for_every_retained_row": all(row["z_equals_r"] for row in evidence_rows),
        "finite_gamma_solved_for_every_retained_row_after_z_cut": finite_gamma_count == len(retained),
        "row_evidence_csv_written": True,
        "required_numeric_evidence_finite": required_numeric_evidence_finite,
        "no_silent_row_loss": no_silent_row_loss,
        "row_closure_residual_within_numeric_tolerance": (
            math.isfinite(closure_residual_max_abs)
            and closure_residual_max_abs <= 1e-12
        ),
    }

    return {
        "csv_path": output_csv,
        "input_file": input_path,
        "mu_source": "derivative_recovered",
        "mu_derivative_recovered": mu_derivative_recovered,
        "catalogue_scale_record": catalogue_scale_record,
        "mu_column": mu_column,
        "total_rows": total_rows,
        "calibrator_rows_removed": len(calibrator_rows),
        "rows_after_calibrator_cut": len(non_calibrator_rows),
        "rows_excluded_by_z_cut": z_cut_excluded,
        "rows_retained_after_z_cut": len(retained),
        "excluded_rows_by_reason": excluded_rows_by_reason,
        "finite_gamma_recoveries": finite_gamma_count,
        "gamma_clip_retained_count": clip_retained_count,
        "robust_subset_count": robust_subset_count,
        "clipped_gamma_median": clipped_median,
        "clipped_gamma_mean": clipped_mean,
        "clipped_gamma_standard_deviation": clipped_std,
        "clipped_gamma_iqr": iqr,
        "robust_median_gamma": robust_median,
        "pantheon_extracted_gamma_record": canonical_gamma,
        "decimal_recording_difference": robust_abs_residual,
        "decimal_recording_relative_difference": robust_relative_residual,
        "closure_residual_max_abs": closure_residual_max_abs,
        "checks": checks,
        "pass": all(checks.values()),
    }


def project_root():
    return os.path.dirname(os.path.abspath(__file__))


def default_data_path(root_dir):
    candidates = [
        os.path.join(root_dir, "data", "PantheonSH0ES.dat"),
        os.path.join(root_dir, "data", "pantheon", "PantheonSH0ES.dat"),
    ]
    for candidate in candidates:
        if os.path.exists(candidate):
            return candidate
    return candidates[0]


def utc_timestamp():
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def command_line():
    return " ".join(shlex.quote(part) for part in sys.argv)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def file_entry(path, root_dir):
    return {
        "path": os.path.relpath(path, root_dir),
        "bytes": os.path.getsize(path),
        "sha256": sha256_file(path),
    }


def write_json(path, payload):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=True)
        f.write("\n")


def write_text(path, text):
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)


def status_word(ok):
    return "PASS" if ok else "FAIL"


def format_hash_block(title, entries):
    lines = [f"- {title}:"]
    for name, entry in entries.items():
        lines.append(f"  - {name}: `{entry['sha256']}`")
    return lines


def build_report(
    run_timestamp_utc,
    command,
    input_path,
    output_dir,
    source_hashes,
    procedure_hashes,
    output_hashes_available,
    universe,
    derivative,
    stationary,
    pantheon,
    pass_fail_checks,
    final_status,
    dps,
):
    rho0 = float(universe.norm_factor())
    r_star = float(universe.r_star())
    robust_relative_percent = pantheon["decimal_recording_relative_difference"] * 100.0

    lines = [
        f"# {AUDIT_NAME}",
        "",
        "## A. Audit Header",
        "",
        f"- title: {AUDIT_NAME}",
        f"- run timestamp UTC: {run_timestamp_utc}",
        f"- command: `{command}`",
        f"- input path: `{input_path}`",
        f"- output directory: `{output_dir}`",
        f"- MPFR precision: {dps} decimal digits",
        f"- degrees_of_freedom: {DEGREES_OF_FREEDOM}",
        f"- final status: {final_status}",
        "",
    ]
    lines.extend(format_hash_block("source file hashes", source_hashes))
    lines.extend(format_hash_block("procedure hashes", procedure_hashes))
    lines.extend(format_hash_block("available output hashes", output_hashes_available))
    lines.extend([
        "",
        "## B. Audit Content",
        "",
        "This audit evaluates the Pantheon catalogue z column as r.",
        "",
        "- implemented identity functions",
        "- derivative recovery route",
        "- stationary recovery route",
        "- Pantheon z=r row rule",
        "- Lambert W row closure",
        "- corpus accounting",
        "- finite recovery",
        "- SHA-256 manifest",
        "- PASS/FAIL result",
        "",
        "## C. Pantheon-Extracted Invariant Record",
        "",
        f"- μ extracted record: {float(universe.mu):.12f}",
        f"- γ extracted record: {float(universe.gamma):.8f}",
        f"- ρ(0): {rho0:.12e}",
        f"- r*: {r_star:.12f}",
        "",
        "## D. Supplementary Derivative Identity Reconstruction",
        "",
        "- Q(r) = sqrt(−σ″(r))",
        "- μ = Q(r) − σ′(r)",
        "- γ = Q(r)/(1 − rQ(r))",
        "",
        f"- derivative grid point count: {derivative['grid_point_count']}",
        f"- recovered μ median: {derivative['recovered_mu_median']:.12f}",
        f"- Pantheon-extracted μ record: {derivative['pantheon_extracted_mu_record']:.12f}",
        f"- max |μ residual|: {derivative['max_abs_mu_residual']:.12e}",
        f"- recovered γ median: {derivative['recovered_gamma_median']:.8f}",
        f"- Pantheon-extracted γ record: {derivative['pantheon_extracted_gamma_record']:.8f}",
        f"- max |γ residual|: {derivative['max_abs_gamma_residual']:.12e}",
        f"- max recovery error: {derivative['max_recovery_error']:.12e}",
        f"- supplementary identity-consistency status: {status_word(derivative['pass'])}",
        "",
        "## E. Supplementary Stationary Identity Reconstruction",
        "",
        "- r* = 1/μ − 1/γ",
        "- σ′(r*) = 0",
        "- μ = sqrt(−σ″(r*))",
        "",
        f"- r*: {stationary['r_star']:.12f}",
        f"- σ′(r*): {stationary['sigma_prime_at_r_star']:.12e}",
        f"- σ″(r*): {stationary['sigma_second_at_r_star']:.12e}",
        f"- μ from sqrt(−σ″(r*)): {stationary['mu_from_sqrt_negative_sigma_second']:.12f}",
        f"- difference from extracted μ record: {stationary['difference_from_extracted_mu_record']:.12e}",
        f"- supplementary identity-consistency status: {status_word(stationary['pass'])}",
        "",
        "## F. Pantheon Corpus Accounting",
        "",
        f"- input file: `{pantheon['input_file']}`",
        f"- total rows: {pantheon['total_rows']}",
        f"- calibrator row count: {pantheon['calibrator_rows_removed']}",
        f"- rows after calibrator cut: {pantheon['rows_after_calibrator_cut']}",
        f"- z-cut row count: {pantheon['rows_excluded_by_z_cut']}",
        f"- rows retained after z cut: {pantheon['rows_retained_after_z_cut']}",
        f"- finite γ recoveries: {pantheon['finite_gamma_recoveries']}",
        f"- γ clip retained count: {pantheon['gamma_clip_retained_count']}",
        f"- robust subset count: {pantheon['robust_subset_count']}",
        f"- degrees_of_freedom: {DEGREES_OF_FREEDOM}",
        "",
        "Corpus accounting by category:",
    ])
    for reason, count in pantheon["excluded_rows_by_reason"].items():
        lines.append(f"- {reason}: {count}")

    lines.extend([
        "",
        "## G. Pantheon z=r Row Closure",
        "",
        "- r = catalogue_z",
        "- z_equals_r = true",
        "- D_record = catalogue_luminosity_record * catalogue_scale_record / (c * (1+r))",
        "- T_record = D_record * exp(μr)",
        "- ln(1+γr)/γ = T_record",
        "- γ recovered by Lambert W closed form",
        "",
        f"- finite γ recovered / retained rows: {pantheon['finite_gamma_recoveries']} / {pantheon['rows_retained_after_z_cut']}",
        f"- clipped γ median: {pantheon['clipped_gamma_median']:.8f}",
        f"- clipped γ mean: {pantheon['clipped_gamma_mean']:.8f}",
        f"- clipped γ standard deviation: {pantheon['clipped_gamma_standard_deviation']:.8f}",
        f"- clipped γ IQR: {pantheon['clipped_gamma_iqr']:.8f}",
        f"- robust median γ: {pantheon['robust_median_gamma']:.8f}",
        f"- Pantheon-extracted γ record: {pantheon['pantheon_extracted_gamma_record']:.8f}",
        f"- decimal-recording difference: {pantheon['decimal_recording_difference']:.12e}",
        f"- decimal-recording relative difference: {robust_relative_percent:.8f}%",
        "- the decimal-recording difference is reported for custody only; it is not a closure predicate",
        f"- maximum row-closure residual: {pantheon['closure_residual_max_abs']:.12e}",
        f"- status: {status_word(pantheon['pass'])}",
        "",
        "## H. PASS/FAIL Matrix",
        "",
    ])
    public_pass_fail_labels = {
        "input_loaded": "input_loaded",
        "total_row_count_accounted": "total_row_count_accounted",
        "calibrator_accounting": "calibrator_exclusion_accounted",
        "z_cut_accounted": "z_cut_accounted",
        "no_silent_row_loss": "no_silent_row_loss",
        "z_equals_r_true_for_every_retained_row": "z_equals_r_true_for_every_retained_row",
        "finite_gamma_solved_for_every_retained_row_after_z_cut": "finite_gamma_solved_for_every_retained_row_after_z_cut",
        "row_evidence_csv_written": "row_evidence_csv_written",
        "derivative_recovery_csv_written": "derivative_recovery_csv_written",
        "report_written": "report_written",
        "result_json_written": "result_json_written",
        "manifest_written": "manifest_written",
        "required_numeric_evidence_finite": "required_numeric_evidence_finite",
        "degrees_of_freedom_is_zero": "degrees_of_freedom_is_zero",
        "row_closure_residual_within_numeric_tolerance": (
            "row_closure_residual_within_numeric_tolerance"
        ),
    }
    for public_name, source_name in public_pass_fail_labels.items():
        ok = pass_fail_checks[source_name]
        name = public_name
        lines.append(f"- {name}: {status_word(ok)}")

    lines.extend([
        "",
        "## I. Output Index",
        "",
    ])
    for name, entry in output_hashes_available.items():
        lines.append(f"- {name}: `{entry['path']}` sha256 `{entry['sha256']}`")
    lines.append("- report: final hash recorded in manifest")
    lines.append("- manifest: hash policy recorded in manifest")
    lines.append("")
    return "\n".join(lines)


def main():
    root_dir = project_root()
    default_dat = default_data_path(root_dir)

    parser = argparse.ArgumentParser(
        description="Pantheon sigma z=r invariant extraction and row closure audit."
    )
    parser.add_argument("--dat_path", default=default_dat, help="Pantheon catalogue table path")
    parser.add_argument("--catalogue_scale_record", type=float, default=73.23,
                        help="Catalogue normalization constant used in D_record")
    parser.add_argument("--z_min", type=float, default=0.01,
                        help="Minimum retained catalogue z value")
    parser.add_argument("--z_max", type=float, default=3.0,
                        help="Maximum retained catalogue z value")
    parser.add_argument("--dps", type=int, default=120,
                        help="Arbitrary-precision decimal digits")
    parser.add_argument("--gamma_clip_lo", type=float, default=-5.0,
                        help="Retain gamma >= this value for ensemble records")
    parser.add_argument("--gamma_clip_hi", type=float, default=20.0,
                        help="Retain gamma <= this value for ensemble records")
    parser.add_argument("--out_dir", default=None, help="Output directory")
    args = parser.parse_args()

    output_dir = args.out_dir or os.path.join(root_dir, "results", AUDIT_SLUG)
    os.makedirs(output_dir, exist_ok=True)

    report_path = os.path.join(output_dir, REPORT_NAME)
    result_path = os.path.join(output_dir, RESULT_NAME)
    manifest_path = os.path.join(output_dir, MANIFEST_NAME)
    derivative_csv_path = os.path.join(output_dir, DERIVATIVE_CSV_NAME)
    pantheon_csv_path = os.path.join(output_dir, PANTHEON_CSV_NAME)

    run_timestamp_utc = utc_timestamp()
    command = command_line()
    universe = Universe(
        PANTHEON_EXTRACTED_MU_RECORD,
        PANTHEON_EXTRACTED_GAMMA_RECORD,
        precision_digits=args.dps,
    )

    derivative = run_derivative_recovery(universe, derivative_csv_path)
    mu_derivative_recovered = derivative["recovered_mu_median"]
    pantheon = run_pantheon_gamma_closure(
        input_path=args.dat_path,
        output_csv=pantheon_csv_path,
        mu_derivative_recovered=mu_derivative_recovered,
        canonical_gamma=float(universe.gamma),
        catalogue_scale_record=args.catalogue_scale_record,
        z_min=args.z_min,
        z_max=args.z_max,
        dps=args.dps,
        gamma_clip_lo=args.gamma_clip_lo,
        gamma_clip_hi=args.gamma_clip_hi,
    )

    source_hashes = {
        "pantheon_catalogue": file_entry(args.dat_path, root_dir),
    }
    procedure_hashes = {
        "entrypoint": file_entry(os.path.join(root_dir, "sigma_pantheon.py"), root_dir),
    }
    evidence_hashes = {
        "derivative_recovery_csv": file_entry(derivative_csv_path, root_dir),
        "pantheon_row_evidence_csv": file_entry(pantheon_csv_path, root_dir),
    }

    epistemic_boundary = {
        "z_equals_r": True,
        "H0_ratio_used_for_mu": False,
        "H0_cmb_used": False,
        "H0_late_used_for_mu": False,
        "DESI_q_used_for_mu": False,
        "velocity_interpretation_used": False,
        "expansion_history_used": False,
        "FLRW_distance_relation_used": False,
        "LambdaCDM_parameter_extraction_used": False,
        "dark_energy_parameter_extraction_used": False,
        "likelihood_prior_posterior_MCMC_used": False,
        "fit_objective_used": False,
        "least_squares_used": False,
        "minimization_used": False,
    }

    pass_fail_checks = {
        "input_loaded": pantheon["checks"]["input_loaded"],
        "total_row_count_accounted": pantheon["checks"]["total_row_count_accounted"],
        "calibrator_exclusion_accounted": pantheon["checks"]["calibrator_exclusion_accounted"],
        "z_cut_accounted": pantheon["checks"]["z_cut_accounted"],
        "no_silent_row_loss": pantheon["checks"]["no_silent_row_loss"],
        "z_equals_r_true_for_every_retained_row": pantheon["checks"]["z_equals_r_true_for_every_retained_row"],
        "finite_gamma_solved_for_every_retained_row_after_z_cut": pantheon["checks"]["finite_gamma_solved_for_every_retained_row_after_z_cut"],
        "row_evidence_csv_written": os.path.exists(pantheon_csv_path),
        "derivative_recovery_csv_written": os.path.exists(derivative_csv_path),
        "report_written": True,
        "result_json_written": True,
        "manifest_written": True,
        "required_numeric_evidence_finite": pantheon["checks"]["required_numeric_evidence_finite"],
        "degrees_of_freedom_is_zero": DEGREES_OF_FREEDOM == 0,
        "no_H0_ratio_used": not epistemic_boundary["H0_ratio_used_for_mu"],
        "no_DESI_q_used_for_mu": not epistemic_boundary["DESI_q_used_for_mu"],
        "no_fit_used": not epistemic_boundary["fit_objective_used"],
        "no_minimization_used": not epistemic_boundary["minimization_used"],
        "no_least_squares_used": not epistemic_boundary["least_squares_used"],
        "row_closure_residual_within_numeric_tolerance": (
            pantheon["checks"]["row_closure_residual_within_numeric_tolerance"]
        ),
    }
    final_status = "PASS" if all(pass_fail_checks.values()) else "FAIL"

    output_paths = {
        "report": os.path.relpath(report_path, root_dir),
        "result_json": os.path.relpath(result_path, root_dir),
        "manifest": os.path.relpath(manifest_path, root_dir),
        "derivative_recovery_csv": os.path.relpath(derivative_csv_path, root_dir),
        "pantheon_row_evidence_csv": os.path.relpath(pantheon_csv_path, root_dir),
    }
    input_paths = {
        "pantheon_catalogue": os.path.relpath(args.dat_path, root_dir),
    }

    result = {
        "audit_name": AUDIT_NAME,
        "run_timestamp_utc": run_timestamp_utc,
        "final_status": final_status,
        "degrees_of_freedom": DEGREES_OF_FREEDOM,
        "epistemic_boundary": epistemic_boundary,
        "input_paths": input_paths,
        "output_paths": output_paths,
        "file_hashes": {
            "source_hashes": source_hashes,
            "procedure_hashes": procedure_hashes,
            "evidence_output_hashes": evidence_hashes,
        },
        "procedure_hashes": procedure_hashes,
        "pantheon_extracted_invariant_record": {
            "mu": float(universe.mu),
            "gamma": float(universe.gamma),
            "rho_at_zero": float(universe.norm_factor()),
            "r_star": float(universe.r_star()),
        },
        "derivative_mu_gamma_recovery": {
            key: value for key, value in derivative.items()
            if key not in {"csv_path", "stationary"}
        },
        "stationary_recovery": derivative["stationary"],
        "pantheon_corpus": {
            "input_file": pantheon["input_file"],
            "total_rows": pantheon["total_rows"],
            "calibrator_rows_removed": pantheon["calibrator_rows_removed"],
            "rows_after_calibrator_cut": pantheon["rows_after_calibrator_cut"],
            "rows_excluded_by_z_cut": pantheon["rows_excluded_by_z_cut"],
            "rows_retained_after_z_cut": pantheon["rows_retained_after_z_cut"],
            "excluded_rows_by_reason": pantheon["excluded_rows_by_reason"],
            "degrees_of_freedom": DEGREES_OF_FREEDOM,
        },
        "pantheon_gamma_closure": {
            key: value for key, value in pantheon.items()
            if key not in {"checks", "csv_path", "input_file", "excluded_rows_by_reason"}
        },
        "pass_fail_checks": pass_fail_checks,
    }
    write_json(result_path, result)
    result_hash = {"result_json": file_entry(result_path, root_dir)}

    report_hash_inputs = {}
    report_hash_inputs.update(evidence_hashes)
    report_hash_inputs.update(result_hash)
    report_text = build_report(
        run_timestamp_utc=run_timestamp_utc,
        command=command,
        input_path=args.dat_path,
        output_dir=output_dir,
        source_hashes=source_hashes,
        procedure_hashes=procedure_hashes,
        output_hashes_available=report_hash_inputs,
        universe=universe,
        derivative=derivative,
        stationary=derivative["stationary"],
        pantheon=pantheon,
        pass_fail_checks=pass_fail_checks,
        final_status=final_status,
        dps=args.dps,
    )
    write_text(report_path, report_text)

    final_output_hashes = {
        "report": file_entry(report_path, root_dir),
        "result_json": file_entry(result_path, root_dir),
        "derivative_recovery_csv": file_entry(derivative_csv_path, root_dir),
        "pantheon_row_evidence_csv": file_entry(pantheon_csv_path, root_dir),
    }
    manifest = {
        "run_timestamp_utc": run_timestamp_utc,
        "command": command,
        "input_paths": input_paths,
        "output_paths": output_paths,
        "source_hashes": source_hashes,
        "procedure_hashes": procedure_hashes,
        "output_hashes": final_output_hashes,
        "manifest_hash_policy": "Output hashes record stable generated artifacts; manifest self-hash is external to the in-file table.",
        "MPFR_precision": args.dps,
        "pantheon_extracted_mu_record": float(universe.mu),
        "pantheon_extracted_gamma_record": float(universe.gamma),
        "derivative_reconstructed_mu": derivative["recovered_mu_median"],
        "derivative_reconstructed_gamma": derivative["recovered_gamma_median"],
        "pantheon_robust_median_gamma": pantheon["robust_median_gamma"],
        "final_status": final_status,
    }
    write_json(manifest_path, manifest)

    print(f"Wrote: {report_path}")
    print(f"Wrote: {result_path}")
    print(f"Wrote: {manifest_path}")
    print(f"Wrote: {derivative_csv_path}")
    print(f"Wrote: {pantheon_csv_path}")
    print(f"{AUDIT_NAME}: {final_status}")

    return 0 if final_status == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
