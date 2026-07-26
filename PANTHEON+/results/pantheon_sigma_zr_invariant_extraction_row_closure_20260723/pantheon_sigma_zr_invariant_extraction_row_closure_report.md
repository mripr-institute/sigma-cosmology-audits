# Pantheon σ z=r Invariant Extraction and Row Closure Audit

## A. Audit Header

- title: Pantheon σ z=r Invariant Extraction and Row Closure Audit
- run timestamp UTC: 2026-07-23T17:05:45Z
- command: `sigma_pantheon.py --out_dir audit_records/pantheon_sigma_zr_invariant_extraction_row_closure_20260723`
- input path: `/home/alex-albert/Documents/GIPR/PANTHEON+/data/PantheonSH0ES.dat`
- output directory: `audit_records/pantheon_sigma_zr_invariant_extraction_row_closure_20260723`
- MPFR precision: 120 decimal digits
- degrees_of_freedom: 0
- final status: PASS

- source file hashes:
  - pantheon_catalogue: `1cb0fc379ef066afdc2ffd1857681cc478024570d8a3eba284fb645775198cf8`
- procedure hashes:
  - entrypoint: `2ac5cede8bf407d19d78bcd1ad32aa98fe22b572bbfdb34abefd3b35d66691a2`
- available output hashes:
  - derivative_recovery_csv: `bfad106189002a3b744895922d2c4d73c225ad93fd4b8c5a76b6c35da149c87e`
  - pantheon_row_evidence_csv: `c76688840f049527a3c0559aab34d8da901ec6440cd90a1ea2d0e5a079fee217`
  - result_json: `d40ce696fda3f6ba5607456fe4bce12953b3fabb1fdaad326f852db05a55d2be`

## B. Audit Content

This audit evaluates the Pantheon catalogue z column as r.

- implemented identity functions
- derivative recovery route
- stationary recovery route
- Pantheon z=r row rule
- Lambert W row closure
- corpus accounting
- finite recovery
- SHA-256 manifest
- PASS/FAIL result

## C. Pantheon-Extracted Invariant Record

- μ extracted record: 0.082912607552
- γ extracted record: 0.38603416
- ρ(0): 1.465944743997e-02
- r*: 9.470447610694

## D. Supplementary Derivative Identity Reconstruction

- Q(r) = sqrt(−σ″(r))
- μ = Q(r) − σ′(r)
- γ = Q(r)/(1 − rQ(r))

- derivative grid point count: 1299
- recovered μ median: 0.082912607552
- Pantheon-extracted μ record: 0.082912607552
- max |μ residual|: 6.050924866952e-122
- recovered γ median: 0.38603416
- Pantheon-extracted γ record: 0.38603416
- max |γ residual|: 4.937554691433e-120
- max recovery error: 4.937554691433e-120
- supplementary identity-consistency status: PASS

## E. Supplementary Stationary Identity Reconstruction

- r* = 1/μ − 1/γ
- σ′(r*) = 0
- μ = sqrt(−σ″(r*))

- r*: 9.470447610694
- σ′(r*): -2.420369946781e-122
- σ″(r*): -6.874500491072e-03
- μ from sqrt(−σ″(r*)): 0.082912607552
- difference from extracted μ record: 2.420369946781e-122
- supplementary identity-consistency status: PASS

## F. Pantheon Corpus Accounting

- input file: `/home/alex-albert/Documents/GIPR/PANTHEON+/data/PantheonSH0ES.dat`
- total rows: 1701
- calibrator row count: 77
- rows after calibrator cut: 1624
- z-cut row count: 46
- rows retained after z cut: 1578
- finite γ recoveries: 1578
- γ clip retained count: 1463
- robust subset count: 1134
- degrees_of_freedom: 0

Corpus accounting by category:
- calibrator: 77
- z_below_min: 46
- z_above_max: 0
- z_nonfinite: 0

## G. Pantheon z=r Row Closure

- r = catalogue_z
- z_equals_r = true
- D_record = catalogue_luminosity_record * catalogue_scale_record / (c * (1+r))
- T_record = D_record * exp(μr)
- ln(1+γr)/γ = T_record
- γ recovered by Lambert W closed form

- finite γ recovered / retained rows: 1578 / 1578
- clipped γ median: 0.38608764
- clipped γ mean: 0.65622081
- clipped γ standard deviation: 2.60806695
- clipped γ IQR: 0.88805588
- robust median γ: 0.38603416
- Pantheon-extracted γ record: 0.38603416
- decimal-recording difference: 2.242115049178e-09
- decimal-recording relative difference: 0.00000058%
- the decimal-recording difference is reported for custody only; it is not a closure predicate
- maximum row-closure residual: 2.220446049250e-16
- status: PASS

## H. PASS/FAIL Matrix

- input_loaded: PASS
- total_row_count_accounted: PASS
- calibrator_accounting: PASS
- z_cut_accounted: PASS
- no_silent_row_loss: PASS
- z_equals_r_true_for_every_retained_row: PASS
- finite_gamma_solved_for_every_retained_row_after_z_cut: PASS
- row_evidence_csv_written: PASS
- derivative_recovery_csv_written: PASS
- report_written: PASS
- result_json_written: PASS
- manifest_written: PASS
- required_numeric_evidence_finite: PASS
- degrees_of_freedom_is_zero: PASS
- row_closure_residual_within_numeric_tolerance: PASS

## I. Output Index

- derivative_recovery_csv: `audit_records/pantheon_sigma_zr_invariant_extraction_row_closure_20260723/supplementary_identity_reconstruction.csv` sha256 `bfad106189002a3b744895922d2c4d73c225ad93fd4b8c5a76b6c35da149c87e`
- pantheon_row_evidence_csv: `audit_records/pantheon_sigma_zr_invariant_extraction_row_closure_20260723/pantheon_zr_row_extraction_closure.csv` sha256 `c76688840f049527a3c0559aab34d8da901ec6440cd90a1ea2d0e5a079fee217`
- result_json: `audit_records/pantheon_sigma_zr_invariant_extraction_row_closure_20260723/pantheon_sigma_zr_invariant_extraction_row_closure_result.json` sha256 `d40ce696fda3f6ba5607456fe4bce12953b3fabb1fdaad326f852db05a55d2be`
- report: final hash recorded in manifest
- manifest: hash policy recorded in manifest
