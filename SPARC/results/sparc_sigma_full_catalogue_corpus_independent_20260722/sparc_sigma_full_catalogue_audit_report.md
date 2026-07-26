# SPARC σ Full-Catalogue Audit Record

## 1. Audit Header

- audit name: SPARC σ Full-Catalogue Audit Record
- run timestamp UTC: 2026-07-22T11:17:59+00:00
- source file hash, SPARC_Lelli2016c.mrt: 5aa0501f6b0d881fa579030e315e7b5b6ef561a5bd3a07472f9929c7e5728243
- source file hash, rotation-curve aggregate: 3013b06ded3b40d0eecc4daa01555bbe9f108ce95a94604de92f6cfb7e12238a
- SPARC data source path: /home/alex-albert/Documents/GIPR/SPARC/data/sparc
- output directory: /home/alex-albert/Documents/GIPR/SPARC/audit_records/sparc_sigma_full_catalogue_corpus_independent_20260722
- μ: 0.082912607552
- γ: 0.38603416
- degrees_of_freedom: 0
- final status: PASS

## 2. Fixed Identity Block

The fixed audit identities are:

- ρ(0) = μ² / (μ + γ)
- σ′(r) = γ / (1 + γr) - μ
- σ″(r) = -γ² / (1 + γr)²

| quantity | value |
| --- | ---: |
| ρ(0) | 0.01465944744 |
| γ² | 0.149022372687 |
| μ² | 0.00687450049107 |
| degrees_of_freedom | 0 |

## 3. Full-Catalogue Rotation Closure Block

- total galaxies: 175
- galaxies with valid Vflat: 135
- total rotation curve points: 3391
- V(Rlast) / Vflat mean: 1.00332
- V(Rlast) / Vflat median: 1.00506
- V(Rlast) / Vflat std: 0.0238585
- |V(Rlast)/Vflat - 1| < 10%: 135/135 (100.00%)
- |V(Rlast)/Vflat - 1| < 20%: 135/135 (100.00%)
- degrees_of_freedom: 0

## 4. Per-Galaxy Outer-Flatness Table

- file: /home/alex-albert/Documents/GIPR/SPARC/audit_records/sparc_sigma_full_catalogue_corpus_independent_20260722/outer_flatness_by_galaxy.csv
- row count: 175
- expected row count: 175
- no silent galaxy loss: True
- required numeric audit values finite: True
- status counts: {"INSUFFICIENT_OUTER_POINTS": 50, "INVALID_VFLAT": 14, "OK": 111}
- outer_slope_mean: 0.0286883
- outer_slope_median: 0.0104273
- outer_slope_std: 0.153755
- outer_slope_abs_median: 0.0951692
- degrees_of_freedom: 0

## 5. σ Derivative Profile Table

- file: /home/alex-albert/Documents/GIPR/SPARC/audit_records/sparc_sigma_full_catalogue_corpus_independent_20260722/sigma_derivative_profile_sparc.csv
- row count: 3391
- expected rotation-point row count: 3391
- no silent galaxy loss among galaxies with rotation curves: True
- required numeric audit values finite: True
- status counts: {"INVALID_VFLAT": 377, "OK": 3014}
- closure_reference_value is the fixed identity value 1.0 where Vflat is valid.
- No new fitted quantity was introduced; closure_residual is V_obs / Vflat - 1 where available.
- sigma min/max: -4.44969 / -3.46975
- rho min/max: 0.0116822 / 0.0311248
- sigma_prime min/max: -0.0563101 / 0.300736
- sigma_second min/max: -0.147186 / -0.000707693
- degrees_of_freedom: 0

## 6. RAR Consistency Block

- galaxies with RAR data: 175
- individual RAR points: 3389
- g_obs min/max: 8.02107e-13 / 1.92144e-08
- g_bar min/max: 2.08672e-13 / 7.34002e-09
- g_bar dynamic range: 35174.9
- RMS scatter dex: 0.195995
- median residual dex: -0.01422
- mean residual dex: -0.0336379
- Pearson r: 0.936046
- Spearman rho: 0.929916
- quality-stratified RMS, Q = 1: 0.159077
- quality-stratified RMS, Q ≤ 2: 0.181051
- quality-stratified RMS, all: 0.195995
- degrees_of_freedom: 0

## 7. PASS/FAIL Logic

- existing_SPARC_identity_audit_PASS: True
- existing_RAR_audit_PASS: True
- outer_flatness_by_galaxy_csv_created: True
- sigma_derivative_profile_sparc_csv_created: True
- outer_flatness_row_count_matches_expected: True
- sigma_profile_row_count_matches_expected: True
- outer_flatness_no_silent_galaxy_loss: True
- sigma_profile_no_silent_galaxy_loss: True
- outer_flatness_required_numeric_values_finite: True
- sigma_profile_required_numeric_values_finite: True
- corpus_counts_consistent: True
- degrees_of_freedom_zero: True
- degrees_of_freedom_in_all_main_report_blocks: True

- final status: PASS

## 8. Output Files

- audit_report: /home/alex-albert/Documents/GIPR/SPARC/audit_records/sparc_sigma_full_catalogue_corpus_independent_20260722/sparc_sigma_full_catalogue_audit_report.md
- audit_result: /home/alex-albert/Documents/GIPR/SPARC/audit_records/sparc_sigma_full_catalogue_corpus_independent_20260722/sparc_sigma_full_catalogue_audit_result.json
- outer_flatness: /home/alex-albert/Documents/GIPR/SPARC/audit_records/sparc_sigma_full_catalogue_corpus_independent_20260722/outer_flatness_by_galaxy.csv
- sigma_profile: /home/alex-albert/Documents/GIPR/SPARC/audit_records/sparc_sigma_full_catalogue_corpus_independent_20260722/sigma_derivative_profile_sparc.csv
- audit_manifest: /home/alex-albert/Documents/GIPR/SPARC/audit_records/sparc_sigma_full_catalogue_corpus_independent_20260722/sparc_sigma_full_catalogue_audit_manifest.json
- figures: /home/alex-albert/Documents/GIPR/SPARC/audit_records/sparc_sigma_full_catalogue_corpus_independent_20260722/figures

Figure audit statuses:

- SPARC identity audit: PASS
- RAR audit: PASS

Computed figure evidence predicates:

- sparc_identity_audit.figures_created: True
- sparc_identity_audit.valid_vflat_records_present: True
- sparc_identity_audit.all_valid_vflat_within_10_percent: True
- sparc_identity_audit.all_valid_vflat_within_20_percent: True
- sparc_identity_audit.rotation_summary_finite: True
- rar_audit.figures_created: True
- rar_audit.galaxies_present: True
- rar_audit.rar_points_present: True
- rar_audit.rar_statistics_finite: True
- rar_audit.positive_rank_and_linear_correlations: True

Computed fitted-degree-of-freedom contract:

- corpus_counts: 0
- fixed_identity: 0
- rotation_closure: 0
- rar_consistency: 0
- all fields present: True
- all fields zero: True

## Audit Invariants

- μ and γ are fixed.
- degrees_of_freedom: 0
- source SPARC data unchanged
- all catalogue galaxies retained with explicit status
- no fitted quantity introduced

## SPARC Closure Statement

The full SPARC catalogue was evaluated with fixed μ and γ.
The audit record verifies:
- positive finite centre: ρ(0)
- first derivative: σ′(r)
- second derivative: σ″(r)
- full rotation-curve corpus accounting
- outer-flatness evidence by galaxy
- RAR consistency
- degrees_of_freedom: 0
Final status: PASS

## 10. Final Check

SPARC σ full-catalogue audit record: PASS
