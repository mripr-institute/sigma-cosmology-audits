# DESI σ Squared-Radial Closure Audit

Official reproducibility repository of the Mathematical Research Institute of Physical Reality (MRIPR) for the DESI σ closure and full-corpus audits.

## Identity

The audit evaluates the fixed identity

\[
\sigma(r)=\ln\!\left[\frac{\mu^2(1+\gamma r)}{\mu+\gamma}\right]-\mu r,
\qquad
\rho(r)=e^{\sigma(r)},
\]

with

\[
\mu=0.082912607552,
\qquad
\gamma=0.38603416,
\qquad
r_*=\frac1\mu-\frac1\gamma.
\]

The DESI spatial audit evaluates the squared-radial closure predicate recorded by the procedure:

\[
\|Z\hat n-C\|^2=r_*^2.
\]

The identity has zero DESI-fitted degrees of freedom. The geometric offset scan does not change \(\mu\), \(\gamma\), or σ and is repeated under the null procedures.

## Authoritative results

- Closure and optimized null audit: `results/sigma_desi_lss_closure_squared_radial_fastnull_20260724/`
- Full 13,097,304-object audit: `results/sigma_desi_full_rho_squared_radial_final/`
- Pooled invariant-support coverage: `0.9708459316265691`
- Full-corpus fast/direct arithmetic maximum difference: `2.842170943040401e-14`
- Optimized null ensemble: 200 trials in each of three families; observed φ-excess exceeded all 600 null realizations.

The Markdown reports are the human-readable entry points. JSON manifests and CSV/NPZ products are the machine-readable evidence.

## Repository layout

```text
.
├── README.md
├── AUTHORS.md
├── CITATION.cff
├── LICENSE
├── REPRODUCIBILITY.md
├── RESULTS.md
├── requirements.txt
├── download_data.py
├── data_sources.json
├── analysis/
│   ├── sigma_desi_lss_closure_squared_radial_audit.py
│   └── sigma_desi_full_rho_squared_radial_audit.py
├── Makefile
├── data/desi/                 # public catalogues; never modified by the audit
└── results/                   # immutable timestamped and final audit records
```

## Reproduce

Create an isolated Python environment and install the pinned minimum dependencies:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

Download missing public catalogues or verify existing files:

```bash
.venv/bin/python download_data.py
.venv/bin/python download_data.py --verify-only
```

Run the closure audit with the complete optimized null ensemble:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
.venv/bin/python analysis/sigma_desi_lss_closure_squared_radial_audit.py \
  --data-root data/desi \
  --output-dir results/sigma_desi_lss_closure_squared_radial_release \
  --seed 12345 --null-trials 200 --null-sub 500000 --null-workers 12
```

Run the full-corpus audit using the completed closure manifest:

```bash
.venv/bin/python analysis/sigma_desi_full_rho_squared_radial_audit.py \
  --data-root data/desi \
  --completed-closure results/sigma_desi_lss_closure_squared_radial_release/sigma_desi_lss_closure_squared_radial_manifest.json \
  --output-dir results/sigma_desi_full_rho_squared_radial_release \
  --mode full --seed 12345 --q-route compare --log-level normal
```

## Data custody

`data_sources.json` pins authoritative public URLs and cryptographic digests. `download_data.py` downloads to a temporary file, verifies the pinned digest, fsyncs it, and only then atomically installs it. Existing mismatching files are never overwritten. Scientific data are third-party public records and are not altered.

## Integrity

Each authoritative audit directory contains its own source hash, input hashes, output hashes, manifest, report and detailed evidence tables. Absolute paths inside historical records identify the original execution environment; cryptographic digests provide portable custody.

## Citation and authorship

Use `CITATION.cff` for machine-readable citation metadata and `AUTHORS.md` for institutional authorship. Cite DESI DR1 independently when using its public catalogues.

## Licence

MRIPR-authored code and documentation are released under CC0 1.0. DESI catalogue files retain their original ownership, attribution requirements and terms.
