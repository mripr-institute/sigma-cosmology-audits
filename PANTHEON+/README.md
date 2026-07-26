# Pantheon+ σ Invariant Extraction and Row-Closure Audit

Official reproducibility repository of the Mathematical Research Institute of Physical Reality (MRIPR) for the Pantheon+ `z = r` invariant extraction and row-closure audit.

## Identity

\[
\sigma(r)=\ln\!\left[\frac{\mu^2(1+\gamma r)}{\mu+\gamma}\right]-\mu r,
\qquad
\rho(r)=\frac{\mu^2}{\mu+\gamma}(1+\gamma r)e^{-\mu r}.
\]

The Pantheon+ audit records

\[
\mu=0.082912607552,
\qquad
\gamma=0.38603416,
\qquad
r_*=\frac1\mu-\frac1\gamma=9.470447610694.
\]

The catalogue row rule is `r = catalogue_z`. γ is recovered row-by-row through the closed-form Lambert-W route. The identity has zero fitted degrees of freedom.

## Authoritative result

The publication record is:

`results/pantheon_sigma_zr_invariant_extraction_row_closure_20260723/`

Recorded result: **PASS**.

- 1,701 total catalogue rows accounted for
- 1,578 retained rows after explicit calibrator and redshift cuts
- 1,578/1,578 finite Lambert-W γ recoveries
- maximum row-closure residual: `2.220446049250e-16`
- derivative reconstruction maximum error below `5e-120`
- fitted degrees of freedom: `0`

The derivative and stationary reconstructions are supplementary identity-consistency evidence. The corpus PASS is based on complete accounting, finite row extraction, exact `z = r`, evidence-product creation, zero degrees of freedom and numerical row closure.

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
│   └── sigma_pantheon.py        # audit entry point
├── data/PantheonSH0ES.dat      # public catalogue; never altered
└── results/                    # immutable audit evidence
```

## Reproduce

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python download_data.py
.venv/bin/python download_data.py --verify-only
.venv/bin/python analysis/sigma_pantheon.py \
  --out_dir results/pantheon_sigma_zr_invariant_extraction_row_closure_release
```

## Data custody

The repository-local manifest pins the authoritative Pantheon+SH0ES public URL and SHA-256. The downloader verifies before atomic installation and refuses to overwrite any mismatching existing file. The catalogue is public third-party evidence and is never rewritten by the audit.

## Evidence products

The authoritative record contains a Markdown report, JSON result, JSON manifest, complete row-evidence CSV and supplementary reconstruction CSV. File and procedure hashes are embedded in the custody record.

## Citation and licence

Use `CITATION.cff` for the MRIPR audit and cite the Pantheon+SH0ES data release independently. MRIPR-authored code and documentation are CC0 1.0; Pantheon+ data retain their original terms and attribution requirements.
