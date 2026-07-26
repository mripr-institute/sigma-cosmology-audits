# SPARC σ Full-Catalogue Audit

Official reproducibility repository of the Mathematical Research Institute of Physical Reality (MRIPR) for the zero-degree-of-freedom SPARC σ audit.

## Identity and audit contract

The procedure evaluates the fixed identity

\[
\sigma(r)=\ln\!\left[\frac{\mu^2(1+\gamma r)}{\mu+\gamma}\right]-\mu r,
\qquad \rho(r)=e^{\sigma(r)},
\]

using \(\mu=0.082912607552\) and \(\gamma=0.38603416\). No SPARC-fitted value changes either invariant or the identity.

## Authoritative result

The publication record is:

`results/sparc_sigma_full_catalogue_corpus_independent_20260722/`

Recorded result: **PASS**.

- 175 galaxies accounted for with no silent loss
- 3,391 rotation-curve points
- 135/135 galaxies with valid `Vflat` close within 10% at the last measured radius
- 3,389 RAR points
- RAR RMS scatter: `0.195995` dex
- Pearson correlation: `0.936046`
- Spearman correlation: `0.929916`
- fitted degrees of freedom: `0`

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
│   ├── audit_procedures/      # institutional audit entry point
│   ├── shared/                # catalogue loader and publication style
│   └── sparc/                 # supporting identity procedures
├── data/sparc/                # public SPARC records; never altered
└── results/                   # immutable audit evidence
```

## Reproduce

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python download_data.py
.venv/bin/python download_data.py --verify-only
.venv/bin/python analysis/audit_procedures/sparc_sigma__audit.py \
  --output-directory results/sparc_sigma_full_catalogue_release
```

Use `--help` on the audit entry point if the installed procedure exposes additional execution controls.

## Evidence products

The authoritative record contains:

- a human-readable audit report;
- a machine-readable result and manifest;
- complete per-galaxy outer-flatness accounting;
- the σ derivative profile for every retained rotation-curve row;
- figure audit products and their predicates;
- input, procedure and output hashes.

## Data custody

`data_sources.json` pins the SPARC public property table and rotation-curve archive. The downloader verifies the property-table SHA-256 and the upstream archive MD5, requires exactly 175 extracted rotation-curve members, refuses to overwrite mismatching files, and never modifies scientific records.

## Citation and licence

Use `CITATION.cff` for this audit and cite the original SPARC release independently. MRIPR-authored code and documentation are CC0 1.0; SPARC data retain their original terms and required attribution.
