# GIPR σ Audit Publication Bundle

Publication-ready institutional release of the Mathematical Research Institute of Physical Reality (MRIPR), prepared 25 July 2026.

## Included independent audits

1. **Pantheon+ σ invariant extraction and row closure**
   - 1,701 catalogue rows accounted for
   - 1,578 retained rows
   - 1,578/1,578 finite row extractions
   - maximum row-closure residual: `2.220446049250e-16`
   - fitted degrees of freedom: `0`

2. **SPARC σ full-catalogue audit**
   - 175 galaxies
   - 3,391 rotation-curve rows
   - 3,389 RAR rows
   - complete per-galaxy status accounting
   - fitted degrees of freedom: `0`

3. **DESI σ squared-radial closure and full-corpus audit**
   - 13,097,304 retained physical objects
   - pooled invariant-support coverage: `0.9708459316265691`
   - full-corpus independent arithmetic-route agreement
   - 200 trials in each of three null families
   - observed φ-excess exceeded all 600 null realizations
   - fitted identity degrees of freedom: `0`

## Bundle structure

```text
GIPR_Publication_2026-07-25/
├── README.md
├── AUTHOR
├── AUTHORS.md
├── RELEASE_INFO.json
├── SHA256SUMS
├── PUBLICATION_REPOSITORY_TEMPLATE.md
├── DESI/
├── SPARC/
└── PANTHEON+/
```

Each audit repository is independently reproducible and contains:

- a comprehensive README;
- author and institutional metadata;
- citation and Zenodo metadata;
- licence and third-party-data boundary;
- exact results and reproducibility documentation;
- a repository-local checksum-verifying downloader;
- a scoped public-data source manifest;
- its authoritative audit evidence.

Every repository uses the same top-level contract:

```text
REPOSITORY/
├── analysis/   # scientific audit procedures and supporting methods
├── data/       # public inputs or retrieval instructions
└── results/    # authoritative immutable audit evidence
```

Raw third-party catalogues are intentionally excluded from this publication bundle. They are reproducibly obtained from authoritative public sources with the included downloaders and pinned cryptographic checksums.

## Integrity

`SHA256SUMS` covers every regular file in the bundle except `SHA256SUMS` itself. Verify from the bundle root:

```bash
sha256sum --check SHA256SUMS
```

Do not edit an audit record after checksums are generated. Any revised execution must be written to a new record directory and released with a new bundle checksum.

## Institution

Mathematical Research Institute of Physical Reality  
[https://www.mripr.org/](https://www.mripr.org/)  
[contact@mripr.org](mailto:contact@mripr.org)
