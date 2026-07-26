# MRIPR σ Cosmology Audits

**Author:** Alex Albert\
**ORCID:** https://orcid.org/0009-0005-6981-2087 \
**Institution:** Mathematical Research Institute of Physical Reality\
**Release date:** 25 July 2026

This bundle contains the public audit records for the Pantheon+, SPARC,
and DESI evaluations of σ. Each audit is independently reproducible and
preserves the code, procedures, numerical results, and integrity records
required to verify the published execution.

## Included audits

### Pantheon+ σ invariant extraction and row closure

-   Catalogue rows accounted for: `1,701`
-   Rows retained: `1,578`
-   Finite row extractions: `1,578 / 1,578`
-   Maximum row-closure residual: `2.220446049250e-16`
-   Fitted degrees of freedom: `0`

### SPARC σ full-catalogue audit

-   Galaxies accounted for: `175`
-   Rotation-curve rows: `3,391`
-   RAR rows: `3,389`
-   Per-galaxy status accounting: complete
-   Fitted degrees of freedom: `0`

### DESI σ squared-radial closure and full-corpus audit

-   Retained physical objects: `13,097,304`
-   Pooled invariant-support coverage: `0.9708459316265691`
-   Independent arithmetic routes: full-corpus agreement
-   Null families: `3`
-   Trials per null family: `200`
-   Total null realizations: `600`
-   Observed φ-excess: greater than every null realization
-   Fitted identity degrees of freedom: `0`

## Bundle contents

``` text
MRIPR_Publication_2026-07-25/
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

The three audit repositories contain:

-   the audit code and supporting procedures;
-   the published numerical results;
-   author, institutional, citation, and release metadata;
-   licensing information and third-party data boundaries;
-   exact reproduction instructions;
-   checksum-verifying data downloaders;
-   manifests identifying the public source data;
-   the authoritative evidence produced by the audit.

Each repository follows the same top-level layout:

``` text
REPOSITORY/
├── analysis/
├── data/
└── results/
```

## Third-party data

Raw third-party catalogues are not redistributed in this bundle.

The required inputs are obtained from their authoritative public sources
by the included downloaders. Expected files are identified by pinned
cryptographic checksums so that the retrieved inputs can be verified
before an audit is executed.

## Reproduction

Each audit repository includes its own instructions and verification
procedure. The repositories can therefore be reproduced separately; the
enclosing bundle records the institutional release in which their
results were published together.

The published files under `results/` are audit records rather than
working directories. A new execution must be written to a new record
directory. Existing records must not be overwritten or edited in place.

## Bundle integrity

`SHA256SUMS` covers every regular file in the bundle except the checksum
file itself.

``` bash
sha256sum --check SHA256SUMS
```

A successful verification confirms that the local files match the files
included in this release.

Any revision requires:

1.  a new record directory;
2.  a new release identifier;
3.  regenerated checksums;
4.  a newly issued publication bundle.

## Author

**Dr. Alex Albert**\
Mathematical Research Institute of Physical Reality

## Institution

**Mathematical Research Institute of Physical Reality**

Website: https://www.mripr.org/\
Email: contact@mripr.org
