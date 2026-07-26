# MRIPR publication repository template

Every independently publishable MRIPR audit repository must expose the same
top-level interface: scientific audit procedures under `analysis/`, public
inputs and retrieval instructions under `data/`, and immutable evidence under
`results/`.

```text
.
├── README.md              # identity, scope, result, layout and quick start
├── AUTHOR                 # compact author identity
├── AUTHORS.md             # institutional authorship and third-party boundary
├── CITATION.cff           # GitHub/archival citation metadata
├── .zenodo.json           # Zenodo deposition metadata
├── LICENSE                # MRIPR code/docs licence; third-party exclusion
├── CHANGELOG.md           # release history
├── CONTRIBUTING.md        # correction and reproduction protocol
├── SECURITY.md            # private vulnerability reporting
├── REPRODUCIBILITY.md     # exact custody and rerun sequence
├── RESULTS.md             # authoritative record and exact values
├── requirements.txt       # runtime dependencies
├── Makefile               # verify-data, check and audit targets
├── download_data.py       # atomic checksum-verifying downloader
├── data_sources.json      # authoritative URLs and pinned digests
├── analysis/              # scientific procedures and support modules
├── data/
│   └── README.md          # immutable third-party-data boundary
└── results/               # reports, manifests, tables, figures and hashes
```

## Required rules

1. Public scientific data are never edited or silently replaced.
2. Every source has an authoritative URL and pinned checksum or verified extracted aggregate.
3. Downloaders write temporary files, verify them, fsync them and install atomically.
4. Existing mismatching files are never overwritten.
5. Every input row is retained or assigned an explicit exclusion/status reason.
6. Identity constants and degrees of freedom are stated in the README and machine record.
7. The authoritative result directory is named explicitly; historical records remain preserved.
8. Reports are human-readable; JSON/CSV/NPZ and manifests carry machine evidence.
9. Source, input and output hashes are recorded before packaging.
10. Third-party data attribution and licensing remain separate from MRIPR code/documentation licensing.
11. `make verify-data`, `make check` and `make audit` are the uniform public entry points.
12. The release archive checksum is recorded before simultaneous website, arXiv and Zenodo deposition.

## Release gate

- [ ] `python download_data.py --verify-only` passes.
- [ ] `make check` passes in a clean environment.
- [ ] The authoritative audit rerun passes.
- [ ] README values match the authoritative machine manifest.
- [ ] `CITATION.cff` and `.zenodo.json` contain the release title, author, ORCID, institution, version and date.
- [ ] No secrets, caches, macOS sidecars, temporary files or accidental symlinks are present.
- [ ] Archive SHA-256 and publication timestamp are recorded.
- [ ] Website, arXiv and Zenodo artifacts are byte-identical or separately checksummed.
