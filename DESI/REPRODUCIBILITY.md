# Reproducibility and custody protocol

1. Use a clean Python environment and install `requirements.txt`.
2. Run `python download_data.py --verify-only`. If files are absent, run without `--verify-only`, then verify again.
3. Confirm that every source reports `OK`; do not substitute catalogue editions.
4. Run the closure audit with seed `12345`, 41 scan positions, 200 trials per null family, a 500,000-object null subsample and the optimized exact-score null engine.
5. Supply the resulting closure manifest directly to the full-corpus audit.
6. Use `--q-route compare` for publication custody so both arithmetic routes are evaluated across the full retained corpus.
7. Preserve the complete output directories. Do not edit generated manifests, reports, CSV, NPZ or figures after execution.
8. Record archive-level SHA-256 checksums after packaging and before upload to the website, arXiv and Zenodo.

The scientific catalogues are immutable inputs. The downloader refuses to overwrite a mismatching existing file. Generated records contain source, input and output hashes sufficient to distinguish code changes, data changes and post-run artifact changes.

The identity parameters are fixed before the DESI execution. The offset scan is a geometric-location scan and is repeated identically in the null families; it does not fit μ, γ or σ.
