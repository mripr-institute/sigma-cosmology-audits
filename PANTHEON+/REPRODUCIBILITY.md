# Reproducibility and custody protocol

1. Create a clean Python environment and install `requirements.txt`.
2. Run `python download_data.py --verify-only`; download the pinned public catalogue if absent.
3. Require the catalogue SHA-256 to equal `1cb0fc379ef066afdc2ffd1857681cc478024570d8a3eba284fb645775198cf8`.
4. Execute `analysis/sigma_pantheon.py` into a new, uniquely named results directory.
5. Confirm all 1,701 input rows are accounted for, 1,578 rows are retained after explicit cuts, and all retained rows have finite Lambert-W extraction.
6. Confirm the maximum row-closure residual remains within the recorded numerical tolerance and fitted degrees of freedom remain zero.
7. Preserve the complete output directory and compute an archive SHA-256 before publication.

The catalogue is read-only evidence. The downloader refuses to overwrite a mismatching existing file. Historical audit records are preserved rather than rewritten.
