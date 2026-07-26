# Reproducibility and custody protocol

1. Create a clean Python environment and install `requirements.txt`.
2. Run `python download_data.py --verify-only`; download missing public records first when required.
3. Require `OK` for both pinned sources and exactly 175 extracted rotation-curve members.
4. Execute `analysis/audit_procedures/sparc_sigma__audit.py` into a new, uniquely named results directory.
5. Confirm the final status, zero-degree-of-freedom contract, corpus counts, finite-value checks and no-silent-loss predicates.
6. Compare the newly generated result and evidence tables with the authoritative record.
7. Preserve the complete record and compute an archive SHA-256 before publication.

The data loader reads the public records; it does not rewrite them. μ and γ are fixed before SPARC evaluation. Generated CSV tables retain an explicit status for unavailable or invalid catalogue quantities rather than silently removing their galaxies.
