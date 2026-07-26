PYTHON ?= python3
OUTPUT_DIR ?= results/pantheon_sigma_zr_invariant_extraction_row_closure_release_20260725

.PHONY: verify-data check audit

verify-data:
	$(PYTHON) download_data.py --verify-only

check:
	$(PYTHON) -m py_compile download_data.py analysis/sigma_pantheon.py

audit: verify-data check
	$(PYTHON) analysis/sigma_pantheon.py --out_dir $(OUTPUT_DIR)
