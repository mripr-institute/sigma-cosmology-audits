PYTHON ?= python3
OUTPUT_DIR ?= results/sparc_sigma_full_catalogue_release_20260725

.PHONY: verify-data check audit

verify-data:
	$(PYTHON) download_data.py --verify-only

check:
	$(PYTHON) -m py_compile download_data.py analysis/audit_procedures/sparc_sigma__audit.py analysis/shared/sparc_data.py analysis/shared/pub_style.py analysis/sparc/sigma_rar_identity.py analysis/sparc/sigma_sparc_identity.py

audit: verify-data check
	$(PYTHON) analysis/audit_procedures/sparc_sigma__audit.py --output-directory $(OUTPUT_DIR)
