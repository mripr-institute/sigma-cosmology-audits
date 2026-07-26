PYTHON ?= python3
DESI_DATA_ROOT ?= data/desi
DESI_FULL_RHO_SQUARED_OUTPUT ?= results/sigma_desi_full_rho_squared_radial_release_20260725
DESI_CLOSURE_SQUARED_OUTPUT ?= results/sigma_desi_lss_closure_squared_radial_release_20260725
SQUARED_CLOSURE_MANIFEST ?= $(DESI_CLOSURE_SQUARED_OUTPUT)/sigma_desi_lss_closure_squared_radial_manifest.json
DESI_FULL_RHO_SQUARED_FINAL_OUTPUT ?= results/sigma_desi_full_rho_squared_radial_release_20260725
AUDIT_PYTHONPATH ?=
AUDIT_MPLCONFIGDIR ?= results/.mplconfig
NULL_WORKERS ?= 12

.PHONY: verify-data check audit desi-full-rho-squared desi-closure-squared desi-squared-all
verify-data:
	$(PYTHON) download_data.py --verify-only

check:
	$(PYTHON) -m py_compile download_data.py analysis/sigma_desi_lss_closure_squared_radial_audit.py analysis/sigma_desi_full_rho_squared_radial_audit.py
	$(PYTHON) analysis/sigma_desi_lss_closure_squared_radial_audit.py --synthetic-test --output-dir $(DESI_CLOSURE_SQUARED_OUTPUT)
	$(PYTHON) analysis/sigma_desi_full_rho_squared_radial_audit.py --synthetic-test --output-dir $(DESI_FULL_RHO_SQUARED_FINAL_OUTPUT)

audit: desi-squared-all

desi-full-rho-squared:
	PYTHONPATH=$(AUDIT_PYTHONPATH) MPLCONFIGDIR=$(AUDIT_MPLCONFIGDIR) $(PYTHON) analysis/sigma_desi_full_rho_squared_radial_audit.py --data-root $(DESI_DATA_ROOT) --completed-closure $(SQUARED_CLOSURE_MANIFEST) --output-dir $(DESI_FULL_RHO_SQUARED_OUTPUT) --mode full --seed 12345 --q-route fast --log-level normal

desi-closure-squared:
	PYTHONPATH=$(AUDIT_PYTHONPATH) MPLCONFIGDIR=$(AUDIT_MPLCONFIGDIR) $(PYTHON) analysis/sigma_desi_lss_closure_squared_radial_audit.py --data-root $(DESI_DATA_ROOT) --output-dir $(DESI_CLOSURE_SQUARED_OUTPUT) --seed 12345 --null-trials 200 --null-sub 500000 --null-workers $(NULL_WORKERS)

desi-squared-all:
	$(PYTHON) -m py_compile analysis/sigma_desi_lss_closure_squared_radial_audit.py analysis/sigma_desi_full_rho_squared_radial_audit.py
	PYTHONPATH=$(AUDIT_PYTHONPATH) MPLCONFIGDIR=$(AUDIT_MPLCONFIGDIR) $(PYTHON) analysis/sigma_desi_lss_closure_squared_radial_audit.py --synthetic-test --output-dir $(DESI_CLOSURE_SQUARED_OUTPUT)
	PYTHONPATH=$(AUDIT_PYTHONPATH) MPLCONFIGDIR=$(AUDIT_MPLCONFIGDIR) $(PYTHON) analysis/sigma_desi_full_rho_squared_radial_audit.py --synthetic-test --output-dir $(DESI_FULL_RHO_SQUARED_FINAL_OUTPUT)
	PYTHONPATH=$(AUDIT_PYTHONPATH) MPLCONFIGDIR=$(AUDIT_MPLCONFIGDIR) $(PYTHON) analysis/sigma_desi_lss_closure_squared_radial_audit.py --data-root $(DESI_DATA_ROOT) --output-dir $(DESI_CLOSURE_SQUARED_OUTPUT) --seed 12345 --null-trials 200 --null-sub 500000 --null-workers $(NULL_WORKERS)
	PYTHONPATH=$(AUDIT_PYTHONPATH) MPLCONFIGDIR=$(AUDIT_MPLCONFIGDIR) $(PYTHON) analysis/sigma_desi_full_rho_squared_radial_audit.py --data-root $(DESI_DATA_ROOT) --completed-closure $(SQUARED_CLOSURE_MANIFEST) --output-dir $(DESI_FULL_RHO_SQUARED_FINAL_OUTPUT) --mode full --seed 12345 --q-route compare --log-level normal
