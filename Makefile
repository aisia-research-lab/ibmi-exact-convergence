PYTHON ?= python3

.PHONY: install test smoke synthetic public figures verify env

install:
	$(PYTHON) -m pip install -r requirements.txt

test:
	$(PYTHON) -m pytest -q

smoke:
	$(PYTHON) run_reproduction.py --mode smoke

synthetic:
	$(PYTHON) experiments/covariance_validation.py
	$(PYTHON) experiments/ordering_stress.py

public:
	$(PYTHON) scripts/download_suitesparse.py
	$(PYTHON) experiments/public_spd_validation.py

figures:
	$(PYTHON) scripts/make_figures.py

verify:
	$(PYTHON) scripts/verify_paper_outputs.py

env:
	$(PYTHON) scripts/capture_environment.py
