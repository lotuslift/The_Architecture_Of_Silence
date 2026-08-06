.PHONY: verify checks all

verify:
	python scripts/verify_release.py

checks:
	python scripts/run_all_checks.py

all: verify checks
