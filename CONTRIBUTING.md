# Contributing (scientific bug reports welcome)

Setup: `python3 -m venv .venv && ./.venv/bin/pip install numpy scipy pytest`
(tested with Python 3.14.7; `requires-python = >=3.10`).
Run: `PYTHONPATH=src ./.venv/bin/python -m pytest tests/ -q`.
Full check: `sh scripts/reproduce.sh`.

## Adding a constitutive law safely

1. Keep the saturated isothermal Darcy core untouched; add the law as a new
   permeability/exchange option behind validation, not a silent default.
2. State units, admissible ranges, and provenance (measured / literature /
   illustrative / fitted) in `catalog.json`.
3. Add an isolated-limit test (analytic or manufactured) before coupling,
   plus ledger/positivity checks for transport changes.
4. Never fit illustrative coefficients to taste or relabel synthetic data
   as measured. Separate verification, calibration, and validation in docs.

## Issue reports (see .github/ISSUE_TEMPLATE)

Include: configuration files/overrides, environment (`python --version`,
`pip list`), expected vs actual result, conservation/balance numbers, and
whether inputs are measured, synthetic, or illustrative.
