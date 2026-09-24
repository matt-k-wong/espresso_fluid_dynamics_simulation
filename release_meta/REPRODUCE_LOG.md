# Reproduce log (2026-09-23 runs; rebuilt 2026-09-24 for ROADMAP sync, no code changes)

Environment: macOS, Python 3.14.7 (venv), numpy 2.5.3, scipy 1.18.1,
pytest 9.1.1 (see pinned_frozen.txt). Supported range per pyproject:
`requires-python = >=3.10`; this log evidences 3.14.7 only.

## Workspace evidence

- `PYTHONPATH=src ./.venv/bin/python -m pytest tests/ -q` → `112 passed`.
- `PYTHONPATH=src ./.venv/bin/python review/probe_submission.py
  --out /tmp/probe_m6.json` → `17/17 review checks passed`.
- Manifest check from inside candidate: `shasum -a 256 -c MANIFEST.sha256`
  → all listed files OK.

## Fresh wheel-install smoke (candidate wheel, clean target)

- `python -m pip wheel ./release_candidate --no-deps -w /tmp/rcdist` →
  `espresso_m1-0.1.0-py3-none-any.whl` (includes `catalog.json`,
  `shots.py`, `experiment.py`).
- Install with `--no-deps -t /tmp/rctarget`; run with `env -u PYTHONPATH
  PYTHONPATH=/tmp/rctarget` from `/tmp` (source tree unreachable):
  - catalog load → `catalog ok: 34`
  - `cli run` (uniform, Nr=8 Nz=16) → converged, `Q_out = 1.98e-06 m³/s`
  - `cli transport` (tracer, small grid, `t_end=2`) → converged, 7 steps
  - `cli shot-validate` (synthetic_linear) → valid, n=31
  - `cli experiment-save` (linear+irregular, t0=0) → `q=1.2000 g/s`,
    held-out RMSE `0.0000 g`
  - `cli experiment-reproduce` → `PASS` (q, RMSE, mean, hash all ok)

Slow validation (`sh scripts/reproduce.sh`, sweeps + `--wheel` probe) was
not re-run in this task; prior evidence retained in `review/`.
