#!/bin/sh
# Easy startup: build the demo experiment bundle + static viewer, then open it.
# Usage: sh scripts/view.sh
# Output goes to outputs/site_demo/ (ignored by .gitignore). No server, no uploads.
set -eu
cd "$(dirname "$0")/.."

PYTHONPATH=src ./.venv/bin/python -m espresso_m1.cli experiment-save \
  --shots examples/shots/synthetic_linear.json examples/shots/synthetic_irregular.json \
  --already-wet-t0 0.0 --label demo --out-dir outputs/site_demo/exp_demo

PYTHONPATH=src ./.venv/bin/python -m espresso_m1.cli experiment-site \
  --bundles outputs/site_demo/exp_demo \
  --out outputs/site_demo/index.html --title "Espresso experiments (demo)"

echo "viewer: outputs/site_demo/index.html"
if command -v open >/dev/null 2>&1; then
  open outputs/site_demo/index.html
elif command -v xdg-open >/dev/null 2>&1; then
  xdg-open outputs/site_demo/index.html
else
  echo "open outputs/site_demo/index.html in a browser"
fi
