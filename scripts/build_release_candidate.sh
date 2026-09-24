#!/bin/sh
# Build an allowlisted local release candidate (no publishing, no history rewrite).
# Run from the repository root: sh scripts/build_release_candidate.sh
set -eu
cd "$(dirname "$0")/.."
DEST="release_candidate"
rm -rf "$DEST"
mkdir -p "$DEST"
copy() { mkdir -p "$DEST/$(dirname "$1")"; cp "$1" "$DEST/$1"; }
copy LICENSE
copy README.md
copy ARCHITECTURE.md
copy REPORT.md
copy REPORT_TRANSPORT.md
copy ROADMAP.md
copy SHOT_DATA.md
copy REFERENCES.md
copy CITATION.cff
copy PROVENANCE_NOTE.md
copy LICENSE_CHOICE_MEMO.md
copy CONTRIBUTING.md
copy pyproject.toml
mkdir -p "$DEST/configs"; cp configs/*.json "$DEST/configs/"
mkdir -p "$DEST/src/espresso_m1"
for f in src/espresso_m1/*.py src/espresso_m1/catalog.json; do cp "$f" "$DEST/$f"; done
mkdir -p "$DEST/examples/shots"; cp examples/shots/synthetic_* "$DEST/examples/shots/"
mkdir -p "$DEST/templates"; cp templates/shot_template.json templates/series_template.csv templates/README.md "$DEST/templates/"
mkdir -p "$DEST/docs/figures"; cp docs/figures/*.svg "$DEST/docs/figures/"
mkdir -p "$DEST/tests"
for f in tests/*.py; do cp "$f" "$DEST/$f"; done
if [ -d tests/fixtures ]; then mkdir -p "$DEST/tests/fixtures"; cp -R tests/fixtures/* "$DEST/tests/fixtures/"; fi
mkdir -p "$DEST/scripts"; cp scripts/reproduce.sh scripts/build_release_candidate.sh "$DEST/scripts/"
mkdir -p "$DEST/review"; cp review/probe_submission.py review/ASTRA_REVIEW.md review/M2_REVIEW_AND_NEXT.md review/GS3_RESULTS_NOTE.md review/submission_results.json review/repaired_results.json review/m2_hydraulics_recheck.json review/gs3_hardening_recheck.json "$DEST/review/"
mkdir -p "$DEST/.github/ISSUE_TEMPLATE"; cp .github/ISSUE_TEMPLATE/scientific-bug-report.md "$DEST/.github/ISSUE_TEMPLATE/"
cp release_meta/PRIVACY_AUDIT.md release_meta/RELEASE_NOTES.md release_meta/DECISION_MEMO.md release_meta/REPRODUCE_LOG.md "$DEST/"
./.venv/bin/pip freeze --exclude-editable > "$DEST/pinned_frozen.txt" 2>/dev/null || ./.venv/bin/pip freeze > "$DEST/pinned_frozen.txt"
cd "$DEST" && find . -type f ! -name "MANIFEST.sha256" -exec shasum -a 256 {} + | sort > MANIFEST.sha256
echo "candidate: $DEST ($(wc -l < MANIFEST.sha256) files)"
