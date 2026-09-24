# Decision memo (owner approval required, nothing published)

## Candidate

Local directory `release_candidate/` (77 allowlisted files + this memo set;
see MANIFEST.sha256). No public repository created, no data uploaded, no
package published, no authors contacted.

## Checks completed

- Privacy audit clean (redacted pattern report, no values).
- Provenance/references/citation files present; AI use disclosed.
- Fresh wheel-install smoke test passed from the candidate wheel.
- Suite 112 passed; probe 17/17; demo numbers match reports.

## Blockers / still required

1. [done 2026-09-24] Owner confirms ownership/contributor authority: all
   code written in-session, no third-party code pasted in.
2. [done 2026-09-24] Copyright holder + year: `Copyright 2026 matt-k-wong`.
3. [done 2026-09-24] License chosen: Apache-2.0; LICENSE text installed;
   CITATION.cff `license` set.
4. Re-run `sh scripts/build_release_candidate.sh` + smoke test after (3)
   (in progress below).
5. Explicit owner approval of the concrete candidate + license before any
   publication step (separate action, not taken here).
