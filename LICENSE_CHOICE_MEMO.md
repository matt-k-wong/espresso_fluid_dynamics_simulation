# License-choice memo (decided 2026-09-24: Apache-2.0, holder matt-k-wong)

LICENSE installed with `Copyright 2026 matt-k-wong`; CITATION.cff set to
`Apache-2.0`. Publication itself remains a separate owner decision; nothing
here publishes code, uploads data, or grants rights beyond the license text.

## Options (per prompt 90)

- MIT: simplest, permissive; no explicit patent grant language.
- Apache-2.0: permissive with an explicit contributor patent grant; preferred
  if the owner wants that grant spelled out.
- Reciprocal (e.g. GPL-family): only if the owner affirmatively wants the
  copyleft tradeoff. Note SciPy binary bundles already carry notices
  (GCC runtime exception, LGPL libquadmath) relevant to redistribution.

Open source permits commercial reuse; it does not mean competitors must ask
permission for every use.

## Required owner inputs before any release grant

1. Confirm ownership/contributor authority for all files in the candidate.
2. Fill in the copyright holder name(s) and year (do not invent).
3. Choose one license, add its full text as LICENSE, and set the
   CITATION.cff `license` field accordingly.
4. Re-run the clean-export + fresh-install checks after adding LICENSE.
