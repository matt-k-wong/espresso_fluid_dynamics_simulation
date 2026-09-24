# Privacy audit (2026-09-23, local candidate only)

Method: `rg` scans of the workspace excluding `.venv/` and `.pytest_cache/`
for credential/secret/account patterns and for the absolute builder path;
manual listing of `outputs/`, `review/`, `examples/`, `templates/`;
confirmation that no git repository exists.

## Findings (values redacted by reporting pattern classes, not contents)

- No credentials, API keys, tokens, passwords, account IDs, email addresses,
  model transcripts, or private shot logs found in the candidate tree.
- No exact machine serials or personal shot metadata: all example shots are
  `source: synthetic` with placeholder profile IDs; templates carry
  `REPLACE-...` placeholders and an explicit do-not-share-serials note.
- No absolute builder paths (home directory and drive-letter forms) in
  candidate files; the release test asserts this dynamically via
  `os.path.expanduser("~")` so no personal path is written down.
- "serial" matches are guidance text (keep serials private) and JSON
  serialization wording — no serial values present.
- Marbling workspace: not present inside this workspace or candidate
  (kept completely outside per prompt 90).
- `outputs/` (generated run artifacts) is excluded from the candidate via
  allowlist; private originals were not deleted and no history was rewritten
  (no git repository exists: `fatal: not a git repository`).

Redaction note: this report describes pattern classes and counts, never
secret values (none were found to redact).
