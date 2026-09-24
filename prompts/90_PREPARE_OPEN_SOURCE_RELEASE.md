# Future optional task: prepare a reviewable public-release candidate

This task is queued, not an instruction to publish now. The owner has reserved
the final choice to open-source. Execute after the selected useful milestone
is complete. Prepare all local artifacts before requesting the final license
choice or publication approval. Do not create a public repository, upload data,
send messages to paper authors, or publish a package without explicit authority.

## 1. Define the release honestly

Write a concise README with the actual equations, scope, install/run commands,
example outputs, and scientific limits. Separate code verification, parameter
calibration, and independent experimental validation. The release must not
promise optimal pressure, flavor prediction, a full GS3 digital twin, or 3D
channel formation unless new evidence actually supports those claims.
Publish known failures and unresolved spatial accuracy alongside benchmarks.

## 2. Audit provenance and rights

Inventory copied/adapted code, dependencies, data, figures, documentation, and
AI-assisted contributions. Check every dependency's actual license and bundled
notices. Link papers and cite specific ideas/equations; attribution is not
permission to redistribute their PDF, figures, or datasets. Use our own plots
from released example data where possible. Avoid author/company endorsement
claims and confusing product branding. Do not assume AI-generated output is
automatically free of third-party obligations or guarantees exclusive rights.

Create REFERENCES.md, third-party notices where needed, and CITATION.cff with
accurate human/software attribution. Research-paper authors are inspirations
or sources unless they actually contributed to this software. Describe AI use
transparently in a provenance note without assigning the paper authors
responsibility for our implementation.

Prepare a license-choice memo: MIT for simplicity; Apache-2.0 for an explicit
contributor patent grant; reciprocal licensing only if the owner wants that
tradeoff. Open source permits commercial reuse; it does not mean competitors
must request permission for every use. Do not install a license with an invented
copyright-holder name. Confirm ownership/contributor authority and the final
choice before the release grant. Relevant references:
- https://opensource.org/osd
- https://choosealicense.com/licenses/apache-2.0/
- https://www.copyright.gov/ai/

## 3. Prepare a clean export

Scan tracked/untracked files and any available history for credentials, local
absolute paths, private shot logs, exact serials, account IDs, model transcripts,
and personal metadata. Report findings with redaction, not secret values.
Do not delete private originals or rewrite history casually. Build an explicit
allowlisted candidate export with public synthetic/measured-permissioned data.
Keep the separate marbling workspace completely outside this release.

There was no Git repository at the 2026-09-23 review. If still absent, prepare
a local export/version-control plan; do not fabricate historical commits.
Archive review findings truthfully, sanitize personal filesystem paths in public
copies, and retain private originals for provenance.

## 4. Make reproduction practical

Provide supported Python versions, pinned reproducibility environment or lock,
a fresh wheel-install smoke test, console commands, and deterministic examples.
Verify package data in an installed wheel outside the source tree. Include
lightweight CI for tests and scientific references, then a separate documented
slow validation command. Provide a numerical-method note, parameter schema,
unit conventions, and instructions for adding a new constitutive law safely.
Add CONTRIBUTING.md and issue templates suited to scientific bug reports
(configuration, environment, measured vs simulated result, conservation check).

## 5. Release gate

Deliver a candidate directory/archive, a file manifest, clean-install logs,
test and benchmark evidence, known limitations, attribution/license audit,
release notes, and a short decision memo. Check all documentation commands in
a clean environment. State anything still blocking release. Owner approval of
the concrete candidate and license is the final step before any publication.
