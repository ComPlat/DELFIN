# Licence Guard Policy — LICENSE_GUARD.md

Single source of truth for the forbidden patterns:
`scripts/license_guard.py`. Enforced by three layers (defense in depth):
`tests/test_license_guard.py` (CI), the pre-commit hook, and the CLI
audit (`python scripts/license_guard.py --package|--staged|PATH...`).

## Principle

This is a public LGPL-3.0 repository, usable by anyone including
commercially. Data produced by licensed software, or derived from a
licensed database, must never be committed here.

## Rule families

### 1. CCDC/CSD (existing)

CCDC/CSD data is licensed and used ONLY INTERNALLY (validation,
calibration, metrics); it lives outside this repo (agent_workspace/).
See `scripts/license_guard.py` FORBIDDEN_* definitions for the exact
patterns. Open sources (COD, Pyykkö/Cordero radii, VSEPR) stay
shippable.

### 2. Licensed QM program output (LC1, new)

Orca EULA_ORCA_2025 section 3 restricts "DATA" — data generated with
ORCA (SOFTWARE):
- c) DATA may be shared with third parties only for ACADEMIC PURPOSES
  (journal publication permitted);
- d) DATA may go to a database only non-commercially, with a copyright
  notice, an explicit EULA reference and the warranty disclaimer;
- e) any software using DATA must itself be used subject to the EULA;
- g) DATA may feed machine learning only for ACADEMIC PURPOSES.

A public LGPL repository is a commercial-usable third-party
distribution, so ORCA output (.out, .json, .property.json, .hess,
.engrad, scan .dat, .gbw, ORCA-written .xyz/.allxyz) must not be
committed. The guard recognises ORCA output by line signatures
(banner, TERMINATED NORMALLY, version line, run-header lines, INPUT
FILE echo, THERMOCHEMISTRY block) and by artifact paths
(.gbw/.hess/.engrad).

The same rule family covers the other licensed QM programs DELFIN may
drive — Gaussian (ban on redistribution of output data files),
Turbomole, Q-Chem, VASP (output distribution restricted by its own
licence), Molpro (output distribution restricted) — no such fixtures
exist today; the path and signature rules block them if one appears.
xtb and CREST are LGPL-3.0 and their output may ship.

**Sequencing note (phase 2 review, F5)**: the QM-signature content scan
applies to every non-.py tracked file from the moment the class-A
ORCA fixtures leave the tree. Before that removal the committed
fixture inventory (tests/fixtures/lc1_synthetic/INVENTORY.md, and
while they exist the orca611 .out files themselves) is classified, not
hidden; the removal commit in this package switches the guard live.
CI green with the guard on requires that order — guard live before
removal would flag the very files this package documents and removes.

**Exception (the only one besides `# license-guard: allow`)**: a
fixture that declares itself on its FIRST line

```
# SYNTHETIC test fixture written for DELFIN -- not ORCA output
```

is hand-written in the program's text format, invented but physically
consistent, and shippable. The marker is checked as the first line by
the guard; content without the marker is blocked.

### 3. Where licensed data lives instead

The private store outside the repository (operator-managed path).
Tests that genuinely need real program output (format drift between
versions) read it via an environment variable and SKIP when the
variable is absent — they never fail in CI. The variable is set by the
developer/operator locally, never in CI.

The store is managed with `delfin-dev fixtures` (`delfin/
fixtures_store.py`): `init <dir>` creates the layout, `add <src> --as
<rel>` copies a file into a repository-mirroring relative path and
records sha256, ORCA version and provenance in `manifest.json`,
`verify` re-checks the hashes and the no-personal-data rule, `where`
prints the resolved store root. The root comes from the
`DELFIN_PRIVATE_FIXTURES` environment variable, then from
`paths.fixtures_dir` in the DELFIN settings file — never a hard-coded
path. The test helper (`tests/test_lc1_fixtures_helper.py`) refuses a
file whose sha256 does not match the manifest. Tests that use the
store carry the `private_fixtures` marker (auto-applied to every
module importing the helper), so a local run can select them with
`pytest -m private_fixtures`. See docs/PRIVATE_FIXTURE_STORE.md.

## Escape hatch

`# license-guard: allow` on the offending line — for a deliberate,
documented reference. Use sparingly.

## History note

Git history is not rewritten by this package. Files removed here stay
in older commits; deciding on a history rewrite is the person's
decision.
