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

### 0. Share-alike reference databases (LC2)

A published reference set is INPUT to this software, never content of
it. The case that prompted the rule: QUESTDB, the obvious source of
reference vertical excitation energies (GitHub pfloos/QUESTDB, Zenodo
10.5281/zenodo.15671385; Loos, Boggio-Pasqua, Blondel, Lipparini and
Jacquemin, J. Chem. Theory Comput. 2025, 21, 8010) is CC-BY-SA-4.0.
Share-alike plus this repository's LGPL-3.0 is the conflict the
principle above forbids, and EU database rights apply on top.

So `delfin/agent/reference_oracle.py` takes its references as an
ARGUMENT -- a path to a file the user supplies -- and
`tests/test_a_reference_the_agent_cannot_author.py` asserts the module
embeds no data of its own. The forbidden pattern is a PATH rule and
deliberately not a content rule: the comparator CITES the database, and
a citation is not a copy.

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

### 4. Third-party ML potentials and model weights (family 4)

The MLIP registry (`delfin/mlp_tools/registry.py`) can call optional
third-party packages (orca-mlips, fairchem, mace, orb-models,
aimnetcentral, gpuma, torchani, ...) and load model weights. CODE licence
and WEIGHT licence are separate things and both are declared per backend
(`code_licence` / `weight_licence`, enforced by
`tests/test_mlip_licence_records.py`).

- **No weights or checkpoints in the repo.** Never commit, vendor or
  download model files into this tree: `*.pt`, `*.pth`, `*.ckpt`,
  `*.safetensors`, `*.npz`, or any trained-model file. Weights stay in the
  user's own cache (HF cache, backend download cache); a test blocks
  tracked weight files.
- **No tokens.** DELFIN never reads, prints, logs or stores Hugging Face
  or other access tokens. Gated checkpoints (UMA) are the user's own
  acceptance + their own login, outside DELFIN.
- **Gated weights** need the user's own acceptance of the upstream terms.
  DELFIN's job is to say so, with the URL, not to work around the gate.
- **Non-commercial weights** (e.g. MACE foundation models under the
  Academic Software License) are flagged `non_commercial=True`; the
  option `refuse_noncommercial_weights=True` (on `mlp_preoptimize` and
  `registry.create_calculator`) makes DELFIN refuse such backends so
  commercial users cannot trip over them.
- **Every optional MLIP dependency DELFIN can call must have a cited
  upstream licence** in the registry; anything not yet verified against
  the upstream text is marked `verified=False`, never guessed.
- **What an installer may do:** hint only (`pip install <pkg>` in an
  install_hint string). DELFIN never installs an MLIP package or weights
  on its own. The guided set-up (`delfin/mlp_tools/guided_setup.py`,
  `delfin mlp_setup`) shows the exact commands with the backend's code
  and weight licence records and runs a command only after the user
  approved that command (terminal y/N, or the agent's shell gate, which
  asks for every installer command in every mode). It installs packages
  into a separate environment; weights are fetched by the backend into
  the user's own cache on first use (the verify job), never into the
  repository. For a gated checkpoint the plan tells the user to run
  `hf auth login` themselves; the token never passes through DELFIN.
- Every backend's `weight_licence.url` names the model card / licence page
  the user must check before depending on the checkpoint.

## Escape hatch

`# license-guard: allow` on the offending line — for a deliberate,
documented reference. Use sparingly.

## History note

Git history is not rewritten by this package. Files removed here stay
in older commits; deciding on a history rewrite is the person's
decision.

## Bring-your-own-licence tools

A small number of scripts under `tools/` exist for users who hold their own CSD licence. They read
the user's local CSD installation and produce, on the user's machine, a file DELFIN can be pointed
at. The rules for every such tool:

- It lives in `tools/` only. No module under `delfin/` imports it (checked by
  `tests/test_byol_tools.py`), and DELFIN runs without it.
- Its single CSD import carries `# license-guard: allow` and is executed only under the CSD
  Python API interpreter (`run_csd_python_api`).
- What it produces stays local. The derived table is never committed to this repository and is
  not redistributed; shipping it would need the written consent of the CCDC.
- The shipped default stays open data: the metal-donor target table delivered with MANTA is the
  COD (CC0) table `delfin/manta/data/ml_targets_cod.tsv`.

Current tools: `tools/export_csd_crystals.py` (exports single-complex crystal coordinates for
`tools/derive_ml_targets.py`, which derives the bond-valence metal-donor table and is licence-free
itself).

## Stronger copyleft than this repository: the AGPL process boundary

DELFIN is **LGPL-3.0** and this repository is public and usable
commercially. The Chemotion ELN work touches two **AGPL-3.0** projects:

- `chemotion_api` (`github.com/Chemotion/chemotionapi`), the Python
  client for the ELN's REST/Grape API;
- `ChemAssistent` (`github.com/StarmanMartin/ChemAssistent`), a
  colleague's assistant built on it.

AGPL-3.0 is stronger copyleft than LGPL-3.0. Linking either into
DELFIN would make the combined work AGPL — which this repository
cannot become by importing a dependency, and which would also reach
everyone who uses DELFIN **over a network**. It is a licence
incompatibility for distribution, not a preference.

**The rule.** `delfin/` never imports `chemotion_api`, and
`pyproject.toml` never declares it in `dependencies` or in any extra.
The ELN is reached through a **separate process** that DELFIN talks to
over a local boundary (subprocess or loopback IPC). That process is its
own work under AGPL-3.0, with its own source offered as the AGPL
requires; DELFIN ships a client for the boundary, not the library.
Enforced by `tests/test_the_agpl_boundary_holds.py`.

**Prior art, not source.** Ideas learned from reading `ChemAssistent`
(confirm-before-act, graded autonomy, the credential living outside the
model's reach, failing closed on anything unrecognised) are design
decisions and may be used freely. Its **code** may not be copied,
adapted or paraphrased into this repository under any arrangement,
including one meant not to be noticed: the obligation follows the work,
not the resemblance. Where its solution is the right one, the honest
routes are to depend on it across the process boundary, to contribute
upstream, or to ask the author.

**If this is ever revisited**, the decision is the user's and the entry
here changes first. Three options were weighed on 2026-11-19: accept
AGPL for DELFIN (refused — it is LGPL and public); the separate process
(chosen); a thin independent HTTP layer against the ELN's routes
(cheapest on licence, most expensive on safety, since the client's own
guards against whole-record overwrites and restricted elements would
have to be rebuilt).


## Test-suite dependencies: optional, permissively licensed

Three ordinary installable packages gate parts of the test suite. None is
one of the heavy optional builders (MACE, molSimplify, Molassembler,
architector, rxembed, metallogen); each is small, permissively licensed and
declared as an **optional extra**, so a machine without it runs the rest of
the suite and the gated tests skip by naming that package. The licence of
each is recorded here so that the compatibility with this repository's
LGPL-3.0 is a cited decision, not an assumption (the family-4 rule for
optional third-party dependencies, applied to the test stack):

- `hypothesis` — **MIT**. Property-based tests of the outside-input parsers
  (`tests/test_outside_input_parsers_hold_on_any_text.py`). Declared
  pinned in the `dev` extra (`hypothesis==6.168.5`); the exact pin is
  deliberate, because a new release changes the inputs a derandomized run
  generates.
- `morfeus` — **MIT**. The steric-descriptor package, imported as `morfeus`
  and shipped on PyPI as `morfeus-ml`. Gates the steric wrapper and adapter
  tests. Declared in the `analysis` extra.
- `playwright` — **Apache-2.0**. The headless-browser driver for the page
  and viewer layout tests, and for the turn-end page check
  (`delfin/agent/artifact_check.py`). Declared in the `browser` extra
  (`playwright>=1.40,<2`); the browser itself is a one-time
  `playwright install chromium`, separate from the package.

MIT and Apache-2.0 are both permissive and compatible with the
repository's LGPL-3.0, so depending on any of them does not move the
combined work's licence. A test that needs one skips by naming that
package; the skip is legitimate on a machine that does not have it, and the
test must not silently depend on anything else.
