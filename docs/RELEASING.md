# Releasing DELFIN

How a version gets from a tag to PyPI and a GitHub release, what is
verifiable about the result, and the repository and PyPI settings the
release workflow depends on. Branch protection and the `v*` tag ruleset
(a pushed release tag cannot be moved or deleted) are in
`docs/MAIN_PROTECTION.md`.

## What happens on a tag

`.github/workflows/release.yml` runs on a pushed `v*` tag (and on dispatch,
see below). Five jobs; the workflow grants no permission, each job names its
own:

| job | permissions | does |
|---|---|---|
| `build (sdist, wheel)` | `contents: read` | checks the tag equals `v` + `project.version` in `pyproject.toml`; `python -m build` (pip, build, setuptools and wheel installed with `--require-hashes` from `.github/ci/release-requirements.txt`, no pip cache); `scripts/check_distribution.py --release dist/*`, the checker packaging.yml runs on every pull request: exactly one wheel and one sdist; no file in an installer's output directory (`qm_tools/bin/`, `qm_tools/share/`, …) and no file named like a program the installers download; no ELF, Mach-O or PE image under any name and no native-library suffix; `LICENSE`, `COPYING`, `COPYING.LESSER` and `THIRD_PARTY_NOTICES.md` in the wheel's `dist-info/licenses/` and the sdist's root, with `License-Expression: LGPL-3.0-or-later` and their `License-File` lines in the metadata |
| `SBOM (CycloneDX)` | `contents: read` | installs the built wheel into an empty Python 3.11 environment as a user would and writes a CycloneDX 1.6 JSON SBOM of it (`cyclonedx-py environment`, generator installed with `--require-hashes` from `.github/ci/sbom-requirements.txt`, reproducible: timestamp from the commit). Separate from `build` because the dependency install runs code from PyPI |
| `provenance` | `contents: read`, `id-token: write`, `attestations: write` | `actions/attest-build-provenance`: a SLSA provenance attestation for each file; `actions/attest` with `sbom-path`: an SBOM attestation for each file. Both signed through Sigstore with the workflow's OIDC identity, stored in the repository. Runs no project code |
| `publish to PyPI` | `id-token: write` | in the `pypi` environment; `pypa/gh-action-pypi-publish` exchanges the OIDC token for a short-lived PyPI token (Trusted Publishing: no stored secret), runs `twine check`, uploads, and uploads PEP 740 attestations |
| `GitHub release` | `contents: write` | the same files, the SBOM (`delfin_complat-<version>.cdx.json`) and `SHA256SUMS.txt` (covering both) on the tag's release; creates the release from the annotated tag's message if there is none. Never replaces an existing asset |

`publish to PyPI` and `GitHub release` both depend on `build` and
`provenance`, not on each other: before PyPI is configured, a tag still gets
its release assets, and each job's own result says which half succeeded.

A dispatch (Actions, Release, Run workflow) builds, checks and attests the
selected ref. It publishes only with **publish** ticked and a `v*` tag as
the ref; on a branch it is a dry run.

## Cutting a release

1. Set `version` in `pyproject.toml` (and the citation metadata) and merge.
2. Tag the merge commit with an annotated tag; its message becomes the
   release notes: `git tag -a v1.3.4 -m "..."`, then
   `git push origin v1.3.4`. Once pushed, the tag cannot be moved or
   deleted (`protect-release-tags`, `docs/MAIN_PROTECTION.md`).
3. Approve the `pypi` deployment when the run asks for it.

## Verifying a release

```sh
pip download delfin-complat==1.3.4 --no-deps -d dl/
gh attestation verify dl/delfin_complat-1.3.4-py3-none-any.whl -R ComPlat/DELFIN
gh attestation verify dl/delfin_complat-1.3.4-py3-none-any.whl -R ComPlat/DELFIN \
    --predicate-type https://cyclonedx.org/bom     # the SBOM attestation
sha256sum -c SHA256SUMS.txt     # from the GitHub release
```

`gh attestation verify` checks that the file was built by
`release.yml` in ComPlat/DELFIN, from which commit and ref; with
`--predicate-type https://cyclonedx.org/bom` it checks the SBOM attested
for that file instead. PyPI shows the PEP 740 provenance on each file's
page. The SBOM lists the dependency versions a plain install resolved on
the release day; a later install can resolve newer ones inside the
ranges of `pyproject.toml`.

## One-time settings (maintainer)

None of this is in the repository; until it is set, `publish to PyPI` fails
and the other three jobs still run.

1. **PyPI trusted publisher.** pypi.org, project `delfin-complat`, Manage,
   Publishing, Add a new publisher, GitHub: owner `ComPlat`, repository
   `DELFIN`, workflow `release.yml`, environment `pypi`. After the first
   publish through it, remove any API token kept for uploads.
2. **The `pypi` environment.** Settings, Environments, New environment
   `pypi`: Required reviewers (at least one maintainer other than the one
   who tags; tick "Prevent self-review"), and Deployment branches and tags:
   Selected, tag pattern `v*`.
3. **Allow the publish action.** The repository allows actions owned by
   GitHub only, and `pypa/gh-action-pypi-publish` is not. Settings, Actions,
   General, Actions permissions: "Allow ComPlat, and select non-ComPlat,
   actions and reusable workflows", keep "Allow actions created by GitHub",
   and add `pypa/gh-action-pypi-publish@dc37677b2e1c63e2034f94d8a5b11f265b73ba33`
   (v1.14.2). An entry naming the SHA must move with each Dependabot bump of
   the action; `pypa/gh-action-pypi-publish@*` together with item 4 avoids
   that and still runs only pinned commits.
4. **Require SHA pins.** Same page: "Require actions to be pinned to a
   full-length commit SHA". Every `uses:` already is
   (`tests/test_every_action_is_pinned_to_a_commit.py`).
5. **Who can start a release.** A `v*` tag push publishes. A tag ruleset
   (Settings, Rules, Rulesets, New tag ruleset, target `v*`, restrict
   creations, updates and deletions, bypass for maintainers) limits that to
   maintainers; the environment's reviewers are the second barrier.
6. **Optional: immutable releases** (Settings, General, Releases). Assets
   cannot be added to a published release then, so let the workflow create
   the release: push the tag, do not create the release in the web
   interface first.

## Supply chain: what is pinned and where it is checked

| what | pinned by | updated by | checked by |
|---|---|---|---|
| workflow actions | full commit SHA + `# vX.Y.Z` | Dependabot `github-actions`, 7-day cooldown | `tests/test_every_action_is_pinned_to_a_commit.py` |
| CI tools (pip, pytest, ruff, pip-audit, build, coverage) | `.github/ci/constraints-ci.txt` | by hand | `tests/test_ci_installs_what_it_names.py` |
| release build tools (pip, build, setuptools, wheel + closure), SBOM generator (cyclonedx-bom + closure) | version + SHA-256 of every PyPI file, `.github/ci/release-requirements.txt`, `.github/ci/sbom-requirements.txt` (`pip install --require-hashes`) | by hand: `.github/ci/hashed_requirements.py`, command in each file's header (7-day cooldown) | `tests/test_the_release_tools_are_installed_by_hash.py` |
| micromamba in CI | version 2.9.0-0 + SHA-256 | by hand | `sha256sum -c` in the job |
| xtb release binary (slow tests) | version 6.7.1 + SHA-256 | by hand | `sha256sum -c` in the job |
| Bun | version + SHA-256 | by hand | `tui/scripts/install-bun.sh` |
| TUI packages | `tui/bun.lock` (sha512 per entry) | Dependabot `bun` (minor/patch, not `@opentui/*`) | `bun install --frozen-lockfile` |
| xtb4stda, stda, their parameter files | tag v1.0 / its commit + SHA-256 per file | by hand | `fetch_verified` in `delfin/qm_tools/install_qm_tools.sh` |
| Ketcher drawing editor (fetched on first use, served by the dashboard) | version 3.18.0 + size + SHA-256 in `delfin/dashboard/ketcher.py` | by hand | `tests/test_the_editor_download_is_pinned_and_unpacked_inside.py`; `install()` refuses a mismatch and any archive member outside its directory |
| agent container image (`agent.bash_isolation = "container"`) | `python@sha256:` digest in `delfin/agent/container_exec.py` | Dependabot `docker` on `.github/ci/agent-container/Dockerfile`, 7-day cooldown; the PR stays red until the constant matches | `tests/test_what_the_agent_runs_on_is_named_exactly.py` |
| Python dependencies | ranges in `pyproject.toml`; voila, jupyter-server and tornado floored at the releases that fix their known advisories | Dependabot `pip` (grouped) | `dependency-review.yml` on PRs; `audit.yml` nightly |

`audit.yml` fails on a known advisory in the closure of `.[agent,docs]`; an
advisory accepted on purpose goes in `.github/ci/audit-ignore.txt` with its
reason. CodeQL analyses Python, the workflows, the TypeScript window
(`tui/src`) and the C# launcher.

Not yet pinned: the micromamba bootstraps of the installers under
`delfin/` and the Settings tab (`micro.mamba.pm/.../latest`), and the
conda-forge packages the toolchains job installs (versions float; conda
checks each package's digest against the channel's repodata).

## Releases already on PyPI

`delfin_complat` 1.1.12 and earlier wheels contain the xtb4stda and stda
binaries without their LGPL-3.0 notice (compliance audit item B2). Yanking
them, or publishing a corrected release with a notice, is a maintainer
decision; `scripts/check_distribution.py`, run by packaging.yml on every
pull request and by release.yml before anything is published, keeps it
from happening again.
