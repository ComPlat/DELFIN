"""The release workflow publishes attested files with the least it needs.

Input: .github/workflows/release.yml.

Output: failures for
- a trigger other than a v* tag push and a dispatch;
- a workflow-level permission (every job names its own), an OIDC token
  (id-token: write) in any job but `provenance` and `pypi`, a write
  permission in `build`, or `pypi` holding anything but id-token;
- a publish outside the `pypi` environment, or one not limited to v* tags;
- a pip cache in any job of the release, or an asset replaced with
  --clobber;
- a build job that does not run scripts/check_distribution.py --release
  over dist/*, or a second distribution checker beside it.

The checker itself (installer output, program names, ELF/Mach-O/PE images
under any name, licence files, the one-wheel-one-sdist release set) is
tested in tests/test_a_release_ships_no_third_party_binary.py.

Measured 2026-10-15: `python -m build` on a clean export of this tree
gives a wheel of 979 files and an sdist of 2834, both passed; with
delfin/qm_tools/bin/xtb4stda present (as on a machine where the installer
ran) both carry it and the check fails, as for the 1.1.12 wheel on PyPI.
"""

from __future__ import annotations

from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
RELEASE = ROOT / ".github" / "workflows" / "release.yml"


def _wf() -> dict:
    return yaml.safe_load(RELEASE.read_text(encoding="utf-8"))


def test_it_runs_on_a_version_tag_and_on_demand_only():
    on = _wf()[True]
    assert set(on) == {"push", "workflow_dispatch"}
    assert on["push"] == {"tags": ["v*"]}
    assert _wf()["concurrency"]["cancel-in-progress"] is False


def test_each_job_holds_only_what_it_uses():
    wf = _wf()
    assert wf["permissions"] == {}
    jobs = wf["jobs"]
    for name, job in jobs.items():
        perms = job.get("permissions")
        assert isinstance(perms, dict), f"{name} inherits the workflow's permissions"
        if name not in ("provenance", "pypi"):
            assert "id-token" not in perms, f"{name} can request an OIDC token"
    assert jobs["build"]["permissions"] == {"contents": "read"}
    assert jobs["provenance"]["permissions"] == {
        "contents": "read", "id-token": "write", "attestations": "write"}
    assert jobs["pypi"]["permissions"] == {"id-token": "write"}
    assert jobs["github"]["permissions"] == {"contents": "write"}


def test_publishing_is_gated_on_a_tag_and_the_pypi_environment():
    jobs = _wf()["jobs"]
    for name in ("pypi", "github"):
        assert "startsWith(github.ref, 'refs/tags/v')" in jobs[name]["if"], name
        assert "provenance" in jobs[name]["needs"], name
    assert jobs["pypi"]["environment"]["name"] == "pypi"
    uses = [s.get("uses", "") for s in jobs["pypi"]["steps"]]
    assert any(u.startswith("pypa/gh-action-pypi-publish@") for u in uses), uses
    attest = [s.get("uses", "") for s in jobs["provenance"]["steps"]]
    assert any(u.startswith("actions/attest-build-provenance@") for u in attest), attest


def test_the_build_is_checked_and_nothing_comes_from_a_cache_or_is_replaced():
    jobs = _wf()["jobs"]
    runs = " ".join(s.get("run", "") for s in jobs["build"]["steps"])
    assert "python -m build --no-isolation --sdist --wheel" in runs, (
        "an isolated build fetches an unpinned backend at release time")
    assert "pip install --require-hashes -r .github/ci/release-requirements.txt" in runs
    assert "python scripts/check_distribution.py --release dist/*" in runs
    assert '[ "$TAG" = "v$version" ]' in runs
    for job in jobs.values():
        for step in job["steps"]:
            assert "cache" not in (step.get("with") or {}), step
            assert "--clobber" not in step.get("run", "")


def test_there_is_one_distribution_checker():
    """Two checkers of one distribution drift: each workflow ran its own,
    with different rules, and a rule added to one did not reach the other."""
    assert not (ROOT / ".github" / "ci" / "check_dists.py").exists()
    for path in sorted((ROOT / ".github" / "workflows").glob("*.yml")):
        text = path.read_text(encoding="utf-8")
        assert "check_dists" not in text, path.name
    packaging = (ROOT / ".github" / "workflows" / "packaging.yml").read_text(encoding="utf-8")
    assert "scripts/check_distribution.py" in packaging
