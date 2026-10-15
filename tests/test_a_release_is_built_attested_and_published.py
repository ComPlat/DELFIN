"""The release workflow publishes attested files with the least it needs.

Input: .github/workflows/release.yml and .github/ci/check_dists.py.

Output: failures for
- a trigger other than a v* tag push and a dispatch;
- a workflow-level permission (every job names its own), an OIDC token
  (id-token: write) in any job but `provenance` and `pypi`, a write
  permission in `build`, or `pypi` holding anything but id-token;
- a publish outside the `pypi` environment, or one not limited to v* tags;
- a pip cache in any job of the release, or an asset replaced with
  --clobber;
- check_dists.py passing an ELF, Mach-O or PE image under any name, a file
  under qm_tools/bin/ or qm_tools/share/xtb4stda/, or a dist set that is
  not one wheel and one sdist; or refusing a clean pair.

Measured 2026-10-15: `python -m build` on a clean export of this tree
gives a wheel of 979 files and an sdist of 2834, both passed; with
delfin/qm_tools/bin/xtb4stda present (as on a machine where the installer
ran) both carry it and the check fails, as for the 1.1.12 wheel on PyPI.
"""

from __future__ import annotations

import importlib.util
import io
import tarfile
import zipfile
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
RELEASE = ROOT / ".github" / "workflows" / "release.yml"


def _wf() -> dict:
    return yaml.safe_load(RELEASE.read_text(encoding="utf-8"))


def _check_dists():
    spec = importlib.util.spec_from_file_location("check_dists", ROOT / ".github" / "ci" / "check_dists.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


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
    assert "python -m build --sdist --wheel" in runs
    assert "check_dists.py dist/*" in runs
    assert '[ "$TAG" = "v$version" ]' in runs
    for job in jobs.values():
        for step in job["steps"]:
            assert "cache" not in (step.get("with") or {}), step
            assert "--clobber" not in step.get("run", "")


# --------------------------------------------------------------- check_dists

def _elf() -> bytes:
    return b"\x7fELF\x02\x01\x01" + b"\0" * 64


def _pe() -> bytes:
    head = bytearray(b"MZ" + b"\0" * 0x7E)
    head[0x3C:0x40] = (0x80).to_bytes(4, "little")
    return bytes(head) + b"PE\0\0" + b"\0" * 32


def _dists(tmp_path: Path, members: dict[str, bytes]) -> list[Path]:
    wheel = tmp_path / "pkg-1.0-py3-none-any.whl"
    with zipfile.ZipFile(wheel, "w") as z:
        for name, data in members.items():
            z.writestr(name, data)
    sdist = tmp_path / "pkg-1.0.tar.gz"
    with tarfile.open(sdist, "w:gz") as t:
        for name, data in members.items():
            info = tarfile.TarInfo("pkg-1.0/" + name)
            info.size = len(data)
            t.addfile(info, io.BytesIO(data))
    return [wheel, sdist]


def test_a_clean_pair_passes(tmp_path):
    problems, summary = _check_dists().check(
        _dists(tmp_path, {"delfin/__init__.py": b"x = 1\n", "delfin/data.json": b"{}"}))
    assert problems == [] and len(summary) == 2


@pytest.mark.parametrize("name, data", [
    ("delfin/qm_tools/bin/xtb4stda", b"#!/bin/sh\n"),
    ("delfin/qm_tools/share/xtb4stda/.xtb4stdarc", b"params\n"),
    ("delfin/renamed.dat", _elf()),
    ("delfin/helper.py", bytes.fromhex("cffaedfe") + b"\0" * 28),
    ("delfin/Starter.txt", _pe()),
])
def test_a_program_binary_is_refused_in_both_distributions(tmp_path, name, data):
    problems, _ = _check_dists().check(_dists(tmp_path, {"delfin/__init__.py": b"", name: data}))
    assert len(problems) >= 2, problems
    assert any(p.startswith("pkg-1.0-py3-none-any.whl") for p in problems)
    assert any(p.startswith("pkg-1.0.tar.gz") for p in problems)


def test_a_missing_distribution_is_refused(tmp_path):
    wheel, _ = _dists(tmp_path, {"delfin/__init__.py": b""})
    problems, _ = _check_dists().check([wheel])
    assert problems and "one wheel and one sdist" in problems[0]


def test_text_that_starts_like_a_dos_header_is_not_an_image():
    assert _check_dists().native_image(b"MZ is a postcode prefix" + b" " * 100) is None
