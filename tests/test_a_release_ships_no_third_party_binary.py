"""A built distribution carries the licence texts and no third-party program.

Input: synthetic wheels and sdists written into ``tmp_path``, and the
repository's ``pyproject.toml``, ``MANIFEST.in`` and
``.github/workflows/packaging.yml``. Output: assertions on
``scripts/check_distribution.py`` and on the packaging configuration.

The published 1.1.8 to 1.1.12 wheels and sdists carry
``delfin/qm_tools/bin/xtb4stda`` (25.6 MB) and ``stda_v1.6.1`` (43.4 MB),
LGPL-3.0 programs of grimme-lab/xtb4stda, plus their parameter files, and
no licence text or source for them; the package-data entries that put them
there selected whatever ``install_qm_tools.sh`` had downloaded on the
building machine. The checker is driven here on archives built in the
test, so the suite stays fast; the real build runs in packaging.yml.
"""

from __future__ import annotations

import importlib.util
import io
import re
import sys
import tarfile
import zipfile
from pathlib import Path

import pytest
import tomllib

REPO = Path(__file__).resolve().parents[1]

_spec = importlib.util.spec_from_file_location(
    "check_distribution", REPO / "scripts" / "check_distribution.py")
CD = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = CD      # dataclasses resolve annotations through it
_spec.loader.exec_module(CD)

ELF = b"\x7fELF\x02\x01\x01" + b"\x00" * 64
DIST = "delfin_complat-9.9.9"


def _metadata(expression: str = CD.LICENSE_EXPRESSION,
              files: tuple[str, ...] = CD.REQUIRED_LICENCE_FILES) -> str:
    lines = ["Metadata-Version: 2.4", "Name: delfin-complat", "Version: 9.9.9",
             f"License-Expression: {expression}"]
    lines += [f"License-File: {f}" for f in files]
    return "\n".join(lines) + "\n\n"


def _wheel(tmp_path: Path, extra: dict[str, bytes] | None = None, *,
           licences: tuple[str, ...] = CD.REQUIRED_LICENCE_FILES,
           metadata: str | None = None) -> Path:
    path = tmp_path / f"{DIST}-py3-none-any.whl"
    with zipfile.ZipFile(path, "w") as zf:
        zf.writestr("delfin/__init__.py", b"")
        zf.writestr("delfin/qm_tools/install_qm_tools.sh", b"#!/usr/bin/env bash\n")
        zf.writestr(f"{DIST}.dist-info/METADATA", metadata or _metadata())
        for name in licences:
            zf.writestr(f"{DIST}.dist-info/licenses/{name}", b"text")
        for name, data in (extra or {}).items():
            zf.writestr(name, data)
    return path


def _sdist(tmp_path: Path, extra: dict[str, bytes] | None = None, *,
           licences: tuple[str, ...] = CD.REQUIRED_LICENCE_FILES) -> Path:
    path = tmp_path / f"{DIST}.tar.gz"
    files = {"delfin/__init__.py": b"", "PKG-INFO": _metadata().encode()}
    files.update({name: b"text" for name in licences})
    files.update(extra or {})
    with tarfile.open(path, "w:gz") as tf:
        for name, data in files.items():
            info = tarfile.TarInfo(f"{DIST}/{name}")
            info.size = len(data)
            tf.addfile(info, io.BytesIO(data))
    return path


# ---------------------------------------------------------------------------
# The checker, driven both ways
# ---------------------------------------------------------------------------

def test_a_clean_wheel_and_sdist_pass(tmp_path):
    assert CD.check(_wheel(tmp_path)) == []
    assert CD.check(_sdist(tmp_path)) == []


@pytest.mark.parametrize("member", [
    "delfin/qm_tools/bin/xtb4stda",
    "delfin/qm_tools/bin/stda_v1.6.1",
    "delfin/qm_tools/share/xtb4stda/.param_stda1.xtb",
    "delfin/csp_tools/bin/gnrs",
    "delfin/qm_tools/third_party/orca/orca",
])
def test_an_installer_output_is_refused_in_both_formats(tmp_path, member):
    for build in (_wheel, _sdist):
        findings = CD.check(build(tmp_path, {member: b"payload"}))
        assert any(member in f for f in findings), (build.__name__, findings)


def test_a_program_name_is_refused_outside_the_installer_directories(tmp_path):
    findings = CD.check(_wheel(tmp_path, {"delfin/tools/xtb4stda": b"#!/bin/sh\n"}))
    assert any("third-party program" in f for f in findings)


@pytest.mark.parametrize("member", [
    "delfin/somewhere/helper",              # no suffix, ELF header
    "delfin/somewhere/helper.dat",
])
def test_native_code_is_refused_wherever_it_sits(tmp_path, member):
    findings = CD.check(_wheel(tmp_path, {member: ELF}))
    assert any("native executable code (ELF)" in f for f in findings), findings


def test_a_native_suffix_is_refused_without_a_header(tmp_path):
    findings = CD.check(_wheel(tmp_path, {"delfin/_speedups.so": b"x"}))
    assert any("native library suffix" in f for f in findings)


def test_the_executable_formats_are_recognised():
    pe = bytearray(b"MZ" + b"\x00" * 62 + b"PE\x00\x00")
    pe[60:64] = (64).to_bytes(4, "little")
    assert CD.native_kind(ELF) == "ELF"
    assert CD.native_kind(b"\xcf\xfa\xed\xfe....") == "Mach-O"
    assert CD.native_kind(bytes(pe)) == "PE"
    # A text file that happens to start with MZ is not a PE image.
    assert CD.native_kind(b"MZ is a postcode" + b" " * 60) == ""
    assert CD.native_kind(b"#!/usr/bin/env bash\n") == ""


@pytest.mark.parametrize("missing", CD.REQUIRED_LICENCE_FILES)
def test_a_missing_licence_file_is_named(tmp_path, missing):
    kept = tuple(f for f in CD.REQUIRED_LICENCE_FILES if f != missing)
    for path in (_wheel(tmp_path, licences=kept), _sdist(tmp_path, licences=kept)):
        findings = CD.check(path)
        assert any(f"/{missing} is missing" in f for f in findings), findings


def test_the_metadata_must_state_the_spdx_expression(tmp_path):
    findings = CD.check(_wheel(tmp_path, metadata=_metadata("LGPL-3.0")))
    assert any("License-Expression is 'LGPL-3.0'" in f for f in findings)
    findings = CD.check(_wheel(tmp_path, metadata=_metadata(files=("LICENSE",))))
    assert any("no License-File entry for COPYING" in f for f in findings)


def _pe() -> bytes:
    head = bytearray(b"MZ" + b"\x00" * 0x7E)
    head[0x3C:0x40] = (0x80).to_bytes(4, "little")
    return bytes(head) + b"PE\x00\x00" + b"\x00" * 32


@pytest.mark.parametrize("member, data, kind", [
    ("delfin/renamed.dat", ELF, "ELF"),
    ("delfin/helper.py", bytes.fromhex("cffaedfe") + b"\x00" * 28, "Mach-O"),
    ("delfin/Starter.txt", _pe(), "PE"),
])
def test_an_executable_image_is_refused_under_any_name_in_both_formats(
        tmp_path, member, data, kind):
    """The header is read, not the name: a renamed program is still one."""
    for build in (_wheel, _sdist):
        findings = CD.check(build(tmp_path, {member: data}))
        assert any(f"native executable code ({kind})" in f for f in findings), (
            build.__name__, findings)


def test_installer_output_outside_the_package_is_refused(tmp_path):
    findings = CD.check(_sdist(tmp_path, {"qm_tools/bin/xtb4stda": b"#!/bin/sh\n"}))
    assert any("installer's output directory" in f for f in findings), findings


def test_the_release_set_is_one_wheel_and_one_sdist(tmp_path, capsys):
    wheel, sdist = _wheel(tmp_path), _sdist(tmp_path)
    assert CD.release_set_findings([wheel, sdist]) == []
    assert CD.main(["--release", str(wheel), str(sdist)]) == 0
    capsys.readouterr()

    missing = CD.release_set_findings([wheel])
    assert missing and "one wheel and one sdist" in missing[0]
    assert CD.main(["--release", str(wheel)]) == 1
    assert "one wheel and one sdist" in capsys.readouterr().out

    second = tmp_path / "other"
    second.mkdir()
    assert CD.release_set_findings([wheel, _wheel(second), sdist])
    stray = tmp_path / "notes.txt"
    stray.write_text("x")
    assert any("not a wheel or an sdist" in f
               for f in CD.release_set_findings([wheel, sdist, stray]))
    # Without --release any number of archives is checked one by one, as
    # packaging.yml does with its two wheels.
    assert CD.main([str(wheel), str(_wheel(second)), str(sdist)]) == 0


def test_the_command_line_exits_one_on_a_finding(tmp_path, capsys):
    clean = _wheel(tmp_path)
    assert CD.main([str(clean)]) == 0
    dirty = _sdist(tmp_path, {"delfin/qm_tools/bin/xtb4stda": ELF})
    assert CD.main([str(dirty)]) == 1
    assert "xtb4stda" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# The configuration that decides what a build selects
# ---------------------------------------------------------------------------

_INSTALLER_DIRS = ("bin", "share", "downloads", "build", "third_party", "logs")


def _pyproject() -> dict:
    return tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))


def test_package_data_selects_no_installer_output():
    data = _pyproject()["tool"]["setuptools"]["package-data"]
    for package, patterns in data.items():
        for pattern in patterns:
            assert not re.search(
                r"_tools/({})/".format("|".join(_INSTALLER_DIRS)), pattern), (
                f"{package}: {pattern} selects files an installer downloads")


def test_exclude_package_data_covers_what_the_qm_installer_writes():
    excluded = _pyproject()["tool"]["setuptools"]["exclude-package-data"]["delfin"]
    for sub in _INSTALLER_DIRS:
        assert f"qm_tools/{sub}/*" in excluded, sub


def test_manifest_prunes_every_installer_directory_after_its_includes():
    """MANIFEST.in applies its lines in order: an include after a prune
    re-adds what the prune removed, which is how 1.1.8 to 1.1.12 shipped
    the binaries although ``prune delfin/qm_tools/bin`` was there."""
    lines = [ln.split() for ln in (REPO / "MANIFEST.in").read_text().splitlines()
             if ln.strip() and not ln.lstrip().startswith("#")]
    prunes = {ln[1]: i for i, ln in enumerate(lines) if ln[0] == "prune"}
    for sub in _INSTALLER_DIRS:
        target = f"delfin/qm_tools/{sub}"
        assert target in prunes, f"MANIFEST.in does not prune {target}"
        for i, ln in enumerate(lines):
            if ln[0] in {"include", "recursive-include", "graft"} and i > prunes[target]:
                assert not any(arg.startswith(target) for arg in ln[1:]), (
                    f"line {' '.join(ln)!r} re-adds {target} after its prune")


def test_the_checker_requires_exactly_the_declared_licence_files():
    """Two copies of one list drift; this pins the pair."""
    declared = tuple(_pyproject()["project"]["license-files"])
    assert set(declared) == set(CD.REQUIRED_LICENCE_FILES)
    assert _pyproject()["project"]["license"] == CD.LICENSE_EXPRESSION


def test_every_required_licence_file_exists_and_is_in_the_sdist():
    manifest = (REPO / "MANIFEST.in").read_text()
    included = {arg for ln in manifest.splitlines()
                if ln.startswith("include ") for arg in ln.split()[1:]}
    for name in CD.REQUIRED_LICENCE_FILES:
        assert (REPO / name).is_file(), name
        assert name in included, f"MANIFEST.in does not include {name}"


def test_the_packaging_workflow_builds_and_runs_the_checker():
    text = (REPO / ".github" / "workflows" / "packaging.yml").read_text()
    assert "permissions:\n  contents: read" in text
    assert "python -m build" in text
    assert "scripts/check_distribution.py" in text
    # GitHub-owned actions only, the repository's own allow-list setting.
    for owner in re.findall(r"uses:\s*([\w.-]+)/", text):
        assert owner == "actions", owner
    # The build tools are pinned, so the same commit builds the same set:
    # installed through the CI constraints file, which names exact versions
    # of the frontend and of the backend --no-isolation builds against.
    assert "pip install -c .github/ci/constraints-ci.txt build setuptools wheel" in text
    pins = (REPO / ".github" / "ci" / "constraints-ci.txt").read_text()
    for tool in ("build", "setuptools", "wheel"):
        assert re.search(rf"^{tool}==\d", pins, re.M), tool
