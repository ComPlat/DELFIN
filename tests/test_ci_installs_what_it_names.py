"""CI installs named versions of its tools and checks what it downloads.

Input: the ``run:`` blocks of .github/workflows/*.yml,
.github/ci/constraints-ci.txt, .github/ci/audit_closure.py and
.github/ci/audit-ignore.txt.

Output: failures for
- a constraint that is not ``name==version``;
- a ``pip install`` that does not install the project (``-e``, ``.``,
  ``.[extras]``) and does not pass the constraints file: the tool version
  would be whatever PyPI served that morning (OpenSSF Scorecard,
  Pinned-Dependencies);
- a download piped straight into ``tar``, ``sh``, ``bash`` or ``python``:
  nothing can check bytes that never touch the disk;
- a dependency audit that cannot fail: ``pip-audit`` followed by
  ``|| true``, or living in ci.yml, whose conclusion is the merge verdict;
- an audit closure that is empty, which pip-audit reports as clean.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
WORKFLOWS = ROOT / ".github" / "workflows"
CONSTRAINTS = ROOT / ".github" / "ci" / "constraints-ci.txt"
_FILES = sorted(WORKFLOWS.glob("*.yml"))


def _runs(path: Path) -> list[str]:
    wf = yaml.safe_load(path.read_text(encoding="utf-8"))
    out = []
    for job in (wf.get("jobs") or {}).values():
        for step in job.get("steps", []):
            if "run" in step:
                out.append(step["run"])
    return out


def _commands(run: str) -> list[str]:
    """Shell lines with backslash continuations joined."""
    return [c.strip() for c in run.replace("\\\n", " ").splitlines() if c.strip()]


def _constraints() -> dict[str, str]:
    pins = {}
    for line in CONSTRAINTS.read_text(encoding="utf-8").splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            name, sep, version = line.partition("==")
            assert sep and version and re.fullmatch(r"[0-9][0-9A-Za-z.+!-]*", version), line
            pins[name.lower().replace("_", "-")] = version
    return pins


def test_every_constraint_is_an_exact_version():
    pins = _constraints()
    for tool in ("pip", "pytest", "pytest-timeout", "ruff", "pip-audit", "build", "coverage"):
        assert tool in pins, f"{tool} is installed by CI and not pinned"


def _is_project_install(cmd: str) -> bool:
    args = cmd.split()
    return any(a in (".", "-e") or a.startswith(('".[', "'.[", ".[")) for a in args)


@pytest.mark.parametrize("path", _FILES, ids=[p.name for p in _FILES])
def test_a_tool_install_names_its_version(path):
    bad = []
    for run in _runs(path):
        for cmd in _commands(run):
            if not re.search(r"\bpip install\b", cmd) or "--dry-run" in cmd:
                continue
            if _is_project_install(cmd):
                continue
            if "-c .github/ci/constraints-ci.txt" not in cmd:
                bad.append(cmd)
    assert not bad, f"{path.name}: unpinned tool installs: {bad}"


@pytest.mark.parametrize("path", _FILES, ids=[p.name for p in _FILES])
def test_no_download_is_piped_into_a_program(path):
    pipe = re.compile(r"\b(curl|wget)\b[^|]*\|\s*(sudo\s+)?(tar|sh|bash|python3?)\b")
    bad = [cmd for run in _runs(path) for cmd in _commands(run) if pipe.search(cmd)]
    assert not bad, f"{path.name}: {bad}"


def test_the_audit_can_fail_and_is_not_in_the_merge_verdict():
    for path in _FILES:
        for cmd in (c for run in _runs(path) for c in _commands(run)):
            assert not re.search(r"pip-audit.*\|\|\s*true", cmd), (path.name, cmd)
    ci = "\n".join(_runs(WORKFLOWS / "ci.yml"))
    assert "pip-audit" not in ci, "an audit in ci.yml makes every PR red on upstream news"

    audit = yaml.safe_load((WORKFLOWS / "audit.yml").read_text(encoding="utf-8"))
    on = audit[True]
    assert "schedule" in on and "workflow_dispatch" in on
    for job in audit["jobs"].values():
        assert job.get("continue-on-error") not in (True, "true")
    run = "\n".join(_runs(WORKFLOWS / "audit.yml"))
    assert 'exit "$code"' in run, "the audit's exit status must be the step's"
    assert ".github/ci/audit-ignore.txt" in run


def _closure(tmp_path: Path, report: dict, *args: str) -> subprocess.CompletedProcess:
    path = tmp_path / "report.json"
    path.write_text(json.dumps(report), encoding="utf-8")
    return subprocess.run(
        [sys.executable, str(ROOT / ".github" / "ci" / "audit_closure.py"), str(path), *args],
        capture_output=True, text=True, timeout=60)


def _item(name: str, version: str, url: str = "https://files.example/x.whl") -> dict:
    return {"metadata": {"name": name, "version": version}, "download_info": {"url": url}}


def test_an_empty_closure_is_refused(tmp_path):
    done = _closure(tmp_path, {"install": []}, "--minimum", "1")
    assert done.returncode == 1 and done.stdout == ""


def test_the_closure_pins_what_pip_would_download_and_leaves_the_project_out(tmp_path):
    report = {"install": [_item("numpy", "1.26.4"), _item("RDKit", "2026.9.1"),
                          _item("delfin-complat", "1.3.3", "file:///work/DELFIN")]}
    done = _closure(tmp_path, report, "--minimum", "2", "--expect", "numpy", "--expect", "rdkit")
    assert done.returncode == 0, done.stderr
    assert done.stdout.splitlines() == ["numpy==1.26.4", "RDKit==2026.9.1"]

    missing = _closure(tmp_path, report, "--expect", "scipy")
    assert missing.returncode == 1 and "scipy" in missing.stderr
