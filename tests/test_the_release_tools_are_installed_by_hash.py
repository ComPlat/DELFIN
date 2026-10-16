"""The tools that build a release are installed by hash, not by version alone.

Input: .github/workflows/release.yml, .github/workflows/packaging.yml, the
hashed requirement files .github/ci/*-requirements.txt, the version pins
in .github/ci/constraints-ci.txt, and .github/ci/hashed_requirements.py,
the script that writes the hashed files.

Output: failures for
- a build tool (pip, build, setuptools, wheel) installed in release.yml's
  ``build`` job or packaging.yml's ``distribution`` job other than by
  ``pip install --require-hashes -r .github/ci/release-requirements.txt``;
  with ``-c constraints-ci.txt`` a re-uploaded or substituted file of the
  same version is installed without notice (OpenSSF Scorecard
  Pinned-Dependencies, SLSA build dependencies);
- a requirement in a hashed file that is not ``name==version`` or carries
  no ``--hash=sha256:`` of 64 hex digits (``--require-hashes`` would
  refuse the whole install, but only on the runner);
- a version that differs from the one constraints-ci.txt names for the
  same distribution: one tool, one version across CI;
- a generator that takes the hash of one file only, or resolves without
  the 7-day cooldown.
"""

from __future__ import annotations

import importlib.util
import json
import re
import subprocess
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
CI = ROOT / ".github" / "ci"
WORKFLOWS = ROOT / ".github" / "workflows"
RELEASE_REQS = CI / "release-requirements.txt"
_HASH = re.compile(r"--hash=sha256:([0-9a-f]{64})\b")


def _job_runs(workflow: str, job: str) -> list[str]:
    wf = yaml.safe_load((WORKFLOWS / workflow).read_text(encoding="utf-8"))
    return [s["run"] for s in wf["jobs"][job]["steps"] if "run" in s]


def _commands(run: str) -> list[str]:
    return [c.strip() for c in run.replace("\\\n", " ").splitlines() if c.strip()]


def parse_hashed(path: Path) -> dict[str, tuple[str, list[str]]]:
    """{name: (version, [sha256, ...])} of a pip requirements file with hashes."""
    text = re.sub(r"\\\n", " ", path.read_text(encoding="utf-8"))
    pins: dict[str, tuple[str, list[str]]] = {}
    for line in text.splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        spec = line.split()[0]
        name, sep, version = spec.partition("==")
        assert sep and re.fullmatch(r"[0-9][0-9A-Za-z.+!-]*", version), (
            f"{path.name}: {spec} is not name==version")
        hashes = _HASH.findall(line)
        assert hashes, f"{path.name}: {spec} has no sha256 hash"
        assert len(hashes) == line.count("--hash="), f"{path.name}: {spec}: a malformed hash"
        pins[name.lower().replace("_", "-")] = (version, hashes)
    return pins


def _constraints() -> dict[str, str]:
    pins = {}
    for line in (CI / "constraints-ci.txt").read_text(encoding="utf-8").splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            name, _, version = line.partition("==")
            pins[name.lower().replace("_", "-")] = version
    return pins


@pytest.mark.parametrize("workflow,job", [("release.yml", "build"),
                                          ("packaging.yml", "distribution")])
def test_the_build_tools_are_installed_from_the_hashed_file(workflow, job):
    cmds = [c for run in _job_runs(workflow, job) for c in _commands(run)
            if re.search(r"\bpip install\b", c)]
    hashed = [c for c in cmds
              if "--require-hashes" in c and "-r .github/ci/release-requirements.txt" in c]
    assert hashed, f"{workflow} {job}: the build tools are not installed by hash: {cmds}"
    for c in cmds:
        if c in hashed:
            continue
        assert not re.search(r"\b(build|setuptools|wheel|--upgrade pip)\b", c), (
            f"{workflow} {job}: {c!r} installs a build tool by version only")


def test_the_release_file_holds_the_build_tools_and_their_closure():
    pins = parse_hashed(RELEASE_REQS)
    for name in ("pip", "build", "setuptools", "wheel", "packaging", "pyproject-hooks"):
        assert name in pins, f"{name} is installed for the release and not hashed"


def test_every_hashed_file_is_pinned_and_agrees_with_the_constraints():
    files = sorted(CI.glob("*-requirements.txt"))
    assert RELEASE_REQS in files
    constraints = _constraints()
    for path in files:
        for name, (version, _) in parse_hashed(path).items():
            if name in constraints:
                assert constraints[name] == version, (
                    f"{path.name}: {name}=={version}, constraints-ci.txt names "
                    f"{constraints[name]}; regenerate one of them")


def _generator():
    spec = importlib.util.spec_from_file_location(
        "hashed_requirements", CI / "hashed_requirements.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_the_generator_hashes_every_file_of_a_release():
    gen = _generator()
    payload = {"urls": [{"digests": {"sha256": "b" * 64}},
                        {"digests": {"sha256": "a" * 64}}]}
    asked = []

    def fetch(url):
        asked.append(url)
        return payload

    assert gen.release_hashes("wheel", "0.48.0", fetch=fetch) == ["a" * 64, "b" * 64]
    assert asked == ["https://pypi.org/pypi/wheel/0.48.0/json"]
    with pytest.raises(SystemExit):
        gen.release_hashes("wheel", "0.48.0", fetch=lambda url: {"urls": []})


def test_the_generator_output_reads_back(tmp_path):
    gen = _generator()
    out = tmp_path / "x-requirements.txt"
    out.write_text(gen.render([("build", "1.6.1", ["a" * 64, "b" * 64]),
                               ("wheel", "0.48.0", ["c" * 64])], "cmd"),
                   encoding="utf-8")
    assert parse_hashed(out) == {"build": ("1.6.1", ["a" * 64, "b" * 64]),
                                 "wheel": ("0.48.0", ["c" * 64])}


def test_the_generator_resolves_with_the_cooldown(monkeypatch):
    gen = _generator()
    seen = {}

    def run(argv, check):
        seen["argv"] = argv
        report = Path(argv[argv.index("--report") + 1])
        report.write_text(json.dumps({"install": [
            {"metadata": {"name": "Pyproject_Hooks", "version": "1.3.3"}}]}),
            encoding="utf-8")
        return subprocess.CompletedProcess(argv, 0)

    monkeypatch.setattr(gen.subprocess, "run", run)
    assert gen.resolve(["build==1.6.1"]) == [("pyproject-hooks", "1.3.3")]
    argv = seen["argv"]
    assert argv[argv.index("--uploaded-prior-to") + 1] == "P7D"
    assert "--dry-run" in argv and "--ignore-installed" in argv
