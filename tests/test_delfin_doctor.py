"""Tests for delfin/doctor.py: installation self-check.

No real ORCA, no network: binaries are mocked as executable scripts in a
tmp directory prepended to PATH; scratch dirs live in tmp_path; the KIT
key is removed from the environment.
"""

import json
import os
import stat

import pytest

from delfin import doctor
from delfin.doctor import (
    BROKEN,
    MISSING,
    OK,
    CheckResult,
    check_docs_index,
    check_kit_toolbox_key,
    check_openmpi,
    check_orca,
    check_scratch_dir,
    check_slurm,
    check_xtb,
    exit_code,
    run_all,
)

# ---------------------------------------------------------------------------
# Helpers / fixtures
# ---------------------------------------------------------------------------

BINARIES = ["orca", "xtb", "mpirun", "sinfo"]


def make_binary(bindir, name, behavior="ok"):
    """Write an executable script into bindir.

    behavior:
      "ok"     -> prints a version line, exit 0
      "broken" -> exits 1
    """
    path = bindir / name
    if behavior == "ok":
        body = "#!/bin/sh\necho 'mock 1.0'\nexit 0\n"
    else:
        body = "#!/bin/sh\necho 'boom' >&2\nexit 1\n"
    path.write_text(body)
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return path


@pytest.fixture
def clean_env(monkeypatch, tmp_path):
    """PATH with only an empty mock bin dir; no KIT key, no scratch vars.

    ORCA and xtb are found through DELFIN's resolver (qm_health), which
    also looks outside PATH; tests steer it through ``tool_health``."""
    bindir = tmp_path / "bin"
    bindir.mkdir()
    monkeypatch.setenv("PATH", str(bindir))
    monkeypatch.delenv("KIT_TOOLBOX_API_KEY", raising=False)
    monkeypatch.delenv("DELFIN_SCRATCH", raising=False)
    import delfin.agent.credentials as creds
    monkeypatch.setattr(creds, "load_credential",
                        lambda name, **kw: os.environ.get(name, ""))
    monkeypatch.setattr(doctor, "_tool_health", lambda name: _health(name, "absent"))
    return bindir


def _health(name, level, path="/opt/x/bin/" , why=""):
    from delfin.qm_health import ToolHealth
    return ToolHealth(name=name, level=level, present=level != "absent",
                      path=(path + name) if level != "absent" else "",
                      why=why or ("not found" if level == "absent" else ""))


@pytest.fixture
def tool_health(monkeypatch):
    """Set what the resolver reports for orca/xtb: tool_health("orca", "fail")."""
    levels = {}
    monkeypatch.setattr(doctor, "_tool_health",
                        lambda name: _health(name, levels.get(name, "absent")))

    def set_level(name, level):
        levels[name] = level
    return set_level


@pytest.fixture
def scratch_dir(tmp_path):
    d = tmp_path / "scratch"
    d.mkdir()
    return d


@pytest.fixture
def docs_index_file(monkeypatch, tmp_path):
    """Redirect get_default_index_path to a file inside tmp_path.

    Returns a factory: call with a dict to write valid JSON, with a str
    to write that raw text (invalid JSON), or with None to remove the file.
    """
    import delfin.doc_server.indexer as indexer

    index_file = tmp_path / "doc_index.json"

    def set_default_index_path():
        return index_file

    monkeypatch.setattr(indexer, "get_default_index_path", set_default_index_path)

    def write(payload):
        if payload is None:
            index_file.unlink(missing_ok=True)
        elif isinstance(payload, str):
            index_file.write_text(payload)
        else:
            index_file.write_text(json.dumps(payload))
        return index_file

    return write


# ---------------------------------------------------------------------------
# Binary checks: ok / missing / broken
# ---------------------------------------------------------------------------

BINARY_CHECKS = [
    ("mpirun", check_openmpi),
    ("sinfo", check_slurm),
]

QM_TOOL_CHECKS = [
    ("orca", check_orca),
    ("xtb", check_xtb),
]


@pytest.mark.parametrize("tool,func", QM_TOOL_CHECKS, ids=[t for t, _ in QM_TOOL_CHECKS])
class TestQmToolChecks:
    """ORCA and xtb go through DELFIN's resolver, not a --version probe:
    ORCA has no --version and exits 2, which reported every working ORCA
    as broken in the first run."""

    def test_ok(self, clean_env, tool_health, tool, func):
        tool_health(tool, "ok")
        result = func()
        assert result.status == OK and tool in result.detail

    def test_missing(self, clean_env, tool_health, tool, func):
        assert func().status == MISSING

    def test_broken(self, clean_env, tool_health, tool, func):
        tool_health(tool, "fail")
        assert func().status == BROKEN


def test_orca_is_not_probed_with_version(monkeypatch, clean_env):
    """A mock `orca` on PATH that fails on --version must not make ORCA broken:
    the resolver, not --version, decides."""
    make_binary(clean_env, "orca", "broken")
    monkeypatch.setattr(doctor, "_tool_health", lambda name: _health(name, "ok"))
    assert check_orca().status == OK


@pytest.mark.parametrize("binary,func", BINARY_CHECKS, ids=[b for b, _ in BINARY_CHECKS])
class TestBinaryChecks:
    def test_ok(self, clean_env, binary, func):
        make_binary(clean_env, binary, "ok")
        result = func()
        assert isinstance(result, CheckResult)
        assert result.status == OK

    def test_missing(self, clean_env, binary, func):
        result = func()
        assert result.status == MISSING

    def test_broken(self, clean_env, binary, func):
        make_binary(clean_env, binary, "broken")
        result = func()
        assert result.status == BROKEN

    def test_result_has_name(self, clean_env, binary, func):
        result = func()
        assert result.name


# ---------------------------------------------------------------------------
# Scratch dir check
# ---------------------------------------------------------------------------

class TestScratchDir:
    def test_ok_explicit_arg(self, clean_env, scratch_dir):
        assert check_scratch_dir(scratch_dir).status == OK

    def test_ok_via_env(self, clean_env, scratch_dir, monkeypatch):
        monkeypatch.setenv("DELFIN_SCRATCH", str(scratch_dir))
        assert check_scratch_dir().status == OK

    def test_ok_fallback_tmpdir(self, clean_env, tmp_path, monkeypatch):
        monkeypatch.setenv("TMPDIR", str(tmp_path))
        assert check_scratch_dir().status == OK

    def test_broken_readonly(self, clean_env, scratch_dir):
        scratch_dir.chmod(0o500)
        if os.access(scratch_dir, os.W_OK):
            pytest.skip("running as root: chmod does not revoke write access")
        try:
            assert check_scratch_dir(scratch_dir).status == BROKEN
        finally:
            scratch_dir.chmod(0o700)

    def test_missing_nonexistent_dir(self, clean_env, tmp_path):
        assert check_scratch_dir(tmp_path / "does-not-exist").status == MISSING

    def test_leaves_no_probe_file(self, clean_env, scratch_dir):
        check_scratch_dir(scratch_dir)
        assert not (scratch_dir / ".delfin_doctor_probe").exists()


# ---------------------------------------------------------------------------
# KIT key check
# ---------------------------------------------------------------------------

class TestKitKey:
    def test_ok(self, clean_env, monkeypatch):
        monkeypatch.setenv("KIT_TOOLBOX_API_KEY", "dummy")
        assert check_kit_toolbox_key().status == OK

    def test_missing(self, clean_env):
        assert check_kit_toolbox_key().status == MISSING

    def test_a_key_in_the_credential_store_counts(self, clean_env, monkeypatch):
        import delfin.agent.credentials as creds
        monkeypatch.setattr(creds, "load_credential",
                            lambda name, **kw: "stored" if name == "KIT_TOOLBOX_API_KEY" else "")
        result = check_kit_toolbox_key()
        assert result.status == OK and "stored" not in result.detail


# ---------------------------------------------------------------------------
# Docs index check
# ---------------------------------------------------------------------------

class TestDocsIndex:
    def test_ok(self, clean_env, docs_index_file):
        docs_index_file({"documents": []})
        assert check_docs_index().status == OK

    def test_missing(self, clean_env, docs_index_file):
        docs_index_file(None)
        assert check_docs_index().status == MISSING

    def test_broken_invalid_json(self, clean_env, docs_index_file):
        docs_index_file("{not json!")
        assert check_docs_index().status == BROKEN


# ---------------------------------------------------------------------------
# run_all and exit_code
# ---------------------------------------------------------------------------

class TestRunAll:
    def test_returns_every_check_in_order(self, clean_env, scratch_dir, docs_index_file):
        docs_index_file(None)
        results = run_all(scratch_dir=scratch_dir)
        names = [r.name for r in results]
        assert len(results) == 9
        assert len(set(names)) == 9  # distinct, no duplicates
        assert all(isinstance(r, CheckResult) for r in results)

    def test_all_ok(self, clean_env, scratch_dir, docs_index_file, monkeypatch, tool_health):
        for b in ("mpirun", "sinfo"):
            make_binary(clean_env, b, "ok")
        tool_health("orca", "ok")
        tool_health("xtb", "ok")
        # Stored, not exported: an exported key is itself a finding.
        monkeypatch.delenv("KIT_TOOLBOX_API_KEY", raising=False)
        monkeypatch.setattr("delfin.agent.credentials.load_credential",
                            lambda name, **kw: "dummy")
        monkeypatch.setattr("delfin.doctor.check_command_isolation",
                            lambda: CheckResult("command_isolation", OK, "supplied"))
        docs_index_file({"documents": []})
        results = run_all(scratch_dir=scratch_dir)
        assert all(r.status == OK for r in results), [
            (r.name, r.status, r.detail) for r in results
        ]

    def test_never_raises_on_crash(self, clean_env, monkeypatch, scratch_dir):
        def boom():
            raise RuntimeError("unexpected")

        monkeypatch.setattr("delfin.doctor.check_orca", boom)
        results = run_all(scratch_dir=scratch_dir)
        crashed = [r for r in results if r.status == BROKEN]
        assert crashed, "a crashing check must surface as broken"


class TestExitCode:
    def test_zero_when_all_ok(self):
        results = [
            CheckResult("a", OK),
            CheckResult("b", OK),
        ]
        assert exit_code(results) == 0

    def test_zero_when_only_missing(self):
        results = [
            CheckResult("a", OK),
            CheckResult("b", MISSING),
        ]
        assert exit_code(results) == 0

    def test_one_when_broken(self):
        results = [
            CheckResult("a", OK),
            CheckResult("b", BROKEN),
        ]
        assert exit_code(results) == 1

    def test_one_beats_missing(self):
        results = [
            CheckResult("a", MISSING),
            CheckResult("b", BROKEN),
        ]
        assert exit_code(results) == 1

    def test_end_to_end_broken_binary(self, clean_env, scratch_dir, docs_index_file):
        docs_index_file(None)
        make_binary(clean_env, "mpirun", "broken")
        make_binary(clean_env, "sinfo", "ok")
        results = run_all(scratch_dir=scratch_dir)
        assert exit_code(results) == 1

    def test_end_to_end_missing_only(self, clean_env, scratch_dir, docs_index_file):
        docs_index_file(None)
        results = run_all(scratch_dir=scratch_dir)
        assert exit_code(results) == 0
        # no binary, no key, no index: the only "ok" is the writable scratch dir
        assert all(r.status in (MISSING, OK) for r in results)
        assert not any(r.status == BROKEN for r in results)


# ---------------------------------------------------------------------------
# CheckResult dataclass shape
# ---------------------------------------------------------------------------

class TestCheckResult:
    def test_defaults(self):
        r = CheckResult("x", "ok")
        assert r.detail == ""
        assert r.fix_hint == ""

    def test_status_constants(self):
        assert (OK, MISSING, BROKEN) == ("ok", "missing", "broken")
