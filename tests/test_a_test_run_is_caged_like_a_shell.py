"""run_tests runs a test file the way bash runs a command.

A test file is code the agent may have written itself. run_tests handed
it the full environment (provider keys included) and ran it outside the
process cage that confines the agent's shell -- in a locked session too,
where bash is held to one folder by bubblewrap (review 2026-09-16).
"""
import json

import pytest

from delfin.agent import api_client as A
from delfin.agent import test_runner as TR
from delfin.agent.api_client import KitToolPermissions


def _perms(ws, **kw):
    return KitToolPermissions(workspace=str(ws), **kw)


def test_the_test_process_gets_the_shells_scrubbed_environment(tmp_path, monkeypatch):
    monkeypatch.setenv("KIT_TOOLBOX_API_KEY", "kit-secret-value")
    monkeypatch.setenv("SOME_SERVICE_TOKEN", "tok")
    monkeypatch.setenv("DELFIN_HARMLESS", "1")
    seen = {}
    monkeypatch.setattr(TR, "run_tests", lambda **kw: seen.update(kw) or {"status": "ok"})
    A._doc_executor._execute_run_tests({"target": "tests"}, _perms(tmp_path))
    env = seen["env"]
    assert "KIT_TOOLBOX_API_KEY" not in env and "SOME_SERVICE_TOKEN" not in env
    assert env["DELFIN_HARMLESS"] == "1"
    assert callable(seen["wrap"])


def test_the_report_directory_stays_writable_inside_the_cage(tmp_path, monkeypatch):
    seen = {}

    def fake(cmd, cwd, perms, mode=None, extra_write=()):
        seen.update(cmd=cmd, extra_write=tuple(extra_write))
        return ["/bin/bash", "-c", cmd]

    monkeypatch.setattr(A, "_bash_isolation_argv", fake)
    out = A._test_run_argv(["python", "-m", "pytest", "a b.py"], _perms(tmp_path), tmp_path)
    assert out == ["/bin/bash", "-c", "python -m pytest 'a b.py'"]
    assert seen["extra_write"] == (tmp_path,)


def test_a_bwrap_wrap_binds_the_extra_directory_after_the_fresh_tmp(tmp_path, monkeypatch):
    monkeypatch.setattr(A, "_bwrap_functional", lambda: True)
    monkeypatch.setattr(A.shutil, "which", lambda _x: "/usr/bin/bwrap")
    perms = KitToolPermissions(workspace=str(tmp_path), lock_workspace=True)
    report = tmp_path / "report"
    out = A._bash_isolation_argv("true", tmp_path, perms, extra_write=(report,))
    rd = str(report.resolve())
    i = out.index(rd)
    assert out[i - 1] == "--bind" and i > out.index("--tmpfs")


@pytest.mark.skipif(not A._bwrap_functional(), reason="bubblewrap does not work here")
@pytest.mark.parametrize("locked", [False, True])
def test_a_real_run_in_the_cage_still_reports_its_result(tmp_path, monkeypatch, locked):
    monkeypatch.setenv("KIT_TOOLBOX_API_KEY", "kit-secret-value")
    ws = tmp_path / "ws"
    ws.mkdir()
    (ws / "test_env.py").write_text(
        "import os\n"
        "def test_no_key():\n"
        "    assert 'KIT_TOOLBOX_API_KEY' not in os.environ\n"
        "def test_second():\n"
        "    assert True\n")
    perms = _perms(ws, lock_workspace=locked)
    out = json.loads(A._doc_executor._execute_run_tests(
        {"target": "test_env.py", "timeout_s": 120}, perms))
    assert out["status"] == "ok", out
    assert out["summary"].get("passed") == 2, out
