"""Security: calc/archive/delfin are REACHABLE for reading but not freely
destroyable (user 2026-06-26: "er kann jetzt nicht einfach calc zerstören ohne
nachzufragen … archive sind fix … in agent_workspace kann er sich austoben").

Tiers:
  agent_workspace + launch workspace → writable (no prompt)
  calc                               → readable; WRITE needs an explicit confirm
  archive, delfin checkout           → readable; WRITE hard-denied
Defense in depth is double-checked: the write EXECUTOR resolves writable-only,
so it can never touch a read-only dir even if the gate were bypassed.
"""
from pathlib import Path
import json
import pytest
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor


@pytest.fixture
def env(tmp_path):
    ws = tmp_path / "proj"; ws.mkdir()
    aw = tmp_path / "agent_workspace"; aw.mkdir()
    calc = tmp_path / "calc"; calc.mkdir()
    arch = tmp_path / "archive"; arch.mkdir()
    delfin = tmp_path / "software" / "delfin"; delfin.mkdir(parents=True)
    confirms = []
    perms = KitToolPermissions(
        workspace=ws, mode="acceptEdits",
        extra_workspace_dirs=(aw, calc),
        confirm_write_dirs=(calc,),
        read_only_workspace_dirs=(arch, delfin),
        confirm_callback=lambda n, a, p: (confirms.append(p) or True),
    )
    return _DocToolExecutor(), perms, locals()


def _w(ex, perms, p):
    return ex._run_permission_gate("write_file", {"path": str(p)}, perms)


def test_write_to_archive_is_hard_denied(env):
    ex, perms, L = env
    err = _w(ex, perms, L["arch"] / "result.out")
    assert err and "READ-ONLY" in err and "copy" in err.lower()
    assert L["confirms"] == []                      # never even asked


def test_write_to_delfin_checkout_is_hard_denied(env):
    ex, perms, L = env
    err = _w(ex, perms, L["delfin"] / "delfin" / "x.py")
    assert err and "READ-ONLY" in err


def test_write_to_calc_requires_confirm(env):
    ex, perms, L = env
    err = _w(ex, perms, L["calc"] / "job.inp")
    assert err is None                              # allowed AFTER confirm
    assert len(L["confirms"]) == 1
    assert "CALC EDIT" in L["confirms"][0]


def test_calc_write_denied_when_user_says_no(tmp_path):
    ws = tmp_path / "p"; ws.mkdir(); calc = tmp_path / "calc"; calc.mkdir()
    perms = KitToolPermissions(
        workspace=ws, mode="acceptEdits", extra_workspace_dirs=(calc,),
        confirm_write_dirs=(calc,), confirm_callback=lambda *a: False)
    err = _DocToolExecutor()._run_permission_gate(
        "write_file", {"path": str(calc / "j.inp")}, perms)
    assert err and "denied" in err.lower()


def test_write_to_agent_workspace_is_free(env):
    ex, perms, L = env
    assert _w(ex, perms, L["aw"] / "app.py") is None
    assert _w(ex, perms, L["ws"] / "main.py") is None     # primary workspace
    assert L["confirms"] == []                            # no prompt


def test_reads_reach_calc_archive_delfin_without_prompt(env):
    ex, perms, L = env
    for d in (L["calc"], L["arch"], L["delfin"]):
        f = d / "data.txt"; f.write_text("hi")
        out = ex._execute_read_file({"path": str(f)}, perms)
        assert "hi" in out and "error" not in out.lower()
    assert L["confirms"] == []


def test_defense_in_depth_executor_refuses_readonly_write(env):
    """Even bypassing the gate, the write EXECUTOR resolves writable-only, so a
    read-only path can't be written."""
    ex, perms, L = env
    out = ex._execute_write_file(
        {"path": str(L["arch"] / "x.out"), "content": "evil"}, perms)
    assert "error" in out.lower() and "outside" in out.lower()
    assert not (L["arch"] / "x.out").exists()


def test_secret_deny_still_hard_blocks_in_reachable_dir(env):
    ex, perms, L = env
    (L["arch"] / ".env").write_text("SECRET=1")
    out = ex._execute_read_file({"path": str(L["arch"] / ".env")}, perms)
    assert "error" in out.lower() and "deny" in out.lower()
