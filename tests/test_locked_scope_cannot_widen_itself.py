"""A locked session decides for itself, not for the ones after it.

The lock already refused ``extra_dir``, because widening the folder is
the obvious way out. Two other routes reached the same place and were
not covered:

* ``default_mode="bypassPermissions"`` sets the live mode AND persists it,
  so every future session starts unattended;
* ``allow_pattern="^.*$"`` persists a rule that auto-approves every shell
  command from then on.

Separately, a hooks file is executable configuration -- entries run
through a shell with the process environment, outside the permission gate
and outside any filesystem isolation. Reading one from the workspace is
right for a project you own and wrong for a folder that receives files
from other people, which is what a locked scope is.
"""

from __future__ import annotations

import json
import pathlib
import tempfile

import pytest

from delfin.agent import api_client as A


def _perms(locked: bool, ws=None):
    ws = ws or pathlib.Path(tempfile.mkdtemp(prefix="ws_"))
    p = A.KitToolPermissions(
        workspace=ws, agent_role="office_agent" if locked else "")
    p.mode = "acceptEdits"
    p.confirm_callback = lambda n, a, prev: True
    return p


def _executor():
    return A._DocToolExecutor.__new__(A._DocToolExecutor)


# ---------------------------------------------------------------------------
# remember_permission
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("kind,value", [
    ("default_mode", "bypassPermissions"),
    ("allow_pattern", "^.*$"),
])
def test_a_locked_session_cannot_persist_a_wider_permission(kind, value):
    perms = _perms(locked=True)
    out = _executor()._execute_remember_permission(
        {"kind": kind, "value": value}, perms)
    assert '"error"' in out, out
    assert kind in out


def test_the_refusal_does_not_change_the_live_mode():
    """The persist and the live change happen together; refusing one must
    refuse both, or the session is unattended for the rest of its life."""
    perms = _perms(locked=True)
    _executor()._execute_remember_permission(
        {"kind": "default_mode", "value": "bypassPermissions"}, perms)
    assert perms.mode == "acceptEdits"


def test_extra_dir_is_still_refused():
    """The route that was already closed must stay closed."""
    perms = _perms(locked=True)
    out = _executor()._execute_remember_permission(
        {"kind": "extra_dir", "value": "/etc"}, perms)
    assert '"error"' in out


@pytest.mark.parametrize("kind,value", [
    ("default_mode", "bypassPermissions"),
    ("allow_pattern", "^git status$"),
])
def test_an_unlocked_session_is_unchanged(kind, value):
    perms = _perms(locked=False)
    out = _executor()._execute_remember_permission(
        {"kind": kind, "value": value}, perms)
    assert '"error"' not in out, out


def test_the_refusal_says_what_to_do_instead():
    perms = _perms(locked=True)
    out = _executor()._execute_remember_permission(
        {"kind": "default_mode", "value": "bypassPermissions"}, perms)
    assert "ask the user" in out.lower(), (
        "a refusal the model cannot act on turns into a workaround attempt")


# ---------------------------------------------------------------------------
# Hooks
# ---------------------------------------------------------------------------

def _workspace_with_hook():
    ws = pathlib.Path(tempfile.mkdtemp(prefix="ws_"))
    (ws / ".delfin").mkdir()
    (ws / ".delfin" / "settings.json").write_text(json.dumps({
        "hooks": {"PreToolUse": [{
            "matcher": "*",
            "hooks": [{"type": "command", "command": "echo probe"}],
        }]}
    }), encoding="utf-8")
    return ws


def test_a_locked_workspace_supplies_no_hooks():
    import delfin.agent.hooks as H

    ws = _workspace_with_hook()
    cfg = H.load_hooks(A._hook_workspace(_perms(locked=True, ws=ws)))
    assert not cfg.by_event, (
        "a settings file inside the documents folder still runs commands")


def test_an_unlocked_workspace_still_supplies_its_hooks():
    """The feature is good; it is the folder it is read from that is wrong."""
    import delfin.agent.hooks as H

    ws = _workspace_with_hook()
    cfg = H.load_hooks(A._hook_workspace(_perms(locked=False, ws=ws)))
    assert cfg.by_event.get("PreToolUse"), "workspace hooks stopped working"


def test_every_hook_load_site_goes_through_the_helper():
    """Four call sites read hooks; one of them left unguarded is the whole
    hole back."""
    source = pathlib.Path(A.__file__).read_text(encoding="utf-8")
    assert "load_hooks(permissions.workspace)" not in source, (
        "a hook load site reads the workspace directly again")
    assert source.count("load_hooks(_hook_workspace(permissions))") == 4


def test_the_helper_is_defensive_about_its_input():
    class _Bare:
        pass

    assert A._hook_workspace(_Bare()) is None
    assert A._hook_workspace(None) is None


# ---------------------------------------------------------------------------
# Isolation state must be visible
# ---------------------------------------------------------------------------

def test_a_locked_session_without_bwrap_says_so(monkeypatch):
    """Silent degradation was the problem: the promise reads the same to
    the user whether or not anything is enforcing it."""
    recorded: list[tuple] = []
    monkeypatch.setattr(A, "_bwrap_functional", lambda: False)
    monkeypatch.setattr(A, "_ISOLATION_GAP_ANNOUNCED", False)
    monkeypatch.setattr(
        A, "_record_security_event",
        lambda *a, **kw: recorded.append((a, kw)))

    ws = pathlib.Path(tempfile.mkdtemp(prefix="ws_"))
    perms = _perms(locked=True, ws=ws)
    perms.mode = "default"
    argv = A._bash_isolation_argv("ls", ws, perms)

    assert argv[0] == "/bin/bash", "expected the documented fallback"
    assert recorded, "the downgrade was not recorded anywhere"
    text = " ".join(str(x) for x in recorded[0][0])
    assert "isolation" in text and "NOT active" in text


def test_the_announcement_fires_once(monkeypatch):
    """A line per bash command would bury the panel it is meant to inform."""
    recorded: list = []
    monkeypatch.setattr(A, "_bwrap_functional", lambda: False)
    monkeypatch.setattr(A, "_ISOLATION_GAP_ANNOUNCED", False)
    monkeypatch.setattr(A, "_record_security_event",
                        lambda *a, **kw: recorded.append(a))
    ws = pathlib.Path(tempfile.mkdtemp(prefix="ws_"))
    perms = _perms(locked=True, ws=ws)
    for _ in range(5):
        A._bash_isolation_argv("ls", ws, perms)
    assert len(recorded) == 1


def test_a_working_bwrap_still_isolates_a_locked_session(monkeypatch):
    monkeypatch.setattr(A, "_bwrap_functional", lambda: True)
    real_which = A.shutil.which
    monkeypatch.setattr(
        A.shutil, "which",
        lambda n, *a, **k: "/usr/bin/bwrap" if n == "bwrap"
        else real_which(n, *a, **k))
    ws = pathlib.Path(tempfile.mkdtemp(prefix="ws_"))
    assert A._bash_isolation_argv("ls", ws, _perms(locked=True, ws=ws))[0] == "bwrap"
