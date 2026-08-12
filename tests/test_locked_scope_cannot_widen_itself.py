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
    """The feature is good; it is the folder it is read from that is wrong.

    Unlocked is now necessary and not sufficient: the workspace also has
    to be one the USER trusted, because "not a registered office folder"
    was true of every directory including a fresh clone. See
    tests/test_a_checked_out_repository_cannot_run_commands.py.
    """
    import delfin.agent.hooks as H
    from delfin.agent import workspace_trust as WT

    ws = _workspace_with_hook()
    WT.trust_workspace(ws, [WT.KIND_HOOKS], actor=WT.ACTOR_USER)
    cfg = H.load_hooks(A._hook_workspace(_perms(locked=False, ws=ws)))
    assert cfg.by_event.get("PreToolUse"), "workspace hooks stopped working"


def test_an_unlocked_but_untrusted_workspace_supplies_none():
    """The lock was never the only reason to refuse."""
    import delfin.agent.hooks as H

    ws = _workspace_with_hook()
    cfg = H.load_hooks(A._hook_workspace(_perms(locked=False, ws=ws)))
    assert not cfg.by_event


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


# ---------------------------------------------------------------------------
# One child's write must not satisfy the guard for the next
# ---------------------------------------------------------------------------

def test_each_subagent_gets_its_own_stale_write_tracker():
    """dataclasses.replace copies field VALUES, so a dict is shared by
    reference. read_tracker is the stale-write guard's memory of "this file
    was read at this mtime, so overwriting it is safe" -- and that question
    is per agent. Shared, sibling A's write bumped the entry and sibling B's
    overwrite then saw an unchanged mtime, passed, and clobbered it."""
    from delfin.agent.subagents import _derive_perms

    parent = A.KitToolPermissions(
        workspace=pathlib.Path(tempfile.mkdtemp(prefix="ws_")))
    parent.read_tracker["shared.txt"] = 111.0

    first = _derive_perms(parent, "acceptEdits")
    second = _derive_perms(parent, "acceptEdits")

    assert first.read_tracker is not parent.read_tracker
    assert first.read_tracker is not second.read_tracker


def test_a_child_still_inherits_what_the_parent_read():
    """Copying, not clearing: a child that has to re-read everything the
    parent read would turn one bug into a pile of wasted rounds."""
    from delfin.agent.subagents import _derive_perms

    parent = A.KitToolPermissions(
        workspace=pathlib.Path(tempfile.mkdtemp(prefix="ws_")))
    parent.read_tracker["shared.txt"] = 111.0
    child = _derive_perms(parent, "acceptEdits")
    assert child.read_tracker.get("shared.txt") == 111.0


def test_a_childs_write_does_not_reach_its_sibling_or_its_parent():
    from delfin.agent.subagents import _derive_perms

    parent = A.KitToolPermissions(
        workspace=pathlib.Path(tempfile.mkdtemp(prefix="ws_")))
    parent.read_tracker["shared.txt"] = 111.0
    first = _derive_perms(parent, "acceptEdits")
    second = _derive_perms(parent, "acceptEdits")

    first.read_tracker["shared.txt"] = 999.0
    assert second.read_tracker["shared.txt"] == 111.0
    assert parent.read_tracker["shared.txt"] == 111.0


def test_the_child_still_inherits_the_lock():
    """The isolation must not have cost the containment."""
    from delfin.agent.subagents import _derive_perms

    parent = _perms(locked=True)
    child = _derive_perms(parent, "bypassPermissions")
    assert child.scope_locked is True
    assert child.workspace == parent.workspace


# ---------------------------------------------------------------------------
# An MCP shell is still a shell
# ---------------------------------------------------------------------------

def _mcp_scene():
    ws = pathlib.Path(tempfile.mkdtemp(prefix="office_ws_"))
    outside = pathlib.Path(tempfile.mkdtemp(prefix="outside_"))
    (outside / "secret.txt").write_text("s", encoding="utf-8")
    perms = _perms(locked=True, ws=ws)
    return _executor(), perms, outside


@pytest.mark.parametrize("template,label", [
    ("cat {out}/secret.txt", "read outside"),
    ("echo leak > {out}/x.txt", "write outside"),
    ("cd .. && cat secret.txt", "climb out"),
])
def test_an_mcp_shell_is_bounded_by_the_folder(template, label):
    """_run_permission_gate carries the deny-list and the secret scan; the
    WORKSPACE boundary lives in two separate gates that only _execute_bash
    called. So an MCP shell arrived with everything applied except the one
    thing a locked scope is for -- routed around by naming another server."""
    executor, perms, outside = _mcp_scene()
    cmd = template.format(out=outside)
    assert executor._gate_mcp_tool(
        "mcp__kit-coding__bash", {"command": cmd}, perms) is not None, label


def test_an_ordinary_mcp_shell_command_still_runs():
    executor, perms, _ = _mcp_scene()
    assert executor._gate_mcp_tool(
        "mcp__kit-coding__bash", {"command": "ls -la"}, perms) is None


def test_unmapped_mcp_tools_are_refused_to_a_scoped_role():
    """The office allow-list judges an MCP tool by its base name, so a
    filesystem server's read tools -- which the MCP gate has no family for
    -- never reach dispatch."""
    from delfin.agent.api_client import _tool_denied_for_role

    for tool in ("mcp__fs__read_text_file", "mcp__fs__list_directory",
                 "mcp__other__shell", "mcp__delfin-tools__run_application"):
        assert _tool_denied_for_role("office_agent", tool), tool


# ---------------------------------------------------------------------------
# Network egress: gated where the user has not opted out, and only there
# ---------------------------------------------------------------------------

def _net_perms(mode):
    perms = _perms(locked=True)
    perms.mode = mode
    return perms


@pytest.mark.parametrize("mode", ["default", "acceptEdits"])
def test_a_locked_session_asks_before_data_leaves(mode):
    """The folder lock bounds where data may be WRITTEN. A record leaves
    through a fetched URL without any path crossing the boundary, and the
    egress scanner reads shell commands -- this is not a shell command."""
    asked: list[str] = []
    perms = _net_perms(mode)
    perms.confirm_callback = lambda n, a, prev: (asked.append(n), True)[1]
    _executor()._run_permission_gate("web_fetch", {"url": "https://x/y"}, perms)
    assert asked == ["web_fetch"]


def test_bypass_asks_nothing_because_that_is_what_it_says():
    """A gate that overrides the profile turns the setting into a
    suggestion. The first field report after this shipped was a user asking
    why Bypass still prompted."""
    asked: list[str] = []
    perms = _net_perms("bypassPermissions")
    perms.confirm_callback = lambda n, a, prev: (asked.append(n), True)[1]
    assert _executor()._run_permission_gate(
        "web_fetch", {"url": "https://x/y"}, perms) is None
    assert asked == []


@pytest.mark.parametrize("mode", ["default", "acceptEdits", "bypassPermissions"])
def test_the_mcp_variant_behaves_exactly_like_the_native_one(mode):
    """Reaching the internet by naming a different server is the same
    bypass-by-another-name the shell branch had."""
    native: list[str] = []
    through_mcp: list[str] = []

    p1 = _net_perms(mode)
    p1.confirm_callback = lambda n, a, prev: (native.append(n), True)[1]
    _executor()._run_permission_gate("web_search", {"query": "x"}, p1)

    p2 = _net_perms(mode)
    p2.confirm_callback = lambda n, a, prev: (through_mcp.append(n), True)[1]
    _executor()._gate_mcp_tool(
        "mcp__kit-coding__web_search", {"query": "x"}, p2)

    assert bool(native) == bool(through_mcp), mode


def test_an_unlocked_session_never_asks_about_the_network():
    asked: list[str] = []
    perms = _perms(locked=False)
    perms.confirm_callback = lambda n, a, prev: (asked.append(n), True)[1]
    assert _executor()._run_permission_gate(
        "web_fetch", {"url": "https://x"}, perms) is None
    assert asked == []
