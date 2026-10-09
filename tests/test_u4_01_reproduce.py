"""U4 phase 1 — reproduce the three reported problems.

  1. The new-session picker offers a session-worktree directory as a plain
     selectable folder whether or not a live session holds it. Reproduced via
     ``agent_sessions.workspace_choices`` (the public path the picker calls).
     The start path refuses a held directory only later (``_on_start``,
     agent_sessions.py:1445-1452); the choices list already offered it, with
     no busy/holder marker (comment at agent_sessions.py:1433-1439).

  2. Root-folder download (#86): FIXED on main
     (``calc_update_download_btn``, tab_calculations_browser.py:6186, plus
     the every-selection re-evaluation at 3158-3170; commits 1250d446 +
     f49caf46). The logic sits inside the nested ``create_tab`` builder, so
     the tab-level behaviour cannot be driven through a module-level public
     path in a unit test; it is recorded here as NOT independently driven.

  3. The agent "step budget" lives in two PROTECTED resolvers, not the
     dashboard: ``_resolve_max_tool_rounds`` (api_client.py:9497, per-turn
     round cap; precedence settings -> model profile -> 500) and
     ``_subagent_limits`` (subagents.py:109; ``max_tool_calls`` default 120
     for subagents). Called directly to pin current defaults/precedence.
     Review (QR) confirmed item 3 HOLDS on main: the dashboard writes both
     fields and save/load round-trips, so there is no "not saved" defect
     left -- only the two pins below, recorded so a regression that takes a
     default too low again is caught.

The picker assertions in the first group are the reproduction for phase 2:
they are red now (the held worktree is offered, unbounded) and flip green
after the fix.
"""
from __future__ import annotations

import json
import os
from types import SimpleNamespace

from delfin.dashboard import agent_sessions as asm


def _fake_ctx(project):
    """A minimal dashboard context rooted at ``project`` (a tmp_path)."""
    return SimpleNamespace(repo_dir=project, agent_dir=project,
                           calc_dir=project)


def _make_session_worktree(project, name="sb") -> "object":
    """A directory that looks like a session worktree: it carries the
    ``.delfin/session_worktree.json`` sidecar agent_sessions writes."""
    wt = project / name
    wt.mkdir(parents=True, exist_ok=True)
    side = wt / ".delfin"
    side.mkdir(parents=True, exist_ok=True)
    (side / "session_worktree.json").write_text(json.dumps({
        "path": str(wt), "branch": "session/abc123", "repo_dir": str(project),
        "base_ref": "deadbeef", "created_at": 0.0, "host": "somehost",
        "pid": -1, "key": "k",
    }), encoding="utf-8")
    return wt


# ---------------------------------------------------------------------------
# Item 1 — a session worktree is offered unboundedly by the picker
# ---------------------------------------------------------------------------

def test_a_session_worktree_dir_is_offered_by_the_picker(tmp_path):
    """Green (documentation): the picker lists a session-worktree directory.

    This proves the folder-picker path includes session worktrees at all —
    the precondition of the bug. A here ``sb`` (sidecar present) appears in
    the choices a new session would get.
    """
    project = tmp_path / "proj"
    project.mkdir()
    wt = _make_session_worktree(project)
    choices = asm.workspace_choices(_fake_ctx(project), typed=str(project) + "/")
    assert str(wt) in choices


def test_a_held_worktree_is_still_selectable(tmp_path, monkeypatch):
    """RED: a worktree a live session holds is offered, not refused in-band.

    When the picker knows the chosen directory is held, it must not offer it
    as a start folder. Today ``workspace_choices`` has no holding concept, so
    the held ``sb`` is offered. This assertion drives the phase-2 fix.
    """
    project = tmp_path / "proj"
    project.mkdir()
    wt = _make_session_worktree(project)
    holder = "My Session on kitzbuhel"
    monkeypatch.setattr(asm, "_live_session_in", lambda ws: (
        holder if str(workspace := str(ws)) == str(wt) else ""))
    choices = asm.workspace_choices(_fake_ctx(project), typed=str(project) + "/")
    assert str(wt) not in choices, (
        "a held worktree must not be offered as a start folder")


def test_a_held_worktree_state_names_who_holds_it(tmp_path, monkeypatch):
    """RED: a classifier must report the holder of a held worktree.

    There is no such classifier today; phase 2 adds one so the picker can
    show "busy by <holder>". Asserting its result here drives its contract.
    """
    project = tmp_path / "proj"
    project.mkdir()
    wt = _make_session_worktree(project)
    holder = "My Session on kitzbuhel"
    monkeypatch.setattr(asm, "_live_session_in", lambda ws: (
        holder if str(workspace := str(ws)) == str(wt) else ""))
    state = asm.session_worktree_state(str(wt))
    assert state["is_worktree"] is True
    assert state["holder"] == holder


def test_an_orphaned_worktree_is_releasable(tmp_path, monkeypatch):
    """RED: a worktree whose session is gone is offered for release.

    Phase 2 adds a releasable flag (sidecar present, no live holder) so the
    picker can offer the existing ``release_session_worktree`` for it. Today
    no such flag exists and ``_live_session_in`` returns "" (session gone)
    even for a live session absent from the presence file.
    """
    project = tmp_path / "proj"
    project.mkdir()
    wt = _make_session_worktree(project)
    monkeypatch.setattr(asm, "_live_session_in", lambda ws: "")
    state = asm.session_worktree_state(str(wt))
    assert state["is_worktree"] is True
    assert state["holder"] == ""
    assert state["releasable"] is True


# ---------------------------------------------------------------------------
# Item 3 — step-budget resolvers (protected files), pinned read-side
# ---------------------------------------------------------------------------

def test_main_loop_default_round_cap_is_the_fallback(tmp_path, monkeypatch):
    """Pin the per-turn round-cap default with nothing configured.

    ``_resolve_max_tool_rounds('')`` falls back to the documented 500. This
    is the main-loop half of #72/#87; changing it is a protected-file patch.
    """
    monkeypatch.setattr("delfin.user_settings.load_settings", lambda: {})
    from delfin.agent.api_client import _resolve_max_tool_rounds
    assert _resolve_max_tool_rounds("") == 500


def test_subagent_default_tool_call_cap(tmp_path, monkeypatch):
    """Pin the default subagent tool-call budget (120 calls).

    ``agent.subagents.max_tool_calls`` was already raised once (subagents.py
    comment 93-96); this pins the current floor so a regression that takes it
    "too low" again is caught. This is the "for subagents too" half of #72/#87.
    """
    monkeypatch.setattr("delfin.user_settings.load_settings", lambda: {})
    from delfin.agent.subagents import _subagent_limits
    limits = _subagent_limits()
    assert limits["max_tool_calls"] >= 120
def test_a_worktree_with_alive_owning_pid_is_not_offered(tmp_path, monkeypatch):
    """RED: a worktree whose owning process is still alive is not offered,
    even when its session has left the presence file (stale presence).

    This is the reviewer's adversarial case: judging by presence alone would
    offer it (`_live_session_in` says ""). The exclusion must consult the
    sidecar's owning pid. The sidecar is written with the current process's
    pid, which is definitely alive, and ``_live_session_in`` is forced to ""
    (session gone) -- yet the tree must not appear as a start folder.
    """
    project = tmp_path / "proj"
    project.mkdir()
    wt = project / "sb"
    wt.mkdir(parents=True, exist_ok=True)
    side = wt / ".delfin"
    side.mkdir(parents=True, exist_ok=True)
    (side / "session_worktree.json").write_text(json.dumps({
        "path": str(wt), "branch": "session/abc123", "repo_dir": str(project),
        "base_ref": "deadbeef", "created_at": 0.0, "host": "somehost",
        "pid": os.getpid(), "key": "k",
    }), encoding="utf-8")
    monkeypatch.setattr(asm, "_live_session_in", lambda ws: "")
    choices = asm.workspace_choices(_fake_ctx(project), typed=str(project) + "/")
    assert str(wt) not in choices, (
        "a worktree whose owning process is alive must not be a start folder")


def test_a_plain_directory_is_still_offered(tmp_path, monkeypatch):
    """Green (specificity): excluding a held worktree must not sweep up plain
    directories -- a normal folder in the tree stays a valid start folder."""
    project = tmp_path / "proj"
    project.mkdir()
    plain = project / "plain"
    plain.mkdir()
    monkeypatch.setattr(asm, "_live_session_in", lambda ws: "")
    choices = asm.workspace_choices(_fake_ctx(project), typed=str(project) + "/")
    assert str(plain) in choices
