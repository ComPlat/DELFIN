"""One rung per thing that can be turned off.

    Plan          read only
    Ask All       asks before a file changes AND before a shell command
    Accept Edits  file changes go through, shell still asks
    Bypass        asks nothing

Until the write confirmation existed, nothing ever asked before a write
inside the workspace. Ask All and Accept Edits therefore decided every
case identically -- measured through the real gate, for writes, edits,
calc-confirm paths and bash alike -- which left the middle rung inert and
the top one describing something it did not do.

The lock and the ladder are independent: the lock decides WHERE the agent
may write, the ladder how much it asks first. A locked office session
still gets all four rungs.
"""

from __future__ import annotations

import pathlib
import tempfile

import pytest

from delfin.dashboard.tab_agent import _perm_options_for_mode

_SOURCE = pathlib.Path(
    __import__("delfin.dashboard.tab_agent", fromlist=["x"]).__file__
).read_text(encoding="utf-8")


def _gate():
    from delfin.agent import api_client as A
    ex = A._DocToolExecutor.__new__(A._DocToolExecutor)
    return A, ex._run_permission_gate


def _scene(mode, answer=True, locked=True):
    A, gate = _gate()
    ws = pathlib.Path(tempfile.mkdtemp(prefix="ws_"))
    (ws / "a.txt").write_text("x")
    seen: list[str] = []
    p = A.KitToolPermissions(
        workspace=ws, agent_role="office_agent" if locked else "")
    p.mode = mode
    p.confirm_callback = lambda n, a, prev: (seen.append(n), answer)[1]
    return gate, ws, p, seen


# ---------------------------------------------------------------------------
# The ladder
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("mode,asks_write,asks_bash", [
    ("default", True, True),
    ("acceptEdits", False, True),
    ("bypassPermissions", False, False),
])
def test_each_rung_turns_off_exactly_one_more_thing(mode, asks_write, asks_bash):
    gate, ws, perms, seen = _scene(mode)

    seen.clear()
    assert gate("write_file", {"path": str(ws / "n.txt"), "content": "y"},
                perms) is None
    assert bool(seen) is asks_write, f"{mode}: write prompt should be {asks_write}"

    seen.clear()
    assert gate("bash", {"command": "curl http://example.invalid"},
                perms) is None
    assert bool(seen) is asks_bash, f"{mode}: bash prompt should be {asks_bash}"


def test_plan_changes_nothing_at_all():
    gate, ws, perms, _seen = _scene("plan")
    assert gate("write_file", {"path": str(ws / "n.txt"), "content": "y"},
                perms) is not None
    assert gate("bash", {"command": "ls"}, perms) is not None


@pytest.mark.parametrize("tool,args_for", [
    ("write_file", lambda ws: {"path": str(ws / "n.txt"), "content": "y"}),
    ("edit_file", lambda ws: {"path": str(ws / "a.txt"),
                              "old_string": "x", "new_string": "y"}),
    ("multi_edit", lambda ws: {"path": str(ws / "a.txt"),
                               "edits": [{"old_string": "x",
                                          "new_string": "y"}]}),
])
def test_every_write_tool_goes_through_the_same_prompt(tool, args_for):
    """One gate for all of them -- a tool that skipped it would be the way
    around the setting."""
    gate, ws, perms, seen = _scene("default")
    seen.clear()
    assert gate(tool, args_for(ws), perms) is None
    assert seen, f"{tool} wrote without asking under Ask All"


def test_a_refusal_is_reported_so_the_model_stops():
    gate, ws, perms, _ = _scene("default", answer=False)
    err = gate("write_file", {"path": str(ws / "n.txt"), "content": "y"}, perms)
    assert err is not None
    assert "denied" in err
    assert "Do NOT retry" in err, (
        "a denial the model reads as a transient failure gets worked around")


# ---------------------------------------------------------------------------
# Head-less
# ---------------------------------------------------------------------------

def test_a_head_less_run_can_still_write():
    """Bash may refuse without a callback because its auto-allow list is the
    escape hatch. Writing has none, so refusing would leave every unattended
    run unable to produce anything."""
    A, gate = _gate()
    ws = pathlib.Path(tempfile.mkdtemp(prefix="ws_"))
    p = A.KitToolPermissions(workspace=ws, agent_role="office_agent")
    p.mode = "default"
    p.confirm_callback = None
    assert gate("write_file", {"path": str(ws / "n.txt"), "content": "y"},
                p) is None


# ---------------------------------------------------------------------------
# The lock is a separate axis
# ---------------------------------------------------------------------------

def test_the_prompt_does_not_replace_the_lock():
    """Approving cannot put a file outside the folder: outside is refused
    before any dialog, and no answer can grant it."""
    gate, ws, perms, seen = _scene("default", answer=True)
    seen.clear()
    err = gate("write_file", {"path": "/etc/shadow", "content": "y"}, perms)
    assert err is not None
    assert not seen, "a path outside a locked scope was offered for approval"


def test_a_locked_mode_still_offers_the_whole_ladder():
    assert [v for _, v in _perm_options_for_mode("office")] == [
        "plan", "ask_all", "repo_free", "all_free"]


def test_every_mode_gets_the_same_ladder():
    for mode in ("office", "code", "solo", "dashboard", "reviewed", ""):
        assert _perm_options_for_mode(mode) == _perm_options_for_mode("code"), mode


# ---------------------------------------------------------------------------
# Labels
# ---------------------------------------------------------------------------

def test_the_labels_say_what_the_rung_does():
    labels = [lbl for lbl, _ in _perm_options_for_mode("office")]
    assert labels == ["Plan", "Ask All", "Accept Edits", "Bypass"]


def test_no_rung_is_called_auto():
    """"Auto" says a mode is convenient, not what it stops asking."""
    for label, _ in _perm_options_for_mode("code"):
        assert "auto" not in label.lower(), label


def test_the_stored_values_are_unchanged():
    """Saved sessions and settings carry these; only labels are cosmetic."""
    assert [v for _, v in _perm_options_for_mode("code")] == [
        "plan", "ask_all", "repo_free", "all_free"]
