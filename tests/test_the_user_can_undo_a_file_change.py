"""The one command that restores content has to be reachable by typing it.

`/undo-file` is advertised in the palette as the command that "restores
content, unlike /undo which only drops context". It was not in the set of
tokens the send path recognises as built-in commands, so typing it fell
through to the command-template/skill lookup, matched nothing, and the
literal string went to the model as a chat message. Picking it from the
palette only pastes the same text into the box. Its handler was live code
reachable only by the MODEL, through the action dispatcher.

The previous version of this file grepped the dashboard source for the
string `"/undo-file"` and for `.revert(` — both were there the whole time
the command did nothing. These tests drive the routing predicate the send
path actually applies, and the renderers the handler actually calls, with
output from the real journal.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import change_journal as cj
from delfin.dashboard import tab_agent


@pytest.fixture
def home(monkeypatch, tmp_path):
    monkeypatch.setattr(Path, "home", lambda: tmp_path / "home")
    (tmp_path / "home").mkdir()
    return tmp_path / "home"


@pytest.fixture
def ws(tmp_path):
    d = tmp_path / "ws"
    d.mkdir()
    return d


# ---------------------------------------------------------------------------
# Routing: typing the command reaches the dispatcher
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("typed", [
    "/undo-file",
    "/undo-file list",
    "/undo-file last",
    "/UNDO-FILE session",
])
def test_typing_the_command_reaches_the_dispatcher(typed):
    assert tab_agent.slash_command_routes_to_builtin(typed), (
        f"{typed!r} is sent to the model as chat instead of being executed")


def test_the_context_undo_still_routes_and_is_a_different_command():
    assert tab_agent.slash_command_routes_to_builtin("/undo")
    assert "/undo" in tab_agent._BUILTIN_SLASH_PREFIXES
    assert "/undo-file" in tab_agent._BUILTIN_SLASH_PREFIXES


def test_a_file_path_is_not_a_slash_command():
    """The prefix set exists to keep /home/... out of the dispatcher."""
    assert not tab_agent.slash_command_routes_to_builtin("/home/user/x.txt")
    assert not tab_agent.slash_command_routes_to_builtin("/etc/passwd")


def test_an_unknown_token_still_goes_to_the_template_or_skill_route():
    assert not tab_agent.slash_command_routes_to_builtin("/my-custom-command")


def test_every_command_the_dispatcher_handles_can_be_typed():
    """The general form of the bug: a handler nobody can reach.

    Any ``cmd ==``/``cmd.startswith(`` branch in the slash dispatcher is
    a command the user is meant to be able to run. If its token is not
    in the prefix set, typing it produces a chat message instead.
    """
    import re

    src = (Path(tab_agent.__file__)).read_text(encoding="utf-8")
    body = src[src.index("def _handle_slash_command"):]
    handled = set()
    for m in re.finditer(
            r'cmd(?:\.strip\(\))?\s*(?:==|\.startswith\(|in \()\s*\(?'
            r'["\'](/[^"\']+)["\']', body):
        handled.add(m.group(1).split()[0])
    assert "/undo-file" in handled, "the handler itself disappeared"
    unreachable = sorted(
        t for t in handled
        if not tab_agent.slash_command_routes_to_builtin(t + " x"))
    assert not unreachable, (
        f"handled but not typeable — these reach the model as chat: "
        f"{unreachable}")


def test_every_advertised_palette_command_is_typeable():
    """What the palette offers must be what the box accepts.

    The delegation commands are the deliberate exception: they are
    expanded into a prompt by the template/skill route rather than
    executed by the dispatcher.
    """
    expansion_route = {"/explore", "/review", "/delegate"}
    broken = []
    for _cat, cmd, _desc, _args in tab_agent._SLASH_COMMANDS:
        token = cmd.split()[0]
        if token in expansion_route:
            continue
        if not tab_agent.slash_command_routes_to_builtin(cmd):
            broken.append(cmd)
    assert not broken, f"palette offers commands that are not executed: {broken}"


# ---------------------------------------------------------------------------
# What the handler renders — driven with real journal output
# ---------------------------------------------------------------------------

def test_a_restore_is_reported_with_the_file_it_restored(home, ws):
    p = ws / "a.txt"
    p.write_text("agent\n", encoding="utf-8")
    cj.record_change("s-ok", tool="write_file", path=p,
                     old_text="user\n", new_text="agent\n")

    out = tab_agent.format_undo_result(
        "last", cj.revert("s-ok", scope="last", workspace=ws))
    assert "Restored:" in out
    assert str(p) in out
    assert p.read_text(encoding="utf-8") == "user\n"


def test_a_conflict_says_why_the_file_was_not_touched(home, ws):
    p = ws / "b.txt"
    p.write_text("agent\n", encoding="utf-8")
    cj.record_change("s-conf", tool="write_file", path=p,
                     old_text="user\n", new_text="agent\n")
    p.write_text("the user edited it afterwards\n", encoding="utf-8")

    out = tab_agent.format_undo_result(
        "last", cj.revert("s-conf", scope="last", workspace=ws))
    assert "Conflicts" in out
    assert "content changed since the agent's edit" in out, (
        "the reason was dropped — a bare path cannot be acted on")


def test_a_skip_says_why_it_was_skipped(home, ws):
    big = "x" * (cj.MAX_PRE_IMAGE_BYTES + 10)
    p = ws / "big.txt"
    p.write_text("small\n", encoding="utf-8")
    cj.record_change("s-skip", tool="write_file", path=p,
                     old_text=big, new_text="small\n")

    out = tab_agent.format_undo_result(
        "last", cj.revert("s-skip", scope="last", workspace=ws))
    assert "Skipped:" in out
    assert "truncated" in out
    assert str(p) in out


def test_the_two_refusals_do_not_read_the_same(home, ws):
    """A refusal to corrupt the file and a refusal to overwrite the
    user's later edit are different situations with different answers."""
    a = ws / "trunc.txt"
    a.write_text("small\n", encoding="utf-8")
    cj.record_change("s-two", tool="write_file", path=a,
                     old_text="x" * (cj.MAX_PRE_IMAGE_BYTES + 10),
                     new_text="small\n")
    b = ws / "moved.txt"
    b.write_text("agent\n", encoding="utf-8")
    cj.record_change("s-two", tool="write_file", path=b,
                     old_text="user\n", new_text="agent\n")
    b.write_text("user came back\n", encoding="utf-8")

    out = tab_agent.format_undo_result(
        "session", cj.revert("s-two", scope="session", workspace=ws))
    trunc_line = next(ln for ln in out.splitlines() if "trunc.txt" in ln)
    moved_line = next(ln for ln in out.splitlines() if "moved.txt" in ln)
    assert trunc_line != moved_line
    assert "truncated" in trunc_line
    assert "changed since" in moved_line


def test_nothing_to_undo_is_still_said(home, ws):
    out = tab_agent.format_undo_result(
        "last", cj.revert("s-empty", scope="last", workspace=ws))
    assert "Nothing to undo." in out


# ---------------------------------------------------------------------------
# The listing the user reads before acting
# ---------------------------------------------------------------------------

def test_the_listing_marks_what_cannot_be_restored(home, ws, monkeypatch):
    monkeypatch.setattr(cj, "MAX_ENTRIES_PER_SESSION", 2)
    p = ws / "c.txt"
    for i in range(4):
        p.write_text(f"v{i + 1}", encoding="utf-8")
        cj.record_change("s-list", tool="edit_file", path=p,
                         old_text=f"v{i}", new_text=f"v{i + 1}")

    out = tab_agent.format_undo_listing(cj.list_changes("s-list"))
    assert "dropped at the session cap" in out
    assert out.count(str(p)) == 4       # every change is still listed


def test_the_listing_shows_an_undo_as_an_undo(home, ws):
    p = ws / "d.txt"
    p.write_text("agent\n", encoding="utf-8")
    cj.record_change("s-undo", tool="write_file", path=p,
                     old_text="user\n", new_text="agent\n")
    cj.revert("s-undo", scope="last", workspace=ws)

    out = tab_agent.format_undo_listing(cj.list_changes("s-undo"))
    assert "(undo)" in out
    assert "already undone" in out


# ---------------------------------------------------------------------------
# The journal side still behaves
# ---------------------------------------------------------------------------

def test_the_journal_functions_are_what_the_handler_expects(home, ws):
    target = ws / "a.txt"
    target.write_text("changed\n", encoding="utf-8")
    cj.record_change("s-undo-test", tool="write_file", path=str(target),
                     old_text="original\n", new_text="changed\n")

    listed = cj.list_changes("s-undo-test")
    assert listed and listed[-1].get("path", "").endswith("a.txt")

    out = cj.revert("s-undo-test", scope="last", workspace=ws)
    assert set(out) >= {"reverted", "conflicts", "skipped"}
    assert target.read_text(encoding="utf-8") == "original\n"


def test_a_file_changed_since_the_edit_is_a_conflict_not_an_overwrite(home, ws):
    target = ws / "b.txt"
    target.write_text("agent wrote this\n", encoding="utf-8")
    cj.record_change("s-undo-conflict", tool="write_file", path=str(target),
                     old_text="original\n", new_text="agent wrote this\n")
    target.write_text("the user then edited it\n", encoding="utf-8")

    out = cj.revert("s-undo-conflict", scope="last", workspace=ws)
    assert out["conflicts"], out
    assert target.read_text(encoding="utf-8") == "the user then edited it\n"
