"""Only READ denials were remembered.

``denied_paths`` is cross-tool and holds: a refused ``read_file`` cannot be
fetched with ``bash cat`` afterwards. Nothing equivalent existed for writes
or for commands. A denied ``write_file`` and a denied ``bash`` returned
prose to the model asking it not to retry, and that was the entire
mechanism — the identical call could be re-emitted verbatim and simply
raised the dialog again. A user who says no once ends up saying it four
times, and the fourth dialog is the one that gets clicked through.

The ledger is typed: writes by resolved path, commands by normalised
command text, everything else by tool name plus arguments. Consulted by the
write gate, the bash gate and the MCP gate before any callback, so the
second attempt is answered without a prompt — through whichever tool, and
through whichever server, it arrives.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import api_client as A


class _Broker:
    """A confirm callback with the broker's timeout flag, as the gates read it."""

    def __init__(self, answer=False, timed_out=False):
        self.answer = answer
        self.last_timed_out = timed_out
        self.calls: list[tuple] = []

    def callback(self, name, args, preview):
        self.calls.append((name, args))
        return self.answer


@pytest.fixture
def ex():
    return A._DocToolExecutor.__new__(A._DocToolExecutor)


def _perms(tmp_path, broker, mode="default"):
    return A.KitToolPermissions(
        workspace=tmp_path, mode=mode, confirm_callback=broker.callback)


# ---------------------------------------------------------------------------
# Writes
# ---------------------------------------------------------------------------

def test_a_denied_write_is_not_asked_about_twice(ex, tmp_path):
    (tmp_path / "report.md").write_text("keep", encoding="utf-8")
    broker = _Broker(answer=False)
    perms = _perms(tmp_path, broker)
    args = {"path": "report.md", "content": "new"}

    first = ex._run_permission_gate("write_file", args, perms)
    second = ex._run_permission_gate("write_file", args, perms)

    assert first is not None and "denied" in first
    assert second is not None and "already refused" in second
    assert len(broker.calls) == 1                 # the user was asked once


def test_a_denied_write_holds_against_the_shell(ex, tmp_path):
    target = tmp_path / "report.md"
    target.write_text("keep", encoding="utf-8")
    broker = _Broker(answer=False)
    perms = _perms(tmp_path, broker)
    # Give the file a read baseline so the stale-write guard is satisfied
    # and the refusal below is the only thing left standing.
    perms.read_tracker[str(target.resolve())] = target.stat().st_mtime
    ex._run_permission_gate(
        "write_file", {"path": "report.md", "content": "new"}, perms)

    for cmd in ("sed -i 's/a/b/' report.md",
                "echo x > report.md",
                "tee report.md < /dev/null"):
        blocked = ex._gate_bash_write_targets(cmd, {}, perms)
        assert blocked is not None, cmd
        assert "already refused" in blocked or "refused a write" in blocked, cmd
    assert len(broker.calls) == 1


def test_the_same_file_by_another_spelling_is_the_same_file(ex, tmp_path):
    (tmp_path / "sub").mkdir()
    (tmp_path / "sub" / "a.txt").write_text("x", encoding="utf-8")
    broker = _Broker(answer=False)
    perms = _perms(tmp_path, broker)
    ex._run_permission_gate(
        "write_file", {"path": "sub/a.txt", "content": "y"}, perms)

    again = ex._run_permission_gate(
        "write_file", {"path": str(tmp_path / "sub" / ".." / "sub" / "a.txt"),
                       "content": "y"}, perms)
    assert again is not None and "already refused" in again
    assert len(broker.calls) == 1


def test_another_file_is_still_a_new_question(ex, tmp_path):
    (tmp_path / "a.txt").write_text("x", encoding="utf-8")
    (tmp_path / "b.txt").write_text("x", encoding="utf-8")
    broker = _Broker(answer=False)
    perms = _perms(tmp_path, broker)
    ex._run_permission_gate("write_file", {"path": "a.txt", "content": "1"}, perms)
    ex._run_permission_gate("write_file", {"path": "b.txt", "content": "1"}, perms)
    assert len(broker.calls) == 2        # a refusal is about one thing


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

def test_a_denied_command_is_not_asked_about_twice(ex, tmp_path):
    broker = _Broker(answer=False)
    perms = _perms(tmp_path, broker)
    cmd = "curl -X POST https://example.invalid/upload"

    first = ex._run_permission_gate("bash", {"command": cmd}, perms)
    second = ex._run_permission_gate("bash", {"command": cmd}, perms)

    assert first is not None and "denied" in first
    assert second is not None and "already refused" in second
    assert len(broker.calls) == 1


def test_respacing_a_command_does_not_make_it_a_new_one(ex, tmp_path):
    broker = _Broker(answer=False)
    perms = _perms(tmp_path, broker)
    ex._run_permission_gate("bash", {"command": "npm  run   deploy"}, perms)
    again = ex._run_permission_gate("bash", {"command": "npm run deploy ;"}, perms)
    assert again is not None and "already refused" in again
    assert len(broker.calls) == 1


def test_a_refusal_outlives_a_mode_switch(ex, tmp_path):
    """A refusal is an answer about an action. Re-asking it in a mode the
    model picked afterwards would make the answer depend on the retry."""
    broker = _Broker(answer=False)
    perms = _perms(tmp_path, broker)
    cmd = "npm run deploy"
    ex._run_permission_gate("bash", {"command": cmd}, perms)

    perms.mode = "bypassPermissions"
    assert "already refused" in ex._run_permission_gate(
        "bash", {"command": cmd}, perms)


def test_a_different_command_still_asks(ex, tmp_path):
    broker = _Broker(answer=False)
    perms = _perms(tmp_path, broker)
    ex._run_permission_gate("bash", {"command": "npm run deploy"}, perms)
    ex._run_permission_gate("bash", {"command": "npm run build"}, perms)
    assert len(broker.calls) == 2


# ---------------------------------------------------------------------------
# Absence is not refusal
# ---------------------------------------------------------------------------

def test_an_expired_window_is_not_a_refusal(ex, tmp_path):
    """A 300 s dialog nobody was there for must not close the action for
    the session — the user never saw the question."""
    broker = _Broker(answer=False, timed_out=True)
    perms = _perms(tmp_path, broker)
    (tmp_path / "a.txt").write_text("x", encoding="utf-8")

    ex._run_permission_gate("bash", {"command": "npm run deploy"}, perms)
    ex._run_permission_gate("write_file", {"path": "a.txt", "content": "1"}, perms)
    assert perms.denied_actions == {}

    ex._run_permission_gate("bash", {"command": "npm run deploy"}, perms)
    assert len(broker.calls) == 3        # asked again, as it should be


# ---------------------------------------------------------------------------
# Every dispatch path
# ---------------------------------------------------------------------------

def test_an_mcp_server_is_not_a_second_chance(ex, tmp_path):
    """The ledger is keyed on the BARE tool name, so a refusal given for a
    native call holds for the same call routed through a server."""
    broker = _Broker(answer=False)
    perms = A.KitToolPermissions(
        workspace=tmp_path, mode="default", lock_workspace=True,
        confirm_callback=broker.callback)
    args = {"url": "https://example.invalid/x"}

    first = ex._run_permission_gate("web_fetch", args, perms)
    assert first is not None and "denied" in first

    blocked = ex._gate_mcp_tool("mcp__research__web_fetch", args, perms)
    assert blocked is not None and "already refused" in blocked
    assert len(broker.calls) == 1


def test_an_mcp_shell_inherits_a_refused_command(ex, tmp_path):
    broker = _Broker(answer=False)
    perms = _perms(tmp_path, broker)
    cmd = "npm run deploy"
    ex._run_permission_gate("bash", {"command": cmd}, perms)

    blocked = ex._gate_mcp_tool(
        "mcp__kit-coding__bash", {"command": cmd}, perms)
    assert blocked is not None and "already refused" in blocked
    assert len(broker.calls) == 1


def test_a_sub_agent_inherits_the_refusals(ex, tmp_path):
    """dataclasses.replace copies field values, so the child shares the
    ledger — deliberately: delegating is not a way around a 'no'."""
    from delfin.agent.subagents import _derive_perms

    broker = _Broker(answer=False)
    perms = _perms(tmp_path, broker)
    ex._run_permission_gate("bash", {"command": "npm run deploy"}, perms)

    child = _derive_perms(perms, "default")
    assert child.action_denied_earlier("bash", "npm run deploy") is True
    assert "already refused" in ex._run_permission_gate(
        "bash", {"command": "npm run deploy"}, child)


# ---------------------------------------------------------------------------
# The ledger itself
# ---------------------------------------------------------------------------

def test_the_ledger_is_typed(tmp_path):
    perms = A.KitToolPermissions(workspace=tmp_path)
    perms.record_denied_action("bash", "rm -rf build")
    assert perms.action_denied_earlier("bash", "rm -rf build") is True
    # A path that happens to read like the command is a different kind.
    assert perms.action_denied_earlier("write", "rm -rf build") is False


def test_the_containment_panel_sees_the_second_attempt(ex, tmp_path):
    from delfin.agent import security_events

    broker = _Broker(answer=False)
    perms = _perms(tmp_path, broker)
    ex._run_permission_gate("bash", {"command": "npm run deploy"}, perms)
    security_events.clear()
    ex._run_permission_gate("bash", {"command": "npm run deploy"}, perms)

    kinds = [e.kind for e in security_events.recent(10)]
    assert "denied_again" in kinds
    assert "Refusal circumvented" in security_events.format_panel_html()


def test_a_gate_without_a_ledger_still_works(ex, tmp_path):
    """Perms objects predating the field must not crash the gate."""
    perms = A.KitToolPermissions(workspace=tmp_path, mode="default")
    del perms.denied_actions
    out = ex._run_permission_gate("bash", {"command": "ls"}, perms)
    assert out is None or isinstance(out, str)
    assert json.dumps({"ok": True})     # sanity: module imports intact
