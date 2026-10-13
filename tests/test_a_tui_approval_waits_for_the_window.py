"""An approval asked over the TUI protocol waits for a window's answer.

The worker asks through the same broker the terminal REPL uses, so the
questions, their options and what each answer does must be the REPL's.
Nothing is answered on the user's behalf: a question with nobody to answer
it waits, and a key the question did not offer changes nothing.
"""

from __future__ import annotations

import time

import pytest

from delfin.agent import terminal_confirm as tc
from delfin.agent.service import protocol
from tui_worker_client import WorkerClient, WorkerError


def _ask(w: WorkerClient, text: str) -> tuple[dict, int]:
    start = len(w.events)
    w.request("session.prompt", {"text": text})
    return w.wait_for("permission.asked", after=start)


def test_an_edit_is_asked_with_the_terminals_options_and_waits(tmp_path):
    with WorkerClient(tmp_path, "--permission-mode", "default") as w:
        asked, at = _ask(w, "please edit the parser")
        req = tc.ConfirmRequest(kind=tc.CONFIRM, tool="edit_file",
                                args={"path": "delfin/xtb.py"},
                                preview=asked["preview"])
        assert [o["key"] for o in asked["options"]] == \
            [o.key for o in tc.options_for(req)]
        assert asked["flags"]["write"] is True
        w.wait_for(lambda e: e["type"] == "session.state"
                   and e["state"] == "needs_input", after=at)
        time.sleep(1.0)
        assert not any(e["type"] in ("permission.resolved", "turn.ended")
                       for e in w.events[at:])
        w.request("permission.reply", {"requestId": asked["requestId"], "choice": "y"})
        resolved, _ = w.wait_for("permission.resolved", after=at)
        assert (resolved["decision"], resolved["by"]) == ("approved", "client")
        _, end = w.wait_for("turn.ended", after=at)
        assert any(e["type"] == "tool.ended" and e["name"] == "edit_file"
                   for e in w.events[at:end])


def test_a_key_that_was_not_offered_changes_nothing(tmp_path):
    with WorkerClient(tmp_path, "--permission-mode", "default") as w:
        asked, at = _ask(w, "please edit the parser")
        with pytest.raises(WorkerError) as bad:
            w.request("permission.reply", {"requestId": asked["requestId"],
                                           "choice": "A"})   # shell-only key
        assert bad.value.code == protocol.INVALID_PARAMS
        with pytest.raises(WorkerError) as unknown:
            w.request("permission.reply", {"requestId": "p999", "choice": "y"})
        assert unknown.value.code == protocol.UNKNOWN_REQUEST
        assert not any(e["type"] == "permission.resolved" for e in w.events[at:])
        w.request("permission.reply", {"requestId": asked["requestId"], "choice": "y"})
        w.wait_for("turn.ended", after=at)


def test_a_refusal_with_a_reason_refuses_the_tool(tmp_path):
    with WorkerClient(tmp_path, "--permission-mode", "default") as w:
        asked, at = _ask(w, "edit it")
        w.request("permission.reply", {"requestId": asked["requestId"],
                                       "choice": "n", "reason": "not today"})
        resolved, _ = w.wait_for("permission.resolved", after=at)
        assert resolved["decision"] == "denied"
        denied, _ = w.wait_for("tool.denied", after=at)
        assert denied["name"] == "edit_file"
        w.wait_for("turn.ended", after=at)


def test_abort_ends_the_turn_and_the_next_turn_asks_again(tmp_path):
    with WorkerClient(tmp_path, "--permission-mode", "default") as w:
        asked, at = _ask(w, "edit it")
        w.request("permission.reply", {"requestId": asked["requestId"], "choice": "a"})
        w.wait_for("turn.ended", after=at)
        asked2, at2 = _ask(w, "edit it again")
        assert asked2["requestId"] != asked["requestId"]
        w.request("permission.reply", {"requestId": asked2["requestId"], "choice": "y"})
        w.wait_for("turn.ended", after=at2)


def test_a_shell_command_offers_the_patterns_the_terminal_offers(tmp_path):
    with WorkerClient(tmp_path, "--permission-mode", "default") as w:
        asked, at = _ask(w, "run the bash tests")
        keys = [o["key"] for o in asked["options"]]
        assert keys[:3] == ["y", "n", "A"] and "a" in keys
        assert asked["flags"]["shell"] is True
        w.request("permission.reply", {"requestId": asked["requestId"], "choice": "y"})
        w.wait_for("turn.ended", after=at)


def test_a_plan_and_a_question_are_answered_like_in_the_terminal(tmp_path):
    with WorkerClient(tmp_path) as w:
        asked, at = _ask(w, "make a plan")
        assert asked["kind"] == "plan"
        assert [o["key"] for o in asked["options"]] == ["d", "e", "n"]
        w.request("permission.reply", {"requestId": asked["requestId"], "choice": "e"})
        _, end = w.wait_for("turn.ended", after=at)
        text = "".join(e["text"] for e in w.events[at:end] if e["type"] == "text.ended")
        assert "acceptEdits" in text

        asked, at = _ask(w, "ask me something")
        assert asked["kind"] == "ask" and asked["choices"] == ["CPCM", "SMD", "none"]
        with pytest.raises(WorkerError):
            w.request("permission.reply", {"requestId": asked["requestId"],
                                           "answers": ["COSMO"]})
        w.request("permission.reply", {"requestId": asked["requestId"],
                                       "answers": ["SMD"]})
        _, end = w.wait_for("turn.ended", after=at)
        text = "".join(e["text"] for e in w.events[at:end] if e["type"] == "text.ended")
        assert "Using SMD." in text


def test_a_window_attaching_later_is_told_what_is_waiting(tmp_path):
    with WorkerClient(tmp_path, "--permission-mode", "default") as w:
        asked, at = _ask(w, "edit it")
        again = w.request("initialize", {"protocolVersion": 1})
        assert again["state"] == "needs_input"
        assert [p["requestId"] for p in again["pending"]] == [asked["requestId"]]
        w.request("permission.reply", {"requestId": asked["requestId"], "choice": "y"})
        w.wait_for("turn.ended", after=at)


def test_closing_the_session_refuses_what_is_still_waiting(tmp_path):
    w = WorkerClient(tmp_path, "--permission-mode", "default")
    asked, at = _ask(w, "edit it")
    assert w.close() == 0
    resolved = [e for e in w.events[at:] if e["type"] == "permission.resolved"]
    assert resolved and resolved[0]["by"] == "abort"
    assert w.events[-1]["type"] == "session.ended"
