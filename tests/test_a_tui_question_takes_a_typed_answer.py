"""A question asked over the TUI protocol takes a line typed in the window.

``permission.reply`` with ``text`` answers an ask_user_question the way the
REPL's dialog and the dashboard's chat input answer it
(``approval_answers.typed_answer``): a line that names an option (its label
or its number) is that option, any other line is the answer as typed. An
approval (confirm) and a plan are never answered with typed text: a
sentence is not a yes, and the request stays open.
"""

from __future__ import annotations

import time

import pytest

from delfin.agent.service import protocol
from tui_worker_client import WorkerClient, WorkerError


def _ask(w: WorkerClient, text: str) -> tuple[dict, int]:
    start = len(w.events)
    w.request("session.prompt", {"text": text})
    return w.wait_for("permission.asked", after=start)


def test_a_typed_line_is_the_answer_as_typed(tmp_path):
    with WorkerClient(tmp_path) as w:
        asked, at = _ask(w, "ask me something")
        assert asked["kind"] == "ask"
        got = w.request("permission.reply", {"requestId": asked["requestId"],
                                             "text": "toluene, implicit"})
        assert got == {"applied": True, "message": "answered"}
        resolved, _ = w.wait_for("permission.resolved", after=at)
        assert (resolved["decision"], resolved["by"]) == ("answered", "client")
        _, end = w.wait_for("turn.ended", after=at)
        assert "Using toluene, implicit." in w.text_of(at, end)


def test_a_typed_line_that_names_an_option_is_that_option(tmp_path):
    with WorkerClient(tmp_path) as w:
        asked, at = _ask(w, "ask me something")
        w.request("permission.reply", {"requestId": asked["requestId"], "text": "smd"})
        _, end = w.wait_for("turn.ended", after=at)
        assert "Using SMD." in w.text_of(at, end)
        asked, at = _ask(w, "ask me again")
        w.request("permission.reply", {"requestId": asked["requestId"], "text": "3"})
        _, end = w.wait_for("turn.ended", after=at)
        assert "Using none." in w.text_of(at, end)


def test_typed_text_never_answers_an_approval_or_a_plan(tmp_path):
    with WorkerClient(tmp_path, "--permission-mode", "default") as w:
        for prompt in ("please edit the parser", "run the bash tests"):
            asked, at = _ask(w, prompt)
            assert asked["kind"] == "confirm"
            for wrong in ({"text": "y"}, {"text": "yes please"},
                          {"choice": "y", "text": "yes"}):
                with pytest.raises(WorkerError) as bad:
                    w.request("permission.reply",
                              {"requestId": asked["requestId"], **wrong})
                assert bad.value.code == protocol.INVALID_PARAMS
            time.sleep(0.3)
            assert not any(e["type"] == "permission.resolved" for e in w.events[at:])
            w.request("permission.reply", {"requestId": asked["requestId"], "choice": "n"})
            w.wait_for("turn.ended", after=at)
    with WorkerClient(tmp_path) as w:
        asked, at = _ask(w, "make a plan")
        assert asked["kind"] == "plan"
        with pytest.raises(WorkerError) as bad:
            w.request("permission.reply", {"requestId": asked["requestId"], "text": "d"})
        assert bad.value.code == protocol.INVALID_PARAMS
        assert not any(e["type"] == "permission.resolved" for e in w.events[at:])
        w.request("permission.reply", {"requestId": asked["requestId"], "choice": "n"})
        w.wait_for("turn.ended", after=at)


def test_a_typed_line_and_picked_answers_together_are_refused(tmp_path):
    with WorkerClient(tmp_path) as w:
        asked, at = _ask(w, "ask me something")
        with pytest.raises(WorkerError) as bad:
            w.request("permission.reply", {"requestId": asked["requestId"],
                                           "answers": ["SMD"], "text": "CPCM"})
        assert bad.value.code == protocol.INVALID_PARAMS
        w.request("permission.reply", {"requestId": asked["requestId"], "answers": ["SMD"]})
        _, end = w.wait_for("turn.ended", after=at)
        assert "Using SMD." in w.text_of(at, end)
