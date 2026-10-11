"""U7 phase 3 - RED control for the mark_read delivery wiring.

Package U7 (wave 14) finding 2: ``session_messages.mark_read`` has no
production caller, so a message's status never reaches ``"read"`` and the
sender cannot tell it was seen. ``take()`` records only ``delivered``
(receipt ``taken_at``); ``mark_read`` (session_messages.py:357) is the only
thing that flips the receipt to ``read_at`` and no delivery point calls it.

The task's contract: ``mark_read`` is called where a message is handed to
the model - terminal and dashboard both, in the SUCCESS branch only. The
terminal handoffs are the idle prompt (``_operator_messages``, take at
repl.py:3251, returned as the next prompt) and the running-turn steer
(``_steer_operator_mail``, take at repl.py:2931, ``engine.steer`` at :2937,
put-back at :2940-2952); the dashboard handoff is ``tab_agent.py``
``_deliver_session_messages`` (take at :7862, push_steer/send_on_its_own at
:7871-7874).

RED on the unpatched tree (no delivery point calls mark_read); GREEN once
the operator applies the repl.py/tab_agent.py patches. A red control for
protected files, so it lives under tests/test_zz_* for the operator to run
when building - never pytest inside a script.

Two behavioural through the real public handoffs (idle prompt and running
turn), one source assertion for the dashboard (driving its delivery closure
needs the whole dashboard state machine). ``mark_read`` must fire only for
messages actually handed over: a message whose steer failed and was put back
at repl.py:2940-2952 stays QUEUED by the put-back re-send and is never
marked read.
"""

from __future__ import annotations

import pytest

from delfin.agent import session_messages as msgs


class _Theme:
    def dim(self, t): return t


class _Transcript:
    theme = _Theme()
    def __init__(self):
        self.lines = []
    def chrome(self, line):
        self.lines.append(line)
    #: _steer_operator_mail shows the delivery line on the transcript.
    blank = ""


class _Opts:
    def __init__(self, name="target-s1"):
        self.session_name = name


class _Engine:
    """A steerable terminal engine: ``steer_ok`` decides whether the backend
    accepts the steered mail (True = into the turn; False = put back)."""
    def __init__(self, steer_ok=True):
        self.steer_ok = steer_ok
        self.steered = []

    def steer(self, text: str):
        self.steered.append(text)
        return bool(self.steer_ok)


def _agent(name="target-s1", steer_ok=True):
    from delfin.agent.repl import TerminalAgent
    a = object.__new__(TerminalAgent)
    a.opts = _Opts(name)
    a.engine = _Engine(steer_ok)
    a.transcript = _Transcript()
    a._mail_steer_last = 0.0
    return a


@pytest.fixture
def inbox(tmp_path, monkeypatch):
    monkeypatch.setattr(msgs, "_DIR", tmp_path / "session_inbox")
    return tmp_path


def test_a_message_delivered_to_the_idle_prompt_reaches_read(inbox):
    """The sender must be able to tell the message was seen: after the idle
    prompt hands a message to the model, its status is ``read``, not the
    delivered-with-no-eyes-on-it ``delivered``."""
    sent = msgs.send("target-s1", "look at the failing test", from_key="ops")
    a = _agent()
    got = a._operator_messages("")
    assert "look at the failing test" in got, "the message must reach the prompt"
    assert msgs.status(sent["id"], "ops") == "read", (
        "handed to the model but never marked read")


def test_a_message_steered_into_a_running_turn_is_read(inbox):
    """The running-turn handoff: after engine.steer() succeeds, the message
    is marked read."""
    sent = msgs.send("target-s1", "steer into this turn", from_key="ops")
    a = _agent(steer_ok=True)
    a._steer_operator_mail()
    assert a.engine.steered, "the mail must be steered into the turn"
    assert msgs.status(sent["id"], "ops") == "read", (
        "steered into the turn but never marked read")


def test_a_message_whose_steer_failed_is_put_back_and_not_read(inbox):
    """A backend that takes no steer returns False; the message goes back to
    the inbox, unread, so the prompt path can deliver it. It must never be
    marked read - that would lie to the sender that it was seen."""
    sent = msgs.send("target-s1", "mine for the prompt path", from_key="ops")
    a = _agent(steer_ok=False)
    a._steer_operator_mail()
    # The put-back ends with the message back in the sender's view as queued
    # (it was re-sent to the same inbox), and never as read.
    assert msgs.status(sent["id"], "ops") != "read", (
        "a message that could not be handed over must not read as seen")


def test_the_dashboard_delivery_calls_mark_read_on_success():
    """Source assertion for the dashboard handoff (driving it needs the whole
    dashboard state machine): _deliver_session_messages marks each handed
    message read in its success branch."""
    import pathlib
    from delfin.dashboard import tab_agent as TA
    src = pathlib.Path(TA.__file__).read_text(encoding="utf-8")
    i = src.index("def _deliver_session_messages")
    window = src[i:i + 2200]
    assert "mark_read" in window, (
        "the dashboard hands messages to the model and never marks them read")
    assert "mark_read(" in window


def test_the_dashboard_marks_read_only_what_it_delivered():
    """A held dashboard session (paused, stopped, waiting for the user's
    answer) leaves the text unsent: that message was not read. The read
    receipt therefore hangs on what ``_send_on_its_own`` returned."""
    from pathlib import Path
    import delfin.dashboard.tab_agent as tab
    src = Path(tab.__file__).read_text(encoding="utf-8")
    at = src.index("delivered = _send_on_its_own(text, leave_text_when_held=True)")
    tail = src[at:at + 900]
    assert "if delivered:" in tail
    assert tail.index("if delivered:") < tail.index("_msgs.mark_read(")
