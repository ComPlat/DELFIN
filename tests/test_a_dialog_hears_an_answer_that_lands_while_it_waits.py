"""A dialog hears an outside answer that lands while it is waiting.

The dialog checked for an outside answer (approvals approve/deny) only
before it started to wait for a key. A supervisor answers a dialog that
is already on screen, so the answer lands DURING the wait -- and the
dialog read on: the session never came back to its prompt, mail waiting
there was never delivered, and the next key typed into the pane went to
the question long since answered, where an ``a`` aborted the session.
Measured in a supervised run: a session stood 4 minutes with unread
mail behind a dialog that had been denied from outside.
"""

import io

from delfin.agent import repl
from delfin.agent import repl_render as rr
from delfin.agent import terminal_confirm as tc

PLAIN = rr.Theme(enabled=False)


class _AnsweredMidWait:
    """A reader whose question is answered from outside while it waits.

    The first reads are silent, as a waiting dialog sees them; on the
    third the answer lands, and the keys after it are the start of a
    line typed for the NEXT prompt -- an ``a`` among them.
    """

    active = True

    def __init__(self, req, later_keys, decision=False):
        self._req = req
        self._later = list(later_keys)
        self._decision = decision
        self.reads = 0

    def read_ready(self, timeout):
        self.reads += 1
        if self.reads < 3:
            return ""
        if self.reads == 3:
            self._req.decision = self._decision
            self._req.resolved = True
            return ""
        return self._later.pop(0) if self._later else ""

    @property
    def untouched(self):
        return list(self._later)


class _Engine:
    def __init__(self):
        self.token_usage = {"input": 0, "output": 0}
        self.stopped = False

    def get_status(self):
        return {}

    def request_stop(self):
        self.stopped = True


def _agent(broker):
    engine = _Engine()
    err = io.StringIO()
    agent = repl.TerminalAgent(
        engine, out=io.StringIO(), err=err, broker=broker)
    agent.transcript.theme = PLAIN
    return agent, engine, err


def test_a_confirm_answered_during_the_wait_takes_no_later_key():
    broker = tc.TerminalConfirmBroker(timeout_s=5.0)
    agent, engine, err = _agent(broker)
    req = tc.ConfirmRequest(
        kind=tc.CONFIRM, tool="bash",
        args={"command": "pytest -x tests/"}, preview="$ pytest -x tests/")
    keys = _AnsweredMidWait(req, ["a", "l", "s"])

    agent._answer_request(req, keys)

    assert keys.untouched == ["a", "l", "s"], (
        "the keys after the answer belong to the next prompt line")
    assert broker.aborted is False, "a later 'a' must not abort the session"
    assert engine.stopped is False
    assert req.decision is False, "the outside answer stands"
    assert "answered elsewhere" in err.getvalue()


def test_a_question_answered_during_the_wait_takes_no_later_key():
    broker = tc.TerminalConfirmBroker(timeout_s=5.0)
    agent, engine, err = _agent(broker)
    req = tc.ConfirmRequest(
        kind=tc.ASK, payload={"question": "Which one?",
                              "options": ["first", "second"]})
    outside = {"answers": ["second"]}
    keys = _AnsweredMidWait(req, ["1", "2"], decision=outside)

    agent._answer_question(req, keys)

    assert keys.untouched == ["1", "2"]
    assert req.decision == outside, "the outside answer stands"
    assert "answered elsewhere" in err.getvalue()
    assert "too late" not in err.getvalue(), (
        "nobody pressed a late key -- nothing arrived too late")


def test_a_dialog_nobody_answered_still_takes_its_key():
    """Control: without an outside answer the key is the answer."""
    broker = tc.TerminalConfirmBroker(timeout_s=5.0)
    agent, engine, err = _agent(broker)
    req = tc.ConfirmRequest(
        kind=tc.CONFIRM, tool="bash",
        args={"command": "ls"}, preview="$ ls")

    class _Keys:
        active = True

        def __init__(self):
            self._k = ["", "n"]

        def read_ready(self, timeout):
            return self._k.pop(0) if self._k else ""

    agent._answer_request(req, _Keys())

    assert req.resolved is True
    assert req.decision is False
    assert "refused" in err.getvalue()
