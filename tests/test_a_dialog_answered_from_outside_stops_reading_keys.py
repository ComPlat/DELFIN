"""A dialog ends when its question is answered, not when a key arrives.

Answered from outside (approvals approve/deny), the dialog in the pane
kept reading keys. A text typed later -- a command for the NEXT turn --
was swallowed as [a] abort of a question long since answered, and the
turn it belonged to never ran.
"""

import io

from delfin.agent import repl
from delfin.agent import repl_render as rr
from delfin.agent import terminal_confirm as tc

PLAIN = rr.Theme(enabled=False)


class _Keys:
    """Stands in for the raw-mode reader: hands over scripted keystrokes.

    read_ready never blocks long (0.1 s slices, like the real one), so
    a dialog that keeps polling after its question was answered keeps
    arriving at the same scripted keys -- which is exactly what the
    test wants to see it NOT take.
    """

    active = True

    def __init__(self, keys):
        self._keys = list(keys)
        self.reads = 0

    def read_ready(self, timeout):
        self.reads += 1
        return self._keys.pop(0) if self._keys else ""


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
    agent = repl.TerminalAgent(
        engine, out=io.StringIO(), err=io.StringIO(), broker=broker)
    agent.transcript.theme = PLAIN
    return agent


def test_a_dialog_answered_from_outside_stops_reading_keys():
    """THE case: the question is resolved while the dialog waits.

    The dialog must return as resolved -- without reading the keys that
    arrive later, which belong to the next prompt line, not to it.
    """
    broker = tc.TerminalConfirmBroker(timeout_s=5.0)
    agent = _agent(broker)
    req = tc.ConfirmRequest(
        kind=tc.CONFIRM, tool="bash",
        args={"command": "pytest -x tests/"}, preview="$ pytest -x tests/")

    # The question is answered from outside, immediately.
    req.resolved = True
    req.decision = False

    # Keys typed LATER -- a whole command line for the next turn. None
    # of these belongs to the dialog any more.
    keys = _Keys(["p", "y", "t", "e", "s", "t", "\r"])

    agent._answer_request(req, keys)

    assert req.resolved is True
    # The dialog must not have consumed the later keystrokes: they stay
    # in the reader for the prompt line that follows.
    assert keys._keys == ["p", "y", "t", "e", "s", "t", "\r"], (
        "a dialog whose question was answered from outside consumed "
        "keystrokes typed later, and a typed text acted as [a] abort "
        "of a question that was already answered")
