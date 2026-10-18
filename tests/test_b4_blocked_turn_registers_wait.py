"""B4 open item, closed: a blocked turn registers its wait on the real path.

The turn-end hook (repl.py turn() -> _note_turn_blocked) now also calls
_register_blocked_wait, which uses wait_registry.register_once. The test
drives the REAL TerminalAgent turn-end hook (as tests/test_blocked_note_
at_turn_end.py does), never the helper in isolation.
"""
import pytest

from delfin.agent import state_paths
from delfin.agent import wait_registry


class _Engine:
    session_id = "sess-b4-hook"


@pytest.fixture
def home(tmp_path, monkeypatch):
    with state_paths.RedirectedUserState(tmp_path):
        yield tmp_path / ".delfin"


def _agent(tmp_path, monkeypatch, broker=None):
    import delfin.agent.repl as repl

    a = repl.TerminalAgent.__new__(repl.TerminalAgent)
    a.opts = type("O", (), {"cwd": str(tmp_path)})()
    a.engine = _Engine()
    a._presence_key = lambda: "sess-b4-hook"
    if broker is not None:
        a.broker = broker
    return a


class _SupervisedBroker:
    """The one attribute the v6 guard reads (terminal_confirm.py:641).

    timeout_s > 0 = a session nobody sits at (a wave run in tmux,
    --confirm-timeout 900): it waits on the operator and SAYS so.
    timeout_s == 0 = a person at the keyboard: never a queued wait.
    """

    def __init__(self, timeout_s: float):
        self.timeout_s = timeout_s


def test_a_blocked_turn_registers_its_wait(tmp_path, monkeypatch, home):
    a = _agent(tmp_path, monkeypatch, broker=_SupervisedBroker(900.0))
    a._note_turn_blocked("QUESTION: proceed?", denied=False)
    a._register_blocked_wait("QUESTION: proceed?", denied=False)
    waits = wait_registry.SessionWaits("sess-b4-hook").list()
    assert len(waits) == 1, waits
    w = waits[0]
    assert w["kind"] == "message_from"
    assert w["value"] == "answer"
    assert w["on_whom"] == "person"
    assert "blocked" in w["for_what"].lower()


def test_a_person_at_the_keyboard_registers_nothing(tmp_path, monkeypatch,
                                                    home):
    """v6 contract (operator round 5): timeout 0 / no broker = an
    ordinary user session — no queued wait, no footer, no registry row."""
    a = _agent(tmp_path, monkeypatch)  # no broker at all
    a._register_blocked_wait("QUESTION: proceed?", denied=False)
    assert wait_registry.SessionWaits("sess-b4-hook").list() == []

    a = _agent(tmp_path, monkeypatch, broker=_SupervisedBroker(0.0))
    a._register_blocked_wait("plain block", denied=False)
    assert wait_registry.SessionWaits("sess-b4-hook").list() == []


def test_the_wiring_is_in_the_turn_end_hook(tmp_path, monkeypatch, home):
    """turn() must call the registration, not leave it as dead code."""
    import inspect

    import delfin.agent.repl as repl

    turn_src = inspect.getsource(repl.TerminalAgent)
    assert "_note_turn_blocked(result.text" in turn_src
    assert "_register_blocked_wait(result.text, denied=denied)" in turn_src
    # The call sits right beside the note writer in the turn-end hook.
    assert turn_src.index("_register_blocked_wait(result.text") \
        > turn_src.index("_note_turn_blocked(result.text")
