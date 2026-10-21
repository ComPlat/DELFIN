"""The turn-end wait registration, on the REAL repl hook.

The hook lives in repl.py and arrived as a protected-file patch built by
the operator; before that build this suite is red by design (no
_register_blocked_wait on TerminalAgent). It runs against the real
``delfin.agent.repl`` module — never a .gate copy — so what it proves is
what every session runs.

Contract (package B4 target 1 + operator rulings rounds 2-5):
- a supervised turn that ended blocked registers its wait with the
  on_whom the turn end actually knows (question/denial -> person, else
  operator);
- a person at the keyboard (no broker, or --confirm-timeout 0)
  registers nothing — the status footer and ops health stay unchanged
  for every ordinary user;
- a session nobody watches (a wave run in tmux, --confirm-timeout > 0)
  registers, so mutual waits are findable;
- the hook's own previous wait is cleared before the new turn is
  judged — the footer cannot stack;
- a registry failure (ValueError, OSError) is recorded in the session's
  failure log, never printed to the person's terminal; anything else
  propagates loudly.
"""
import pytest


@pytest.fixture
def agent(tmp_path):
    """A TerminalAgent on the real repl, with the hook's inputs stubbed.

    TerminalAgent.__new__ skips __init__: the hook needs only opts.cwd,
    the presence key and the broker. The RedirectedUserState keeps every
    registry/failure write inside the test's own throwaway state.
    """
    from delfin.agent import repl, state_paths

    class _Broker:
        """The one attribute the hook reads (terminal_confirm.py:641)."""

        def __init__(self, timeout_s: float):
            self.timeout_s = timeout_s

    with state_paths.RedirectedUserState(tmp_path):
        def make(presence="sess-hook", broker=None):
            a = repl.TerminalAgent.__new__(repl.TerminalAgent)
            a.opts = type("O", (), {"cwd": str(tmp_path)})()
            a._presence_key = lambda: presence
            a.broker = broker
            return a
        make.broker_cls = _Broker
        yield make


def test_no_marker_at_all_registers_nothing(agent):
    """A session with no broker at all registers nothing (the plain
    terminal, the pty test's shape)."""
    from delfin.agent import wait_registry
    a = agent(broker=None)
    a._register_blocked_wait("QUESTION: proceed?", denied=False)
    a._register_blocked_wait("plain block", denied=False)
    assert wait_registry.SessionWaits("sess-hook").list() == []


def test_a_person_at_the_keyboard_registers_nothing(agent):
    """cli.py:2114 builds a broker for EVERY person at a tty, but a
    person passes --confirm-timeout 0 — the documented marker
    (cli.py:2116). timeout 0: no wait, no 'waiting:' footer."""
    from delfin.agent import wait_registry
    a = agent(broker=agent.broker_cls(0.0))
    a._register_blocked_wait("QUESTION: proceed?", denied=False)
    assert wait_registry.SessionWaits("sess-hook").list() == []


def test_a_watched_unattended_session_registers(agent):
    """--confirm-timeout > 0: a session nobody watches (a wave run in
    tmux) registers, so mutual waits and ops health see it."""
    from delfin.agent import wait_registry
    a = agent(broker=agent.broker_cls(900.0))
    a._register_blocked_wait("QUESTION: proceed?", denied=False)
    w = wait_registry.SessionWaits("sess-hook").list()
    assert w and w[0]["on_whom"] == "person", w


def test_question_waits_on_the_person(agent):
    from delfin.agent import wait_registry
    a = agent(broker=agent.broker_cls(900.0))
    a._register_blocked_wait("QUESTION: which isomer?", denied=False)
    w = wait_registry.SessionWaits("sess-hook").list()[0]
    assert w["on_whom"] == "person"
    assert w["value"] == "answer"


def test_denial_waits_on_the_person(agent):
    from delfin.agent import wait_registry
    a = agent(broker=agent.broker_cls(900.0))
    a._register_blocked_wait("stopped.", denied=True)
    w = wait_registry.SessionWaits("sess-hook").list()[0]
    assert w["on_whom"] == "person"


def test_plain_block_waits_on_the_operator(agent):
    from delfin.agent import wait_registry
    a = agent(broker=agent.broker_cls(900.0))
    a._register_blocked_wait("waiting for the operator to answer",
                             denied=False)
    w = wait_registry.SessionWaits("sess-hook").list()[0]
    assert w["on_whom"] == "operator"


def test_the_previous_own_wait_is_cleared_not_stacked(agent):
    """The hook's own wait from the PREVIOUS turn is gone once this turn
    was judged: at most THIS turn's wait remains, on one key."""
    from delfin.agent import wait_registry
    a = agent(broker=agent.broker_cls(900.0))
    a._register_blocked_wait("QUESTION: proceed?", denied=False)
    assert wait_registry.SessionWaits("sess-hook").list()
    a._register_blocked_wait("done, all tasks closed", denied=False)
    rows = wait_registry.SessionWaits("sess-hook").list()
    assert not [w for w in rows if str(w.get("value")) == "answer"], rows
    assert len(rows) == 1 and str(rows[0]["value"]) == "operator", rows


def test_other_sessions_waits_are_not_touched(agent, tmp_path):
    """The clear is scoped to the hook's own footprint: a wait on another
    session's reply survives the hook on THIS session."""
    from delfin.agent import wait_registry
    wait_registry.register_once("message_from", "nacht-s20",
                                on_whom="session",
                                for_what="cross review", key="sess-hook")
    a = agent(broker=agent.broker_cls(900.0))
    a._register_blocked_wait("QUESTION: proceed?", denied=False)
    a._register_blocked_wait("done", denied=False)
    rows = wait_registry.SessionWaits("sess-hook").list()
    assert any(str(w.get("value")) == "nacht-s20" for w in rows), rows


def test_a_registry_error_goes_to_the_failure_log_not_the_terminal(
        agent, caplog):
    """A full disk / lost lock costs one failure-log record, not the
    turn loop — and never a line on the person's screen (the pty
    regression the logging.warning version caused)."""
    import delfin.agent.wait_registry as wr
    from delfin.agent import failure_log
    a = agent(broker=agent.broker_cls(900.0))
    orig = wr.register_once

    def boom(*args, **kwargs):
        raise OSError("disk full")

    wr.register_once = boom
    try:
        with caplog.at_level("WARNING"):
            a._register_blocked_wait("plain", denied=False)
    finally:
        wr.register_once = orig
    assert not caplog.records, caplog.records
    rows = failure_log.read_failures(last_n=5)
    assert any(r.get("tool") == "wait_registration" for r in rows), rows


def test_a_registry_value_error_is_also_recorded(agent, caplog):
    import delfin.agent.wait_registry as wr
    from delfin.agent import failure_log
    a = agent(broker=agent.broker_cls(900.0))
    orig = wr.register_once

    def boom(*args, **kwargs):
        raise ValueError("forced")

    wr.register_once = boom
    try:
        with caplog.at_level("WARNING"):
            a._register_blocked_wait("plain", denied=False)
    finally:
        wr.register_once = orig
    assert not caplog.records, caplog.records
    rows = failure_log.read_failures(last_n=5)
    assert any(r.get("tool") == "wait_registration" for r in rows), rows


def test_anything_else_propagates(agent):
    """The bug class we want loud: NameError/ImportError/TypeError are
    NOT swallowed by the hook's narrow except."""
    import delfin.agent.wait_registry as wr
    a = agent(broker=agent.broker_cls(900.0))
    orig = wr.register_once

    def boom(*args, **kwargs):
        raise NameError("a real bug")

    wr.register_once = boom
    try:
        with pytest.raises(NameError):
            a._register_blocked_wait("plain", denied=False)
    finally:
        wr.register_once = orig
