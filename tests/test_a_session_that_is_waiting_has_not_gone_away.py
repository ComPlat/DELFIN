"""A session waiting on an approval is still there, and still says so.

Presence was refreshed from the idle poll, which is where the job
wake-up and the operator's inbox are read. A session that goes from one
approval straight into the next never reaches it, and after
``_STALE_S`` -- fifteen minutes -- it drops out of ``open_sessions``.

Seen while supervising, 2026-09-20: five sessions had been running
twenty minutes, two were blocked on questions of mine at that very
moment, and `delfin-agent sessions` printed no "Open now" block at all.
The overview went quiet exactly when a supervisor most needs it, and it
went quiet about the sessions that were WAITING FOR THE SUPERVISOR.

Asking is proof of life, so asking is where the heartbeat belongs. The
broker does not learn about presence for it: it calls back, and the REPL
decides what that means.
"""

from __future__ import annotations

import pytest

from delfin.agent import terminal_confirm as tc


@pytest.fixture
def room(tmp_path, monkeypatch):
    monkeypatch.setattr(tc, "_PENDING_DIR", tmp_path / "pending")
    return tmp_path


def _req():
    return tc.ConfirmRequest(kind=tc.CONFIRM, tool="bash",
                             args={"command": "ls"}, preview="x")


class TestTheCallback:
    def test_asking_reports_activity(self, room):
        beats = []
        b = tc.TerminalConfirmBroker(session_key="s1",
                                     on_activity=lambda: beats.append(1))
        b._enqueue(_req())
        assert beats == [1]

    def test_every_question_is_a_beat(self, room):
        beats = []
        b = tc.TerminalConfirmBroker(session_key="s1",
                                     on_activity=lambda: beats.append(1))
        b._enqueue(_req())
        b._enqueue(_req())
        assert len(beats) == 2

    def test_a_broker_without_one_still_works(self, room):
        b = tc.TerminalConfirmBroker(session_key="s1")
        req = b._enqueue(_req())
        assert req.resolved is False

    def test_a_throwing_callback_does_not_cost_the_question(self, room):
        def boom():
            raise OSError("presence is on a dead mount")
        b = tc.TerminalConfirmBroker(session_key="s1", on_activity=boom)
        b._enqueue(_req())                      # must not raise
        assert len(tc.pending_at_terminals()) == 1


class TestTheReplWiresIt:
    def test_the_agent_hands_the_broker_its_heartbeat(self):
        from delfin.agent.repl import TerminalAgent
        agent = object.__new__(TerminalAgent)
        broker = tc.TerminalConfirmBroker(session_key="s1")
        agent.broker = broker
        TerminalAgent._arm_presence_heartbeat(agent)
        assert broker.on_activity == agent._announce_presence

    def test_no_broker_is_not_an_error(self):
        from delfin.agent.repl import TerminalAgent
        agent = object.__new__(TerminalAgent)
        agent.broker = None
        TerminalAgent._arm_presence_heartbeat(agent)   # must not raise
