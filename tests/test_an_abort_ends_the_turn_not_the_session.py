"""The abort of one turn must not silence the questions of the next.

[a] abort refuses everything in flight and ends the turn. That refusal
posture is for THAT turn: a session that kept it for hours answered
every later question with a silent "user denied" no one had decided.
"""

import io
import threading

from delfin.agent import repl
from delfin.agent import terminal_confirm as tc


class EngineThatAsks:
    """stream_response asks one confirm question, like a gated tool."""

    def __init__(self, perms):
        self._perms = perms
        self.stop_requested = threading.Event()
        self.busy = False

    def stream_response(self, **kwargs):
        if self.busy:
            return "busy"
        self.busy = True
        try:
            # The tool call the model chose, gated by the broker.
            import delfin.agent.api_client as _A
            err = _A._doc_executor._run_permission_gate(
                "bash", {"command": "echo hello"}, self._perms)
            if err:
                return f"REFUSED: {err}"
            return "done"
        finally:
            self.busy = False

    def request_stop(self):
        self.stop_requested.set()

    def clear_stop(self):
        pass

    def get_status(self):
        return {}


def _agent_with_broker(tmp_path):
    import delfin.agent.api_client as _A
    broker = tc.TerminalConfirmBroker(timeout_s=5.0)
    perms = _A.KitToolPermissions(
        workspace=str(tmp_path), mode="default",
        confirm_callback=broker.callback)
    engine = EngineThatAsks(perms)
    agent = repl.TerminalAgent(engine, out=io.StringIO(), err=io.StringIO(),
                               broker=broker)
    return agent, broker


class TestTheAbortOfOneTurn:
    def test_a_question_in_the_next_turn_is_not_still_aborted(self, tmp_path):
        """THE case: one [a] abort, then a fresh turn asks again.

        Before the fix every question after the abort was refused in
        _enqueue with nobody deciding it.
        """
        agent, broker = _agent_with_broker(tmp_path)

        # One abort, as the terminal dialog would do it.
        broker.abort_all()
        assert broker.aborted is True

        # A question that arrives NOW -- same turn -- must still be
        # refused: the abort means "refuse everything that arrives
        # after" for THIS turn.
        inflight = broker._enqueue(tc.ConfirmRequest(kind=tc.CONFIRM))
        assert inflight.resolved is True, (
            "the abort must refuse what arrives in its own turn")
        # resolve it away so the broker is not left holding it
        broker.resolve(inflight, {"answers": []})

        # The NEXT turn resets the posture: its questions are asked.
        agent._begin_turn_resets_the_abort()

        req = broker._enqueue(tc.ConfirmRequest(kind=tc.CONFIRM))
        assert req.resolved is False, (
            "a question of the NEXT turn was silently refused by the "
            "abort of the PREVIOUS one")
        broker.resolve(req, {"answers": []})
