"""The restored read-ledger was thrown away before any guard read it.

``restore_state`` deliberately leaves the backend client alone -- only
conversation state is restored. The client mints
``_observed_files_session`` fresh at the start of its first turn. And
after every stream the engine did, unconditionally:

    self._last_observed_files = set(_obs_ledger or ())

A REPLACE. So on the first resumed turn the set that
``_restore_evidence`` had just brought back was overwritten by the
client's empty one, and nothing had read it in between. The grounding
guard then ran in its enforcing branch (the availability flag is
re-derived as True the moment the client has a ledger at all) against a
record of nothing, and flagged every ``file:line`` the answer restated
from the restored history as unsupported.

That is the exact false correction the evidence export was written to
prevent. The export was correct; the value it restored survived about
one turn's worth of code.

Two changes, because either alone leaves a hole: the client's session
ledger is SEEDED from the restored set on restore, so the client's own
readers see the history too; and the post-stream update is a UNION, so
what the turn observed is added to what the session already knew instead
of replacing it.
"""

from __future__ import annotations

from delfin.agent.engine import AgentEngine


class _Client:
    """Stands in for a backend client: keeps a session ledger, mints it
    lazily at the start of a turn exactly as the real one does."""

    def __init__(self, *, keeps_a_ledger: bool = True):
        self._keeps = keeps_a_ledger

    def begin_turn(self, observed=()):
        if self._keeps:
            if not hasattr(self, "_observed_files_session"):
                self._observed_files_session = set()
            self._observed_files_session |= set(observed)


def _engine(client=None) -> AgentEngine:
    eng = AgentEngine.__new__(AgentEngine)
    for spec in AgentEngine._SESSION_FIELDS:
        setattr(eng, spec.attr, spec.reset())
    eng.mode = "solo"
    eng.route = []
    eng.messages = []
    eng.session_id = ""
    eng.client = client
    return eng


def _saved_state() -> dict:
    src = _engine()
    src.session_id = "s1"
    src._last_observed_files = {"delfin/agent/engine.py", "README.md"}
    src._observed_ledger_available = True
    return src.export_state()


def _post_stream_update(eng: AgentEngine) -> None:
    """The engine's own post-stream ledger update, isolated. Mirrors the
    block in stream_response that runs after every turn."""
    _obs = getattr(eng.client, "_observed_files_session", None)
    _restored = bool(getattr(eng, "_observed_ledger_available", False))
    eng._observed_ledger_available = (_obs is not None) or _restored
    eng._last_observed_files = (
        set(getattr(eng, "_last_observed_files", None) or ()) | set(_obs or ()))


# ---------------------------------------------------------------------------
# The restored value must still be there when the guard looks
# ---------------------------------------------------------------------------

def test_the_first_resumed_turn_does_not_discard_the_restored_reads():
    client = _Client()
    eng = _engine(client)
    eng.restore_state(_saved_state())

    client.begin_turn(observed={"delfin/agent/cli.py"})
    _post_stream_update(eng)

    assert eng._last_observed_files == {
        "delfin/agent/engine.py", "README.md", "delfin/agent/cli.py"}, (
        "the turn's own reads replaced the session's instead of joining "
        "them -- every file:line restated from restored history is now "
        "unsupported")


def test_the_client_is_seeded_so_its_own_readers_see_the_history():
    client = _Client()
    eng = _engine(client)
    eng.restore_state(_saved_state())
    assert getattr(client, "_observed_files_session") == {
        "delfin/agent/engine.py", "README.md"}


def test_a_turn_that_reads_nothing_leaves_the_restored_record_intact():
    client = _Client()
    eng = _engine(client)
    eng.restore_state(_saved_state())
    client.begin_turn()
    _post_stream_update(eng)
    assert eng._last_observed_files == {"delfin/agent/engine.py", "README.md"}
    assert eng._observed_ledger_available is True


# ---------------------------------------------------------------------------
# ...without inventing a ledger that never existed
# ---------------------------------------------------------------------------

def test_a_save_with_no_ledger_does_not_give_the_client_one():
    """"the backend keeps no ledger" is a different fact from "the ledger
    is empty", and only the second one enforces. Seeding must not turn
    the first into the second."""
    src = _engine()
    src._observed_ledger_available = False
    client = _Client(keeps_a_ledger=False)
    eng = _engine(client)
    eng.restore_state(src.export_state())
    assert not hasattr(client, "_observed_files_session")
    assert eng._observed_ledger_available is False


def test_restoring_without_a_client_does_not_raise():
    eng = _engine(None)
    eng.restore_state(_saved_state())
    assert eng._last_observed_files == {"delfin/agent/engine.py", "README.md"}


def test_the_engine_updates_the_ledger_by_union_not_by_replacement():
    """Read off the source, because the failure was one operator: the
    post-stream update is the only place the restored value can be lost,
    and it is inside a method too large to exercise end to end here."""
    import inspect
    src = inspect.getsource(AgentEngine.stream_response)
    assert "self._last_observed_files = set(_obs_ledger or ())" not in src, (
        "the post-stream update replaces the session ledger again")
    assert "_carried | set(_obs_ledger or ())" in src
