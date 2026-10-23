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


# ---------------------------------------------------------------------------
# ...through the real turn, not through a copy of it
# ---------------------------------------------------------------------------
# The union was asserted twice over source text -- a NEGATIVE substring
# ("...= set(_obs_ledger or ())" not in src) and a positive one -- with
# the reason "the post-stream update is inside a method too large to
# exercise end to end here". A negative substring is satisfied forever by
# renaming the local; and ``_post_stream_update`` above is a COPY of the
# engine's block living in this file, so the behavioural tests around it
# would keep passing against an engine that no longer does any of it.
#
# A real turn is not too large to drive: a stub client and the same
# minimal pack the other engine tests build.

def _real_engine(tmp_path, observed_during_turn):
    """An engine on a stub client that observes files while it streams."""
    import textwrap
    from unittest.mock import MagicMock, patch

    from delfin.agent.api_client import StreamEvent

    lite = tmp_path / "pack_lite"
    (lite / "modes").mkdir(parents=True)
    (lite / "modes" / "solo.md").write_text("# quick mode")
    (lite / "manifest.yaml").write_text(textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - session_manager
    """))

    client = MagicMock()

    def _stream(*a, **k):
        # The real backend opens every turn with this, and on the first
        # turn it BINDS a fresh session ledger -- discarding whatever
        # ``_seed_client_observed_files`` put there. Called for real
        # rather than modelled, so the test cannot drift away from it.
        from delfin.agent.api_client import _begin_observed_ledgers
        _begin_observed_ledgers(client)
        client._observed_files_session |= set(observed_during_turn)
        yield StreamEvent(type="text_delta", text="fertig")
        yield StreamEvent(type="message_delta", output_tokens=3, cost_usd=0.0)

    client.stream_message = MagicMock(side_effect=_stream)
    with patch("delfin.agent.engine.create_client", return_value=client):
        from delfin.agent.engine import AgentEngine as _AE
        eng = _AE(repo_dir=tmp_path, backend="cli", mode="quick",
                  pack_dir=tmp_path)
    return eng


def test_a_real_resumed_turn_keeps_what_the_session_had_read(tmp_path):
    """The failure was one operator, on the FIRST resumed turn, before any
    guard had read the restored set."""
    eng = _real_engine(tmp_path, {"delfin/agent/cli.py"})
    eng.restore_state(_saved_state())
    assert eng._last_observed_files == {
        "delfin/agent/engine.py", "README.md"}

    eng.stream_response("weiter wie besprochen")

    assert eng._last_observed_files == {
        "delfin/agent/engine.py", "README.md", "delfin/agent/cli.py"}, (
        "the turn's own reads replaced the session's instead of joining "
        "them -- every file:line restated from restored history becomes "
        "unsupported")


def test_a_real_turn_that_reads_nothing_keeps_the_restored_record(tmp_path):
    eng = _real_engine(tmp_path, set())
    eng.restore_state(_saved_state())
    eng.stream_response("nur eine Rückfrage")
    assert eng._last_observed_files == {
        "delfin/agent/engine.py", "README.md"}
    assert eng._observed_ledger_available is True


def test_a_fresh_session_still_records_what_the_turn_read(tmp_path):
    """The union must not turn into "never update"."""
    eng = _real_engine(tmp_path, {"delfin/agent/cli.py"})
    eng.stream_response("lies mal die cli")
    assert eng._last_observed_files == {"delfin/agent/cli.py"}
