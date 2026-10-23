"""A resumed session must estimate its context the way the live one did.

The window bar reads ``max(local character estimate, provider floor)``.
The local estimate can only see what is in the message list plus the
system prompt; the tool schemas are sent as a separate argument and exist
only inside the provider's own count, which is why the floor matters --
the schema budget alone allows several thousand tokens the character
count can never see.

Both inputs used to be dropped on resume, so a loaded session read far
emptier than the conversation it held. That is not only cosmetic: the
same estimate drives the sliding-window trim and auto-compaction, so an
under-reading session does not trim when it should.

The prompt TEXT is deliberately not persisted -- it carries the injected
memory. Only its length is, because only its length is what the estimate
needs.
"""

from __future__ import annotations

import inspect
import pathlib

import pytest

from delfin.agent.engine import AgentEngine


def _bare_engine():
    eng = AgentEngine.__new__(AgentEngine)
    eng.mode = "solo"
    eng.route = ["solo_agent"]
    eng.current_role_index = 0
    eng.role_outputs = {}
    eng.compaction_summaries = {}
    eng.messages = []
    eng.token_usage = {"input": 0, "output": 0}
    eng.cost_usd = 0.0
    eng.session_id = ""
    eng._project_dir = ""
    eng._last_input_tokens = 0
    eng._system_prompt_chars = 0
    return eng


def test_the_prompt_size_survives_a_round_trip():
    eng = _bare_engine()
    eng._system_prompt_chars = 24_000
    state = eng.export_state()
    assert state["system_prompt_chars"] == 24_000

    eng2 = _bare_engine()
    eng2.restore_state(state)
    assert eng2._system_prompt_chars == 24_000


def test_the_prompt_text_itself_is_not_written_into_the_session():
    """It carries the injected memory and the estimate does not need it."""
    eng = _bare_engine()
    eng.last_system_prompt = "SECRET-MEMORY-CONTENT"
    eng._system_prompt_chars = len(eng.last_system_prompt)
    state = eng.export_state()
    assert "SECRET-MEMORY-CONTENT" not in repr(state)
    assert "last_system_prompt" not in state


def test_a_legacy_session_restores_clean():
    """Sessions saved before the field existed must not inherit a stale size."""
    eng = _bare_engine()
    eng._system_prompt_chars = 999
    eng.restore_state({"mode": "solo", "engine_messages": [], "session_id": "x"})
    assert eng._system_prompt_chars == 0


def test_export_survives_an_engine_that_never_ran_init():
    eng = AgentEngine.__new__(AgentEngine)
    for name, value in (("mode", "solo"), ("route", []), ("current_role_index", 0),
                        ("role_outputs", {}), ("compaction_summaries", {}),
                        ("messages", []), ("token_usage", {}), ("cost_usd", 0.0),
                        ("session_id", ""), ("_project_dir", ""),
                        ("_last_input_tokens", 0)):
        setattr(eng, name, value)
    assert eng.export_state()["system_prompt_chars"] == 0


def test_the_estimate_counts_the_restored_size():
    """The point of persisting it: the bar after a load matches the bar before."""
    eng = _bare_engine()
    eng.messages = [{"role": "user", "content": "x" * 4_000}]
    eng.last_system_prompt = "p" * 20_000
    eng._system_prompt_chars = 20_000
    live = eng._estimate_context_tokens()

    resumed = _bare_engine()
    resumed.restore_state(eng.export_state())
    # The text is gone by design -- the size is what carries the estimate.
    assert not getattr(resumed, "last_system_prompt", "")
    assert resumed._system_prompt_chars == 20_000
    assert resumed._estimate_context_tokens() == live, (
        'a resumed session estimates its context differently from the live '
        'one it was saved from')


def test_the_estimate_reads_the_persisted_size_not_the_text():
    """If the estimate went back to len(last_system_prompt), every resumed
    session would silently lose the prompt from its count again."""
    source = pathlib.Path(AgentEngine.__module__.replace(".", "/") + ".py")
    if not source.exists():
        source = pathlib.Path(inspect.getfile(AgentEngine))
    text = source.read_text(encoding="utf-8")
    assert '_sp_live or int(getattr(self, "_system_prompt_chars", 0) or 0)' in text, (
        'the context estimate no longer falls back to the persisted prompt '
        'size, so a resumed session counts the system prompt as nothing')


# ---------------------------------------------------------------------------
# The enumeration trap
# ---------------------------------------------------------------------------

def test_the_saver_accepts_every_field_the_exporter_produces():
    """export_state gaining a field is useless if save_session cannot take it.

    A named parameter is one way; a catch-all that persists the key
    verbatim is the other, and it is the one that cannot fall behind. What
    must never happen is a TypeError inside a best-effort save, which is a
    silently unwritten session.
    """
    from delfin.agent import session_store

    exported = set(_bare_engine().export_state())
    params = inspect.signature(session_store.save_session).parameters
    accepted = set(params)
    catch_all = any(p.kind is inspect.Parameter.VAR_KEYWORD
                    for p in params.values())
    missing = {k for k in exported if k not in accepted}
    assert catch_all or not missing, (
        f'save_session cannot store {sorted(missing)} -- the exporter '
        'produces it and the saver would drop it on the floor')


def test_headless_resume_forwards_the_whole_stored_session():
    """The CLI path listed its keys by hand and fell behind twice already."""
    source = pathlib.Path(inspect.getfile(
        __import__("delfin.agent.cli", fromlist=["x"]))).read_text(encoding="utf-8")
    assert "**data," in source, (
        'delfin/agent/cli.py enumerates restore fields again, so every new '
        'piece of session state has to be remembered in a second place')


# ---------------------------------------------------------------------------
# The evidence a resumed session judges itself against
# ---------------------------------------------------------------------------

def _engine_with_evidence():
    eng = _bare_engine()
    eng._last_observed_files = {"a.py", "b.xlsx"}
    eng._exec_commands_session = ["pytest -q", "python build.py"]
    eng._session_tool_names = {"subagent", "read_file"}
    eng._delegation_satisfied = True
    eng.role_verdicts = {"test_agent": {"verdict": "PASS"}}
    eng._trimmed_chars_since_floor = 40_000
    return eng


def test_the_ledgers_survive_a_round_trip():
    """A resumed session got the whole conversation and none of the
    evidence, while the "a ledger exists" flags came back True -- which is
    the ENFORCING branch, not the silent one."""
    resumed = _bare_engine()
    resumed.restore_state(_engine_with_evidence().export_state())

    assert resumed._last_observed_files == {"a.py", "b.xlsx"}
    assert resumed._exec_commands_session == ["pytest -q", "python build.py"]
    assert resumed._session_tool_names == {"subagent", "read_file"}
    assert resumed._delegation_satisfied is True
    assert resumed.role_verdicts == {"test_agent": {"verdict": "PASS"}}


def test_the_trim_credit_survives():
    """The sharpest of them: restore brings back _last_input_tokens, a
    PRE-trim provider snapshot, and the credit that offsets it reset to
    zero -- so a resumed session read its context as larger than it was and
    compacted early."""
    resumed = _bare_engine()
    resumed.restore_state(_engine_with_evidence().export_state())
    assert resumed._trimmed_chars_since_floor == 40_000


def test_a_session_saved_before_this_existed_still_loads():
    eng = _engine_with_evidence()
    eng.restore_state({"mode": "solo", "engine_messages": [], "session_id": "x"})
    assert eng._last_observed_files == set()
    assert eng._exec_commands_session == []
    assert eng._trimmed_chars_since_floor == 0
    assert eng._delegation_satisfied is False


def test_a_malformed_evidence_block_does_not_raise():
    for junk in ("not a dict", 42, [], None):
        eng = _bare_engine()
        eng.restore_state({"mode": "solo", "engine_messages": [],
                           "session_id": "x", "evidence": junk})
        assert eng._last_observed_files == set()


def test_the_saver_takes_the_new_field():
    """export_state gaining a key is useless if save_session cannot store
    it -- the exact break this file already guards against for the others."""
    import inspect
    from delfin.agent import session_store

    assert "evidence" in inspect.signature(session_store.save_session).parameters
