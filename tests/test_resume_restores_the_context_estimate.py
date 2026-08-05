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
    """export_state gaining a field is useless if save_session cannot take it."""
    from delfin.agent import session_store

    exported = set(_bare_engine().export_state())
    accepted = set(inspect.signature(session_store.save_session).parameters)
    # route/mode/role_index and friends are named identically on both sides.
    missing = {k for k in exported if k not in accepted}
    assert not missing, (
        f'save_session cannot store {sorted(missing)} -- the exporter '
        'produces it and the saver would drop it on the floor')


def test_headless_resume_forwards_the_whole_stored_session():
    """The CLI path listed its keys by hand and fell behind twice already."""
    source = pathlib.Path(inspect.getfile(
        __import__("delfin.agent.cli", fromlist=["x"]))).read_text(encoding="utf-8")
    assert "**data," in source, (
        'delfin/agent/cli.py enumerates restore fields again, so every new '
        'piece of session state has to be remembered in a second place')
