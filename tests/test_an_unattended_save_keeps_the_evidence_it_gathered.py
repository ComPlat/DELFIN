"""The headless save wrote ``evidence: {}`` and nobody was watching.

``export_state`` produces the evidence ledgers, the compaction summaries,
the project-directory pin, the last input-token count and the
system-prompt size. ``delfin/agent/cli.py::_save_session`` listed the
fields it forwarded BY HAND -- eight of them, none of those -- so
``save_session`` defaulted every one of them away.

The restore side had already been fixed for this exact hazard, with a
comment in the same file saying so: "Listing the fields by hand was the
bug ... headless --resume came back thinner than the session that was
saved." The writer eighty lines below still listed them by hand.

What the resumed session then believed: it had read no files, run no
commands, used no tools, held no role verdict and no test evidence --
while ``_observed_ledger_available`` is re-derived as True on the first
turn, which is the ENFORCING branch of the grounding guard. It judged as
if it had checked. Every scheduled unattended run saves through this
path.

The fix is not "add five more keyword arguments". It is that the
exported state is the single source: the caller forwards it whole, and
``save_session`` has a catch-all so a key the exporter learns tomorrow
reaches disk today instead of raising TypeError inside a best-effort
save. These tests go through the CALLER, because the previous tests all
went through ``save_session`` directly -- which is why the call site
could be wrong for as long as it was.
"""

from __future__ import annotations

import json
import types

import pytest

from delfin.agent import cli as agent_cli
from delfin.agent import session_store as ss
from delfin.agent.engine import AgentEngine


@pytest.fixture()
def sessions_dir(tmp_path, monkeypatch):
    d = tmp_path / "sessions"
    d.mkdir()
    monkeypatch.setattr(ss, "_SESSIONS_DIR", d)
    monkeypatch.setattr(ss, "_ensure_dir", lambda: d)
    monkeypatch.setattr(ss, "acquire_session_lock", lambda *a, **k: True)
    return d


def _engine_that_did_some_work() -> AgentEngine:
    eng = AgentEngine.__new__(AgentEngine)
    for spec in AgentEngine._SESSION_FIELDS:
        setattr(eng, spec.attr, spec.reset())
    eng.mode = "solo"
    eng.route = ["solo_agent"]
    eng.session_id = "unattended-1"
    eng.messages = [
        {"role": "user", "content": "check the parser"},
        {"role": "assistant", "content": "delfin/agent/engine.py:2915 is it"},
    ]
    eng._last_observed_files = {"delfin/agent/engine.py"}
    eng._observed_ledger_available = True
    eng._exec_commands_session = ["python -m pytest -q"]
    eng._session_tool_names = {"read_file", "bash", "subagent"}
    eng._delegation_satisfied = True
    eng.role_verdicts = {"test_agent": {"verdict": "PASS"}}
    eng.role_test_evidence = {"test_agent": ["pytest -q: 40 passed"]}
    eng.compaction_summaries = {"solo_agent": "earlier: parser rewritten"}
    eng._project_dir = "/work/parser"
    eng._last_input_tokens = 31_000
    eng._system_prompt_chars = 22_500
    eng._trimmed_chars_since_floor = 4_096
    return eng


def _saved(sessions_dir, engine) -> dict:
    sid = agent_cli._save_session(engine, sessions_dir)
    return json.loads((sessions_dir / f"{sid}.json").read_text())


# ---------------------------------------------------------------------------
# What the unattended save actually wrote
# ---------------------------------------------------------------------------

def test_the_evidence_reaches_the_file(sessions_dir):
    """It was written as ``{}`` -- the whole record of what the session
    had checked, replaced by the parameter default."""
    ev = _saved(sessions_dir, _engine_that_did_some_work())["evidence"]
    assert ev["observed_files"] == ["delfin/agent/engine.py"]
    assert ev["exec_commands"] == ["python -m pytest -q"]
    assert ev["observed_ledger_available"] is True
    assert ev["role_test_evidence"] == {"test_agent": ["pytest -q: 40 passed"]}


@pytest.mark.parametrize("key,expected", [
    ("project_dir", "/work/parser"),
    ("last_input_tokens", 31_000),
    ("system_prompt_chars", 22_500),
    ("compaction_summaries", {"solo_agent": "earlier: parser rewritten"}),
])
def test_the_other_dropped_fields_reach_the_file(sessions_dir, key, expected):
    assert _saved(sessions_dir, _engine_that_did_some_work())[key] == expected


def test_every_exported_key_reaches_the_file(sessions_dir):
    """The mechanism, not the symptom: no key of the exported state may be
    lost between the exporter and the file, whatever it is called."""
    engine = _engine_that_did_some_work()
    exported = engine.export_state()
    on_disk = _saved(sessions_dir, engine)
    lost = [k for k in exported if k not in on_disk]
    assert not lost, f"the save path dropped {lost}"


def test_a_key_the_saver_has_never_heard_of_still_reaches_the_file(
        sessions_dir, monkeypatch):
    """The guarantee that makes the enumeration unnecessary. A field added
    to the exporter must not need a matching keyword argument before it
    can be persisted -- and must never raise TypeError inside a
    best-effort save, which would silently write no session at all."""
    engine = _engine_that_did_some_work()
    real_export = engine.export_state

    def _export_with_a_new_field():
        state = real_export()
        state["a_field_invented_after_the_saver_was_written"] = 7
        return state

    engine.export_state = _export_with_a_new_field
    data = _saved(sessions_dir, engine)
    assert data["a_field_invented_after_the_saver_was_written"] == 7


# ---------------------------------------------------------------------------
# ...and back in through the caller, which is where the round trip broke
# ---------------------------------------------------------------------------

def test_the_headless_round_trip_restores_what_it_saved(sessions_dir):
    """Through _save_session and _resume_or_create, not through
    save_session/restore_state. Both halves were individually tested and
    the pair was still broken."""
    saved = _engine_that_did_some_work()
    sid = agent_cli._save_session(saved, sessions_dir)

    resumed = AgentEngine.__new__(AgentEngine)
    for spec in AgentEngine._SESSION_FIELDS:
        setattr(resumed, spec.attr, spec.reset())
    resumed.mode = "solo"
    resumed.route = []
    resumed.messages = []
    resumed.session_id = "throwaway"
    resumed.client = None

    args = types.SimpleNamespace(session=sid, mode="solo", verbose=False)
    assert agent_cli._resume_or_create(resumed, args) == sid

    assert resumed._last_observed_files == {"delfin/agent/engine.py"}
    assert resumed._exec_commands_session == ["python -m pytest -q"]
    assert resumed._observed_ledger_available is True
    assert resumed.role_test_evidence == {
        "test_agent": ["pytest -q: 40 passed"]}
    assert resumed.compaction_summaries == {
        "solo_agent": "earlier: parser rewritten"}
    assert resumed._project_dir == "/work/parser"
    assert resumed._last_input_tokens == 31_000
    assert resumed._system_prompt_chars == 22_500
    assert resumed._trimmed_chars_since_floor == 4_096


def test_the_enforcing_flag_is_never_true_over_an_empty_record(sessions_dir):
    """The precise harm. A ledger that exists but is empty means 'nothing
    was read', and the guard enforces on it. Saving the flag while
    dropping the files is the one combination that must be impossible."""
    ev = _saved(sessions_dir, _engine_that_did_some_work())["evidence"]
    if ev.get("observed_ledger_available"):
        assert ev.get("observed_files"), (
            "the save says a read-ledger existed and carries nothing in "
            "it -- the resumed session would enforce against an empty "
            "record it did not earn")
