"""Tests for auto-memory distillation (Roadmap B1).

Contract: opt-in (default off), LLM injectable, robust fact parsing,
dedupe against existing memories, short sessions skipped, never raises.
"""

from __future__ import annotations

import pytest

from delfin.agent import memory_distill as md


_CHAT = [
    {"role": "user", "content": "bitte immer englische strings im code"},
    {"role": "assistant", "content": "Verstanden — alle Strings englisch."},
    {"role": "user", "content": "und nutze azure.gpt-5.4, nicht gpt-5.5"},
    {"role": "assistant", "content": "OK."},
    {"role": "user", "content": "porpoise liegt unter /pfs/.../Porpoise"},
]


def test_default_settings_disabled():
    from delfin.user_settings import DEFAULT_SETTINGS
    cfg = DEFAULT_SETTINGS["agent"]["auto_memory"]
    assert cfg["enabled"] is False
    assert cfg["max_facts"] >= 1


def test_disabled_means_zero_llm_calls():
    calls = []
    n = md.distill_and_save(
        _CHAT,
        settings={"agent": {"auto_memory": {"enabled": False}}},
        llm_fn=lambda *a: calls.append(1) or "fact one two three",
    )
    assert n == 0 and calls == []           # opt-in respected, no tokens


def test_short_sessions_are_skipped():
    calls = []
    n = md.distill_and_save(
        [{"role": "user", "content": "hi"}],
        settings={"agent": {"auto_memory": {"enabled": True}}},
        llm_fn=lambda *a: calls.append(1) or "x",
    )
    assert n == 0 and calls == []


def test_parse_facts_filters_noise():
    raw = ("- User prefers English strings in code\n"
           "* azure.gpt-5.4 is the working model; gpt-5.5 is broken\n"
           "NONE\n"
           "ok\n"                                   # too short
           "Project Porpoise lives at /pfs/.../Porpoise\n")
    facts = md.parse_facts(raw, max_facts=5)
    assert len(facts) == 3
    assert facts[0].startswith("User prefers English")


def test_parse_none_returns_empty():
    assert md.parse_facts("NONE") == []
    assert md.parse_facts("") == []


def test_distill_saves_deduped_facts(monkeypatch):
    saved = []
    monkeypatch.setattr("delfin.agent.memory_store.load_memories",
                        lambda: [{"text": "already known fact here"}])
    monkeypatch.setattr("delfin.agent.memory_store.save_memory",
                        lambda text, source="": saved.append(text))
    n = md.distill_and_save(
        _CHAT,
        settings={"agent": {"auto_memory": {"enabled": True}}},
        llm_fn=lambda *a: ("already known fact here\n"
                           "User prefers English strings in code\n"),
    )
    assert n == 1                            # duplicate skipped
    assert saved == ["User prefers English strings in code"]


def test_force_overrides_optin(monkeypatch):
    monkeypatch.setattr("delfin.agent.memory_store.load_memories", lambda: [])
    saved = []
    monkeypatch.setattr("delfin.agent.memory_store.save_memory",
                        lambda text, source="": saved.append(text))
    n = md.distill_and_save(
        _CHAT,
        settings={"agent": {"auto_memory": {"enabled": False}}},
        llm_fn=lambda *a: "Porpoise project path is /pfs/.../Porpoise",
        force=True,                           # the manual /memorize
    )
    assert n == 1 and saved


def test_llm_failure_is_contained():
    def _boom(*a):
        raise RuntimeError("no api key")
    n = md.distill_and_save(
        _CHAT,
        settings={"agent": {"auto_memory": {"enabled": True}}},
        llm_fn=_boom,
    )
    assert n == 0                            # never raises


def test_excerpt_drops_tool_noise():
    chat = _CHAT + [{"role": "tool", "content": "HUGE TOOL OUTPUT " * 500}]
    ex = md._transcript_excerpt(chat)
    assert "HUGE TOOL OUTPUT" not in ex
    assert "USER:" in ex and "ASSISTANT:" in ex
