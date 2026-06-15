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


# ---------------------------------------------------------------------------
# Typed project-memory store (repo_root) — "memory like Claude Code"
# ---------------------------------------------------------------------------


def test_save_facts_writes_typed_store(monkeypatch, tmp_path):
    from pathlib import Path
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    from delfin.agent import memory_store as ms
    repo = tmp_path / "repo"; repo.mkdir()
    n = md.save_facts(
        ["feedback: never add a Co-Authored-By trailer to commits",
         "project: shipping the memory layer this week"],
        repo_root=repo,
    )
    assert n == 2
    mdir = ms._claude_memory_dir(repo)
    files = sorted(p.name for p in mdir.glob("*.md"))
    assert "MEMORY.md" in files
    assert any(f.startswith("feedback_") for f in files)
    assert any(f.startswith("project_") for f in files)
    index = (mdir / "MEMORY.md").read_text(encoding="utf-8")
    assert "Feedback" in index and "Project" in index


def test_save_facts_typed_dedups_on_body(monkeypatch, tmp_path):
    from pathlib import Path
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    repo = tmp_path / "repo"; repo.mkdir()
    assert md.save_facts(["user: Max is a quantum chemist at KIT"],
                         repo_root=repo) == 1
    # Same fact again (even without the prefix) must be skipped.
    assert md.save_facts(["Max is a quantum chemist at KIT"],
                         repo_root=repo) == 0


def test_distill_and_save_threads_repo_root(monkeypatch, tmp_path):
    from pathlib import Path
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    from delfin.agent import memory_store as ms
    repo = tmp_path / "repo"; repo.mkdir()
    n = md.distill_and_save(
        _CHAT,
        settings={"agent": {"auto_memory": {"enabled": True}}},
        llm_fn=lambda *a: "feedback: always use English strings in code",
        repo_root=repo,
    )
    assert n == 1
    mdir = ms._claude_memory_dir(repo)
    assert any(p.name.startswith("feedback_") for p in mdir.glob("*.md"))
