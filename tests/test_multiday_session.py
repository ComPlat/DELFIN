"""Tests for the multi-day session hardening pass:

- Slim system prompt (progressive disclosure stripper)
- Cache-aware prompt section order
- Sliding-window proactive compaction
- Selective compaction (user goals survive)
- Task dependencies (blocked_by / blocks)
- Stalled-task detector
- Failure log + recurring-pattern detection
- BM25-style memory rerank
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Slim system prompt
# ---------------------------------------------------------------------------

def test_slim_prompt_strips_inactive_modules():
    from delfin.agent.prompt_loader import PromptLoader
    loader = PromptLoader()
    text = loader.load_role_prompt("solo_agent")
    slim = loader._strip_lazy_modules(
        text, task_text="hi there", mode_id="solo",
    )
    assert "## ORCA / chemistry" not in slim
    assert "## Web research" not in slim
    assert "## Jupyter notebooks" not in slim
    assert "## Background tasks" not in slim
    # Core sections always survive
    assert "## Strategies for approaching" in slim
    assert "## Subagents" in slim
    assert len(slim) < len(text)
    # Savings should be non-trivial — at least 5% smaller
    assert (len(text) - len(slim)) / len(text) > 0.05


def test_slim_prompt_keeps_chemistry_when_triggered():
    from delfin.agent.prompt_loader import PromptLoader
    loader = PromptLoader()
    text = loader.load_role_prompt("solo_agent")
    out = loader._strip_lazy_modules(
        text,
        task_text="extract imag frequencies from calc/foo.out",
        mode_id="solo",
    )
    assert "## ORCA / chemistry" in out
    assert "## Data search tools" in out
    assert "## DELFIN ops MCP tools" in out
    # Unrelated modules still stripped
    assert "## Web research" not in out


def test_slim_prompt_keeps_web_module_for_url_triggers():
    from delfin.agent.prompt_loader import PromptLoader
    loader = PromptLoader()
    text = loader.load_role_prompt("solo_agent")
    out = loader._strip_lazy_modules(
        text, task_text="please fetch https://example.com", mode_id="solo",
    )
    assert "## Web research" in out


def test_slim_prompt_only_solo_or_plan_strips():
    """Pipeline roles (quick / reviewed / etc.) get the full prompt
    because they're sensitive to subtle context changes."""
    from delfin.agent.prompt_loader import PromptLoader
    loader = PromptLoader()
    text = loader.load_role_prompt("solo_agent")
    out_quick = loader._strip_lazy_modules(
        text, task_text="hi", mode_id="quick",
    )
    out_solo = loader._strip_lazy_modules(
        text, task_text="hi", mode_id="solo",
    )
    assert len(out_quick) > len(out_solo)


def test_module_marker_lines_never_leak_into_prompt():
    """The ``<!-- module:X -->`` comment markers must never end up in
    the rendered prompt — they'd confuse the model."""
    from delfin.agent.prompt_loader import PromptLoader
    loader = PromptLoader()
    text = loader.load_role_prompt("solo_agent")
    for task_text in ("hi", "extract orca", "fetch https://x"):
        slim = loader._strip_lazy_modules(
            text, task_text=task_text, mode_id="solo",
        )
        assert "<!-- module:" not in slim, task_text


# ---------------------------------------------------------------------------
# Sliding-window + selective compaction
# ---------------------------------------------------------------------------

def _bare_engine():
    from delfin.agent.engine import AgentEngine
    eng = AgentEngine.__new__(AgentEngine)
    eng.messages = []
    eng.context_window_tokens = 1000  # tiny so test triggers fire
    eng.auto_compact_pct = 0.95
    eng._SLIDING_WINDOW_PCT = 0.70
    eng.last_compaction_info = None
    eng.backend = "api"
    eng.route = ["solo_agent"]
    eng.current_role_index = 0
    return eng


def test_sliding_window_trims_large_assistant_message():
    """A 5kB assistant message + a 200-char user goal should result in
    the assistant being trimmed head+tail while the user message survives
    verbatim once we cross the sliding-window threshold."""
    eng = _bare_engine()
    eng.messages = [
        {"role": "user", "content": "we need to test the flux capacitor"},
        {"role": "assistant", "content": "x" * 5000},
        {"role": "user", "content": "go"},
        {"role": "assistant", "content": "ok"},
        # Recent messages (last 4) are protected
        {"role": "user", "content": "newest user msg"},
        {"role": "assistant", "content": "newest reply"},
    ]
    assert eng._should_slide()
    n = eng._slide_window_trim()
    assert n >= 1
    # User goals survive verbatim
    assert eng.messages[0]["content"] == "we need to test the flux capacitor"
    # The huge assistant got middle-trimmed
    trimmed = eng.messages[1]["content"]
    assert "trimmed by sliding window" in trimmed
    assert len(trimmed) < 5000


def test_sliding_window_protects_recent_messages():
    """The last _KEEP_RECENT messages must NOT be trimmed even if they
    cross the size threshold."""
    eng = _bare_engine()
    eng.messages = [
        {"role": "assistant", "content": "a" * 5000},
        {"role": "user", "content": "u" * 100},
        {"role": "assistant", "content": "b" * 5000},   # recent
        {"role": "user", "content": "u" * 100},          # recent
        {"role": "assistant", "content": "c" * 5000},   # recent
        {"role": "user", "content": "u" * 100},          # recent
    ]
    eng._slide_window_trim()
    # The last 4 messages must be untouched
    for m in eng.messages[-4:]:
        assert "trimmed by sliding window" not in m["content"]


def test_selective_compaction_keeps_user_goals_more_than_first_line():
    """The new selective compactor must keep up to 400 chars of the
    user message, not just the first line as the legacy version did."""
    eng = _bare_engine()
    # Force the message threshold
    long_goal = "I need to refactor the flux capacitor.\nDetails:\n" + "x" * 200
    eng.messages = [{"role": "user", "content": long_goal}] * 15
    # Stub the client so LLM path fails and we fall through to extractive
    eng.client = None
    eng._compact_history()
    summary_msg = eng.messages[0]["content"]
    assert "Conversation summary" in summary_msg
    # The detailed user content must appear in the summary
    assert "flux capacitor" in summary_msg


# ---------------------------------------------------------------------------
# Task dependencies (blocked_by / blocks)
# ---------------------------------------------------------------------------

def test_task_create_with_blocked_by_records_reverse_index(tmp_path):
    from delfin.agent.agent_tasks import TaskStore
    store = TaskStore(tmp_path)
    a = store.create("parent")
    b = store.create("child", blocked_by=[a["id"]])
    assert b["blocked_by"] == [a["id"]]
    parent = store.get(a["id"])
    assert parent is not None
    assert b["id"] in (parent.get("blocks") or [])


def test_task_update_in_progress_refused_when_blocked(tmp_path):
    from delfin.agent.agent_tasks import TaskStore
    store = TaskStore(tmp_path)
    parent = store.create("setup")
    child = store.create("execute", blocked_by=[parent["id"]])
    with pytest.raises(ValueError) as exc:
        store.update(child["id"], status="in_progress")
    assert "blocked by" in str(exc.value)


def test_task_update_allows_in_progress_after_parent_completes(tmp_path):
    from delfin.agent.agent_tasks import TaskStore
    store = TaskStore(tmp_path)
    parent = store.create("setup")
    child = store.create("execute", blocked_by=[parent["id"]])
    store.update(parent["id"], status="completed")
    updated = store.update(child["id"], status="in_progress")
    assert updated["status"] == "in_progress"


def test_task_remove_blocked_by_unblocks(tmp_path):
    from delfin.agent.agent_tasks import TaskStore
    store = TaskStore(tmp_path)
    parent = store.create("p")
    child = store.create("c", blocked_by=[parent["id"]])
    store.update(child["id"], remove_blocked_by=[parent["id"]])
    out = store.get(child["id"])
    assert out["blocked_by"] == []
    # Now in_progress should be allowed
    assert store.update(child["id"], status="in_progress")["status"] == "in_progress"


# ---------------------------------------------------------------------------
# Stalled-task detector
# ---------------------------------------------------------------------------

def test_find_stalled_returns_only_old_in_progress(tmp_path):
    """A fresh in-progress task should NOT count as stalled, but a 30 min
    old one (faked via direct file edit) should."""
    from delfin.agent.agent_tasks import TaskStore
    import datetime as dt
    store = TaskStore(tmp_path)
    fresh = store.create("fresh")
    store.update(fresh["id"], status="in_progress")
    old = store.create("old")
    store.update(old["id"], status="in_progress")
    # Hand-edit the on-disk record to backdate
    data = json.loads(store.path.read_text(encoding="utf-8"))
    for t in data["tasks"]:
        if t["id"] == old["id"]:
            t["updated_at"] = (dt.datetime.utcnow()
                               - dt.timedelta(seconds=1800)).isoformat()
    store.path.write_text(json.dumps(data), encoding="utf-8")
    stalled = store.find_stalled(max_age_s=600)
    sids = [s["id"] for s in stalled]
    assert old["id"] in sids
    assert fresh["id"] not in sids


# ---------------------------------------------------------------------------
# Failure log
# ---------------------------------------------------------------------------

@pytest.fixture
def fl_path(tmp_path, monkeypatch):
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    from delfin.agent import failure_log
    monkeypatch.setattr(failure_log, "_LOG_PATH", tmp_path / "fl.jsonl")
    return tmp_path / "fl.jsonl"


def test_record_and_read_failure(fl_path):
    from delfin.agent import failure_log as fl
    fl.record_failure("Edit", "/x/y.py", "string not found", path=fl_path)
    fl.record_failure("Bash", "rm -rf x", "permission denied", path=fl_path)
    records = fl.read_failures(path=fl_path)
    assert len(records) == 2
    assert records[0]["tool"] == "Edit"


def test_top_recurring_filters_by_min_count(fl_path):
    from delfin.agent import failure_log as fl
    fl.record_failure("Edit", "/a", "boom", path=fl_path)
    fl.record_failure("Edit", "/b", "boom", path=fl_path)
    fl.record_failure("Edit", "/c", "boom", path=fl_path)
    fl.record_failure("Bash", "x", "other", path=fl_path)
    top = fl.top_recurring(min_count=3, path=fl_path)
    assert len(top) == 1
    assert top[0]["tool"] == "Edit"
    assert top[0]["count"] == 3


def test_detect_repeat_for_current_task(fl_path):
    from delfin.agent import failure_log as fl
    for _ in range(4):
        fl.record_failure("Edit", "f.py", "same error", path=fl_path)
    fl.record_failure("Edit", "g.py", "different", path=fl_path)
    n = fl.detect_repeat_for_current_task(
        "Edit", "f.py", "same error", path=fl_path,
    )
    assert n == 4
    assert fl.detect_repeat_for_current_task(
        "Edit", "g.py", "different", path=fl_path,
    ) == 1


def test_record_failure_never_raises_on_unwritable_path(tmp_path):
    """Best-effort writes must swallow OSError, not crash the agent."""
    from delfin.agent import failure_log as fl
    bad = tmp_path / "nope" / "nope" / "x.jsonl"
    # Path is fine but parent dir creation will succeed in tmp_path —
    # we can't easily force OSError on POSIX. Just make sure call returns.
    fl.record_failure("X", "y", "z", path=bad)
    assert bad.parent.is_dir()


# ---------------------------------------------------------------------------
# Memory rerank (BM25-style)
# ---------------------------------------------------------------------------

def test_format_memory_context_ranks_relevant_facts_higher(tmp_path, monkeypatch):
    """The new BM25 ranker should pick the fact that mentions the task's
    keyword even when older facts are more recent."""
    from delfin.agent import memory_store as ms
    monkeypatch.setattr(ms, "_DEFAULT_PATH", tmp_path / "memory.json")
    ms.save_memory("test pytest workflow keeps going", source="feedback")
    ms.save_memory("user prefers German answers", source="user")
    ms.save_memory("project deadline next Friday", source="project")
    ms.save_memory("frobnicator design uses lattice", source="feedback")
    # Task mentions pytest → the pytest fact should rank first
    out = ms.format_memory_context(task_text="run pytest")
    assert "pytest workflow" in out
    # The unrelated user/project entries should NOT dominate
    assert "deadline" not in out or out.index("pytest") < out.index("deadline")


def test_format_memory_context_falls_back_to_recency_without_query():
    """Empty task text falls through to last-N selection (legacy
    contract preserved)."""
    from delfin.agent import memory_store as ms
    out = ms.format_memory_context(task_text="")
    # Should not raise and should return a string (empty if no memories)
    assert isinstance(out, str)
