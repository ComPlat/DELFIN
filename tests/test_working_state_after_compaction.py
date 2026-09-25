"""A compacted agent must still know where it stands.

The compaction path guarantees (test_compaction_continuity.py) that the
user's GOAL survives verbatim and that the full transcript is archived.
Nothing guarantees the WORKING STATE -- what a person wants first after
a pause:

* which files this session changed (and whether anything was committed),
* the last test outcome per test file (red/green, when, after which change),
* open tasks from the task tool,
* the operator's instructions and denial REASONS ("not X, because Y"),
* the file/function names the session was last working on.

These tests build an artificial session with all of the above, run
``_compact_history`` with the model call stubbed (CLI/extractive path),
and assert the surviving context still answers those questions.
They are the control for the working-state block: red before it exists.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent.engine import AgentEngine


def _bare_engine(tmp_path: Path) -> AgentEngine:
    """An engine via __new__, exactly like the re-compaction tests in
    test_engine_compaction.py: no pack tree, just the attributes the
    compaction path touches."""
    eng = AgentEngine.__new__(AgentEngine)
    eng.messages = []
    eng.context_window_tokens = 10      # tiny window -> compaction fires
    eng.auto_compact_pct = 0.95
    eng.backend = "cli"                 # extractive path, no LLM call
    eng.client = None
    eng.last_compaction_info = {}
    eng.session_id = ""
    eng.repo_dir = tmp_path
    return eng


def _history_with_working_state() -> list[dict]:
    """A session that did real work: reads, an edit, a test run, a denial
    with a reason, and an operator instruction."""
    msgs: list[dict] = []
    msgs.append({"role": "user", "content":
                 "Fix the gate in module foo so lint stops failing."})
    # A machine turn: the agent edited a file and ran its tests.
    msgs.append({"role": "user", "content":
                 "[Command results]\n"
                 "edit_file delfin/foo.py: replaced the gate body\n"
                 'bash: gate tests/test_foo.py -> 12 passed'})
    msgs.append({"role": "assistant", "content":
                 "Edited delfin/foo.py (_check_path) and ran the tests; "
                 "all 12 green."})
    # A denial with a reason -- the shape the audit log feeds back.
    msgs.append({"role": "user", "content":
                 '[Command results]\nbash: git push origin main -> '
                 '{"error": "denied: push to the default branch is '
                 'not allowed for agent sessions"}'})
    msgs.append({"role": "assistant", "content":
                 "Understood, no push -- I commit on the branch instead."})
    # An operator instruction the agent must keep obeying after compaction.
    msgs.append({"role": "user", "content":
                 "[System] Operator: never edit api_client.py, it is "
                 "security code; requests go through session_message."})
    # Enough filler that the compaction floor lets these go.
    for i in range(6):
        msgs.append({"role": "user", "content": f"step {i}: " + "work " * 20})
        msgs.append({"role": "assistant", "content": f"did {i}: " + "done " * 20})
    return msgs


def test_working_state_survives_compaction(tmp_path):
    """Files changed, test outcome, denial reason, operator instruction and
    last-touched names must all still be readable in the post-compaction
    context."""
    eng = _bare_engine(tmp_path)
    eng.messages = _history_with_working_state()
    eng._compact_history()
    context = "\n".join(
        m.get("content", "") for m in eng.messages if isinstance(m, dict))
    missing = [
        name for name, needle in [
            ("changed file", "delfin/foo.py"),
            ("touched symbol", "_check_path"),
            ("test outcome", "tests/test_foo.py"),
            ("denial reason", "default branch"),
            ("operator instruction", "api_client.py"),
        ] if needle not in context
    ]
    assert not missing, (
        f"after compaction the context no longer mentions: {missing}")


def test_open_tasks_survive_compaction(tmp_path, monkeypatch):
    """Outstanding task-tool work is durable state on disk; the compacted
    context must still name the open tasks."""
    from delfin.agent import agent_tasks

    store_dir = tmp_path / "ws"
    store_dir.mkdir()
    store = agent_tasks.get_store(store_dir)
    t = store.create("Wire the working-state block into engine",
                     description="compaction part only")
    store.update(t["id"], status="in_progress")

    eng = _bare_engine(tmp_path)
    eng.repo_dir = store_dir
    eng.messages = _history_with_working_state()
    eng._compact_history()
    context = "\n".join(
        m.get("content", "") for m in eng.messages if isinstance(m, dict))
    assert "Wire the working-state block into engine" in context, (
        "the open task vanished from the compacted context")


def test_working_state_block_is_bounded(tmp_path):
    """The working-state block must have a hard character ceiling -- an
    unbounded recap can itself blow the window it is meant to save."""
    eng = _bare_engine(tmp_path)
    eng.messages = _history_with_working_state()
    # A pathological machine turn: enormous tool output.
    eng.messages.insert(3, {"role": "user", "content":
                            "[Command results]\n" + "x" * 400_000})
    eng._compact_history()
    total = sum(len(m.get("content", "")) for m in eng.messages
                if isinstance(m, dict))
    assert total < 200_000, (
        f"compacted context is {total} chars; the recap is unbounded")


def test_working_state_is_redacted(tmp_path):
    """No credential may ride into the compacted context, not even one that
    arrived inside a tool result."""
    eng = _bare_engine(tmp_path)
    eng.messages = _history_with_working_state()
    eng.messages.insert(3, {"role": "user", "content":
                            "[Command results]\n"
                            "export GITHUB_TOKEN=ghp_0123456789abcdefABCD"})
    eng._compact_history()
    context = "\n".join(
        m.get("content", "") for m in eng.messages if isinstance(m, dict))
    assert "ghp_0123456789abcdefABCD" not in context, (
        "a credential from a tool result survived compaction verbatim")


def test_deterministic_not_model_authored(tmp_path):
    """The working-state block is built from tool results, not model prose:
    identical history -> identical block. A model-authored summary would
    differ between runs (and cost a call)."""
    blocks = []
    for _ in range(2):
        eng = _bare_engine(tmp_path)
        eng.messages = _history_with_working_state()
        eng._compact_history()
        blocks.append(eng.messages[0]["content"])
    assert blocks[0] == blocks[1], (
        "the post-compaction context differs between identical runs")
