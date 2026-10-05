"""The engine's one-shot fresh start keeps only the current, fused prompt.

When the terminal finds the current context over budget it fuses a small
fresh block into THIS turn's user message and sets engine.start_fresh. The
engine then archives the earlier history, sends only that one message, and
clears the flag, so the next turn is judged on its own small context and is
not a fresh restart again.
"""

from __future__ import annotations

import textwrap
from types import SimpleNamespace
from unittest.mock import patch

from delfin.agent.engine import AgentEngine


def _pack_tree(tmp_path):
    """Minimal prompt-pack tree so AgentEngine boots (as in
    tests/test_engine_compaction.py:agent_tree)."""
    agent_dir = tmp_path / "pack"
    shared = agent_dir / "shared"
    shared.mkdir(parents=True)
    agents = agent_dir / "agents"
    agents.mkdir()
    (shared / "delfin_context.md").write_text("# DELFIN Context\nTest.")
    (shared / "work_cycle_rules.md").write_text("# Work Cycle Rules\nRule 1.")
    (shared / "universal_input_template.md").write_text("# Input Template")
    (shared / "minimal_final_verdict.md").write_text("# Verdict")
    (agents / "solo_agent.md").write_text("# Solo Agent\nYou are solo.")
    lite = tmp_path / "pack_lite"
    modes = lite / "modes"
    modes.mkdir(parents=True)
    (modes / "solo.md").write_text("# Mode: quick\nQuick mode.")
    (lite / "manifest.yaml").write_text(textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        recommended_default_mode: solo
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - solo_agent
    """))
    return tmp_path


class _CaptureClient:
    """Scripted streaming client that records whatever message list the
    engine sends to the transport."""

    def __init__(self):
        self.last_messages = None

    def stream_message(self, *, messages, system="", max_tokens=0, **kw):
        self.last_messages = list(messages or [])
        yield SimpleNamespace(type="text_delta", text="ok")
        yield SimpleNamespace(type="message_delta", input_tokens=30,
                              output_tokens=3, cost_usd=0.0,
                              stop_reason="end_turn", text="")


def _make_engine(tmp_path):
    client = _CaptureClient()
    with patch("delfin.agent.engine.create_client", return_value=client):
        engine = AgentEngine(
            repo_dir=tmp_path, backend="cli", mode="quick", pack_dir=tmp_path,
        )
    engine.client = client
    return engine, client


def test_a_fresh_turn_sends_only_the_fused_prompt_and_archives_the_rest(
        tmp_path, monkeypatch):
    tmp_path = _pack_tree(tmp_path)
    engine, client = _make_engine(tmp_path)
    archived = {}
    import delfin.agent.session_store as store
    monkeypatch.setattr(store, "archive_pre_compaction_transcript",
                        lambda sid, msgs, info=None: archived.update(
                            n=len(msgs), info=info))
    engine.messages = [{"role": "user" if i % 2 == 0 else "assistant",
                        "content": f"turn {i} " + "x" * 40} for i in range(30)]
    engine.start_fresh = True
    fused = ("[Fresh context - full history cut at the token budget]\n"
             "task: build task_state.py\n\ncontinue with phase 3")
    engine.stream_response(user_message=fused)
    sent = list(client.last_messages or [])
    assert len(sent) == 1, len(sent)
    assert "task_state.py" in sent[0]["content"]
    assert "continue with phase 3" in sent[0]["content"]
    assert engine.start_fresh is False
    assert archived.get("n") == 30
    assert (archived.get("info") or {}).get("kind") == "fresh_restart"


def test_the_turn_after_a_fresh_start_is_not_fresh_again(tmp_path, monkeypatch):
    tmp_path = _pack_tree(tmp_path)
    engine, client = _make_engine(tmp_path)
    import delfin.agent.session_store as store
    monkeypatch.setattr(store, "archive_pre_compaction_transcript",
                        lambda *a, **k: None)
    engine.messages = [{"role": "user", "content": "old " * 50},
                       {"role": "assistant", "content": "ok"}]
    engine.start_fresh = True
    engine.stream_response(user_message="fused prompt")
    engine.stream_response(user_message="next prompt")
    sent = list(client.last_messages or [])
    assert len(sent) == 3, [m["role"] for m in sent]
    assert sent[-1]["content"] == "next prompt"


def test_without_the_flag_history_is_sent_whole(tmp_path):
    tmp_path = _pack_tree(tmp_path)
    engine, client = _make_engine(tmp_path)
    engine.messages = [{"role": "user", "content": "a"},
                       {"role": "assistant", "content": "b"}]
    engine.stream_response(user_message="c")
    assert len(client.last_messages) == 3
