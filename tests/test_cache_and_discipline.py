"""Prompt-cache breakpoints, cache-hit visibility, and the git/orchestration
discipline addendum.

The prompt loader keeps the system prompt byte-stable for caching, but the
Anthropic API path never emitted a cache_control breakpoint and nothing
surfaced the cache-hit ratio — a prefix-stability regression would have
been invisible while doubling prefill cost.
"""

from __future__ import annotations

import textwrap
from unittest.mock import MagicMock, patch

import pytest


# ---------------------------------------------------------------------------
# Anthropic path emits a cache_control breakpoint on the system prompt
# ---------------------------------------------------------------------------

class _FakeStream:
    def __enter__(self):
        return iter(())

    def __exit__(self, *a):
        return False


def _api_client_with_capture(captured):
    from delfin.agent.api_client import APIClient
    client = APIClient.__new__(APIClient)
    client.model = "test-model"
    fake = MagicMock()

    def _stream(**kwargs):
        captured.update(kwargs)
        return _FakeStream()

    fake.messages.stream = _stream
    client.client = fake
    return client


def test_system_prompt_gets_cache_breakpoint():
    captured: dict = {}
    client = _api_client_with_capture(captured)
    list(client.stream_message("SYS PROMPT", [
        {"role": "user", "content": "hi"}], max_tokens=16))
    sys_param = captured["system"]
    assert isinstance(sys_param, list)
    assert sys_param[0]["text"] == "SYS PROMPT"
    assert sys_param[0]["cache_control"] == {"type": "ephemeral"}


def test_empty_system_prompt_passes_through():
    captured: dict = {}
    client = _api_client_with_capture(captured)
    list(client.stream_message("", [
        {"role": "user", "content": "hi"}], max_tokens=16))
    assert captured["system"] == ""


# ---------------------------------------------------------------------------
# Cache-hit ratio surfaces in the context-status block
# ---------------------------------------------------------------------------

@pytest.fixture
def agent_tree(tmp_path):
    lite_dir = tmp_path / "pack_lite"
    modes = lite_dir / "modes"
    modes.mkdir(parents=True)
    (modes / "quick.md").write_text("# quick mode")
    manifest = textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: quick
            file: modes/quick.md
            route:
              - session_manager
    """)
    (lite_dir / "manifest.yaml").write_text(manifest)
    return tmp_path


def _engine(agent_tree):
    from delfin.agent.engine import AgentEngine
    fake = MagicMock()
    fake.stream_message = MagicMock(side_effect=lambda *a, **k: iter(()))
    with patch("delfin.agent.engine.create_client", return_value=fake):
        return AgentEngine(repo_dir=agent_tree, backend="cli",
                          mode="quick", pack_dir=agent_tree)


def test_context_status_shows_cache_ratio(agent_tree):
    engine = _engine(agent_tree)
    engine.context_window_tokens = 100_000
    engine.token_usage["input"] = 10_000
    engine.token_usage["cached"] = 7_500
    block = engine._build_context_status_block()
    assert "Prompt cache" in block
    assert "75%" in block


def test_context_status_omits_cache_line_without_data(agent_tree):
    engine = _engine(agent_tree)
    engine.context_window_tokens = 100_000
    block = engine._build_context_status_block()
    assert "Prompt cache" not in block


# ---------------------------------------------------------------------------
# Git/orchestration discipline addendum is injected
# ---------------------------------------------------------------------------

def test_git_workflow_addendum_ships_in_the_pack():
    from pathlib import Path
    pack = (Path(__file__).resolve().parent.parent / "delfin" / "agent"
            / "pack" / "shared" / "git_workflow_addendum.md")
    text = pack.read_text(encoding="utf-8")
    assert "Push the branch before any merge" in text
    assert "Disjoint ownership" in text
    assert "watch_job" in text


def test_git_workflow_addendum_injected_when_present(agent_tree, tmp_path):
    from delfin.agent.prompt_loader import PromptLoader
    shared = agent_tree / "pack" / "shared"
    shared.mkdir(parents=True, exist_ok=True)
    (agent_tree / "pack" / "agents").mkdir(exist_ok=True)
    (agent_tree / "pack" / "agents" / "solo_agent.md").write_text(
        "# Solo\nYou are solo.")
    (shared / "git_workflow_addendum.md").write_text(
        "# Git discipline\nGITWF-MARKER")
    loader = PromptLoader(agent_tree)
    prompt = loader.build_system_prompt(
        role_id="solo_agent", mode_id="quick", mode_description="solo",
        route=["solo_agent"], role_index=0)
    assert "GITWF-MARKER" in prompt


def test_git_rules_do_not_contradict_each_other():
    """The addendum demanded a commit after each work unit while the role
    prompt forbade committing unprompted — a contradiction the model
    cannot follow. One rule now: commit on YOUR branch, never on the
    user's; push and merge always wait for the user."""
    from pathlib import Path
    pack = Path(__file__).resolve().parent.parent / "delfin" / "agent" / "pack"
    addendum = (pack / "shared" / "git_workflow_addendum.md").read_text(
        encoding="utf-8")
    solo = (pack / "agents" / "solo_agent.md").read_text(encoding="utf-8")

    # The blanket prohibition is gone from both places ...
    assert "never commit unprompted" not in solo
    # ... and both state the branch-scoped resolution.
    assert "on YOUR branch" in addendum
    assert "Where you commit decides whether you may" in solo
    for text in (addendum, solo):
        assert "default branch" in text
    # Push / merge stay user-gated everywhere.
    assert "Merge to main only" in addendum
    assert "wait for the user" in solo
