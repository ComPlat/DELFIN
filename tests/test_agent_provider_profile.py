"""Regression tests for provider-specific DELFIN agent profiles."""

import json
import textwrap
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest


@pytest.fixture
def agent_tree(tmp_path):
    """Create a minimal agent prompt tree for engine tests."""
    agent_dir = tmp_path / "pack"
    shared = agent_dir / "shared"
    shared.mkdir(parents=True)
    agents = agent_dir / "agents"
    agents.mkdir()

    (shared / "delfin_context.md").write_text("# Context")
    (shared / "work_cycle_rules.md").write_text("# Rules")
    (shared / "goal_decomposition_rules.md").write_text("# Goal Decomposition")
    (shared / "universal_input_template.md").write_text("")
    (shared / "minimal_final_verdict.md").write_text("")
    (agents / "session_manager.md").write_text("# Session Manager")
    (agents / "builder_agent.md").write_text("# Builder Agent")
    (agents / "test_agent.md").write_text("# Test Agent")

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
              - builder_agent
              - test_agent
    """)
    (lite_dir / "manifest.yaml").write_text(manifest)
    return tmp_path


def test_load_provider_profile_uses_codex_alias(tmp_path):
    from delfin.agent.provider_profile import load_provider_profile

    profile_path = tmp_path / "profiles.json"
    profile_path.write_text(
        '{"openai": {"thinking_budget_mult": 1.4, "success_rate": {"solo": 0.9}}}',
        encoding="utf-8",
    )

    profile = load_provider_profile("codex", profile_path)

    assert profile["thinking_budget_mult"] == 1.4
    assert profile["success_rate"]["solo"] == 0.9


def test_load_provider_profile_merges_shared_and_provider_layers(tmp_path):
    from delfin.agent.provider_profile import load_provider_profile

    profile_path = tmp_path / "profiles.json"
    profile_path.write_text(
        textwrap.dedent(
            """\
            {
              "shared": {
                "common_failures": ["grep-first before large reads"],
                "tool_usage": {"rules": ["run pytest after Python edits"]}
              },
              "openai": {
                "thinking_budget_mult": 1.4,
                "common_failures": ["verify provider wiring before edits"],
                "tool_usage": {"rules": ["keep explanations terse"]}
              }
            }
            """
        ),
        encoding="utf-8",
    )

    profile = load_provider_profile("codex", profile_path)

    assert profile["thinking_budget_mult"] == 1.4
    assert profile["common_failures"] == [
        "grep-first before large reads",
        "verify provider wiring before edits",
    ]
    assert profile["tool_usage"]["rules"] == [
        "run pytest after Python edits",
        "keep explanations terse",
    ]


def test_load_task_profile_prefers_task_performance_block(tmp_path):
    from delfin.agent.provider_profile import load_task_profile

    profile_path = tmp_path / "profiles.json"
    profile_path.write_text(
        textwrap.dedent(
            """\
            {
              "shared": {
                "task_performance": {
                  "coding": {
                    "success_rate": 0.82,
                    "thinking_budget_mult": 1.15
                  }
                }
              },
              "openai": {
                "task_performance": {
                  "coding": {
                    "success_rate": 0.91
                  }
                }
              }
            }
            """
        ),
        encoding="utf-8",
    )

    task_profile = load_task_profile("codex", "coding", profile_path)

    assert task_profile["success_rate"] == 0.91
    assert task_profile["thinking_budget_mult"] == 1.15


def test_save_provider_profile_normalizes_codex_alias(tmp_path):
    from delfin.agent.provider_profile import save_provider_profile

    profile_path = tmp_path / "profiles.json"
    save_provider_profile(
        "codex",
        {"thinking_budget_mult": 1.2, "success_rate": {"reviewed": 0.8}},
        profile_path,
    )

    saved = profile_path.read_text(encoding="utf-8")
    assert '"openai"' in saved
    assert '"codex"' not in saved


def test_save_provider_profile_keeps_shared_rules_in_shared_block(tmp_path):
    from delfin.agent.provider_profile import save_provider_profile

    profile_path = tmp_path / "profiles.json"
    profile_path.write_text(
        textwrap.dedent(
            """\
            {
              "shared": {
                "common_failures": ["grep-first before large reads"],
                "tool_usage": {"rules": ["run pytest after Python edits"]}
              }
            }
            """
        ),
        encoding="utf-8",
    )

    save_provider_profile(
        "codex",
        {
            "thinking_budget_mult": 1.2,
            "common_failures": [
                "grep-first before large reads",
                "verify provider wiring before edits",
            ],
            "tool_usage": {
                "rules": [
                    "run pytest after Python edits",
                    "keep explanations terse",
                ]
            },
        },
        profile_path,
    )

    data = json.loads(profile_path.read_text(encoding="utf-8"))
    assert data["shared"]["common_failures"] == ["grep-first before large reads"]
    assert data["openai"]["common_failures"] == ["verify provider wiring before edits"]
    assert data["openai"]["tool_usage"]["rules"] == ["keep explanations terse"]


def test_format_profile_context_includes_shared_rules_without_cycles(tmp_path):
    from delfin.agent.provider_profile import format_profile_context

    profile_path = tmp_path / "profiles.json"
    profile_path.write_text(
        textwrap.dedent(
            """\
            {
              "shared": {
                "common_failures": ["generic QM answers without docs"],
                "tool_usage": {"rules": ["search_docs before chemistry answers"]}
              },
              "openai": {
                "thinking_budget_mult": 1.1
              }
            }
            """
        ),
        encoding="utf-8",
    )

    context = format_profile_context("openai", profile_path)

    assert "Shared DELFIN failure patterns" in context
    assert "search_docs before chemistry answers" in context
    assert "Provider: openai (0 cycles)" in context


def test_format_profile_context_includes_persisted_next_steps(tmp_path):
    from delfin.agent.provider_profile import format_profile_context

    profile_path = tmp_path / "profiles.json"
    profile_path.write_text(
        textwrap.dedent(
            """\
            {
              "openai": {
                "next_steps": [
                  {
                    "task": "Track FAIL solo outcomes",
                    "status": "pending",
                    "why": "PASS-only tracking is implemented first"
                  },
                  "Expand multi-module playbook selection"
                ]
              }
            }
            """
        ),
        encoding="utf-8",
    )

    context = format_profile_context("openai", profile_path)

    assert "Persisted next steps" in context
    assert "[pending] Track FAIL solo outcomes" in context
    assert "Expand multi-module playbook selection" in context


def test_format_profile_context_includes_provider_rules(tmp_path):
    from delfin.agent.provider_profile import format_profile_context

    profile_path = tmp_path / "profiles.json"
    profile_path.write_text(
        textwrap.dedent(
            """\
            {
              "openai": {
                "communication": {
                  "rules": [
                    "I'll leave the commit to the user"
                  ]
                },
                "tool_usage": {
                  "rules": [
                    "Max 40 tool calls per task",
                    "If a command fails twice, stop retrying"
                  ]
                }
              }
            }
            """
        ),
        encoding="utf-8",
    )

    context = format_profile_context("openai", profile_path)

    assert "Provider communication rules" in context
    assert "I'll leave the commit to the user" in context
    assert "Provider tool rules" in context
    assert "Max 40 tool calls per task" in context
    assert "If a command fails twice, stop retrying" in context


def test_update_from_outcome_writes_task_performance_overlay(tmp_path):
    from delfin.agent.outcome_tracker import CycleOutcome
    from delfin.agent.provider_profile import update_from_outcome

    profile_path = tmp_path / "profiles.json"
    profile_path.write_text(
        textwrap.dedent(
            """\
            {
              "shared": {
                "task_performance": {
                  "coding": {
                    "success_rate": 0.84,
                    "thinking_budget_mult": 1.1
                  }
                }
              }
            }
            """
        ),
        encoding="utf-8",
    )

    changes = update_from_outcome(
        "openai",
        CycleOutcome(
            provider="openai",
            mode="solo",
            verdict="PASS",
            task_class="coding",
            task="Fix engine task profile routing",
        ),
        profile_path,
    )

    data = json.loads(profile_path.read_text(encoding="utf-8"))
    assert "task_success_rate" in changes
    assert data["openai"]["task_performance"]["coding"]["success_rate"] > 0.84


def test_engine_passes_provider_to_prompt_loader(agent_tree):
    from delfin.agent.engine import AgentEngine

    mock_client = MagicMock()
    with patch("delfin.agent.engine.create_client", return_value=mock_client):
        engine = AgentEngine(
            repo_dir=agent_tree,
            backend="cli",
            provider="openai",
            mode="quick",
            pack_dir=agent_tree,
        )

    assert engine.loader._active_provider == "openai"


def test_engine_passes_task_text_to_prompt_loader(agent_tree):
    from delfin.agent.engine import AgentEngine

    mock_client = MagicMock()
    with patch("delfin.agent.engine.create_client", return_value=mock_client):
        engine = AgentEngine(
            repo_dir=agent_tree,
            backend="cli",
            provider="openai",
            mode="quick",
            pack_dir=agent_tree,
        )

    with patch.object(engine.loader, "build_system_prompt", return_value="PROMPT") as build_prompt:
        prompt = engine._build_current_system_prompt(
            memory_context="mem",
            task_text="Fix delfin/build_up_complex.py",
        )

    assert prompt == "PROMPT"
    build_prompt.assert_called_once()
    assert build_prompt.call_args.kwargs["task_text"] == "Fix delfin/build_up_complex.py"
