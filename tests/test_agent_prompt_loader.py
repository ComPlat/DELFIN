"""Tests for delfin.agent.prompt_loader."""

import json
import textwrap
from pathlib import Path

import pytest


@pytest.fixture
def agent_tree(tmp_path):
    """Create a minimal agent directory tree for testing."""
    agent_dir = tmp_path / "pack"
    shared = agent_dir / "shared"
    shared.mkdir(parents=True)
    agents = agent_dir / "agents"
    agents.mkdir()
    routing = agent_dir / "routing"
    routing.mkdir()

    (shared / "delfin_context.md").write_text("# DELFIN Context\nTest context.")
    (shared / "work_cycle_rules.md").write_text("# Work Cycle Rules\nRule 1.")
    (shared / "universal_input_template.md").write_text("# Input Template")
    (shared / "minimal_final_verdict.md").write_text("# Verdict Template")
    (agents / "session_manager.md").write_text("# Session Manager\nYou are the SM.")
    (agents / "builder_agent.md").write_text("# Builder Agent\nYou are the Builder.")
    (agents / "test_agent.md").write_text("# Test Agent\nYou are the Test agent.")
    (routing / "minimal_workflow_routing.md").write_text("# Routing\nRoute rules.")

    lite_dir = tmp_path / "pack_lite"
    modes = lite_dir / "modes"
    modes.mkdir(parents=True)
    (modes / "quick.md").write_text("# Mode: quick\nDaily mode.")
    (modes / "reviewed.md").write_text("# Mode: reviewed\nReviewed mode.")

    manifest = textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        recommended_default_mode: quick
        modes:
          - id: quick
            file: modes/quick.md
            route:
              - session_manager
              - builder_agent
              - test_agent
          - id: reviewed
            file: modes/reviewed.md
            route:
              - session_manager
              - critic_agent
              - builder_agent
              - test_agent
    """)
    (lite_dir / "manifest.yaml").write_text(manifest)

    return tmp_path


def test_load_shared_context(agent_tree):
    from delfin.agent.prompt_loader import PromptLoader

    loader = PromptLoader(agent_tree)
    ctx = loader.load_shared_context()
    assert "DELFIN Context" in ctx
    assert "Work Cycle Rules" in ctx


def test_load_role_prompt(agent_tree):
    from delfin.agent.prompt_loader import PromptLoader

    loader = PromptLoader(agent_tree)
    sm = loader.load_role_prompt("session_manager")
    assert "Session Manager" in sm
    builder = loader.load_role_prompt("builder_agent")
    assert "Builder Agent" in builder


def test_load_role_prompt_missing(agent_tree):
    from delfin.agent.prompt_loader import PromptLoader

    loader = PromptLoader(agent_tree)
    result = loader.load_role_prompt("nonexistent_agent")
    assert result == ""


def test_load_routing_rules(agent_tree):
    from delfin.agent.prompt_loader import PromptLoader

    loader = PromptLoader(agent_tree)
    routing = loader.load_routing_rules()
    assert "Routing" in routing


def test_load_mode_default(agent_tree):
    from delfin.agent.prompt_loader import PromptLoader

    loader = PromptLoader(agent_tree)
    mode = loader.load_mode("quick")
    assert mode["route"] == ["session_manager", "builder_agent", "test_agent"]
    assert "quick" in mode["description"].lower()


def test_load_mode_reviewed(agent_tree):
    from delfin.agent.prompt_loader import PromptLoader

    loader = PromptLoader(agent_tree)
    mode = loader.load_mode("reviewed")
    assert "critic_agent" in mode["route"]


def test_load_mode_unknown(agent_tree):
    from delfin.agent.prompt_loader import PromptLoader

    loader = PromptLoader(agent_tree)
    with pytest.raises(ValueError, match="Unknown agent mode"):
        loader.load_mode("nonexistent")


def test_available_modes(agent_tree):
    from delfin.agent.prompt_loader import PromptLoader

    loader = PromptLoader(agent_tree)
    modes = loader.available_modes()
    assert "quick" in modes
    assert "reviewed" in modes


def test_build_system_prompt_basic(agent_tree):
    from delfin.agent.prompt_loader import PromptLoader

    loader = PromptLoader(agent_tree)
    prompt = loader.build_system_prompt(
        role_id="session_manager",
        mode_id="quick",
        mode_description="Daily mode.",
        route=["session_manager", "builder_agent", "test_agent"],
        role_index=0,
    )
    assert "Session Manager" in prompt
    assert "DELFIN Context" in prompt
    assert "Daily mode" in prompt
    assert "step 1 of 3" in prompt
    # Session manager should get routing rules
    assert "Routing" in prompt


def test_build_system_prompt_builder_no_routing(agent_tree):
    from delfin.agent.prompt_loader import PromptLoader

    loader = PromptLoader(agent_tree)
    prompt = loader.build_system_prompt(
        role_id="builder_agent",
        mode_id="quick",
        route=["session_manager", "builder_agent", "test_agent"],
        role_index=1,
    )
    assert "Builder Agent" in prompt
    # Builder should NOT get routing rules
    assert "Routing" not in prompt


def test_build_system_prompt_with_prior_outputs(agent_tree):
    from delfin.agent.prompt_loader import PromptLoader

    loader = PromptLoader(agent_tree)
    prompt = loader.build_system_prompt(
        role_id="builder_agent",
        mode_id="quick",
        route=["session_manager", "builder_agent"],
        role_index=1,
        prior_outputs={"session_manager": "Plan: fix the bug in cli.py"},
    )
    assert "Prior Role Outputs" in prompt
    assert "fix the bug in cli.py" in prompt


def test_build_system_prompt_with_memory(agent_tree):
    from delfin.agent.prompt_loader import PromptLoader

    loader = PromptLoader(agent_tree)
    prompt = loader.build_system_prompt(
        role_id="session_manager",
        mode_id="quick",
        memory_context="Last session: refactored config.py",
    )
    assert "Project Memory" in prompt
    assert "refactored config.py" in prompt


def test_build_system_prompt_injects_relevant_profile_playbook(agent_tree):
    from delfin.agent.prompt_loader import PromptLoader

    profile_path = agent_tree / "profiles.json"
    profile_path.write_text(
        json.dumps(
            {
                "shared": {
                    "playbooks": {
                        "build_up_complex": {
                            "description": "Editing the metal complex builder",
                            "steps": [
                                "1. Grep for the target function name in build_up_complex.py",
                            ],
                            "key_invariants": [
                                "PSO places ligands as rigid bodies",
                            ],
                        }
                    },
                    "codebase_map": {
                        "modules": {
                            "build_up_complex.py": {"lines": 1531},
                        },
                        "test_mapping": {
                            "build_up_complex.py": ["test_build_up_complex.py"],
                        },
                    },
                }
            }
        ),
        encoding="utf-8",
    )

    loader = PromptLoader(agent_tree)
    loader._active_provider = "openai"
    loader._profile_path = profile_path
    prompt = loader.build_system_prompt(
        role_id="builder_agent",
        mode_id="quick",
        task_text="Fix a regression in delfin/build_up_complex.py",
    )

    assert "Relevant Playbook" in prompt
    assert "Target module: build_up_complex.py" in prompt
    assert "Editing the metal complex builder" in prompt
    assert "PSO places ligands as rigid bodies" in prompt
    assert "test_build_up_complex.py" in prompt


def test_build_system_prompt_skips_profile_playbook_for_unknown_module(agent_tree):
    from delfin.agent.prompt_loader import PromptLoader

    profile_path = agent_tree / "profiles.json"
    profile_path.write_text(
        json.dumps(
            {
                "shared": {
                    "playbooks": {
                        "build_up_complex": {
                            "description": "Editing the metal complex builder",
                            "steps": ["1. Grep first"],
                        }
                    },
                    "codebase_map": {
                        "modules": {
                            "build_up_complex.py": {"lines": 1531},
                        }
                    },
                }
            }
        ),
        encoding="utf-8",
    )

    loader = PromptLoader(agent_tree)
    loader._active_provider = "openai"
    loader._profile_path = profile_path
    prompt = loader.build_system_prompt(
        role_id="builder_agent",
        mode_id="quick",
        task_text="Fix a regression in delfin/unknown_module.py",
    )

    assert "Target module: build_up_complex.py" not in prompt


def test_build_system_prompt_injects_repo_map(agent_tree):
    from delfin.agent.prompt_loader import PromptLoader

    src_dir = agent_tree / "delfin"
    src_dir.mkdir()
    (src_dir / "build_up_complex.py").write_text(
        "def _pso_fitness():\n    return 1\n\nclass SwarmBuilder:\n    pass\n",
        encoding="utf-8",
    )
    tests_dir = agent_tree / "tests"
    tests_dir.mkdir()
    (tests_dir / "test_build_up_complex.py").write_text(
        "def test_placeholder():\n    assert True\n",
        encoding="utf-8",
    )

    profile_path = agent_tree / "profiles.json"
    profile_path.write_text(
        json.dumps(
            {
                "shared": {
                    "playbooks": {
                        "build_up_complex": {
                            "description": "Editing the metal complex builder",
                            "steps": ["1. Grep first"],
                        }
                    },
                    "codebase_map": {
                        "modules": {
                            "build_up_complex.py": {},
                        },
                        "test_mapping": {
                            "build_up_complex.py": ["test_build_up_complex.py"],
                        },
                    },
                }
            }
        ),
        encoding="utf-8",
    )

    loader = PromptLoader(agent_tree)
    loader._active_provider = "openai"
    loader._profile_path = profile_path
    prompt = loader.build_system_prompt(
        role_id="builder_agent",
        mode_id="quick",
        task_text="Fix regression in delfin/build_up_complex.py around _pso_fitness",
    )

    assert "Repo Map" in prompt
    assert "delfin/build_up_complex.py" in prompt
    assert "_pso_fitness" in prompt
    assert "tests/test_build_up_complex.py" in prompt


def test_caching(agent_tree):
    from delfin.agent.prompt_loader import PromptLoader

    loader = PromptLoader(agent_tree)
    # Load twice, should use cache
    ctx1 = loader.load_shared_context()
    ctx2 = loader.load_shared_context()
    assert ctx1 == ctx2
    # Verify cache is populated
    assert len(loader._cache) > 0


def test_build_system_prompt_skips_duplicate_profile_in_same_session(agent_tree):
    from delfin.agent.prompt_loader import PromptLoader

    profile_path = agent_tree / "profiles.json"
    profile_path.write_text(
        json.dumps(
            {
                "shared": {
                    "common_failures": ["grep-first before large reads"],
                    "tool_usage": {"rules": ["run pytest after Python edits"]},
                },
                "openai": {
                    "success_rate": {"solo": 0.91},
                },
            }
        ),
        encoding="utf-8",
    )

    loader = PromptLoader(agent_tree)
    loader._active_provider = "openai"
    loader._profile_path = profile_path

    prompt1 = loader.build_system_prompt(
        role_id="builder_agent",
        mode_id="quick",
        session_key="s1",
    )
    prompt2 = loader.build_system_prompt(
        role_id="builder_agent",
        mode_id="quick",
        session_key="s1",
    )

    assert "Provider Profile" in prompt1
    assert "Provider Profile" not in prompt2
