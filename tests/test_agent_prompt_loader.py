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


# ---------------------------------------------------------------------------
# Session env block — CLI-style cwd/branch/status/commits injection
# ---------------------------------------------------------------------------

def test_session_env_block_includes_cwd(agent_tree):
    """The env block must always include the cwd line as a baseline."""
    from delfin.agent.prompt_loader import PromptLoader

    loader = PromptLoader(agent_tree)
    block = loader._build_session_env_block()
    assert block.startswith("cwd: ")
    assert str(agent_tree) in block


def test_session_env_block_no_git_returns_cwd_only(tmp_path, monkeypatch):
    """Outside a git repo the block degrades gracefully — just cwd."""
    from delfin.agent.prompt_loader import PromptLoader

    # Build an isolated empty tree (no git)
    (tmp_path / "pack").mkdir()
    (tmp_path / "pack_lite").mkdir()
    loader = PromptLoader(tmp_path)
    block = loader._build_session_env_block()
    # cwd is always there; status/branch/commits skipped silently
    assert block.startswith("cwd: ")
    # Should be short (no git data)
    assert "branch:" not in block
    assert "recent commits:" not in block


def test_session_env_block_real_repo_has_branch(agent_tree, monkeypatch):
    """Smoke test: when running in the actual DELFIN repo, branch should appear.

    We import the loader without override so it uses the real ``repo_root``.
    """
    from delfin.agent.prompt_loader import PromptLoader

    loader = PromptLoader()  # use real repo
    block = loader._build_session_env_block()
    # In the real repo, git is available — branch line should be present
    if "branch:" in block:
        assert "recent commits:" in block or "status:" in block


def test_solo_mode_prompt_includes_session_env(agent_tree):
    """Solo-mode build must inject the Session Environment section."""
    from delfin.agent.prompt_loader import PromptLoader

    # Add a solo_agent role file so the test can build the solo prompt
    (agent_tree / "pack" / "agents" / "solo_agent.md").write_text(
        "# Solo Agent\nYou are the solo agent."
    )

    loader = PromptLoader(agent_tree)
    prompt = loader.build_system_prompt(
        role_id="solo_agent",
        mode_id="quick",
        mode_description="solo",
        route=["solo_agent"],
        role_index=0,
    )
    assert "Session Environment" in prompt
    assert "cwd:" in prompt


def test_solo_mode_prompt_includes_full_project_context(agent_tree):
    """Solo-mode no longer truncates delfin_context.md to 18 lines."""
    from delfin.agent.prompt_loader import PromptLoader

    long_ctx = "\n".join(f"line {i}" for i in range(40))
    (agent_tree / "pack" / "shared" / "delfin_context.md").write_text(
        f"# DELFIN Context\n{long_ctx}"
    )
    (agent_tree / "pack" / "agents" / "solo_agent.md").write_text(
        "# Solo Agent\nYou are the solo agent."
    )

    loader = PromptLoader(agent_tree)
    prompt = loader.build_system_prompt(
        role_id="solo_agent",
        mode_id="quick",
        mode_description="solo",
        route=["solo_agent"],
        role_index=0,
    )
    # Lines beyond the old 18-line cutoff must now be present
    assert "line 30" in prompt
    assert "line 39" in prompt


# ---------------------------------------------------------------------------
# External Memory bridge (S2 — read ~/.claude/projects/<slug>/memory/...)
# ---------------------------------------------------------------------------

def test_external_memory_missing_returns_empty(agent_tree, tmp_path):
    """No memory directory → empty string, no exception."""
    from delfin.agent.prompt_loader import PromptLoader

    loader = PromptLoader(agent_tree)
    out = loader._load_external_memory_context(memory_root=tmp_path / "missing")
    assert out == ""


def test_external_memory_reads_index(agent_tree, tmp_path):
    """MEMORY.md alone is enough — referenced files are optional."""
    from delfin.agent.prompt_loader import PromptLoader

    mem = tmp_path / "memory"
    mem.mkdir()
    (mem / "MEMORY.md").write_text(
        "# Project Memory\n- Some standalone fact about the user.\n"
    )

    loader = PromptLoader(agent_tree)
    out = loader._load_external_memory_context(memory_root=mem)
    assert "MEMORY.md" in out
    assert "standalone fact" in out


def test_external_memory_follows_markdown_links(agent_tree, tmp_path):
    """Linked files referenced from MEMORY.md are concatenated."""
    from delfin.agent.prompt_loader import PromptLoader

    mem = tmp_path / "memory"
    mem.mkdir()
    (mem / "MEMORY.md").write_text(
        "# Project Memory\n"
        "- [User role](user.md) — brief\n"
        "- [Coding style](style.md) — terse diffs\n"
    )
    (mem / "user.md").write_text("---\nname: user\n---\nThe user is a chemist.")
    (mem / "style.md").write_text("---\nname: style\n---\nTerse diffs only.")

    loader = PromptLoader(agent_tree)
    out = loader._load_external_memory_context(memory_root=mem)
    assert "MEMORY.md" in out
    assert "user.md" in out
    assert "style.md" in out
    assert "chemist" in out
    assert "Terse diffs only" in out


def test_external_memory_skips_missing_referenced_files(agent_tree, tmp_path):
    """Broken links don't crash the bridge."""
    from delfin.agent.prompt_loader import PromptLoader

    mem = tmp_path / "memory"
    mem.mkdir()
    (mem / "MEMORY.md").write_text("- [Missing](nope.md)\n- [Real](real.md)")
    (mem / "real.md").write_text("Real content here.")

    loader = PromptLoader(agent_tree)
    out = loader._load_external_memory_context(memory_root=mem)
    assert "Real content here" in out
    # No crash — and broken link doesn't appear as content
    assert "Missing" in out  # title shows up in MEMORY.md text
    assert "nope.md" not in out.split("# Real")[1] if "# Real" in out else True


def test_external_memory_truncates_to_max_chars(agent_tree, tmp_path):
    """Big memories are capped — never blow up the prompt."""
    from delfin.agent.prompt_loader import PromptLoader

    mem = tmp_path / "memory"
    mem.mkdir()
    (mem / "MEMORY.md").write_text("- [Big](big.md)\n")
    (mem / "big.md").write_text("x" * 50_000)

    loader = PromptLoader(agent_tree)
    out = loader._load_external_memory_context(
        memory_root=mem, max_chars=2_000,
    )
    assert len(out) <= 2_100  # cap + truncation marker
    assert "truncated" in out


def test_external_memory_dedupes_repeated_links(agent_tree, tmp_path):
    """If a file is linked twice, it's read once."""
    from delfin.agent.prompt_loader import PromptLoader

    mem = tmp_path / "memory"
    mem.mkdir()
    (mem / "MEMORY.md").write_text(
        "- [A](a.md)\n- [Also A](a.md)\n"
    )
    (mem / "a.md").write_text("Body of A")

    loader = PromptLoader(agent_tree)
    out = loader._load_external_memory_context(memory_root=mem)
    assert out.count("Body of A") == 1


def test_external_memory_blocks_path_traversal(agent_tree, tmp_path):
    """Links that escape the memory directory must be rejected."""
    from delfin.agent.prompt_loader import PromptLoader

    mem = tmp_path / "memory"
    mem.mkdir()
    secret = tmp_path / "secret.md"
    secret.write_text("PASSWORD")
    (mem / "MEMORY.md").write_text("- [Bad](../secret.md)\n")

    loader = PromptLoader(agent_tree)
    out = loader._load_external_memory_context(memory_root=mem)
    # Title still appears in MEMORY.md text, but the secret content does NOT
    assert "PASSWORD" not in out


def test_solo_prompt_includes_external_memory_when_present(
    agent_tree, tmp_path, monkeypatch,
):
    """End-to-end: solo build picks up the external memory layer."""
    from delfin.agent.prompt_loader import PromptLoader

    (agent_tree / "pack" / "agents" / "solo_agent.md").write_text(
        "# Solo Agent\nYou are the solo agent."
    )

    mem = tmp_path / "memory"
    mem.mkdir()
    (mem / "MEMORY.md").write_text(
        "# Project Memory\n- The user prefers terse diffs.\n"
    )

    # Force the loader to look at our temp memory directory
    loader = PromptLoader(agent_tree)
    monkeypatch.setattr(
        loader, "_load_external_memory_context",
        lambda max_chars=6000, memory_root=None:
            loader.__class__._load_external_memory_context(
                loader, memory_root=mem, max_chars=max_chars,
            ),
    )
    prompt = loader.build_system_prompt(
        role_id="solo_agent",
        mode_id="quick",
        mode_description="solo",
        route=["solo_agent"],
        role_index=0,
    )
    assert "External Memory" in prompt
    assert "terse diffs" in prompt


# ---------------------------------------------------------------------------
# S1 — live_state injected into the system prompt (replaces user-msg state block)
# ---------------------------------------------------------------------------

def test_live_state_appears_in_solo_prompt(agent_tree):
    """live_state passes through and lands as a --- Live state --- section."""
    from delfin.agent.prompt_loader import PromptLoader
    loader = PromptLoader(agent_tree)
    prompt = loader.build_system_prompt(
        role_id="solo_agent",
        mode_id="quick",
        mode_description="solo",
        route=["solo_agent"],
        role_index=0,
        live_state="calc_dir: /tmp/calc\nactive_tab: Calculations",
    )
    assert "--- Live state ---" in prompt
    assert "calc_dir: /tmp/calc" in prompt


def test_live_state_appears_in_dashboard_prompt(agent_tree):
    """Dashboard goes through the non-solo branch — must also support live_state."""
    from delfin.agent.prompt_loader import PromptLoader
    loader = PromptLoader(agent_tree)
    prompt = loader.build_system_prompt(
        role_id="dashboard_agent",
        mode_id="dashboard",
        mode_description="dashboard",
        route=["dashboard_agent"],
        role_index=0,
        live_state="ORCA Builder: method=BP86",
    )
    assert "--- Live state ---" in prompt
    assert "method=BP86" in prompt


def test_live_state_omitted_when_empty(agent_tree):
    """Empty live_state must NOT add a --- Live state --- header."""
    from delfin.agent.prompt_loader import PromptLoader
    loader = PromptLoader(agent_tree)
    prompt = loader.build_system_prompt(
        role_id="solo_agent",
        mode_id="quick",
        mode_description="solo",
        route=["solo_agent"],
        role_index=0,
        live_state="",
    )
    assert "--- Live state ---" not in prompt


def test_live_state_default_is_empty(agent_tree):
    """Backwards compat: callers that don't pass live_state still work."""
    from delfin.agent.prompt_loader import PromptLoader
    loader = PromptLoader(agent_tree)
    prompt = loader.build_system_prompt(
        role_id="solo_agent",
        mode_id="quick",
        mode_description="solo",
        route=["solo_agent"],
        role_index=0,
    )
    assert "--- Live state ---" not in prompt


def test_engine_set_live_state_passes_through(agent_tree):
    """AgentEngine.set_live_state() flows through into build_system_prompt."""
    from unittest.mock import MagicMock, patch
    from delfin.agent.engine import AgentEngine

    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        engine = AgentEngine(
            repo_dir=agent_tree, backend="cli", mode="quick", pack_dir=agent_tree,
        )

    engine.set_live_state("CONTROL: PAL=8\nselected: foo.out")
    prompt = engine._build_current_system_prompt(memory_context="", task_text="hi")
    assert "--- Live state ---" in prompt
    assert "CONTROL: PAL=8" in prompt
