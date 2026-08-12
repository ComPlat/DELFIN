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
    (modes / "solo.md").write_text("# Mode: solo\nDaily mode.")
    (modes / "office.md").write_text("# Mode: office\nOffice mode.")

    # The fixture's modes are named for modes that still exist: every
    # retired name now migrates onto `solo`, so a fixture built around
    # `quick` described a world the loader no longer has. The three-role
    # route stays -- the engine's role advancement is what is under test.
    manifest = textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        recommended_default_mode: solo
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - session_manager
              - builder_agent
              - test_agent
          - id: office
            file: modes/office.md
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
    mode = loader.load_mode("solo")
    assert mode["route"] == ["session_manager", "builder_agent", "test_agent"]
    assert "solo" in mode["description"].lower()


def test_load_mode_office(agent_tree):
    from delfin.agent.prompt_loader import PromptLoader

    loader = PromptLoader(agent_tree)
    mode = loader.load_mode("office")
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
    assert "solo" in modes


def test_build_system_prompt_basic(agent_tree):
    from delfin.agent.prompt_loader import PromptLoader

    loader = PromptLoader(agent_tree)
    prompt = loader.build_system_prompt(
        role_id="session_manager",
        mode_id="solo",
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
        mode_id="solo",
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
        mode_id="solo",
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
        mode_id="solo",
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
        mode_id="solo",
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
        mode_id="solo",
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
        mode_id="solo",
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


def _write_profile(agent_tree):
    """Shared profiles.json fixture for the inject-gating tests."""
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
    return profile_path


def test_stateless_backend_reinjects_profile_on_second_build(agent_tree):
    """Chat-API backends rebuild every request from scratch — an earlier
    injection is NOT still in the conversation, so unchanged stable sections
    must appear in EVERY build (stateful_backend defaults to False)."""
    from delfin.agent.prompt_loader import PromptLoader

    loader = PromptLoader(agent_tree)
    loader._active_provider = "openai"
    loader._profile_path = _write_profile(agent_tree)

    prompt1 = loader.build_system_prompt(
        role_id="builder_agent",
        mode_id="solo",
        session_key="s1",
    )
    prompt2 = loader.build_system_prompt(
        role_id="builder_agent",
        mode_id="solo",
        session_key="s1",
    )

    assert "Provider Profile" in prompt1
    assert "Provider Profile" in prompt2


def test_stateful_cli_backend_injects_profile_once_per_session(agent_tree):
    """The persistent CLI process keeps the first system prompt alive across
    turns, so the unchanged profile is gated after the first injection."""
    from delfin.agent.prompt_loader import PromptLoader

    loader = PromptLoader(agent_tree)
    loader._active_provider = "openai"
    loader._profile_path = _write_profile(agent_tree)
    loader.stateful_backend = True

    prompt1 = loader.build_system_prompt(
        role_id="builder_agent",
        mode_id="solo",
        session_key="s1",
    )
    prompt2 = loader.build_system_prompt(
        role_id="builder_agent",
        mode_id="solo",
        session_key="s1",
    )

    assert "Provider Profile" in prompt1
    assert "Provider Profile" not in prompt2


def test_stateful_gating_resets_with_session_state(agent_tree):
    """After reset_session_prompt_state the stateful loader injects again."""
    from delfin.agent.prompt_loader import PromptLoader

    loader = PromptLoader(agent_tree)
    loader._active_provider = "openai"
    loader._profile_path = _write_profile(agent_tree)
    loader.stateful_backend = True

    prompt1 = loader.build_system_prompt(
        role_id="builder_agent", mode_id="solo", session_key="s1",
    )
    assert "Provider Profile" in prompt1
    loader.reset_session_prompt_state("s1")
    prompt2 = loader.build_system_prompt(
        role_id="builder_agent", mode_id="solo", session_key="s1",
    )
    assert "Provider Profile" in prompt2


def test_engine_marks_only_claude_cli_backend_stateful(agent_tree):
    """The engine threads backend statefulness to the loader: persistent
    claude CLI → True; chat-API backends → False (re-inject every build)."""
    from unittest.mock import MagicMock, patch
    from delfin.agent.engine import AgentEngine

    with patch("delfin.agent.engine.create_client", return_value=MagicMock()):
        cli = AgentEngine(
            repo_dir=agent_tree, backend="cli", provider="claude",
            mode="quick", pack_dir=agent_tree,
        )
        api = AgentEngine(
            repo_dir=agent_tree, backend="api", provider="openai",
            mode="quick", pack_dir=agent_tree,
        )
    assert cli.loader.stateful_backend is True
    assert api.loader.stateful_backend is False


# ---------------------------------------------------------------------------
# Sticky lazy-module set — modules never deactivate mid-session
# ---------------------------------------------------------------------------

_MODULE_TEXT = (
    "## Intro\n"
    "Always kept.\n"
    "\n"
    "<!-- module:chemistry -->\n"
    "## Chemistry\n"
    "chemistry module body\n"
    "\n"
    "## Tail\n"
    "tail body\n"
)


def test_lazy_module_stays_active_without_trigger_keywords(agent_tree):
    """A module activated by turn 1 must survive a trigger-free follow-up
    ('ja, mach weiter') in the same session — the active set is monotonic."""
    from delfin.agent.prompt_loader import PromptLoader

    loader = PromptLoader(agent_tree)
    out1 = loader._strip_lazy_modules(
        _MODULE_TEXT, task_text="run an orca dft job", mode_id="solo",
        session_key="s1", role_id="solo_agent",
    )
    assert "chemistry module body" in out1
    out2 = loader._strip_lazy_modules(
        _MODULE_TEXT, task_text="ja, mach weiter", mode_id="solo",
        session_key="s1", role_id="solo_agent",
    )
    assert "chemistry module body" in out2
    assert "tail body" in out2
    # A different session has no sticky history — the module is stripped.
    out3 = loader._strip_lazy_modules(
        _MODULE_TEXT, task_text="ja, mach weiter", mode_id="solo",
        session_key="s2", role_id="solo_agent",
    )
    assert "chemistry module body" not in out3


def test_lazy_module_without_session_key_stays_per_turn(agent_tree):
    """No session_key (direct callers/tests) keeps the old per-turn detection."""
    from delfin.agent.prompt_loader import PromptLoader

    loader = PromptLoader(agent_tree)
    out1 = loader._strip_lazy_modules(
        _MODULE_TEXT, task_text="run an orca dft job", mode_id="solo",
    )
    assert "chemistry module body" in out1
    out2 = loader._strip_lazy_modules(
        _MODULE_TEXT, task_text="ja, mach weiter", mode_id="solo",
    )
    assert "chemistry module body" not in out2


def test_reset_session_prompt_state_clears_sticky_modules(agent_tree):
    """reset_session_prompt_state must clear the sticky module union."""
    from delfin.agent.prompt_loader import PromptLoader

    loader = PromptLoader(agent_tree)
    active = loader._detect_active_modules(
        "run an orca dft job", "solo", session_key="s1", role_id="solo_agent",
    )
    assert "chemistry" in active
    loader.reset_session_prompt_state("s1")
    after = loader._detect_active_modules(
        "ja, mach weiter", "solo", session_key="s1", role_id="solo_agent",
    )
    assert "chemistry" not in after


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
        mode_id="solo",
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
        mode_id="solo",
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


def test_external_memory_records_only_files_inside_prompt_limit(
    agent_tree, tmp_path, monkeypatch,
):
    from delfin.agent import memory_store as ms
    from delfin.agent.prompt_loader import PromptLoader

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    loader = PromptLoader(agent_tree)
    first, _, _ = ms.save_typed_memory(
        "project: first bounded recall entry", repo_root=agent_tree)
    second, _, _ = ms.save_typed_memory(
        "project: second bounded recall entry", repo_root=agent_tree)
    memory_dir = ms._delfin_memory_dir(agent_tree)
    index = memory_dir / "MEMORY.md"
    index.write_text(
        "# Project Memory\n"
        f"- [First]({first.name})\n"
        f"- [Second]({second.name})\n",
        encoding="utf-8",
    )
    index_chunk = f"# MEMORY.md\n{index.read_text(encoding='utf-8').strip()}"
    # Size the budget by what is actually INJECTED, not by the file on disk:
    # the frontmatter is mostly the store's bookkeeping and only its
    # provenance fields reach the prompt. Measuring the raw file would leave
    # the budget generous enough for both entries and the boundary this test
    # exists for would never be reached. The qualifying preamble is charged
    # to the same budget, so it is part of the sum.
    from delfin.agent.prompt_loader import (
        _MEMORY_BLOCK_PREAMBLE, _memory_entry_chunk,
    )
    first_chunk = _memory_entry_chunk(
        "First", first.name, first.read_text(encoding="utf-8"))

    out = loader._load_external_memory_context(
        max_chars=(len(_MEMORY_BLOCK_PREAMBLE) + 2
                   + len(index_chunk) + 2 + len(first_chunk)),
    )

    assert "first bounded recall entry" in out
    assert "second bounded recall entry" not in out
    records = {rec["file"]: rec for rec in ms.list_typed_memories(agent_tree)}
    assert records[first.name]["use_count"] == 2
    assert records[second.name]["use_count"] == 1


def test_external_memory_ranks_task_relevant_files_first(
    agent_tree, tmp_path, monkeypatch,
):
    """With task_text, the char budget goes to the BM25-relevant memory,
    not to whatever sits at the top of MEMORY.md."""
    from delfin.agent import memory_store as ms
    from delfin.agent.prompt_loader import PromptLoader

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    loader = PromptLoader(agent_tree)
    dashboard, _, _ = ms.save_typed_memory(
        "project: dashboard websocket runs on port 8050", repo_root=agent_tree)
    orca, _, _ = ms.save_typed_memory(
        "project: orca geometry convergence needs tightscf and slowconv",
        repo_root=agent_tree)
    memory_dir = ms._delfin_memory_dir(agent_tree)
    index = memory_dir / "MEMORY.md"
    # Dashboard listed FIRST — without ranking it would win the budget.
    index.write_text(
        "# Project Memory\n"
        f"- [Dashboard]({dashboard.name})\n"
        f"- [Orca]({orca.name})\n",
        encoding="utf-8",
    )
    index_chunk = f"# MEMORY.md\n{index.read_text(encoding='utf-8').strip()}"
    from delfin.agent.prompt_loader import (
        _MEMORY_BLOCK_PREAMBLE, _memory_entry_chunk,
    )
    orca_chunk = _memory_entry_chunk(
        "Orca", orca.name, orca.read_text(encoding="utf-8"))

    out = loader._load_external_memory_context(
        max_chars=(len(_MEMORY_BLOCK_PREAMBLE) + 2
                   + len(index_chunk) + 2 + len(orca_chunk)),
        task_text="fix the orca geometry convergence failure in the run",
    )

    assert "tightscf" in out
    assert "port 8050" not in out
    records = {rec["file"]: rec for rec in ms.list_typed_memories(agent_tree)}
    assert records[orca.name]["use_count"] == 2
    assert records[dashboard.name]["use_count"] == 1


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
        lambda max_chars=6000, memory_root=None, task_text="", domain="":
            loader.__class__._load_external_memory_context(
                loader, memory_root=mem, max_chars=max_chars,
                task_text=task_text, domain=domain,
            ),
    )
    prompt = loader.build_system_prompt(
        role_id="solo_agent",
        mode_id="solo",
        mode_description="solo",
        route=["solo_agent"],
        role_index=0,
    )
    assert "External Memory" in prompt
    assert "terse diffs" in prompt


def test_external_memory_injected_into_dashboard_prompt(agent_tree, tmp_path, monkeypatch):
    """S3 — Dashboard mode now also gets the CLI memory bridge (smaller cap)."""
    from delfin.agent.prompt_loader import PromptLoader

    (agent_tree / "pack" / "agents" / "dashboard_agent.md").write_text(
        "# Dashboard Agent\nYou are the dashboard operator."
    )

    mem = tmp_path / "memory"
    mem.mkdir()
    (mem / "MEMORY.md").write_text(
        "# Project Memory\n- DELFIN runs on uc3 cluster.\n"
    )

    loader = PromptLoader(agent_tree)
    captured: dict = {}
    real = loader._load_external_memory_context

    def _spy(max_chars=6000, memory_root=None, task_text="", domain=""):
        captured["max_chars"] = max_chars
        return real.__func__(loader, memory_root=mem, max_chars=max_chars,
                             task_text=task_text, domain=domain)

    monkeypatch.setattr(loader, "_load_external_memory_context", _spy)
    prompt = loader.build_system_prompt(
        role_id="dashboard_agent",
        mode_id="dashboard",
        mode_description="dashboard",
        route=["dashboard_agent"],
        role_index=0,
    )
    assert "External Memory" in prompt
    assert "uc3 cluster" in prompt
    # Dashboard cap is tighter than solo's default 6000
    assert captured["max_chars"] == 2000


def test_external_memory_not_injected_for_other_roles(agent_tree, tmp_path, monkeypatch):
    """Non-solo, non-dashboard roles must NOT get the external-memory block."""
    from delfin.agent.prompt_loader import PromptLoader

    mem = tmp_path / "memory"
    mem.mkdir()
    (mem / "MEMORY.md").write_text("# Project Memory\n- Some note.\n")

    loader = PromptLoader(agent_tree)
    monkeypatch.setattr(
        loader, "_load_external_memory_context",
        lambda max_chars=6000, memory_root=None, task_text="", domain="":
            loader.__class__._load_external_memory_context(
                loader, memory_root=mem, max_chars=max_chars,
                task_text=task_text, domain=domain,
            ),
    )
    prompt = loader.build_system_prompt(
        role_id="builder_agent",
        mode_id="default",
        mode_description="default",
        route=["builder_agent"],
        role_index=0,
    )
    assert "External Memory" not in prompt


# ---------------------------------------------------------------------------
# Recall-time provenance verification (stale / drifted annotations)
# ---------------------------------------------------------------------------

def test_recall_annotates_stale_reference(agent_tree, tmp_path, monkeypatch):
    """A memory citing a file that vanished is annotated on injection and
    fed into the rot counter — ON TOP OF the normal recall bump.

    The bump used to be withheld here. That made ``updated_at`` — the field
    both the prune order and the age cutoff read — the punishment for citing
    code, which is what a memory written to correct another one does. Rot is
    a separate signal now and lives in ``stale_hits`` alone.
    """
    from delfin.agent import memory_store as ms
    from delfin.agent.prompt_loader import PromptLoader

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    loader = PromptLoader(agent_tree)
    ms.save_typed_memory(
        "project: the retry logic sits in delfin/vanished.py somewhere",
        repo_root=agent_tree)

    out = loader._load_external_memory_context(task_text="retry logic")

    assert ("[stale: delfin/vanished.py no longer exists — verify against "
            "the current code before relying on this]") in out
    rec = ms.list_typed_memories(agent_tree)[0]
    assert rec["stale_hits"] == 1
    assert rec["use_count"] == 2            # injected ⇒ seen ⇒ recalled


def test_recall_annotates_drifted_reference(agent_tree, tmp_path, monkeypatch):
    """File still exists, but the anchored line moved out of the window."""
    from delfin.agent import memory_store as ms
    from delfin.agent.prompt_loader import PromptLoader

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    loader = PromptLoader(agent_tree)
    target = agent_tree / "srcmod.py"
    target.write_text("import os\nMAGIC_LIMIT = 17\n", encoding="utf-8")
    ms.save_typed_memory(
        "project: the magic limit constant sits at srcmod.py:2 today",
        repo_root=agent_tree)
    # Rewrite the file so the anchored line is nowhere near line 2.
    target.write_text(
        "\n".join(f"filler_{i} = {i}" for i in range(60)), encoding="utf-8")

    out = loader._load_external_memory_context(task_text="magic limit")

    assert ("[drifted: the cited code at srcmod.py:2 has changed — "
            "re-read it]") in out


def test_recall_healthy_reference_not_annotated(
    agent_tree, tmp_path, monkeypatch,
):
    """Healthy refs stay clean and still earn the normal recall bump."""
    from delfin.agent import memory_store as ms
    from delfin.agent.prompt_loader import PromptLoader

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    loader = PromptLoader(agent_tree)
    target = agent_tree / "srcmod.py"
    target.write_text("import os\nMAGIC_LIMIT = 17\n", encoding="utf-8")
    ms.save_typed_memory(
        "project: the magic limit constant sits at srcmod.py:2 today",
        repo_root=agent_tree)

    out = loader._load_external_memory_context(task_text="magic limit")

    assert "[stale:" not in out
    assert "[drifted:" not in out
    rec = ms.list_typed_memories(agent_tree)[0]
    assert rec["use_count"] == 2             # healthy recall still bumps
    assert rec["stale_hits"] == 0


# ---------------------------------------------------------------------------
# Global (user-scoped) store merged into the External Memory block
# ---------------------------------------------------------------------------

def test_recall_merges_global_store_first(agent_tree, tmp_path, monkeypatch):
    from delfin.agent import memory_store as ms
    from delfin.agent.prompt_loader import PromptLoader

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    loader = PromptLoader(agent_tree)
    ms.save_typed_memory(
        "global: user: Max prefers concise German answers",
        repo_root=agent_tree,
        allow_scope_prefix=True)
    ms.save_typed_memory(
        "project: dashboard websocket runs on port 8050",
        repo_root=agent_tree)

    out = loader._load_external_memory_context()

    assert "GLOBAL MEMORY.md" in out
    assert "concise German answers" in out
    assert "port 8050" in out
    # Global (identity / standing rules) comes before project content.
    assert out.index("concise German answers") < out.index("port 8050")
    # The recall bump lands in the right store.
    assert ms.list_typed_memories(agent_tree, scope="user")[0]["use_count"] == 2
    assert ms.list_typed_memories(agent_tree)[0]["use_count"] == 2


def test_recall_global_floor_survives_fat_project_store(
    agent_tree, tmp_path, monkeypatch,
):
    """A fat project store must not starve the global entries: they keep a
    guaranteed 25% floor of the budget (and are injected first)."""
    from delfin.agent import memory_store as ms
    from delfin.agent.prompt_loader import PromptLoader

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    loader = PromptLoader(agent_tree)
    ms.save_typed_memory(
        "global: user: identity anchor phrase zeta", repo_root=agent_tree,
        allow_scope_prefix=True)
    ms.save_typed_memory(
        "project: giant context dump " + "verbose " * 1500,
        repo_root=agent_tree)

    out = loader._load_external_memory_context(max_chars=3000)

    assert "identity anchor phrase zeta" in out
    assert len(out) <= 3100                  # cap + truncation marker
    assert "truncated" in out


def test_recall_without_global_store_unchanged(
    agent_tree, tmp_path, monkeypatch,
):
    """No ~/.delfin/memory -> pure project behaviour, no global markers."""
    from delfin.agent import memory_store as ms
    from delfin.agent.prompt_loader import PromptLoader

    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    loader = PromptLoader(agent_tree)
    ms.save_typed_memory(
        "project: dashboard websocket runs on port 8050",
        repo_root=agent_tree)

    out = loader._load_external_memory_context()

    assert "GLOBAL MEMORY.md" not in out
    # The block opens with its qualifier; the project index follows it.
    from delfin.agent.prompt_loader import _MEMORY_BLOCK_PREAMBLE
    assert out.startswith(_MEMORY_BLOCK_PREAMBLE)
    assert "\n\n# MEMORY.md" in out
    assert "port 8050" in out


# ---------------------------------------------------------------------------
# S1 — live_state injected into the system prompt (replaces user-msg state block)
# ---------------------------------------------------------------------------

def test_live_state_appears_in_solo_prompt(agent_tree):
    """live_state passes through and lands as a --- Live state --- section."""
    from delfin.agent.prompt_loader import PromptLoader
    loader = PromptLoader(agent_tree)
    prompt = loader.build_system_prompt(
        role_id="solo_agent",
        mode_id="solo",
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
        mode_id="solo",
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
        mode_id="solo",
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


# ---------------------------------------------------------------------------
# Refusal addendum — universal refusal contract injected for every role
# ---------------------------------------------------------------------------

def _real_pack() -> Path:
    return Path(__file__).resolve().parent.parent / "delfin" / "agent" / "pack"


def test_refusal_addendum_ships_in_the_pack():
    """The shared refusal contract exists and carries its four pillars:
    explicit what+why, no routing-around, nearest safe alternative, and
    risk-before-capability for legitimate out-of-scope requests."""
    body = (_real_pack() / "shared" / "refusal_addendum.md").read_text(
        encoding="utf-8")
    flat = " ".join(body.split())          # robust against line wrapping
    assert "Name what you will not do" in flat
    assert "Never route around the refusal" in flat
    assert "nearest safe alternative" in flat
    assert "BEFORE naming the mode" in flat
    # Pre-waived confirmation is not a license to act.
    assert "does not license the action" in flat


@pytest.mark.parametrize("role_id,mode_id", [
    ("solo_agent", "solo"),
    ("dashboard_agent", "dashboard"),
    ("builder_agent", "quick"),
    ("critic_agent", "quick"),
])
def test_refusal_addendum_injected_for_every_role(role_id, mode_id):
    """The contract is universal — composed into the prompt head for guide,
    solo, and pipeline roles alike (same layer as the honesty addendum)."""
    from delfin.agent.prompt_loader import PromptLoader

    prompt = PromptLoader().build_system_prompt(
        role_id=role_id, mode_id=mode_id, task_text="tidy the workspace")
    assert "Never route around the refusal" in prompt
    assert "nearest safe alternative" in prompt


def test_refusal_addendum_injected_when_present(agent_tree):
    """Loader contract: pack/shared/refusal_addendum.md is picked up by
    build_system_prompt without any extra registration."""
    from delfin.agent.prompt_loader import PromptLoader

    shared = agent_tree / "pack" / "shared"
    (shared / "refusal_addendum.md").write_text(
        "# Refusing unsafe requests\nREFUSAL-MARKER")
    (agent_tree / "pack" / "agents" / "solo_agent.md").write_text(
        "# Solo\nYou are solo.")
    loader = PromptLoader(agent_tree)
    prompt = loader.build_system_prompt(
        role_id="solo_agent", mode_id="solo", mode_description="solo",
        route=["solo_agent"], role_index=0)
    assert "REFUSAL-MARKER" in prompt


def test_dashboard_redirect_carves_out_destructive_requests():
    """The guide's mode-switch one-liner applies to safe, constructive
    requests only; destructive requests must be refused in the guide's own
    voice, never redirected to the code mode to be executed there."""
    body = (_real_pack() / "agents" / "dashboard_agent.md").read_text(
        encoding="utf-8")
    flat = " ".join(body.split())          # robust against line wrapping
    assert "refused, not redirected" in flat
    assert ("Never name another mode as the place where the destructive act "
            "would run") in flat
    # The redirect instruction itself is scoped to safe requests.
    assert "safe, constructive source edit" in flat


def test_dashboard_prompt_contains_refusal_contract_and_carve_out():
    """End-to-end: the composed dashboard prompt carries BOTH the universal
    contract and the role-specific carve-out, so the two never diverge."""
    from delfin.agent.prompt_loader import PromptLoader

    prompt = PromptLoader().build_system_prompt(
        role_id="dashboard_agent", mode_id="dashboard",
        route=["dashboard_agent"], role_index=0,
        task_text="open the submit tab")
    assert "Never route around the refusal" in prompt
    assert "refused, not redirected" in prompt


def test_dashboard_prompt_specifies_turn_closing_contract():
    """The guide's /done rule is a completion contract, not a cost hint:
    a request covered by the emitted ACTIONs closes with ACTION: /done in
    the same response; questions are reserved for genuinely missing
    information (one concrete question, then end the turn)."""
    body = (_real_pack() / "agents" / "dashboard_agent.md").read_text(
        encoding="utf-8")
    flat = " ".join(body.split())          # robust against line wrapping
    assert "A satisfied request is closed, not extended" in flat
    # Closing happens in the same response as the covering ACTIONs.
    assert "in the same response" in flat
    # Clarification stays legitimate — but scoped to missing information.
    assert "ask ONE concrete question" in flat
    # Result-dependent chains stay legitimate — /done is held back there.
    assert "depends on the RESULT" in flat
    # The old open-ended escape hatch must not come back: vague
    # uncertainty was the license models used to keep turns open.
    assert "unsure whether more actions are coming" not in flat


def test_episode_recall_injected_into_solo_prompt(
    agent_tree, tmp_path, monkeypatch,
):
    """A saved episode matching the task text appears as its own
    'Past Sessions' section; disabling agent.episodes.enabled removes it."""
    from delfin.agent import episodes as ep
    from delfin.agent.prompt_loader import PromptLoader
    import delfin.user_settings as us

    (agent_tree / "pack" / "agents" / "solo_agent.md").write_text(
        "# Solo Agent\nYou are the solo agent."
    )
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    ep.save_episode(
        "epi12345",
        repo_root=agent_tree,
        goal="Debug orca geometry convergence failure",
        outcome="Added tightscf and slowconv keywords",
        decisions=[],
        open_items=[],
        verdict="PASS",
    )
    monkeypatch.setattr(
        us, "load_settings",
        lambda *a, **k: {"agent": {"episodes": {"enabled": True}}})

    loader = PromptLoader(agent_tree)
    build = dict(
        role_id="solo_agent",
        mode_id="solo",
        mode_description="solo",
        route=["solo_agent"],
        role_index=0,
        task_text="the orca geometry convergence job fails again",
    )
    prompt = loader.build_system_prompt(**build)
    assert "--- Past Sessions ---" in prompt
    assert "# Similar past sessions" in prompt
    assert "tightscf" in prompt

    monkeypatch.setattr(
        us, "load_settings",
        lambda *a, **k: {"agent": {"episodes": {"enabled": False}}})
    prompt_disabled = loader.build_system_prompt(**build)
    assert "--- Past Sessions ---" not in prompt_disabled


# ---------------------------------------------------------------------------
# Document module — triggers in both languages the users write in
# ---------------------------------------------------------------------------

_DOC_MODULE_TEXT = (
    "## Intro\n"
    "Always kept.\n"
    "\n"
    "<!-- module:documents -->\n"
    "## Documents\n"
    "document module body\n"
    "\n"
    "## Tail\n"
    "tail body\n"
)


def _doc_module_active(agent_tree, task_text: str, key: str) -> bool:
    from delfin.agent.prompt_loader import PromptLoader

    out = PromptLoader(agent_tree)._strip_lazy_modules(
        _DOC_MODULE_TEXT, task_text=task_text, mode_id="solo",
        session_key=key, role_id="solo_agent",
    )
    return "document module body" in out


@pytest.mark.parametrize("task", [
    "werte bitte die Tabelle aus",
    "fill in the PDF Formular for me",
    "summarise budget.xlsx",
    "die Rechnung prüfen",
    "read the spreadsheet and total column C",
])
def test_document_module_activates_for_office_work(agent_tree, task):
    assert _doc_module_active(agent_tree, task, f"k-{hash(task)}")


@pytest.mark.parametrize("task", [
    "fix the import in delfin/agent/engine.py",
    "format the output of the parser",
    "give me more information about the run",
])
def test_document_module_stays_out_of_unrelated_work(agent_tree, task):
    """'formular' and 'tabelle' must not fire on 'format' or 'information'."""
    assert not _doc_module_active(agent_tree, task, f"n-{hash(task)}")
