"""What the agent records about itself becomes what it reads about itself.

Two entries found in the user's real outcome history on 2026-09-11:
a FAILED task whose text was "[Verify] The following physical quantities
were stated without any evidence ..." -- the verifier's own correction
turn, recorded as if the user had asked it, and then quoted back into
later prompts by the briefing as a lesson from the user's past.

And a rule the same prompts carried for every solo session, under
"Directory permissions": "Never run real ORCA/xTB/SLURM -- only pytest."
Both models, interviewed as operators, named it as a contradiction with
the submit and pipeline tools in the same prompt. It was written for
work on DELFIN's own code, and that is where it now stands; what a
calculation may not do is start from bash instead of the sanctioned path.
"""

from __future__ import annotations

import textwrap
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def agent_tree(tmp_path):
    """The minimal pack tree the engine tests use."""
    agent_dir = tmp_path / "pack"
    shared = agent_dir / "shared"
    shared.mkdir(parents=True)
    agents = agent_dir / "agents"
    agents.mkdir()
    for name in ("delfin_context.md", "work_cycle_rules.md",
                 "goal_decomposition_rules.md", "universal_input_template.md",
                 "minimal_final_verdict.md"):
        (shared / name).write_text("# x")
    for name in ("session_manager.md", "builder_agent.md", "test_agent.md"):
        (agents / name).write_text("# role")
    lite = tmp_path / "pack_lite"
    (lite / "modes").mkdir(parents=True)
    (lite / "modes" / "solo.md").write_text("# solo mode")
    (lite / "manifest.yaml").write_text(textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - session_manager
    """))
    return tmp_path


@pytest.fixture
def engine(agent_tree):
    from delfin.agent.engine import AgentEngine
    client = MagicMock()
    client.model = "test-model"
    with patch("delfin.agent.engine.create_client", return_value=client):
        eng = AgentEngine(repo_dir=agent_tree, backend="cli", mode="solo",
                          pack_dir=agent_tree)
    return eng


def _history_lines():
    from delfin.agent import outcome_tracker
    p = outcome_tracker._DEFAULT_PATH
    return p.read_text().splitlines() if p.exists() else []


def test_a_verify_turn_writes_no_outcome(engine):
    before = len(_history_lines())
    engine.record_cycle_outcome(
        "FAIL", "[Verify] The following physical quantities were stated "
                "without any evidence act this turn: '-25 Eh'.",
        error_type="empty_turn")
    assert len(_history_lines()) == before


def test_a_users_task_still_writes_one(engine):
    before = len(_history_lines())
    engine.record_cycle_outcome("PASS", "Welche Rechnung hat die niedrigste Energie?")
    lines = _history_lines()
    assert len(lines) == before + 1
    assert "niedrigste Energie" in lines[-1]


def _built_solo_prompt() -> str:
    from delfin.agent.prompt_loader import PromptLoader
    loader = PromptLoader()
    loader.workspace_root = _ROOT
    return loader.build_system_prompt(
        role_id="solo_agent", mode_id="solo",
        task_text="Welche meiner Rechnungen ist noch nicht fertig?",
        session_key="internal-turn-prompt", model="kit.deepseek-v4-flash")


def test_the_solo_prompt_no_longer_forbids_its_own_tools():
    prompt = _built_solo_prompt()
    assert "only pytest" not in prompt
    assert "never by launching ORCA, xTB or SLURM from bash" in prompt
    # The rule for DELFIN's own code survives, where it belongs.
    assert "never launch a real ORCA, xTB or SLURM job to test a change" in prompt


def test_the_rule_moved_rather_than_vanished():
    text = (_ROOT / "delfin" / "agent" / "pack" / "agents" / "solo_agent.md").read_text()
    perms = text.split("## Directory permissions", 1)[1].split("\n## ", 1)[0]
    assert "pytest" not in perms
    edits = text.split("## After every code edit", 1)[1].split("\n## ", 1)[0]
    assert "never launch a real ORCA" in edits
