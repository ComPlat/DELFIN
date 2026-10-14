"""AGENTS.md, CONTRIBUTING.md and CLAUDE.md reach every entry path once.

Input: a workspace with the three files and an engine (or a bare prompt
loader) rooted in it. Output: the built system prompt.

Before: the instruction-file loader had two callers, the dashboard turn
and the sub-agent spawn. The terminal CLI and ``-p`` runs -- which build
their prompt through ``AgentEngine._build_current_system_prompt`` --
never read AGENTS.md, although ``delfin-agent init`` writes one.
CONTRIBUTING.md and CLAUDE.md were read by nobody. Measured with the
foreign-repo probe: the agentsrepo prompt held none of its AGENTS.md.
"""
from __future__ import annotations

import pathlib
from pathlib import Path

from delfin.agent import untrusted

AGENTS = "AGENTS-RULE: run make check before each commit"
CONTRIB = "CONTRIB-RULE: open pull requests against develop"
OTHER = "OTHER-AGENT-RULE: never edit the gen directory"


def _ws(tmp_path: Path) -> Path:
    ws = tmp_path / "repo"
    ws.mkdir()
    (ws / "AGENTS.md").write_text(f"# A\n{AGENTS}\n", encoding="utf-8")
    (ws / "CONTRIBUTING.md").write_text(f"# C\n{CONTRIB}\n", encoding="utf-8")
    (ws / "CLAUDE.md").write_text(f"# O\n{OTHER}\n", encoding="utf-8")
    return ws


def _engine(ws: Path, mode: str = "solo"):
    from delfin.agent.engine import AgentEngine
    eng = AgentEngine(repo_dir=ws, backend="api", provider="kit",
                      api_key="stub-key", model="kit.glm-5.3", mode=mode,
                      permission_mode="default")
    eng.messages = [{"role": "user", "content": "fix the failing test"}]
    return eng


def test_the_terminal_engine_prompt_carries_all_three(tmp_path):
    prompt = _engine(_ws(tmp_path))._build_current_system_prompt(
        task_text="fix the failing test")
    for rule in (AGENTS, CONTRIB, OTHER):
        assert prompt.count(rule) == 1, rule


def test_they_arrive_fenced(tmp_path):
    prompt = _engine(_ws(tmp_path))._build_current_system_prompt(
        task_text="fix the failing test")
    start = prompt.index("UNTRUSTED EXTERNAL CONTENT — repository text")
    inner = untrusted.unwrap(prompt[start - 1:])
    for rule in (AGENTS, CONTRIB, OTHER):
        assert rule in inner


def test_the_section_is_byte_stable(tmp_path):
    eng = _engine(_ws(tmp_path))
    a = eng._build_current_system_prompt(task_text="fix the failing test")
    b = eng._build_current_system_prompt(task_text="fix the failing test")
    block = a[a.index("--- Project Instructions ---"):]
    block = block[:block.index("[END UNTRUSTED EXTERNAL CONTENT")]
    assert block in b


def test_a_pipeline_role_gets_them_too(tmp_path):
    from delfin.agent.prompt_loader import PromptLoader
    loader = PromptLoader()
    loader.workspace_root = _ws(tmp_path)
    loader.is_delfin_workspace = False
    prompt = loader.build_system_prompt(
        role_id="builder_agent", mode_id="solo", task_text="fix it")
    assert AGENTS in prompt and CONTRIB in prompt


def test_contributing_is_capped(tmp_path):
    from delfin.agent.project_memory import load_project_memory
    ws = _ws(tmp_path)
    (ws / "CONTRIBUTING.md").write_text(
        CONTRIB + "\n" + ("long text " * 2000), encoding="utf-8")
    out = load_project_memory(ws)
    assert CONTRIB in out
    assert len(untrusted.unwrap(out)) < 8000


def test_contributing_and_claude_md_are_read_from_the_workspace_only(
        tmp_path):
    """A CONTRIBUTING.md or CLAUDE.md in a PARENT directory belongs to some
    other project; only AGENTS.md / DELFIN.MD walk up."""
    from delfin.agent.project_memory import load_project_memory
    parent = tmp_path / "outer"
    ws = parent / "repo"
    ws.mkdir(parents=True)
    (parent / "CONTRIBUTING.md").write_text(CONTRIB, encoding="utf-8")
    (parent / "CLAUDE.md").write_text(OTHER, encoding="utf-8")
    (parent / "AGENTS.md").write_text(AGENTS, encoding="utf-8")
    out = load_project_memory(ws)
    assert AGENTS in out
    assert CONTRIB not in out and OTHER not in out


def test_the_sub_agent_inherits_contributing_too(tmp_path):
    from delfin.agent.subagents import _inherited_context
    text = _inherited_context(_ws(tmp_path), "do the task")
    assert CONTRIB in text


def test_one_loader_two_callers():
    """The dashboard turn no longer loads the files itself: the engine's
    prompt does, so a dashboard prompt does not carry them twice."""
    root = pathlib.Path(__file__).resolve().parents[1]
    callers = set()
    for path in sorted(root.joinpath("delfin").rglob("*.py")):
        text = path.read_text(encoding="utf-8", errors="replace")
        if "load_project_memory(" in text:
            callers.add(str(path.relative_to(root)))
    assert callers == {"delfin/agent/project_memory.py",
                       "delfin/agent/prompt_loader.py",
                       "delfin/agent/subagents.py"}, sorted(callers)


def test_no_workspace_no_section():
    """A loader nobody rooted in a workspace reads nothing (the budget
    tests build exactly that)."""
    from delfin.agent.prompt_loader import PromptLoader
    assert PromptLoader()._load_project_instructions() == ""
