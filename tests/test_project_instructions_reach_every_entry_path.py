"""AGENTS.md, CONTRIBUTING.md and CLAUDE.md reach every entry path once.

Input: a workspace with the three files and a real api/kit engine rooted
in it (no request sent). Output: what ``turn_context`` -- the one turn
composition the terminal session, ``-p``/run and the dashboard share --
puts into the system prompt, and what a sub-agent inherits.

Before: CONTRIBUTING.md and CLAUDE.md were read by nobody, so a
repository's branch, sign-off and test conventions never reached the
model; the instruction loader read a file whole whatever its size.
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


def _engine(ws: Path):
    from delfin.agent.engine import AgentEngine
    eng = AgentEngine(repo_dir=ws, backend="api", provider="kit",
                      api_key="stub-key", model="kit.glm-5.3", mode="solo",
                      permission_mode="default")
    eng.messages = [{"role": "user", "content": "fix the failing test"}]
    return eng


def _turn_prompt(eng) -> str:
    from delfin.agent import turn_context
    memory = turn_context.memory_context(eng, "fix the failing test")
    return eng._build_current_system_prompt(
        memory_context=memory, task_text="fix the failing test")


def test_the_terminal_turn_carries_all_three_once(tmp_path):
    prompt = _turn_prompt(_engine(_ws(tmp_path)))
    for rule in (AGENTS, CONTRIB, OTHER):
        assert prompt.count(rule) == 1, rule


def test_they_arrive_fenced(tmp_path):
    prompt = _turn_prompt(_engine(_ws(tmp_path)))
    start = prompt.index("[UNTRUSTED EXTERNAL CONTENT — repository text")
    inner = untrusted.unwrap(prompt[start:])
    for rule in (AGENTS, CONTRIB, OTHER):
        assert rule in inner


def test_the_block_is_byte_stable_between_turns(tmp_path):
    from delfin.agent import turn_context
    eng = _engine(_ws(tmp_path))
    assert (turn_context.project_instructions(eng)
            == turn_context.project_instructions(eng))


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


def test_a_file_linking_out_of_the_workspace_is_not_read(tmp_path):
    from delfin.agent.project_memory import load_project_memory
    ws = tmp_path / "repo"
    ws.mkdir()
    secret = tmp_path / "key.txt"
    secret.write_text("PRIVATE-KEY-MATERIAL", encoding="utf-8")
    (ws / "AGENTS.md").symlink_to(secret)
    assert "PRIVATE-KEY-MATERIAL" not in load_project_memory(ws)


def test_the_sub_agent_inherits_contributing_too(tmp_path):
    from delfin.agent.subagents import _inherited_context
    text = _inherited_context(_ws(tmp_path), "do the task")
    assert CONTRIB in text


def test_one_loader_two_callers():
    """turn_context for the three entry paths, the sub-agent spawn for a
    child; nothing else reads the files on its own."""
    root = pathlib.Path(__file__).resolve().parents[1]
    callers = set()
    for path in sorted(root.joinpath("delfin").rglob("*.py")):
        text = path.read_text(encoding="utf-8", errors="replace")
        if "load_project_memory(" in text:
            callers.add(str(path.relative_to(root)))
    assert callers == {"delfin/agent/project_memory.py",
                       "delfin/agent/turn_context.py",
                       "delfin/agent/subagents.py"}, sorted(callers)
