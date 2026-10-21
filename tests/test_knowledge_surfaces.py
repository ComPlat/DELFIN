"""The knowledge commands reach every surface, and only the person can run
the ones that share or curate.
"""
from __future__ import annotations

import inspect
from pathlib import Path
from types import SimpleNamespace

import pytest

from delfin.agent import knowledge_private as kp
from delfin.agent import repl_commands as rc

NEW = ("/memory", "/share", "/unshare", "/lessons", "/swarm", "/knowledge")
LESSON = ("Situation: a parallel pytest run reports failures. Lesson: "
          "re-run them serially first. Why: shared state.")


@pytest.fixture(autouse=True)
def _home(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.delenv("DELFIN_AGENT_SESSION", raising=False)
    monkeypatch.setenv("DELFIN_KNOWLEDGE_QUEUE", str(tmp_path / "q"))


def test_the_terminal_registers_every_knowledge_command():
    for name in NEW:
        assert name in rc.BUILTINS, name
        assert rc.BUILTINS[name].summary.strip()
    rows = {r[1] for r in rc.palette_rows()}
    assert set(NEW) <= rows


def test_memory_and_share_dispatch_through_the_repl_router(tmp_path):
    ctx = rc.ReplContext(engine=None, workspace=tmp_path, ask=None,
                         say=lambda _s: None)
    out = rc.dispatch("/memory", ctx)
    assert out.handled and "layers:" in out.output
    les = kp.propose(LESSON, evidence=["test:t"], source="agent",
                     topic="tools-and-sandbox", title="Serial first")
    out = rc.dispatch(f"/share {les.id} global as maxh", ctx)
    assert out.handled and out.prompt == ""          # nothing to the model
    assert "/share confirm" in out.output or "refused" in out.output


def test_the_window_lists_their_arguments():
    from delfin.agent.service import choices
    for name in NEW:
        assert name in choices._USAGE or name in choices._WORDS, name


def test_the_dashboard_routes_them_and_refuses_them_from_the_agent():
    from delfin.dashboard import tab_agent as T
    for name in NEW:
        assert name in T._BUILTIN_SLASH_PREFIXES, name
        assert any(row[1].split()[0] == name for row in T._SLASH_COMMANDS)
    for cmd in ("/share 3 global", "/share confirm abc ABC-DEF",
                "/unshare entry", "/lessons edit 1 lesson=x",
                "/swarm close wave-1", "/knowledge retire x because",
                "/knowledge accept q1"):
        assert T._slash_is_user_only(cmd), cmd
    assert not T._slash_is_user_only("/memory")


def test_the_dashboard_hands_the_original_text_to_the_handler():
    """The router lower-cases its own copy; a lesson edit or a code must
    reach the handler as typed."""
    from delfin.dashboard import tab_agent as T
    src = Path(inspect.getfile(T)).read_text(encoding="utf-8")
    block = src[src.index("# Layered memory (knowledge_commands)"):]
    block = block[:block.index("return True")]
    assert "text.strip()" in block and "knowledge_commands" in block


def test_the_repl_shows_the_share_hint_once_after_a_turn():
    from delfin.agent import repl
    kp.propose(LESSON, evidence=["test:t"], source="agent", title="Serial")
    lines: list = []
    agent = repl.TerminalAgent.__new__(repl.TerminalAgent)
    agent.transcript = SimpleNamespace(
        chrome=lines.append,
        theme=SimpleNamespace(dim=lambda s: s, yellow=lambda s: s))
    agent._report_lessons()
    agent._report_lessons()
    hints = [ln for ln in lines if "could help others" in ln]
    assert len(hints) == 1 and "/share 1" in hints[0]


def test_the_cli_has_a_knowledge_command(tmp_path, monkeypatch, capsys):
    from delfin.agent import cli
    from delfin.agent import knowledge as kn
    gdir = tmp_path / "k"
    kn.write_topics(gdir, {t: "" for t in kn.DEFAULT_TOPICS})
    kn.write_index(gdir)
    monkeypatch.setattr(kn, "global_dir", lambda: gdir)
    rc_ = cli.main(["knowledge", "lint"])
    assert rc_ == 0 and "lint: clean" in capsys.readouterr().out
