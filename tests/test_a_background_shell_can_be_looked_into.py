"""A background shell in the panel can be looked into: its last output lines
land in the chat, the way a background task's output can be opened.

Asked for on 2026-09-16: the row showed only the command and how long it had
run; the output was reachable through the agent alone (bash_output).
"""
import inspect
import json
from pathlib import Path

from delfin.agent import background_view as BV
from delfin.agent import bash_jobs as bj


def test_a_shell_shows_the_tail_of_its_output(tmp_path, monkeypatch):
    out = tmp_path / "kit_bg_job.stdout"; out.write_text("\n".join(f"line {i}" for i in range(100)) + "\n")
    err = tmp_path / "kit_bg_job.stderr"; err.write_text("warning: x\n")
    reg = {"jobs": {"bg-1": {"command": "pytest -q", "stdout_path": str(out), "stderr_path": str(err)}}}
    monkeypatch.setattr(bj, "_load_registry_file", lambda ws: reg)
    text = BV.peek(tmp_path, "shells", "bg-1", lines=5)
    assert text.startswith("Shell bg-1 · pytest -q")
    assert "line 99" in text and "line 94" not in text
    assert "[stderr, last 1 line(s)]" in text and "warning: x" in text


def test_an_unknown_shell_says_so(tmp_path, monkeypatch):
    monkeypatch.setattr(bj, "_load_registry_file", lambda ws: {"jobs": {}})
    assert BV.peek(tmp_path, "shells", "nope") == "Shell nope: no record of it."


def test_the_other_groups_say_where_to_look(tmp_path):
    assert "result arrives in the chat" in BV.peek(tmp_path, "watches", "ci:x@y")
    assert "fires as a message" in BV.peek(tmp_path, "wakeups", "w1")


def test_the_panel_row_has_the_look_inside_button():
    src = Path(inspect.getfile(__import__("delfin.dashboard.tab_agent", fromlist=["x"]))).read_text()
    i = src.index("def _refresh_subagent_panel(")
    body = src[i:i + 6000]
    assert 'if _row["group"] in ("shells", "agents"):' in body
    assert "_peek_background(_w, *_k)" in body
    assert "def _peek_background(" in src
