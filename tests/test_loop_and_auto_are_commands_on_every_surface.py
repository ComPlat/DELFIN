"""/loop and /auto: one handler, reached from the terminal table.

Input: a typed line through ``repl_commands.dispatch`` with a
``ReplContext`` whose ``auto`` is the session's AutoController. Output:
the CommandResult. The window runs the same table (its worker is the
terminal loop), and the dashboard calls the same ``loop_commands``
functions, so the words are the same everywhere.
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

from delfin.agent import auto_mode, loop_commands, loops
from delfin.agent import repl_commands as rc


class Engine:
    def __init__(self, mode="acceptEdits"):
        self.kit_permissions = SimpleNamespace(mode=mode, auto_veto=None)
        self.token_usage = {"input": 0, "output": 0}
        self.messages = [{"role": "user", "content": "evaluate job 4711"},
                         {"role": "assistant", "content": "done"}]
        self.session_id = "s1"

    def request_stop(self):
        pass


@pytest.fixture()
def ws(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
    w = tmp_path / "ws"
    w.mkdir()
    for args in (["init", "-q", "-b", "main"], ["config", "user.name", "t"],
                 ["config", "user.email", "t@example.org"]):
        subprocess.run(["git", "-C", str(w), *args], check=True)
    (w / "f.txt").write_text("x\n", encoding="utf-8")
    subprocess.run(["git", "-C", str(w), "add", "-A"], check=True)
    subprocess.run(["git", "-C", str(w), "commit", "-q", "-m", "i"], check=True)
    return w


def _ctx(ws, answers=()):
    eng = Engine()
    ctl = auto_mode.AutoController(ws, session_id="s1",
                                   run_check=lambda e, c: (1, ""),
                                   notify=lambda *a: None,
                                   settings={"idle_seconds": 60})
    said = []
    queue = list(answers)
    ctx = rc.ReplContext(engine=eng, workspace=ws, session_id="s1",
                         ask=(lambda p: queue.pop(0) if queue else ""),
                         say=said.append, auto=ctl)
    return ctx, ctl, said


def _run(ctx, line):
    res = rc.dispatch(line, ctx)
    assert res.handled, line
    return res


def test_both_are_in_the_table_and_in_help():
    assert "/loop" in rc.BUILTINS and "/auto" in rc.BUILTINS
    names = [r[1] for r in rc.palette_rows()]
    assert "/loop" in names and "/auto" in names


def test_auto_on_off_and_status(ws):
    ctx, ctl, _ = _ctx(ws)
    assert "off" in _run(ctx, "/auto").output
    out = _run(ctx, "/auto on").output
    assert ctl.enabled and "auto mode on" in out
    assert "auto" in _run(ctx, "/auto").output
    assert "stopped" in _run(ctx, "/auto off").output and not ctl.enabled


def test_auto_on_refuses_bypass(ws):
    ctx, ctl, _ = _ctx(ws)
    ctx.engine.kit_permissions.mode = "bypassPermissions"
    assert "bypass" in _run(ctx, "/auto on").output.lower()
    assert not ctl.enabled


def test_a_surface_without_a_controller_says_so(ws):
    ctx = rc.ReplContext(engine=Engine(), workspace=ws)
    assert "not available" in _run(ctx, "/auto on").output
    assert "not available" in _run(ctx, "/loop run x").output


def test_loop_list_show_confirm_run_stop(ws):
    ctx, ctl, _ = _ctx(ws)
    assert "no loops" in _run(ctx, "/loop list").output.lower()
    d = ws / ".delfin" / "loops"
    d.mkdir(parents=True)
    (d / "tidy.md").write_text("---\ngoal: tidy notes\n---\n1. Tidy one.\n",
                               encoding="utf-8")
    listing = _run(ctx, "/loop").output
    assert "tidy" in listing and "not confirmed" in listing
    shown = _run(ctx, "/loop show tidy").output
    assert "1. Tidy one." in shown and "measure" in shown
    # run without confirmation and nobody to ask: refused, nothing started
    ctx.ask = None
    assert "confirm" in _run(ctx, "/loop run tidy").output
    assert ctl.loop is None
    assert "confirmed" in _run(ctx, "/loop confirm tidy").output
    assert loops.is_confirmed(loops.load("tidy", ws))
    out = _run(ctx, "/loop run tidy").output
    assert "started" in out and ctl.loop is not None
    assert "tidy" in _run(ctx, "/loop").output
    assert "stopped" in _run(ctx, "/loop stop").output and ctl.loop is None


def test_run_asks_to_confirm_what_it_shows(ws):
    d = ws / ".delfin" / "loops"
    d.mkdir(parents=True)
    (d / "tidy.md").write_text("---\ngoal: tidy notes\n---\n1. Tidy one.\n",
                               encoding="utf-8")
    ctx, ctl, said = _ctx(ws, answers=["n"])
    out = _run(ctx, "/loop run tidy").output
    assert ctl.loop is None and "not started" in out
    assert any("1. Tidy one." in s for s in said)
    ctx, ctl, said = _ctx(ws, answers=["y"])
    out = _run(ctx, "/loop run tidy").output
    assert ctl.loop is not None and "started" in out


def test_an_invalid_loop_says_why(ws):
    d = ws / ".delfin" / "loops"
    d.mkdir(parents=True)
    (d / "bad.md").write_text("---\ngoal: x\nbranch: main\n---\n1. a\n",
                              encoding="utf-8")
    ctx, ctl, _ = _ctx(ws)
    assert "branch" in _run(ctx, "/loop show bad").output
    assert "branch" in _run(ctx, "/loop confirm bad").output
    assert not loops.is_confirmed(loops.load("bad", ws))


def test_new_from_the_conversation_asks_the_model_to_draft(ws):
    ctx, ctl, _ = _ctx(ws)
    res = _run(ctx, "/loop new job-table")
    assert res.prompt and ".delfin/loops/job-table.md" in res.prompt
    assert "loops" in res.prompt          # the skill it follows
    assert "confirm" in res.prompt.lower()
    assert not (ws / ".delfin" / "loops" / "job-table.md").exists()


def test_new_blank_writes_an_unconfirmed_template(ws):
    ctx, ctl, _ = _ctx(ws)
    out = _run(ctx, "/loop new weekly --blank").output
    assert (ws / ".delfin" / "loops" / "weekly.md").is_file()
    assert "/loop show weekly" in out
    assert not loops.is_confirmed(loops.load("weekly", ws))


def test_an_interval_line_is_an_ad_hoc_loop(ws):
    ctx, ctl, _ = _ctx(ws)
    out = _run(ctx, "/loop 30m run the suite and fix the first failure").output
    assert ctl.loop is not None and ctl.loop.spec.every_seconds == 1800
    assert "started" in out


def test_needs_you_lists_the_queue(ws):
    ctx, ctl, _ = _ctx(ws)
    assert "nothing" in _run(ctx, "/auto needs").output.lower()
    ctl._queue("outward", "bash: git push", "pushes to a remote")
    out = _run(ctx, "/auto needs").output
    assert "#1" in out and "git push" in out
    _run(ctx, "/auto needs clear")
    assert ctl.needs_you == []


def test_the_skill_ships_and_names_the_format():
    path = Path(loop_commands.__file__).parent / "pack" / "skills" / "loops.md"
    text = path.read_text(encoding="utf-8")
    assert text.startswith("---\nname: loops\n")
    for word in ("trigger", "budget", "measure", "LOOP:", "/loop confirm"):
        assert word in text
