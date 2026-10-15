"""The dashboard's agent tab drives auto mode and shows its state.

Input: a real agent tab (no browser; ``create_tab`` with a scratch
context), a stand-in engine with a permission object, typed commands sent
through the input box and the Send button, and the tab's own tick, turn-end
hook and buttons. Output: chat messages, the auto panel, the toggle, the
scheduler entries and the controller's state.

Judged: /auto and the /loop sub-commands reach the shared handler; the old
``/loop <interval> <prompt>`` stays a scheduler entry; both commands are
the person's only (registered user-only, refused from a turn nobody
typed); the tick hands the controller's prompt to the turn starter and
takes the iteration back when the turn cannot start; a message typed
while auto mode runs pauses it; the panel shows the state and the needs-you
queue, and its button stops auto mode.
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

widgets = pytest.importorskip("ipywidgets")

from delfin.agent import agent_tasks  # noqa: E402
from delfin.agent import scheduler as S  # noqa: E402
from delfin.dashboard import tab_agent  # noqa: E402


class _Engine:
    session_id = "dash-s1"
    mode = "solo"

    def __init__(self, ws):
        self.kit_permissions = SimpleNamespace(mode="acceptEdits",
                                               auto_veto=None,
                                               workspace=str(ws))
        self.token_usage = {"input": 0, "output": 0}
        self.messages = []
        self.pending_turn_continuation = ""
        # A live process: the send path's stale-turn check finds one.
        self.client = SimpleNamespace(
            model="m", _proc=SimpleNamespace(poll=lambda: None))

    def request_stop(self):
        pass

    def note_typed(self, text):
        pass


def _git(ws, *a):
    subprocess.run(["git", "-C", str(ws), *a], check=True,
                   capture_output=True)


@pytest.fixture()
def tab(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))
    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir()
    ws = tmp_path / "ws"
    ws.mkdir()
    _git(ws, "init", "-q", "-b", "main")
    _git(ws, "config", "user.name", "t")
    _git(ws, "config", "user.email", "t@example.org")
    (ws / "f.txt").write_text("x\n", encoding="utf-8")
    _git(ws, "add", "-A")
    _git(ws, "commit", "-q", "-m", "i")
    from delfin.dashboard.context import DashboardContext
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script, **kw: None
    ctx.agent_workspace = str(ws)
    t, refs = tab_agent.create_tab(ctx)
    eng = _Engine(ws)
    refs["state"]["engine"] = eng
    yield t, refs, eng, ws
    refs["shutdown"]()


def _walk(root):
    stack = [root]
    while stack:
        node = stack.pop()
        yield node
        stack.extend(getattr(node, "children", ()) or ())


def _type(t, text):
    box = next(w for w in _walk(t) if isinstance(w, widgets.Textarea)
               and "delfin-agent-input" in w._dom_classes)
    send = next(w for w in _walk(t) if isinstance(w, widgets.Button)
                and w.description == "Send")
    box.value = text
    send.click()


def _chat(refs):
    return "\n".join(str(m.get("content", ""))
                     for m in refs["state"]["chat_messages"])


def _button(t, description):
    return next(w for w in _walk(t) if isinstance(w, widgets.Button)
                and w.description == description)


def _panel(t):
    return next(w for w in _walk(t) if isinstance(w, widgets.HTML)
                and "delfin-auto-panel" in w._dom_classes)


def test_the_registers_know_both_commands():
    assert "/auto" in tab_agent._BUILTIN_SLASH_PREFIXES
    assert "/auto" in tab_agent._USER_ONLY_SLASH_COMMANDS
    palette = {row[1] for row in tab_agent._SLASH_COMMANDS}
    assert {"/auto", "/loop"} <= palette


def test_auto_on_shows_the_panel_and_its_button_stops_it(tab):
    t, refs, eng, ws = tab
    _type(t, "/auto on")
    ctl = refs["auto"]()
    assert ctl is not None and ctl.enabled
    assert "auto mode on" in _chat(refs)
    refs["refresh_loop"]()
    panel = _panel(t)
    assert "auto 0/" in panel.value
    assert _button(t, "Auto: on")
    _button(t, "Stop auto").click()
    assert not ctl.enabled and "auto stopped" in _chat(refs)
    refs["refresh_loop"]()
    assert _panel(t).value == ""
    assert _button(t, "Auto: off")


def test_the_toggle_turns_it_on(tab):
    t, refs, eng, ws = tab
    _button(t, "Auto: off").click()
    assert refs["auto"]().enabled


def test_loop_commands_reach_the_shared_handler(tab):
    t, refs, eng, ws = tab
    _type(t, "/loop")
    assert "no loops defined" in _chat(refs)
    d = ws / ".delfin" / "loops"
    d.mkdir(parents=True)
    (d / "tidy.md").write_text("---\ngoal: tidy notes\n---\n1. Tidy one.\n",
                               encoding="utf-8")
    _type(t, "/loop show tidy")
    assert "1. Tidy one." in _chat(refs)
    _type(t, "/loop run tidy")
    assert "not confirmed" in _chat(refs)
    _type(t, "/loop confirm tidy")
    _type(t, "/loop run tidy")
    assert refs["auto"]().loop is not None
    _type(t, "/loop stop")
    assert refs["auto"]().loop is None


def test_the_interval_form_stays_a_scheduler_entry(tab):
    t, refs, eng, ws = tab
    _type(t, "/loop 30m /check")
    entries = [e for e in S._GLOBAL.list_entries()
               if (e.reason or "").startswith("loop:")]
    assert len(entries) == 1 and entries[0].every_seconds == 1800
    _type(t, "/loop list")
    assert "Active loops" in _chat(refs)


def test_a_turn_nobody_typed_cannot_switch_it_on(tab):
    t, refs, eng, ws = tab
    refs["state"]["_on_its_own"] = True
    try:
        _type(t, "/auto on")
    finally:
        refs["state"]["_on_its_own"] = False
    ctl = refs["auto"]()
    assert ctl is None or not ctl.enabled
    assert "only from the person" in _chat(refs)


def test_the_tick_starts_an_iteration_and_takes_it_back_when_held(tab):
    t, refs, eng, ws = tab
    store = agent_tasks.get_store(ws)
    task = store.create("write the parser", session_id=eng.session_id)
    store.update(task["id"], status="in_progress")
    _type(t, "/auto on")
    ctl = refs["auto"]()
    ctl.settings["idle_seconds"] = 0
    sent = []
    prompt = refs["auto_tick"](send=lambda text: sent.append(text) or True)
    assert prompt and sent and "write the parser" in sent[0]
    assert sent[0].startswith("[auto")
    assert ctl.iteration is not None
    # The turn ended: judged, digest in the chat.
    refs["state"]["_last_answer_text"] = "LOOP: progress parser started"
    refs["auto_after_turn"]()
    assert ctl.iteration is None and "iteration 1" in _chat(refs)
    # A held turn takes the iteration back.
    refs["auto_tick"](send=lambda text: False)
    assert ctl.iteration is None and ctl.state()["iteration"] == 1


def test_a_message_typed_while_auto_runs_pauses_it(tab):
    t, refs, eng, ws = tab
    _type(t, "/auto on")
    ctl = refs["auto"]()
    refs["state"]["streaming"] = True      # a turn runs: the line is queued
    try:
        _type(t, "what are you doing?")
    finally:
        refs["state"]["streaming"] = False
    assert ctl.paused
    assert refs["state"]["message_queue"] == ["what are you doing?"]


def test_the_panel_lists_what_needs_the_user(tab):
    t, refs, eng, ws = tab
    _type(t, "/auto on")
    ctl = refs["auto"]()
    ctl._queue("outward", "bash: git push origin HEAD", "pushes to a remote")
    refs["refresh_loop"]()
    value = _panel(t).value
    assert "needs you" in value and "git push" in value
    assert "<script" not in value
