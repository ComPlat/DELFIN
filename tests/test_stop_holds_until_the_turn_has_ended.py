"""Stop stays a button until the turn has really ended.

Pressing Stop while the model sent nothing yet (a queued KIT endpoint)
set the tab to "not streaming" at once, which disabled Stop -- while the
turn's thread still ran. A stop that missed (the flag cleared by the turn
starting, a request that did not read it) could then not be repeated, and
the model carried on with nothing on screen to stop it.

Now Stop is enabled for as long as the turn's thread is alive, every press
stops again, and one press repeats the stop until that thread has ended.
"""

from __future__ import annotations

import threading
import time

import pytest


class _Client:
    def __init__(self):
        self.signals = 0

    def signal_stop(self):
        self.signals += 1


class _Engine:
    def __init__(self):
        self.stops = 0
        self.client = _Client()

    def request_stop(self):
        self.stops += 1


def _tab(tmp_path):
    from delfin.agent import scheduler as S
    from delfin.dashboard import tab_agent
    from delfin.dashboard.context import DashboardContext

    S._GLOBAL = S.Scheduler(path=tmp_path / "cron.json")
    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir(exist_ok=True)
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script, **kw: None
    root, refs = tab_agent.create_tab(ctx)
    return root, refs


def _stop_button(root):
    import ipywidgets as widgets
    stack = [root]
    while stack:
        n = stack.pop()
        if (isinstance(n, widgets.Button) and n.description == "Stop"
                and n.button_style == "danger"):
            return n
        stack.extend(getattr(n, "children", ()) or ())
    raise AssertionError("no Stop button")


def _wait(pred, seconds=5.0):
    end = time.monotonic() + seconds
    while time.monotonic() < end:
        if pred():
            return True
        time.sleep(0.02)
    return pred()


@pytest.fixture
def running_turn(tmp_path, monkeypatch):
    from delfin.dashboard import tab_agent
    monkeypatch.setattr(tab_agent, "_STOP_REPEAT_S", 0.05, raising=False)
    root, refs = _tab(tmp_path)
    state = refs["state"]
    release = threading.Event()
    thread = threading.Thread(target=release.wait, args=(30,), daemon=True)
    thread.start()
    engine = _Engine()
    state["engine"] = engine
    state["streaming"] = True
    state["_turn_thread"] = thread
    yield root, state, engine, release, thread
    release.set()
    try:
        refs["shutdown"]()
    except Exception:
        pass


def test_stop_stays_pressable_while_the_turn_still_runs(running_turn):
    root, state, engine, release, thread = running_turn
    stop = _stop_button(root)
    stop.disabled = False
    stop.click()
    assert state["streaming"] is False
    assert engine.stops >= 1 and engine.client.signals >= 1
    # The thread has not ended: the button is still there to press.
    assert stop.disabled is False
    before = engine.stops
    stop.click()
    assert engine.stops > before


def test_one_press_repeats_the_stop_until_the_turn_ends(running_turn):
    root, state, engine, release, thread = running_turn
    stop = _stop_button(root)
    stop.disabled = False
    stop.click()
    assert _wait(lambda: engine.stops >= 4), engine.stops
    release.set()
    thread.join(5)
    # The turn is over: the repeating ends and Stop is greyed out again.
    assert _wait(lambda: stop.disabled is True)
    settled = engine.stops
    time.sleep(0.3)
    assert engine.stops == settled


@pytest.mark.parametrize("content,kept", [("", False), ("  \n", False),
                                          ("Erst lese ich die Datei", True)])
def test_a_stop_before_any_text_leaves_no_empty_answer(tmp_path, content, kept):
    """Field report 20261020-140059: Stop while the model had sent nothing
    left an empty "Code Agent" bubble in the chat; a partial answer stays."""
    root, refs = _tab(tmp_path)
    state = refs["state"]
    try:
        state["engine"] = _Engine()
        state["streaming"] = True
        state["chat_messages"].append({"role": "user", "content": "go"})
        state["chat_messages"].append(
            {"role": "assistant", "content": content, "_streaming": True})
        stop = _stop_button(root)
        stop.disabled = False
        stop.click()
        answers = [m for m in state["chat_messages"] if m.get("role") == "assistant"]
        assert bool(answers) is kept
        assert not any(m.get("_streaming") for m in state["chat_messages"])
    finally:
        refs["shutdown"]()
