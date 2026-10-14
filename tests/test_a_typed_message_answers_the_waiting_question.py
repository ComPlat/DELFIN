"""A message typed into the chat answers the choice question that waits.

ask_user_question blocks the tool thread until an option is picked. A
line typed into the normal chat input instead went past it: the
dashboard steered it into the running turn (which could not read it --
the turn was blocked inside the question) and left the question on
screen until its 300 s expiry; the terminal dialog read single digit
keys and dropped every other key. Now, on both surfaces, a line sent
while a question waits IS the answer: the question is resolved with the
typed text (an option when it names one, free text otherwise), its
widget is closed, and the line starts no second turn.
"""

from __future__ import annotations

import io
import json
import tempfile
import threading
import time
from pathlib import Path
from types import SimpleNamespace

import pytest

from delfin.agent import approval_answers as aa
from delfin.agent import repl
from delfin.agent import repl_render as rr
from delfin.agent import terminal_confirm as tc
from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor

PAYLOAD = {
    "question": "Aufgabe — Was genau soll gebaut werden?",
    "header": "Aufgabe",
    "options": [{"label": "Ein CLI-Werkzeug", "description": ""},
                {"label": "Eine Bibliothek", "description": ""},
                {"label": "Both", "description": ""}],
    "multiSelect": False,
}
TYPED = "Ein Dashboard-Panel für Konformere"


# ---------------------------------------------------------------------------
# The shared answer shape and the tool result
# ---------------------------------------------------------------------------

def test_free_text_is_kept_and_marked_as_typed():
    out = aa.typed_answer(TYPED, PAYLOAD)
    assert out["answers"] == [TYPED]
    assert out["answered_in_chat"] is True
    assert out["free_text"] is True


@pytest.mark.parametrize("typed,label", [("both", "Both"), ("2", "Eine Bibliothek"),
                                         (" Ein CLI-Werkzeug ", "Ein CLI-Werkzeug")])
def test_a_typed_option_or_number_is_that_option(typed, label):
    out = aa.typed_answer(typed, PAYLOAD)
    assert out["answers"] == [label]
    assert out["answered_in_chat"] is True
    assert not out.get("free_text")


def test_the_tool_result_tells_the_model_it_was_typed():
    perms = KitToolPermissions(
        workspace=Path(tempfile.mkdtemp(prefix="askchat_")), mode="default",
        ask_user_callback=lambda args: aa.typed_answer(TYPED, args))
    out = json.loads(_DocToolExecutor().execute(
        "ask_user_question", dict(PAYLOAD), permissions=perms))
    assert out["answers"] == [TYPED]
    assert out["answered_in_chat"] is True
    assert out["free_text"] is True
    assert "chat" in out["note"] and "outside the options" in out["note"]
    assert not out.get("timed_out")


# ---------------------------------------------------------------------------
# Terminal: the broker and the dialog
# ---------------------------------------------------------------------------

@pytest.fixture
def rooms(tmp_path, monkeypatch):
    monkeypatch.setattr(tc, "_PENDING_DIR", tmp_path / "terminal")
    return tmp_path


def test_the_terminal_broker_takes_a_typed_answer(rooms):
    broker = tc.TerminalConfirmBroker(timeout_s=10.0, poll_s=0.02)
    out: dict = {}
    t = threading.Thread(
        target=lambda: out.setdefault("r", broker.ask_user(dict(PAYLOAD))),
        daemon=True)
    t.start()
    for _ in range(500):
        if broker.take() is not None:
            break
        time.sleep(0.01)
    assert broker.take() is not None, "the question never arrived"
    assert broker.answer_question_with_text(TYPED) is True
    t.join(5)
    assert not t.is_alive(), "the asking thread still waits"
    assert out["r"]["answers"] == [TYPED]
    assert out["r"]["answered_in_chat"] is True
    assert broker.take() is None, "the question is still pending"
    # Nothing left to answer: a second line is a message, not an answer.
    assert broker.answer_question_with_text("noch etwas") is False


class _Keys:
    active = True

    def __init__(self, chunks):
        self._chunks = list(chunks)

    def read_ready(self, timeout):
        return self._chunks.pop(0) if self._chunks else ""


class _Engine:
    def __init__(self):
        self.token_usage = {"input": 0, "output": 0}
        self.steered: list = []

    def get_status(self):
        return {}

    def request_stop(self):
        pass

    def steer(self, text, **kw):
        self.steered.append(text)
        return True


def _agent(broker):
    engine = _Engine()
    err = io.StringIO()
    agent = repl.TerminalAgent(engine, out=io.StringIO(), err=err,
                               broker=broker)
    agent.transcript.theme = rr.Theme(enabled=False)
    return agent, engine, err


def test_a_line_typed_at_the_terminal_dialog_is_the_answer(rooms):
    broker = tc.TerminalConfirmBroker(timeout_s=5.0)
    agent, engine, err = _agent(broker)
    req = broker._enqueue(tc.ConfirmRequest(
        kind=tc.ASK, tool="ask_user_question", payload=dict(PAYLOAD)))
    # Typed in pieces, with a digit inside the text and one correction.
    agent._answer_question(req, _Keys(["", "Ein Dashboard-Panel", " v2x",
                                       "\x7f", " für Konformere", "\r"]))
    assert req.resolved is True
    assert req.decision["answers"] == ["Ein Dashboard-Panel v2 für Konformere"]
    assert req.decision["answered_in_chat"] is True
    assert broker.take() is None
    assert engine.steered == [] and agent.queued == [], (
        "the answer must not also go out as a message")


def test_the_terminal_dialog_shows_labels_not_dicts(rooms):
    broker = tc.TerminalConfirmBroker(timeout_s=5.0)
    agent, engine, err = _agent(broker)
    req = broker._enqueue(tc.ConfirmRequest(
        kind=tc.ASK, tool="ask_user_question", payload=dict(PAYLOAD)))
    agent._answer_question(req, _Keys(["3"]))
    assert req.decision == {"answers": ["Both"]}
    screen = err.getvalue()
    assert "3. Both" in screen and "{'label'" not in screen
    assert "type" in screen.lower(), "the dialog says an answer can be typed"


# ---------------------------------------------------------------------------
# Dashboard: the real tab, its widget tree
# ---------------------------------------------------------------------------

def _walk(root):
    stack = [root]
    while stack:
        node = stack.pop()
        yield node
        stack.extend(getattr(node, "children", ()) or ())


def test_a_chat_message_answers_the_dashboard_question(tmp_path, monkeypatch):
    widgets = pytest.importorskip("ipywidgets")
    from delfin.agent import scheduler as S
    from delfin.dashboard import tab_agent
    from delfin.dashboard.context import DashboardContext

    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))
    for name in ("calc", "archive", "office", "ws"):
        (tmp_path / name).mkdir(exist_ok=True)
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script: None
    ctx.agent_workspace = str(tmp_path / "ws")
    tab, refs = tab_agent.create_tab(ctx)
    try:
        state = refs["state"]
        engine = SimpleNamespace(kit_permissions=SimpleNamespace(),
                                 steered=[], typed=[])
        engine.steer = lambda text, **kw: engine.steered.append(text) or True
        engine.note_typed = lambda text: engine.typed.append(text)
        state["engine"] = engine
        state["_wire_phase5_callbacks"](engine)
        ask = engine.kit_permissions.ask_user_callback
        # The turn is running: it is blocked inside the question.
        state["streaming"] = True

        box = next(w for w in _walk(tab) if isinstance(w, widgets.VBox)
                   and "delfin-agent-request" in w._dom_classes
                   and w.layout.border == "2px solid #0a84ff")
        out: dict = {}
        t = threading.Thread(target=lambda: out.setdefault(
            "r", ask(dict(PAYLOAD))), daemon=True)
        t.start()
        for _ in range(500):
            if box.layout.display == "":
                break
            time.sleep(0.01)
        assert box.layout.display == "", "the question was never shown"
        assert any("Was genau soll gebaut werden" in getattr(w, "value", "")
                   for w in _walk(box))

        inp = next(w for w in _walk(tab) if isinstance(w, widgets.Textarea)
                   and "delfin-agent-input" in w._dom_classes)
        send = next(w for w in _walk(tab) if isinstance(w, widgets.Button)
                    and w.description == "Send")
        inp.value = TYPED
        send.click()
        # Closed by the click itself, before the asking thread wakes.
        assert box.layout.display == "none"
        t.join(5)
        assert not t.is_alive(), "the turn still waits on the question"
        assert out["r"]["answers"] == [TYPED]
        assert out["r"]["answered_in_chat"] is True
        assert box.layout.display == "none"
        assert not any(isinstance(w, widgets.Button) for w in _walk(box)), (
            "the option buttons are still in the tree")
        assert state.get("_ask_user_event") is None
        assert inp.value == ""
        assert engine.steered == [], "the answer also went into the turn"
        assert list(state["message_queue"]) == [], "the answer was queued"
        users = [m for m in state["chat_messages"] if m.get("role") == "user"]
        assert [m["content"] for m in users] == [TYPED], (
            "the typed answer appears in the chat exactly once")
    finally:
        refs["shutdown"]()
