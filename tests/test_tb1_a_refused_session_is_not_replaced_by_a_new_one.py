"""TB1: a tab that was refused a session does not quietly start a new one.

Seen 2026-11-17: a kept session is held by another live process; opening it
in a tab shows "This conversation is held by a live session ... start a new
session from it." and returns. The tab holds no session, so the NEXT line
the person types starts a brand-new, empty session — silently. The person
believes they write into the held conversation; the text lands in a stranger
session and the held one never sees it.

Target behaviour, while the refused state stands:
  - a send does NOT start a session (nothing new in the session store),
  - the text stays in the input box (never lost, never sent anywhere),
  - one line shows the two choices: "Start a new session from it" and "Wait".
  - choosing a session in the list, "New session", or the holder ending
    clears the refused state.

The widgets are the real ones (tab_agent.create_tab, as
tests/test_a_reload_takes_over_from_its_dying_kernel.py drives them); the
holder is a real second process, as in the field.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time

import pytest

pytest.importorskip("ipywidgets")

from delfin.agent import session_store as ss                       # noqa: E402

SID = "b1c2d3e4f5a6b7c8d9e0f1a2b3c4d5e6"

_HOLD = r"""
import sys, time
from pathlib import Path
from delfin.agent import session_store as ss
ss._SESSIONS_DIR = Path(sys.argv[1])
ss.set_lock_identity(kernel="holder-proc")
ss.acquire_session_lock(sys.argv[2])
print("HELD", flush=True)
for line in sys.stdin:
    if line.strip() == "release":
        ss.release_session_lock(sys.argv[2])
        print("RELEASED", flush=True)
"""


def _chat_text(state) -> str:
    return "\n".join(str(m.get("content") or "")
                     for m in state.get("chat_messages") or [])


@pytest.fixture
def store(tmp_path, monkeypatch):
    path = tmp_path / "sessions"
    path.mkdir()
    monkeypatch.setattr(ss, "_SESSIONS_DIR", path)
    # A saved conversation with history, so a fork/new-session-from-it has
    # something real to start from.
    history = []
    for i in range(4):
        history += [{"role": "user", "content": f"task {i}"},
                    {"role": "assistant", "content": f"answer {i}"}]
    ss.save_session(SID, mode="solo", engine_messages=history,
                    chat_messages=[dict(m, role_label="") for m in history],
                    workspace=str(tmp_path / "ws"))
    ss.release_session_lock(SID)
    return path


@pytest.fixture
def holder(store, tmp_path):
    """A second live process holding SID's writer lock."""
    from delfin.agent import session_store as _ss_child
    env = dict(os.environ)
    # A private HOME for the child too (wave 33: tests redirect every
    # ~/.delfin sink; the child resolves its own at import time).
    env["HOME"] = str(tmp_path / "child_home")
    (tmp_path / "child_home").mkdir(exist_ok=True)
    env.pop("JUPYTER_RUNTIME_DIR", None)
    env["PYTHON_KEYRING_BACKEND"] = "keyring.backends.null.Keyring"
    env["PYTHONPATH"] = os.pathsep.join(p for p in sys.path if p)
    proc = subprocess.Popen(
        [sys.executable, "-c", _HOLD, str(store), SID],
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, env=env)
    try:
        assert proc.stdout.readline().strip() == "HELD"
        yield proc
    finally:
        proc.stdin.write("release\n")
        proc.stdin.flush()
        try:
            proc.stdout.readline()
        except Exception:
            pass
        if proc.poll() is None:
            proc.kill()
            proc.wait(10)


def _open_tab(tmp_path, monkeypatch):
    from delfin import user_settings
    from delfin.agent import scheduler as S
    from delfin.dashboard import tab_agent
    from delfin.dashboard.context import DashboardContext
    # Redirected HOME first: _open_tab writes a settings file next to the
    # real home otherwise (user_settings.py:518-519); a developer's file
    # must never be touched by a test (reviewer finding F3).
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))
    path = user_settings.get_settings_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"agent": {
        "provider": "ollama", "model": "fake-model",
        "permission_profile": "all_free", "probe_endpoint": False}}))
    for name in ("calc", "archive", "office", "ws"):
        (tmp_path / name).mkdir(exist_ok=True)
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script: None
    ctx.agent_workspace = str(tmp_path / "ws")
    ctx.sessions_managed = True
    # OFFLINE: no real engine is ever built. tab_agent imports AgentEngine
    # inside _ensure_engine (tab_agent.py:12533), so the module attribute
    # is what the patch replaces — a stub with the keys the tab reads and
    # a client with kill. A real kit engine would probe ki-toolbox
    # (operator refusal 2026-11; _probe_endpoint_in_background :12850).
    created_engines = []
    session_state = {}
    saved_states = {}

    class _StubClient:
        def __init__(self):
            self.calls = []

        def kill(self):
            self.calls.append("kill")

    class _StubEngine:
        def __init__(self, sid=""):
            self.session_id = sid
            self.client = _StubClient()
            self.messages = []
            self.role_outputs = {}
            self.route = ["solo_agent"]
            self.current_role = ""
            self.mode = "solo"
            self.backend = "stub"
            self.created = created_engines
            self.is_cycle_complete = True
            self.current_role_index = 0
            self.request_stop = lambda: None
            self.reset_cycle = lambda mode=None: None

        def get_status(self):
            return {"mode": self.mode, "backend": self.backend, "role": "",
                    "role_index": 0, "role_total": 1, "input_tokens": 0,
                    "output_tokens": 0, "cost_usd": 0.0}

        def restore_state(self, data):
            session_state.update(data if isinstance(data, dict) else {})
            return type("R", (), {"restored": 1, "failed": 0})()

        def export_state(self):
            saved_states.setdefault("state", self.messages)
            return {}

    def _factory(*args, **kwargs):
        engine = _StubEngine(kwargs.get("session_id", "") or "")
        created_engines.append(engine)
        return engine

    from delfin.agent import engine as engine_mod
    monkeypatch.setattr(engine_mod, "AgentEngine", _factory)
    tab, refs = tab_agent.create_tab(ctx)
    refs["_stub"] = {"engines": created_engines, "state": session_state,
                     "exported": saved_states}
    return tab, refs


def _find_send_button(widget):
    """The REAL chat Send button: the one with the delfin-agent-send class
    (delfin/dashboard/tab_agent.py:10024). Two buttons are described
    'Send'; the hidden bug-report one (tab_agent.py:9045) is built first,
    so the description alone picks the wrong one (reviewer finding F2)."""
    if "delfin-agent-send" in set(getattr(widget, "_dom_classes", ()) or ()):
        return widget
    for child in getattr(widget, "children", ()) or ():
        found = _find_send_button(child)
        if found is not None:
            return found
    return None


def _host(refs) -> dict:
    """The real host adapter dict: refs["host_adapter"] is a FACTORY
    (tab_agent.py:25791), not the dict itself."""
    return refs["host_adapter"]()


def _press_send(tab, refs, text):
    """The gesture the person makes: text in the box, the REAL chat Send
    pressed. _host_send restores the draft itself
    (tab_agent.py:25734-25746), so 'text kept' through it proves nothing
    (reviewer F2)."""
    _host(refs)["input"].value = text
    btn = _find_send_button(tab)
    assert btn is not None, "the real chat Send button was not found"
    btn.click()


def _writable_sessions():
    return sorted(r["session_id"] for r in ss.list_sessions()
                  if not r.get("deleted"))


def test_a_send_after_the_refusal_starts_no_session_and_keeps_the_text(
        tmp_path, monkeypatch, store, holder):
    """THE red test. Today: the refusal stores nothing, the send mints a new
    empty session and the typed text is wiped."""
    tab, refs = _open_tab(tmp_path, monkeypatch)
    try:
        state = refs["state"]
        assert state.get("active_session_id") != SID
        assert _writable_sessions() == [SID]
        # The person opens the held conversation; the open is refused.
        refs["load_session"](SID)
        assert "held by a live session" in _chat_text(state)
        # The person types their next line into the box and presses the
        # REAL chat Send button -- the gesture the person makes.
        # THE DEFECT, observable on main (first, so the red line names the
        # behaviour, not a missing key -- reviewer F1): the send mints a
        # fresh engine + a brand-new empty session, and the text goes to a
        # stranger conversation.
        _press_send(tab, refs, "please continue the analysis")
        assert len(refs["_stub"]["engines"]) == 0, (
            "a send while refused must not mint a new engine: "
            "the text went to a stranger")
        assert _writable_sessions() == [SID], (
            "a send while refused must not start a new session: "
            "the text went to a stranger")
        # One line names both choices.
        chat = _chat_text(state)
        assert "Start a new session from it" in chat
        assert "Wait" in chat
        # The typed text was never lost and never sent anywhere.
        assert _host(refs)["input"].value == "please continue the analysis"
        assert not state.get("_last_user_message")
        # And the refusal remembered WHAT it refused (the key the patch
        # sets; asserted LAST so the red line above names the defect, not
        # the missing variable -- reviewer finding F1).
        assert state.get("_held_refused", {}).get("session_id") == SID
    finally:
        refs["shutdown"](save=False)


def test_choosing_start_makes_exactly_one_session_from_the_held_one(
        tmp_path, monkeypatch, store, holder):
    tab, refs = _open_tab(tmp_path, monkeypatch)
    try:
        state = refs["state"]
        refs["load_session"](SID)
        # The refusal precondition must actually be reached: on main the
        # refusal stores nothing, so the observable fact is the refusal
        # line in the chat.
        assert "held by a live session" in _chat_text(state)
        n_before = len(_writable_sessions())
        _press_send(tab, refs, "Start a new session from it")
        sessions = _writable_sessions()
        # On main the choice line does not exist and the send mints a
        # stranger session with empty history; the first assertion below
        # (exactly one new session, history from the held one) is what
        # distinguishes the feature from today's silent defect.
        assert len(sessions) == n_before + 1, (
            "exactly one new session, the held one untouched")
        assert SID in sessions
        # The new session starts from the held one's history.
        new_sid = [s for s in sessions if s != SID][0]
        data = ss.load_session(new_sid)
        history = [m.get("content") for m in data.get("engine_messages", [])]
        assert "task 3" in history and "answer 3" in history
        assert state.get("_held_refused") is None, (
            "the choice cleared the refused state")
    finally:
        refs["shutdown"](save=False)


def test_wait_then_the_holder_ends_and_the_tab_opens_the_conversation(
        tmp_path, monkeypatch, store, holder):
    tab, refs = _open_tab(tmp_path, monkeypatch)
    try:
        state = refs["state"]
        refs["load_session"](SID)
        assert "held by a live session" in _chat_text(state)
        assert state.get("active_session_id") != SID
        # The person answers "Wait": the state stays, the tab watches the
        # holder on its own (the Wait branch arms the refused-tab poller,
        # its own bounded timer — no engine exists to arm a background
        # tick, reviewer F6).
        _press_send(tab, refs, "Wait")
        assert "Waiting" in _chat_text(state), (
            "the Wait reply is acknowledged and the text stays")
        # The holder ends; nothing else changes. The tab notices on its own.
        holder.stdin.write("release\n")
        holder.stdin.flush()
        holder.stdout.readline()
        deadline = time.time() + 30
        while time.time() < deadline:
            if state.get("active_session_id") == SID:
                break
            try:
                refs["refresh_chat"]()
            except Exception:
                pass
            time.sleep(0.2)
        assert state.get("active_session_id") == SID, (
            "when the holder ends, the tab opens the held conversation")
        assert state.get("_held_refused") is None
        assert _writable_sessions() == [SID]
        # And now a send goes to THAT conversation, not to a stranger.
        _press_send(tab, refs, "now it is mine")
        assert _writable_sessions() == [SID]
    finally:
        refs["shutdown"](save=False)


def test_a_new_session_choice_clears_the_refused_state(
        tmp_path, monkeypatch, store, holder):
    tab, refs = _open_tab(tmp_path, monkeypatch)
    try:
        state = refs["state"]
        refs["load_session"](SID)
        assert "held by a live session" in _chat_text(state)
        # The real "New session" gesture: "(+ New Session", "") selected
        # in the Session dropdown, the Load button pressed
        # (_on_load_session, tab_agent.py:25125-25131 -> _on_new_cycle).
        # Calling _load_saved_session("") would bypass it.
        session_dd = None

        def _find_dropdown(widget):
            nonlocal session_dd
            if (session_dd is None
                    and isinstance(widget, __import__(
                        "ipywidgets").Dropdown)
                    and widget.description == "Session:"):
                session_dd = widget
            for child in getattr(widget, "children", ()) or ():
                _find_dropdown(child)

        _find_dropdown(tab)
        assert session_dd is not None, "the Session dropdown was not found"
        session_dd.value = ""
        load_btn = None

        def _find_load(widget):
            nonlocal load_btn
            if (load_btn is None
                    and getattr(widget, "description", "") == "Load"
                    and widget.tooltip == "Load the selected session"):
                load_btn = widget
            for child in getattr(widget, "children", ()) or ():
                _find_load(child)

        _find_load(tab)
        assert load_btn is not None, "the Load button was not found"
        load_btn.click()
        assert state.get("_held_refused") is None, (
            "choosing a different session clears the refused state")
        assert _writable_sessions() == [SID]
    finally:
        refs["shutdown"](save=False)


def test_an_ordinary_send_with_no_refusal_still_starts_its_session(
        tmp_path, monkeypatch, store):
    """The green direction (QS condition 4, second half): with NO refused
    state, a send must behave exactly as today — exactly one session
    started, the guard must not block the normal path."""
    tab, refs = _open_tab(tmp_path, monkeypatch)
    try:
        state = refs["state"]
        assert state.get("_held_refused") is None
        n_before = len(_writable_sessions())
        _press_send(tab, refs, "an ordinary first message")
        assert len(refs["_stub"]["engines"]) >= 1, (
            "an ordinary send with no refusal must mint an engine, "
            f"minted {len(refs['_stub']['engines'])}")
        assert len(_writable_sessions()) >= n_before + 1, (
            "an ordinary send with no refusal must start a session")
    finally:
        refs["shutdown"](save=False)
