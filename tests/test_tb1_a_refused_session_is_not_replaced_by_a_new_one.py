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
    import os
    from tests.conftest import child_env
    env = child_env(tmp_path)
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
    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))
    monkeypatch.setenv("KIT_TOOLBOX_API_KEY", "sk-test")
    path = user_settings.get_settings_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"agent": {
        "provider": "kit", "model": "kit.deepseek-v4-flash",
        "permission_profile": "all_free", "probe_endpoint": False}}))
    for name in ("calc", "archive", "office", "ws"):
        (tmp_path / name).mkdir(exist_ok=True)
    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script: None
    ctx.agent_workspace = str(tmp_path / "ws")
    ctx.sessions_managed = True
    _tab, refs = tab_agent.create_tab(ctx)
    return refs


def _writable_sessions():
    return sorted(r["session_id"] for r in ss.list_sessions()
                  if not r.get("deleted"))


def test_a_send_after_the_refusal_starts_no_session_and_keeps_the_text(
        tmp_path, monkeypatch, store, holder):
    """THE red test. Today: the refusal stores nothing, the send mints a new
    empty session and the typed text is wiped."""
    refs = _open_tab(tmp_path, monkeypatch)
    try:
        state = refs["state"]
        assert state.get("active_session_id") != SID
        assert _writable_sessions() == [SID]
        # The person opens the held conversation; the open is refused.
        refs["load_session"](SID)
        assert "held by a live session" in _chat_text(state)
        # The refusal remembered WHAT it refused: the session id is known.
        assert state.get("_held_refused", {}).get("session_id") == SID
        # The person types their next line into the box.
        refs["input"].value = "please continue the analysis"
        # ... and sends it (the real path a line takes, as from a host too).
        refs["host_adapter"]["send"]("please continue the analysis")
        # Nothing new in the store: the send did NOT start a session.
        assert _writable_sessions() == [SID], (
            "a send while refused must not start a new session")
        # One line names both choices.
        chat = _chat_text(state)
        assert "Start a new session from it" in chat
        assert "Wait" in chat
        # The typed text was never lost and never sent anywhere.
        assert refs["input"].value == "please continue the analysis"
        assert not state.get("_last_user_message")
    finally:
        refs["shutdown"](save=False)


def test_choosing_start_makes_exactly_one_session_from_the_held_one(
        tmp_path, monkeypatch, store, holder):
    state = refs = None
    refs = _open_tab(tmp_path, monkeypatch)
    try:
        state = refs["state"]
        refs["load_session"](SID)
        assert state.get("_held_refused", {}).get("session_id") == SID
        n_before = len(_writable_sessions())
        refs["host_adapter"]["send"](
            "Start a new session from it")
        sessions = _writable_sessions()
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
    refs = _open_tab(tmp_path, monkeypatch)
    try:
        state = refs["state"]
        refs["load_session"](SID)
        assert state.get("_held_refused", {}).get("session_id") == SID
        assert state.get("active_session_id") != SID
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
        refs["host_adapter"]["send"]("now it is mine")
        assert _writable_sessions() == [SID]
    finally:
        refs["shutdown"](save=False)


def test_a_new_session_choice_clears_the_refused_state(
        tmp_path, monkeypatch, store, holder):
    refs = _open_tab(tmp_path, monkeypatch)
    try:
        state = refs["state"]
        refs["load_session"](SID)
        assert state.get("_held_refused", {}).get("session_id") == SID
        refs["load_session"]("")           # "New session" selected
        assert state.get("_held_refused") is None, (
            "choosing a different session clears the refused state")
        assert _writable_sessions() == [SID]
    finally:
        refs["shutdown"](save=False)
