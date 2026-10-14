"""A reload takes its conversation over from the kernel it replaces.

Input: a dashboard kernel holding a session's writer lock, whose window has
gone (a reload: the server keeps the old kernel for the grace so a dropped
connection can come back), and a new kernel opening the same session.

Output: the new kernel opens the conversation and holds the lock; the old
kernel's next save is refused. Before this, the new kernel was refused with
"held by a live session", the chat came back empty, and the session list
was saved without the conversation, so the next start did not reopen it.

Semantics: the lock is handed over only when its holder is a kernel on this
host whose window the server has recorded as gone and which runs no turn.
A holder with a window, or one on another host, still refuses: two live
views would save over each other's turns.

The holder is a real second process, as in the field: a lock held by this
test's own pid would read as its own and prove nothing.
"""
from __future__ import annotations

import asyncio
import json
import subprocess
import sys
import time

import pytest

pytest.importorskip("ipywidgets")

from delfin.agent import session_store as ss                    # noqa: E402
from delfin.dashboard import resume_server as R                  # noqa: E402

from tests.conftest import child_env                             # noqa: E402
from tests.test_the_server_ends_what_nobody_watches import _Base  # noqa: E402

SID = "b1c2d3e4f5a6b7c8d9e0f1a2b3c4d5e6"
OLD_KERNEL = "0ld0ld00-0000-4000-8000-000000000001"

_HOLDER = r"""
import sys
from delfin.agent import session_store as ss
from pathlib import Path
ss._SESSIONS_DIR = Path(sys.argv[1])
_identify = getattr(ss, "set_lock_identity", None)
if _identify:
    _identify(kernel=sys.argv[3])
ss.acquire_session_lock(sys.argv[2])
print("HELD", flush=True)
for line in sys.stdin:
    if line.strip() == "save":
        try:
            ss.acquire_session_lock(sys.argv[2])
            print("SAVED", flush=True)
        except ss.SessionLockedError as exc:
            print(f"REFUSED {exc.holder_pid}", flush=True)
"""


@pytest.fixture
def holder(tmp_path, monkeypatch):
    """A second process holding SID's lock as kernel OLD_KERNEL."""
    store = tmp_path / "child_home" / ".delfin" / "agent_sessions"
    store.mkdir(parents=True)
    monkeypatch.setattr(ss, "_SESSIONS_DIR", store)
    history = []
    for i in range(6):
        history += [{"role": "user", "content": f"task {i}"},
                    {"role": "assistant", "content": f"answer {i}"}]
    ss.save_session(SID, mode="solo", engine_messages=history,
                    chat_messages=[dict(m, role_label="") for m in history],
                    workspace=str(tmp_path / "ws"))
    ss.release_session_lock(SID)
    import os
    env = child_env(tmp_path)
    env["PYTHONPATH"] = os.pathsep.join(p for p in sys.path if p)
    proc = subprocess.Popen(
        [sys.executable, "-c", _HOLDER, str(store), SID, OLD_KERNEL],
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, env=env)
    try:
        assert proc.stdout.readline().strip() == "HELD"
        yield proc, history
    finally:
        proc.kill()
        proc.wait(10)


def _ask_holder_to_save(proc) -> str:
    proc.stdin.write("save\n")
    proc.stdin.flush()
    return proc.stdout.readline().strip()


def _server(windows_left: int):
    """The real server-side manager over a stand-in base: the old kernel
    had one window and now has ``windows_left``."""
    manager = R.resume_kernel_manager_class(_Base)()
    asyncio.run(manager.start_kernel(kernel_id=OLD_KERNEL))
    manager.notify_connect(OLD_KERNEL)
    if windows_left == 0:
        manager.notify_disconnect(OLD_KERNEL)
    return manager


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
    ctx.initial_session_id = SID
    ctx.sessions_managed = True
    ctx.presence_key = "reload"
    _tab, refs = tab_agent.create_tab(ctx)
    return refs


def _chat_text(state) -> str:
    return "\n".join(str(m.get("content") or "")
                     for m in state.get("chat_messages") or [])


def test_a_reload_opens_the_conversation_its_windowless_kernel_held(
        tmp_path, monkeypatch, holder):
    proc, history = holder
    _server(windows_left=0)
    refs = _open_tab(tmp_path, monkeypatch)
    try:
        state = refs["state"]
        assert "held by a live session" not in _chat_text(state)
        assert state.get("active_session_id") == SID
        assert len(state["engine"].messages) == len(history)
        assert ss.session_lock_holder(SID) is None, "the lock is ours now"
        ss.acquire_session_lock(SID)          # a save here is not refused
        assert _ask_holder_to_save(proc) == f"REFUSED {__import__('os').getpid()}"
    finally:
        refs["shutdown"](save=False)


def test_a_holder_that_still_has_its_window_keeps_the_conversation(
        tmp_path, monkeypatch, holder):
    proc, _history = holder
    _server(windows_left=1)
    refs = _open_tab(tmp_path, monkeypatch)
    try:
        state = refs["state"]
        assert "held by a live session" in _chat_text(state)
        assert state.get("active_session_id") != SID
        assert _ask_holder_to_save(proc) == "SAVED"
    finally:
        refs["shutdown"](save=False)


def test_a_conversation_refused_on_start_stays_in_the_open_list(
        tmp_path, monkeypatch):
    """The session list remembers what was open. A conversation the new
    page could not open because a live view holds it is still open there,
    and dropping it from the list meant the next start forgot it."""
    import ipywidgets as widgets
    from delfin.agent import scheduler as S
    from delfin.dashboard import agent_sessions as AS
    from delfin.dashboard.context import DashboardContext

    store = tmp_path / "sessions"
    store.mkdir()
    monkeypatch.setattr(ss, "_SESSIONS_DIR", store)
    monkeypatch.setattr(AS, "_OPEN_SESSIONS_PATH", tmp_path / "open.json")
    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))
    monkeypatch.delenv("DELFIN_RESUME_SESSION", raising=False)
    for name in ("launch", "calc", "agent_workspace", "archive", "office"):
        (tmp_path / name).mkdir()
    monkeypatch.setenv("DELFIN_LAUNCH_CWD", str(tmp_path / "launch"))
    monkeypatch.setattr(AS, "_saved_session",
                        lambda sid: {"session_id": sid} if sid == SID else None)
    # pid 1 is alive on every POSIX host and is not this process.
    (store / f"{SID}.lock").write_text(json.dumps(
        {"pid": 1, "ts": time.time(), "host": ss._this_host()}))
    AS.save_open_sessions([{"session_id": SID, "workspace": ""}])

    def _build(ctx):
        state = {"active_session_id": "", "chat_messages": [],
                 "streaming": False, "engine": None}
        ctx.agent_state = state
        return widgets.VBox(), {"state": state, "shutdown": lambda: None,
                                "workspace": lambda: ctx.agent_workspace}

    ctx = DashboardContext(calc_dir=tmp_path / "calc",
                           agent_dir=tmp_path / "agent_workspace",
                           archive_dir=tmp_path / "archive",
                           office_dir=tmp_path / "office")
    ctx.run_js = lambda script: None
    _widget, refs = AS.create_tab(ctx, build=_build)
    refs["refresh"]()
    assert [r["session_id"] for r in AS.load_open_sessions()] == [SID]
