"""A dashboard opened again takes its kept session back into the full chat.

Input: a session with Keep session on, open when its dashboard closed, and
the dashboard opened again while the closed window's kernel still holds
the session (the window grace) -- or after that kernel ended.

Before (reported 2026-11-16): the new dashboard showed the agent
service's view first (the old kernel hosted the session there), then,
once the old kernel had ended, the hand-over resumed the session in the
service, and only the next take-over stopped it there and opened the full
chat -- the person saw the wrong view and had to wait.
Now the row says the session is being handed over, the new dashboard
waits for the old kernel, and takes the session from its record straight
into its own full chat, kept again; the hand-over leaves a record alone
while a dashboard waits for it. A record nobody waits for any more, a
stop, a closed row: as before.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time

import pytest

pytest.importorskip("ipywidgets")

import ipywidgets as widgets                                     # noqa: E402

from delfin.dashboard import agent_sessions as AS                # noqa: E402
from delfin.dashboard import kept_agent_sessions as K            # noqa: E402

SID = "ks2a000000000001"


class _Build:
    """A stand-in agent tab that records what it was asked to open."""

    def __init__(self):
        self.contexts = []
        self.inputs = []

    def __call__(self, ctx):
        state = {"active_session_id": ctx.initial_session_id or "",
                 "chat_messages": [], "streaming": False, "engine": None}
        ctx.agent_state = state
        self.contexts.append(ctx)
        stop = widgets.Button(description="Stop")
        row = widgets.HBox([stop])
        box = widgets.Textarea()
        self.inputs.append(box)
        return widgets.VBox([row]), {"state": state, "shutdown": lambda: None,
                                     "control_row": row, "stop_button": stop,
                                     "host_adapter": lambda: {"input": box}}

    def opened(self):
        return [c.initial_session_id for c in self.contexts]


@pytest.fixture
def home(tmp_path, monkeypatch):
    from delfin.agent import scheduler as S

    monkeypatch.setattr(AS, "_OPEN_SESSIONS_PATH", tmp_path / "open.json")
    monkeypatch.setattr(S, "_GLOBAL", S.Scheduler(path=tmp_path / "cron.json"))
    monkeypatch.delenv("DELFIN_RESUME_SESSION", raising=False)
    monkeypatch.setenv(K.DIR_ENV, str(tmp_path / "kept"))
    for name in ("launch", "calc", "agent_workspace", "archive", "office"):
        (tmp_path / name).mkdir(exist_ok=True)
    monkeypatch.setenv("DELFIN_LAUNCH_CWD", str(tmp_path / "launch"))
    monkeypatch.setattr(AS, "_saved_session",
                        lambda sid: {"session_id": sid, "title": "kept job"})
    monkeypatch.setattr(AS, "keep_session_refusal", lambda: "")
    # Never the person's real stop file.
    from delfin.agent import stop_all
    monkeypatch.setattr(stop_all, "wakes_allowed", lambda since: True)
    return tmp_path


@pytest.fixture
def old_kernel():
    """The closed window's kernel: a process of its own that holds the record."""
    proc = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(120)"])
    try:
        yield proc
    finally:
        proc.kill()
        proc.wait(timeout=30)


def _ctx(tmp_path):
    from delfin.dashboard.context import DashboardContext

    ctx = DashboardContext(calc_dir=tmp_path / "calc", agent_dir=tmp_path / "agent_workspace",
                           archive_dir=tmp_path / "archive", office_dir=tmp_path / "office")
    ctx.run_js = lambda script: None
    return ctx


def _kept_by(home, proc):
    K.arm(SID, str(home / "launch"), {"model": "m"}, pid=proc.pid, title="kept job")
    AS.save_open_sessions([{"session_id": SID, "workspace": str(home / "launch"),
                            "front": True}])


def _end(proc):
    proc.kill()
    proc.wait(timeout=30)


def _record():
    return json.loads(K._path(SID).read_text())


# -- where a kept session is, and who takes it --------------------------------------

def test_where_kept_tells_held_free_moving_and_nothing(home, old_kernel):
    assert K.where_kept(SID) == ""
    K.arm(SID, "/w", pid=old_kernel.pid)
    assert K.where_kept(SID) == "held"
    _end(old_kernel)
    assert K.where_kept(SID) == "free"
    claimed = K._path(SID).with_name(f".{K._path(SID).name}.claimed-{os.getpid()}")
    os.rename(K._path(SID), claimed)
    assert K.where_kept(SID) == "moving"
    claimed.rename(claimed.with_name(f".{K._path(SID).name}.claimed-999999999"))
    assert K.where_kept(SID) == "", "a claim whose process is gone is not a move"
    K.arm(SID, "/w")                                  # held by this process
    assert K.where_kept(SID) == ""


def test_the_hand_over_leaves_a_wanted_record_to_its_dashboard(home, old_kernel):
    K.arm(SID, "/w", pid=old_kernel.pid)
    _end(old_kernel)
    K.want(SID)
    rows = K.hand_over(connect=lambda: pytest.fail("not to the service"))
    assert rows == [{"session_id": SID, "result": "wanted: a dashboard takes it back"}]
    assert K.is_armed(SID)
    K.unwant(SID)
    assert not K.wanted_here(SID)


def test_a_wanted_marker_of_a_process_gone_is_not_wanted(home, old_kernel):
    K.arm(SID, "/w", pid=old_kernel.pid)
    K._write(K._wanted_path(SID), {"session_id": SID, "pid": old_kernel.pid,
                                   "start_ticks": K._start_ticks(old_kernel.pid),
                                   "host": K._host()})
    assert K.wanted_here(SID)
    _end(old_kernel)
    assert not K.wanted_here(SID)


def test_claim_here_takes_only_a_free_record_once(home, old_kernel):
    K.arm(SID, "/w", {"model": "m"}, pid=old_kernel.pid)
    assert K.claim_here(SID) == (None, "")            # its kernel still runs
    _end(old_kernel)
    record, why = K.claim_here(SID)
    assert why == "" and record["session_id"] == SID and record["settings"] == {"model": "m"}
    assert not K.is_armed(SID)
    assert K.claim_here(SID) == (None, "")            # nothing left to take
    K.put_back(record)
    assert K.is_armed(SID) and _record()["pid"] == old_kernel.pid


def test_a_stop_after_keeping_drops_the_claim(home, old_kernel, monkeypatch):
    from delfin.agent import stop_all
    K.arm(SID, "/w", pid=old_kernel.pid)
    _end(old_kernel)
    monkeypatch.setattr(stop_all, "wakes_allowed", lambda since: False)
    record, why = K.claim_here(SID)
    assert record is None and "emergency stop" in why
    assert not K.is_armed(SID)


def test_a_kernel_that_ends_while_waiting_hands_over_again(home, monkeypatch):
    spawned = []
    monkeypatch.setattr(K, "spawn_hand_over", lambda **kw: spawned.append(kw) or "p")
    K.want(SID)
    assert K.on_kernel_end() == "p"
    assert spawned == [{"wait_pid": os.getpid()}]
    assert not K._wanted_path(SID).exists()


# -- the dashboard opened again ------------------------------------------------------

def test_reopened_in_the_grace_it_says_handed_over_then_opens_the_full_chat(home, old_kernel):
    _kept_by(home, old_kernel)
    build = _Build()
    _w, refs = AS.create_tab(_ctx(home), build=build, service=False)
    try:
        (row,) = refs["sessions"]()
        assert row.get("handing_over") and refs["active"]() is row
        assert build.opened() == [], "nothing resumed while the old kernel holds it"
        text = row["tab"].children[0].value
        assert "handed over" in text and "kept job" in text
        assert "delfin-state-running" in row["row"]._dom_classes
        assert K.wanted_here(SID)
        refs["refresh"]()
        assert refs["sessions"]()[0] is row             # still waiting
        _end(old_kernel)
        refs["refresh"]()
        (rec,) = refs["sessions"]()
        assert not rec.get("dormant") and not rec.get("live")
        assert build.opened() == [SID]
        assert refs["active"]() is rec
        assert rec["kept"] and rec["keep_switch"].value is True
        assert _record()["pid"] == os.getpid(), "kept again, by this kernel"
        assert not K.wanted_here(SID)
    finally:
        refs["stop_refresh"]()


def test_reopened_after_the_old_kernel_ended_it_opens_at_once(home, old_kernel):
    _kept_by(home, old_kernel)
    _end(old_kernel)
    build = _Build()
    _w, refs = AS.create_tab(_ctx(home), build=build, service=False)
    try:
        (rec,) = refs["sessions"]()
        assert not rec.get("dormant") and build.opened() == [SID]
        assert refs["active"]() is rec and rec["kept"]
        assert _record()["pid"] == os.getpid()
    finally:
        refs["stop_refresh"]()


def test_closing_the_waiting_row_leaves_the_session_kept(home, old_kernel, monkeypatch):
    spawned = []
    monkeypatch.setattr(K, "spawn_hand_over", lambda **kw: spawned.append(kw))
    _kept_by(home, old_kernel)
    _w, refs = AS.create_tab(_ctx(home), build=_Build(), service=False)
    try:
        (row,) = refs["sessions"]()
        refs["close"](row["key"])
        assert K.is_armed(SID) and _record()["pid"] == old_kernel.pid
        assert not K.wanted_here(SID)
        assert spawned, "the hand-over runs again for it"
    finally:
        refs["stop_refresh"]()


def test_a_session_no_longer_kept_reads_not_running(home, old_kernel, monkeypatch):
    _kept_by(home, old_kernel)
    build = _Build()
    _w, refs = AS.create_tab(_ctx(home), build=build, service=False)
    try:
        (row,) = refs["sessions"]()
        K.disarm(SID)                                   # switched off elsewhere
        monkeypatch.setattr(AS, "_HAND_OVER_SETTLE_S", 0.0)
        row["settle_until"] = time.monotonic() - 1
        refs["refresh"]()
        # In front, so loaded from its record, as any open session that did
        # not keep working.
        assert build.opened() == [SID]
        assert not [r for r in refs["sessions"]() if r.get("handing_over")]
    finally:
        refs["stop_refresh"]()


def test_a_stop_given_after_keeping_starts_nothing(home, old_kernel, monkeypatch):
    from delfin.agent import stop_all
    _kept_by(home, old_kernel)
    _end(old_kernel)
    monkeypatch.setattr(stop_all, "wakes_allowed", lambda since: False)
    AS.save_open_sessions([{"session_id": SID, "workspace": str(home / "launch")},
                           {"session_id": "ks2b000000000002", "workspace": str(home / "launch"),
                            "front": True}])
    build = _Build()
    _w, refs = AS.create_tab(_ctx(home), build=build, service=False)
    try:
        row = refs["sessions"]()[0]
        assert row.get("dormant") and not row.get("handing_over")
        assert "emergency stop" in row["tab"].children[0].value
        assert SID not in build.opened()
    finally:
        refs["stop_refresh"]()


def test_a_chat_still_loading_its_session_is_one_row_not_two(home):
    """A tab asked to open a conversation has its id only once loaded;
    meanwhile the service's row of the same session was listed beside it."""
    class _Loading(_Build):
        def __call__(self, ctx):
            tab, refs = super().__call__(ctx)
            refs["state"]["active_session_id"] = ""         # not loaded yet
            return tab, refs

    _w, refs = AS.create_tab(_ctx(home), build=_Loading(), service=False)
    try:
        assert not refs["running_panel"].is_open(SID)
        refs["open"](str(home / "launch"), SID, in_kernel=True)
        assert refs["running_panel"].is_open(SID)
        assert not refs["running_panel"].is_open("")
    finally:
        refs["stop_refresh"]()


# -- with the agent service: the closed window's kernel still hosts it ------------------

from fake_host import FakeHost                                   # noqa: E402
from tui_service import ServiceProcess                           # noqa: E402

HSID = "ks2c0ffee0000000000000000000003"
CHAT = ('<template class="delfin-chat-src"><div class="delfin-agent-chat">'
        '<div class="delfin-chat-msg delfin-chat-user">hi</div></div></template>')


def _until(cond, timeout: float = 30.0):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        got = cond()
        if got:
            return got
        time.sleep(0.05)
    return cond()


@pytest.fixture
def svc(tmp_path):
    (tmp_path / "svcroot").mkdir()
    service = ServiceProcess(tmp_path / "svcroot")
    try:
        yield service
    finally:
        service.stop()


@pytest.fixture
def kernel_mode(home, svc, monkeypatch):
    from delfin.dashboard import live_sessions as LS
    dash = LS.LiveClient(svc.socket)
    svc.clients.append(dash)
    monkeypatch.setattr(LS, "connect", lambda **kw: dash)
    monkeypatch.delenv(AS.SESSIONS_ENV, raising=False)
    return dash


def _rows_for(refs, sid):
    mine = [r for r in refs["sessions"]() if r.get("session_id") == sid
            or r["ctx"].initial_session_id == sid]
    listed = [r for r in refs["running_panel"].ordered() if r["sid"] == sid]
    return mine, listed


def _classes_of(w):
    return tuple(getattr(w, "_dom_classes", ()))


def _all(w):
    yield w
    for child in getattr(w, "children", ()):
        yield from _all(child)


def test_reopened_while_held_it_looks_like_the_chat_says_so_takes_input_then_is_the_chat(
        home, svc, kernel_mode, old_kernel):
    host = FakeHost(svc.socket, HSID, str(home / "launch"), keep=True)
    host.emit("dashboard.chat", v=1, html=CHAT)
    host.emit("turn.started", turnId="t1", source="user")      # mid-step there
    K.arm(HSID, str(home / "launch"), pid=old_kernel.pid, title="kept job")
    AS.save_open_sessions([{"session_id": HSID, "workspace": str(home / "launch"),
                            "front": True}])
    build = _Build()
    _w, refs = AS.create_tab(_ctx(home), build=build, service=lambda: kernel_mode,
                              start_in="kernel")
    try:
        assert _until(lambda: refs["running"].summary(HSID))
        refs["refresh"]()
        mine, listed = _rows_for(refs, HSID)
        assert len(mine) == 1 and not listed, "one row for it, never two"
        view = mine[0]
        assert view.get("live") and view.get("into_chat") and view.get("kept_here")
        assert K.wanted_here(HSID)
        # Laid out as the full chat: its root, its composer, its chips.
        classes = {c for w in _all(view["tab"]) for c in _classes_of(w)}
        for cls in ("delfin-agent-root", "delfin-agent-composer", "delfin-agent-input",
                    "delfin-agent-chip-perm", "delfin-agent-send", "delfin-agent-submit"):
            assert cls in classes, cls
        # The header as the full chat's: Mode, Provider, Auto, Stop, then
        # Keep session right after Stop; the beetle; "Details" at the foot.
        row = view["refs"]["controls"].children
        words = [getattr(w, "description", "") for w in row]
        assert words[:2] == ["Mode:", "Provider:"] and "Auto: off" in words
        stop = view["refs"]["cancel_btn"]
        assert "delfin-keep-session" in _classes_of(row[list(row).index(stop) + 1])
        assert "🐞" in words
        assert view["refs"]["details"].get_title(0) == "Details"
        refs["activate"](view["key"])
        line = view["refs"]["transit_html"].value
        assert "Coming over" in line and "type" in line, line
        # Typing works meanwhile: sent to where it runs.
        view["refs"]["submit"].value = '{"seq": 1, "text": "from the reopened window"}'
        assert host.wait_asked("session.prompt")["text"] == "from the reopened window"
        view["refs"]["input"].value = "a draft not sent yet"
        # The closed window's kernel finishes and ends.
        host.emit("turn.ended", turnId="t1", durationMs=5, inputTokens=1, outputTokens=1,
                  stopReason="end_turn")
        host.close()
        _end(old_kernel)
        refs["refresh"]()
        mine, listed = _rows_for(refs, HSID)
        assert len(mine) == 1 and not listed
        chat = mine[0]
        assert not chat.get("live") and not chat.get("dormant")
        assert build.opened() == [HSID] and refs["active"]() is chat
        assert chat["kept"]
        assert json.loads(K._path(HSID).read_text())["pid"] == os.getpid()
        assert build.inputs[-1].value == "a draft not sent yet", "the draft came along"
        assert not K.wanted_here(HSID)
    finally:
        refs["stop_refresh"]()
        host.close()


def test_picked_from_the_list_while_held_it_is_the_same_view_and_then_the_chat(
        home, svc, kernel_mode, old_kernel):
    host = FakeHost(svc.socket, HSID, str(home / "launch"), keep=True)
    host.emit("dashboard.chat", v=1, html=CHAT)
    K.arm(HSID, str(home / "launch"), pid=old_kernel.pid, title="kept job")
    build = _Build()
    _w, refs = AS.create_tab(_ctx(home), build=build, service=lambda: kernel_mode,
                              start_in="kernel")
    try:
        assert _until(lambda: refs["running"].summary(HSID))
        view = refs["open_from_list"](HSID)
        assert view.get("live") and view.get("kept_here")
        assert "Coming over" in view["refs"]["transit_html"].value
        mine, listed = _rows_for(refs, HSID)
        assert len(mine) == 1 and not listed
        host.close()
        _end(old_kernel)
        refs["refresh"]()
        mine, listed = _rows_for(refs, HSID)
        assert len(mine) == 1 and not mine[0].get("live") and not listed
        # (the empty first session of a dashboard opened with nothing to restore)
        assert build.opened().count(HSID) == 1 and mine[0]["kept"]
    finally:
        refs["stop_refresh"]()
        host.close()


def test_a_conversation_another_process_writes_is_not_opened_until_it_lets_go(
        home, old_kernel, tmp_path):
    """Seen 2026-11-16: the service woke the session once the closed
    window's kernel was gone, and the dashboard opened it as well -- its
    chat stayed empty ("held by a live session"). Now it waits for the
    writer to let go, keeps the record, then opens it."""
    from delfin.agent import session_store
    _kept_by(home, old_kernel)
    build = _Build()
    _w, refs = AS.create_tab(_ctx(home), build=build, service=False)
    writer = subprocess.Popen([sys.executable, "-c", (
        "import sys, time; from delfin.agent import session_store as S; "
        "assert S.acquire_session_lock(sys.argv[1]); print('held', flush=True); "
        "time.sleep(120)"), SID], stdout=subprocess.PIPE, text=True,
        env=dict(os.environ, PYTHONPATH=os.pathsep.join(sys.path)))
    try:
        assert writer.stdout.readline().strip() == "held"
        assert session_store.session_lock_holder(SID) is not None
        _end(old_kernel)
        refs["refresh"]()
        assert build.opened() == [], "not opened while another process writes it"
        assert K.is_armed(SID), "its record is kept for the next try"
        (row,) = refs["sessions"]()
        assert row.get("handing_over")
        writer.kill()
        writer.wait(timeout=30)
        refs["refresh"]()
        (rec,) = refs["sessions"]()
        assert build.opened() == [SID] and not rec.get("dormant") and rec["kept"]
    finally:
        if writer.poll() is None:
            writer.kill()
            writer.wait(timeout=30)
        refs["stop_refresh"]()
