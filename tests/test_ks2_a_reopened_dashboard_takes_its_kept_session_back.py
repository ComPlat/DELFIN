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

    def __call__(self, ctx):
        state = {"active_session_id": ctx.initial_session_id or "",
                 "chat_messages": [], "streaming": False, "engine": None}
        ctx.agent_state = state
        self.contexts.append(ctx)
        stop = widgets.Button(description="Stop")
        row = widgets.HBox([stop])
        return widgets.VBox([row]), {"state": state, "shutdown": lambda: None,
                                     "control_row": row, "stop_button": stop}

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


def test_picked_from_the_list_while_held_it_is_handed_over_not_viewed(home, old_kernel):
    """Not among the remembered open sessions (another dashboard rewrote
    that list; seen 2026-11-16): picked from the session list while the
    closed window's kernel still hosts it, it is handed over as well."""
    K.arm(SID, str(home / "launch"), pid=old_kernel.pid, title="kept job")
    build = _Build()
    _w, refs = AS.create_tab(_ctx(home), build=build, service=False)
    try:
        refs["running"].summary = lambda sid: ({"name": "kept job", "host": "dashboard",
                                                "state": "idle", "keep": True,
                                                "workspace": str(home / "launch")}
                                               if sid == SID else None)
        rec = refs["open_from_list"](SID)
        assert rec.get("handing_over") and not rec.get("live")
        assert refs["active"]() is rec
        _end(old_kernel)
        refs["refresh"]()
        mine = [r for r in refs["sessions"]() if r["ctx"].initial_session_id == SID]
        assert mine and not mine[0].get("dormant") and not mine[0].get("live")
        assert refs["active"]() is mine[0] and mine[0]["kept"]
    finally:
        refs["stop_refresh"]()


def test_a_kept_session_working_there_is_still_followed_live(home, old_kernel):
    K.arm(SID, str(home / "launch"), pid=old_kernel.pid)
    _w, refs = AS.create_tab(_ctx(home), build=_Build(), service=False)
    try:
        refs["running"].summary = lambda sid: ({"host": "dashboard", "state": "working"}
                                               if sid == SID else None)
        rec = refs["open_from_list"](SID)
        assert rec is None or not rec.get("handing_over")
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
