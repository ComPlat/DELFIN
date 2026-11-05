"""One live session, open in several windows at once, under one name.

The agent service holds every session; the TUI, the dashboard and the
phone bridge are its windows. These tests drive a private service with the
scripted engine and check what all windows must agree on:

* the name: one title stored in the service, trimmed, at most
  ``protocol.TITLE_MAX`` (50, SimpleX's group-name limit) characters, no
  control characters; empty means the default (folder name and short id);
  every summary carries the display name (``name``); a rename in one window
  reaches every other live; it outlives a service restart;
* an approval: the first answer wins, the second window is told, and
  ``permission.resolved`` names the window (``client``) that answered --
  stamped by the service from that window's hello, never taken from what
  the window claims;
* people mode: each person has a service of their own (own HOME and
  DELFIN_RUNTIME_DIR); a window of one person neither lists, attaches to,
  nor renames another person's session.
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path

import pytest

from delfin.agent.service import protocol
from delfin.agent.service.client import ServiceClient, ServiceError
from tui_service import ServiceProcess, turn


@pytest.fixture
def service(tmp_path):
    svc = ServiceProcess(tmp_path)
    try:
        yield svc
    finally:
        svc.stop()


def _summary(c: ServiceClient, sid: str) -> dict:
    return next(s for s in c.request("sessions.list")["sessions"] if s["sessionId"] == sid)


def _until(cond, timeout: float = 20.0):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        got = cond()
        if got:
            return got
        time.sleep(0.05)
    return cond()


def test_the_name_limit_is_simplex_s_group_name_limit():
    assert protocol.TITLE_MAX == 50
    assert protocol.clean_title("  " + "y" * 70 + "  ") == "y" * 50
    assert protocol.clean_title("a\x00b\x1b[31mc​ d") == "a b c d"


def test_the_default_name_is_the_folder_and_a_short_id():
    sid = "0123456789abcdef"
    assert protocol.default_name("/work/scan-a", sid) == "scan-a · 01234567"
    assert protocol.display_name("", "/work/scan-a", sid) == "scan-a · 01234567"
    assert protocol.display_name("my scan", "/work/scan-a", sid) == "my scan"
    assert protocol.display_name(None, "", sid) == "session · 01234567"


def test_every_summary_carries_the_display_name(service):
    a = service.client()
    sid = service.session(a, cwd=str(service.workspace("scan-a")))
    assert _summary(a, sid)["name"] == f"scan-a · {sid[:8]}"
    a.request("session.setTitle", {"sessionId": sid, "title": "solvent scan"})
    assert _until(lambda: _summary(a, sid)["name"] == "solvent scan")


def test_a_rename_is_trimmed_cut_and_cleaned(service):
    a = service.client()
    sid = service.session(a)
    got = a.request("session.setTitle", {"sessionId": sid,
                                         "title": "  \x1b[31mred\x07 " + "n" * 80})
    assert got["title"] == ("red " + "n" * 80)[:50].rstrip()
    assert len(got["title"]) <= 50
    assert _until(lambda: _summary(a, sid)["title"] == got["title"])


def test_an_empty_name_goes_back_to_the_default(service):
    a, b = service.client(), service.client()
    sid = service.session(a, cwd=str(service.workspace("calc")), title="first name")
    b.request("session.attach", {"sessionId": sid})
    assert _summary(b, sid)["name"] == "first name"
    mark = len(b.events)
    assert a.request("session.setTitle", {"sessionId": sid, "title": "  \x1b[0m "}) == {"title": ""}
    b.wait_for(lambda e: e.get("type") == "session.updated" and e.get("session") == sid
               and e.get("title") == "", after=mark)
    summary, _ = b.wait_for(lambda e: e.get("type") == "session.summary"
                            and e.get("session") == sid
                            and e["summary"].get("name") == f"calc · {sid[:8]}", after=mark)
    assert summary["summary"]["title"] is None


def test_two_sessions_may_share_a_name_and_the_short_id_tells_them_apart(service):
    a = service.client()
    one = service.session(a, cwd=str(service.workspace("w1")), title="scan")
    two = service.session(a, cwd=str(service.workspace("w2")), title="scan")
    names = {s["sessionId"]: s["name"] for s in a.request("sessions.list")["sessions"]}
    assert names[one] == names[two] == "scan" and one[:8] != two[:8]


def test_a_rename_in_one_window_reaches_every_other_live(service):
    tui, dash, phone = service.client(), service.client(), service.client()
    sid = service.session(tui)
    for renamer, name in ((tui, "from the terminal"), (dash, "from the dashboard"),
                          (phone, "from the phone")):
        marks = {c: len(c.events) for c in (tui, dash, phone)}
        renamer.request("session.setTitle", {"sessionId": sid, "title": name})
        for c, mark in marks.items():
            c.wait_for(lambda e: e.get("type") == "session.summary" and e.get("session") == sid
                       and e["summary"].get("name") == name, after=mark)


def test_a_name_outlives_a_service_restart(service):
    a = service.client()
    sid = service.session(a, fakeSave=True, cwd=str(service.workspace("kept")))
    turn(a, sid, "hello there")
    a.request("session.setTitle", {"sessionId": sid, "title": "kept name"})
    _until(lambda: _summary(a, sid)["title"] == "kept name")
    a.request("service.shutdown", timeout=30)
    service.proc.wait(timeout=60)
    service.start()
    b = service.client()
    asleep = _summary(b, sid)
    assert (asleep["state"], asleep["title"], asleep["name"]) == ("sleeping", "kept name", "kept name")
    attached = b.request("session.attach", {"sessionId": sid, "turns": 0}, timeout=120)
    assert attached["session"]["title"] == "kept name"
    assert attached["summary"]["name"] == "kept name"


def test_the_first_answer_wins_and_every_window_hears_who_gave_it(service):
    tui = ServiceClient(service.socket, client="delfin-tui")
    dash = ServiceClient(service.socket, client="delfin-dashboard")
    service.clients += [tui, dash]
    sid = service.session(tui, permissionMode="default")
    dash.request("session.attach", {"sessionId": sid})
    marks = {tui: len(tui.events), dash: len(dash.events)}
    tui.request("session.prompt", {"sessionId": sid, "text": "edit it"})
    asked, _ = dash.wait_for("permission.asked", after=marks[dash])
    tui.wait_for("permission.asked", after=marks[tui])
    # A window may not say it is another: the service stamps who answered.
    dash.request("permission.reply", {"sessionId": sid, "requestId": asked["requestId"],
                                      "choice": "y", "client": "delfin-tui"})
    with pytest.raises(ServiceError) as late:
        tui.request("permission.reply", {"sessionId": sid,
                                         "requestId": asked["requestId"], "choice": "n"})
    assert late.value.code == protocol.UNKNOWN_REQUEST
    for c, mark in marks.items():
        resolved, _ = c.wait_for("permission.resolved", after=mark)
        assert (resolved["decision"], resolved["by"], resolved["client"]) == \
            ("approved", "client", "delfin-dashboard")


def test_a_late_reply_never_answers_the_next_question(service):
    a, b = service.client(), service.client()
    sid = service.session(a, permissionMode="default")
    b.request("session.attach", {"sessionId": sid})
    mark = len(a.events)
    a.request("session.prompt", {"sessionId": sid, "text": "edit it"})
    first, _ = a.wait_for("permission.asked", after=mark)
    a.request("permission.reply", {"sessionId": sid, "requestId": first["requestId"], "choice": "n"})
    a.wait_for(lambda e: e.get("type") == "turn.ended" and e.get("session") == sid, after=mark)
    mark = len(a.events)
    a.request("session.prompt", {"sessionId": sid, "text": "edit it again"})
    second, _ = a.wait_for("permission.asked", after=mark)
    assert second["requestId"] != first["requestId"]
    with pytest.raises(ServiceError) as late:
        b.request("permission.reply", {"sessionId": sid, "requestId": first["requestId"], "choice": "y"})
    assert late.value.code == protocol.UNKNOWN_REQUEST
    assert sid in {s["sessionId"] for s in a.request("sessions.list")["sessions"]
                   if s["pending"] == 1}


# -- people mode ----------------------------------------------------------------------

def _person(tmp_path: Path, name: str) -> ServiceProcess:
    root = tmp_path / name
    root.mkdir()
    return ServiceProcess(root)


def test_a_person_sees_and_reaches_only_their_own_service(tmp_path, monkeypatch):
    from delfin.agent.service import rundir

    anna, ben = _person(tmp_path, "anna"), _person(tmp_path, "ben")
    try:
        a = anna.client()
        b = ben.client()
        anna_sid = anna.session(a, title="anna's scan")
        ben_sid = ben.session(b, title="ben's scan")
        # Anna's window: her runtime directory, as the launcher's people mode
        # exports it (DELFIN_RUNTIME_DIR=<her home>/.run).
        monkeypatch.setenv("DELFIN_RUNTIME_DIR", str(anna.runtime))
        assert rundir.socket_path() == anna.socket
        with ServiceClient(client="delfin-dashboard") as mine:
            listed = mine.request("sessions.list")
            ids = {s["sessionId"] for s in listed["sessions"]}
            assert ids == {anna_sid}
            assert ben_sid not in {e.get("sessionId") for e in listed["elsewhere"]}
            with pytest.raises(ServiceError) as attach:
                mine.request("session.attach", {"sessionId": ben_sid})
            assert attach.value.code == protocol.UNKNOWN_SESSION
            with pytest.raises(ServiceError) as rename:
                mine.request("session.setTitle", {"sessionId": ben_sid, "title": "taken over"})
            assert rename.value.code == protocol.UNKNOWN_SESSION
            with pytest.raises(ServiceError):
                mine.request("permission.reply", {"sessionId": ben_sid, "requestId": "p1-1",
                                                  "choice": "y"})
        assert _summary(b, ben_sid)["title"] == "ben's scan"
        # The socket is the user's alone (0600, in a 0700 directory).
        assert oct(os.stat(anna.socket).st_mode & 0o777) == "0o600"
        assert oct(os.stat(anna.socket.parent).st_mode & 0o777) == "0o700"
        assert json.loads(json.dumps(listed))      # plain data, nothing else
    finally:
        anna.stop()
        ben.stop()
