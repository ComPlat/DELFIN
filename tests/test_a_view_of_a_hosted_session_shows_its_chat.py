"""A dashboard's view of a session another dashboard hosts shows that
dashboard's chat, not a second rendering of it.

The hosting tab's chat HTML (dashboard.chat: whole on attach, then diffs)
is what the view's chat widget shows -- same markup, so the page's chat
scripts draw the same streaming text, tool rows, cards and 3D viewers.
Open questions still come as permission.asked and get their buttons; the
composer still sends to the session. A session a worker runs (no host)
is still rendered from its events.
"""
from __future__ import annotations

import time

import pytest

pytest.importorskip("ipywidgets")

from delfin.dashboard import live_sessions as LS                 # noqa: E402
from delfin.dashboard import live_tab                            # noqa: E402
from fake_host import FakeHost                                   # noqa: E402
from tui_service import ServiceProcess, turn                     # noqa: E402

SID = "c0ffee00c0ffee00c0ffee00c0ffee00"
FIRST = ('<template class="delfin-chat-src"><div class="delfin-agent-chat">'
         '<div class="delfin-chat-msg delfin-chat-user">hi</div></div></template>')
CARD = ("<div class='delfin-chat-msg delfin-chat-tool'>"
        "<div class='delfin-mol-card' data-fmt='xyz'>3D</div></div>")


def _until(cond, timeout: float = 20.0):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        got = cond()
        if got:
            return got
        time.sleep(0.05)
    return cond()


@pytest.fixture
def service(tmp_path):
    svc = ServiceProcess(tmp_path)
    try:
        yield svc
    finally:
        svc.stop()


def test_the_view_shows_the_hosting_tabs_chat_and_follows_it(service):
    host = FakeHost(service.socket, SID, str(service.workspace()), keep=True)
    host.emit("dashboard.chat", v=1, html=FIRST)
    client = LS.LiveClient(service.socket)
    service.clients.append(client)
    _w, refs = live_tab.build_live_tab(client, SID, page_scripts=False)
    try:
        assert _until(lambda: refs["chat"].value == FIRST), refs["chat"].value[:200]
        at = FIRST.index("</div></template>")
        new = FIRST[:at] + CARD + FIRST[at:]
        host.emit("dashboard.chat", v=2, base=1, keep=at, tail=new[at:])
        assert _until(lambda: refs["chat"].value == new)
        host.emit("text.ended", turnId="t", partId="p", text="said by the protocol")
        time.sleep(0.5)
        assert refs["chat"].value == new, "the mirror is what is shown, not a second rendering"
        host.emit("permission.asked", requestId="r9", kind="confirm", tool="Bash",
                  preview="rm -rf build", options=[{"key": "y", "label": "Approve"},
                                                   {"key": "n", "label": "Deny"}], flags={})
        assert _until(lambda: len(refs["approval_box"].children) == 1)
        refs["input"].value = "from the view"
        refs["send"]()
        assert host.wait_asked("session.prompt")["text"] == "from the view"
    finally:
        refs["close"]()
        host.close()


def test_a_session_a_worker_runs_is_still_rendered_from_its_events(service):
    client = LS.LiveClient(service.socket)
    service.clients.append(client)
    tui = service.client()
    sid = service.session(tui)
    turn(tui, sid, "hello there")
    _w, refs = live_tab.build_live_tab(client, sid, page_scripts=False)
    try:
        assert _until(lambda: "hello there" in refs["chat"].value)
    finally:
        refs["close"]()
