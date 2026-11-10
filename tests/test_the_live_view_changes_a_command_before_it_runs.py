"""A shell approval in the dashboard's live view: change the command, and
the second click of a lasting rule.

"change the command" (key c) approves once the person's own version of the
command, never the one asked (approvals.answer: command with choice c).
The live view's button opens the command to edit and sends what was
edited; without an editable command (it holds a secret) it is not offered.
A lasting rule (A, k) needs a second click within CONFIRM_S; after that the
button says again what it does.

Real path: a private service, the scripted engine's shell request, the
live view's own buttons. Control: "yes" approves as asked.
"""

from __future__ import annotations

import time

import pytest

pytest.importorskip("ipywidgets")

from delfin.agent.service.client import ServiceClient          # noqa: E402
from delfin.dashboard import live_sessions as LS                # noqa: E402
from tui_service import ServiceProcess                          # noqa: E402


def _until(cond, timeout: float = 30.0):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        got = cond()
        if got:
            return got
        time.sleep(0.05)
    return cond()


@pytest.fixture
def setup(tmp_path):
    svc = ServiceProcess(tmp_path)
    tui = ServiceClient(svc.socket, client="delfin-tui")
    svc.clients.append(tui)
    dash = LS.LiveClient(svc.socket)
    try:
        yield svc, tui, dash
    finally:
        dash.close()
        svc.stop()


def _asked(refs):
    return _until(lambda: refs["approval_box"].children[0] if refs["approval_box"].children else None)


def _button(card, key):
    return next((b for row in card.children if hasattr(row, "children")
                 for b in row.children if getattr(b, "tooltip", "").startswith(f"{key}:")), None)


def test_change_the_command_sends_the_edited_command(setup):
    from delfin.dashboard import live_tab

    svc, tui, dash = setup
    sid = svc.session(tui, permissionMode="default")
    tui.request("session.attach", {"sessionId": sid})
    _w, refs = live_tab.build_live_tab(dash, sid)
    mark = len(tui.events)
    tui.request("session.prompt", {"sessionId": sid, "text": "bash please"})
    asked, _ = tui.wait_for("permission.asked", after=mark)
    assert asked.get("command")
    card = _asked(refs)
    change = _button(card, "c")
    assert change is not None
    change.click()
    edit = refs["change_input"]
    assert refs["change_box"].layout.display != "none" and edit.value == asked["command"]
    edit.value = "echo changed by the person"
    refs["change_run"].click()
    resolved, _ = tui.wait_for("permission.resolved", after=mark)
    assert resolved["decision"] == "changed" and resolved["client"] == "delfin-dashboard"
    refs["close"]()


def test_yes_is_the_control(setup):
    from delfin.dashboard import live_tab

    svc, tui, dash = setup
    sid = svc.session(tui, permissionMode="default")
    tui.request("session.attach", {"sessionId": sid})
    _w, refs = live_tab.build_live_tab(dash, sid)
    mark = len(tui.events)
    tui.request("session.prompt", {"sessionId": sid, "text": "bash please"})
    tui.wait_for("permission.asked", after=mark)
    _button(_asked(refs), "y").click()
    resolved, _ = tui.wait_for("permission.resolved", after=mark)
    assert resolved["decision"] == "approved"
    refs["close"]()


def test_a_lasting_rule_button_says_again_what_it_does_after_its_window(setup, monkeypatch):
    from delfin.dashboard import live_tab

    monkeypatch.setattr(live_tab, "CONFIRM_S", 0.4)
    svc, tui, dash = setup
    sid = svc.session(tui, permissionMode="default")
    tui.request("session.attach", {"sessionId": sid})
    _w, refs = live_tab.build_live_tab(dash, sid)
    mark = len(tui.events)
    tui.request("session.prompt", {"sessionId": sid, "text": "bash please"})
    tui.wait_for("permission.asked", after=mark)
    always = _button(_asked(refs), "A")
    label = always.description
    always.click()
    assert always.description != label                       # armed: asks again
    assert _until(lambda: always.description == label, 5), always.description
    # Past its window the next click arms again; it does not save a rule.
    always.click()
    time.sleep(0.2)
    assert not [e for e in tui.events[mark:] if e.get("type") == "permission.resolved"]
    refs["close"]()
