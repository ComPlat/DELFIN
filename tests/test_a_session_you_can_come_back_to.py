"""Closing the window ends the session — unless somebody said not to.

Voila gives every browser connection a fresh kernel and re-executes the
notebook, so a reload is a restart and closing the window ends the run.
That stays the default: it is what people expect and nothing here
changes it.

What was missing is a way to say "not this one". The state people lose
is not the agent's alone — 2054 widget constructions across nineteen
tabs, of which the agent is 198 — plus the Python objects beside them:
parsed calculations, scanned archives, a half-filled submit form. None
of it has to be SAVED if the kernel simply does not go away, and that is
why the mechanism is an inversion rather than a snapshot: the kernel
watches its own page and ends itself, and the opt-in says only "do not".

The property this file cares about most is the third one. A watchdog
that fires before the page has ever spoken would tear down the default
path for anybody whose frontend did not load the heartbeat — so it arms
on the FIRST beat and never before.
"""

from __future__ import annotations

import json
import os
import threading
import time
from pathlib import Path

import pytest

from delfin.dashboard import session as S


@pytest.fixture(autouse=True)
def clean(tmp_path, monkeypatch):
    # conftest redirects RECORD_DIR for the whole suite; this pins it per
    # test as well, so one test's armed session cannot appear in the
    # next one's landing page.
    monkeypatch.setattr(S, "RECORD_DIR", str(tmp_path / "kept"))
    S._reset_for_tests()
    yield
    S._reset_for_tests()


# ---------------------------------------------------------------------------
# The default: the window closes, the session ends
# ---------------------------------------------------------------------------

def test_nothing_is_kept_alive_by_default():
    assert S.is_kept_alive() is False
    assert S.describe()["kept_alive"] is False


def test_a_page_that_stops_beating_ends_the_kernel():
    ended = threading.Event()
    S.beat()
    S.GRACE_SECONDS, old = 0.05, S.GRACE_SECONDS
    try:
        S.start_watchdog(shutdown=ended.set, poll_seconds=0.01)
        assert ended.wait(3.0), "the kernel outlived its page"
    finally:
        S.GRACE_SECONDS = old


def test_a_beating_page_is_left_alone():
    ended = threading.Event()
    S.beat()
    S.GRACE_SECONDS, old = 0.4, S.GRACE_SECONDS
    try:
        S.start_watchdog(shutdown=ended.set, poll_seconds=0.02)
        for _ in range(20):          # keep beating for ~1s
            S.beat()
            time.sleep(0.05)
        assert not ended.is_set(), "a live page was torn down"
    finally:
        S.GRACE_SECONDS = old


# ---------------------------------------------------------------------------
# The safety property: never fire before the page has spoken
# ---------------------------------------------------------------------------

def test_a_page_that_never_beats_is_not_a_page_that_left():
    """A frontend whose scripts did not load cannot send heartbeats. It
    must not be read as a window that closed, or the default path breaks
    for everyone it happens to."""
    ended = threading.Event()
    S.GRACE_SECONDS, old = 0.05, S.GRACE_SECONDS
    try:
        S.start_watchdog(shutdown=ended.set, poll_seconds=0.01)
        time.sleep(0.6)
        assert not ended.is_set(), (
            "the watchdog fired before the page had ever spoken")
        assert S.describe()["watchdog_armed"] is False
    finally:
        S.GRACE_SECONDS = old


def test_the_watchdog_arms_on_the_first_beat():
    assert S.describe()["watchdog_armed"] is False
    S.beat()
    assert S.describe()["watchdog_armed"] is True


# ---------------------------------------------------------------------------
# The opt-in
# ---------------------------------------------------------------------------

def test_an_armed_session_survives_its_page():
    ended = threading.Event()
    S.beat()
    S.keep_alive(True, session_name="probe-1")
    S.GRACE_SECONDS, old = 0.05, S.GRACE_SECONDS
    try:
        S.start_watchdog(shutdown=ended.set, poll_seconds=0.01)
        time.sleep(0.6)
        assert not ended.is_set(), "an armed session was torn down"
    finally:
        S.GRACE_SECONDS = old


def test_disarming_lets_it_end_again():
    ended = threading.Event()
    S.beat()
    S.keep_alive(True, session_name="probe-2")
    S.GRACE_SECONDS, old = 0.05, S.GRACE_SECONDS
    try:
        S.start_watchdog(shutdown=ended.set, poll_seconds=0.01)
        time.sleep(0.3)
        assert not ended.is_set()
        S.keep_alive(False)
        assert ended.wait(3.0), "disarming did not let the session end"
    finally:
        S.GRACE_SECONDS = old


def test_an_armed_session_says_so():
    """Kept alive and invisible is how a machine ends up holding ten of
    them. Everything the strip and the terminal line need is in one
    call."""
    S.keep_alive(True, session_name="uc3n990-2")
    d = S.describe()
    assert d["kept_alive"] is True
    assert d["session_name"] == "uc3n990-2"
    assert d["pid"] > 0
    assert d["grace_seconds"] == S.GRACE_SECONDS


# ---------------------------------------------------------------------------
# Coming back: nothing is restored, because nothing left
# ---------------------------------------------------------------------------

def test_resume_shows_the_widgets_that_are_still_there():
    shown = []

    class _W:
        def __init__(self, name): self.name = name

    a, b = _W("header"), _W("tabs")
    S.register_root(a, b)

    import IPython.display as ipd
    real = ipd.display
    try:
        ipd.display = lambda *args, **kw: shown.extend(args)
        assert S.resume() is True
    finally:
        ipd.display = real
    assert shown == [a, b], "resume must re-show the same objects, in order"


def test_resume_says_no_when_there_is_nothing_to_show():
    """A kernel that never built a dashboard is a caller error, not a
    state to paper over."""
    assert S.resume() is False


def test_registering_again_replaces_rather_than_appends():
    class _W: pass
    first, second = _W(), _W()
    S.register_root(first)
    S.register_root(second)
    assert S.roots() == [second]


def test_a_none_root_is_dropped():
    class _W: pass
    w = _W()
    S.register_root(w, None)
    assert S.roots() == [w]


# ---------------------------------------------------------------------------
# The watchdog itself
# ---------------------------------------------------------------------------

def test_starting_twice_does_not_race_itself():
    """Two watchdogs would both call shutdown, and the second would do it
    to a kernel that is already going."""
    first = S.start_watchdog(shutdown=lambda: None, poll_seconds=5)
    second = S.start_watchdog(shutdown=lambda: None, poll_seconds=5)
    assert first is second


# ---------------------------------------------------------------------------
# Wired into the dashboard, not just available
# ---------------------------------------------------------------------------
#
# A module nobody calls is a module that does nothing. These read the
# source rather than building a dashboard, because create_dashboard
# wants a live kernel, a calc directory and nineteen tabs — and what
# has to be true here is a property of the wiring, not of a run.

def _dashboard_source() -> str:
    import inspect

    from delfin import dashboard as d
    return inspect.getsource(d)


def test_the_dashboard_registers_what_it_displays():
    src = _dashboard_source()
    assert "_session.register_root(" in src
    # Before the display calls, so a resume can never find a half-built
    # registry.
    assert src.index("_session.register_root(") < src.index("display(_header_root)")


def test_the_dashboard_starts_the_watchdog():
    """Without it the default teardown does not happen at all — the
    kernel would simply live on, which is the opposite of the default."""
    assert "_session.start_watchdog()" in _dashboard_source()


def test_both_roots_are_registered_and_both_are_shown():
    src = _dashboard_source()
    assert "_session.register_root(_header_root, _body_root)" in src
    assert "display(_header_root)" in src
    assert "display(_body_root)" in src


# ---------------------------------------------------------------------------
# The page's end of the heartbeat
# ---------------------------------------------------------------------------

def test_a_write_from_the_page_is_a_beat():
    field = S.build_heartbeat_widget()
    assert S.describe()["watchdog_armed"] is False
    field.value = "1700000000"
    assert S.describe()["watchdog_armed"] is True


def test_the_field_is_hidden_and_findable():
    """Hidden because it is plumbing; classed because the script has to
    find it without knowing the widget's generated ids."""
    field = S.build_heartbeat_widget()
    assert field.layout.display == "none"
    assert S._HEARTBEAT_CLASS in field._dom_classes


def test_the_script_writes_through_the_native_setter():
    """A plain `el.value = x` does not notify a widget: its own handler
    listens for the events a real edit produces. This is the idiom the
    dashboard already uses (see molecule_viewer)."""
    js = S.heartbeat_js()
    assert "getOwnPropertyDescriptor" in js
    assert "new Event('input'" in js and "new Event('change'" in js
    assert S._HEARTBEAT_CLASS in js


def test_the_script_beats_when_a_hidden_tab_comes_back():
    """Browsers throttle timers in a background tab to once a minute or
    worse, which looks exactly like a window that closed."""
    js = S.heartbeat_js()
    assert "visibilitychange" in js
    assert "document.hidden" in js


def test_several_beats_fit_inside_the_grace():
    """A slow link may drop beats; it must not drop the session."""
    assert S.GRACE_SECONDS >= 4 * S.BEAT_SECONDS


def test_the_script_replaces_its_own_timer():
    """Voila re-runs page scripts on some navigations; two intervals
    would double the traffic and outlive each other."""
    js = S.heartbeat_js()
    assert "clearInterval(window.__delfinHeartbeat)" in js


# ---------------------------------------------------------------------------
# Where to come back to
# ---------------------------------------------------------------------------

_RESUME_PATH = "/voila/render/delfin_voila_runtime/delfin_resume.ipynb"


def test_the_resume_url_keeps_host_port_and_token():
    url = S.resume_url(
        "uc3n990-ab12",
        request_url=("http://uc3n990:8866/voila/render/"
                     "delfin_voila_runtime/delfin_dashboard.ipynb?token=abc"),
        resume_path=_RESUME_PATH)
    assert url == ("http://uc3n990:8866" + _RESUME_PATH
                   + "?token=abc&session=uc3n990-ab12")


def test_the_resume_url_is_not_the_dashboard_url():
    """The dashboard's own address renders and EXECUTES the whole
    notebook. Against a live kernel that would run all nineteen tabs a
    second time — the opposite of coming back."""
    url = S.resume_url(
        "n", request_url="http://h:1/voila/render/delfin_dashboard.ipynb",
        resume_path=_RESUME_PATH)
    assert "delfin_dashboard.ipynb" not in url
    assert "delfin_resume.ipynb" in url


def test_a_stale_session_in_the_address_is_replaced_not_doubled():
    """Coming back from a resume page and arming again must not leave two
    session keys for the manager to choose between."""
    url = S.resume_url(
        "second", resume_path=_RESUME_PATH,
        request_url="http://h:1/voila/render/r.ipynb?token=t&session=first")
    assert url.count("session=") == 1
    assert "session=second" in url


@pytest.mark.parametrize("name,req", [
    ("", "http://h:1/x"),
    ("n", ""),
    ("n", "not a url"),
    ("n", "///"),
])
def test_an_unknown_address_is_empty_not_a_guess(name, req):
    assert S.resume_url(name, request_url=req, resume_path=_RESUME_PATH) == ""


def test_no_staged_resume_notebook_means_no_address(monkeypatch):
    """The launcher sets the path. Without it there is nowhere to come
    back to, and a guess would send somebody to a 404."""
    monkeypatch.delenv(S.RESUME_PATH_ENV, raising=False)
    assert S.resume_url("n", request_url="http://h:1/x", resume_path="") == ""


def test_two_sessions_started_together_do_not_collide():
    assert S.default_session_name() != S.default_session_name()


# ---------------------------------------------------------------------------
# The strip
# ---------------------------------------------------------------------------

def test_an_unarmed_strip_says_the_session_ends():
    html = S._strip_html(False, "", "")
    assert "endet" in html
    assert "delfin-session-dot" in html and "dot on" not in html


def test_an_armed_strip_shows_the_name_and_where_to_come_back():
    html = S._strip_html(True, "uc3n990-ab12", "http://h:8866/delfin/resume/x")
    assert "dot on" in html
    assert "uc3n990-ab12" in html
    assert "delfin/resume" in html


def test_the_toggle_arms_and_disarms_the_session():
    strip = S.build_status_strip()
    toggle = next(c for c in strip.children if hasattr(c, "value")
                  and isinstance(getattr(c, "value"), bool))
    assert S.is_kept_alive() is False
    toggle.value = True
    assert S.is_kept_alive() is True
    assert S.session_name(), "an armed session needs a name to come back to"
    toggle.value = False
    assert S.is_kept_alive() is False


def test_the_tooltip_says_what_does_not_come_back():
    """A control that promises more than it delivers is worse than none.
    Scroll positions live in the browser, not the kernel."""
    strip = S.build_status_strip()
    toggle = next(c for c in strip.children
                  if isinstance(getattr(c, "tooltip", None), str)
                  and getattr(c, "tooltip"))
    assert "Scroll" in toggle.tooltip


def test_arming_prints_where_to_come_back(capsys):
    """The server's terminal is usually inside tmux, which is where
    somebody looks tomorrow after `tmux attach`."""
    S.keep_alive(True, session_name="probe-9")
    S.announce()
    out = capsys.readouterr().out
    assert "probe-9" in out
    assert "Beenden" in out


# ---------------------------------------------------------------------------
# Landing while a session is still running
# ---------------------------------------------------------------------------
#
# Offered where people land rather than behind a route of its own: they
# open the address they always open, and it tells them. That is also why
# none of this needs a server extension.

def _kept(name, kid, *, hours_ago=1.0, url="http://h:8866/voila/render/x?token=t"):
    import json
    import os
    import time as _t

    path = S.write_record(name, kid=kid)
    with open(path, encoding="utf-8") as h:
        rec = json.load(h)
    rec["request_url"] = url
    rec["started_at"] = _t.time() - hours_ago * 3600
    with open(path, "w", encoding="utf-8") as h:
        json.dump(rec, h)
    return path


@pytest.fixture(autouse=True)
def _resume_path(monkeypatch):
    monkeypatch.setenv(
        S.RESUME_PATH_ENV,
        "/voila/render/delfin_voila_runtime/delfin_resume.ipynb")


def test_a_new_window_is_told_about_a_running_session():
    _kept("uc3n990-ab12", "kernel-aaa")
    others = S.other_sessions(exclude_kernel="kernel-bbb")
    assert [r["session_name"] for r in others] == ["uc3n990-ab12"]


def test_the_session_you_are_looking_at_is_not_offered_back():
    """A resume renders INTO the kept kernel, so the page you are on IS
    that session — offering to go back to it would be a loop."""
    _kept("uc3n990-ab12", "kernel-aaa")
    assert S.other_sessions(exclude_kernel="kernel-aaa") == []


def test_nothing_is_shown_when_nothing_is_running():
    assert S._banner_html([]) == ""
    assert S.build_returning_banner(exclude_kernel="kernel-x") is None


def test_the_banner_carries_a_link_and_an_age():
    _kept("uc3n990-ab12", "kernel-aaa", hours_ago=5)
    html = S._banner_html(S.other_sessions(exclude_kernel="kernel-bbb"))
    assert "uc3n990-ab12" in html
    assert "delfin_resume.ipynb" in html
    assert "seit 5" in html


def test_a_young_session_is_counted_in_minutes():
    """A dashboard opened three minutes ago should not read "seit 0 h"."""
    import re

    _kept("fresh-1", "kernel-aaa", hours_ago=0.05)
    html = S._banner_html(S.other_sessions(exclude_kernel="kernel-bbb"))
    age = re.search(r"seit[^&]*&nbsp;(min|h)", html)
    assert age and age.group(1) == "min", html


def test_a_record_without_an_address_says_so_rather_than_linking_nowhere():
    _kept("no-url", "kernel-aaa", url="")
    html = S._banner_html(S.other_sessions(exclude_kernel="kernel-bbb"))
    assert "Adresse unbekannt" in html
    assert "<a href" not in html


def test_a_broken_timestamp_does_not_take_the_banner_down():
    import json

    path = _kept("odd-1", "kernel-aaa")
    with open(path, encoding="utf-8") as h:
        rec = json.load(h)
    rec["started_at"] = "not a number"
    with open(path, "w", encoding="utf-8") as h:
        json.dump(rec, h)
    html = S._banner_html(S.other_sessions(exclude_kernel="kernel-bbb"))
    assert "odd-1" in html


def test_the_dashboard_shows_it_above_everything_else():
    import inspect

    from delfin import dashboard as d
    src = inspect.getsource(d)
    assert "_session.build_returning_banner()" in src
    assert "([_returning] if _returning else [])" in src, (
        "no running session must not leave an empty box in the header")


def _header_root_segment() -> str:
    """The header the dashboard actually puts on the page, as source."""
    import ast
    import inspect

    from delfin import dashboard as d

    src = inspect.getsource(d.create_dashboard)
    tree = ast.parse(src)
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and any(
            isinstance(t, ast.Name) and t.id == "_header_root"
            for t in node.targets
        ):
            return ast.get_source_segment(src, node.value) or ""
    raise AssertionError("create_dashboard no longer builds a _header_root")


def test_the_control_that_keeps_a_session_is_on_the_page():
    """Built, tested, and never mounted — which a unit test cannot see.

    The module had the strip, the hidden field and the script, and every
    one of them had a test. The dashboard displayed none of them, so a
    browser got a page with no way to arm a session at all. Reaching the
    page is a property of the assembler, so it is pinned here.
    """
    import re

    header = _header_root_segment()
    assert "_session_strip" in header, (
        "the session strip is built but not in the header the dashboard shows"
    )
    assert re.search(r"\b_heartbeat\b", header), (
        "the heartbeat field is not displayed; a page that cannot beat "
        "reads as a window that closed"
    )
    assert re.search(r"\b_heartbeat_js\b", header), (
        "the heartbeat script is not sent, so nothing writes to the field"
    )


def test_the_dashboard_builds_all_three_pieces():
    src = _dashboard_source()
    for call in (
        "_session.build_status_strip()",
        "_session.build_heartbeat_widget()",
        "_session.heartbeat_js()",
    ):
        assert call in src, f"create_dashboard never calls {call}"


def test_the_heartbeat_script_is_not_sent_through_run_js():
    """``ctx.run_js`` clears its output before writing.

    The beat has to keep running for the life of the page, so it gets an
    Output of its own; sending it through the shared one would let the
    next startup script wipe it.
    """
    src = _dashboard_source()
    assert "run_js(_session.heartbeat_js" not in src
    assert "display(Javascript(_session.heartbeat_js()))" in src


# ---------------------------------------------------------------------------
# Ending a session must end it, not hand it back
# ---------------------------------------------------------------------------
#
# A kernel that exits on its own reads to the server as a crash, and the
# restarter puts a replacement in its place under the same id. Under
# Voila that replacement is blank -- the notebook is not re-executed --
# and a blank page never beats, so the watchdog never arms in it and it
# stays for the life of the server. Driving this in a browser is what
# showed it: the kernel id was still listed after the teardown, under a
# process that had started seconds earlier.

def test_the_shutdown_address_is_built_from_what_the_server_gave_us(monkeypatch):
    from delfin.dashboard import session as s

    monkeypatch.setenv(
        "VOILA_REQUEST_URL",
        "http://127.0.0.1:8890/voila/render/x.ipynb?token=abc",
    )
    monkeypatch.setenv("JUPYTER_TOKEN", "abc")
    monkeypatch.setattr(s, "kernel_id", lambda: "kid-1")

    url = s.server_shutdown_url()
    assert url == "http://127.0.0.1:8890/api/kernels/kid-1?token=abc"


def test_the_port_stands_in_when_no_request_was_recorded(monkeypatch):
    from delfin.dashboard import session as s

    monkeypatch.delenv("VOILA_REQUEST_URL", raising=False)
    monkeypatch.setenv("VOILA_APP_PORT", "8899")
    monkeypatch.setenv("JUPYTER_TOKEN", "t")
    monkeypatch.setattr(s, "kernel_id", lambda: "kid-2")

    assert s.server_shutdown_url() == (
        "http://127.0.0.1:8899/api/kernels/kid-2?token=t"
    )


def test_outside_a_kernel_there_is_nobody_to_ask(monkeypatch):
    """Every test and every CLI call is outside a kernel."""
    from delfin.dashboard import session as s

    monkeypatch.setattr(s, "kernel_id", lambda: "")
    assert s.server_shutdown_url() == ""
    assert s._ask_server_to_end_this_kernel() is False


def test_a_nonsense_port_is_not_an_address(monkeypatch):
    from delfin.dashboard import session as s

    monkeypatch.delenv("VOILA_REQUEST_URL", raising=False)
    monkeypatch.setenv("VOILA_APP_PORT", "not-a-port")
    monkeypatch.setattr(s, "kernel_id", lambda: "kid-3")
    assert s.server_shutdown_url() == ""


def test_the_delete_carries_the_token_in_the_header_too(monkeypatch):
    """Query token and Authorization header, because a server may be
    configured to accept only one of them."""
    from delfin.dashboard import session as s

    seen = {}

    class _Resp:
        status = 204

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    def _urlopen(req, timeout=0):
        seen["method"] = req.get_method()
        seen["url"] = req.full_url
        seen["auth"] = req.get_header("Authorization")
        return _Resp()

    monkeypatch.setenv("VOILA_REQUEST_URL", "http://h:1/x?token=tok")
    monkeypatch.setenv("JUPYTER_TOKEN", "tok")
    monkeypatch.setattr(s, "kernel_id", lambda: "k9")
    import urllib.request

    monkeypatch.setattr(urllib.request, "urlopen", _urlopen)

    assert s._ask_server_to_end_this_kernel() is True
    assert seen["method"] == "DELETE"
    assert seen["url"] == "http://h:1/api/kernels/k9?token=tok"
    assert seen["auth"] == "token tok"


def test_a_server_that_refuses_is_not_taken_for_a_shutdown(monkeypatch):
    from delfin.dashboard import session as s

    def _boom(req, timeout=0):
        raise OSError("connection refused")

    monkeypatch.setenv("VOILA_REQUEST_URL", "http://h:1/x")
    monkeypatch.setattr(s, "kernel_id", lambda: "k9")
    import urllib.request

    monkeypatch.setattr(urllib.request, "urlopen", _boom)
    assert s._ask_server_to_end_this_kernel() is False


# ---------------------------------------------------------------------------
# Three things a browser found after the mechanism "worked"
# ---------------------------------------------------------------------------

def test_the_first_beat_does_not_wait_for_the_interval():
    """A window closed inside its first ten seconds never beat, so its
    kernel was never torn down. The script now retries the first beat
    quickly until the field it writes to is there."""
    from delfin.dashboard import session as s

    js = s.heartbeat_js()
    assert "return true" in js and "return false" in js
    assert "if (!beat())" in js
    assert "250" in js, "the retry must be far quicker than the interval"


def test_disarming_restarts_the_grace(monkeypatch):
    """The last beat on record may be hours old when the option is switched
    off on a resumed page; judging by it would end the kernel under the
    person looking at it."""
    from delfin.dashboard import session as s

    s._reset_for_tests()
    now = {"t": 1000.0}
    monkeypatch.setattr(s.time, "monotonic", lambda: now["t"])

    s.beat()
    s.keep_alive(True, session_name="x")
    now["t"] += 8 * 3600                      # overnight
    assert s.seconds_since_beat() > s.GRACE_SECONDS

    s.keep_alive(False)
    assert s.seconds_since_beat() == 0.0, "disarming must restart the grace"
    s._reset_for_tests()


def test_disarming_a_page_that_never_beat_does_not_arm_the_watchdog():
    from delfin.dashboard import session as s

    s._reset_for_tests()
    s.keep_alive(True, session_name="x")
    s.keep_alive(False)
    assert s.seconds_since_beat() is None
    s._reset_for_tests()


def test_a_record_whose_process_is_gone_is_not_offered(tmp_path, monkeypatch):
    from delfin.dashboard import session as s

    monkeypatch.setattr(s, "_hostname", lambda: "thishost")
    # A pid that cannot exist on Linux: above pid_max's ceiling.
    dead = {"session_name": "gone", "kernel_id": "k1", "pid": 4194305,
            "host": "thishost", "started_at": 1.0}
    live = {"session_name": "here", "kernel_id": "k2", "pid": os.getpid(),
            "host": "thishost", "started_at": 2.0}
    for rec in (dead, live):
        (tmp_path / f"{rec['session_name']}.json").write_text(json.dumps(rec))

    names = [r["session_name"] for r in s.other_sessions(root=str(tmp_path))]
    assert names == ["here"]
    assert not (tmp_path / "gone.json").exists(), "the dead record is dropped"
    assert (tmp_path / "here.json").exists()


def test_a_record_from_another_host_is_left_alone(tmp_path, monkeypatch):
    """A home directory may be shared; a pid means nothing elsewhere."""
    from delfin.dashboard import session as s

    monkeypatch.setattr(s, "_hostname", lambda: "thishost")
    rec = {"session_name": "far", "kernel_id": "k3", "pid": 4194305,
           "host": "otherhost", "started_at": 1.0}
    (tmp_path / "far.json").write_text(json.dumps(rec))
    assert [r["session_name"] for r in s.other_sessions(root=str(tmp_path))] == ["far"]
    assert (tmp_path / "far.json").exists()


def test_a_record_without_a_host_is_left_alone(tmp_path, monkeypatch):
    """Written by an older DELFIN; nothing to judge it by."""
    from delfin.dashboard import session as s

    monkeypatch.setattr(s, "_hostname", lambda: "thishost")
    rec = {"session_name": "old", "kernel_id": "k4", "pid": 4194305,
           "started_at": 1.0}
    (tmp_path / "old.json").write_text(json.dumps(rec))
    assert [r["session_name"] for r in s.other_sessions(root=str(tmp_path))] == ["old"]


def test_a_new_record_names_its_host(tmp_path, monkeypatch):
    from delfin.dashboard import session as s

    monkeypatch.setattr(s, "kernel_id", lambda: "kid")
    monkeypatch.setattr(s, "_hostname", lambda: "thishost")
    path = s.write_record("named", root=str(tmp_path))
    assert json.loads(Path(path).read_text())["host"] == "thishost"


def test_the_strip_cannot_be_shrunk_by_the_header():
    """Seen in a browser: a 130px control rendered 20px wide the moment
    the header row was full, and the address beside it wrapped into a
    70px column. Clickable, unreadable."""
    from delfin.dashboard import session as s

    strip = s.build_status_strip()
    toggle = next(w for w in strip.children
                  if getattr(w, "description", "") == "Offen halten")
    assert toggle.layout.flex == "0 0 auto"
    assert strip.layout.flex == "0 0 auto"
    assert "white-space:nowrap" in s._STRIP_CSS


def test_inside_a_kernel_the_announcement_reaches_the_server(monkeypatch, capsys):
    """A kernel's print goes to the frontend, not the terminal: ipykernel
    captures it. Inside a kernel the line is written to the server's own
    stdout instead -- and a test is not inside a kernel, which is why the
    test above still reads it from print."""
    import io

    class _Out(io.StringIO):
        def close(self):            # `with out:` must not discard the text
            pass

    out = _Out()
    monkeypatch.setattr(S, "kernel_id", lambda: "k1")
    monkeypatch.setattr(S, "_server_stdout", lambda: out)
    S.keep_alive(True, session_name="probe-10")
    S.announce()
    assert "probe-10" in out.getvalue() and "Beenden" in out.getvalue()
    assert capsys.readouterr().out == ""
