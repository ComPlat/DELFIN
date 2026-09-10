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

import threading
import time

import pytest

from delfin.dashboard import session as S


@pytest.fixture(autouse=True)
def clean():
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

def test_the_resume_url_keeps_host_port_and_token():
    url = S.resume_url(
        "uc3n990-ab12",
        request_url="http://uc3n990:8866/voila/render/dash.ipynb?token=abc123")
    assert url == "http://uc3n990:8866/delfin/resume/uc3n990-ab12?token=abc123"


def test_the_resume_url_is_not_the_dashboard_url():
    """The normal address goes through Voila's renderer, which EXECUTES
    the notebook — against a live kernel that would run everything a
    second time instead of showing what is there."""
    assert "/voila/render/" not in S.RESUME_PATH
    url = S.resume_url("n", request_url="http://h:1/voila/render/d.ipynb")
    assert "/voila/render/" not in url


@pytest.mark.parametrize("name,req", [
    ("", "http://h:1/x"),
    ("n", ""),
    ("n", "not a url"),
    ("n", "///"),
])
def test_an_unknown_address_is_empty_not_a_guess(name, req):
    assert S.resume_url(name, request_url=req) == ""


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
