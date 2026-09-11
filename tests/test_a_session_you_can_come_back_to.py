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

import re

import json
import os
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


# ---------------------------------------------------------------------------
# The safety property: never fire before the page has spoken
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# The opt-in
# ---------------------------------------------------------------------------


def test_an_armed_session_says_so():
    """describe() is what the strip and the terminal are built from, so an
    armed session cannot be invisible."""
    S.keep_alive(True, session_name="night-run")
    told = S.describe()
    assert told["kept_alive"] is True
    assert told["session_name"] == "night-run"
    assert told["pid"] == __import__("os").getpid()
    S.keep_alive(False)
    assert S.describe()["kept_alive"] is False


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


def test_both_roots_are_registered_and_both_are_shown():
    src = _dashboard_source()
    assert "_session.register_root(_header_root, _body_root)" in src
    assert "display(_header_root)" in src
    assert "display(_body_root)" in src


# ---------------------------------------------------------------------------
# The page's end of the heartbeat
# ---------------------------------------------------------------------------


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

def test_an_unarmed_strip_is_a_dot_that_says_the_session_ends():
    html = S._strip_html(False, "", "")
    assert "delfin-session-dot" in html and "dot on" not in html
    assert 'title="Session ends when the window closes"' in html
    assert re.sub(r"<[^>]+>", "", html).strip() == "", "the strip shows text beside the dot"


def test_an_armed_strip_is_a_dot_that_links_to_the_return_address():
    html = S._strip_html(True, "uc3n990-ab12", "http://h:8866/delfin/resume/x")
    assert "dot on" in html
    assert 'href="http://h:8866/delfin/resume/x"' in html
    assert "Kept as uc3n990-ab12" in html
    assert re.sub(r"<[^>]+>", "", html).strip() == "", "the strip shows text beside the dot"
    assert not re.search(r"[äöüÄÖÜß]", html)


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
    assert "End:" in out


# ---------------------------------------------------------------------------
# Landing while a session is still running
# ---------------------------------------------------------------------------
#
# Offered where people land rather than behind a route of its own: they
# open the address they always open, and it tells them. That is also why
# none of this needs a server extension.

def _kept(name, kid, *, hours_ago=1.0, url="http://h:8866/voila/render/x?token=t"):
    import json
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
    assert "for 5" in html


def test_a_young_session_is_counted_in_minutes():
    """A dashboard opened three minutes ago should not read "for 0 h"."""
    import re

    _kept("fresh-1", "kernel-aaa", hours_ago=0.05)
    html = S._banner_html(S.other_sessions(exclude_kernel="kernel-bbb"))
    age = re.search(r"for[^&]*&nbsp;(min|h)", html)
    assert age and age.group(1) == "min", html


def test_a_record_without_an_address_says_so_rather_than_linking_nowhere():
    _kept("no-url", "kernel-aaa", url="")
    html = S._banner_html(S.other_sessions(exclude_kernel="kernel-bbb"))
    assert "address unknown" in html
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

    header = _header_root_segment()
    assert "_session_strip" in header, (
        "the session strip is built but not in the header the dashboard shows"
    )


def test_the_dashboard_builds_the_strip_and_nothing_of_the_old_mechanism():
    src = _dashboard_source()
    assert "_session.build_status_strip()" in src
    for gone in ("build_heartbeat_widget", "heartbeat_js", "start_watchdog"):
        assert gone not in src, f"{gone} is back; the server ends kernels now"


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


# ---------------------------------------------------------------------------
# Three things a browser found after the mechanism "worked"
# ---------------------------------------------------------------------------


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
                  if getattr(w, "description", "") == "Keep session")
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
    assert "probe-10" in out.getvalue() and "End:" in out.getvalue()
    assert capsys.readouterr().out == ""


# ---------------------------------------------------------------------------
# A return that found nothing says why
# ---------------------------------------------------------------------------

def test_a_return_without_an_address_names_voila(monkeypatch):
    monkeypatch.delenv("VOILA_REQUEST_URL", raising=False)
    why = S.why_not_resumed()
    assert "VOILA_REQUEST_URL" in why


def test_a_return_without_a_record_names_the_record_dir(monkeypatch, tmp_path):
    monkeypatch.setattr(S, "RECORD_DIR", str(tmp_path))
    why = S.why_not_resumed(request_url="http://h:8866/voila/render/x.ipynb?session=uc3n990-ab12")
    assert "No record" in why and "uc3n990-ab12" in why and str(tmp_path) in why


def test_a_return_whose_kernel_is_gone_says_so(monkeypatch, tmp_path):
    monkeypatch.setattr(S, "RECORD_DIR", str(tmp_path))
    S.write_record("uc3n990-ab12", kid="aaaa1111-0000-4000-8000-000000000001")
    monkeypatch.setattr(S, "kernel_id", lambda: "bbbb2222-0000-4000-8000-000000000002")
    why = S.why_not_resumed(request_url="http://h:8866/voila/render/x.ipynb?session=uc3n990-ab12")
    assert "exists" in why and "fresh kernel" in why
    assert "aaaa1111" in why and "bbbb2222" in why


def test_every_word_the_session_shows_is_english():
    """Asked for on 2026-09-11: the dashboard speaks English. The strip,
    the toggle, the banner and the terminal line."""
    import inspect
    src = inspect.getsource(S)
    for word in ("Sitzung", "Offen halten", "Läuft weiter", "zurück über",
                 "Adresse unbekannt", "wieder hineingehen", "bleibt bestehen"):
        assert word not in src, word


# ---------------------------------------------------------------------------
# A record that could not be written says so
# ---------------------------------------------------------------------------

def test_a_record_that_cannot_be_written_names_the_reason(monkeypatch, tmp_path):
    """Seen on a cluster: the terminal said the session was kept and gave
    a return address, the server found no record. write_record had
    swallowed the OSError."""
    blocker = tmp_path / "not-a-dir"
    blocker.write_text("x")                      # a file where the directory should be
    monkeypatch.setattr(S, "RECORD_DIR", str(blocker / "kept_sessions"))
    monkeypatch.setattr(S, "kernel_id", lambda: "aaaa1111-0000-4000-8000-000000000001")
    assert S.write_record("uc3n990-ab12") == ""
    why = S.last_write_error()
    assert "Error" in why and str(blocker) in why


def test_a_failed_chmod_does_not_lose_the_record(monkeypatch, tmp_path):
    monkeypatch.setattr(S, "RECORD_DIR", str(tmp_path))
    monkeypatch.setattr(S, "kernel_id", lambda: "aaaa1111-0000-4000-8000-000000000001")
    monkeypatch.setattr(S.os, "chmod", lambda *a, **k: (_ for _ in ()).throw(OSError("no chmod here")))
    path = S.write_record("uc3n990-ab12")
    assert path and S.last_write_error() == ""
    assert any(r["session_name"] == "uc3n990-ab12" for r in S.list_records(root=str(tmp_path)))


def test_the_toggle_reports_a_failed_write_and_does_not_stay_armed(monkeypatch, tmp_path, capsys):
    blocker = tmp_path / "not-a-dir"
    blocker.write_text("x")
    monkeypatch.setattr(S, "RECORD_DIR", str(blocker / "kept_sessions"))
    monkeypatch.setattr(S, "kernel_id", lambda: "")        # outside a kernel: print, not the server fd
    monkeypatch.setattr(S, "write_record", lambda name: "")
    monkeypatch.setattr(S, "last_write_error", lambda: "PermissionError: [Errno 13] (record dir /x)")
    strip = S.build_status_strip()
    toggle = next(w for w in strip.children if getattr(w, "description", "") == "Keep session")
    note = next(w for w in strip.children if w.__class__.__name__ == "HTML" and 'class="delfin-session-strip"' in (w.value or ""))
    toggle.value = True
    assert toggle.value is False, "the control stayed armed with nothing on disk"
    assert "delfin-session-dot failed" in note.value
    assert "Could not keep the session" in note.value and "PermissionError" in note.value
    assert not S.is_kept_alive()
    out = capsys.readouterr().out
    assert "could NOT be kept" in out and "PermissionError" in out


def test_a_missing_record_report_lists_what_is_on_disk(monkeypatch, tmp_path):
    monkeypatch.setattr(S, "RECORD_DIR", str(tmp_path))
    monkeypatch.setattr(S, "kernel_id", lambda: "bbbb2222-0000-4000-8000-000000000002")
    (tmp_path / "other-1.json").write_text('{"session_name": "other-1", "kernel_id": "x"}')
    why = S.why_not_resumed(request_url="http://h:8866/voila/render/x.ipynb?session=uc3n990-ab12")
    assert "No record" in why and "holds 1 record" in why and "other-1.json" in why
    assert "HOME=" in why and "bbbb2222" in why
    monkeypatch.setattr(S, "RECORD_DIR", str(tmp_path / "nowhere"))
    why = S.why_not_resumed(request_url="http://h:8866/voila/render/x.ipynb?session=uc3n990-ab12")
    assert "does not exist" in why


def test_removing_a_record_says_who_asked(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(S, "RECORD_DIR", str(tmp_path))
    monkeypatch.setattr(S, "kernel_id", lambda: "")
    S.write_record("uc3n990-ab12", kid="aaaa1111-0000-4000-8000-000000000001")
    assert S.drop_record("uc3n990-ab12")
    out = capsys.readouterr().out
    assert 'record of session "uc3n990-ab12" removed' in out
    assert "test_removing_a_record_says_who_asked" in out
    assert not S.drop_record("uc3n990-ab12")            # gone already: nothing said twice
