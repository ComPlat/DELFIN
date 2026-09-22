"""A page that closed said so; a connection that dropped did not.

The server sees one fact -- no websocket -- where the person at the
browser lived two: they closed the page, or their link went quiet. The
short grace was written for the first (it is really about a reload) and
was being applied to the second, so a sleeping laptop or a ping timeout
ended a session nobody had closed.

A page that is unloaded sends a beacon; silence is a drop. What is
pinned:

  a closed window            the beacon arrived -> the short grace
  a dropped connection       nothing arrived -> far longer
  the long one still ends    an unwatched session is not immortal
  a beacon is about ITS window  an older one does not shorten a later drop
  coming back clears it      a reconnect forgets the last closing
  the beacon can only shorten no route to a longer stay
  the page carries the script pagehide + sendBeacon, naming the kernel
  the script survives the widget  an HTML widget would strip it
"""

from __future__ import annotations

import asyncio
import logging

import pytest

from delfin.dashboard import resume_server as R
from delfin.dashboard import session as S
from delfin.dashboard import turn_record as T
from delfin.dashboard import window_close as W


class _Base:
    def __init__(self):
        self._kernel_connections = {}
        self._kernels = {}
        self.log = logging.getLogger("test")
        self.ended = []

    def __contains__(self, kid):
        return kid in self._kernels

    async def start_kernel(self, **kwargs):
        kid = kwargs.get("kernel_id") or f"k{len(self._kernels) + 1}"
        self._kernels[kid] = object()
        self._kernel_connections[kid] = 0
        return kid

    def notify_connect(self, kid):
        self._kernel_connections[kid] += 1

    def notify_disconnect(self, kid):
        self._kernel_connections[kid] -= 1

    async def shutdown_kernel(self, kid, *a, **k):
        self.ended.append(kid)
        self._kernels.pop(kid, None)
        self._kernel_connections.pop(kid, None)

    async def shutdown_all(self, *a, **k):
        for kid in list(self._kernels):
            await self.shutdown_kernel(kid)


def _run(coro):
    return asyncio.run(coro)


@pytest.fixture
def manager(tmp_path, monkeypatch):
    monkeypatch.setattr(S, "RECORD_DIR", str(tmp_path / "kept"))
    monkeypatch.setattr(T, "RECORD_DIR", str(tmp_path / "turns"))
    monkeypatch.delenv(R.GRACE_ENV, raising=False)
    monkeypatch.delenv(R.DROPPED_ENV, raising=False)
    cls = R.resume_kernel_manager_class(_Base)
    return cls()


def _leave(manager, kid):
    manager.notify_connect(kid)
    manager.notify_disconnect(kid)


def _age(manager, kid, seconds):
    since, had = manager._delfin_unwatched()[kid]
    manager._delfin_unwatched()[kid] = (since - seconds, had)
    for key, when in list(manager._delfin_closed().items()):
        if key == kid:
            manager._delfin_closed()[key] = when - seconds


# -- the two kinds of quiet -------------------------------------------------

def test_a_window_that_closed_gets_the_short_grace(manager):
    kid = _run(manager.start_kernel())
    _leave(manager, kid)
    manager.delfin_window_closed(kid)
    _age(manager, kid, R.GRACE_SECONDS + 1)
    _run(manager.cull_kernel_if_idle(kid))
    assert manager.ended == [kid]


def test_a_connection_that_dropped_is_given_far_longer(manager):
    kid = _run(manager.start_kernel())
    _leave(manager, kid)
    # No beacon: the page is still open somewhere.
    _age(manager, kid, R.GRACE_SECONDS + 1)
    _run(manager.cull_kernel_if_idle(kid))
    assert manager.ended == []


def test_even_a_dropped_connection_does_not_last_for_ever(manager):
    kid = _run(manager.start_kernel())
    _leave(manager, kid)
    _age(manager, kid, R.DROPPED_SECONDS + 1)
    _run(manager.cull_kernel_if_idle(kid))
    assert manager.ended == [kid]


def test_a_beacon_from_an_earlier_window_does_not_shorten_a_later_drop(manager):
    """Two tabs. One is closed; much later the other one's link dies."""
    kid = _run(manager.start_kernel())
    manager.delfin_window_closed(kid)
    manager._delfin_closed()[kid] -= R.DROPPED_SECONDS  # an hour ago, say
    _leave(manager, kid)
    _age(manager, kid, R.GRACE_SECONDS + 1)
    _run(manager.cull_kernel_if_idle(kid))
    assert manager.ended == []


def test_coming_back_forgets_the_last_closing(manager):
    kid = _run(manager.start_kernel())
    _leave(manager, kid)
    manager.delfin_window_closed(kid)
    manager.notify_connect(kid)
    assert kid not in manager._delfin_closed()

    manager.notify_disconnect(kid)
    _age(manager, kid, R.GRACE_SECONDS + 1)
    _run(manager.cull_kernel_if_idle(kid))
    assert manager.ended == []


def test_a_page_still_loading_is_not_a_closed_window(manager):
    """A kernel that never had a window keeps its own allowance."""
    kid = _run(manager.start_kernel())
    manager._delfin_unwatched()[kid] = (0.0, False)
    manager.delfin_window_closed(kid)
    assert manager._delfin_allowance(kid, False) == R.NEVER_CONNECTED_SECONDS


def test_a_running_turn_outranks_a_closed_window(manager):
    """Closing the tab on a working agent does not throw the run away."""
    kid = _run(manager.start_kernel())
    T.mark(True, kid=kid)
    try:
        _leave(manager, kid)
        manager.delfin_window_closed(kid)
        _age(manager, kid, R.GRACE_SECONDS + 1)
        _run(manager.cull_kernel_if_idle(kid))
        assert manager.ended == []
    finally:
        T.mark(False, kid=kid)


def test_the_override_is_read_from_the_environment(manager, monkeypatch):
    monkeypatch.setenv(R.DROPPED_ENV, "30")
    assert R.dropped_seconds() == 30
    kid = _run(manager.start_kernel())
    _leave(manager, kid)
    _age(manager, kid, 31)
    _run(manager.cull_kernel_if_idle(kid))
    assert manager.ended == [kid]


# -- the beacon itself ------------------------------------------------------

def test_the_script_names_the_kernel_and_sends_on_pagehide():
    script = W.beacon_script("kernel-42")
    assert "pagehide" in script
    assert "sendBeacon" in script
    assert '"kernel-42"' in script
    assert W.ROUTE in script
    # A token-authenticated Voila render does not necessarily create the
    # cookie its own beforeunload code expects.  Prime it on the authenticated
    # route, then carry it in the beacon instead of disabling XSRF.
    assert "fetch(" in script and "credentials: 'same-origin'" in script
    assert "_xsrf" in script and "FormData" in script
    # A page restored from the browser's back/forward cache must be able to
    # announce its later, real close as well.
    assert "pageshow" in script and "event.persisted" in script


def test_the_script_is_not_put_where_a_widget_would_strip_it():
    """An HTML widget sanitises its value and would strip the script,
    leaving every close looking like a drop. An Output does not."""
    widgets = pytest.importorskip("ipywidgets")
    shown = []
    import IPython.display as _d

    real = _d.display
    _d.display = lambda *args, **kw: shown.extend(args)
    try:
        out = W.beacon_widget("kernel-42")
    finally:
        _d.display = real
    assert isinstance(out, widgets.Output)
    assert not isinstance(out, widgets.HTML)
    # Displayed INTO the Output, script intact. (An Output captures
    # nothing outside a kernel, so what was displayed is read here
    # rather than what the widget holds.)
    assert shown and "sendBeacon" in str(getattr(shown[0], "data", ""))
    assert "kernel-42" in str(getattr(shown[0], "data", ""))


def test_outside_a_kernel_there_is_no_beacon(monkeypatch):
    monkeypatch.setattr(S, "kernel_id", lambda: "")
    assert W.beacon_widget() is None


def test_the_route_is_added_to_the_server():
    added = []

    class _App:
        settings = {"base_url": "/prefix/"}

        def add_handlers(self, host, handlers):
            added.append((host, handlers))

    class _Server:
        web_app = _App()

    route = W.register(_Server())
    assert route == "/prefix/" + W.ROUTE
    assert added and added[0][1][0][0] == route


def test_a_server_that_cannot_take_the_route_is_not_a_failure():
    class _Server:
        pass

    assert W.register(_Server()) is None


def test_the_handler_tells_the_manager_which_window_closed():
    """The one thing the route does, and the only thing it can do."""
    from tornado.web import Application

    told = []

    class _Manager:
        def __contains__(self, kid):
            return kid == "kernel-42"

        def delfin_window_closed(self, kid):
            told.append(kid)

    handler_cls = W.handler_class()
    app = Application(kernel_manager=_Manager())
    handler = handler_cls.__new__(handler_cls)
    handler.application = app
    handler.get_argument = lambda name, default=None: (
        "kernel-42" if name == W.KERNEL_ARG else default)
    handler.set_status = lambda code: told.append(code)

    handler.finish = lambda: None
    # The decorated post is the authenticated wrapper; the body under it
    # is what runs once a request is authenticated.
    _run(handler_cls.post.__wrapped__(handler))
    assert told == ["kernel-42", 204]


def test_the_close_route_primes_xsrf_without_bypassing_it():
    """The GET mints a cookie; the POST inherits Jupyter's normal check."""
    from tornado.web import Application

    events = []
    handler_cls = W.handler_class()
    assert "check_xsrf_cookie" not in handler_cls.__dict__
    handler_cls.xsrf_token = property(lambda self: events.append("xsrf") or b"token")

    handler = handler_cls.__new__(handler_cls)
    handler.application = Application(kernel_manager=object())
    handler.set_header = lambda name, value: events.append((name, value))
    handler.set_status = lambda code: events.append(code)
    handler.finish = lambda: events.append("finished")

    _run(handler_cls.get.__wrapped__(handler))
    assert events == ["xsrf", ("Cache-Control", "no-store"), 204, "finished"]


def test_the_close_route_rejects_a_post_without_authentication():
    """Knowing the local port and a kernel id is still not access."""
    from types import SimpleNamespace
    from tornado import web

    handler_cls = W.handler_class()
    handler = handler_cls.__new__(handler_cls)
    handler._current_user = None
    handler.request = SimpleNamespace(method="POST")

    with pytest.raises(web.HTTPError) as raised:
        handler_cls.post(handler)
    assert raised.value.status_code == 403


def test_the_close_route_rejects_a_cookie_post_without_xsrf():
    """Authentication alone does not turn a cross-site POST into authority."""
    from types import SimpleNamespace
    from tornado import web
    from tornado.web import Application

    handler_cls = W.handler_class()
    # Exercise JupyterHandler's inherited check deterministically: this is a
    # cookie-authenticated request, not the explicit-token exemption.
    handler_cls.token_authenticated = property(lambda self: False)
    handler_cls.check_origin = lambda self: True
    handler = handler_cls.__new__(handler_cls)
    handler.application = Application(disable_check_xsrf=False)
    handler._jupyter_current_user = object()
    handler.request = SimpleNamespace(method="POST", headers={})
    handler.get_argument = lambda name, default=None: default

    with pytest.raises(web.HTTPError) as raised:
        handler.check_xsrf_cookie()
    assert raised.value.status_code == 403
    assert "_xsrf" in str(raised.value)


def test_the_close_route_discloses_nothing_about_an_unknown_kernel():
    from tornado import web
    from tornado.web import Application

    class _Manager:
        def __contains__(self, kid):
            return False

        def delfin_window_closed(self, kid):
            raise AssertionError("an unknown kernel must not be marked")

    handler_cls = W.handler_class()
    handler = handler_cls.__new__(handler_cls)
    handler.application = Application(kernel_manager=_Manager())
    handler.get_argument = lambda name, default=None: "somebody-elses-kernel"

    with pytest.raises(web.HTTPError) as raised:
        _run(handler_cls.post.__wrapped__(handler))
    assert raised.value.status_code == 404


def test_voilas_shutdown_keeps_auth_and_xsrf_but_uses_the_safe_grace():
    """Patch only the action, never the API handler's security checks."""
    from jupyter_server.base.handlers import APIHandler
    from tornado.web import Application

    events = []

    class _Manager:
        def __contains__(self, kid):
            return kid == "kernel-42"

        def delfin_window_closed(self, kid):
            events.append(("closed", kid))

        async def shutdown_kernel(self, kid):
            events.append(("shutdown", kid))

    class _VoilaHandler(APIHandler):
        async def post(self, kernel_id):
            await self.kernel_manager.shutdown_kernel(kernel_id)

    xsrf_check = _VoilaHandler.check_xsrf_cookie
    assert W.bridge_voila_shutdown(_VoilaHandler)
    patched = _VoilaHandler.post
    assert W.bridge_voila_shutdown(_VoilaHandler)
    assert _VoilaHandler.post is patched
    assert _VoilaHandler.check_xsrf_cookie is xsrf_check
    assert "check_xsrf_cookie" not in _VoilaHandler.__dict__

    handler = _VoilaHandler.__new__(_VoilaHandler)
    handler.application = Application(kernel_manager=_Manager())
    handler.set_status = lambda code: events.append(("status", code))
    handler.finish = lambda: events.append(("finish", None))
    _run(_VoilaHandler.post.__wrapped__(handler, "kernel-42"))

    assert events == [
        ("closed", "kernel-42"),
        ("status", 204),
        ("finish", None),
    ]


@pytest.mark.slow
def test_the_real_server_requires_token_and_xsrf_for_a_close_notice(tmp_path):
    """The complete Jupyter/Voila route, without a browser or auth shortcut."""
    import http.cookiejar
    import json
    import os
    import socket
    import subprocess
    import sys
    import time
    import urllib.error
    import urllib.parse
    import urllib.request
    from pathlib import Path

    token = "delfin-close-route-integration-test-token"
    with socket.socket() as candidate:
        candidate.bind(("127.0.0.1", 0))
        port = int(candidate.getsockname()[1])

    records = tmp_path / "records"
    runtime = tmp_path / "runtime"
    root_dir = tmp_path / "root"
    for directory in (records, runtime, root_dir):
        directory.mkdir()
        directory.chmod(0o700)
    log = tmp_path / "server.log"
    repo = Path(__file__).resolve().parents[1]
    env = dict(os.environ)
    env.update({
        "PYTHONPATH": str(repo),
        "XDG_RUNTIME_DIR": str(runtime),
        "DELFIN_SESSION_RECORD_DIR": str(records),
        "DELFIN_VOILA_ROOT_DIR": str(root_dir),
        "DELFIN_DASHBOARD_STAY_UP": "1",
    })

    with log.open("wb") as output:
        proc = subprocess.Popen(
            [sys.executable, "-m", "delfin.cli_voila", "--port", str(port),
             "--ip", "127.0.0.1", "--token", token, "--no-browser"],
            cwd=str(repo), env=env, stdout=output, stderr=subprocess.STDOUT,
        )

    root = f"http://127.0.0.1:{port}"
    auth = {"Authorization": f"token {token}"}
    kid = ""
    try:
        deadline = time.time() + 90
        while time.time() < deadline:
            if proc.poll() is not None:
                pytest.fail("server exited while starting:\n" +
                            log.read_text(errors="replace")[-2000:])
            try:
                request = urllib.request.Request(root + "/api/kernels", headers=auth)
                with urllib.request.urlopen(request, timeout=5):
                    break
            except Exception:
                time.sleep(0.5)
        else:
            pytest.fail("server did not answer:\n" +
                        log.read_text(errors="replace")[-2000:])

        create = urllib.request.Request(
            root + "/api/kernels", data=b"{}", method="POST",
            headers={**auth, "Content-Type": "application/json"})
        with urllib.request.urlopen(create, timeout=30) as response:
            kid = str(json.load(response)["id"])

        class _NoRedirect(urllib.request.HTTPRedirectHandler):
            def redirect_request(self, *args, **kwargs):
                return None

        anonymous = urllib.request.build_opener(_NoRedirect)
        try:
            anonymous.open(root + "/delfin-api/window-closing", timeout=10)
            anonymous_status = 200
        except urllib.error.HTTPError as exc:
            anonymous_status = exc.code
        assert anonymous_status in (302, 403), anonymous_status

        cookies = http.cookiejar.CookieJar()
        browser = urllib.request.build_opener(
            urllib.request.HTTPCookieProcessor(cookies))
        armed_url = (root + "/delfin-api/window-closing?token=" +
                     urllib.parse.quote(token))
        with browser.open(armed_url, timeout=10) as response:
            assert response.status == 204
        xsrf = next((cookie.value for cookie in cookies
                     if cookie.name == "_xsrf"), "")
        assert xsrf, "the authenticated primer did not mint an XSRF cookie"

        missing = urllib.request.Request(
            root + "/voila/api/shutdown/" + kid, data=b"", method="POST")
        with pytest.raises(urllib.error.HTTPError) as missing_error:
            browser.open(missing, timeout=10)
        assert missing_error.value.code == 403

        body = urllib.parse.urlencode({"_xsrf": xsrf}).encode("ascii")
        accepted = urllib.request.Request(
            root + "/voila/api/shutdown/" + kid, data=body, method="POST",
            headers={"Content-Type": "application/x-www-form-urlencoded"})
        with browser.open(accepted, timeout=10) as response:
            assert response.status == 204

        listed = urllib.request.Request(root + "/api/kernels", headers=auth)
        with urllib.request.urlopen(listed, timeout=10) as response:
            kernels = json.load(response)
        assert any(str(row.get("id")) == kid for row in kernels), (
            "Voila's unload route killed the kernel instead of applying "
            "DELFIN's reload grace")
    finally:
        if kid and proc.poll() is None:
            try:
                delete = urllib.request.Request(
                    root + "/api/kernels/" + kid, method="DELETE", headers=auth)
                urllib.request.urlopen(delete, timeout=10).read()
            except Exception:
                pass
        if proc.poll() is None:
            proc.terminate()
            try:
                proc.wait(timeout=20)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait(timeout=10)
