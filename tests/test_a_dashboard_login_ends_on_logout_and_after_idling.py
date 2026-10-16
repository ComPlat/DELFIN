"""A dashboard login ends when its owner logs out or stops using it.

Input: the login cookie jupyter-server sets for a token login. Output: a
request that carries the cookie is authenticated, or it is not.

Properties asserted:

* logging out (``GET /logout``) ends the login on the server, not only in
  the browser: the same cookie value, sent again, no longer authenticates.
  jupyter-server's cookie is a signed, self-contained value, so without a
  server-side session a copied cookie stayed valid for its whole lifetime
  after the logout;
* a login nobody uses for longer than the idle timeout
  (``DELFIN_DASHBOARD_IDLE_HOURS``, default 12 h) ends; one in use does
  not; the absolute lifetime (``DELFIN_DASHBOARD_LOGIN_DAYS``) still holds;
* the session store is owner-only, holds no session id in clear, and
  survives a server restart (a browser that was logged in stays logged
  in, a logged-out cookie stays refused);
* a token login cannot set a permanent password (jupyter-server's
  ``allow_password_change``);
* the page has a logout control and reports user activity, so an open
  page that is being used does not idle out.

Server-level tests run one ``delfin.cli_voila`` server (``--no-browser``,
HOME in tmp_path); unit tests drive the session store with a fake clock.
"""

from __future__ import annotations

import http.client
import json
import os
import socket
import stat
import subprocess
import sys
import time
import urllib.parse
from pathlib import Path

import pytest
from conftest import child_env

from delfin.dashboard import login_session

TOKEN = "logout-test-token-0123456789abcdef0123456789"
REPO = Path(__file__).resolve().parents[1]


# -- the session store, with a fake clock ---------------------------------------

class Clock:
    def __init__(self):
        self.now = 1_000_000.0

    def __call__(self):
        return self.now


def _store(tmp_path, clock, idle=3600.0, lifetime=7 * 86400.0):
    return login_session.Sessions(tmp_path / "sessions.json", idle=idle,
                                  lifetime=lifetime, clock=clock)


def test_a_login_in_use_stays_and_an_idle_one_ends(tmp_path):
    clock = Clock()
    store = _store(tmp_path, clock)
    sid = store.open()
    for _ in range(5):               # used every 50 minutes for 4 hours
        clock.now += 3000
        assert store.touch(sid)
    clock.now += 3601                # then left alone for an hour
    assert not store.touch(sid)
    assert not store.touch(sid)      # and it stays ended


def test_logging_out_ends_the_login(tmp_path):
    store = _store(tmp_path, Clock())
    sid = store.open()
    other = store.open()
    store.close(sid)
    assert not store.touch(sid)
    assert store.touch(other)


def test_the_absolute_lifetime_holds_for_a_login_in_use(tmp_path):
    clock = Clock()
    store = _store(tmp_path, clock, idle=3600.0, lifetime=10_000.0)
    sid = store.open()
    while clock.now < 1_000_000.0 + 9_000:
        clock.now += 1000
        assert store.touch(sid)
    clock.now += 2000
    assert not store.touch(sid)


def test_an_unknown_or_malformed_session_is_refused(tmp_path):
    store = _store(tmp_path, Clock())
    for sid in ("", "x", None, 5, "a" * 10_000):
        assert not store.touch(sid)


def test_the_store_is_owner_only_and_holds_no_session_id(tmp_path):
    store = _store(tmp_path, Clock())
    sid = store.open()
    path = tmp_path / "sessions.json"
    assert stat.S_IMODE(path.stat().st_mode) == 0o600
    assert sid not in path.read_text(encoding="utf-8")


def test_a_restart_keeps_a_login_and_keeps_a_logout(tmp_path):
    clock = Clock()
    first = _store(tmp_path, clock)
    kept, ended = first.open(), first.open()
    first.close(ended)
    again = _store(tmp_path, clock)
    assert again.touch(kept)
    assert not again.touch(ended)


def test_the_store_does_not_grow_without_bound(tmp_path):
    clock = Clock()
    store = _store(tmp_path, clock)
    for _ in range(login_session.MAX_SESSIONS + 50):
        store.open()
        clock.now += 1
    assert len(store) <= login_session.MAX_SESSIONS


def test_one_off_logins_do_not_push_out_a_browser_that_is_logged_in(tmp_path):
    clock = Clock()
    store = _store(tmp_path, clock)
    browser = store.open()
    clock.now += 10
    assert store.touch(browser)
    for _ in range(login_session.MAX_SESSIONS * 2):   # a script, no cookie jar
        clock.now += 1
        store.open()
    assert store.touch(browser)


def test_the_idle_timeout_is_bounded(monkeypatch):
    monkeypatch.delenv(login_session.IDLE_ENV, raising=False)
    assert login_session.idle_seconds() == 12 * 3600
    monkeypatch.setenv(login_session.IDLE_ENV, "2")
    assert login_session.idle_seconds() == 2 * 3600
    monkeypatch.setenv(login_session.IDLE_ENV, "0.001")
    assert login_session.idle_seconds() == 15 * 60
    monkeypatch.setenv(login_session.IDLE_ENV, "100000")
    from delfin.dashboard import web_security
    assert login_session.idle_seconds() == web_security.login_days() * 86400
    monkeypatch.setenv(login_session.IDLE_ENV, "not a number")
    assert login_session.idle_seconds() == 12 * 3600


# -- the page --------------------------------------------------------------------

def test_the_header_has_a_logout_control_and_the_page_reports_activity():
    html = login_session.logout_control_html()
    assert "logout" in html and "Log out" in html
    js = login_session.keepalive_js()
    assert login_session.ALIVE_ROUTE.lstrip("/") in js
    for event in ("pointerdown", "keydown"):
        assert event in js
    dashboard = (REPO / "delfin" / "dashboard" / "__init__.py").read_text(
        encoding="utf-8")
    assert "logout_control_html()" in dashboard
    assert "keepalive_js()" in dashboard


# -- a real server ---------------------------------------------------------------

def _free_port() -> int:
    with socket.socket() as candidate:
        candidate.bind(("127.0.0.1", 0))
        return int(candidate.getsockname()[1])


def _request(port, method, path, headers=None, body=None):
    conn = http.client.HTTPConnection("127.0.0.1", port, timeout=15)
    try:
        conn.request(method, path, body=body, headers=headers or {})
        response = conn.getresponse()
        data = response.read()
        return (response.status,
                [(k.lower(), v) for k, v in response.getheaders()], data)
    finally:
        conn.close()


def _set_cookies(headers) -> dict:
    jar = {}
    for name, value in headers:
        if name == "set-cookie":
            pair = value.split(";", 1)[0]
            key, _, val = pair.partition("=")
            jar[key.strip()] = val.strip()
    return jar


def _start(base: Path, port: int):
    runtime = base / "runtime"
    root = base / "root"
    for directory in (runtime, root):
        directory.mkdir(exist_ok=True)
        directory.chmod(0o700)
    env = child_env(base)
    for name in list(env):
        if name.startswith("VSCODE_") or name in (
                "BROWSER", "TERM_PROGRAM", "DISPLAY", "WAYLAND_DISPLAY",
                "DBUS_SESSION_BUS_ADDRESS"):
            env.pop(name)
    env.update({
        "PYTHONPATH": str(REPO),
        "XDG_RUNTIME_DIR": str(runtime),
        "DELFIN_VOILA_ROOT_DIR": str(root),
        "DELFIN_DASHBOARD_STAY_UP": "1",
        "DELFIN_NO_BROWSER": "1",
        "JUPYTER_CONFIG_DIR": str(base / "jupyter-config"),
    })
    log = base / f"server-{port}.log"
    output = log.open("ab")
    proc = subprocess.Popen(
        [sys.executable, "-m", "delfin.cli_voila", "--port", str(port),
         "--ip", "127.0.0.1", "--token", TOKEN, "--no-browser",
         "--strict-port"],
        cwd=str(REPO), env=env, stdout=output, stderr=subprocess.STDOUT,
        start_new_session=True)
    deadline = time.time() + 90
    while time.time() < deadline:
        if proc.poll() is not None:
            pytest.fail("server exited:\n" + log.read_text(errors="replace")[-2000:])
        try:
            status, _, _ = _request(port, "GET", "/api/status",
                                    {"Authorization": f"token {TOKEN}"})
            if status == 200:
                return proc, log
        except OSError:
            pass
        time.sleep(0.3)
    proc.terminate()
    pytest.fail("server did not answer:\n" + log.read_text(errors="replace")[-2000:])


def _stop(proc):
    if proc.poll() is None:
        proc.terminate()
        try:
            proc.wait(timeout=20)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait()


@pytest.fixture(scope="module")
def server(tmp_path_factory):
    base = tmp_path_factory.mktemp("login-session")
    port = _free_port()
    proc, log = _start(base, port)
    state = {"base": base, "port": port, "proc": proc, "log": log}
    try:
        yield state
    finally:
        _stop(state["proc"])


def _login(port) -> str:
    """A token login; returns the ``Cookie`` header it earned."""
    status, headers, _ = _request(port, "GET", f"/api/status?token={TOKEN}")
    assert status == 200, status
    jar = _set_cookies(headers)
    login = {k: v for k, v in jar.items() if k.startswith("username-")}
    assert login, headers
    return "; ".join(f"{k}={v}" for k, v in login.items())


def _status_with(port, cookie) -> int:
    status, _, _ = _request(port, "GET", "/api/status", {"Cookie": cookie})
    return status


def test_logging_out_ends_the_login_on_the_server(server):
    port = server["port"]
    cookie = _login(port)
    assert _status_with(port, cookie) == 200
    status, _, _ = _request(port, "GET", "/logout", {"Cookie": cookie})
    assert status == 200, status
    # The browser has dropped the cookie; a copy of it must not work either.
    assert _status_with(port, cookie) == 403


def test_another_login_survives_a_logout(server):
    port = server["port"]
    mine, other = _login(port), _login(port)
    _request(port, "GET", "/logout", {"Cookie": mine})
    assert _status_with(port, other) == 200


def test_the_page_reports_activity_on_its_own_route(server):
    port = server["port"]
    cookie = _login(port)
    status, _, body = _request(port, "GET", login_session.ALIVE_ROUTE,
                               {"Cookie": cookie})
    assert status == 204, (status, body)
    status, _, _ = _request(port, "GET", login_session.ALIVE_ROUTE)
    assert status in (302, 403), status


def test_a_token_login_cannot_set_a_permanent_password(server):
    port = server["port"]
    status, headers, _ = _request(port, "GET", "/login")
    xsrf = _set_cookies(headers).get("_xsrf", "")
    body = urllib.parse.urlencode({"password": TOKEN,
                                   "new_password": "a-password-of-my-own",
                                   "_xsrf": xsrf})
    _request(port, "POST", "/login", {
        "Content-Type": "application/x-www-form-urlencoded",
        "Cookie": f"_xsrf={xsrf}"}, body)
    config = server["base"] / "jupyter-config" / "jupyter_server_config.json"
    if config.exists():
        data = json.loads(config.read_text(encoding="utf-8") or "{}")
        assert not (data.get("IdentityProvider") or {}).get("hashed_password")
        assert not (data.get("ServerApp") or {}).get("password")


def test_a_restart_keeps_a_login_and_keeps_a_logout(server):
    port = server["port"]
    kept, ended = _login(port), _login(port)
    _request(port, "GET", "/logout", {"Cookie": ended})
    _stop(server["proc"])
    server["proc"], server["log"] = _start(server["base"], port)
    assert _status_with(port, kept) == 200
    assert _status_with(port, ended) == 403


def test_the_session_store_on_the_host_is_owner_only(server):
    _login(server["port"])
    # jupyter-server's runtime directory, below the child's HOME.
    stores = list(server["base"].rglob(login_session.STORE_PREFIX + "*"))
    assert stores, sorted(str(p) for p in server["base"].rglob("*runtime*"))
    for path in stores:
        assert stat.S_IMODE(path.stat().st_mode) == 0o600, path
        assert os.stat(path).st_uid == os.getuid()
