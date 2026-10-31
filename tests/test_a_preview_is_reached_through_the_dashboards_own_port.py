"""A preview is reached through the main dashboard's port, never its own.

Field report: through VS Code port forwarding or the Windows launcher's one
SSH tunnel only the main dashboard's port reaches the browser, so a link to
``http://localhost:<preview port>/...`` was refused. The main server now
proxies ``<base>/delfin-preview/<slug>/...`` (HTTP and websocket) to the
preview, which serves under that same base URL on 127.0.0.1.

Input: a jupyter server started here with only the
``delfin.dashboard.preview_proxy`` extension and a token; a stand-in preview
(a tornado app on 127.0.0.1 that answers only ``Authorization: token
<preview token>``, echoes what it received, and echoes websocket messages);
its record in ``~/.delfin/branch_previews`` under a temporary HOME.

Output: without the main server's authentication the proxy answers neither
HTTP nor websocket and the preview sees nothing; with it, the request
reaches the preview carrying the preview's token (added by the proxy), the
main token and the main login cookie are not passed on; another slug, a
path that climbs out with ``..`` and a closed preview answer 404; the link
the dashboard shows is a path on the page's own origin with no port and no
token.

Universal: skipped where jupyter_server or tornado is not installed.
"""

from __future__ import annotations

import json
import os
import socket
import subprocess
import sys
import textwrap
import time
import urllib.error
import urllib.request
from pathlib import Path

import pytest

from conftest import child_env

pytest.importorskip("jupyter_server")
pytest.importorskip("tornado")

import delfin  # noqa: E402
from delfin.dashboard import branch_preview as BP  # noqa: E402

ROOT = Path(delfin.__file__).resolve().parent.parent
MAIN_TOKEN = "m" * 43
PREV_TOKEN = "p" * 43

_UPSTREAM = textwrap.dedent('''
    import json, sys
    from tornado import ioloop, web, websocket
    port, token, log = int(sys.argv[1]), sys.argv[2], sys.argv[3]

    def note(kind, h):
        with open(log, "a") as fh:
            fh.write(json.dumps({"kind": kind, "path": h.request.path,
                                 "query": h.request.query,
                                 "auth": h.request.headers.get("Authorization", ""),
                                 "cookie": h.request.headers.get("Cookie", "")}) + "\\n")

    class Echo(web.RequestHandler):
        def get(self, _p):
            note("http", self)
            if self.request.headers.get("Authorization") != "token " + token:
                raise web.HTTPError(403)
            self.set_header("Content-Type", "application/json")
            self.write({"path": self.request.path, "query": self.request.query})

        post = get

    class WS(websocket.WebSocketHandler):
        def check_origin(self, origin):
            return True
        def get(self, *a):
            note("ws", self)
            if self.request.headers.get("Authorization") != "token " + token:
                raise web.HTTPError(403)
            return super().get(*a)
        def on_message(self, m):
            self.write_message("echo:" + m)

    app = web.Application([(r"(.*/channels)", WS), (r"(.*)", Echo)])
    app.listen(port, address="127.0.0.1")
    ioloop.IOLoop.current().start()
''')


def _free_port() -> int:
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def _wait_port(port: int, proc, timeout=60.0) -> None:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        with socket.socket() as s:
            if s.connect_ex(("127.0.0.1", port)) == 0:
                return
        assert proc.poll() is None, "server exited"
        time.sleep(0.2)
    raise AssertionError(f"port {port} never opened")


@pytest.fixture
def servers(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    # child_env: no provider credential reaches the servers, and their own
    # ~/.delfin sinks resolve outside the real home; HOME is then set to
    # the same private home the test process uses.
    env = {k: v for k, v in child_env(tmp_path).items()
           if not k.startswith(("VSCODE_", "JPY_", "JUPYTER_"))
           and k not in ("BROWSER", "PYTHONPATH")}
    env.update({"HOME": str(home), "PYTHONPATH": str(ROOT),
                "JUPYTER_TOKEN": MAIN_TOKEN, "DELFIN_NO_BROWSER": "1",
                "JUPYTER_CONFIG_DIR": str(tmp_path / "jcfg"),
                "JUPYTER_DATA_DIR": str(tmp_path / "jdata"),
                "JUPYTER_RUNTIME_DIR": str(tmp_path / "jrun")})
    up_port, main_port = _free_port(), _free_port()
    log = tmp_path / "upstream.jsonl"
    script = tmp_path / "upstream.py"
    script.write_text(_UPSTREAM)
    up = subprocess.Popen([sys.executable, str(script), str(up_port), PREV_TOKEN,
                           str(log)], env=env, start_new_session=True,
                          stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    main = subprocess.Popen(
        [sys.executable, "-m", "jupyter", "server", "--no-browser",
         f"--port={main_port}", "--ServerApp.port_retries=0",
         "--ServerApp.ip=127.0.0.1", f"--ServerApp.root_dir={tmp_path}",
         "--ServerApp.jpserver_extensions={'delfin.dashboard.preview_proxy': True}"],
        env=env, start_new_session=True, stdout=open(tmp_path / "main.log", "w"),
        stderr=subprocess.STDOUT)
    try:
        _wait_port(up_port, up)
        _wait_port(main_port, main)
        from delfin.agent import proc_identity
        p = BP.Preview(
            slug="feat", ref="feat", sha="0" * 40, branch="feat", repo=str(tmp_path),
            worktree=str(tmp_path / "wt"), port=up_port, token=PREV_TOKEN,
            url_path="/delfin-preview/feat/voila/render/x.ipynb", pid=up.pid,
            proc_start=proc_identity.process_start(up.pid), ticks=None,
            host=socket.gethostname(), started_at=time.time(), delfin_file="",
            log="", base_url="/delfin-preview/feat/")
        BP._save(p)
        yield {"main": main_port, "up": up_port, "log": log, "preview": p,
               "tmp": tmp_path}
    finally:
        for proc in (main, up):
            try:
                os.killpg(proc.pid, 9)
            except OSError:
                pass
            proc.wait(timeout=10)


def _get(port, path, *, auth=True, headers=None, method="GET"):
    h = dict(headers or {})
    if auth:
        h["Authorization"] = f"token {MAIN_TOKEN}"
    req = urllib.request.Request(f"http://127.0.0.1:{port}{path}", headers=h, method=method)

    class NoRedirect(urllib.request.HTTPRedirectHandler):
        def redirect_request(self, *a, **k):
            return None
    opener = urllib.request.build_opener(NoRedirect)
    try:
        r = opener.open(req, timeout=20)
        return r.status, r.read().decode(), dict(r.headers)
    except urllib.error.HTTPError as e:
        return e.code, e.read().decode(errors="replace"), dict(e.headers)


def _seen(log: Path) -> list:
    if not log.exists():
        return []
    return [json.loads(x) for x in log.read_text().splitlines()]


def test_without_the_main_login_nothing_reaches_the_preview(servers):
    status, _body, _h = _get(servers["main"], "/delfin-preview/feat/voila/x", auth=False)
    assert status in (302, 401, 403), status
    assert _seen(servers["log"]) == []


def test_an_authenticated_request_reaches_the_preview_with_its_own_token(servers):
    status, body, headers = _get(
        servers["main"], f"/delfin-preview/feat/voila/x?token={MAIN_TOKEN}&a=1",
        headers={"Cookie": f"username-127-0-0-1-{servers['main']}=secret; _xsrf=q1"})
    assert status == 200, body
    seen = _seen(servers["log"])[-1]
    assert seen["path"] == "/delfin-preview/feat/voila/x"
    assert seen["auth"] == f"token {PREV_TOKEN}"          # added by the proxy
    assert MAIN_TOKEN not in json.dumps(seen)             # main token stays here
    assert "a=1" in seen["query"]
    assert "username-" not in seen["cookie"] and "_xsrf=q1" in seen["cookie"]
    assert PREV_TOKEN not in json.dumps(headers)


@pytest.mark.parametrize("path", [
    "/delfin-preview/other/voila/x",
    "/delfin-preview/feat/../other/x",
    "/delfin-preview/feat/%2e%2e/%2e%2e/api/contents",
    "/delfin-preview/..%2ffeat/x",
])
def test_another_slug_or_a_climbing_path_is_not_found(servers, path):
    before = len(_seen(servers["log"]))
    status, _b, _h = _get(servers["main"], path)
    assert status in (400, 404), (path, status)
    assert len(_seen(servers["log"])) == before


def test_a_closed_preview_is_not_found(servers):
    BP._record_path("feat").unlink()
    status, _b, _h = _get(servers["main"], "/delfin-preview/feat/voila/x")
    assert status == 404


def test_the_websocket_goes_through_only_when_logged_in(servers):
    import asyncio
    from tornado import httpclient, websocket

    async def talk(auth):
        h = {"Authorization": f"token {MAIN_TOKEN}"} if auth else {}
        req = httpclient.HTTPRequest(
            f"ws://127.0.0.1:{servers['main']}/delfin-preview/feat/api/kernels/k1/channels",
            headers=h)
        conn = await websocket.websocket_connect(req)
        await conn.write_message("hello")
        msg = await conn.read_message()
        conn.close()
        return msg

    assert asyncio.run(talk(True)) == "echo:hello"
    with pytest.raises(Exception):
        asyncio.run(talk(False))
    ws_seen = [s for s in _seen(servers["log"]) if s["kind"] == "ws"]
    assert len(ws_seen) == 1 and ws_seen[0]["auth"] == f"token {PREV_TOKEN}"


def test_the_link_is_a_path_on_the_pages_origin(monkeypatch):
    from delfin.dashboard.branch_panel import preview_link_html
    monkeypatch.setenv(BP.PROXY_ENV, "/")
    p = BP.Preview(
        slug="feat", ref="feat", sha="0" * 40, branch="feat", repo="/r", worktree="/w",
        port=45999, token=PREV_TOKEN, url_path="/delfin-preview/feat/voila/render/x.ipynb",
        pid=1, proc_start="", ticks=None, host="h", started_at=0.0, delfin_file="",
        log="", base_url="/delfin-preview/feat/")
    link = preview_link_html(p)
    assert 'href="/delfin-preview/feat/voila/render/x.ipynb"' in link
    assert PREV_TOKEN not in link and "45999" not in link and "localhost" not in link


def test_a_preview_is_started_under_the_proxy_base_url(tmp_path, monkeypatch):
    """The base URL reaches the preview server through a config file on
    JUPYTER_CONFIG_PATH, so a branch whose launcher has no such flag is
    served under it too."""
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv(BP.PROXY_ENV, "/")
    env = BP.proxy_env({}, "feat")
    cfg_dirs = env["JUPYTER_CONFIG_PATH"].split(os.pathsep)
    cfg = json.loads((Path(cfg_dirs[0]) / "jupyter_server_config.json").read_text())
    assert cfg["ServerApp"]["base_url"] == "/delfin-preview/feat/"
    monkeypatch.delenv(BP.PROXY_ENV)
    assert BP.proxy_env({"A": "1"}, "feat") == {"A": "1"}


@pytest.mark.parametrize("auth, expected", [(True, 200), (False, 403)])
def test_post_authentication_is_ready_before_xsrf_check(servers, auth, expected):
    status, body, _ = _get(servers["main"],
                          "/delfin-preview/feat/voila/api/shutdown/kernel",
                          auth=auth, method="POST")
    assert status == expected, body
    assert bool(_seen(servers["log"])) == auth


def test_cookie_login_still_requires_a_matching_xsrf_token(servers):
    import http.cookiejar
    jar = http.cookiejar.CookieJar()
    opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(jar))
    with opener.open(f"http://127.0.0.1:{servers['main']}/api?token={MAIN_TOKEN}",
                     timeout=20) as response:
        assert response.status == 200
    login = "; ".join(f"{c.name}={c.value}" for c in jar)
    assert "username-" in login
    path = "/delfin-preview/feat/delfin-api/window-closing"
    token = "a" * 32
    headers = {"Cookie": login + "; _xsrf=" + token}
    for sent, expected in ((None, 403), ("b" * 32, 403), (token, 200)):
        if sent:
            headers["X-Xsrftoken"] = sent
        status, body, _ = _get(servers["main"], path, auth=False,
                              headers=headers, method="POST")
        assert status == expected, body
    assert len(_seen(servers["log"])) == 1


def test_an_exited_preview_is_not_proxied_while_waiting_to_be_reaped(servers):
    import signal
    pid = servers["preview"].pid
    stat_path = Path(f"/proc/{pid}/stat")
    if not stat_path.exists():
        pytest.skip("Linux process state required")
    os.kill(pid, signal.SIGTERM)
    deadline = time.monotonic() + 5
    while time.monotonic() < deadline:
        stat = stat_path.read_text()
        if stat[stat.rindex(")") + 1:].split()[0] == "Z":
            break
        time.sleep(.01)
    else:
        pytest.fail("preview never exited")
    status, body, _ = _get(servers["main"], "/delfin-preview/feat/api/kernels")
    assert status == 404, body
