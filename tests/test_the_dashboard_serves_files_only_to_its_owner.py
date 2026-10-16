"""What the dashboard server hands out over HTTP, measured on a real server.

One ``delfin.cli_voila`` server per module, rooted in a tmp directory that
holds harmless sample files: a PNG, an HTML page, an SVG, a PDF, a file in a
hidden directory, a link that points out of the root, and the drawing
editor's launch copy under ``.delfin/ketcher``. Every request is a plain
``http.client`` GET against ``127.0.0.1``; nothing outside ``tmp_path`` is
read or written, and no browser is started (``--no-browser``).

Properties asserted:

* ``/voila/files/`` answers only an authenticated request (it used to be a
  bare tornado ``StaticFileHandler`` that served every allow-listed file
  below the root to anybody who could reach the port).
* Active content (HTML, SVG) from user data is never rendered on the
  dashboard's origin: it comes back as a ``text/plain`` attachment under a
  ``sandbox`` CSP.
* A link out of the root and a hidden path are refused; PDFs and images
  stay viewable in place, and the editor's own copy still loads as HTML.
"""

from __future__ import annotations

import http.client
import socket
import subprocess
import sys
import time
from pathlib import Path

import pytest
from conftest import child_env

TOKEN = "files-route-test-token-0123456789abcdef0123"
PNG = b"\x89PNG\r\n\x1a\n" + b"\0" * 16


def _free_port() -> int:
    with socket.socket() as candidate:
        candidate.bind(("127.0.0.1", 0))
        return int(candidate.getsockname()[1])


@pytest.fixture(scope="module")
def server(tmp_path_factory):
    base = tmp_path_factory.mktemp("files-route")
    root = base / "root"
    runtime = base / "runtime"
    records = base / "records"
    outside = base / "outside"
    for directory in (root, runtime, records, outside):
        directory.mkdir()
        directory.chmod(0o700)
    (root / "pic.png").write_bytes(PNG)
    (root / "page.html").write_text("<p>sample</p>", encoding="utf-8")
    (root / "pic.svg").write_text(
        '<svg xmlns="http://www.w3.org/2000/svg"></svg>', encoding="utf-8")
    (root / "doc.pdf").write_bytes(b"%PDF-1.4\n%%EOF\n")
    (root / ".hidden").mkdir()
    (root / ".hidden" / "pic.png").write_bytes(PNG)
    (outside / "real.png").write_bytes(PNG)
    (root / "link.png").symlink_to(outside / "real.png")
    editor = root / ".delfin" / "ketcher"
    editor.mkdir(parents=True)
    (editor / "index.html").write_text("<html></html>", encoding="utf-8")

    port = _free_port()
    repo = Path(__file__).resolve().parents[1]
    env = child_env(base)
    for name in list(env):
        if name.startswith("VSCODE_") or name in ("BROWSER", "TERM_PROGRAM"):
            env.pop(name)
    env.update({
        "PYTHONPATH": str(repo),
        "XDG_RUNTIME_DIR": str(runtime),
        "DELFIN_SESSION_RECORD_DIR": str(records),
        "DELFIN_VOILA_ROOT_DIR": str(root),
        "DELFIN_DASHBOARD_STAY_UP": "1",
        "DELFIN_NO_BROWSER": "1",
    })
    log = base / "server.log"
    with log.open("wb") as output:
        proc = subprocess.Popen(
            [sys.executable, "-m", "delfin.cli_voila", "--port", str(port),
             "--ip", "127.0.0.1", "--token", TOKEN, "--no-browser"],
            cwd=str(repo), env=env, stdout=output, stderr=subprocess.STDOUT,
            start_new_session=True,
        )
    try:
        deadline = time.time() + 90
        while time.time() < deadline:
            if proc.poll() is not None:
                pytest.fail("server exited while starting:\n"
                            + log.read_text(errors="replace")[-2000:])
            try:
                status, _, _ = _get(port, "/api/status",
                                    {"Authorization": f"token {TOKEN}"})
                if status == 200:
                    break
            except OSError:
                pass
            time.sleep(0.3)
        else:
            pytest.fail("server did not answer:\n"
                        + log.read_text(errors="replace")[-2000:])
        yield {"port": port, "root": root, "log": log}
    finally:
        if proc.poll() is None:
            proc.terminate()
            try:
                proc.wait(timeout=20)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait()


def _get(port: int, path: str, headers: dict | None = None):
    conn = http.client.HTTPConnection("127.0.0.1", port, timeout=15)
    try:
        conn.request("GET", path, headers=headers or {})
        response = conn.getresponse()
        body = response.read()
        return response.status, {k.lower(): v for k, v in response.getheaders()}, body
    finally:
        conn.close()


def _owner(port: int, path: str):
    return _get(port, path, {"Authorization": f"token {TOKEN}"})


def test_a_file_is_not_served_without_the_token(server):
    status, headers, body = _get(server["port"], "/voila/files/pic.png")
    assert status in (302, 403), (status, headers)
    assert body != PNG


def test_the_owner_still_gets_an_image_inline(server):
    status, headers, body = _owner(server["port"], "/voila/files/pic.png")
    assert status == 200
    assert headers["content-type"].startswith("image/png")
    assert "attachment" not in headers.get("content-disposition", "")
    assert headers.get("x-content-type-options") == "nosniff"
    assert body == PNG


def test_a_pdf_stays_viewable_in_place(server):
    status, headers, _ = _owner(server["port"], "/voila/files/doc.pdf")
    assert status == 200
    assert headers["content-type"].startswith("application/pdf")
    assert "attachment" not in headers.get("content-disposition", "")
    assert headers.get("x-content-type-options") == "nosniff"


@pytest.mark.parametrize("name", ["page.html", "pic.svg"])
def test_active_content_from_user_data_is_a_sandboxed_download(server, name):
    status, headers, _ = _owner(server["port"], f"/voila/files/{name}")
    if status == 403:
        return        # refused outright by the allow-list: equally safe
    assert status == 200, status
    assert headers["content-type"].startswith("text/plain"), headers
    assert headers.get("content-disposition", "").startswith("attachment")
    assert "sandbox" in headers.get("content-security-policy", "")
    assert headers.get("x-content-type-options") == "nosniff"


def test_a_link_out_of_the_root_is_refused(server):
    status, _, body = _owner(server["port"], "/voila/files/link.png")
    assert status in (403, 404), status
    assert body != PNG


def test_a_hidden_path_is_refused(server):
    status, _, body = _owner(server["port"], "/voila/files/.hidden/pic.png")
    assert status in (403, 404), status
    assert body != PNG


def test_the_editor_launch_copy_still_loads_as_a_page(server):
    status, headers, _ = _owner(
        server["port"], "/voila/files/.delfin/ketcher/index.html")
    assert status == 200
    assert headers["content-type"].startswith("text/html")
    assert "attachment" not in headers.get("content-disposition", "")
    assert "sandbox" not in headers.get("content-security-policy", "")


# -- the root the server is pointed at ----------------------------------------

def _root_for(monkeypatch, tmp_path, value, **extra):
    from delfin import cli_voila

    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    monkeypatch.setenv("DELFIN_VOILA_ROOT_DIR", value)
    for key, val in extra.items():
        monkeypatch.setenv(key, val)
    return cli_voila._default_voila_root()


@pytest.mark.parametrize("broad", ["/", "/home", "HOME"])
def test_a_root_override_that_spans_the_machine_is_not_taken(
        monkeypatch, tmp_path, broad):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    value = str(home) if broad == "HOME" else broad
    root = _root_for(monkeypatch, tmp_path, value)
    assert root != Path(value).resolve()
    assert root.is_dir()


def test_a_broad_root_is_taken_when_asked_for_explicitly(monkeypatch, tmp_path):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    root = _root_for(monkeypatch, tmp_path, str(home),
                     DELFIN_VOILA_ALLOW_BROAD_ROOT="1")
    assert root == home.resolve()


def test_a_root_override_that_is_not_a_directory_is_not_taken(
        monkeypatch, tmp_path):
    missing = tmp_path / "missing"
    a_file = tmp_path / "a_file"
    a_file.write_text("x", encoding="utf-8")
    assert _root_for(monkeypatch, tmp_path, str(missing)) != missing.resolve()
    assert _root_for(monkeypatch, tmp_path, str(a_file)) != a_file.resolve()


def test_a_project_directory_override_is_taken(monkeypatch, tmp_path):
    project = tmp_path / "project"
    project.mkdir()
    assert _root_for(monkeypatch, tmp_path, str(project)) == project.resolve()


# -- the handler itself, with the allow-list widened ----------------------------
#
# The server above refuses HTML/SVG at the allow-list already. The handler's
# own rule -- user data that is not an image or a PDF goes out as a sandboxed
# text/plain attachment -- is what still holds when a config widens the
# allow-list, so it is driven here directly with ``.*``. Authentication is
# covered by the live server; this subclass skips only Jupyter's request
# preparation and the auth decorators.

def test_the_handler_sends_user_html_as_a_sandboxed_download(tmp_path):
    import asyncio

    from tornado import web
    from tornado.httpclient import AsyncHTTPClient
    from tornado.httpserver import HTTPServer
    from tornado.netutil import bind_sockets

    from delfin.dashboard import web_security

    root = tmp_path / "root"
    root.mkdir()
    (root / "page.html").write_text("<p>sample</p>", encoding="utf-8")
    (root / "pic.svg").write_text("<svg/>", encoding="utf-8")
    (root / "pic.png").write_bytes(PNG)
    (root / "doc.pdf").write_bytes(b"%PDF-1.4\n%%EOF\n")
    editor = root / ".delfin" / "ketcher"
    editor.mkdir(parents=True)
    (editor / "main.js").write_text("void 0;", encoding="utf-8")

    base = web_security.handler_class()

    class Unauthenticated(base):
        async def prepare(self):
            return None

        def get(self, path, include_body=True):
            return web.StaticFileHandler.get(self, path, include_body)

    async def scenario():
        app = web.Application(
            [(r"/voila/files/(.*)", Unauthenticated,
              {"path": str(root), "allowlist": [".*"], "denylist": []})],
            base_url="/", headers=dict(web_security.SECURITY_HEADERS))
        sockets = bind_sockets(0, "127.0.0.1")
        port = sockets[0].getsockname()[1]
        server = HTTPServer(app)
        server.add_sockets(sockets)
        client = AsyncHTTPClient()
        out = {}
        try:
            for name in ("page.html", "pic.svg", "pic.png", "doc.pdf",
                         ".delfin/ketcher/main.js"):
                response = await client.fetch(
                    f"http://127.0.0.1:{port}/voila/files/{name}",
                    raise_error=False)
                out[name] = response
        finally:
            server.stop()
        return out

    got = asyncio.run(scenario())
    for name in ("page.html", "pic.svg"):
        response = got[name]
        assert response.code == 200, (name, response.code)
        assert response.headers["Content-Type"].startswith("text/plain")
        assert response.headers["Content-Disposition"].startswith("attachment")
        assert "sandbox" in response.headers["Content-Security-Policy"]
        assert response.headers["X-Content-Type-Options"] == "nosniff"
    assert got["pic.png"].headers["Content-Type"] == "image/png"
    assert "Content-Disposition" not in got["pic.png"].headers
    assert got["doc.pdf"].headers["Content-Type"] == "application/pdf"
    assert "sandbox" not in got["doc.pdf"].headers["Content-Security-Policy"]
    script = got[".delfin/ketcher/main.js"]
    assert script.code == 200
    assert "javascript" in script.headers["Content-Type"]
    assert "Content-Disposition" not in script.headers


def test_a_failed_handler_build_leaves_the_route_refusing(monkeypatch):
    from tornado import web
    from voila.static_file_handler import AllowListFileHandler

    from delfin.dashboard import web_security

    app = web.Application([(r"/voila/files/(.*)", AllowListFileHandler,
                            {"path": "/nonexistent", "allowlist": [".*"]})])

    def broken():
        raise ImportError("handler unavailable")

    monkeypatch.setattr(web_security, "handler_class", broken)
    assert web_security.install_files_route(app) == 1
    targets = [rule.target for rule in web_security._rules(app.wildcard_router)]
    assert AllowListFileHandler not in targets
    assert any(t.__name__ == "RefuseFiles" for t in targets if isinstance(t, type))
