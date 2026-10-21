"""Shared pieces of the egress proxy tests.

Used by test_a_child_process_reaches_the_network_only_through_the_proxy.py
(the proxy's decisions, host-independent) and by
test_a_sandboxed_command_reaches_the_network_through_its_proxy.py (the
kernel layers, which need a host that has them). No ``test_`` prefix: this
module holds no tests and is not scanned by the machine-binding ratchet.

``serve_web`` and ``session_proxy`` are generators; each test file wraps
them in its own fixture so pyflakes sees every name it imports used.
"""
from __future__ import annotations

import base64
import hashlib
import http.server
import io
import socket
import threading
import zipfile
from pathlib import Path

from delfin.agent import api_client as A
from delfin.agent import egress_proxy as EP

HELPER = str(Path(A.__file__).with_name("landlock_exec.py"))


def wheel() -> tuple[str, bytes]:
    """A minimal valid wheel ``tinypkg`` 0.1 that pip can install."""
    name = "tinypkg-0.1-py3-none-any.whl"
    files = {
        "tinypkg/__init__.py": b"VALUE = 42\n",
        "tinypkg-0.1.dist-info/METADATA":
            b"Metadata-Version: 2.1\nName: tinypkg\nVersion: 0.1\n",
        "tinypkg-0.1.dist-info/WHEEL":
            b"Wheel-Version: 1.0\nGenerator: test\nRoot-Is-Purelib: true\n"
            b"Tag: py3-none-any\n",
    }
    record = []
    for path, data in files.items():
        digest = base64.urlsafe_b64encode(
            hashlib.sha256(data).digest()).rstrip(b"=").decode()
        record.append(f"{path},sha256={digest},{len(data)}")
    record.append("tinypkg-0.1.dist-info/RECORD,,")
    files["tinypkg-0.1.dist-info/RECORD"] = ("\n".join(record) + "\n").encode()
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        for path, data in files.items():
            zf.writestr(path, data)
    return name, buf.getvalue()


def serve_web():
    """A local web on 127.0.0.1: a page, an upload sink, a package index.

    Yields the server; ``server.seen`` lists (method, Host, path)."""
    wheel_name, wheel_bytes = wheel()
    seen = []

    class H(http.server.BaseHTTPRequestHandler):
        protocol_version = "HTTP/1.0"

        def _send(self, body: bytes, ctype="text/plain"):
            self.send_response(200)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def do_GET(self):
            seen.append(("GET", self.headers.get("Host"), self.path))
            if self.path.rstrip("/") == "/simple/tinypkg":
                self._send(f'<a href="/files/{wheel_name}">{wheel_name}</a>'.encode(),
                           "text/html")
            elif self.path == f"/files/{wheel_name}":
                self._send(wheel_bytes, "application/octet-stream")
            else:
                self._send(f"hello {self.headers.get('Host')}".encode())

        def do_POST(self):
            n = int(self.headers.get("Content-Length") or 0)
            self.rfile.read(n)
            seen.append(("POST", self.headers.get("Host"), self.path))
            self._send(f"got {n}".encode())

        def log_message(self, *a):
            pass

    srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), H)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    srv.seen = seen
    yield srv
    srv.shutdown()


def session_proxy(web, monkeypatch):
    """The process-wide egress proxy in proxy mode, allowing ``allowed.test``
    on the web's port; every host name connects to 127.0.0.1.

    Yields the proxy; ``proxy.refused`` lists its refusals."""
    refused = []
    p = EP.EgressProxy(
        allowed_domains=["allowed.test"], allowed_ports=[web.server_address[1]],
        connector=lambda h, port, t: socket.create_connection(("127.0.0.1", port), timeout=t),
        on_refusal=refused.append, upload_ask_bytes=1024 * 1024)
    p.refused = refused
    monkeypatch.setattr(EP, "_PROXY", p)
    monkeypatch.setattr(EP, "_settings", lambda: {"mode": "proxy"})
    monkeypatch.setattr(EP, "_MODE_AT_START", None)
    monkeypatch.setattr(EP, "_ENDPOINT_HOSTS", set())
    yield p
    p.close()


def perms_in(tmp_path, mode="default", asker=None):
    """Session permissions with the workspace ``tmp_path/ws``."""
    ws = tmp_path / "ws"
    ws.mkdir(exist_ok=True)
    return A.KitToolPermissions(workspace=str(ws), mode=mode,
                                confirm_callback=asker)
