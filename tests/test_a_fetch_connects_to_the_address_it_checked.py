"""web_fetch connects to the address it checked, not to a second answer.

The URL check resolved the host and refused private addresses; urllib then
resolved the host again to connect. A name that answers with a public
address the first time and 127.0.0.1 or 169.254.169.254 the second (DNS
rebinding) passed the check and reached the private service (audit
2026-09-16).
"""
import http.server
import socket
import threading

import pytest

from delfin.agent import web_tools as W


@pytest.fixture
def local_service():
    class H(http.server.BaseHTTPRequestHandler):
        def do_GET(self):
            body = b"INTERNAL-SERVICE-REACHED"
            self.send_response(200)
            self.send_header("Content-Type", "text/plain")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def log_message(self, *a):
            pass

    srv = http.server.HTTPServer(("127.0.0.1", 0), H)
    t = threading.Thread(target=srv.serve_forever, daemon=True)
    t.start()
    yield srv.server_address[1]
    srv.shutdown()


def test_a_name_that_rebinds_to_loopback_is_refused_at_connect(local_service, monkeypatch):
    real = socket.getaddrinfo
    answers = {"n": 0}

    def rebinding(host, port, *a, **kw):
        if host == "rebind.example":
            answers["n"] += 1
            ip = "93.184.216.34" if answers["n"] == 1 else "127.0.0.1"
            return [(socket.AF_INET, socket.SOCK_STREAM, 6, "", (ip, port or 0))]
        return real(host, port, *a, **kw)

    monkeypatch.setattr(socket, "getaddrinfo", rebinding)
    monkeypatch.delenv("http_proxy", raising=False)
    monkeypatch.delenv("HTTP_PROXY", raising=False)
    monkeypatch.setenv("no_proxy", "*")
    W._GUARDED_OPENER = None
    out = W.web_fetch(f"http://rebind.example:{local_service}/")
    assert "INTERNAL-SERVICE-REACHED" not in str(out)
    assert answers["n"] >= 2, "the connect did not resolve again; the test measured nothing"


def test_a_public_address_still_connects(monkeypatch):
    calls = []

    def fake_connect(self, sockaddr):
        calls.append(sockaddr)
        raise ConnectionRefusedError("test stops here")

    monkeypatch.setattr(socket.socket, "connect", fake_connect)
    monkeypatch.setattr(socket, "getaddrinfo", lambda h, p, *a, **k: [
        (socket.AF_INET, socket.SOCK_STREAM, 6, "", ("93.184.216.34", p))])
    with pytest.raises(ConnectionRefusedError):
        W._checked_create_connection(("example.org", 80), timeout=2)
    assert calls == [("93.184.216.34", 80)]


def test_a_mapped_loopback_address_is_refused():
    assert W._ip_refusal("::ffff:127.0.0.1", "x")
    assert W._ip_refusal("169.254.169.254", "x")
    assert W._ip_refusal("93.184.216.34", "x") is None
