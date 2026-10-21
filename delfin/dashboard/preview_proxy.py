"""Serve each branch preview through the main dashboard's own port.

Input: a request to ``<base>/delfin-preview/<slug>/<rest>`` on the main
dashboard server (HTTP, or a websocket for ``.../channels``), and the
preview records in ``~/.delfin/branch_previews`` (``branch_preview``).
Output: the same request sent to ``http://127.0.0.1:<preview port>`` --
the preview serves under that same base URL, so the path passes through
unchanged -- and its answer streamed back; websocket frames are relayed
both ways.

Why: through VS Code port forwarding or the Windows launcher's single SSH
tunnel only the main port reaches the browser; a link to the preview's own
port was refused (field report). The preview keeps binding to 127.0.0.1
and keeps its own token.

Rules:

* Authentication is the main server's: ``@web.authenticated`` for HTTP,
  ``current_user`` before the websocket handshake, and the main server's
  websocket origin check. A state-changing request authenticated by the
  login cookie must carry an XSRF token equal to one of the ``_xsrf``
  cookies sent (the preview's page sends its own, set under the preview's
  path; tornado's own check reads only one of the two cookies of that
  name). The preview's token is added server-side as ``Authorization:
  token ...``; it never reaches the browser.
* Not passed on: the main ``token`` query argument, the main login cookie
  (only ``_xsrf`` and the preview's own ``username-*-<port>`` cookie go
  through), ``Host``, ``Origin``, ``Referer`` and hop-by-hop headers.
* 404 for: a slug that is not ``[A-Za-z0-9][A-Za-z0-9._-]*``, a record
  that is missing, from another host, not started for this proxy (its
  ``base_url`` differs) or whose server is gone; a path with a ``.`` or
  ``..`` segment or a percent-encoded dot, slash or backslash.

Rejected: ``jupyter-server-proxy`` (a new dependency for one route, and its
generic ``/proxy/<port>/`` route reaches any local port); asking people to
forward a second port (impossible through the Windows launcher).
"""

from __future__ import annotations

import hmac
import os
import re
import socket
from typing import Optional
from urllib.parse import parse_qsl, urlencode

PREFIX = "delfin-preview"
SLUG_RE = r"[A-Za-z0-9][A-Za-z0-9._-]{0,79}"
_BAD_REST = re.compile(r"%2e|%2f|%5c|\\", re.I)
_HOP = {"connection", "keep-alive", "proxy-authenticate", "proxy-authorization",
        "te", "trailers", "transfer-encoding", "upgrade", "host", "origin",
        "referer", "authorization", "cookie", "content-length"}
_RESP_DROP = {"connection", "keep-alive", "transfer-encoding", "proxy-authenticate",
              "trailers", "upgrade"}


def base_prefix(base_url: str, slug: str) -> str:
    """``<base>/delfin-preview/<slug>/``."""
    base = "/" + (base_url or "/").strip("/")
    return (base.rstrip("/") + f"/{PREFIX}/{slug}/")


def safe_rest(path: str) -> bool:
    if _BAD_REST.search(path or ""):
        return False
    return not any(seg in (".", "..") for seg in (path or "").split("/"))


def target(slug: str, path: str, base_url: str = "/"):
    """The preview record the request may reach, or None."""
    from delfin.dashboard import branch_preview as BP
    if not re.fullmatch(SLUG_RE, slug or "") or not safe_rest(path):
        return None
    prefix = base_prefix(base_url, slug)
    if not (path or "").startswith(prefix):
        return None
    p = BP.load(slug)
    if p is None or p.host != socket.gethostname():
        return None
    if (getattr(p, "base_url", "") or "") != prefix:
        return None
    try:
        from delfin.agent import proc_identity
        if not proc_identity.alive(p.pid, p.proc_start, p.host):
            return None
    except Exception:
        return None
    return p


def forward_query(query: str) -> str:
    return urlencode([(k, v) for k, v in parse_qsl(query or "", keep_blank_values=True)
                      if k != "token"])


def forward_cookie(cookie_header: str, port: int) -> str:
    keep = []
    for part in (cookie_header or "").split(";"):
        name, _, _value = part.strip().partition("=")
        if name == "_xsrf" or (name.startswith("username-") and name.endswith(f"-{int(port)}")):
            keep.append(part.strip())
    return "; ".join(keep)


def forward_headers(headers, p) -> dict:
    out = {}
    for name, value in headers.get_all():
        if name.lower() in _HOP:
            continue
        out[name] = value
    out["Authorization"] = f"token {p.token}"
    cookie = forward_cookie(headers.get("Cookie", ""), p.port)
    if cookie:
        out["Cookie"] = cookie
    return out


def _xsrf_values(cookie_header: str) -> list:
    vals = []
    for part in (cookie_header or "").split(";"):
        name, _, value = part.strip().partition("=")
        if name == "_xsrf" and value:
            vals.append(value)
    return vals


def handlers():
    """``(http_handler, websocket_handler)``, built late so importing this
    module costs nothing outside a server."""
    from jupyter_server.base.handlers import JupyterHandler
    from jupyter_server.base.websocket import WebSocketMixin
    from tornado import httpclient, httputil, web, websocket

    def _missing():
        raise web.HTTPError(404)

    class PreviewProxyHandler(JupyterHandler):
        SUPPORTED_METHODS = ("GET", "HEAD", "POST", "PUT", "PATCH", "DELETE", "OPTIONS")

        def check_xsrf_cookie(self):
            if self.request.method in ("GET", "HEAD", "OPTIONS"):
                return
            if getattr(self, "token_authenticated", False):
                return
            token = (self.request.headers.get("X-Xsrftoken")
                     or self.request.headers.get("X-Csrftoken")
                     or self.get_argument("_xsrf", None))
            if not token:
                raise web.HTTPError(403, "'_xsrf' argument missing")
            _v, want, _t = self._decode_xsrf_token(token)
            for raw in _xsrf_values(self.request.headers.get("Cookie", "")):
                _v2, have, _t2 = self._decode_xsrf_token(raw)
                if want and have and hmac.compare_digest(want, have):
                    return
            raise web.HTTPError(403, "XSRF cookie does not match")

        @web.authenticated
        async def get(self, slug, _rest):
            await self._proxy(slug)

        @web.authenticated
        async def head(self, slug, _rest):
            await self._proxy(slug)

        @web.authenticated
        async def post(self, slug, _rest):
            await self._proxy(slug)

        @web.authenticated
        async def put(self, slug, _rest):
            await self._proxy(slug)

        @web.authenticated
        async def patch(self, slug, _rest):
            await self._proxy(slug)

        @web.authenticated
        async def delete(self, slug, _rest):
            await self._proxy(slug)

        @web.authenticated
        async def options(self, slug, _rest):
            await self._proxy(slug)

        async def _proxy(self, slug):
            base = self.settings.get("base_url", "/")
            p = target(slug, self.request.path, base)
            if p is None:
                _missing()
            query = forward_query(self.request.query)
            url = f"http://127.0.0.1:{int(p.port)}{self.request.path}" + (
                f"?{query}" if query else "")
            body = self.request.body if self.request.method in (
                "POST", "PUT", "PATCH", "DELETE") else None
            state = {"started": False, "code": 502, "reason": "", "headers": None}

            def _header(line):
                line = line.rstrip("\r\n")
                if line.startswith("HTTP/"):
                    parts = line.split(" ", 2)
                    state["code"] = int(parts[1])
                    state["reason"] = parts[2] if len(parts) > 2 else ""
                    state["headers"] = httputil.HTTPHeaders()
                elif line and state["headers"] is not None:
                    state["headers"].parse_line(line)

            def _start(code, reason, headers):
                state["started"] = True
                self.set_status(code, reason or None)
                if headers is None:
                    return
                for name, value in headers.get_all():
                    low = name.lower()
                    if low in _RESP_DROP:
                        continue
                    if low == "location":
                        value = re.sub(r"^https?://127\.0\.0\.1:\d+", "", value)
                    if low == "set-cookie":
                        self.add_header(name, value)
                    else:
                        self.set_header(name, value)

            def _chunk(data):
                if not state["started"]:
                    _start(state["code"], state["reason"], state["headers"])
                if data:
                    self.write(data)
                    self.flush()

            req = httpclient.HTTPRequest(
                url, method=self.request.method, headers=forward_headers(
                    self.request.headers, p), body=body, follow_redirects=False,
                decompress_response=False, allow_nonstandard_methods=True,
                connect_timeout=10, request_timeout=3600,
                header_callback=_header, streaming_callback=_chunk)
            client = httpclient.AsyncHTTPClient()
            try:
                resp = await client.fetch(req, raise_error=False)
            except Exception as exc:
                if state["started"]:
                    self.finish()
                    return
                raise web.HTTPError(502, f"preview not reachable: {exc}")
            if not state["started"]:
                if resp.code == 599:
                    raise web.HTTPError(502, "preview not reachable")
                _start(resp.code, resp.reason, resp.headers)
            self.finish()

    class PreviewProxyWebSocket(WebSocketMixin, websocket.WebSocketHandler, JupyterHandler):
        _up = None
        _chosen = None
        _pending: list = []

        async def get(self, slug, rest):
            if not self.current_user:
                raise web.HTTPError(403)
            base = self.settings.get("base_url", "/")
            p = target(slug, self.request.path, base)
            if p is None:
                _missing()
            query = forward_query(self.request.query)
            url = f"ws://127.0.0.1:{int(p.port)}{self.request.path}" + (
                f"?{query}" if query else "")
            protos = [s.strip() for s in self.request.headers.get(
                "Sec-WebSocket-Protocol", "").split(",") if s.strip()]
            headers = {"Authorization": f"token {p.token}"}
            cookie = forward_cookie(self.request.headers.get("Cookie", ""), p.port)
            if cookie:
                headers["Cookie"] = cookie
            self._pending = []
            try:
                self._up = await websocket.websocket_connect(
                    httpclient.HTTPRequest(url, headers=headers, connect_timeout=10),
                    subprotocols=protos or None, on_message_callback=self._from_up)
            except Exception as exc:
                raise web.HTTPError(502, f"preview not reachable: {exc}")
            self._chosen = self._up.selected_subprotocol
            return await super().get(slug, rest)

        def select_subprotocol(self, subprotocols):
            return self._chosen

        def open(self, *args, **kwargs):
            super().open(*args, **kwargs)
            pending, self._pending = self._pending, None
            for m in pending or []:
                self.write_message(m, binary=isinstance(m, bytes))

        def _from_up(self, message):
            if message is None:
                self.close()
                return
            if self._pending is not None:
                self._pending.append(message)
                return
            try:
                self.write_message(message, binary=isinstance(message, bytes))
            except websocket.WebSocketClosedError:
                pass

        def on_message(self, message):
            if self._up is not None:
                self._up.write_message(message, binary=isinstance(message, bytes))

        def on_close(self):
            if self._up is not None:
                self._up.close()

    return PreviewProxyHandler, PreviewProxyWebSocket


def register(serverapp) -> Optional[str]:
    """Add the routes to a running server and announce them to its kernels
    (``branch_preview.PROXY_ENV`` = the base URL). The route prefix, or None."""
    from delfin.dashboard import branch_preview as BP
    try:
        web_app = serverapp.web_app
        base = web_app.settings.get("base_url", "/") or "/"
        http_h, ws_h = handlers()
        root = "/" + base.strip("/")
        root = root.rstrip("/") + f"/{PREFIX}/({SLUG_RE})"
        web_app.add_handlers(".*$", [
            (root + r"/(.*/channels|api/events/subscribe)", ws_h),
            (root + r"/(.*)", http_h),
        ])
    except Exception:
        return None
    os.environ[BP.PROXY_ENV] = base
    return root


def _jupyter_server_extension_points():
    return [{"module": "delfin.dashboard.preview_proxy"}]


def _load_jupyter_server_extension(serverapp) -> None:
    if register(serverapp) is None:
        serverapp.log.warning("[delfin] the preview proxy could not be installed; "
                              "previews are linked by their own port")


load_jupyter_server_extension = _load_jupyter_server_extension
