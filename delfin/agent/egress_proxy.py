"""The one way out to the network for a sandboxed command.

A command that runs isolated (a locked session, an unattended run, a forced
isolation) reaches the network only through this proxy, which runs in the
agent process, outside every sandbox. It lets through the domains the user
allows and nothing else: a command cannot post a workspace file to an
arbitrary host, reach an internal service, or ask the cloud metadata
endpoint for credentials. The sandbox makes the proxy the only door --
a network namespace under bubblewrap, the socket guard under Landlock --
so a tool that ignores the proxy settings simply has no network.

The proxy listens on a Unix socket in a private (0700) directory, and the
sandbox helper forwards a loopback port to it. Every request must carry the
proxy token (``Proxy-Authorization``); the token reaches the command through
its environment only. Supported: ``CONNECT host:port`` (HTTPS, git over
HTTPS) and plain ``http://`` requests. A refused request answers 403 with
the reason and the setting that would allow it.

Settings (``agent.sandbox_network``): ``allowed_domains`` adds domains
(``example.org`` or ``*.example.org``), ``allowed_ports`` replaces the port
list, ``mode`` is ``"proxy"`` (default), ``"none"`` or ``"open"``.
"""
from __future__ import annotations

import base64
import os
import secrets
import select
import socket
import tempfile
import threading
from typing import Callable, Optional

DEFAULT_ALLOWED_DOMAINS = (
    # Python packages
    "pypi.org", "files.pythonhosted.org",
    # conda / micromamba
    "conda.anaconda.org", "repo.anaconda.com", "repo.prefix.dev",
    # git over HTTPS and release downloads
    "github.com", "api.github.com", "codeload.github.com",
    "objects.githubusercontent.com", "raw.githubusercontent.com",
    "release-assets.githubusercontent.com",
)
DEFAULT_ALLOWED_PORTS = (80, 443)

_MAX_HEADER = 64 * 1024


def _settings() -> dict:
    try:
        from delfin.user_settings import load_settings
        agent = (load_settings() or {}).get("agent") or {}
        net = agent.get("sandbox_network") or {}
        return net if isinstance(net, dict) else {}
    except Exception:
        return {}


def network_mode() -> str:
    """``proxy`` (default), ``none`` or ``open``, from the settings."""
    mode = str(_settings().get("mode") or "proxy").strip().lower()
    return mode if mode in ("proxy", "none", "open") else "proxy"


def domain_allowed(host: str, allowed) -> bool:
    host = (host or "").strip().lower().rstrip(".")
    if not host:
        return False
    for entry in allowed:
        entry = str(entry or "").strip().lower().rstrip(".")
        if not entry:
            continue
        if entry.startswith("*."):
            if host.endswith(entry[1:]) and host != entry[2:]:
                return True
        elif host == entry:
            return True
    return False


class EgressProxy:
    """A minimal forward proxy with a domain allow-list. Thread per client."""

    def __init__(self, allowed_domains=None, allowed_ports=None, *,
                 connector: Optional[Callable[[str, int, float], socket.socket]] = None,
                 on_refusal: Optional[Callable[[str], None]] = None):
        conf = _settings()
        extra = conf.get("allowed_domains") or []
        self.allowed_domains = tuple(allowed_domains if allowed_domains is not None
                                     else list(DEFAULT_ALLOWED_DOMAINS) + list(extra))
        ports = allowed_ports if allowed_ports is not None else conf.get("allowed_ports")
        self.allowed_ports = tuple(int(p) for p in (ports or DEFAULT_ALLOWED_PORTS))
        self.token = secrets.token_urlsafe(32)
        self._connector = connector or self._connect_public
        self._on_refusal = on_refusal
        self._dir = tempfile.mkdtemp(prefix="delfin-egress-")
        os.chmod(self._dir, 0o700)
        self.socket_path = os.path.join(self._dir, "proxy.sock")
        # The sandbox helper reads the token from here (0600), never from
        # its command line, which every user on the machine can read.
        token_path = os.path.join(self._dir, "token")
        fd = os.open(token_path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        with os.fdopen(fd, "w") as fh:
            fh.write(self.token)
        self._server = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        self._server.bind(self.socket_path)
        os.chmod(self.socket_path, 0o600)
        self._server.listen(64)
        self._closed = False
        threading.Thread(target=self._accept_loop, name="delfin-egress-proxy",
                         daemon=True).start()

    # -- lifecycle ----------------------------------------------------------
    def close(self) -> None:
        self._closed = True
        try:
            self._server.close()
        except OSError:
            pass
        for p in (self.socket_path, os.path.join(self._dir, "token")):
            try:
                os.unlink(p)
            except OSError:
                pass
        try:
            os.rmdir(self._dir)
        except OSError:
            pass

    def _accept_loop(self) -> None:
        while not self._closed:
            try:
                conn, _ = self._server.accept()
            except OSError:
                return
            threading.Thread(target=self._serve, args=(conn,), daemon=True).start()

    # -- policy -------------------------------------------------------------
    @staticmethod
    def _connect_public(host: str, port: int, timeout: float) -> socket.socket:
        from .web_tools import _checked_create_connection
        return _checked_create_connection((host, port), timeout=timeout)

    def _refuse(self, conn: socket.socket, code: int, reason: str) -> None:
        body = (reason + "\n").encode("utf-8", "replace")
        head = (f"HTTP/1.1 {code} {'Forbidden' if code == 403 else 'Proxy Authentication Required' if code == 407 else 'Bad Gateway'}\r\n"
                "Content-Type: text/plain; charset=utf-8\r\n"
                f"Content-Length: {len(body)}\r\n"
                + ('Proxy-Authenticate: Basic realm="delfin"\r\n' if code == 407 else "")
                + "Connection: close\r\n\r\n").encode()
        try:
            conn.sendall(head + body)
        except OSError:
            pass
        if code == 403 and self._on_refusal is not None:
            try:
                self._on_refusal(reason)
            except Exception:
                pass

    def _authorised(self, headers: dict) -> bool:
        value = headers.get("proxy-authorization", "")
        if not value.lower().startswith("basic "):
            return False
        try:
            user_pass = base64.b64decode(value[6:].strip()).decode("utf-8", "replace")
        except Exception:
            return False
        _user, _, password = user_pass.partition(":")
        return secrets.compare_digest(password, self.token)

    def _check(self, host: str, port: int) -> Optional[str]:
        if port not in self.allowed_ports:
            return (f"blocked by DELFIN's sandbox: port {port} is not allowed "
                    f"(allowed: {', '.join(map(str, self.allowed_ports))}; "
                    "setting agent.sandbox_network.allowed_ports)")
        if not domain_allowed(host, self.allowed_domains):
            return (f"blocked by DELFIN's sandbox: {host} is not an allowed "
                    "domain. Ask the user to add it to "
                    "agent.sandbox_network.allowed_domains.")
        return None

    # -- serving ------------------------------------------------------------
    def _serve(self, conn: socket.socket) -> None:
        upstream = None
        try:
            conn.settimeout(30)
            data = b""
            while b"\r\n\r\n" not in data:
                chunk = conn.recv(4096)
                if not chunk:
                    return
                data += chunk
                if len(data) > _MAX_HEADER:
                    self._refuse(conn, 400, "request header too large")
                    return
            head, _, rest = data.partition(b"\r\n\r\n")
            lines = head.decode("latin-1").split("\r\n")
            parts = lines[0].split(" ")
            if len(parts) != 3:
                self._refuse(conn, 400, "malformed request line")
                return
            method, target, version = parts
            headers = {}
            for line in lines[1:]:
                name, _, value = line.partition(":")
                headers[name.strip().lower()] = value.strip()
            if not self._authorised(headers):
                self._refuse(conn, 407, "proxy token missing or wrong")
                return
            if method.upper() == "CONNECT":
                host, _, port_s = target.rpartition(":")
                host = host.strip("[]")
                try:
                    port = int(port_s)
                except ValueError:
                    self._refuse(conn, 400, "CONNECT target needs host:port")
                    return
                why = self._check(host, port)
                if why:
                    self._refuse(conn, 403, why)
                    return
                upstream = self._open(conn, host, port)
                if upstream is None:
                    return
                conn.sendall(b"HTTP/1.1 200 Connection established\r\n\r\n")
                if rest:
                    upstream.sendall(rest)
            else:
                from urllib.parse import urlsplit
                url = urlsplit(target)
                if url.scheme != "http" or not url.hostname:
                    self._refuse(conn, 400, "only absolute http:// URLs or CONNECT")
                    return
                port = url.port or 80
                why = self._check(url.hostname, port)
                if why:
                    self._refuse(conn, 403, why)
                    return
                upstream = self._open(conn, url.hostname, port)
                if upstream is None:
                    return
                path = url.path or "/"
                if url.query:
                    path += "?" + url.query
                kept = [line for line in lines[1:]
                        if line.split(":", 1)[0].strip().lower()
                        not in ("proxy-authorization", "proxy-connection", "connection")]
                request = "\r\n".join([f"{method} {path} {version}", *kept,
                                       "Connection: close"]) + "\r\n\r\n"
                upstream.sendall(request.encode("latin-1") + rest)
            self._relay(conn, upstream)
        except OSError:
            pass
        finally:
            for s in (conn, upstream):
                if s is not None:
                    try:
                        s.close()
                    except OSError:
                        pass

    def _open(self, conn, host: str, port: int):
        try:
            return self._connector(host, port, 30.0)
        except OSError as exc:
            reason = str(exc)
            if reason.startswith("refused:"):
                self._refuse(conn, 403, f"blocked by DELFIN's sandbox: {reason[8:].strip()}")
            else:
                self._refuse(conn, 502, f"could not reach {host}:{port}: {exc}")
            return None

    @staticmethod
    def _relay(a: socket.socket, b: socket.socket) -> None:
        a.settimeout(None)
        b.settimeout(None)
        socks = [a, b]
        while True:
            readable, _, broken = select.select(socks, [], socks, 300)
            if broken or not readable:
                return
            for s in readable:
                other = b if s is a else a
                try:
                    chunk = s.recv(65536)
                except OSError:
                    return
                if not chunk:
                    return
                try:
                    other.sendall(chunk)
                except OSError:
                    return


_PROXY: Optional[EgressProxy] = None
_LOCK = threading.Lock()


def get_proxy() -> EgressProxy:
    """The agent process's proxy, started on first use."""
    global _PROXY
    with _LOCK:
        if _PROXY is None or _PROXY._closed:
            def _record(reason: str) -> None:
                try:
                    from .api_client import _record_security_event
                    _record_security_event("egress_blocked", "bash", reason[:200],
                                           blocked=True)
                except Exception:
                    pass
            _PROXY = EgressProxy(on_refusal=_record)
            import atexit
            atexit.register(_PROXY.close)
        return _PROXY


def slurm_controllers() -> list[str]:
    """``host:port`` of the SLURM controller(s) from slurm.conf, so sbatch,
    squeue and scancel keep working inside the sandbox. Empty without SLURM."""
    path = os.environ.get("SLURM_CONF") or "/etc/slurm/slurm.conf"
    hosts: list[str] = []
    port = "6817"
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for raw in fh:
                line = raw.split("#", 1)[0].strip()
                if not line or "=" not in line:
                    continue
                key, _, value = line.partition("=")
                key = key.strip().lower()
                value = value.strip()
                if key in ("slurmctldhost", "controlmachine", "backupcontroller"):
                    name = value.split("(", 1)[0].strip()
                    if name:
                        hosts.append(name)
                elif key == "slurmctldport":
                    port = value.split("-", 1)[0].strip() or port
    except OSError:
        return []
    out = []
    for name in hosts:
        try:
            for info in socket.getaddrinfo(name, int(port), 0, socket.SOCK_STREAM):
                out.append(f"{info[4][0]}:{port}")
        except (OSError, ValueError):
            continue
    return sorted(set(out))
