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
list, ``mode`` is ``"proxy"`` (default), ``"none"`` or ``"open"``,
``upload_ask_bytes`` the request body size above which a request counts
as an upload (default 1 MiB).

Decisions per request (``EgressProxy.decide``), for the session whose
token the request carries (``token_for``):

1. a port outside ``allowed_ports`` is refused;
2. a ``web_fetch(domain:...)`` permission rule (``permission_rules``,
   when that module is present) denies, asks or allows first;
3. a KNOWN host -- the defaults, ``allowed_domains``, and what the
   session itself is configured to use (its model endpoint, its git
   remotes, its pip and conda indexes) -- passes; an upload to it
   (a body above ``upload_ask_bytes``) is asked when someone can be
   asked and passes otherwise, so an unattended push or model call works;
4. an UNKNOWN host, and any POST/PUT/PATCH to one, is asked when the
   session has an approval channel and is not in plan or bypass mode;
   otherwise refused with a message naming the host and the rule that
   would allow it. An answer is kept for the session.

The process-wide token (``self.token``) carries no session: allow-list
only, unknown refused -- the behaviour before sessions existed.

The mode is read when the process first asks; a later change to the
settings file may tighten it at once but loosens it only at the next
start, so a command that edits the file cannot switch the proxy off for
the session it runs in.
"""
from __future__ import annotations

import base64
import ipaddress
import os
import re
import secrets
import select
import socket
import tempfile
import threading
import weakref
from typing import Callable, Optional
from urllib.parse import urlsplit

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
DEFAULT_UPLOAD_ASK_BYTES = 1024 * 1024

_MAX_HEADER = 64 * 1024
_UPLOAD_METHODS = ("POST", "PUT", "PATCH")


def _settings() -> dict:
    try:
        from delfin.user_settings import load_settings
        agent = (load_settings() or {}).get("agent") or {}
        net = agent.get("sandbox_network") or {}
        return net if isinstance(net, dict) else {}
    except Exception:
        return {}


_MODE_ORDER = {"none": 0, "proxy": 1, "open": 2}
_MODE_AT_START: Optional[str] = None


def network_mode() -> str:
    """``proxy`` (default), ``none`` or ``open``, from the settings.

    The first answer in a process is kept as its ceiling: a later value
    that is stricter applies, a looser one waits for the next start."""
    global _MODE_AT_START
    mode = str(_settings().get("mode") or "proxy").strip().lower()
    mode = mode if mode in _MODE_ORDER else "proxy"
    # A managed domain ceiling holds only behind the proxy, so "open" is
    # not available under one (managed_policy.allowed_domains_max).
    if mode == "open":
        try:
            from . import managed_policy as _mp
            if _mp.network_ceiling_active():
                mode = "proxy"
        except Exception:
            mode = "proxy"
    if _MODE_AT_START is None:
        _MODE_AT_START = mode
        return mode
    if _MODE_ORDER[mode] > _MODE_ORDER[_MODE_AT_START]:
        return _MODE_AT_START
    return mode


# -- what the session itself is configured to reach ------------------------

_ENDPOINT_HOSTS: set = set()


def _host_of(url: str) -> str:
    try:
        return (urlsplit(str(url or "").strip()).hostname or "").lower()
    except ValueError:
        return ""


def note_endpoint(url: str) -> None:
    """Register the host of a model endpoint this process talks to, so a
    script the agent starts can call the same endpoint."""
    host = _host_of(url)
    if host:
        _ENDPOINT_HOSTS.add(host)


def _index_urls() -> list[str]:
    """pip and conda index URLs from the environment and their config files."""
    urls: list[str] = []
    for var in ("PIP_INDEX_URL", "PIP_EXTRA_INDEX_URL", "UV_INDEX_URL",
                "UV_EXTRA_INDEX_URL", "CONDA_CHANNEL_ALIAS"):
        urls += (os.environ.get(var) or "").split()
    home = os.path.expanduser("~")
    files = [os.path.join(home, ".config", "pip", "pip.conf"),
             os.path.join(home, ".pip", "pip.conf"), "/etc/pip.conf",
             os.path.join(home, ".condarc"), os.path.join(home, ".mambarc"),
             os.environ.get("CONDARC") or "", os.environ.get("PIP_CONFIG_FILE") or ""]
    for path in files:
        if not path:
            continue
        try:
            with open(path, encoding="utf-8", errors="replace") as fh:
                urls += re.findall(r"https?://[^\s'\",]+", fh.read(65536))
        except OSError:
            continue
    return urls


_GIT_REMOTES: dict = {}


def git_remotes(workspace) -> list[str]:
    """The workspace repository's remote URLs (fetch and push). Cached."""
    key = str(workspace or "")
    if not key:
        return []
    if key not in _GIT_REMOTES:
        out: list[str] = []
        try:
            from . import safe_git
            done = safe_git.run(["remote", "-v"], key, timeout=5)
            if done.returncode == 0:
                for line in done.stdout.splitlines():
                    parts = line.split()
                    if len(parts) >= 2 and parts[1] not in out:
                        out.append(parts[1])
        except Exception:
            out = []
        _GIT_REMOTES[key] = out
    return _GIT_REMOTES[key]


def configured_domains(workspace=None) -> list[str]:
    """Hosts the session is configured to use beyond the defaults: the
    ``allowed_domains`` setting, the model endpoints, the pip/conda indexes
    and the workspace's HTTP(S) git remotes."""
    out = [str(d) for d in (_settings().get("allowed_domains") or [])]
    out += sorted(_ENDPOINT_HOSTS)
    out += [h for h in (_host_of(u) for u in _index_urls()) if h]
    out += [h for h in (_host_of(u) for u in git_remotes(workspace)
                        if u.startswith(("http://", "https://"))) if h]
    return out


def _resolve(host: str, port: int) -> list[str]:
    try:
        return sorted({f"{info[4][0]}:{port}" for info in socket.getaddrinfo(
            host, port, 0, socket.SOCK_STREAM)})
    except (OSError, ValueError, UnicodeError):
        return []


def direct_tcp_endpoints(workspace=None) -> list[str]:
    """``ip:port`` a command may reach WITHOUT the proxy, because the
    protocol is not HTTP: the SLURM controller, git remotes over SSH, and
    license servers named in ``*_LICENSE_FILE`` (``port@host``; port 0 =
    every port of that host, since a FlexLM vendor daemon listens on a
    port of its own). Cached per workspace: it resolves names."""
    key = str(workspace or "")
    if key not in _DIRECT_TCP:
        _DIRECT_TCP[key] = _direct_tcp_endpoints(workspace)
    return list(_DIRECT_TCP[key])


_DIRECT_TCP: dict = {}


def _direct_tcp_endpoints(workspace) -> list[str]:
    out = list(slurm_controllers())
    for url in git_remotes(workspace):
        m = re.match(r"^ssh://(?:[^@/]+@)?(\[[^\]]+\]|[^:/]+)(?::(\d+))?/", url)
        if m:
            out += _resolve(m.group(1).strip("[]"), int(m.group(2) or 22))
            continue
        m = re.match(r"^(?:[^@/:]+@)?([A-Za-z0-9.\-]+):(?!//)", url)
        if m and "://" not in url:
            out += _resolve(m.group(1), 22)
    for var, value in os.environ.items():
        if not var.upper().endswith("LICENSE_FILE"):
            continue
        for entry in re.split(r"[:;]", value or ""):
            m = re.match(r"^\s*(\d*)@([A-Za-z0-9.\-]+)\s*$", entry)
            if m:
                out += [f"{hp.rsplit(':', 1)[0]}:0"
                        for hp in _resolve(m.group(2), int(m.group(1) or 27000))]
    return sorted(set(out))


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


def _short_base() -> Optional[str]:
    """A directory whose paths stay inside a Unix socket's 104-108 bytes.

    TMPDIR on a cluster is often a long job scratch path, where binding
    the proxy socket failed and an isolated command had no network."""
    candidates = [tempfile.gettempdir(), "/tmp", "/var/tmp", "/dev/shm"]
    for base in candidates:
        try:
            if (len(base) <= 60 and os.path.isdir(base)
                    and os.access(base, os.W_OK | os.X_OK)):
                return base
        except OSError:
            continue
    return None


class EgressProxy:
    """A minimal forward proxy with a domain allow-list. Thread per client."""

    def __init__(self, allowed_domains=None, allowed_ports=None, *,
                 connector: Optional[Callable[[str, int, float], socket.socket]] = None,
                 on_refusal: Optional[Callable[[str], None]] = None,
                 upload_ask_bytes: Optional[int] = None):
        conf = _settings()
        extra = conf.get("allowed_domains") or []
        self.allowed_domains = tuple(allowed_domains if allowed_domains is not None
                                     else list(DEFAULT_ALLOWED_DOMAINS) + list(extra))
        ports = allowed_ports if allowed_ports is not None else conf.get("allowed_ports")
        self.allowed_ports = tuple(int(p) for p in (ports or DEFAULT_ALLOWED_PORTS))
        try:
            limit = int(upload_ask_bytes if upload_ask_bytes is not None
                        else conf.get("upload_ask_bytes") or DEFAULT_UPLOAD_ASK_BYTES)
        except (TypeError, ValueError):
            limit = DEFAULT_UPLOAD_ASK_BYTES
        self.upload_ask_bytes = max(0, limit)
        self.token = secrets.token_urlsafe(32)
        self._connector = connector or self._connect_public
        self._on_refusal = on_refusal
        # token -> weakref to the session's permissions; id(perms) -> token.
        self._sessions: dict = {}
        self._session_tokens: dict = {}
        self._token_files: dict = {}
        # (id(perms), host, kind) -> the user's answer, kept for the session.
        self._answers: dict = {}
        self._ask_lock = threading.Lock()
        self._dir = tempfile.mkdtemp(prefix="dlf-egress-", dir=_short_base())
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
        self.tcp_port = 0
        self._tcp = None
        threading.Thread(target=self._accept_loop, args=(self._server,),
                         name="delfin-egress-proxy", daemon=True).start()

    def ensure_tcp(self) -> int:
        """A loopback port as well, for a sandbox without a forwarder of its
        own (macOS Seatbelt). The token is required there as everywhere."""
        if self._tcp is None:
            self._tcp = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            self._tcp.bind(("127.0.0.1", 0))
            self._tcp.listen(64)
            self.tcp_port = self._tcp.getsockname()[1]
            threading.Thread(target=self._accept_loop, args=(self._tcp,),
                             name="delfin-egress-proxy-tcp", daemon=True).start()
        return self.tcp_port

    # -- lifecycle ----------------------------------------------------------
    def close(self) -> None:
        self._closed = True
        for srv in (self._server, self._tcp):
            if srv is None:
                continue
            try:
                srv.close()
            except OSError:
                pass
        for p in (self.socket_path, os.path.join(self._dir, "token"),
                  *self._token_files.values()):
            try:
                os.unlink(p)
            except OSError:
                pass
        try:
            os.rmdir(self._dir)
        except OSError:
            pass

    def _accept_loop(self, server) -> None:
        while not self._closed:
            try:
                conn, _ = server.accept()
            except OSError:
                return
            threading.Thread(target=self._serve, args=(conn,), daemon=True).start()

    # -- policy -------------------------------------------------------------
    @staticmethod
    def _connect_public(host: str, port: int, timeout: float) -> socket.socket:
        """Public addresses only -- except for a host the user configured
        (a model endpoint, an index, a remote, ``allowed_domains``), which
        may be on the private network. Link-local and the cloud metadata
        addresses never."""
        from .web_tools import _checked_create_connection
        if not domain_allowed(host, configured_domains()):
            return _checked_create_connection((host, port), timeout=timeout)
        last: Optional[BaseException] = None
        for af, socktype, proto, _c, sockaddr in socket.getaddrinfo(
                host, port, 0, socket.SOCK_STREAM):
            ip = ipaddress.ip_address(str(sockaddr[0]).split("%", 1)[0])
            if getattr(ip, "ipv4_mapped", None) is not None:
                ip = ip.ipv4_mapped
            if ip.is_link_local or str(ip) in ("100.100.100.200", "168.63.129.16"):
                last = OSError(f"refused: link-local / metadata address {ip}")
                continue
            sock = socket.socket(af, socktype, proto)
            sock.settimeout(timeout)
            try:
                sock.connect(sockaddr)
                return sock
            except OSError as exc:
                last = exc
                sock.close()
        raise last or OSError(f"cannot resolve {host}")

    # -- sessions -----------------------------------------------------------
    def token_for(self, perms) -> str:
        """The token that identifies ``perms``'s commands to the proxy; the
        process token for None."""
        if perms is None:
            return self.token
        key = id(perms)
        token = self._session_tokens.get(key)
        if token is not None and self._sessions.get(token, lambda: None)() is perms:
            return token
        token = secrets.token_urlsafe(32)
        self._session_tokens[key] = token
        try:
            self._sessions[token] = weakref.ref(perms)
            weakref.finalize(perms, self._forget, key, token)
        except TypeError:
            self._sessions[token] = lambda _p=perms: _p
        return token

    def token_file_for(self, perms) -> str:
        """A 0600 file holding ``token_for(perms)``: the sandbox helper reads
        it there, never from its command line."""
        if perms is None:
            return os.path.join(self._dir, "token")
        token = self.token_for(perms)
        path = self._token_files.get(token)
        if path is None:
            path = os.path.join(self._dir, "token-" + secrets.token_hex(8))
            fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
            with os.fdopen(fd, "w") as fh:
                fh.write(token)
            self._token_files[token] = path
        return path

    def _forget(self, key: int, token: str) -> None:
        self._sessions.pop(token, None)
        if self._session_tokens.get(key) == token:
            self._session_tokens.pop(key, None)
        for k in [k for k in self._answers if k[0] == key]:
            self._answers.pop(k, None)
        path = self._token_files.pop(token, None)
        if path:
            try:
                os.unlink(path)
            except OSError:
                pass

    # -- decisions ----------------------------------------------------------
    @staticmethod
    def _rule_verdict(host: str, perms) -> Optional[tuple]:
        """``(verdict, rule text)`` from a ``web_fetch(domain:...)``
        permission rule, or None. Uses ``permission_rules`` when present;
        without it there are no rules to apply."""
        try:
            from . import permission_rules as pr
        except ImportError:
            return None
        try:
            ws = getattr(perms, "workspace", None)
            rs = pr.load(ws)
            args = {"url": f"https://{host}/"}
            hit = pr.friction(rs, "web_fetch", args, ws)
            if hit:
                return hit[0], str(getattr(hit[1], "text", hit[1]))
            rule = pr.allows(rs, "tool", "web_fetch", args, ws)
            if rule is not None:
                return "allow", str(getattr(rule, "text", rule))
        except Exception:
            return None
        return None

    def known(self, host: str, perms=None) -> bool:
        return domain_allowed(host, list(self.allowed_domains) + configured_domains(
            getattr(perms, "workspace", None)))

    def decide(self, host: str, port: int, *, method: str = "CONNECT",
               body_bytes: int = 0, perms=None) -> tuple[str, str, bool]:
        """``(verdict, reason, upload)``; verdict is allow, ask or deny."""
        if port not in self.allowed_ports:
            return "deny", (
                f"blocked by DELFIN's egress proxy: port {port} is not allowed "
                f"(allowed: {', '.join(map(str, self.allowed_ports))}; "
                "setting agent.sandbox_network.allowed_ports)"), False
        # The administrator's ceiling comes before every user rule, known
        # host and question: outside it nothing reaches the host. An
        # unreadable trusted policy is handled (closed) by managed_policy;
        # an error here keeps main's previous answer.
        try:
            from . import managed_policy as _mp
            within = _mp.domain_within_max(host)
        except Exception:
            within = True
        if not within:
            return "deny", (
                f"blocked by DELFIN's egress proxy: {host} is outside the "
                "domains the administrator's managed policy allows "
                "(network.allowed_domains_max); no user setting adds it"), False
        method = (method or "").upper()
        rule = self._rule_verdict(host, perms)
        if rule and rule[0] == "deny":
            return "deny", (f"blocked by DELFIN's egress proxy: the permission "
                            f"rule {rule[1]} denies {host}"), False
        known = self.known(host, perms) or bool(rule and rule[0] == "allow")
        big = body_bytes > self.upload_ask_bytes
        upload = (method in _UPLOAD_METHODS and not known) or big
        if rule and rule[0] == "ask":
            return "ask", f"the permission rule {rule[1]} asks for {host}", upload
        if known and not big:
            return "allow", "", False
        if known:
            return "ask", f"an upload of {body_bytes} bytes to {host}", True
        return "ask", f"{host} is not an allowed domain", upload

    def _no_one_to_ask(self, host: str, why: str) -> str:
        return (f"blocked by DELFIN's egress proxy: {why}, and nobody can be "
                f"asked. To allow {host}, the user adds the permission rule "
                f"web_fetch(domain:{host}) or adds {host} to "
                "agent.sandbox_network.allowed_domains in ~/.delfin_settings.json.")

    def settle(self, host: str, port: int, *, method: str = "CONNECT",
               body_bytes: int = 0, perms=None) -> Optional[str]:
        """Apply ``decide``: None to proceed, else the refusal reason. Asks
        the session's approval channel where there is one and the mode
        allows a question; a known host's upload passes unasked otherwise."""
        verdict, why, upload = self.decide(host, port, method=method,
                                           body_bytes=body_bytes, perms=perms)
        if verdict == "allow":
            return None
        if verdict == "deny":
            return why
        known = self.known(host, perms)
        mode = str(getattr(perms, "mode", "") or "")
        callback = getattr(perms, "confirm_callback", None)
        rule_ask = why.startswith("the permission rule")
        can_ask = callback is not None and (
            rule_ask or mode not in ("plan", "bypassPermissions"))
        if not can_ask:
            if known and upload and not rule_ask:
                _note_event("egress_upload", f"{host}: {why} (passed unasked)", False)
                return None
            return self._no_one_to_ask(host, why)
        kind = "upload" if upload else "domain"
        key = (id(perms), host.lower(), kind)
        with self._ask_lock:
            if key not in self._answers:
                what = (f"UPLOAD {body_bytes or 'an unknown number of'} bytes "
                        f"({method or 'CONNECT'}) to {host}:{port}" if upload
                        else f"connect to {host}:{port}")
                preview = (f"A program started by the bash tool wants to {what} "
                           f"through DELFIN's egress proxy ({why}).")
                args = {"url": f"https://{host}:{port}/", "method": method,
                        "bytes": int(body_bytes), "source": "bash child process"}
                try:
                    ok = bool(callback("web_fetch", args, preview))
                except Exception:
                    ok = False
                self._answers[key] = ok
                _note_event("egress_asked", f"{host} {kind}: "
                            + ("approved" if ok else "refused"), not ok)
            ok = self._answers[key]
        if ok:
            return None
        return (f"blocked by DELFIN's egress proxy: the user refused "
                f"{'the upload to ' if upload else ''}{host}. Do not retry it "
                "or reach the same host another way.")

    def _refuse(self, conn: socket.socket, code: int, reason: str) -> None:
        body = (reason + "\n").encode("utf-8", "replace")
        head = (f"HTTP/1.1 {code} {'Forbidden' if code == 403 else 'Proxy Authentication Required' if code == 407 else 'Bad Gateway'}\r\n"
                "Content-Type: text/plain; charset=utf-8\r\n"
                f"Content-Length: {len(body)}\r\n"
                + ('Proxy-Authenticate: Basic realm="delfin"\r\n' if code == 407 else "")
                + "Connection: close\r\n\r\n").encode()
        # The refusal is recorded BEFORE the answer goes out. Recorded after
        # it, a client could read its 403 and look at the record before the
        # serving thread had written it (SLURM 7210626: proxy.refused == []
        # beside a correct 403) -- a security record trailing the event.
        if code == 403 and self._on_refusal is not None:
            try:
                self._on_refusal(reason)
            except Exception:
                pass
        try:
            conn.sendall(head + body)
        except OSError:
            pass

    def _authorised(self, headers: dict) -> tuple[bool, object]:
        """``(ok, perms)``: perms is None for the process token."""
        value = headers.get("proxy-authorization", "")
        if not value.lower().startswith("basic "):
            return False, None
        try:
            user_pass = base64.b64decode(value[6:].strip()).decode("utf-8", "replace")
        except Exception:
            return False, None
        _user, _, password = user_pass.partition(":")
        if secrets.compare_digest(password, self.token):
            return True, None
        ref = self._sessions.get(password)
        if ref is None:
            return False, None
        perms = ref()
        return perms is not None, perms

    def _check(self, host: str, port: int) -> Optional[str]:
        """The session-less decision (process token)."""
        return self.settle(host, port)

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
            ok, perms = self._authorised(headers)
            if not ok:
                self._refuse(conn, 407, "proxy token missing or wrong")
                return
            meter = None
            if method.upper() == "CONNECT":
                host, _, port_s = target.rpartition(":")
                host = host.strip("[]")
                try:
                    port = int(port_s)
                except ValueError:
                    self._refuse(conn, 400, "CONNECT target needs host:port")
                    return
                why = self.settle(host, port, perms=perms)
                if why:
                    self._refuse(conn, 403, why)
                    return
                if not self.known(host, perms) and self.upload_ask_bytes:
                    # Inside a tunnel the method is not visible; what leaves
                    # is. A host the user admitted by answering a question
                    # is asked again once its upload passes the limit.
                    meter = (lambda sent, _h=host, _p=port, _perms=perms:
                             self.settle(_h, _p, method="CONNECT",
                                         body_bytes=sent, perms=_perms) is None)
                upstream = self._open(conn, host, port)
                if upstream is None:
                    return
                conn.sendall(b"HTTP/1.1 200 Connection established\r\n\r\n")
                if rest:
                    upstream.sendall(rest)
            else:
                url = urlsplit(target)
                if url.scheme != "http" or not url.hostname:
                    self._refuse(conn, 400, "only absolute http:// URLs or CONNECT")
                    return
                port = url.port or 80
                if "chunked" in headers.get("transfer-encoding", "").lower():
                    size = self.upload_ask_bytes + 1      # unknown: counts as large
                else:
                    try:
                        size = int(headers.get("content-length") or 0)
                    except ValueError:
                        size = self.upload_ask_bytes + 1
                why = self.settle(url.hostname, port, method=method,
                                  body_bytes=size, perms=perms)
                if why:
                    self._refuse(conn, 403, why)
                    # Read and drop the body the client is still sending, so
                    # it reaches the 403 instead of a reset connection.
                    try:
                        conn.shutdown(socket.SHUT_WR)
                        conn.settimeout(5)
                        while conn.recv(65536):
                            pass
                    except OSError:
                        pass
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
            self._relay(conn, upstream, meter, self.upload_ask_bytes)
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
    def _relay(a: socket.socket, b: socket.socket, meter=None, limit: int = 0) -> None:
        """Copy both ways. ``meter(sent)`` is consulted once when the bytes
        from ``a`` (the command) pass ``limit``; False ends the tunnel."""
        a.settimeout(None)
        b.settimeout(None)
        socks = [a, b]
        sent = 0
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
                if s is a and meter is not None:
                    sent += len(chunk)
                    if sent > limit:
                        if not meter(sent):
                            return
                        meter = None
                try:
                    other.sendall(chunk)
                except OSError:
                    return


def _note_event(kind: str, detail: str, blocked: bool) -> None:
    try:
        from .api_client import _record_security_event
        _record_security_event(kind, "bash", detail[:200], blocked=blocked)
    except Exception:
        pass


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
