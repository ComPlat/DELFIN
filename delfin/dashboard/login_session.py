"""Server-side login sessions of the dashboard: logout and idle timeout.

jupyter-server's login cookie is a signed, self-contained value (the user's
name and display fields). Logging out asks the browser to drop it, and the
server keeps no record, so a copy of the cookie -- from a backup of the
browser profile, a shared machine, a proxy log -- stays valid for the whole
cookie lifetime. :class:`Sessions` is the record that was missing. Every
login cookie carries a random session id (``delfin_session``); the store
holds, per id, when the login was made and when it was last used:

* a request whose session id is unknown, closed, idle for longer than
  :func:`idle_seconds` or older than the login lifetime is not
  authenticated (jupyter-server then clears the cookie and shows the login
  page, which takes the token and continues to the same address);
* ``GET /logout`` closes the session of the cookie it carries;
* any authenticated request is use; the page itself reports user input on
  :data:`ALIVE_ROUTE` at most every five minutes (:func:`keepalive_js`),
  because a page that is being used talks to its kernel over a websocket
  that is authenticated once, when it opens.

Input: session ids (``secrets.token_urlsafe``), a clock. Output: ``touch``
answers whether a session may continue. The store is a JSON file,
``0600``, written by replacing it (never through a link), holding the
SHA-256 of each id, never the id; it survives a restart of the server, so
a browser that was logged in stays logged in and a logout stays a logout.
A store that cannot be read starts empty (everyone logs in again); one that
cannot be written keeps working from memory and says so in the log.

Rejected: closing open websockets when a session ends. jupyter-server has
no hook for it, and the kernel a websocket talks to may be a kept session
whose run must continue; an open page keeps its connection until it
reloads, and a page whose login ended shows a notice (:func:`keepalive_js`).
Rejected: an idle timeout counted on HTTP requests alone -- an open page
makes none while idle (measured: 0 requests in 60 s), but a page in use
makes none either.
"""

from __future__ import annotations

import hashlib
import json
import os
import secrets
import time
from pathlib import Path
from typing import Any, Callable, Optional

IDLE_ENV = "DELFIN_DASHBOARD_IDLE_HOURS"
DEFAULT_IDLE_HOURS = 12.0
MIN_IDLE_SECONDS = 15 * 60

#: The route the page reports user activity on (authenticated, 204).
ALIVE_ROUTE = "/delfin/session/alive"

#: File-name prefix of the store in jupyter-server's runtime directory.
STORE_PREFIX = "delfin-dashboard-sessions-"

#: Sessions kept at most; the least recently used go first.
MAX_SESSIONS = 512

#: Field of the login cookie that carries the session id.
COOKIE_FIELD = "delfin_session"

#: A use is written to the store at most this often (seconds): the store
#: is a file, and every request is a use.
_SAVE_EVERY = 60.0


def idle_seconds() -> float:
    """How long a login may go unused, in seconds.

    ``DELFIN_DASHBOARD_IDLE_HOURS`` (default 12), at least 15 minutes and at
    most the login lifetime (``DELFIN_DASHBOARD_LOGIN_DAYS``).
    """
    from delfin.dashboard.web_security import login_days

    try:
        hours = float(os.environ.get(IDLE_ENV, "") or DEFAULT_IDLE_HOURS)
    except ValueError:
        hours = DEFAULT_IDLE_HOURS
    if hours != hours:            # NaN
        hours = DEFAULT_IDLE_HOURS
    return min(login_days() * 86400.0, max(float(MIN_IDLE_SECONDS), hours * 3600.0))


def _digest(sid: str) -> str:
    return hashlib.sha256(sid.encode("utf-8")).hexdigest()


class Sessions:
    """The live logins of one dashboard server."""

    def __init__(self, path: Optional[Path], *, idle: float, lifetime: float,
                 clock: Callable[[], float] = time.time, log: Any = None):
        self.path = Path(path) if path else None
        self.idle = float(idle)
        self.lifetime = float(lifetime)
        self.clock = clock
        self.log = log
        self._live: dict[str, list[float]] = {}   # digest -> [created, seen]
        self._saved_at = 0.0
        self._write_failed = False
        self._load()

    def __len__(self) -> int:
        return len(self._live)

    # -- the store -----------------------------------------------------------
    def _load(self) -> None:
        if self.path is None:
            return
        try:
            fd = os.open(str(self.path), os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
        except OSError:
            return
        try:
            info = os.fstat(fd)
            if info.st_uid != os.getuid() or info.st_mode & 0o077:
                # Not this account's, or readable by others: not trusted.
                return
            with os.fdopen(fd, "r", encoding="utf-8") as handle:
                fd = -1
                data = json.load(handle)
        except (OSError, ValueError):
            return
        finally:
            if fd >= 0:
                os.close(fd)
        if not isinstance(data, dict):
            return
        for key, value in data.items():
            if (isinstance(key, str) and len(key) == 64
                    and isinstance(value, list) and len(value) == 2
                    and all(isinstance(v, (int, float)) for v in value)):
                self._live[key] = [float(value[0]), float(value[1])]

    def _save(self) -> None:
        if self.path is None:
            return
        now = self.clock()
        tmp = self.path.with_name(f".{self.path.name}.{secrets.token_hex(4)}")
        try:
            fd = os.open(str(tmp), os.O_WRONLY | os.O_CREAT | os.O_EXCL
                         | getattr(os, "O_NOFOLLOW", 0), 0o600)
            with os.fdopen(fd, "w", encoding="utf-8") as handle:
                json.dump(self._live, handle)
            os.replace(str(tmp), str(self.path))
            self._saved_at = now
            self._write_failed = False
        except OSError as exc:
            try:
                os.unlink(str(tmp))
            except OSError:
                pass
            if not self._write_failed and self.log is not None:
                self.log.warning(
                    "[delfin] dashboard logins are kept in memory only "
                    "(the session store %s could not be written: %s); a "
                    "restart of the server ends them", self.path, exc)
            self._write_failed = True

    def _prune(self, now: float) -> None:
        for key in [k for k, (created, seen) in self._live.items()
                    if now - seen > self.idle or now - created > self.lifetime]:
            del self._live[key]
        if len(self._live) > MAX_SESSIONS:
            # A session never used again after it was made goes first: a
            # client without a cookie jar (curl with the token header) gets
            # a new one per request, and those must not push out the
            # browser that is logged in. Then the least recently used.
            order = sorted(self._live, key=lambda k: (
                self._live[k][1] > self._live[k][0], self._live[k][1]))
            for key in order[:len(self._live) - MAX_SESSIONS]:
                del self._live[key]

    # -- the answers ---------------------------------------------------------
    def open(self) -> str:
        """A new session; returns its id (for the cookie)."""
        now = self.clock()
        sid = secrets.token_urlsafe(24)
        self._live[_digest(sid)] = [now, now]
        self._prune(now)
        self._save()
        return sid

    def touch(self, sid: Any) -> bool:
        """True if *sid* may continue; records the use."""
        if not isinstance(sid, str) or not 16 <= len(sid) <= 128:
            return False
        key = _digest(sid)
        entry = self._live.get(key)
        if entry is None:
            return False
        now = self.clock()
        created, seen = entry
        if now - seen > self.idle or now - created > self.lifetime:
            del self._live[key]
            self._save()
            return False
        entry[1] = now
        if now - self._saved_at >= _SAVE_EVERY:
            self._save()
        return True

    def close(self, sid: Any) -> None:
        """End *sid* (logout). Unknown ids are ignored."""
        if not isinstance(sid, str):
            return
        if self._live.pop(_digest(sid), None) is not None:
            self._save()


def store_path(runtime_dir: str, port: Any) -> Path:
    """Where the server on *port* keeps its sessions."""
    return Path(runtime_dir) / f"{STORE_PREFIX}{int(port)}.json"


# -- the browser side ------------------------------------------------------------

def logout_control_html() -> str:
    """The header's logout control: a link to Jupyter's ``/logout``.

    It asks first, because leaving the page ends a session that is not
    kept (as closing the window does).
    """
    return (
        '<a class="delfin-logout" href="#" role="button" '
        'title="Log out of this dashboard on this server. A kept session '
        'keeps running; any other ends as when the window is closed." '
        'style="font-size:12px;color:#555;text-decoration:none;'
        'padding:2px 6px;white-space:nowrap;" '
        'onclick="(function(){var p=window.location.pathname||\'/\';'
        'var i=p.indexOf(\'/voila/\');var b=i>=0?p.slice(0,i+1):\'/\';'
        'if(window.confirm(\'Log out of the dashboard? A kept session keeps '
        'running; any other ends as when the window is closed.\')){'
        # The close beacon (window_close) goes out while the login still
        # stands; sent on the navigation itself it raced the logout and was
        # refused, so the server took the close for a dropped connection.
        'try{window.dispatchEvent(new Event(\'pagehide\'));}catch(e){}'
        'setTimeout(function(){window.location.href=b+\'logout\';},300);'
        '}})();return false;">Log out</a>'
    )


def keepalive_js() -> str:
    """Report user input on :data:`ALIVE_ROUTE`, at most every five minutes.

    A refused report (the login ended: idle, logged out elsewhere, server
    restarted without its store) shows a notice instead of failing later
    on the first image the page loads.
    """
    route = json.dumps(ALIVE_ROUTE.lstrip("/"))
    return (
        "(function(){\n"
        "  if (window.__delfinKeepAlive) { return; }\n"
        "  window.__delfinKeepAlive = true;\n"
        f"  var route = {route}, every = 300000, last = 0, ended = false;\n"
        "  function base(){\n"
        "    var p = window.location.pathname || '/';\n"
        "    var i = p.indexOf('/voila/');\n"
        "    return i >= 0 ? p.slice(0, i + 1) : '/';\n"
        "  }\n"
        "  function notice(){\n"
        "    ended = true;\n"
        "    try {\n"
        "      var d = document.createElement('div');\n"
        "      d.className = 'delfin-login-ended';\n"
        "      d.textContent = 'Your dashboard login has ended (idle, or "
        "logged out). This page keeps its open connection; reload it and "
        "enter the token to continue.';\n"
        "      d.style.cssText = 'position:fixed;left:0;right:0;top:0;"
        "z-index:2147483647;background:#b71c1c;color:#fff;padding:6px 12px;"
        "font:13px sans-serif;';\n"
        "      document.body.appendChild(d);\n"
        "    } catch (e) {}\n"
        "  }\n"
        "  function used(){\n"
        "    var now = Date.now();\n"
        "    if (ended || now - last < every) { return; }\n"
        "    last = now;\n"
        "    try {\n"
        "      fetch(base() + route, {credentials: 'same-origin',\n"
        "        cache: 'no-store', redirect: 'manual'}).then(function(r){\n"
        "        if (r.status === 401 || r.status === 403 ||\n"
        "            r.type === 'opaqueredirect') { notice(); }\n"
        "      }).catch(function(){});\n"
        "    } catch (e) {}\n"
        "  }\n"
        "  ['pointerdown', 'keydown', 'wheel'].forEach(function(t){\n"
        "    window.addEventListener(t, used, {capture: true, passive: true});\n"
        "  });\n"
        "})();\n"
    )
