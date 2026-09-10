"""Keeping a dashboard session — and letting it end, which stays the default.

Voila gives every browser connection a fresh kernel and re-executes the
notebook, so a reload is a restart and closing the window ends the run.
That is the behaviour people expect here and it does not change.

What did not exist is a way to say "not this one". A turn the agent is
part-way through, a half-filled submit form, an opened calculation and
everything parsed behind it live in the KERNEL — 2054 widget
constructions across nineteen tabs, plus the Python state beside them.
None of that has to be saved if the kernel simply does not go away.

So the mechanism is inverted rather than added to: the kernel watches
its own page and shuts ITSELF down when the page stops talking, and the
opt-in says only "do not". Nothing is serialised, nothing is restored,
and the work does not grow when somebody adds a twentieth tab.

Three properties this module is built around:

* **The default is unchanged.** No opt-in means the session ends when
  the window does, as it does today.
* **The watchdog arms on the FIRST beat, never before.** A page that
  cannot send heartbeats — an older template, a frontend that failed to
  load its scripts — must not be read as a page that went away.
* **Keeping a session is visible or it is a leak.** ``describe()`` is
  what the status strip and the terminal print, so an armed session
  cannot be invisible.
"""

from __future__ import annotations

import os
import threading
import time
from typing import Any, Callable, Optional


# How long a page may stay silent before the kernel calls it gone. The
# page beats every ~10s, so this tolerates several missed beats over a
# slow link rather than tearing a working session down.
GRACE_SECONDS = 90.0

# How often the watchdog looks. Cheap, and the grace above is what
# decides; this only bounds how late the teardown is.
_POLL_SECONDS = 5.0


_lock = threading.RLock()
_roots: list[Any] = []
_keep_alive = False
_last_beat: Optional[float] = None      # None until the page beats once
_watchdog: Optional[threading.Thread] = None
_stop = threading.Event()
_session_name = ""


# ---------------------------------------------------------------------------
# What the dashboard registers, and what a resume re-displays
# ---------------------------------------------------------------------------

def register_root(*widgets: Any) -> None:
    """Remember the widgets ``create_dashboard`` displayed, in order.

    A resume does not rebuild anything: the widget objects still exist in
    this kernel with all their values and callbacks, so re-displaying
    them is the whole restore. They are kept here rather than only on the
    context because the resume path has to find them without knowing how
    the dashboard was built.
    """
    with _lock:
        _roots.clear()
        _roots.extend(w for w in widgets if w is not None)


def roots() -> list[Any]:
    with _lock:
        return list(_roots)


def resume() -> bool:
    """Show the existing dashboard in the frontend that is connected now.

    Returns False when there is nothing to show — a kernel that never
    built a dashboard, which is a caller error rather than a state to
    paper over.
    """
    from IPython.display import display

    current = roots()
    if not current:
        return False
    for widget in current:
        display(widget)
    return True


# ---------------------------------------------------------------------------
# The heartbeat, and the teardown it guards
# ---------------------------------------------------------------------------

def beat() -> None:
    """The page is still there. Called from its heartbeat comm."""
    global _last_beat
    with _lock:
        _last_beat = time.monotonic()


def seconds_since_beat() -> Optional[float]:
    """Seconds since the page last spoke, or None if it never has."""
    with _lock:
        if _last_beat is None:
            return None
        return time.monotonic() - _last_beat


def keep_alive(on: bool, *, session_name: str = "") -> None:
    """Arm or disarm the opt-in. Never called by default."""
    global _keep_alive, _session_name
    with _lock:
        _keep_alive = bool(on)
        if session_name:
            _session_name = session_name


def is_kept_alive() -> bool:
    with _lock:
        return _keep_alive


def session_name() -> str:
    with _lock:
        return _session_name


def describe() -> dict:
    """What the status strip and the terminal line are built from.

    A session that is kept alive without saying so is how a machine ends
    up holding ten of them, so everything a caller needs to make that
    visible is here in one place.
    """
    with _lock:
        return {
            "kept_alive": _keep_alive,
            "session_name": _session_name,
            "seconds_since_beat": seconds_since_beat(),
            "watchdog_armed": _last_beat is not None,
            "grace_seconds": GRACE_SECONDS,
            "pid": os.getpid(),
        }


def _default_shutdown() -> None:                       # pragma: no cover
    """End this kernel. Called only when the page is gone and unpinned."""
    os._exit(0)


def start_watchdog(
    *, shutdown: Optional[Callable[[], None]] = None,
    poll_seconds: float = _POLL_SECONDS,
) -> threading.Thread:
    """Watch the page and end the kernel when it goes, unless pinned.

    Idempotent: a second call returns the running thread rather than
    starting a rival that would race it to the shutdown.
    """
    global _watchdog
    with _lock:
        if _watchdog is not None and _watchdog.is_alive():
            return _watchdog
        _stop.clear()

        def _run() -> None:
            while not _stop.wait(poll_seconds):
                if is_kept_alive():
                    continue
                since = seconds_since_beat()
                # None means the page has never beaten. That is a page
                # which cannot, not a page that left, and tearing the
                # kernel down under it would break the default path for
                # everyone whose frontend did not load the script.
                if since is None or since < GRACE_SECONDS:
                    continue
                (shutdown or _default_shutdown)()
                return

        _watchdog = threading.Thread(
            target=_run, name="delfin-session-watchdog", daemon=True)
        _watchdog.start()
        return _watchdog


def stop_watchdog() -> None:
    """For tests, and for a deliberate 'end this session now'."""
    _stop.set()


def _reset_for_tests() -> None:
    global _keep_alive, _last_beat, _watchdog, _session_name
    stop_watchdog()
    with _lock:
        _roots.clear()
        _keep_alive = False
        _last_beat = None
        _watchdog = None
        _session_name = ""


# ---------------------------------------------------------------------------
# The page's end of it
# ---------------------------------------------------------------------------

#: How often the page reports in. Several of these fit inside
#: ``GRACE_SECONDS``, so a slow link drops beats without dropping the
#: session.
BEAT_SECONDS = 10

_HEARTBEAT_CLASS = "delfin-session-heartbeat"

# Written into a hidden text input with the native setter, then announced
# with input+change, which is how this dashboard already gets values from
# JS into Python (see molecule_viewer). A widget's own event handling
# does the rest, so nothing here has to reach for the kernel's comm API —
# which is not addressable from arbitrary script under Voila.
#
# `visibilitychange` matters more than the interval: a browser throttles
# timers in a hidden tab to once a minute or worse, which would look
# exactly like a window that closed. The beat is therefore sent
# immediately whenever the tab comes back, and the grace above is wide
# enough to cover a throttled interval in between.
_HEARTBEAT_JS = """
(function () {
  var CLASS = '%(cls)s';
  var EVERY = %(every)d * 1000;
  function input() {
    var host = document.querySelector('.' + CLASS);
    return host ? host.querySelector('input, textarea') : null;
  }
  function beat() {
    var el = input();
    if (!el) return;
    var proto = (el.tagName === 'TEXTAREA')
      ? window.HTMLTextAreaElement.prototype
      : window.HTMLInputElement.prototype;
    var setter = Object.getOwnPropertyDescriptor(proto, 'value');
    var next = String(Date.now());
    if (setter && setter.set) setter.set.call(el, next);
    else el.value = next;
    el.dispatchEvent(new Event('input', {bubbles: true}));
    el.dispatchEvent(new Event('change', {bubbles: true}));
  }
  if (window.__delfinHeartbeat) clearInterval(window.__delfinHeartbeat);
  window.__delfinHeartbeat = setInterval(beat, EVERY);
  document.addEventListener('visibilitychange', function () {
    if (!document.hidden) beat();
  });
  beat();
})();
""" % {"cls": _HEARTBEAT_CLASS, "every": BEAT_SECONDS}


def heartbeat_js() -> str:
    """The script the page runs to report in."""
    return _HEARTBEAT_JS


def build_heartbeat_widget():
    """A hidden field the page writes to, wired to :func:`beat`.

    Returns the widget. Every change is a beat; the value itself is a
    timestamp nobody reads, so a beat costs one small widget message.
    """
    import ipywidgets as widgets

    field = widgets.Text(value="", layout=widgets.Layout(display="none"))
    field.add_class(_HEARTBEAT_CLASS)
    field.observe(lambda _change: beat(), names="value")
    return field


# ---------------------------------------------------------------------------
# Where to come back to
# ---------------------------------------------------------------------------

#: The query key that names which session an address comes back to.
#: It rides on the URL because that is the only thing the kernel manager
#: is handed -- the request URL reaches it through the environment Voila
#: builds for the kernel.
SESSION_QUERY_KEY = "session"

#: Set by the launcher: where the one-cell resume notebook was staged,
#: as a path the server serves. Deliberately NOT the dashboard's own
#: address -- that one renders and EXECUTES the full notebook, which
#: against a live kernel would run all nineteen tabs a second time.
RESUME_PATH_ENV = "DELFIN_RESUME_URL_PATH"


def _request_url() -> str:
    """The URL this kernel was opened with, as Voila recorded it."""
    return os.environ.get("VOILA_REQUEST_URL", "")


def resume_url(name: str = "", *, request_url: str = "",
               resume_path: str = "") -> str:
    """The address that comes back to this session, or "" if unknown.

    Built from the request that started the kernel, so it carries the
    host, the port and the server token the user already authenticated
    with. That makes it exactly as sensitive as the address they are
    looking at — it belongs in a terminal or a copy button, not in a
    ticket.
    """
    from urllib.parse import parse_qsl, urlencode, urlsplit, urlunsplit

    base = request_url or _request_url()
    who = name or session_name()
    path = resume_path or os.environ.get(RESUME_PATH_ENV, "")
    if not base or not who or not path:
        return ""
    try:
        parts = urlsplit(base)
        if not parts.netloc:
            return ""
        # Keep whatever the current address carries -- the token above
        # all -- and add the session on top, replacing any stale one.
        query = [(k, v) for k, v in parse_qsl(parts.query)
                 if k != SESSION_QUERY_KEY]
        query.append((SESSION_QUERY_KEY, who))
        return urlunsplit((
            parts.scheme or "http", parts.netloc, path,
            urlencode(query), "",
        ))
    except ValueError:
        return ""


def default_session_name() -> str:
    """A name a person can recognise in a list of running sessions.

    Host first because these are read on a machine that may host several,
    then a short random tail so two dashboards opened the same minute do
    not collide.
    """
    import socket
    import uuid

    host = (socket.gethostname() or "delfin").split(".")[0]
    return f"{host}-{uuid.uuid4().hex[:4]}"


# ---------------------------------------------------------------------------
# The strip that makes it visible
# ---------------------------------------------------------------------------

_STRIP_CSS = """
<style>
.delfin-session-strip { display:flex; align-items:center; gap:8px; }
.delfin-session-dot {
  width:8px; height:8px; border-radius:50%; display:inline-block;
  background:#c8ccd4;
}
.delfin-session-dot.on { background:#4b9e5f; }
.delfin-session-note { font-size:12px; color:#5a6270; }
.delfin-session-note code {
  font-size:11px; background:#f2f4f7; padding:1px 4px; border-radius:3px;
}
</style>
"""


def _strip_html(armed: bool, name: str, url: str) -> str:
    """What the strip says. Split out so it can be read without a kernel."""
    dot = '<span class="delfin-session-dot%s"></span>' % (" on" if armed else "")
    if not armed:
        return (
            f'<div class="delfin-session-strip">{dot}'
            '<span class="delfin-session-note">Sitzung endet beim Schlie&szlig;en'
            '</span></div>'
        )
    where = (f' &middot; zur&uuml;ck &uuml;ber <code>{url}</code>' if url else "")
    return (
        f'<div class="delfin-session-strip">{dot}'
        f'<span class="delfin-session-note">L&auml;uft weiter als '
        f'<code>{name}</code>{where}</span></div>'
    )


def build_status_strip():
    """The control that arms the opt-in, and shows it while it is armed.

    Two things it must do, and the second is the one that keeps this
    honest: arm the session, and never let an armed one be invisible. A
    kept session that says nothing is how a machine ends up holding ten
    of them, so the dot stays in the header rather than hiding in a menu.

    The label states what does NOT come back. Scroll positions and open
    HTML panels live in the browser, not in the kernel, and a control
    that promised otherwise would be worse than none.
    """
    import ipywidgets as widgets

    css = widgets.HTML(_STRIP_CSS)
    note = widgets.HTML(_strip_html(False, "", ""))
    toggle = widgets.ToggleButton(
        value=False,
        description="Offen halten",
        tooltip=(
            "Sitzung nach dem Schließen des Browsers weiterlaufen "
            "lassen. Der Agent arbeitet weiter, Formulare und geöffnete "
            "Rechnungen bleiben. Scroll-Positionen kommen nicht zurück."
        ),
        icon="thumb-tack",
        layout=widgets.Layout(width="130px"),
    )

    def _on_toggle(change):
        if change.get("name") != "value":
            return
        on = bool(change.get("new"))
        name = session_name() or default_session_name()
        keep_alive(on, session_name=name)
        # The record is what a later request follows back to this kernel,
        # and disarming has to remove it: a landing page that offers a
        # session which is gone is worse than one that offers nothing.
        if on:
            write_record(name)
        else:
            drop_record(name)
        note.value = _strip_html(on, name, resume_url(name) if on else "")
        if on:
            announce(name)

    toggle.observe(_on_toggle, names="value")
    strip = widgets.HBox(
        [css, toggle, note],
        layout=widgets.Layout(align_items="center", gap="8px"),
    )
    return strip


def announce(name: str = "") -> str:
    """Print where to come back to, on the terminal that runs the server.

    That terminal is usually inside tmux, which is exactly where somebody
    looks tomorrow after `tmux attach`. Printing to stdout from a kernel
    reaches the server's log, not the browser, which is the point.
    """
    who = name or session_name()
    url = resume_url(who)
    line = (
        f"[delfin] Sitzung \"{who}\" bleibt bestehen.\n"
        f"         Zurück:  {url or '(Adresse unbekannt)'}\n"
        f"         Beenden: im Dashboard, oder Ctrl+C hier"
    )
    print(line, flush=True)
    return line


# ---------------------------------------------------------------------------
# The record the server reads to find this kernel again
# ---------------------------------------------------------------------------

#: Where a kept session announces itself. One small JSON file per armed
#: session, removed when it is disarmed or the kernel ends, so the
#: directory is also the answer to "what is running".
RECORD_DIR = os.path.join(
    os.path.expanduser("~"), ".delfin", "kept_sessions")


def kernel_id() -> str:
    """This kernel's id, as the server knows it.

    Taken from the connection file ipykernel was started with —
    ``.../kernel-<id>.json`` — because a kernel is not told its own id
    any other way, and the server's registry is keyed on exactly that.
    Empty outside a kernel, which is every test and every CLI call.
    """
    try:
        from ipykernel.connect import get_connection_file

        stem = os.path.basename(str(get_connection_file() or ""))
    except Exception:
        return ""
    if not stem.startswith("kernel-") or not stem.endswith(".json"):
        return ""
    return stem[len("kernel-"):-len(".json")]


def record_path(name: str, *, root: str = "") -> str:
    return os.path.join(root or RECORD_DIR, f"{name}.json")


def write_record(name: str = "", *, root: str = "", kid: str = "") -> str:
    """Announce this session so a later request can find it.

    Returns the path written, or "" when there is nothing to announce —
    outside a kernel there is no id, and a record without one would send
    the server looking for a kernel that does not exist.
    """
    import json

    who = name or session_name()
    ident = kid or kernel_id()
    if not who or not ident:
        return ""
    directory = root or RECORD_DIR
    try:
        os.makedirs(directory, exist_ok=True)
        path = record_path(who, root=directory)
        payload = {
            "session_name": who,
            "kernel_id": ident,
            "pid": os.getpid(),
            "started_at": time.time(),
            "request_url": _request_url(),
        }
        tmp = path + ".tmp"
        with open(tmp, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=1)
        os.replace(tmp, path)
        os.chmod(path, 0o600)
        return path
    except OSError:
        return ""


def drop_record(name: str = "", *, root: str = "") -> bool:
    """Take the announcement back. Disarming must leave nothing behind,
    or the landing page offers a session that is gone."""
    who = name or session_name()
    if not who:
        return False
    try:
        os.remove(record_path(who, root=root))
        return True
    except OSError:
        return False


def list_records(*, root: str = "") -> list[dict]:
    """Every announced session, newest first.

    The landing page and the settings list are built from this, which is
    why a stale or unreadable file is skipped rather than raised: one bad
    record must not hide the others.
    """
    import json

    directory = root or RECORD_DIR
    out: list[dict] = []
    try:
        names = sorted(os.listdir(directory))
    except OSError:
        return out
    for entry in names:
        if not entry.endswith(".json"):
            continue
        try:
            with open(os.path.join(directory, entry), encoding="utf-8") as h:
                data = json.load(h)
        except (OSError, ValueError):
            continue
        if data.get("session_name") and data.get("kernel_id"):
            out.append(data)
    out.sort(key=lambda d: float(d.get("started_at") or 0), reverse=True)
    return out
