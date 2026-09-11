"""Keeping a dashboard session — and letting it end, which stays the default.

Voila gives every browser connection a fresh kernel and re-executes the
notebook, so a reload is a restart and closing the window ends the run.
That is the behaviour people expect here and it does not change.

What did not exist is a way to say "not this one". A turn the agent is
part-way through, a half-filled submit form, an opened calculation and
everything parsed behind it live in the KERNEL — 2054 widget
constructions across nineteen tabs, plus the Python state beside them.
None of that has to be saved if the kernel simply does not go away.

So nothing is serialised and nothing is restored. The server ends a
kernel that has had no window for a grace period (see
``resume_server``), and the opt-in here says only "do not": it writes a
record naming this kernel, which the server reads as the exemption and
a later request follows back to the kernel. A resume re-displays the
widget objects that never left.

Two properties this module is built around:

* **The default is unchanged.** No opt-in means the session ends when
  the window does, as it does today.
* **Keeping a session is visible or it is a leak.** The strip in the
  header, the record on disk and the line on the terminal are the same
  fact stated three times, so an armed session cannot be invisible.
"""

from __future__ import annotations

import os
import threading
import time
from typing import Any, Callable, Optional


_lock = threading.RLock()
_roots: list[Any] = []
_bootstrap: Optional[Callable[[], str]] = None   # page scripts
_bootstrap_before: Any = None                     # the root they run ahead of
_hooks: list[Callable[[], None]] = []             # after the roots
_keep_alive = False
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


def register_bootstrap(getter: Optional[Callable[[], str]], *,
                       before: Any = None) -> None:
    """The scripts a page needs before any widget on it is worth looking at.

    The widget objects come back by themselves; the JavaScript around
    them does not. The bundled 3Dmol, every tab's startup script and the
    editor's lazily sent bootstraps were sent once, into an Output widget
    that the next script cleared -- so a resumed page had the widgets and
    none of the code that makes a viewer a viewer. The getter returns
    what has to run FIRST, and resume() emits it as the cell's own output
    ahead of the roots, so it has run by the time any Output widget
    replays a viewer into the page.
    """
    global _bootstrap, _bootstrap_before
    with _lock:
        _bootstrap = getter
        _bootstrap_before = before


def on_resume(hook: Callable[[], None]) -> None:
    """Something to do once the page is back, after the roots are shown.

    For state that lives only in the browser -- a drawing inside an
    editor's frame, a text view built by script -- and has to be pushed
    again from what the kernel remembers of it.
    """
    with _lock:
        if hook not in _hooks:
            _hooks.append(hook)


def why_not_resumed(*, root: str = "", request_url: str = "") -> str:
    """One sentence on why a resume request landed in a fresh kernel.

    Printed by the resume notebook when resume() had nothing to show.
    "This session is gone" was all it said for one night, and the
    person reading it on a cluster could not tell a server restart from
    a missing record from a Voila that never passed the address on.
    """
    url = request_url or _request_url()
    if not url:
        return ("The server did not see this as a return: Voila passed no "
                "request address to the kernel (VOILA_REQUEST_URL is unset), "
                "so the session name never reached the kernel manager.")
    from delfin.dashboard.resume_server import requested_session
    name = requested_session(url)
    if not name:
        return ("The address carries no session= parameter, so this is an "
                "ordinary page, not a return.")
    directory = root or RECORD_DIR
    record = next((r for r in list_records(root=root)
                   if r.get("session_name") == name), None)
    mine = kernel_id()
    if record is None:
        # What THIS process sees on disk, spelled out: the server said
        # "no record" for a session the terminal had just called kept,
        # and nothing told the reader whether the file was ever there.
        try:
            names = sorted(n for n in os.listdir(directory) if n.endswith(".json"))
            listing = (f"the directory holds {len(names)} record(s): "
                       + ", ".join(names[:8]) if names
                       else "the directory exists and holds no record")
        except FileNotFoundError:
            listing = "the directory does not exist"
        except OSError as exc:
            listing = f"the directory cannot be listed ({type(exc).__name__}: {exc})"
        home = os.environ.get("HOME", "")
        return (f"No record for session \"{name}\" under {directory}: it was "
                "ended, or it was kept by another account or on another "
                f"machine. As seen from this kernel: {listing}; HOME={home}; "
                f"kernel {mine[:8] or '?'}.")
    wanted = str(record.get("kernel_id") or "")
    return (f"A record for \"{name}\" exists (kernel {wanted[:8]}), but the "
            f"server started a fresh kernel ({mine[:8]}) for this page: it no "
            "longer runs that kernel -- a server restart, or the kernel ended. "
            "The record is dropped on the server side.")


def resume_hooks() -> list[Callable[[], None]]:
    with _lock:
        return list(_hooks)


def resume() -> bool:
    """Show the existing dashboard in the frontend that is connected now.

    Returns False when there is nothing to show — a kernel that never
    built a dashboard, which is a caller error rather than a state to
    paper over.

    Order matters and is the whole reason this is more than a loop of
    display(): page scripts first, as plain cell output, then the widget
    roots, then the hooks that push browser-only state back. A hook that
    fails is reported and skipped; one tab's script must not keep the
    other eighteen off the page.
    """
    from IPython.display import Javascript, display

    current = roots()
    if not current:
        return False
    with _lock:
        getter, before = _bootstrap, _bootstrap_before
    script = ""
    if getter is not None:
        try:
            script = getter() or ""
        except Exception as exc:                     # noqa: BLE001
            print(f"[delfin] page scripts could not be collected: {exc}")
    # Where the scripts go: ahead of the root they were registered
    # against, or ahead of everything. At startup the bundle sits at the
    # top of the body, below a header that is already on the page, and
    # a script that reaches for the header finds it there. The same
    # order here, so what worked at startup works again.
    at = current.index(before) if before in current else 0
    for index, widget in enumerate(current):
        if index == at and script.strip():
            display(Javascript(script))
        display(widget)
    for hook in resume_hooks():
        try:
            hook()
        except Exception as exc:                     # noqa: BLE001
            print(f"[delfin] a resume hook failed and was skipped: {exc}")
    return True


# ---------------------------------------------------------------------------
# The opt-in
# ---------------------------------------------------------------------------

def keep_alive(on: bool, *, session_name: str = "") -> None:
    """Arm or disarm the opt-in. Never called by default.

    The state here is what the strip and the terminal describe; what the
    SERVER honours is the record on disk, written and dropped beside
    this by the control.
    """
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
            "pid": os.getpid(),
        }


def _reset_for_tests() -> None:
    global _keep_alive, _session_name, _bootstrap, _bootstrap_before
    with _lock:
        _roots.clear()
        _hooks.clear()
        _bootstrap = None
        _bootstrap_before = None
        _keep_alive = False
        _session_name = ""


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
.delfin-session-strip { display:flex; align-items:center; gap:8px; white-space:nowrap; }
.delfin-session-dot {
  width:8px; height:8px; border-radius:50%; display:inline-block;
  background:#c8ccd4;
}
.delfin-session-dot.on { background:#4b9e5f; }
.delfin-session-dot.failed { background:#d64545; }
.delfin-session-link { line-height:0; text-decoration:none; }
</style>
"""


def _strip_html(armed: bool, name: str, url: str) -> str:
    """What the strip shows: a dot. Nothing else.

    The state is in the dot's colour and its tooltip; the return address
    is where the dot links to when the session is kept, and on the
    server's terminal. A sentence beside the dot, and the address spelled
    out in the header, were asked to go (2026-09-11): the header is for
    the work, not for the session's paperwork.
    """
    if not armed:
        return (
            '<div class="delfin-session-strip">'
            '<span class="delfin-session-dot" '
            'title="Session ends when the window closes"></span></div>'
        )
    if url:
        # The title sits on the dot, which is what a pointer hovers; the
        # anchor around it is where a click goes.
        title = f"Kept as {name} &middot; click for the return address"
        return (
            '<div class="delfin-session-strip">'
            f'<a class="delfin-session-link" href="{url}">'
            f'<span class="delfin-session-dot on" title="{title}"></span></a></div>'
        )
    return (
        '<div class="delfin-session-strip">'
        f'<span class="delfin-session-dot on" title="Kept as {name} &middot; '
        'return address unknown"></span></div>'
    )


def _failed_strip_html(reason: str) -> str:
    """The dot when the session could not be kept: red, and the reason
    in its tooltip."""
    safe = (reason or "unknown reason").replace('"', "&quot;")
    return (
        '<div class="delfin-session-strip">'
        f'<span class="delfin-session-dot failed" title="Could not keep the '
        f'session: {safe}"></span></div>'
    )


def announce_failure(name: str, reason: str) -> None:
    """Say on the server's terminal that the session could NOT be kept."""
    line = (
        f"[delfin] Session \"{name}\" could NOT be kept: {reason or 'unknown reason'}\n"
        f"         The return address would lead nowhere; fix the record "
        f"directory and switch the toggle on again."
    )
    out = _server_stdout() if kernel_id() else None
    if out is not None:
        try:
            with out:
                out.write(line + "\n")
                out.flush()
            return
        except OSError:
            pass
    print(line)


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
        description="Keep session",
        tooltip=(
            "Keep the session running after the browser closes. The agent "
            "keeps working, forms and open calculations stay. Scroll "
            "positions do not come back."
        ),
        icon="thumb-tack",
        # flex 0 0 auto: the header is a flex row, and a control that may
        # shrink is one that renders twenty pixels wide the moment the
        # row is full -- seen in a browser, clickable and unreadable.
        layout=widgets.Layout(width="130px", flex="0 0 auto"),
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
            if not write_record(name) and last_write_error():
                # Nothing to come back to: say so where the address
                # would have been, and do not leave the control armed.
                reason = last_write_error()
                keep_alive(False, session_name=name)
                # Disarm first: the observer re-enters with on=False and
                # draws the grey dot, and the red one has to come after.
                toggle.value = False
                note.value = _failed_strip_html(reason)
                announce_failure(name, reason)
                return
        else:
            drop_record(name)
        note.value = _strip_html(on, name, resume_url(name) if on else "")
        if on:
            announce(name)

    toggle.observe(_on_toggle, names="value")
    strip = widgets.HBox(
        [css, toggle, note],
        layout=widgets.Layout(align_items="center", gap="8px",
                              flex="0 0 auto"),
    )
    return strip


def _server_stdout():
    """The server's own stdout: the terminal, or the log the launcher keeps.

    A kernel's print never gets there. ipykernel captures Python-level AND
    fd-level stdout and forwards both to the frontend as stream output,
    which under Voila is nowhere -- a server log stayed at zero bytes
    through two kept sessions before this was measured. The server that
    spawned this kernel is its parent, and on Linux its stdout is a path.
    """
    for fd in (1, 2):
        try:
            return open(f"/proc/{os.getppid()}/fd/{fd}", "w",
                        encoding="utf-8", errors="replace")
        except OSError:
            continue
    return None


def announce(name: str = "") -> str:
    """Say where to come back to, on the terminal that runs the server.

    That terminal is usually inside tmux, which is exactly where somebody
    looks tomorrow after `tmux attach`.
    """
    who = name or session_name()
    url = resume_url(who)
    line = (
        f"[delfin] Session \"{who}\" is kept.\n"
        f"         Return:  {url or '(address unknown)'}\n"
        f"         End:     in the dashboard, or Ctrl+C here"
    )
    # Only a kernel has to reach past its own captured stdout; anywhere
    # else -- a test, a CLI -- the parent is not the server, and print
    # is what the caller expects to see.
    out = _server_stdout() if kernel_id() else None
    if out is not None:
        with out:
            out.write(line + "\n")
            out.flush()
    else:
        print(line, flush=True)
    return line


# ---------------------------------------------------------------------------
# The record the server reads to find this kernel again
# ---------------------------------------------------------------------------

#: Where a kept session announces itself. One small JSON file per armed
#: session, removed when it is disarmed or the kernel ends, so the
#: directory is also the answer to "what is running".
RECORD_DIR = os.environ.get("DELFIN_SESSION_RECORD_DIR") or os.path.join(
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

    global _last_write_error
    who = name or session_name()
    ident = kid or kernel_id()
    if not who or not ident:
        # Outside a kernel there is nothing to announce; that is not a
        # failure of the record directory, and the control may still
        # arm the in-memory state a test or a CLI looks at.
        _last_write_error = ""
        return ""
    directory = root or RECORD_DIR
    try:
        os.makedirs(directory, exist_ok=True)
        path = record_path(who, root=directory)
        payload = {
            "session_name": who,
            "kernel_id": ident,
            "pid": os.getpid(),
            "host": _hostname(),
            "started_at": time.time(),
            "request_url": _request_url(),
        }
        tmp = path + ".tmp"
        with open(tmp, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=1)
        os.replace(tmp, path)
    except OSError as exc:
        # A record that could not be written used to fail in silence,
        # and the control then announced a return address that led
        # nowhere -- seen on a cluster on 2026-09-11, where the server
        # found no record for a session the terminal had just called
        # kept. The reason is kept for the announcement.
        _last_write_error = f"{type(exc).__name__}: {exc} (record dir {directory})"
        return ""
    try:
        # Owner-only, but a file system that refuses the mode still has
        # the record: the write is what counts.
        os.chmod(path, 0o600)
    except OSError:
        pass
    _last_write_error = ""
    return path


_last_write_error = ""


def last_write_error() -> str:
    """Why the last write_record returned "", or "" when it succeeded."""
    return _last_write_error


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

    def _started(record: dict) -> float:
        # A record is written by a kernel and read by a server, possibly
        # of a different vintage. The promise above is that one bad file
        # does not hide the others, and a sort key that raises breaks
        # exactly that -- found by a test that put a string in this
        # field.
        try:
            return float(record.get("started_at") or 0)
        except (TypeError, ValueError):
            return 0.0

    out.sort(key=_started, reverse=True)
    return out


# ---------------------------------------------------------------------------
# Landing on a dashboard while another session is still running
# ---------------------------------------------------------------------------

def _hostname() -> str:
    import socket

    try:
        return (socket.gethostname() or "").split(".")[0]
    except Exception:
        return ""


def _record_is_dead(record: dict) -> bool:
    """A record whose kernel process is gone, as far as this host can tell.

    Only judged on the host that wrote it: a home directory may be shared
    between machines, and a pid means nothing on another one. Elsewhere
    the record is left alone and the resume itself decides.
    """
    host = str(record.get("host") or "")
    if not host or host != _hostname():
        return False
    try:
        pid = int(record.get("pid") or 0)
    except (TypeError, ValueError):
        return False
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return True
    except OSError:
        return False
    return False


def other_sessions(*, root: str = "", exclude_kernel: str = "") -> list[dict]:
    """Kept sessions that are not the one being looked at.

    The current kernel is excluded by id rather than by name: a resume
    renders into the running kernel, so the page you are on IS the kept
    session and offering to go back to it would be a loop.

    A record whose process is gone is dropped rather than offered. After
    a machine restart every kept session is such a record, and a landing
    page that lists them all until each is clicked away is one nobody
    trusts.
    """
    mine = exclude_kernel or kernel_id()
    out = []
    for record in list_records(root=root):
        if str(record.get("kernel_id") or "") == mine:
            continue
        if _record_is_dead(record):
            drop_record(str(record.get("session_name") or ""), root=root)
            continue
        out.append(record)
    return out


def _banner_html(records: list[dict]) -> str:
    """What a returning visitor is told. Split out so it reads without a
    kernel, and so the wording is testable."""
    if not records:
        return ""
    rows = []
    for record in records:
        name = str(record.get("session_name") or "")
        url = resume_url(name, request_url=str(record.get("request_url") or ""))
        started = record.get("started_at")
        try:
            age = max(0.0, time.time() - float(started))
            hours = age / 3600.0
            when = (f"for {hours:.0f}&nbsp;h" if hours >= 1
                    else f"for {age / 60.0:.0f}&nbsp;min")
        except (TypeError, ValueError):
            when = ""
        link = (f'<a href="{url}">re-enter</a>' if url
                else '<span style="color:#8a919e">address unknown</span>')
        rows.append(
            f'<li style="margin:2px 0"><code>{name}</code>'
            f'{" &middot; " + when if when else ""} &middot; {link}</li>')
    return (
        '<div style="border:1px solid #d7dbe2; border-left:3px solid #4b9e5f;'
        ' background:#f7f9fb; padding:8px 12px; margin:0 0 8px 0;'
        ' border-radius:4px; font-size:13px">'
        '<b>A session is still running.</b> '
        '<span style="color:#5a6270">This window is new &mdash; '
        'you can continue where you left off.</span>'
        f'<ul style="margin:6px 0 0 18px; padding:0">{"".join(rows)}</ul>'
        '</div>'
    )


def build_returning_banner(*, root: str = "", exclude_kernel: str = ""):
    """A banner offering the running session, or None when there is none.

    Shown where the user lands rather than behind a route of its own:
    they go to the address they always go to, and it tells them. That is
    also why nothing here needs a server extension.
    """
    import ipywidgets as widgets

    records = other_sessions(root=root, exclude_kernel=exclude_kernel)
    html = _banner_html(records)
    if not html:
        return None
    return widgets.HTML(html)
