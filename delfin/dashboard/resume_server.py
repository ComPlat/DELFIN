"""Answering the address a kept session left behind.

The dashboard's own URL cannot do this. It goes through Voila's
renderer, which EXECUTES the notebook — so opening it against a live
kernel would run all nineteen tabs a second time in a kernel that
already has them, which is the opposite of coming back.

So a resume takes a different route and a different notebook. The
notebook has one cell:

    from delfin.dashboard import session; session.resume()

Voila renders that the way it renders anything, and the one thing that
changes is WHICH kernel it renders into: a kernel manager that,
recognising a resume request, hands back the kernel the session is
already running in instead of starting a fresh one. The cell then
re-displays the widget objects that never left, and the page is the
dashboard again — every value, every callback, the agent's thread still
going.

Nothing is serialised and nothing is restored. That is why this is four
small pieces rather than a persistence layer.
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path
from typing import Any, Optional
from urllib.parse import parse_qs, quote, urlsplit

from delfin.dashboard import session as _session


# ---------------------------------------------------------------------------
# Ending a kernel nobody is looking at
# ---------------------------------------------------------------------------
#
# The server counts the websocket connections to every kernel it runs. A
# dashboard whose last window closed has zero, and that is the whole
# signal: no script in the page, no heartbeat to miss, no restarter to
# fool -- the server shuts the kernel down itself, the way it would for
# any other reason. Keeping a session is the one exemption, read from
# the record the page wrote when the option was switched on.

#: How long a kernel may go without a window before it is ended. The
#: connection drops the instant a tab is closed, so this is entirely
#: about a reload: the new page has to connect before the old kernel is
#: taken away, or every reload would look like a departure.
GRACE_SECONDS = 90.0

#: A kernel that has never had a window is a page still loading. The
#: dashboard takes a while to build; give it longer than a reload.
NEVER_CONNECTED_SECONDS = 600.0

#: How often the server looks. Bounds how late a teardown is.
POLL_SECONDS = 10

GRACE_ENV = "DELFIN_SESSION_GRACE_SECONDS"

#: Set to 1 to keep the server up after its last window has gone and no
#: session is kept. By default it stops then and frees its port: the
#: person who closed the last window has no Ctrl+C left to press.
STAY_UP_ENV = "DELFIN_DASHBOARD_STAY_UP"


def stays_up() -> bool:
    return os.environ.get(STAY_UP_ENV, "").strip().lower() in ("1", "true", "yes", "on")


def grace_seconds(default: float = GRACE_SECONDS) -> float:
    """The grace, with an environment override for a suite that has to
    wait it out for real."""
    try:
        value = float(os.environ.get(GRACE_ENV, ""))
    except (TypeError, ValueError):
        return default
    return value if value > 0 else default


def cull_config_args(grace: Optional[float] = None) -> list[str]:
    """The server flags that switch the culler on.

    The culler runs only when cull_idle_timeout is positive; the value
    itself is not what decides here -- the override below is -- but it
    has to be set, and setting it to the grace keeps the server's own
    log lines truthful.
    """
    seconds = max(1, int(round(grace if grace is not None else grace_seconds())))
    return [
        f"--MappingKernelManager.cull_idle_timeout={seconds}",
        f"--MappingKernelManager.cull_interval={POLL_SECONDS}",
        "--MappingKernelManager.cull_connected=False",
    ]


def kept_kernel_ids(*, root: str = "") -> set[str]:
    """Every kernel a record says is kept."""
    return {str(r.get("kernel_id") or "") for r in _session.list_records(root=root)}


#: The query key that marks a request as a resume, and names which
#: session to resume. It rides on the URL because that is the only thing
#: Voila hands the kernel manager -- the request URL reaches it through
#: the kernel environment it builds.
SESSION_QUERY_KEY = "session"

RESUME_NOTEBOOK_NAME = "delfin_resume.ipynb"


# ---------------------------------------------------------------------------
# The notebook a resume renders
# ---------------------------------------------------------------------------

def resume_notebook_source() -> dict:
    """One cell. Anything more would run in somebody's live kernel."""
    return {
        "cells": [{
            "cell_type": "code",
            "id": "delfin-resume",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": [
                "# Re-display the dashboard this kernel already holds.\n",
                "# Nothing is rebuilt: the widget objects, their values and\n",
                "# their callbacks never went away.\n",
                "from delfin.dashboard import session as _s\n",
                "if not _s.resume():\n",
                "    print('This session is gone.')\n",
                "    print(_s.why_not_resumed())\n",
            ],
        }],
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3", "language": "python",
                "name": "python3",
            },
            "language_info": {"name": "python"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }


def stage_resume_notebook(root_dir: str | os.PathLike) -> str:
    """Write the resume notebook beside the dashboard, and return it.

    Rewritten every start rather than only when missing: it is generated,
    it is one cell, and a stale copy from an older DELFIN would run
    against a kernel from a newer one.
    """
    directory = Path(root_dir).resolve() / "delfin_voila_runtime"
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / RESUME_NOTEBOOK_NAME
    path.write_text(json.dumps(resume_notebook_source(), indent=1),
                    encoding="utf-8")
    # Signed here as well as by the launcher's `jupyter trust`: on a
    # cluster the latter left the server saying "is not trusted" at every
    # return. Trust concerns stored outputs, of which this cell has none,
    # so the warning was noise -- but noise beside a failed return reads
    # as its cause.
    try:
        import nbformat
        from nbformat.sign import NotebookNotary
        nb = nbformat.read(str(path), as_version=4)
        NotebookNotary().sign(nb)
    except Exception:
        pass
    return str(path)


def resume_render_url(root_dir: str | os.PathLike, notebook: str,
                      name: str) -> str:
    """The Voila render URL for the resume notebook, naming the session."""
    rel = Path(notebook).resolve().relative_to(
        Path(root_dir).resolve()).as_posix()
    return f"/voila/render/{quote(rel)}?{SESSION_QUERY_KEY}={quote(name)}"


# ---------------------------------------------------------------------------
# Recognising a resume, from the only thing the kernel manager is given
# ---------------------------------------------------------------------------

def requested_session(request_url: str) -> str:
    """The session a request asks to resume, or "" for an ordinary one."""
    if not request_url:
        return ""
    try:
        query = urlsplit(request_url).query
    except ValueError:
        return ""
    values = parse_qs(query).get(SESSION_QUERY_KEY) or []
    return (values[0] or "").strip() if values else ""


def kernel_for_session(name: str, *, root: str = "") -> str:
    """The kernel id a session announced, or "" if it announced none."""
    if not name:
        return ""
    for record in _session.list_records(root=root):
        if record.get("session_name") == name:
            return str(record.get("kernel_id") or "")
    return ""


def resume_kernel_manager_class(base: type) -> type:
    """A kernel manager that reuses a kept session's kernel.

    Wrapping whatever class is configured rather than a fixed one: Voila
    already builds its own manager on top of the server's, and replacing
    that would take the preheat and pooling behaviour with it.

    The reuse is refused unless the named kernel is still known to this
    manager. A record can outlive its kernel -- a crash, a machine
    restart -- and handing Voila a dead id would render the resume
    notebook into nothing at all.
    """

    class ResumeAwareKernelManager(base):                # type: ignore[misc]

        async def _delfin_existing(self, env: Optional[dict]) -> str:
            # Every branch is logged: the server's terminal is where the
            # person who could not get back is looking, and "this session
            # is gone" on the page told them nothing about which of these
            # it was.
            name = requested_session((env or {}).get("VOILA_REQUEST_URL", ""))
            if not name:
                return ""
            kid = kernel_for_session(name)
            log = getattr(self, "log", None)
            if not kid:
                if log:
                    try:
                        seen = sorted(n for n in os.listdir(_session.RECORD_DIR)
                                      if n.endswith(".json"))
                        listing = ", ".join(seen) or "no record files"
                    except OSError as exc:
                        listing = f"cannot list: {type(exc).__name__}: {exc}"
                    log.warning(
                        "[delfin] return to session %r: no record under %s "
                        "(ended, or kept by another account/machine); the "
                        "server sees there: %s; server HOME=%s. Starting a "
                        "fresh kernel.", name, _session.RECORD_DIR, listing,
                        os.environ.get("HOME", ""))
                return ""
            try:
                known = kid in self
            except Exception:
                known = False
            if not known:
                # The session announced a kernel that is gone. Drop the
                # record here rather than leaving it to mislead the next
                # visitor.
                if log:
                    log.warning(
                        "[delfin] return to session %r: its kernel %s is not "
                        "one this server runs (server restarted, or kernel "
                        "ended); dropping the record and starting a fresh "
                        "kernel.", name, kid[:8])
                _session.drop_record(name)
                return ""
            if log:
                log.info("[delfin] return to session %r: reusing kernel %s.",
                         name, kid[:8])
            return kid

        async def start_kernel(self, *args: Any, **kwargs: Any):
            existing = await self._delfin_existing(kwargs.get("env"))
            if existing:
                return existing
            kid = await super().start_kernel(*args, **kwargs)
            # Unwatched from birth until a window connects. Timed from
            # here so a page that never manages to connect is not kept
            # for the life of the server.
            self._delfin_unwatched()[kid] = (time.monotonic(), False)
            return kid

        # -- who is looking ------------------------------------------------

        def _delfin_unwatched(self) -> dict:
            """kernel id -> (since, had_a_window)."""
            store = getattr(self, "_delfin_unwatched_since", None)
            if store is None:
                store = {}
                self._delfin_unwatched_since = store
            return store

        def notify_connect(self, kernel_id):
            super().notify_connect(kernel_id)
            self._delfin_unwatched().pop(kernel_id, None)

        def notify_disconnect(self, kernel_id):
            super().notify_disconnect(kernel_id)
            try:
                left = int(self._kernel_connections.get(kernel_id, 0))
            except Exception:
                left = 0
            if left <= 0:
                self._delfin_unwatched()[kernel_id] = (time.monotonic(), True)

        def _delfin_seconds_unwatched(self, kernel_id) -> Optional[tuple]:
            """(seconds, had_a_window), or None while a window is there."""
            entry = self._delfin_unwatched().get(kernel_id)
            if entry is None:
                return None
            since, had_window = entry
            return (time.monotonic() - since, had_window)

        # -- the decision ---------------------------------------------------

        async def cull_kernel_if_idle(self, kernel_id):
            """End a kernel that has had no window for the grace, unless
            it is kept. Replaces the server's idle rule outright: a
            dashboard kernel is idle while its agent works in a thread
            and busy with widget traffic while nobody watches, so
            activity says nothing here. Windows do.
            """
            state = self._delfin_seconds_unwatched(kernel_id)
            if state is None:
                return
            seconds, had_window = state
            allowed = grace_seconds() if had_window else NEVER_CONNECTED_SECONDS
            if seconds <= allowed:
                return
            kept = kept_kernel_ids()
            if kernel_id in kept:
                return
            # A warning, not info: the default server log hides info, and
            # a kernel ending is the one event the person who kept a
            # session must be able to see on the terminal.
            self.log.warning(
                "[delfin] ending kernel %s: no window for %ds and not kept "
                "(kept kernels on record: %s).",
                kernel_id[:8], int(seconds),
                ", ".join(k[:8] for k in sorted(kept)) or "none")
            await self.shutdown_kernel(kernel_id)
            self._delfin_stop_if_nothing_is_left()

        # -- the last window -------------------------------------------------

        def _delfin_kernels_alive(self) -> list:
            try:
                return list(self.list_kernel_ids())
            except Exception:
                return list(getattr(self, "_kernels", {}) or {})

        def _delfin_stop_if_nothing_is_left(self) -> None:
            """Stop the server once no kernel runs and no session is kept.

            Closing the last window used to leave the server serving
            nothing, its port taken, until somebody found the terminal
            and pressed Ctrl+C. A kept session is the one reason to stay:
            its kernel is alive and its record says so. STAY_UP_ENV keeps
            the old behaviour for a server meant to outlive its windows.
            """
            if stays_up() or getattr(self, "_delfin_stopping", False):
                return
            if self._delfin_kernels_alive() or kept_kernel_ids():
                return
            log = getattr(self, "log", None)
            if log:
                log.warning(
                    "[delfin] the last window closed and no session is kept: "
                    "stopping the server and freeing its port (keep a session, "
                    "or set %s=1, to have it stay).", STAY_UP_ENV)
            app = getattr(self, "parent", None)
            stop = getattr(app, "stop", None)
            if not callable(stop):
                # A manager without a server above it -- a test double,
                # an embedding that built it by hand -- has nothing to
                # stop; a signal to the process would hit whatever
                # process that is (it hit pytest once).
                return
            self._delfin_stopping = True
            try:
                stop()
            except Exception as exc:
                if log:
                    log.warning("[delfin] could not stop the server: %s", exc)

        async def shutdown_all(self, *args: Any, **kwargs: Any):
            # The server stopping is the one shutdown a kept kernel obeys.
            self._delfin_stopping = True
            return await super().shutdown_all(*args, **kwargs)

        async def shutdown_kernel(self, kernel_id, *args: Any, **kwargs: Any):
            # Voila's page says goodbye when it unloads: a beacon to its
            # shutdown route, which lands here as a plain shutdown of the
            # kernel. Closing the window IS the event a kept session must
            # survive, and for one evening this override honoured the
            # goodbye, dropped the record and ended the kernel -- seven
            # seconds before the return request arrived (cluster,
            # 2026-09-11). A kept kernel ignores the page's goodbye. The
            # server stopping and a restart still go through: the first
            # is the owner's Ctrl+C, the second keeps the id.
            restart = bool(kwargs.get("restart")) or (len(args) >= 2 and bool(args[1]))
            kept_as = next((str(r.get("session_name") or "")
                            for r in _session.list_records()
                            if str(r.get("kernel_id") or "") == str(kernel_id)), "")
            if kept_as and not restart and not getattr(self, "_delfin_stopping", False):
                log = getattr(self, "log", None)
                if log:
                    log.warning("[delfin] kernel %s is kept as session %r; "
                                "ignoring the shutdown request the page sent "
                                "on closing.", str(kernel_id)[:8], kept_as)
                return None
            # Whatever else ends the kernel -- the cull rule, Ctrl+C, a
            # restart -- the record must not outlive it, or the landing
            # page offers a session that is gone.
            for record in _session.list_records():
                if str(record.get("kernel_id") or "") == str(kernel_id):
                    name = str(record.get("session_name") or "")
                    log = getattr(self, "log", None)
                    if log:
                        log.warning("[delfin] kernel %s is shutting down; "
                                    "dropping the record of session %r.",
                                    str(kernel_id)[:8], name)
                    _session.drop_record(name)
            self._delfin_unwatched().pop(kernel_id, None)
            result = await super().shutdown_kernel(kernel_id, *args, **kwargs)
            # The page's goodbye ends an un-kept kernel through here, not
            # through the cull rule -- and for one run the server then sat
            # serving nothing, because the stop lived in the cull path
            # alone. Whatever ended the last kernel, nothing is left.
            if not getattr(self, "_delfin_stopping", False):
                self._delfin_stop_if_nothing_is_left()
            return result

    ResumeAwareKernelManager.__name__ = f"ResumeAware{base.__name__}"
    return ResumeAwareKernelManager


# ---------------------------------------------------------------------------
# The concrete class the server is pointed at
# ---------------------------------------------------------------------------
#
# `--ServerApp.kernel_manager_class` takes an import path, not a factory,
# so the wrapping happens here once against the server's own default. In
# extension mode -- which is how DELFIN runs Voila -- this is the manager
# Voila is handed, so it is the only place a resume can be recognised.


def _server_default_manager() -> type:
    from jupyter_server.services.kernels.kernelmanager import (
        AsyncMappingKernelManager,
    )

    return AsyncMappingKernelManager


ResumeAwareMappingKernelManager = resume_kernel_manager_class(
    _server_default_manager())


#: The environment variable that tells a kernel where its own resume
#: address lives. Set by the launcher, which is the only thing that
#: knows where the notebook was staged relative to the server root.
RESUME_PATH_ENV = "DELFIN_RESUME_URL_PATH"
