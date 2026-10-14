"""Find the agent service, or start it, for ``delfin-agent tui``.

The service is started detached -- its own session (``setsid``), no
terminal, its stderr in the service log -- so it outlives the window that
started it, the way a tmux server outlives the first ``tmux`` command. It
lives under no lifeline: a TUI started from a dashboard's terminal keeps its
sessions when the dashboard ends. An emergency stop still ends it
(``stop_all``), like every agent process.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path

from . import client, eventlog, rundir

#: How long a starting service gets before the window gives up on it.
START_TIMEOUT_S = 20.0

#: How long a service that holds the host but does not listen yet is
#: waited for: it is stopping (its sessions save; the stop is bounded at
#: about 40 s) or starting.
STOP_WAIT_S = 90.0


class ServiceDidNotStart(RuntimeError):
    """The service was started and never answered; the message says why."""


def service_log() -> Path:
    return eventlog.service_dir() / "logs" / "service.log"


def _service_env() -> dict[str, str]:
    from .. import lifeline

    env = dict(os.environ)
    for name in (lifeline.ENV_PID, lifeline.ENV_TICKS):
        env.pop(name, None)
    return env


def start_service(*, idle_exit_s: float | None = None) -> subprocess.Popen:
    from .. import state_paths

    log = service_log()
    state_paths.ensure_dir(log.parent)
    fd = os.open(log, os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o600)
    cmd = [sys.executable, "-m", "delfin.agent.service.supervisor"]
    if idle_exit_s is not None:
        cmd += ["--idle-exit", str(float(idle_exit_s))]
    try:
        return subprocess.Popen(cmd, stdin=subprocess.DEVNULL,
                                stdout=fd, stderr=fd, start_new_session=True,
                                close_fds=True, cwd=str(Path.home()),
                                env=_service_env())
    finally:
        os.close(fd)


def ensure_service(*, timeout: float = START_TIMEOUT_S,
                   idle_exit_s: float | None = None) -> dict:
    """The service's hello: from the one running here, or from a new one.

    A service that holds the host but does not answer yet is waited for:
    one that is starting answers soon, one that is stopping (`service
    stop` a moment ago) is gone soon, and a new one is started then.

    One of an older DELFIN release, or of another protocol version (an
    update), is moved over to this one when nothing in it would be cut --
    no turn running, no question open: it stops (its sessions are saved
    and come back asleep, its windows reconnect) and a new one starts. A
    busy one of another protocol, or one that never answers, is named,
    with how to end it -- never stopped from here; a busy one of the same
    protocol is used as it is, and the note says how to move it over.
    """
    deadline = time.monotonic() + timeout
    hello, why = client.probe()
    # A service being moved over to this release: never the one to use.
    old = 0
    if hello is not None:
        if not _older_release(hello) or not _move_over(hello):
            return hello
        old, hello, why = int(hello.get("pid") or 0), None, "none"
    elif why.startswith("it speaks another protocol"):
        # Its protocol is not this one's: only its roster says whether it
        # is idle. An idle one is moved over; a busy one is named.
        old = rundir.service_pid()
        if not _stop_idle_other_protocol():
            raise ServiceDidNotStart(_running_but(why))
        why = "none"
    proc: subprocess.Popen | None = None
    started = 0.0
    began = time.monotonic()
    while time.monotonic() < deadline:
        holder = rundir.service_pid()
        if proc is None and holder and why == "none":
            # One holds the host and does not listen: it is stopping (its
            # sessions save) or starting. Wait for it, not for a new one.
            deadline = max(deadline, began + STOP_WAIT_S)
        if not holder and (proc is None or proc.poll() is not None) \
                and time.monotonic() - started > 1.0:
            if proc is not None and proc.returncode not in (0, None):
                raise ServiceDidNotStart(
                    f"the agent service exited ({proc.returncode}) before it "
                    f"answered; its log: {service_log()}\n{_log_tail()}")
            # None runs (or the one that did has left): start one. A start
            # that loses the race for the lock exits 0 and is looked past;
            # the winner writes its run file within a second.
            rundir.ensure_runtime_dir()
            proc = start_service(idle_exit_s=idle_exit_s)
            started = time.monotonic()
            # A new one gets its whole start time, also after a wait for
            # an old one to stop.
            deadline = max(deadline, min(started + timeout, began + STOP_WAIT_S + timeout))
        hello, why = client.probe(timeout=1.0)
        if hello is not None and old and int(hello.get("pid") or 0) == old:
            # Still answering on its way out: wait for it to go.
            hello, why = None, "none"
        if hello is not None:
            return hello
        if why.startswith("it speaks another protocol"):
            if not (old and rundir.service_pid() == old):
                raise ServiceDidNotStart(_running_but(why))
            why = "none"                 # the old one, on its way out
        time.sleep(0.1)
    pid = rundir.service_pid()
    if pid:
        raise ServiceDidNotStart(_running_but(why or "it does not answer"))
    raise ServiceDidNotStart(f"the agent service did not answer within "
                             f"{timeout:.0f} s; its log: {service_log()}\n"
                             f"{_log_tail()}")


def _release() -> str:
    try:
        from ... import __version__
        return str(__version__)
    except Exception:
        return ""


def _version_key(text: str):
    try:
        from packaging.version import Version
        return Version(text)
    except Exception:
        return tuple(int(p) if p.isdigit() else 0 for p in str(text).split("."))


def _older_release(hello: dict) -> bool:
    """The running service is an older DELFIN release than this one."""
    theirs, ours = str(hello.get("release") or ""), _release()
    if not ours or theirs == ours:
        return False
    if not theirs:
        return True                     # from before services said so
    try:
        return _version_key(theirs) < _version_key(ours)
    except TypeError:
        return False


def _idle(sessions: list) -> bool:
    """Nothing a restart would cut: no turn running, no question open."""
    return all(isinstance(s, dict) and s.get("state") in ("idle", "sleeping", "ended")
               and not s.get("pending") for s in sessions)


def _move_over(hello: dict) -> bool:
    """An older release that is idle restarts onto this one: its sessions
    are saved and come back asleep; its windows reconnect. A busy one is
    left as it is (said once): its turns must not be cut."""
    pid = int(hello.get("pid") or 0)
    sessions = list(hello.get("sessions") or [])
    if not _idle(sessions):
        print(f"note: the agent service runs DELFIN {hello.get('release') or 'of before'}, "
              f"this is {_release()}; `delfin-agent service restart` moves it over once its "
              "sessions are done", file=sys.stderr)
        return False
    print(f"moving the agent service from DELFIN {hello.get('release') or 'of before'} to "
          f"{_release()} ({len(sessions)} session{'s' * (len(sessions) != 1)} kept)...",
          file=sys.stderr, flush=True)
    try:
        with client.ServiceClient(client="delfin-agent tui") as svc:
            svc.request("service.shutdown", {"restart": True})
    except Exception:
        return False
    # Asked to stop: from here the start loop waits for it to be gone.
    return bool(pid)


def _stop_idle_other_protocol() -> bool:
    """A service of another protocol version (an update): stopped when its
    roster says it is idle -- SIGTERM ends and saves its sessions."""
    import signal

    pid = rundir.service_pid()
    if not pid:
        return False
    try:
        roster = json.loads(rundir.roster_file().read_text())
    except (OSError, ValueError):
        return False
    if not isinstance(roster, dict) or int(roster.get("pid") or 0) != pid \
            or not _idle(list(roster.get("sessions") or [])):
        return False
    print("the agent service is of another DELFIN version and idle: restarting it...",
          file=sys.stderr, flush=True)
    try:
        os.kill(pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    return True


def other_host_notes() -> list[str]:
    """What runs on the other login nodes (their run files and rosters on
    the shared home): a window says where those sessions are."""
    notes = []
    for other in rundir.other_hosts():
        host = str(other.get("host") or "")
        try:
            roster = json.loads(rundir.roster_file(host).read_text())
        except (OSError, ValueError):
            continue
        if not isinstance(roster, dict) or time.time() - float(roster.get("updatedAt") or 0) > 7 * 86400:
            continue
        n = len(roster.get("sessions") or []) + len(roster.get("dormant") or [])
        if n:
            shown = "".join(ch for ch in host if ch.isprintable())[:60]
            notes.append(f"{n} session{'s' * (n != 1)} of yours {'are' if n != 1 else 'is'} on "
                         f"{shown} (another login node): log in there to open "
                         f"{'them' if n != 1 else 'it'}")
    return notes


def main(argv: list[str] | None = None) -> int:
    """`python -m delfin.agent.service.launcher`: make sure the service
    runs (a window whose service vanished starts the next one this way)."""
    from .. import state_paths

    state_paths.scratch_state_from_environment()
    try:
        hello = ensure_service()
    except (ServiceDidNotStart, rundir.UnsafeRuntimeDir) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 3
    print(hello.get("pid", ""))
    return 0


def service_command() -> list[str]:
    """The command that starts the service when none runs: for a window."""
    return [sys.executable, "-m", "delfin.agent.service.launcher"]


def why_in_words(why: str) -> str:
    """probe()'s reason, said: "none" there is a service that holds the
    host and does not listen (yet, or any more)."""
    return "it is not listening: it is starting or stopping" if why in ("", "none") else why


def _running_but(why: str) -> str:
    pid = rundir.service_pid()
    who = f" (pid {pid})" if pid else ""
    return (f"an agent service runs on this host{who}, but {why_in_words(why)}. "
            "`delfin-agent service stop` ends it -- its sessions save and can "
            "be resumed -- then start the window again.")


def _log_tail(lines: int = 12) -> str:
    try:
        text = service_log().read_text(errors="replace")
    except OSError:
        return ""
    return "\n".join(text.splitlines()[-lines:])


def tui_root() -> Path:
    return Path(__file__).resolve().parents[3] / "tui"


def tui_command() -> list[str] | None:
    """How to run the terminal UI from this checkout, or None.

    Until the TUI ships as one binary (phase 7) it runs from source with
    Bun: ``tui/src/main.tsx`` next to the ``delfin`` package, started in
    ``tui/`` -- Bun reads its configuration (the Solid JSX preload) and
    tsconfig from the working directory.
    """
    import shutil

    main = tui_root() / "src" / "main.tsx"
    bun = shutil.which("bun") or str(Path.home() / ".bun" / "bin" / "bun")
    if main.is_file() and Path(bun).is_file():
        return [bun, str(main)]
    return None


if __name__ == "__main__":
    sys.exit(main())
