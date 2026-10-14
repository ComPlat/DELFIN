"""Find the agent service, or start it, for ``delfin-agent tui``.

The service is started detached -- its own session (``setsid``), no
terminal, its stderr in the service log -- so it outlives the window that
started it, the way a tmux server outlives the first ``tmux`` command. It
lives under no lifeline: a TUI started from a dashboard's terminal keeps its
sessions when the dashboard ends. An emergency stop still ends it
(``stop_all``), like every agent process.
"""

from __future__ import annotations

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
    stop` a moment ago) is gone soon, and a new one is started then. One
    that refuses this window's protocol version, or never answers, is
    named, with how to end it -- never stopped from here, its sessions may
    be at work.
    """
    deadline = time.monotonic() + timeout
    hello, why = client.probe()
    if hello is not None:
        return hello
    if why.startswith("it speaks another protocol"):
        raise ServiceDidNotStart(_running_but(why))
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
        if hello is not None:
            return hello
        if why.startswith("it speaks another protocol"):
            raise ServiceDidNotStart(_running_but(why))
        time.sleep(0.1)
    pid = rundir.service_pid()
    if pid:
        raise ServiceDidNotStart(_running_but(why or "it does not answer"))
    raise ServiceDidNotStart(f"the agent service did not answer within "
                             f"{timeout:.0f} s; its log: {service_log()}\n"
                             f"{_log_tail()}")


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
