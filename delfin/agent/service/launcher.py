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
    """The service's hello: from the one running here, or from a new one."""
    hello = client.reachable()
    if hello is not None:
        return hello
    rundir.ensure_runtime_dir()
    proc = start_service(idle_exit_s=idle_exit_s)
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        hello = client.reachable(timeout=1.0)
        if hello is not None:
            return hello
        code = proc.poll()
        if code is not None:
            # Another window may have started one at the same moment, and
            # this one left because the lock was taken: look once more.
            hello = client.reachable(timeout=2.0)
            if hello is not None:
                return hello
            raise ServiceDidNotStart(
                f"the agent service exited ({code}) before it answered; "
                f"its log: {service_log()}\n{_log_tail()}")
        time.sleep(0.1)
    raise ServiceDidNotStart(f"the agent service did not answer within "
                             f"{timeout:.0f} s; its log: {service_log()}\n"
                             f"{_log_tail()}")


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
