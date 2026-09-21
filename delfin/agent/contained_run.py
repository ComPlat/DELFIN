"""Run a command so that nothing it starts outlives it -- on every host.

The process cage (bubblewrap with a PID namespace) gives this where
unprivileged user namespaces work. Many hosts have no such thing: a CI
runner without bubblewrap, an HPC login node with user namespaces turned
off, macOS. There ``subprocess.run`` killed only the direct child at a
timeout, so a pipeline member, an ``&`` job or a ``nohup`` survived the
command and the stop, and the command shared the agent's terminal, where
a TIOCSTI ioctl can type into the user's shell on older kernels.

This is the floor that holds everywhere POSIX process groups exist:

* the command leads a session of its own -- no controlling terminal, so
  nothing can be typed into the agent's or the user's terminal;
* stdin is ``/dev/null`` unless input is given;
* when the command exits, or its time is up, or the caller is interrupted,
  its whole process group is ended (SIGTERM, then SIGKILL).

A process that starts a session of its own escapes the group; the cage,
where it works, catches that too. Never used for the background shells,
which are meant to outlive the call that started them.
"""
from __future__ import annotations

import logging
import os
import signal
import subprocess
import threading
import time
from typing import Any, Callable, Optional, Sequence, Union

_GRACE_S = 2.0
_DRAIN_S = 2.0


def end_group(pgid: int, grace_s: float = _GRACE_S) -> None:
    """End every process in group ``pgid``: SIGTERM, then SIGKILL. Never raises."""
    if not pgid or pgid <= 1 or not hasattr(os, "killpg"):
        return
    if pgid == os.getpgrp():
        return              # never the caller's own group
    for sig in (signal.SIGTERM, signal.SIGKILL):
        try:
            os.killpg(pgid, sig)
        except (ProcessLookupError, PermissionError):
            return
        except OSError:
            return
        deadline = time.monotonic() + (grace_s if sig == signal.SIGTERM else 0.5)
        while time.monotonic() < deadline:
            try:
                os.killpg(pgid, 0)
            except (ProcessLookupError, PermissionError, OSError):
                return
            time.sleep(0.02)


def _reader(stream, sink: list) -> None:
    try:
        while True:
            chunk = stream.read(65536)
            if not chunk:
                break
            sink.append(chunk)
    except (OSError, ValueError):
        pass


#: Returned when the caller asked for the command to stop. 130 is what a
#: shell reports for an interrupted command, which is what this is.
STOPPED_RETURNCODE = 130

#: How often the wait looks at ``should_stop``. Short enough that Stop
#: feels immediate, long enough to cost nothing on a quiet wait.
_POLL_S = 0.25


def run(
    args: Union[str, Sequence[str]],
    *,
    cwd: Optional[str] = None,
    env: Optional[dict] = None,
    timeout: Optional[float] = None,
    shell: bool = False,
    input: Optional[str] = None,
    text: bool = True,
    should_stop: Optional[Callable[[], bool]] = None,
) -> subprocess.CompletedProcess:
    """Like ``subprocess.run(..., capture_output=True)``, contained.

    Returns when the command's own process exits; whatever it left running
    in its group is ended then. Raises ``subprocess.TimeoutExpired`` (with
    the output so far) when ``timeout`` passes, after ending the group.

    ``should_stop`` is asked every quarter second while waiting. When it
    says yes the group is ended and the result comes back with
    ``STOPPED_RETURNCODE`` and the output so far. Without it a command
    could only be waited out: a session that pressed Stop during a
    half-hour test run went on looking frozen until the run finished,
    because the worker was inside this wait and nothing in it ever
    looked (2026-09-17).
    """
    proc = subprocess.Popen(
        args, cwd=cwd, env=env, shell=shell,
        stdin=subprocess.PIPE if input is not None else subprocess.DEVNULL,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        start_new_session=True,
    )
    pgid = proc.pid
    out: list[bytes] = []
    err: list[bytes] = []
    readers = [threading.Thread(target=_reader, args=(proc.stdout, out), daemon=True),
               threading.Thread(target=_reader, args=(proc.stderr, err), daemon=True)]
    for t in readers:
        t.start()
    timed_out = False
    stopped = False
    try:
        if input is not None:
            try:
                proc.stdin.write(input.encode("utf-8") if isinstance(input, str) else input)
                proc.stdin.close()
            except (BrokenPipeError, OSError):
                pass
        try:
            probe = should_stop
            if probe is None:
                proc.wait(timeout=timeout)
            else:
                deadline = (None if timeout is None
                            else time.monotonic() + float(timeout))
                while True:
                    slice_s = _POLL_S
                    if deadline is not None:
                        left = deadline - time.monotonic()
                        if left <= 0:
                            raise subprocess.TimeoutExpired(args, timeout)
                        slice_s = min(_POLL_S, left)
                    try:
                        proc.wait(timeout=slice_s)
                        break
                    except subprocess.TimeoutExpired:
                        pass
                    if probe is None:
                        # The probe is gone (it raised); wait the rest out.
                        continue
                    try:
                        asked = bool(probe())
                    except Exception as exc:
                        # A probe that raises answered "never stop" at
                        # every poll, silently: the one failure a STOP
                        # mechanism must not have is the quiet kind.
                        # Found by a session testing the mechanism from
                        # outside, whose probe defined __bool__ and not
                        # __call__ (2026-09-17). It is not turned into a
                        # stop -- a flaky probe would then end every
                        # command -- but it is said once and not asked
                        # again, so the command runs unstoppable rather
                        # than unstoppable AND unremarked.
                        logging.getLogger(__name__).warning(
                            "stop probe raised %s: %s — this command can no "
                            "longer be stopped on request",
                            type(exc).__name__, exc)
                        probe = None
                        asked = False
                    if asked:
                        stopped = True
                        break
        except subprocess.TimeoutExpired:
            timed_out = True
    except BaseException:
        end_group(pgid)
        try:
            proc.kill()
            proc.wait(timeout=5)
        except Exception:
            pass
        raise
    finally:
        # Whatever the command left behind in its group ends with it.
        end_group(pgid, grace_s=_GRACE_S if timed_out else 0.3)
        if timed_out:
            try:
                proc.kill()
                proc.wait(timeout=5)
            except Exception:
                pass
    # A process that left the group may still hold the pipes open; take
    # what arrived and stop waiting for the rest.
    for t in readers:
        t.join(timeout=_DRAIN_S)
    for stream in (proc.stdout, proc.stderr):
        try:
            stream.close()
        except Exception:
            pass

    def _join(parts: list[bytes]) -> Any:
        data = b"".join(parts)
        return data.decode("utf-8", errors="replace") if text else data

    stdout, stderr = _join(out), _join(err)
    if timed_out:
        raise subprocess.TimeoutExpired(args, timeout, output=stdout, stderr=stderr)
    code = proc.returncode
    if stopped:
        # The command was ended on request, not by its own choice: say so
        # in the one place every caller already looks.
        code = STOPPED_RETURNCODE
    return subprocess.CompletedProcess(args, code, stdout, stderr)
