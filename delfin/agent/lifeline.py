"""Everything an agent starts ends with the terminal that started DELFIN.

Ctrl+C in the terminal running delfin-voila -- or that terminal closing, or
its SSH session dropping -- is the one stop for everything: the dashboard
server and its kernels (kept sessions included), the agent's background
shells, its MCP servers and the daemons. Nothing an agent started may keep
running on its own.

How it holds:

* The root -- delfin-voila, or a CLI command run from a terminal -- names
  itself the lifeline in the environment, by pid AND start time, so a reused
  pid never passes for it. Everything started below inherits that.
* A process started detached (a session or process group of its own, which
  no terminal signal reaches) is written to the lifeline's ledger. When the
  root stops, it ends every process on the ledger.
* A long-lived process -- a dashboard kernel, a daemon -- also watches the
  lifeline and ends itself once it is gone. That covers a root that was
  killed without a chance to clean up.
* The same watch ends the process on an emergency stop (``stop_all``),
  which reaches every machine sharing the home directory -- a lifeline
  only reaches its own.
"""

from __future__ import annotations

import json
import os
import signal
import threading
import time
from pathlib import Path
from typing import Callable, Optional, TypeVar

from delfin.agent import proc_identity

ENV_PID = "DELFIN_LIFELINE_PID"
ENV_TICKS = "DELFIN_LIFELINE_TICKS"
_DIR = Path.home() / ".delfin" / "lifeline"
_POLL_S = 3.0

#: Process start time in clock ticks; None where /proc is unavailable.
#: Read in one place for the whole codebase -- see ``proc_identity``.
_start_ticks = proc_identity.start_ticks


def _alive(pid: int, ticks: Optional[int]) -> bool:
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        pass
    except OSError:
        return False
    if ticks is None:
        return True
    return _start_ticks(pid) == ticks


def _env_for(pid: int, ticks: Optional[int]) -> dict[str, str]:
    return {ENV_PID: str(pid), ENV_TICKS: "" if ticks is None else str(ticks)}


def claim_root() -> dict[str, str]:
    """Make this process the lifeline of everything it starts. Returns the
    environment entries, which are also set in ``os.environ``."""
    pid = os.getpid()
    env = _env_for(pid, _start_ticks(pid))
    os.environ.update(env)
    return env


def current() -> Optional[tuple[int, Optional[int]]]:
    """``(pid, start ticks)`` of the lifeline this process lives under."""
    try:
        pid = int(os.environ.get(ENV_PID) or 0)
    except ValueError:
        return None
    if pid <= 0:
        return None
    raw = os.environ.get(ENV_TICKS, "")
    try:
        ticks = int(raw) if raw else None
    except ValueError:
        ticks = None
    return pid, ticks


def terminal_root() -> tuple[int, Optional[int]]:
    """The session leader of this process's terminal -- the shell. A daemon
    a CLI command starts from a terminal ends when that terminal does."""
    sid = os.getsid(0)
    return sid, _start_ticks(sid)


def child_env(base: Optional[dict] = None) -> dict[str, str]:
    """An environment for a detached child: the current lifeline, or the
    terminal's session leader when this process has none."""
    env = dict(os.environ if base is None else base)
    lifeline = current() or terminal_root()
    env.update(_env_for(*lifeline))
    return env


def _ledger(lifeline: tuple[int, Optional[int]]) -> Path:
    pid, ticks = lifeline
    return _DIR / f"{pid}-{ticks or 0}.jsonl"


def record_child(pid: int, kind: str = "") -> None:
    """Put a detached child on the lifeline's ledger. Never raises."""
    lifeline = current()
    if lifeline is None or pid <= 0:
        return
    try:
        from .state_paths import ensure_dir, open_append
        ensure_dir(_DIR)
        entry = {"pid": pid, "ticks": _start_ticks(pid), "kind": kind,
                 "at": time.time()}
        with open_append(_ledger(lifeline)) as fh:
            fh.write(json.dumps(entry) + "\n")
    except Exception:
        pass


def _signal_group(pid: int, sig: int) -> None:
    """A detached child leads its own process group: signal the group, so a
    shell's children go with it; fall back to the pid alone."""
    try:
        os.killpg(pid, sig)
    except Exception:
        try:
            os.kill(pid, sig)
        except Exception:
            pass


def end_children(lifeline: Optional[tuple[int, Optional[int]]] = None, *,
                 grace_s: float = 5.0) -> list[int]:
    """End every process on the ledger that is still the one recorded:
    SIGTERM, then SIGKILL after ``grace_s``. Returns the pids ended."""
    lifeline = lifeline or current()
    if lifeline is None:
        return []
    path = _ledger(lifeline)
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError:
        return []
    targets: list[tuple[int, Optional[int]]] = []
    for line in lines:
        try:
            entry = json.loads(line)
            pid = int(entry.get("pid") or 0)
        except (ValueError, TypeError):
            continue
        ticks = entry.get("ticks")
        if pid != os.getpid() and (pid, ticks) not in targets and _alive(pid, ticks):
            targets.append((pid, ticks))
    for pid, _ticks in targets:
        _signal_group(pid, signal.SIGTERM)
    deadline = time.monotonic() + grace_s
    while time.monotonic() < deadline and any(_alive(p, t) for p, t in targets):
        time.sleep(0.1)
    for pid, ticks in targets:
        if _alive(pid, ticks):
            _signal_group(pid, signal.SIGKILL)
    try:
        path.unlink()
    except OSError:
        pass
    return [pid for pid, _ticks in targets]


def _stopped() -> bool:
    """An emergency stop was given after this process started -- on this
    machine or any other that shares the home directory."""
    try:
        from .stop_all import stopped_since_start
        return stopped_since_start()
    except Exception:
        return False


def watch(on_gone: Callable[[], None], *, poll_s: float = _POLL_S
          ) -> Optional[threading.Thread]:
    """Call ``on_gone`` once, from a daemon thread, when the lifeline this
    process lives under is gone or an emergency stop was given
    (``stop_all``). A process under no lifeline still obeys the stop."""
    lifeline = current()

    def _run() -> None:
        while not _stopped():
            if lifeline is not None and not _alive(*lifeline):
                break
            time.sleep(poll_s)
        try:
            on_gone()
        except Exception:
            pass

    thread = threading.Thread(target=_run, name="delfin-lifeline", daemon=True)
    thread.start()
    return thread


def exit_now() -> None:
    """The end for a process whose lifeline is gone."""
    os._exit(0)


def guard_daemon() -> Optional[threading.Thread]:
    """For a daemon's entry point: end the process with its lifeline.

    Started by delfin-voila or a CLI command, the lifeline is in the
    environment. Started by hand from a shell (``nohup python -m ...``), it
    is that terminal's session leader -- unless the daemon leads a session
    of its own, which is no terminal to follow.
    """
    try:
        from . import process_guard as _process_guard
        _process_guard.protect("daemon")
    except Exception:
        pass
    if current() is None:
        sid = os.getsid(0)
        if sid != os.getpid():
            os.environ.update(_env_for(sid, _start_ticks(sid)))
    return watch(exit_now)


_T = TypeVar("_T")
_FORKER_LOCK = threading.Lock()
_forker: Optional[tuple[int, threading.Thread, "queue.Queue"]] = None


def _forker_loop(jobs: "queue.Queue") -> None:
    while True:
        start, box, done = jobs.get()
        try:
            box.append((True, start()))
        except BaseException as exc:            # noqa: BLE001
            box.append((False, exc))
        finally:
            done.set()


def start_bound_to_process(start: Callable[[], _T]) -> _T:
    """Run *start* -- one that forks -- on a thread that lives as long as the
    process, and return what it returned.

    PR_SET_PDEATHSIG, which bubblewrap's ``--die-with-parent`` and
    ``parent_death_signal`` set, fires when the THREAD that forked ends,
    not the process. A caged child started from a tool call's thread was
    SIGKILLed as soon as the call returned: every background job, and an
    MCP server started on first use. Forked from here, the child ends when
    the process does, which is what the flag is there for.
    """
    import queue

    global _forker
    with _FORKER_LOCK:
        if _forker is None or _forker[0] != os.getpid() \
                or not _forker[1].is_alive():
            jobs: queue.Queue = queue.Queue()
            thread = threading.Thread(target=_forker_loop, args=(jobs,),
                                      name="delfin-forker", daemon=True)
            thread.start()
            _forker = (os.getpid(), thread, jobs)
        _pid, thread, jobs = _forker
    if threading.current_thread() is thread:
        return start()
    box: list = []
    done = threading.Event()
    jobs.put((start, box, done))
    done.wait()
    ok, value = box[0]
    if ok:
        return value
    raise value


def parent_death_signal() -> None:
    """For ``preexec_fn``: the child receives SIGTERM when its parent dies,
    however the parent died. Linux only; a no-op elsewhere."""
    try:
        import ctypes
        # The running process's own C library: "libc.so.6" is glibc's name
        # and does not exist on musl (Alpine), where this was a silent no-op.
        libc = ctypes.CDLL(None, use_errno=True)
        libc.prctl(1, signal.SIGTERM)      # PR_SET_PDEATHSIG
    except Exception:
        pass
