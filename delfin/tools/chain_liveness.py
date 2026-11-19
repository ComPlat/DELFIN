"""Is a chain's open step still being run? (JC3)

JC2's live run 4 (2026-11-15): after a crash the restarted agent found the
ORCA temporary files of the killed optimisation and told the user the step
was "actively running"; the chain state said "submitted" and nothing more.
Each submission records its runner (:func:`this_runner`, stored on
``ChainStep.runner``); :func:`runner_liveness` answers from it.

Two signs are asked, the strongest first:

* the runner's process, by pid AND start fingerprint -- only where that pid
  means the same process: on the same machine and in the same PID
  namespace.  JC3's live run 5 (2026-11-15): the agent started the chain
  inside its sandbox, which has a PID namespace of its own, so the runner
  was "pid 3"; asked from outside, pid 3 is a kernel thread, and the
  killed step read "still running".  A pid only ever answers RUNNING when
  its start fingerprint matches the recorded one.
* a heartbeat file the running step touches every :data:`BEAT_SECONDS`
  (:func:`beating`), which any process on any node can read.  No touch for
  :data:`STALE_SECONDS` means the step is gone.

A SLURM job is no sign: the job an agent runs in outlives the chains it
starts (live run 5's agent read its own job as the chain's).

Kept apart from :mod:`delfin.tools.chain_state` on purpose: that module
must never launch anything (it holds no resume trigger, pinned by test),
and the heartbeat is a thread.
"""

from __future__ import annotations

import contextlib
import os
import signal
import socket
import threading
import time
from pathlib import Path
from typing import Any, Callable, Dict, Iterator, Optional, Tuple

#: Liveness answers: the step's process still runs / has ended / cannot be
#: checked from here.
RUNNING, GONE, UNKNOWN = "running", "gone", "unknown"

#: The running step touches its heartbeat this often ...
BEAT_SECONDS = 30
#: ... and is taken as gone once the last touch is older than this.  Four
#: beats: a slow shared filesystem may delay one or two.
STALE_SECONDS = 120

#: The heartbeat file in a chain's work dir.
HEARTBEAT_NAME = ".chain_heartbeat"


def _pid_namespace() -> str:
    """This process's PID namespace ("pid:[4026531836]"), "" if unknown."""
    try:
        return os.readlink("/proc/self/ns/pid")
    except OSError:
        return ""


def this_runner(chain_dir: Optional[Path] = None) -> Dict[str, Any]:
    """The identity of the current process, to record on a submitted step.

    With *chain_dir* the record names the heartbeat file there, which
    :func:`beating` keeps fresh while the step runs.
    """
    from delfin.agent import proc_identity

    pid = os.getpid()
    runner: Dict[str, Any] = {
        "pid": pid,
        "start": proc_identity.process_start(pid),
        "host": socket.gethostname(),
        "pid_ns": _pid_namespace(),
    }
    if chain_dir is not None:
        runner["heartbeat"] = str(Path(chain_dir) / HEARTBEAT_NAME)
    return runner


@contextlib.contextmanager
def beating(path: Optional[str]) -> Iterator[None]:
    """Touch *path* now and every :data:`BEAT_SECONDS` until the block ends.

    A failing touch is skipped, never raised: the heartbeat must not stop
    the calculation it reports on.  Without *path* nothing happens.
    """
    if not path:
        yield
        return
    target = Path(path)

    def _touch() -> None:
        try:
            target.touch()
        except OSError:
            pass

    stop = threading.Event()

    def _beat() -> None:
        while not stop.wait(BEAT_SECONDS):
            _touch()

    _touch()
    thread = threading.Thread(target=_beat, name="chain-heartbeat",
                              daemon=True)
    thread.start()
    try:
        yield
    finally:
        stop.set()
        thread.join(timeout=5)


#: The signals a run is told to end with: ``scancel`` and a batch system's
#: time limit send SIGTERM, a closing terminal or an ended shell call SIGHUP.
END_SIGNALS = ("SIGTERM", "SIGHUP")


@contextlib.contextmanager
def ended_on_signal(record: Callable[[str], None]) -> Iterator[None]:
    """Call ``record(<signal name>)`` when the run is ended by a signal, then
    end exactly as before.

    Without it, a step ended by ``scancel`` reads "still running" until
    its heartbeat is :data:`STALE_SECONDS` old, and a resume in that window
    is refused (JC3 live run 6). The previous handler is restored first
    and then given the signal. Only where that ends the process (the
    default action) is the end recorded: a Python handler is called and
    decides for itself, an ignored signal stays ignored. Only the main thread may set handlers; elsewhere (a tool
    server's worker thread) this does nothing and the heartbeat decides.
    """
    if threading.current_thread() is not threading.main_thread():
        yield
        return
    previous: Dict[int, Any] = {}

    def _handler(signum, frame):
        prev = previous.get(signum, signal.SIG_DFL)
        signal.signal(signum, prev)
        if callable(prev):
            # Whoever set it decides whether the process ends; if it raises,
            # the pipeline's own except-clause records the end.
            prev(signum, frame)
            return
        if prev == signal.SIG_IGN:
            return
        # The default action: the process ends now. Only then is "ended"
        # true -- recorded for a signal that does not end the run, it would
        # let a resume start the same step twice.
        try:
            record(signal.Signals(signum).name)
        except Exception:  # noqa: BLE001 - the signal must still end the run
            pass
        os.kill(os.getpid(), signum)

    for name in END_SIGNALS:
        sig = getattr(signal, name, None)
        if sig is None:
            continue
        try:
            previous[sig] = signal.getsignal(sig)
            signal.signal(sig, _handler)
        except (ValueError, OSError):
            previous.pop(sig, None)
    try:
        yield
    finally:
        for sig, prev in previous.items():
            try:
                if signal.getsignal(sig) is _handler:
                    signal.signal(sig, prev)
            except (ValueError, OSError):
                pass


def _process_answer(r: Dict[str, Any]) -> Optional[Tuple[str, str]]:
    """RUNNING/GONE from the recorded process, None when it cannot tell."""
    host = str(r.get("host") or "")
    ns = str(r.get("pid_ns") or "")
    if host != socket.gethostname() or not ns or ns != _pid_namespace():
        return None  # that pid names a different process here, or none
    from delfin.agent import proc_identity

    pid, start = r.get("pid"), str(r.get("start") or "")
    alive = proc_identity.alive(pid, start, host)
    if alive is False:
        return GONE, (f"not running — its process {pid} has ended, the "
                      "step was interrupted")
    now = proc_identity.process_start(pid)
    if start and now == start:
        return RUNNING, f"still running (process {pid} on this machine)"
    if start and now and now != start:
        return GONE, (f"not running — its process {pid} has ended, the "
                      "step was interrupted")
    return None  # pid in use, but whose cannot be read (another user's)


def _heartbeat_answer(r: Dict[str, Any]) -> Optional[Tuple[str, str]]:
    path = r.get("heartbeat")
    if not path:
        return None
    try:
        age = time.time() - Path(path).stat().st_mtime
    except OSError:
        return None
    if age <= STALE_SECONDS:
        return RUNNING, f"still running (last sign of life {max(int(age), 0)} s ago)"
    when = time.strftime("%Y-%m-%d %H:%M", time.localtime(time.time() - age))
    return GONE, (f"not running — no sign of life since {when}, the step "
                  "was interrupted")


def runner_liveness(step: Any) -> Tuple[str, str]:
    """``(RUNNING | GONE | UNKNOWN, one-line reason)`` for an open step.

    Only ``submitted``/``running`` steps have a runner to ask; any other
    state answers ``("", "")``.
    """
    if getattr(step, "state", "") not in ("submitted", "running"):
        return "", ""
    r = getattr(step, "runner", None) or {}
    if not r.get("pid"):
        return UNKNOWN, ("no process recorded (state written by an older "
                         "DELFIN) — check before resuming")
    if r.get("ended"):
        # The run raised and its process went on (an agent service, a
        # tool server, a notebook kernel): the pid is alive, the step is not.
        return GONE, (f"not running — its run ended with {r['ended']} "
                      "before a result was recorded")
    answer = _process_answer(r) or _heartbeat_answer(r)
    if answer is not None:
        return answer
    host = str(r.get("host") or "")
    return UNKNOWN, (f"ran on {host or 'another machine'}; whether it still "
                     "runs cannot be checked from here")


__all__ = ["RUNNING", "GONE", "UNKNOWN", "BEAT_SECONDS", "STALE_SECONDS",
           "HEARTBEAT_NAME", "END_SIGNALS", "this_runner", "beating",
           "ended_on_signal", "runner_liveness"]
