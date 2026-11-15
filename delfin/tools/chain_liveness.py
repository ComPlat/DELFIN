"""Is a chain's open step still being run? (JC3)

JC2's live run 4 (2026-11-15): after a crash the restarted agent found the
ORCA temporary files of the killed optimisation and told the user the step
was "actively running"; the chain state said "submitted" and nothing more.
Each submission records its runner (:func:`this_runner`, stored on
``ChainStep.runner``); :func:`runner_liveness` answers from it.

Kept apart from :mod:`delfin.tools.chain_state` on purpose: that module
must never launch a process (it holds no resume trigger, pinned by test),
and asking SLURM about a job that ran on another machine does.
"""

from __future__ import annotations

import os
import re
import socket
import subprocess
from typing import Any, Dict, Optional, Tuple

#: Liveness answers: the step's process still runs / has ended / cannot be
#: checked from this machine.
RUNNING, GONE, UNKNOWN = "running", "gone", "unknown"


def this_runner() -> Dict[str, Any]:
    """The identity of the current process, to record on a submitted step."""
    from delfin.agent import proc_identity

    pid = os.getpid()
    runner: Dict[str, Any] = {
        "pid": pid,
        "start": proc_identity.process_start(pid),
        "host": socket.gethostname(),
    }
    job = os.environ.get("SLURM_JOB_ID", "").strip()
    if job:
        runner["slurm_job"] = job
    return runner


def _slurm_job_state(job: str) -> Optional[str]:
    """The job's queue state, "" when it left the queue, None when unknown.

    One ``squeue`` call; only asked when the step ran on another machine.
    """
    if not re.fullmatch(r"[0-9_]+", job or ""):
        return None
    try:
        out = subprocess.run(["squeue", "-h", "-j", job, "-o", "%T"],
                             capture_output=True, text=True, timeout=15)
    except Exception:
        return None
    if out.returncode == 0:
        text = (out.stdout or "").strip()
        return text.splitlines()[0] if text else ""
    if "invalid job id" in (out.stderr or "").lower():
        return ""  # slurm forgets finished jobs after a while
    return None


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
    from delfin.agent import proc_identity

    pid, host = r.get("pid"), str(r.get("host") or "")
    alive = proc_identity.alive(pid, str(r.get("start") or ""), host)
    if alive is True:
        return RUNNING, f"still running (process {pid} on this machine)"
    if alive is False:
        return GONE, (f"not running — its process {pid} on "
                      f"{host or 'this machine'} has ended, the step was "
                      "interrupted")
    job = str(r.get("slurm_job") or "")
    if job:
        q = _slurm_job_state(job)
        if q:
            return RUNNING, f"still running (SLURM job {job} is {q})"
        if q == "":
            return GONE, (f"not running — SLURM job {job} has left the "
                          "queue, the step was interrupted")
    return UNKNOWN, (f"ran on {host or 'another machine'}; whether it still "
                     "runs cannot be checked from here")


__all__ = ["RUNNING", "GONE", "UNKNOWN", "this_runner", "runner_liveness"]
