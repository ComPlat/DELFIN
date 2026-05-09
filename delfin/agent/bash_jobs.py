"""Background bash-job registry for the KIT-Toolbox coding agent.

Long-running tasks (Bayesian-opt runs, training, large pytest sessions)
need to survive past a single tool-call turn. This module provides a
thread-safe registry of subprocess-Popen jobs whose stdout / stderr
land in tempfiles (no in-memory buffering, no reader threads needed).

Tools that exercise this module:

* ``bash_background(command, cwd, description, timeout_s)`` — start
  a job, return ``job_id`` immediately.
* ``bash_status(job_id)`` — running? exit code? elapsed?
* ``bash_output(job_id, head_lines, tail_lines)`` — read what's been
  written so far. By default keeps the head + tail with a marker in
  between, mirroring the smart-truncate behavior of the synchronous
  bash tool.
* ``bash_kill(job_id)`` — SIGTERM, then SIGKILL after grace period.

Job IDs are short hex tokens (8 chars) so the agent can paste them
into chat without scrolling. Jobs are kept in the registry until
``bash_kill`` or until the singleton is GC'd at process exit. The
caller — the bash gate in api_client.py — is responsible for
running the same auto-allow / deny / secret-scan checks BEFORE any
job is started; this module trusts its inputs.

The job's stdout / stderr tempfiles are opened with ``unbuffered=False``
in line-buffered mode so the agent can read partial progress while a
script is still running. Tempfiles are removed when the job is
explicitly killed or when the registry shuts down. They survive a
crash so post-mortem inspection is possible.
"""

from __future__ import annotations

import os
import secrets
import signal
import subprocess
import tempfile
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


_DEFAULT_BG_TIMEOUT_S = 24 * 3600    # 24 h hard cap
_OUTPUT_HEAD_DEFAULT = 60            # lines kept from head
_OUTPUT_TAIL_DEFAULT = 200           # lines kept from tail
_KILL_GRACE_S = 3.0                  # SIGTERM → SIGKILL gap


@dataclass
class BashJob:
    """One long-running shell command tracked by the registry."""

    job_id: str
    command: str
    description: str
    cwd: str
    started_at: float
    proc: subprocess.Popen
    stdout_path: Path
    stderr_path: Path
    timeout_s: int
    finished_at: Optional[float] = None
    _watchdog: Optional[threading.Thread] = field(default=None, repr=False)

    def poll(self) -> Optional[int]:
        rc = self.proc.poll()
        if rc is not None and self.finished_at is None:
            self.finished_at = time.monotonic()
        return rc

    def elapsed_s(self) -> float:
        end = self.finished_at if self.finished_at is not None else time.monotonic()
        return round(end - self.started_at, 3)

    def status_dict(self) -> dict:
        rc = self.poll()
        return {
            "job_id": self.job_id,
            "running": rc is None,
            "exit_code": rc,
            "elapsed_s": self.elapsed_s(),
            "command": self.command[:300],
            "description": self.description,
            "cwd": self.cwd,
            "stdout_path": str(self.stdout_path),
            "stderr_path": str(self.stderr_path),
        }


class _Registry:
    """Thread-safe job registry. Singleton via module-level instance."""

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._jobs: dict[str, BashJob] = {}

    def _new_job_id(self) -> str:
        # 8 hex chars — short enough for chat readability, large enough
        # to make accidental collisions vanishingly rare.
        for _ in range(10):
            jid = secrets.token_hex(4)
            if jid not in self._jobs:
                return jid
        # Extremely unlikely; fall back to a longer token.
        return secrets.token_hex(8)

    def start(
        self,
        command: str,
        cwd: str,
        description: str = "",
        timeout_s: int = _DEFAULT_BG_TIMEOUT_S,
        env: Optional[dict] = None,
    ) -> BashJob:
        if not command.strip():
            raise ValueError("command must be non-empty")
        if timeout_s <= 0:
            raise ValueError("timeout_s must be positive")
        timeout_s = min(timeout_s, _DEFAULT_BG_TIMEOUT_S)

        # tempfiles for stdout/stderr — opened append+text so the
        # subprocess can buffer to them and our reader can tail-read.
        stdout_fd, stdout_path = tempfile.mkstemp(prefix="kit_bg_", suffix=".stdout")
        stderr_fd, stderr_path = tempfile.mkstemp(prefix="kit_bg_", suffix=".stderr")
        os.close(stdout_fd)
        os.close(stderr_fd)
        sout = open(stdout_path, "w", buffering=1)   # line-buffered
        serr = open(stderr_path, "w", buffering=1)

        run_env = os.environ.copy()
        run_env.setdefault("LC_ALL", "C.UTF-8")
        run_env.setdefault("LANG", "C.UTF-8")
        # Buffering off in the child Python so output appears live.
        run_env.setdefault("PYTHONUNBUFFERED", "1")
        if env:
            run_env.update(env)

        proc = subprocess.Popen(
            ["/bin/bash", "-c", command],
            cwd=cwd,
            env=run_env,
            stdout=sout,
            stderr=serr,
            stdin=subprocess.DEVNULL,
            # New process group so we can SIGTERM the whole tree.
            preexec_fn=os.setsid,
            text=True,
        )

        with self._lock:
            jid = self._new_job_id()
            job = BashJob(
                job_id=jid,
                command=command,
                description=description,
                cwd=cwd,
                started_at=time.monotonic(),
                proc=proc,
                stdout_path=Path(stdout_path),
                stderr_path=Path(stderr_path),
                timeout_s=timeout_s,
            )
            self._jobs[jid] = job

        # Watchdog: enforces timeout + closes log file handles when done.
        def _watch():
            try:
                rc = proc.wait(timeout=timeout_s)
                _ = rc
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
                    time.sleep(_KILL_GRACE_S)
                    if proc.poll() is None:
                        os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
                except Exception:
                    pass
            finally:
                try:
                    sout.close()
                    serr.close()
                except Exception:
                    pass
                if job.finished_at is None:
                    job.finished_at = time.monotonic()

        t = threading.Thread(target=_watch, daemon=True, name=f"bashjob-{jid}")
        t.start()
        job._watchdog = t
        return job

    def get(self, job_id: str) -> Optional[BashJob]:
        with self._lock:
            return self._jobs.get(job_id)

    def list_jobs(self, include_finished: bool = True) -> list[BashJob]:
        with self._lock:
            jobs = list(self._jobs.values())
        if not include_finished:
            jobs = [j for j in jobs if j.poll() is None]
        return jobs

    def kill(self, job_id: str, *, sig: int = signal.SIGTERM) -> tuple[bool, str]:
        job = self.get(job_id)
        if job is None:
            return False, f"unknown job_id: {job_id}"
        if job.poll() is not None:
            return True, f"job already finished (rc={job.proc.returncode})"
        try:
            os.killpg(os.getpgid(job.proc.pid), sig)
        except ProcessLookupError:
            return True, "process already exited"
        except Exception as exc:
            return False, f"kill failed: {exc}"
        # Give SIGTERM a moment, then escalate to SIGKILL.
        if sig == signal.SIGTERM:
            for _ in range(int(_KILL_GRACE_S * 10)):
                if job.poll() is not None:
                    return True, f"terminated (rc={job.proc.returncode})"
                time.sleep(0.1)
            try:
                os.killpg(os.getpgid(job.proc.pid), signal.SIGKILL)
            except Exception:
                pass
            return True, "SIGTERM ignored, sent SIGKILL"
        return True, f"signal {sig} delivered"


_REGISTRY = _Registry()


def get_registry() -> _Registry:
    """Return the process-wide job registry singleton."""
    return _REGISTRY


def read_output(
    job: BashJob,
    head_lines: int = _OUTPUT_HEAD_DEFAULT,
    tail_lines: int = _OUTPUT_TAIL_DEFAULT,
) -> dict:
    """Read the job's stdout/stderr with smart head+tail truncation.

    Tracebacks live at the END so the tail is always preserved. Empty
    requests (head=0 OR tail=0) suppress that side. Returns a dict with
    ``stdout``, ``stderr``, total line counts, and a ``truncated`` flag.
    """
    out_lines = _read_lines(job.stdout_path)
    err_lines = _read_lines(job.stderr_path)
    return {
        "stdout": _format_segment(out_lines, head_lines, tail_lines),
        "stderr": _format_segment(err_lines, head_lines, tail_lines),
        "stdout_total_lines": len(out_lines),
        "stderr_total_lines": len(err_lines),
        "stdout_truncated": (head_lines + tail_lines) < len(out_lines),
        "stderr_truncated": (head_lines + tail_lines) < len(err_lines),
    }


def _read_lines(path: Path) -> list[str]:
    try:
        with path.open("r", encoding="utf-8", errors="replace") as fh:
            return fh.read().splitlines()
    except FileNotFoundError:
        return []
    except Exception:
        return []


def _format_segment(lines: list[str], head: int, tail: int) -> str:
    n = len(lines)
    head = max(0, head)
    tail = max(0, tail)
    if n == 0:
        return ""
    if head + tail >= n:
        return "\n".join(lines)
    head_part = lines[:head] if head else []
    tail_part = lines[-tail:] if tail else []
    omitted = n - len(head_part) - len(tail_part)
    middle = [
        f"... ({omitted} lines from the middle omitted; "
        f"head and tail preserved so tracebacks survive)"
    ]
    return "\n".join(head_part + middle + tail_part)
