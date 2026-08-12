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

Crash safety: every started job is additionally persisted to
``<workspace>/.delfin/bash_jobs.json`` (atomic temp-file + ``os.replace``
writes, same pattern as ``memory_store._atomic_write``). A dashboard or
kernel restart used to orphan every running calculation — the setsid
child survives, but the new process had no record of job ids, output
paths, or start times. Now an unknown job id re-attaches from that file:
a live pid (guarded against pid reuse via the process start time in
``/proc/<pid>/stat``) reports running with its output files reconnected;
a dead pid reports finished with ``exit_code=None`` (the real status is
unrecoverable once init reaped the orphan). ``drain_finished_events()``
turns the registry into a cheap, restart-safe completion-event source:
the watchdog marks the exit in the file at the moment it happens, and
each finished job is reported exactly once across drains and restarts.

The job's stdout / stderr tempfiles are opened with ``unbuffered=False``
in line-buffered mode so the agent can read partial progress while a
script is still running.

They are NOT removed when the job ends, is killed, or when the process
exits: ``bash_output`` is read after a job finishes, and reading a killed
job's output is how one finds out why it was killed. Their single
cleanup point is the registry record's ~7-day prune, which is the last
moment anything still knows where they are. (This paragraph used to claim
removal on kill and on shutdown; neither existed, and 6740 orphaned
``kit_bg_*`` files had accumulated in ``/tmp`` by the time it was
checked. They are 0600, so the content stays private to the owner.)
"""

from __future__ import annotations

import json
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

# How many background commands one agent may have running at once. There
# was no cap at all, in explicit contrast to background SUB-AGENTS, which
# have had an eight-slot semaphore since they were built. On a laptop an
# unbounded fan-out is untidy; on a shared cluster node it is other
# people's CPU, and the agent is the one participant that can start
# hundreds without noticing.
#
# Four rather than eight: a background command here is usually a real
# calculation, not a helper process.
# Ceiling on the derived cap. Background commands are a convenience, not
# the way to use a compute node -- that is what the batch scheduler is for
# -- so a large allocation must not turn into a large fan-out.
_MAX_BG_JOBS_CEILING = 8


def _affinity_cpus() -> int:
    """CPUs this process may actually run on.

    ``sched_getaffinity`` reflects cgroup cpusets and ``taskset``, which
    ``cpu_count`` does not: inside a container or a scheduler-managed
    cpuset, the machine's core count is not the process's.
    """
    try:
        return len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        return os.cpu_count() or 1


def _available_cpus() -> int:
    """How many CPUs this process was granted.

    A SLURM allocation states it outright and wins: on a shared node the
    affinity mask can still show the whole machine while the job was
    granted a slice of it, and spending the machine's width there is
    taking other people's CPU.
    """
    for var in ("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE"):
        try:
            n = int(str(os.environ.get(var, "")).strip())
            if n > 0:
                return n
        except (TypeError, ValueError):
            continue
    try:
        return max(1, int(_affinity_cpus()))
    except Exception:
        return 1


def _max_background_jobs() -> int:
    """The concurrency cap: an explicit setting, else what the grant allows.

    Read per call rather than cached: a user who raises it mid-session
    should not have to restart to get the slot, and an allocation can
    differ between one turn and the next.

    This used to be the constant 4 -- a sensible number for the
    workstation it was chosen on, and unrelated to the shared node the
    agent actually runs on, where it both oversubscribes a small slice
    and wastes a large one. Half the grant, floored at one and capped,
    still yields 4 on an eight-core box and adapts everywhere else.
    """
    try:
        from delfin.user_settings import load_settings
        raw = ((load_settings() or {}).get("agent") or {}).get(
            "max_background_jobs")
        if raw is not None:
            return max(1, int(raw))
    except Exception:
        pass
    try:
        return max(1, min(_MAX_BG_JOBS_CEILING, _available_cpus() // 2))
    except Exception:
        return 1
_OUTPUT_HEAD_DEFAULT = 60            # lines kept from head
_OUTPUT_TAIL_DEFAULT = 200           # lines kept from tail
_KILL_GRACE_S = 3.0                  # SIGTERM → SIGKILL gap

# --- Persistent registry (crash-safe re-attach + completion events) ------
_REGISTRY_DIRNAME = ".delfin"        # per-workspace state directory
_REGISTRY_FILENAME = "bash_jobs.json"
_REGISTRY_MAX_AGE_S = 7 * 24 * 3600  # prune records older than ~7 days
_EVENT_TAIL_CHARS = 500              # stdout/stderr tail carried per event
_RC_UNKNOWN = -257                   # poll() sentinel: finished, code unknown
# Per-user locator: job_id -> workspace, so a NEW process can find the right
# <workspace>/.delfin/bash_jobs.json from a bare job_id (bash_status /
# bash_output / bash_kill carry no workspace argument).
_INDEX_PATH = Path.home() / ".delfin" / "bash_jobs_index.json"
# Serializes read-modify-write cycles on the registry/index files within
# this process (watchdog threads vs. drain/status calls).
_FILE_LOCK = threading.Lock()


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


# ---------------------------------------------------------------------------
# Persistent registry file (crash-safe re-attach)
# ---------------------------------------------------------------------------


def _registry_path(workspace: str | Path) -> Path:
    return Path(workspace).expanduser() / _REGISTRY_DIRNAME / _REGISTRY_FILENAME


def _atomic_write_json(path: Path, data: dict) -> None:
    """Write JSON atomically (temp file + ``os.replace``).

    Same crash-safe pattern as ``memory_store._atomic_write``: a reader
    racing the writer (dashboard vs. watchdog thread vs. a second CLI
    session) must never observe an empty or torn registry file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent))
    tmp = Path(tmp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            json.dump(data, fh, indent=2, ensure_ascii=False)
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            try:
                tmp.unlink()
            except OSError:
                pass


def _load_registry_file(workspace: str | Path) -> dict:
    try:
        data = json.loads(_registry_path(workspace).read_text(encoding="utf-8"))
    except Exception:
        return {"jobs": {}}
    if not isinstance(data, dict) or not isinstance(data.get("jobs"), dict):
        return {"jobs": {}}
    return data


def _base_child_env() -> dict:
    """The environment a background command starts from.

    Foreground bash already ran on a scrubbed environment; this path took
    ``os.environ.copy()``, so backgrounding a command — a choice the MODEL
    makes, ``bash_background`` being model-callable — handed the child the
    provider key the agent itself is running on. Anything that prints its
    environment then wrote that key into the job's output file and returned
    it as a tool result.

    Falls back to the raw environment only if the scrubber cannot be
    imported at all, which would mean a broken install rather than a
    running agent.
    """
    try:
        from .api_client import _scrubbed_bash_env
        return _scrubbed_bash_env()
    except Exception:
        return os.environ.copy()


def _unlink_job_outputs(rec: dict) -> None:
    """Remove a finished job's stdout/stderr tempfiles.

    Nothing ever deleted these: 6740 orphaned ``kit_bg_*`` files were
    observed in ``/tmp``. They are 0600, so no third party reads them, but
    they hold whatever the command printed and they accumulate without
    limit. Best-effort — a file already gone is the expected case.
    """
    for key in ("stdout_path", "stderr_path"):
        raw = (rec or {}).get(key)
        if not raw:
            continue
        try:
            Path(str(raw)).unlink()
        except (OSError, ValueError):
            pass


def _prune_old_records(jobs: dict, now: float) -> bool:
    """Drop records started more than ~7 days ago. The 24 h hard timeout
    guarantees any such job is long over. Returns True when pruned.

    Dropping the record is also the last moment anything knows where the
    job's output files are, so they are unlinked here — otherwise the only
    map to them is gone and they stay in ``/tmp`` forever.
    """
    cutoff = now - _REGISTRY_MAX_AGE_S
    stale = [jid for jid, rec in jobs.items()
             if float((rec or {}).get("started_at") or 0.0) < cutoff]
    for jid in stale:
        _unlink_job_outputs(jobs.get(jid) or {})
        jobs.pop(jid, None)
    return bool(stale)


def _persist_job_start(workspace: str, record: dict) -> None:
    """Write the just-started job into the workspace registry. Best-effort —
    persistence must never prevent a job from starting."""
    try:
        with _FILE_LOCK:
            data = _load_registry_file(workspace)
            _prune_old_records(data["jobs"], time.time())
            data["jobs"][record["job_id"]] = record
            _atomic_write_json(_registry_path(workspace), data)
    except Exception:
        pass


def _update_job_record(workspace: str | Path, job_id: str, **fields) -> None:
    """Merge ``fields`` into one job's persisted record. Best-effort."""
    try:
        with _FILE_LOCK:
            data = _load_registry_file(workspace)
            rec = data["jobs"].get(job_id)
            if rec is None:
                return
            rec.update(fields)
            _atomic_write_json(_registry_path(workspace), data)
    except Exception:
        pass


def _proc_start_ticks(pid: int) -> Optional[int]:
    """Process start time in clock ticks (``/proc/<pid>/stat`` field 22).

    Recorded at job start and compared on re-attach as a pid-reuse guard:
    a recycled pid carries a different start time. Best-effort — returns
    None off Linux, in which case only the aliveness check applies."""
    try:
        stat = Path(f"/proc/{pid}/stat").read_text()
        # comm (field 2) may contain spaces/parens — split after the LAST ')'.
        return int(stat[stat.rindex(")") + 1:].split()[19])
    except Exception:
        return None


def _pid_alive(pid: int, start_ticks: Optional[int] = None) -> bool:
    if not pid or pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        pass         # exists but owned by another user — ticks check decides
    except Exception:
        return False
    if start_ticks:
        current = _proc_start_ticks(pid)
        if current is not None and current != start_ticks:
            return False   # pid was reused by an unrelated process
    return True


def _note_job_workspace(job_id: str, workspace: str) -> None:
    """Record job_id → workspace in the per-user locator index (pruned to
    the same ~7-day horizon as the registry). Best-effort."""
    try:
        with _FILE_LOCK:
            try:
                data = json.loads(_INDEX_PATH.read_text(encoding="utf-8"))
            except Exception:
                data = {}
            jobs = data.get("jobs") if isinstance(data, dict) else None
            if not isinstance(jobs, dict):
                jobs = {}
            now = time.time()
            jobs = {j: e for j, e in jobs.items()
                    if float((e or {}).get("ts") or 0.0) >= now - _REGISTRY_MAX_AGE_S}
            jobs[job_id] = {"workspace": workspace, "ts": now}
            _atomic_write_json(_INDEX_PATH, {"jobs": jobs})
    except Exception:
        pass


def _lookup_job_workspace(job_id: str) -> Optional[str]:
    try:
        data = json.loads(_INDEX_PATH.read_text(encoding="utf-8"))
        ws = ((data.get("jobs") or {}).get(job_id) or {}).get("workspace")
        return str(ws) if ws else None
    except Exception:
        return None


@dataclass
class _ReattachedProc:
    """Minimal Popen stand-in for a re-attached job — ``pid`` and
    ``returncode`` are all the registry's kill path touches."""

    pid: int
    returncode: Optional[int] = None

    def poll(self) -> Optional[int]:
        return self.returncode


class ReattachedJob:
    """A job recovered from the persistent registry after a restart.

    Duck-types the parts of :class:`BashJob` the tool layer uses (``poll``,
    ``elapsed_s``, ``status_dict``, output paths, ``proc.pid``). The child
    is no longer OUR child, so once it dies its exit status is gone (init
    reaped it) — a finished re-attached job reports ``exit_code=None``
    unless a previous session's watchdog recorded the real code."""

    def __init__(self, workspace: str, record: dict) -> None:
        self.workspace = str(workspace)
        self.job_id = str(record.get("job_id", ""))
        self.command = str(record.get("command", ""))
        self.description = str(record.get("description", ""))
        self.cwd = str(record.get("cwd", ""))
        self.timeout_s = int(record.get("timeout_s") or _DEFAULT_BG_TIMEOUT_S)
        self.stdout_path = Path(record.get("stdout_path") or "")
        self.stderr_path = Path(record.get("stderr_path") or "")
        self.started_at = float(record.get("started_at") or 0.0)   # epoch
        self.finished_at = record.get("finished_at")               # epoch|None
        self.exit_code = record.get("exit_code")
        self._start_ticks = record.get("proc_start_ticks")
        self.proc = _ReattachedProc(int(record.get("pid") or 0), self.exit_code)

    def poll(self) -> Optional[int]:
        if self.finished_at is None:
            if _pid_alive(self.proc.pid, self._start_ticks):
                return None
            # Died while no agent process was attached — record the moment we
            # noticed; the real exit status is unrecoverable.
            self.finished_at = time.time()
            _update_job_record(self.workspace, self.job_id,
                               finished_at=self.finished_at,
                               exit_code=self.exit_code)
        self.proc.returncode = self.exit_code
        return self.exit_code if self.exit_code is not None else _RC_UNKNOWN

    def elapsed_s(self) -> float:
        # Epoch-based (monotonic clocks do not survive a restart).
        end = self.finished_at if self.finished_at is not None else time.time()
        return round(max(0.0, end - self.started_at), 3)

    def status_dict(self) -> dict:
        rc = self.poll()
        status = {
            "job_id": self.job_id,
            "running": rc is None,
            "exit_code": self.exit_code,
            "elapsed_s": self.elapsed_s(),
            "command": self.command[:300],
            "description": self.description,
            "cwd": self.cwd,
            "stdout_path": str(self.stdout_path),
            "stderr_path": str(self.stderr_path),
            "reattached": True,
        }
        if rc is not None and self.exit_code is None:
            status["note"] = (
                "exit code unknown — the job outlived the process that "
                "started it (recovered from the persistent registry)")
        return status


def _reattach(job_id: str, workspace: str | Path | None = None) -> Optional[ReattachedJob]:
    """Rebuild a job view from ``<workspace>/.delfin/bash_jobs.json``.

    Called for job ids the in-memory registry does not know — i.e. after a
    dashboard/kernel restart orphaned a running calculation. The workspace
    comes from the caller when available, else from the locator index."""
    ws = str(workspace) if workspace else _lookup_job_workspace(job_id)
    if not ws:
        return None
    rec = _load_registry_file(ws).get("jobs", {}).get(job_id)
    if not rec:
        return None
    return ReattachedJob(ws, rec)


class _Registry:
    """Thread-safe job registry. Singleton via module-level instance."""

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._jobs: dict[str, BashJob] = {}

    def _new_job_id(self) -> str:
        # 8 hex chars — short enough for chat readability, large enough
        # to make accidental collisions vanishingly rare. An all-digit
        # token (~2 % of them) is indistinguishable from a SLURM job id
        # to any consumer that classifies by shape, so those are redrawn.
        for _ in range(10):
            jid = secrets.token_hex(4)
            if jid not in self._jobs and not jid.isdigit():
                return jid
        # Extremely unlikely; fall back to a longer token.
        jid = secrets.token_hex(8)
        return jid if not jid.isdigit() else "bg" + jid[2:]

    def start(
        self,
        command: str,
        cwd: str,
        description: str = "",
        timeout_s: int = _DEFAULT_BG_TIMEOUT_S,
        env: Optional[dict] = None,
        workspace: str | Path | None = None,
    ) -> BashJob:
        if not command.strip():
            raise ValueError("command must be non-empty")
        if timeout_s <= 0:
            raise ValueError("timeout_s must be positive")
        # Refuse rather than queue: a queued job looks started to the model,
        # which would then wait for output that cannot arrive yet. Refusing
        # with the count tells it what to do instead -- wait for one to
        # finish -- which is a decision it can act on.
        running = self.count_running()
        cap = _max_background_jobs()
        if running >= cap:
            raise ValueError(
                f"{running} background jobs are already running and the cap "
                f"is {cap}. Wait for one to finish (bash_status) before "
                "starting another, or raise agent.max_background_jobs if "
                "this machine can carry more.")
        timeout_s = min(timeout_s, _DEFAULT_BG_TIMEOUT_S)

        # tempfiles for stdout/stderr — opened append+text so the
        # subprocess can buffer to them and our reader can tail-read.
        stdout_fd, stdout_path = tempfile.mkstemp(prefix="kit_bg_", suffix=".stdout")
        stderr_fd, stderr_path = tempfile.mkstemp(prefix="kit_bg_", suffix=".stderr")
        os.close(stdout_fd)
        os.close(stderr_fd)
        sout = open(stdout_path, "w", buffering=1)   # line-buffered
        serr = open(stderr_path, "w", buffering=1)

        run_env = _base_child_env()
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

        # Persist the job BEFORE the watchdog starts, so its exit update can
        # never race the initial write. After a restart this record is the
        # only map back to the (setsid-surviving) child's pid, output files,
        # and start time.
        ws = str(Path(workspace).expanduser()) if workspace else cwd
        _persist_job_start(ws, {
            "job_id": jid,
            "pid": proc.pid,
            "proc_start_ticks": _proc_start_ticks(proc.pid),
            "command": command,
            "description": description,
            "cwd": cwd,
            "workspace": ws,
            "stdout_path": str(stdout_path),
            "stderr_path": str(stderr_path),
            "started_at": time.time(),      # epoch — survives restarts
            "timeout_s": timeout_s,
            "exit_code": None,
            "finished_at": None,
            "acknowledged": False,
        })
        _note_job_workspace(jid, ws)

        # Watchdog: enforces timeout, closes log file handles when done, and
        # marks the exit in the persistent registry the moment it happens —
        # the completion-event source for drain_finished_events().
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
                # Reap (instant on the normal path; after a timeout-kill the
                # child needs one more wait to yield its return code).
                try:
                    rc_final = proc.wait(timeout=_KILL_GRACE_S)
                except Exception:
                    rc_final = proc.poll()
                try:
                    sout.close()
                    serr.close()
                except Exception:
                    pass
                if job.finished_at is None:
                    job.finished_at = time.monotonic()
                _update_job_record(ws, jid,
                                   exit_code=rc_final,
                                   finished_at=time.time())

        t = threading.Thread(target=_watch, daemon=True, name=f"bashjob-{jid}")
        t.start()
        job._watchdog = t
        return job

    def get(
        self, job_id: str, workspace: str | Path | None = None,
    ) -> Optional["BashJob | ReattachedJob"]:
        with self._lock:
            job = self._jobs.get(job_id)
        if job is not None:
            return job
        # Unknown to THIS process — e.g. after a dashboard/kernel restart.
        # Re-attach from the persistent workspace registry.
        return _reattach(job_id, workspace)

    def count_running(self) -> int:
        """How many jobs this registry still has in flight.

        Counted from poll() rather than from a stored flag: a process that
        died without anyone asking would otherwise hold a slot forever, and
        a cap that leaks slots is worse than no cap -- it stops the agent
        working and gives no reason a user can act on.
        """
        try:
            return sum(1 for j in self.list_jobs(include_finished=False)
                       if j.poll() is None)
        except Exception:
            return 0

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


# ---------------------------------------------------------------------------
# Completion events (restart-safe, exactly-once)
# ---------------------------------------------------------------------------


def _tail_chars(path: str | Path, limit: int = _EVENT_TAIL_CHARS) -> str:
    """Last ``limit`` characters of a job output file, best-effort."""
    try:
        p = Path(path)
        size = p.stat().st_size
        with p.open("rb") as fh:
            if size > limit * 4:
                # Over-read to keep a full tail after multi-byte repair.
                fh.seek(-limit * 4, os.SEEK_END)
            text = fh.read().decode("utf-8", errors="replace")
    except Exception:
        return ""
    return text[-limit:]


def drain_finished_events(workspace: str | Path) -> list[dict]:
    """Return jobs that finished since the last drain — exactly once.

    The acknowledged flag lives in the registry file, so the exactly-once
    contract survives restarts: a new process only sees what no previous
    process drained. Jobs whose pid vanished while nobody was attached
    (agent restart during a multi-hour ORCA/xtb run) are folded in with
    ``exit_code=None``. This gives the engine an event-driven "job
    finished" signal to inject into the next turn — no blocking
    ``bash_status(wait_seconds=...)`` poll needed.

    Each event: ``job_id``, ``command`` (truncated), ``description``,
    ``exit_code`` (None when unrecoverable), ``runtime_s``, and ~500-char
    ``stdout_tail`` / ``stderr_tail``."""
    events: list[dict] = []
    try:
        with _FILE_LOCK:
            data = _load_registry_file(workspace)
            jobs = data.get("jobs", {})
            now = time.time()
            changed = _prune_old_records(jobs, now)
            for jid, rec in jobs.items():
                if rec.get("acknowledged"):
                    continue
                finished_at = rec.get("finished_at")
                if finished_at is None:
                    if _pid_alive(int(rec.get("pid") or 0),
                                  rec.get("proc_start_ticks")):
                        continue                    # still running
                    # Orphan died unattached — exit moment approximated now.
                    finished_at = rec["finished_at"] = now
                rec["acknowledged"] = True
                changed = True
                started = float(rec.get("started_at") or finished_at)
                events.append({
                    "job_id": jid,
                    "command": str(rec.get("command") or "")[:300],
                    "description": str(rec.get("description") or ""),
                    "exit_code": rec.get("exit_code"),
                    "runtime_s": round(max(0.0, finished_at - started), 3),
                    "stdout_tail": _tail_chars(rec.get("stdout_path") or ""),
                    "stderr_tail": _tail_chars(rec.get("stderr_path") or ""),
                })
            if changed:
                _atomic_write_json(_registry_path(workspace), data)
    except Exception:
        # Event delivery is best-effort; a broken registry file must never
        # take the tool layer down.
        return events
    return events


def running_jobs_for_workspace(workspace: str | Path) -> list[dict]:
    """Jobs still alive that were registered under ``workspace``.

    Read from the persisted registry rather than from this process's
    in-memory one: the question is who still OWNS the directory, and that
    includes jobs started by a thread whose own bookkeeping has already
    been torn down.

    Each entry: ``job_id``, ``pid``, ``command`` (truncated), ``cwd``.
    Never raises — an unreadable registry answers "nothing running",
    which leaves every caller behaving as it did before this existed."""
    out: list[dict] = []
    try:
        jobs = _load_registry_file(workspace).get("jobs", {})
    except Exception:
        return out
    for jid, rec in (jobs or {}).items():
        try:
            if (rec or {}).get("finished_at") is not None:
                continue
            pid = int((rec or {}).get("pid") or 0)
            if not _pid_alive(pid, (rec or {}).get("proc_start_ticks")):
                continue
            out.append({
                "job_id": jid,
                "pid": pid,
                "command": str((rec or {}).get("command") or "")[:200],
                "cwd": str((rec or {}).get("cwd") or ""),
            })
        except Exception:
            continue
    return out


def _known_workspaces() -> list[str]:
    """Every workspace the locator index has seen a job started in.

    Best-effort and already bounded: the index is pruned to the same
    ~7-day horizon as the registries themselves."""
    try:
        data = json.loads(_INDEX_PATH.read_text(encoding="utf-8"))
        jobs = data.get("jobs") if isinstance(data, dict) else None
        if not isinstance(jobs, dict):
            return []
        seen: list[str] = []
        for entry in jobs.values():
            ws = str((entry or {}).get("workspace") or "")
            if ws and ws not in seen:
                seen.append(ws)
        return seen
    except Exception:
        return []


def drain_all_finished_events(workspace: str | Path) -> list[dict]:
    """Completion events from ``workspace`` AND every other workspace a job
    was started in. Same exactly-once contract as ``drain_finished_events``.

    A job belongs to the workspace it was STARTED in, while the engine only
    knew the workspace the session is in NOW. A subagent working in an
    isolated worktree, a switched office folder, a resumed session with a
    different root -- in each case the job finished and nobody ever drained
    it, so the one event that tells the agent a long calculation is done
    never arrived.

    Ordering: the current workspace first, so its events lead the block.
    The acknowledged flag lives in each record, so a job already drained
    from its own workspace is not repeated here."""
    events = list(drain_finished_events(workspace) or [])
    try:
        primary = Path(workspace).expanduser().resolve()
    except (OSError, ValueError):
        primary = None
    for ws in _known_workspaces():
        try:
            candidate = Path(ws).expanduser().resolve()
        except (OSError, ValueError):
            continue
        if primary is not None and candidate == primary:
            continue
        if not candidate.is_dir():
            continue          # the folder is gone; its registry with it
        events.extend(drain_finished_events(candidate) or [])
    return events


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
