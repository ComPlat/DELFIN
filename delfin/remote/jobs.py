"""DELFIN workflows and single ORCA runs as Slurm jobs on one of DELFIN's clusters.

Input: job folders in a session's workspace -- ``input.txt`` and
``CONTROL.txt`` for a DELFIN workflow, an ORCA ``.inp`` for one run -- with
the resources asked for, a Connection (``delfin.remote.connection``) and the
commit the PC runs (``delfin.remote.install``). Output: a ``Batch``, which is
what an approval dialog shows and what refuses the batch before anyone is
asked (``problems``); once a person approved it, the submitted jobs.

Semantics:

* A folder lies inside the workspace, holds at most ``MAX_FILES`` files and
  ``MAX_BYTES`` bytes and no link (rsync -r copies none), and every file
  passes the caller's ``refuse`` check (secrets, DELFIN's masked state). A
  DELFIN workflow needs input.txt and CONTROL.txt; an ORCA run names its
  .inp, which must be in the folder.
* Resources follow the dashboard's conventions (dashboard/backend_slurm):
  one task, cpus-per-task = the cores (CONTROL's PAL, the input's nprocs),
  memory = cores x maxcore MB; the time limit given, else CONTROL's
  job_timeout_hours, else ``DEFAULT_TIME``; the partition given, else the
  cluster's default_partition, and one its settings allow. The job runs
  submit_delfin.sh of the PC's commit's install, the venv used in place
  (DELFIN_STAGE_VENV=0; staging would tar a 1 GB venv) and TMPDIR the
  node's /tmp where the node names none (the script falls back to
  /scratch, which not every cluster has).
* The budget: the batch's hours as the cluster's budget counts them --
  core-hours when its unit names cores, else job-hours, one time limit per
  job -- against what the session has left. A batch that does not fit, or
  that exceeds the managed ceiling on core-hours, is refused before anyone
  is asked.
* The install: when the PC's commit is not installed on the cluster, the
  batch says so (the dialog shows it) -- only a published commit, and not
  when a managed policy turned automatic installs off. An approved batch
  then waits for it without holding anyone: ``queue_after_install`` starts
  the install (detached, on the login node's share) and keeps the batch,
  exactly as approved, in ``~/.delfin/remote_pending.json``; the job watcher
  (delfin.remote.watch) submits it once the install is complete, after the
  same checks of the files and the budget as ``submit`` makes, and asks
  nobody again. ``submit`` itself still installs first and waits, for a
  caller that wants that.
* Submission: each folder goes to ``<work_dir>/<stamp>-<name>/`` and sbatch
  runs there through ``bash -c``, which ctun's own budget guard reads, with
  ``--export=ALL`` and DELFIN's variables in sbatch's environment (as the
  dashboard does: ``--export=ALL,VAR=...`` made bwUniCluster hold jobs) and
  ``--comment=delfin-commit=<sha>``, which keeps the install while the job
  is queued or runs. Every submitted job is recorded in
  ``~/.delfin/remote_jobs.jsonl`` for the job watcher.
* ``cancel`` cancels only active jobs of the user's on that cluster.
"""

from __future__ import annotations

import contextlib
import json
import os
import re
import secrets
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterable, Optional

from . import install as _install
from .errors import RemoteError

MAX_JOBS = 50
MAX_FILES = 2000
MAX_BYTES = 200 * 1024 * 1024
DEFAULT_TIME = "12:00:00"
#: maxcore (MB per core) when neither CONTROL.txt nor the input names one.
DEFAULT_MAXCORE = 2000
#: How long ``submit`` waits for a first install before it gives up on the
#: batch (the install itself carries on): a new commit's venv took 23 min
#: on kcist, where /home is a network filesystem.
INSTALL_WAIT_S = 45 * 60

_SLUG = re.compile(r"[^A-Za-z0-9._-]+")
_PARTITION = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*")
_JOB_ID = re.compile(r"[0-9]+(?:_[0-9]+)?")

#: Runs on the login node in the job's folder: sbatch with DELFIN's
#: variables in its environment, the commit's own job script. Arguments:
#: the folder, the install, the mode, the ORCA input, the job's name, then
#: sbatch's options. "~/" arrives unexpanded (argv words) and is expanded here.
_SUBMIT = r'''set -e
home_path() { case "$1" in "~") printf '%s' "$HOME" ;; "~/"*) printf '%s/%s' "$HOME" "${1#"~/"}" ;; *) printf '%s' "$1" ;; esac; }
dir="$(home_path "$1")"; install="$(home_path "$2")"; mode="$3"; inp="$4"; name="$5"; shift 5
cd "$dir"
if [ -f "$HOME/.delfin_env.sh" ]; then . "$HOME/.delfin_env.sh"; fi
export DELFIN_MODE="$mode" DELFIN_JOB_NAME="$name" DELFIN_INP_FILE="$inp"
export DELFIN_REPO="$install" DELFIN_VENV="$install/.venv"
export DELFIN_OMPI_HOME="${OPENMPI_PREFIX:-}" DELFIN_ORCA_BASE="${ORCA_BASE:-}"
export DELFIN_STAGE_VENV=0 DELFIN_TMP_FALLBACK=/tmp
unset SBATCH_GET_USER_ENV SBATCH_EXPORT TMPDIR
sbatch --parsable --export=ALL "$@" "$install/delfin/submit_templates/submit_delfin.sh"
'''


class JobRefused(RemoteError):
    """A cancel the cluster's job list does not allow."""

    code = "job_refused"


@dataclass(frozen=True)
class Job:
    """One job of a batch, as it will be submitted."""

    folder: Path      # in the workspace
    kind: str         # "delfin" (a workflow) | "orca" (one input)
    input: str        # the ORCA input's file name; "" for a workflow
    name: str         # sbatch --job-name
    cores: int
    mem_mb: int
    time: str         # in sbatch's form
    hours: float      # the time limit in hours
    files: int
    size: int         # bytes
    #: The files approved, relative to the folder: what is uploaded, and
    #: with their sizes and times what must still be there at the upload.
    paths: tuple[str, ...] = ()
    fingerprint: tuple[tuple[str, int, int], ...] = ()


class PartialSubmit(RemoteError):
    """A batch that failed after some of its jobs were submitted."""

    code = "partly_submitted"

    def __init__(self, message: str, done: list[dict]):
        super().__init__(message)
        self.done = done


#: A folder the job watcher brought an ended job's results into
#: (delfin.remote.watch): never part of a job folder that is uploaded.
RESULTS_MARKER = ".delfin-results"


@dataclass
class Batch:
    """A batch as checked before anyone is asked: what the dialog shows."""

    cluster: str
    target: str
    partition: str
    jobs: list[Job] = field(default_factory=list)
    commit: str = ""
    install: str = ""            # "complete" | "needed" | "running" | ""
    unit: str = ""
    requested: float = 0.0       # in the budget's unit
    core_hours: float = 0.0
    used: Optional[float] = None
    limit: Optional[float] = None
    problems: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.problems

    def preview(self) -> str:
        """What the approval dialog says, one fact a line."""
        n = len(self.jobs)
        lines = [f"submit {n} job{'s' if n != 1 else ''} to {self.cluster} ({self.target})",
                 f"partition {self.partition}"]
        for job in self.jobs:
            what = f"ORCA {job.input}" if job.kind == "orca" else "DELFIN workflow"
            lines.append(f"  {job.name}: {what} in {job.folder.name}/ -- {job.cores} cores, "
                         f"{_memory(job.mem_mb)}, {job.time}")
        lines.append(f"at most {self.core_hours:g} core-hours"
                     + (f" ({self.requested:g} {self.unit} by the budget's count)"
                        if self.unit and "core" not in self.unit else ""))
        if self.limit is not None and self.used is not None:
            lines.append(f"budget: {self.used:g} of {self.limit:g} {self.unit} used, "
                         f"{max(0.0, self.limit - self.used):g} left")
        else:
            lines.append("budget: no session limit set")
        if self.install in ("needed", "running"):
            lines.append(f"first installs DELFIN {self.commit[:12]} on {self.cluster} (once per "
                         "commit, on the login node's share: up to an hour); the batch goes "
                         "in by itself once it is installed")
        lines += self.notes
        return "\n".join(lines)


def plan(conn, workspace: Path, requests: Iterable[dict], *,
         refuse: Optional[Callable[[Path], Optional[str]]] = None,
         max_core_hours: Optional[float] = None, auto_install: bool = True,
         commit: Optional["_install.Commit"] = None) -> Batch:
    """The batch ``requests`` make on ``conn``'s cluster, checked."""
    cluster = conn.cluster
    batch = Batch(cluster=cluster.name, target=cluster.target, partition="")
    requests = list(requests)
    if not requests:
        batch.problems.append("no jobs: give at least one folder")
        return batch
    if len(requests) > MAX_JOBS:
        batch.problems.append(f"{len(requests)} jobs in one batch; at most {MAX_JOBS}")
        return batch
    workspace = Path(workspace).resolve()
    partitions = set()
    for number, request in enumerate(requests, 1):
        job, problem = _job(workspace, request, refuse, number)
        if problem:
            batch.problems.append(problem)
            continue
        partition, problem = _partition(cluster, request.get("partition"))
        if problem:
            batch.problems.append(f"job {number}: {problem}")
            continue
        partitions.add(partition)
        batch.jobs.append(job)
    if len(partitions) > 1:
        batch.problems.append(f"one partition per batch, not {sorted(partitions)}")
    batch.partition = next(iter(partitions), "")
    if batch.problems:
        return batch

    batch.core_hours = round(sum(j.cores * j.hours for j in batch.jobs), 2)
    if max_core_hours is not None and batch.core_hours > max_core_hours:
        batch.problems.append(f"{batch.core_hours:g} core-hours exceed the managed limit of "
                              f"{max_core_hours:g} per batch")
    try:
        decision = conn.budget()
        batch.unit = decision.unit or ""
        batch.used = None if decision.limit is None else float(decision.used or 0.0)
        batch.limit = decision.limit
    except RemoteError as exc:
        batch.problems.append(f"the budget could not be read: {exc}")
        return batch
    batch.requested = round(batch.core_hours if "core" in batch.unit.lower()
                            else sum(j.hours for j in batch.jobs), 2)
    if batch.limit is not None and batch.used + batch.requested > batch.limit:
        batch.problems.append(
            f"the batch may use {batch.requested:g} {batch.unit}, but the session has "
            f"{max(0.0, batch.limit - batch.used):g} of {batch.limit:g} left; ask for less "
            f"time or fewer jobs, or log in again with a larger limit")

    try:
        commit = commit or _install.pc_commit()
    except RemoteError as exc:
        batch.problems.append(str(exc))
        return batch
    batch.commit = commit.sha
    if commit.dirty:
        batch.notes.append(f"the PC's checkout has changes beyond {commit.sha[:12]}; "
                           "the cluster runs the commit without them")
    state = _install.state(conn, commit.sha)
    if state == "complete":
        batch.install = "complete"
    elif state == "running":
        batch.install = "running"
    elif state == "stopped" and auto_install:
        # hpc_guard stopped it; approving this batch starts it again.
        batch.install = "needed"
        batch.notes.append(f"the install of DELFIN {commit.sha[:12]} on {cluster.name} was stopped "
                           "by hpc_guard; it starts again before the batch goes in")
    elif state == "stopped":
        batch.problems.append(f"the install of DELFIN {commit.sha[:12]} on {cluster.name} does not "
                              f"go on: {_install.STOPPED_HINT}; see /clusters (space)")
    elif state == "failed":
        batch.problems.append(f"the install of DELFIN {commit.sha[:12]} on {cluster.name} failed "
                              "before; see /clusters (space) and fix it first")
    elif not auto_install:
        batch.problems.append(f"DELFIN {commit.sha[:12]} is not installed on {cluster.name}, and "
                              "a managed policy turned automatic installs off")
    else:
        sources = _install._clusters.install_sources()
        if not _install.published_source(commit, sources):
            batch.problems.append(
                f"DELFIN {commit.sha[:12]} is not installed on {cluster.name} and on no branch "
                f"of {', '.join(sources) or 'any install source'}; push it there first")
        else:
            batch.install = "needed"
    return batch


def _job(workspace: Path, request: Any, refuse, number: int) -> tuple[Optional[Job], str]:
    if not isinstance(request, dict):
        return None, f"job {number}: not an object"
    raw = str(request.get("folder") or "").strip()
    if not raw:
        return None, f"job {number}: no folder"
    folder = Path(raw) if Path(raw).is_absolute() else workspace / raw
    try:
        folder = folder.resolve(strict=True)
    except (OSError, RuntimeError):
        return None, f"job {number}: no folder {raw}"
    if folder == workspace or workspace not in folder.parents or not folder.is_dir():
        return None, f"job {number}: {raw} is not a folder inside the workspace"
    files, size, problem, paths, fingerprint = _contents(folder, refuse)
    if problem:
        return None, f"job {number}: {problem}"
    kind = str(request.get("kind") or "").strip().lower()
    entry = str(request.get("input") or "").strip()
    if not kind:
        kind = "orca" if entry else "delfin"
    if kind not in ("delfin", "orca"):
        return None, f"job {number}: kind is delfin or orca, not {kind!r}"
    control = folder / "CONTROL.txt"
    if kind == "delfin":
        missing = [f for f in ("input.txt", "CONTROL.txt") if not (folder / f).is_file()]
        if missing:
            return None, f"job {number}: a DELFIN workflow needs {' and '.join(missing)} in {raw}"
        pal, maxcore, hours, named = _control_resources(
            control.read_text(encoding="utf-8", errors="replace"))
        source = "CONTROL.txt" if named else "DELFIN's default (CONTROL.txt names no PAL)"
        entry = ""
    else:
        if not entry or "/" in entry or not entry.endswith(".inp") or not (folder / entry).is_file():
            return None, f"job {number}: an ORCA run names its .inp in {raw} (input), not {entry!r}"
        pal, maxcore = _inp_resources((folder / entry).read_text(encoding="utf-8", errors="replace"))
        source = entry if pal else f"{entry} (no %pal: ORCA runs one process)"
        pal, maxcore, hours = pal or 1, maxcore or DEFAULT_MAXCORE, None
    # The run takes its processes and memory from its own file -- DELFIN
    # from CONTROL.txt (or its defaults), ORCA from the input -- so the
    # allocation follows the file; an override that disagrees would be an
    # allocation the run does not use.
    try:
        needed = pal * maxcore
        cores = int(request["cores"]) if request.get("cores") else pal
        mem_mb = int(float(request["mem_gb"]) * 1024) if request.get("mem_gb") else needed
    except (TypeError, ValueError):
        return None, f"job {number}: cores and mem_gb are numbers"
    if cores != pal:
        return None, (f"job {number}: cores={cores}, but the run uses {pal} processes, as "
                      f"{source} says; change it there, or leave cores out")
    if mem_mb < needed:
        return None, (f"job {number}: mem_gb={mem_mb / 1024:g} is less than the run needs: "
                      f"{pal} processes x maxcore {maxcore} MB = {needed} MB")
    if not 1 <= cores <= 512 or not 100 <= mem_mb <= 4 * 1024 * 1024:
        return None, f"job {number}: {cores} cores and {mem_mb} MB are out of range"
    from delfin.slurm_submit import normalize_time_limit, time_limit_seconds
    given = request.get("time") or (f"{hours:02d}:00:00" if hours else DEFAULT_TIME)
    try:
        limit = normalize_time_limit(given)
        # A batch is judged by its time limits: one without is never asked about.
        seconds = 0 if limit in ("UNLIMITED", "INFINITE") else time_limit_seconds(limit)
    except (ValueError, TypeError) as exc:
        return None, f"job {number}: the time limit {given!r}: {exc}"
    if seconds <= 0:
        return None, f"job {number}: a time limit is needed, not {given!r}"
    name = _SLUG.sub("-", str(request.get("name") or folder.name)).strip("-.")[:40] or "delfin"
    return Job(folder=folder, kind=kind, input=entry, name=name, cores=cores, mem_mb=mem_mb,
               time=limit, hours=round(seconds / 3600, 2), files=files, size=size,
               paths=paths, fingerprint=fingerprint), ""


def _contents(folder: Path, refuse):
    """(files, bytes, problem, paths, fingerprint) of a job folder as it would
    be uploaded: no link, every file past ``refuse``, within the limits; the
    results a watcher brought back into it (RESULTS_MARKER) are no part of it."""
    files = size = 0
    paths, fingerprint = [], []
    for root, dirs, names in os.walk(folder):
        dirs[:] = [d for d in dirs if (Path(root) / d).is_symlink()
                   or not (Path(root) / d / RESULTS_MARKER).is_file()]
        for name in dirs + names:
            path = Path(root) / name
            if path.is_symlink():
                return 0, 0, (f"{path.relative_to(folder)} is a link; uploads copy no links "
                              "(rsync -r): put the file itself there"), (), ()
        for name in names:
            path = Path(root) / name
            if refuse is not None:
                why = refuse(path)
                if why:
                    return 0, 0, f"{path.relative_to(folder)}: {why}", (), ()
            st = path.stat()
            rel = path.relative_to(folder).as_posix()
            files += 1
            size += st.st_size
            paths.append(rel)
            fingerprint.append((rel, st.st_size, st.st_mtime_ns))
            if files > MAX_FILES or size > MAX_BYTES:
                return 0, 0, (f"{folder.name} holds more than {MAX_FILES} files or "
                              f"{MAX_BYTES // (1024 * 1024)} MB; a job folder carries its inputs"), \
                    (), ()
    if not files:
        return 0, 0, f"{folder.name} is empty", (), ()
    return files, size, "", tuple(sorted(paths)), tuple(sorted(fingerprint))


def _control_resources(text: str) -> tuple[int, int, Optional[int], bool]:
    """PAL and maxcore as the run uses them -- CONTROL.txt's, else DELFIN's
    defaults (the template's, as delfin.config.read_control_file fills them)
    -- its job_timeout_hours, and whether CONTROL.txt names PAL."""
    from delfin.agent.job_fix import parse_pal_maxcore
    pal, maxcore = parse_pal_maxcore(text)
    found = re.search(r"^\s*job_timeout_hours\s*=\s*(\d+)", text, flags=re.MULTILINE)
    defaults: dict = {}
    if not (pal and maxcore):
        from delfin.config import _load_template_defaults
        defaults = _load_template_defaults()
    return (int(pal or defaults.get("PAL") or 1),
            int(maxcore or defaults.get("maxcore") or DEFAULT_MAXCORE),
            int(found.group(1)) if found else None, bool(pal))


def _inp_resources(text: str) -> tuple[Optional[int], Optional[int]]:
    from delfin.dashboard.input_processing import parse_inp_resources
    return parse_inp_resources(text)


def _partition(cluster, asked: Any) -> tuple[str, str]:
    name = str(asked or cluster.default_partition or "").strip()
    if not name:
        return "", (f"no partition: name one, or set default_partition for {cluster.name} in "
                    "remote.clusters")
    if not _PARTITION.fullmatch(name):
        return "", f"{name!r} is no partition name"
    if cluster.allowed_partitions and name not in cluster.allowed_partitions:
        return "", f"partition {name} is not one of {list(cluster.allowed_partitions)}"
    return name, ""


def submit(conn, batch: Batch, *, session_id: str = "", workspace: Any = None,
           on_progress: Optional[Callable[[str], None]] = None,
           install_wait_s: float = INSTALL_WAIT_S) -> list[dict]:
    """Submit an approved ``batch``: install first when needed, then each job.
    The submitted jobs (``job_id``, ``name``, ``remote_dir``, ``folder``)."""
    if not batch.ok:
        raise RemoteError("a batch with problems is not submitted: " + "; ".join(batch.problems))
    say = on_progress or (lambda message: None)
    install_dir = f"{conn.cluster.install_root.rstrip('/')}/{batch.commit}"
    if batch.install != "complete":
        commit = _install.Commit(sha=batch.commit, checkout=_install.pc_commit().checkout,
                                 dirty=False)
        outcome = _install.ensure_installed(conn, commit=commit, wait_seconds=install_wait_s,
                                            on_progress=say)
        if outcome.state != "complete":
            raise RemoteError(f"DELFIN {batch.commit[:12]} is still installing on "
                              f"{conn.cluster.name}; nothing was submitted -- submit again once "
                              "it is installed")
    # What was approved is what goes: each folder as the dialog counted it,
    # and the budget as it stands now -- a dialog and an install take time.
    for job in batch.jobs:
        changed = _changed_since_approval(job)
        if changed:
            raise RemoteError(f"{job.folder.name} changed after the batch was approved "
                              f"({changed}); nothing was submitted -- submit it again to be "
                              "asked about it as it is now")
    if batch.limit is not None:
        decision = conn.budget()
        if (decision.limit is not None
                and float(decision.used or 0) + batch.requested > float(decision.limit)):
            raise RemoteError(
                f"the budget changed since the dialog: the batch may use {batch.requested:g} "
                f"{batch.unit}, {max(0.0, float(decision.limit) - float(decision.used or 0)):g} "
                f"of {float(decision.limit):g} are left; nothing was submitted")
    stamp = time.strftime("%Y%m%d-%H%M%S")
    done: list[dict] = []
    for job in batch.jobs:
        try:
            done.append(_submit_one(conn, batch, job, stamp=stamp, install_dir=install_dir,
                                    session_id=session_id, workspace=workspace, say=say))
        except RemoteError as exc:
            if not done:
                raise
            ids = ", ".join(d["job_id"] for d in done)
            raise PartialSubmit(f"{exc} -- before it, job(s) {ids} were submitted and run: "
                                "submit only the rest again", done) from exc
    return done


def _submit_one(conn, batch: Batch, job: Job, *, stamp: str, install_dir: str,
                session_id: str, workspace: Any, say) -> dict:
    """Upload one approved job folder -- exactly its approved files -- and
    sbatch it; its record."""
    remote_dir = (f"{conn.cluster.work_dir.rstrip('/')}/{stamp}-{secrets.token_hex(2)}-"
                  f"{job.name}")
    relative = remote_dir[2:] if remote_dir.startswith("~/") else remote_dir
    say(f"uploading {job.folder.name} to {conn.cluster.name}:{remote_dir}")
    made = conn.run(["mkdir", "-p", relative])
    if not made.ok:
        raise RemoteError(f"cannot create {remote_dir}: {(made.stderr or '').strip()}")
    fd, names = tempfile.mkstemp(prefix="delfin-upload-", suffix=".list")
    try:
        with os.fdopen(fd, "wb") as fh:
            fh.write(b"".join(p.encode("utf-8", "surrogateescape") + b"\0" for p in job.paths))
        sent = conn.upload(f"{job.folder}/", f"{relative}/",
                           extra=["--from0", f"--files-from={names}"])
    finally:
        os.unlink(names)
    if not sent.ok:
        raise RemoteError(f"uploading {job.folder.name} failed: {(sent.stderr or '').strip()}")
    options = [f"--job-name={job.name}", f"--partition={batch.partition}", "--ntasks=1",
               f"--cpus-per-task={job.cores}", f"--mem={job.mem_mb}M", f"--time={job.time}",
               f"--signal=B:USR1@{_usr1_offset(job.hours)}",
               f"--comment={_install.job_comment(batch.commit)}"]
    result = conn.run(["bash", "-c", _SUBMIT, "delfin-submit", remote_dir, install_dir,
                       job.kind, job.input, job.name, *options])
    said = (result.stdout or "").strip().splitlines() if result.ok else []
    job_id = said[-1].split(";")[0] if said else ""
    if not _JOB_ID.fullmatch(job_id):
        raise RemoteError(f"sbatch did not take {job.name}: "
                          f"{(result.stderr or result.stdout or '').strip()}")
    record = {"cluster": conn.cluster.name, "job_id": job_id, "name": job.name,
              "kind": job.kind, "input": job.input, "remote_dir": remote_dir,
              "folder": str(job.folder), "workspace": str(workspace or ""),
              "session_id": session_id, "commit": batch.commit,
              "submitted_at": time.time()}
    _record(record)
    say(f"{job.name} is job {job_id} on {conn.cluster.name}")
    return record


def _changed_since_approval(job: Job) -> str:
    """What differs in a job folder from the files approved; "" for nothing."""
    _, _, problem, _, fingerprint = _contents(job.folder, None)
    if problem:
        return problem
    if fingerprint == job.fingerprint:
        return ""
    before = {p: (s, m) for p, s, m in job.fingerprint}
    now = {p: (s, m) for p, s, m in fingerprint}
    parts = []
    for what, names in (("new", sorted(set(now) - set(before))),
                        ("gone", sorted(set(before) - set(now))),
                        ("changed", sorted(p for p in set(now) & set(before)
                                           if now[p] != before[p]))):
        if names:
            parts.append(f"{what}: {', '.join(names[:3])}"
                         + (f" and {len(names) - 3} more" if len(names) > 3 else ""))
    return "; ".join(parts)


def _memory(mb: int) -> str:
    """A job's memory as the dialog shows it: the MB sbatch gets below 10 GB,
    else GB to one decimal (2000 MB is not "1 GB")."""
    return f"{mb} MB" if mb < 10 * 1024 else f"{mb / 1024:.1f} GB"


def _usr1_offset(hours: float) -> int:
    """The dashboard's: USR1 (a final sync) half the time limit before its
    end, at least 30 s and at most 300 s ahead."""
    return max(30, min(300, int(hours * 3600) // 2))


# --- a batch that waits for its install -------------------------------------------


#: What a kept batch is: waiting for its install, submitted (job ids kept),
#: or not submitted (the reason kept).
PENDING_STATES = ("waiting", "submitted", "not_submitted")


def batch_to_dict(batch: Batch) -> dict:
    """An approved batch as it is kept until its install is done."""
    return {"cluster": batch.cluster, "target": batch.target, "partition": batch.partition,
            "commit": batch.commit, "install": batch.install, "unit": batch.unit,
            "requested": batch.requested, "core_hours": batch.core_hours,
            "used": batch.used, "limit": batch.limit, "notes": list(batch.notes),
            "jobs": [{"folder": str(j.folder), "kind": j.kind, "input": j.input,
                      "name": j.name, "cores": j.cores, "mem_mb": j.mem_mb, "time": j.time,
                      "hours": j.hours, "files": j.files, "size": j.size,
                      "paths": list(j.paths),
                      "fingerprint": [list(f) for f in j.fingerprint]} for j in batch.jobs]}


def batch_from_dict(data: dict) -> Batch:
    """The batch ``batch_to_dict`` kept: the jobs, their approved files and
    the fingerprints the files are checked against before sbatch."""
    batch = Batch(cluster=str(data["cluster"]), target=str(data.get("target") or ""),
                  partition=str(data["partition"]), commit=str(data["commit"]),
                  install=str(data.get("install") or ""), unit=str(data.get("unit") or ""),
                  requested=float(data.get("requested") or 0.0),
                  core_hours=float(data.get("core_hours") or 0.0),
                  used=None if data.get("used") is None else float(data["used"]),
                  limit=None if data.get("limit") is None else float(data["limit"]),
                  notes=list(data.get("notes") or []))
    for j in data.get("jobs") or []:
        batch.jobs.append(Job(
            folder=Path(j["folder"]), kind=str(j["kind"]), input=str(j.get("input") or ""),
            name=str(j["name"]), cores=int(j["cores"]), mem_mb=int(j["mem_mb"]),
            time=str(j["time"]), hours=float(j["hours"]), files=int(j["files"]),
            size=int(j["size"]), paths=tuple(str(p) for p in j.get("paths") or ()),
            fingerprint=tuple((str(p), int(sz), int(mt))
                              for p, sz, mt in j.get("fingerprint") or ())))
    return batch


def queue_after_install(conn, batch: Batch, *, session_id: str = "", workspace: Any = None,
                        on_progress: Optional[Callable[[str], None]] = None) -> dict:
    """An approved ``batch`` whose commit is not installed yet: the install
    started (detached; one already running is left to run) and the batch
    kept for the job watcher to submit once it is complete. Returns the
    kept entry (``id``, ``cluster``, ``commit``, ``jobs``). Raises what the
    install raises (CommitNotPublished, InstallFailed) before keeping
    anything."""
    if not batch.ok:
        raise RemoteError("a batch with problems is not submitted: " + "; ".join(batch.problems))
    if batch.install == "running":
        current = _install.state(conn, batch.commit)
        if current not in ("running", "complete"):
            batch.install = "needed"     # it ended since the dialog: start it again
    if batch.install == "needed":
        commit = _install.Commit(sha=batch.commit, checkout=_install.pc_commit().checkout,
                                 dirty=False)
        _install.ensure_installed(conn, commit=commit, wait=False, on_progress=on_progress)
    entry = {"id": secrets.token_hex(4), "cluster": conn.cluster.name, "commit": batch.commit,
             "session_id": session_id, "workspace": str(workspace or ""),
             "queued_at": time.time(), "state": "waiting",
             "jobs": [j.name for j in batch.jobs], "batch": batch_to_dict(batch)}
    with _pending_locked() as kept:
        kept[entry["id"]] = entry
    return entry


def pending(*, cluster: Optional[str] = None, state: str = "waiting") -> list[dict]:
    """The batches kept for their install, oldest first (``state`` "" for all)."""
    with _pending_locked(write=False) as kept:
        out = [dict(e) for e in kept.values()
               if (not state or e.get("state") == state)
               and (cluster is None or e.get("cluster") == cluster)]
    return sorted(out, key=lambda e: float(e.get("queued_at") or 0))


def settle_pending(entry_id: str, **changes: Any) -> None:
    """Record what became of a kept batch (its state, job ids, reason)."""
    with _pending_locked() as kept:
        if entry_id in kept:
            kept[entry_id].update(changes)


def _pending_path() -> Path:
    # Resolved per call, so the suite's state redirect reaches it (state_paths).
    return Path.home() / ".delfin" / "remote_pending.json"


@contextlib.contextmanager
def _pending_locked(write: bool = True):
    """The kept batches, read and (``write``) saved under one lock: the tool
    that keeps a batch and the watcher that settles it are other processes."""
    import fcntl
    from delfin.agent import state_paths
    path = _pending_path()
    state_paths.ensure_dir(path.parent)
    fd = os.open(str(path.with_suffix(".lock")), os.O_RDWR | os.O_CREAT, 0o600)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX)
        try:
            kept = json.loads(path.read_text(encoding="utf-8"))
            if not isinstance(kept, dict):
                kept = {}
        except (OSError, ValueError):
            kept = {}
        yield kept
        if write:
            fd2, tmp = tempfile.mkstemp(prefix=".remote_pending-", dir=str(path.parent))
            try:
                with os.fdopen(fd2, "w", encoding="utf-8") as fh:
                    json.dump(kept, fh, indent=1, sort_keys=True)
                os.replace(tmp, path)
            except BaseException:
                with contextlib.suppress(OSError):
                    os.unlink(tmp)
                raise
    finally:
        fcntl.flock(fd, fcntl.LOCK_UN)
        os.close(fd)


def cancellable(conn, job_ids: Iterable[str]) -> list[str]:
    """The ids, when each is an active job of the user's on this cluster;
    else JobRefused for the whole call (checked before anyone is asked)."""
    ids = [str(j).strip() for j in job_ids]
    if not ids or any(not _JOB_ID.fullmatch(j) for j in ids):
        raise JobRefused("give Slurm job ids (digits), at least one")
    mine = {str(j.jobid) for j in conn.jobs() if getattr(j, "source", "") == "active"}
    foreign = [j for j in ids if j not in mine]
    if foreign:
        raise JobRefused(f"{', '.join(foreign)}: none of your active jobs on {conn.cluster.name}")
    return ids


def cancel(conn, job_ids: Iterable[str]) -> list[str]:
    """Cancel these of the user's active jobs (``cancellable``, checked again
    now); the ids cancelled."""
    ids = cancellable(conn, job_ids)
    result = conn.run(["scancel", *ids])
    if not result.ok:
        raise RemoteError(f"scancel failed: {(result.stderr or '').strip()}")
    return ids


def records(*, cluster: Optional[str] = None) -> list[dict]:
    """The jobs DELFIN submitted, oldest first (for the job watcher)."""
    path = _records_path()
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError:
        return []
    out = []
    for line in lines:
        try:
            rec = json.loads(line)
        except ValueError:
            continue
        if isinstance(rec, dict) and (cluster is None or rec.get("cluster") == cluster):
            out.append(rec)
    return out


def _records_path() -> Path:
    # Resolved per call, so the suite's state redirect reaches it (state_paths).
    return Path.home() / ".delfin" / "remote_jobs.jsonl"


def _record(record: dict) -> None:
    from delfin.agent import state_paths
    path = _records_path()
    state_paths.ensure_dir(path.parent)
    with state_paths.open_append(path) as fh:
        fh.write(json.dumps(record) + "\n")
