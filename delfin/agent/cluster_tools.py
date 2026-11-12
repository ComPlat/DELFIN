"""The agent's tools for DELFIN's clusters: status, reading, transfers, jobs, a shell.

Input: one tool call's arguments, the session's permissions
(``KitToolPermissions``) and the executor (``api_client._DocToolExecutor``),
whose gates these reuse. Output: the tool result the model reads, JSON;
whatever came from a cluster is fenced as untrusted text.

Semantics, tool by tool:

* ``cluster_status``: every cluster of remote.clusters -- logged in or not,
  why one is refused, budget, jobs -- or one cluster's jobs and whether the
  PC's commit is installed there. Unasked.
* ``cluster_read``: a file's text (its first ``READ_BYTES``, or its last
  ``tail`` lines) or a folder's listing. Unasked inside the cluster's work
  dir, asked anywhere else; never ``~/.ssh``, ``~/.config``, ``~/.gnupg``,
  ``~/.delfin``, the shell histories and the like there, a path the secret
  deny list names, nor another user's home. Paths are judged where they
  really lead (``realpath`` on the cluster).
* ``cluster_fetch``: a file or folder from the cluster into the workspace
  (``dest``, default ``cluster/<name>/<its name>``), at most ``FETCH_BYTES``,
  never over an existing file, and only where write_file may write (the same
  gate: asked in the default mode, and under diff_approval, which stages
  nothing for a fetch). Reading is unasked from the work dir, asked
  elsewhere, with the same refusals as ``cluster_read``. No link comes along
  (rsync -r copies none).
* ``cluster_upload``: workspace files to a folder under the work dir, judged
  where it really leads (a link there would carry them elsewhere). Asked;
  links, secrets and DELFIN's masked state are refused; the transfer is
  recorded as an egress event.
* ``cluster_submit``: DELFIN workflows and single ORCA runs
  (``delfin.remote.jobs``). A batch that does not fit -- budget, partition,
  install, what is in its folders -- is refused before anyone is asked; one
  that fits is always asked, in every mode, one dialog per batch.
* ``cluster_cancel``: jobs of the user's that are queued or running, which
  is checked before the dialog. Asked.
* ``cluster_run``: one command on the login node. Asked every time, in every
  mode, after the bash gate's checks that mean something on another machine
  (deny patterns, secret paths, credential prints, the user's rules for
  ``bash`` and for ``cluster_run``, a refusal earlier this session; an
  outbound transfer is named in the dialog). A command that would submit --
  sbatch, srun, salloc, also inside ``bash -c`` or after env/timeout/xargs,
  as ctun's own guard reads it -- or cancel is refused and pointed to
  cluster_submit / cluster_cancel, which check the budget and ask.

* ``cluster_propose``: a cluster for DELFIN, as a NEEDS YOU item that opens
  the window's add form filled in (``delfin.remote.setup.propose``). It
  writes no setting and no ctun entry -- only the person's form does --
  and is advertised with or without a cluster configured. Unasked: the
  item itself is the question.

Every call: refused for a cluster the managed policy does not allow, and,
when nobody is logged in, answered ``login_required`` with a login item in
NEEDS YOU (the window's masked field). Each call that changes something, or
reads outside the work dir, and each refused cluster_run writes an audit row:
cluster, command or files, job ids, approval and outcome.
"""

from __future__ import annotations

import json
import os
import re
from pathlib import Path, PurePosixPath
from typing import Any, Optional

#: How much of a file cluster_read returns, how many lines ``tail`` may ask
#: for, and the most a fetch, an upload and a shell's output may carry.
READ_BYTES = 64 * 1024
TAIL_MAX = 2000
FETCH_BYTES = 200 * 1024 * 1024
UPLOAD_BYTES = 200 * 1024 * 1024
UPLOAD_FILES = 2000
OUTPUT_BYTES = 32 * 1024
RUN_TIMEOUT = 120
RUN_TIMEOUT_MAX = 600

#: Under the cluster's home: never read or fetched, whoever asks.
_PRIVATE = (".ssh", ".config", ".gnupg", ".delfin", ".delfin_settings.json", ".netrc",
            ".aws", ".kube", ".docker", ".git-credentials", ".pgpass", ".bash_history",
            ".zsh_history", ".python_history")

#: Left out of a fetched folder wherever they are in it, beside _PRIVATE.
_KEYS = ("*.pem", "*.key", "*.p12", "*.pfx", "id_rsa*", "id_dsa*", "id_ecdsa*", "id_ed25519*",
         ".env", ".env.*", ".npmrc", ".pypirc")

#: The submitting commands cluster_run never runs, beside the cluster's own
#: guard_commands, and the one that cancels.
_SUBMIT = ("sbatch", "srun", "salloc")
_CANCEL = ("scancel",)

#: A destination under the work dir, as uploads write it.
_DEST = re.compile(r"[A-Za-z0-9._/-]+")

#: Runs on the cluster: its home, the work dir and the path, each resolved,
#: and what the path is. The path must exist unless the third argument is
#: -m (an upload's destination). "~/" arrives unexpanded and is expanded here.
_RESOLVE = r'''home="$(cd ~ && pwd -P)"
case "$1" in "~") w="$home" ;; "~/"*) w="$home/${1#"~/"}" ;; *) w="$1" ;; esac
case "$2" in "~") p="$home" ;; "~/"*) p="$home/${2#"~/"}" ;; /*) p="$2" ;; *) p="$home/$2" ;; esac
printf '%s\n' "$home"
realpath -m -- "$w"
realpath "${3:--e}" -- "$p" || exit 3
if [ -d "$p" ]; then echo dir; elif [ -e "$p" ]; then echo file; else echo none; fi
'''


# --- the tools ---------------------------------------------------------------------


def status(executor, perms, *, cluster: Any = None) -> str:
    from delfin.agent import managed_policy
    if not cluster:
        from delfin.agent.service import clusters_feed
        rows = [r for r in clusters_feed.clusters()["clusters"]
                if managed_policy.cluster_allowed(r["name"])]
        if not rows:
            return _error("no_clusters", "no cluster is configured: remote.clusters in "
                          "~/.delfin_settings.json names none")
        return _fenced("cluster status", {"clusters": rows, "hint": (
            "logged_out: the user logs in (/clusters in the DELFIN window, or the login item); "
            "refused: the problems say what the user must fix")})
    conn, refusal = _connect(perms, cluster, "cluster_status")
    if refusal:
        return refusal
    from delfin.remote import install
    from delfin.remote.errors import LoginRequired, RemoteError
    try:
        if not conn.is_live():
            return _login_needed(perms, conn.name, "cluster_status")
        listing = conn.jobs()
        decision = conn.budget()
        try:
            commit = install.pc_commit()
            installed = install.state(conn, commit.sha)
        except RemoteError as exc:
            commit, installed = None, f"unknown ({exc})"
    except LoginRequired:
        return _login_needed(perms, conn.name, "cluster_status")
    except RemoteError as exc:
        return _error(exc.code, str(exc))
    active = [{"id": j.jobid, "name": j.name, "state": j.state, "partition": j.partition,
               "elapsed": j.elapsed, "limit": j.timelimit} for j in listing
              if j.source == "active"]
    return _fenced(f"cluster {conn.name} status", {
        "cluster": conn.name, "target": conn.cluster.target, "logged_in": True,
        "budget": None if decision.limit is None else {
            "used": decision.used, "limit": decision.limit, "unit": decision.unit},
        "active_jobs": active, "finished_recently": len(listing) - len(active),
        "work_dir": conn.cluster.work_dir, "default_partition": conn.cluster.default_partition,
        "delfin": {"commit": commit.sha[:12] if commit else "", "install": installed}})


def read_path(executor, perms, *, cluster: Any, path: Any, tail: Any = None) -> str:
    conn, refusal = _connect(perms, cluster, "cluster_read")
    if refusal:
        return refusal
    place, refusal = _place(executor, perms, conn, path, "cluster_read", {"cluster": cluster,
                                                                           "path": path})
    if refusal:
        return refusal
    real, kind = place
    if kind == "dir":
        got = conn.run(["find", real, "-mindepth", "1", "-maxdepth", "1", "-printf",
                        "%y %s %TY-%Tm-%Td %TH:%TM %f\\n"], timeout=60)
        if not got.ok:
            return _error("read_failed", (got.stderr or "").strip()[:300])
        entries = sorted((got.stdout or "").splitlines())[:500]
        return _fenced(f"{conn.name}:{real}", {"cluster": conn.name, "path": real,
                                               "entries": entries})
    try:
        lines = int(tail) if tail not in (None, "") else 0
    except (TypeError, ValueError):
        return _error("bad_argument", "tail is a number of lines")
    if lines:
        lines = max(1, min(lines, TAIL_MAX))
        got = conn.run(["bash", "-c", 'tail -n "$1" -- "$2" | head -c "$3"', "delfin-read",
                        str(lines), real, str(READ_BYTES)], timeout=60)
    else:
        got = conn.run(["head", "-c", str(READ_BYTES + 1), "--", real], timeout=60)
    if not got.ok:
        return _error("read_failed", (got.stderr or "").strip()[:300])
    text = got.stdout or ""
    if "\x00" in text:
        return _error("binary", f"{real} is not text; cluster_fetch copies it into the workspace")
    truncated = not lines and len(text.encode("utf-8", "replace")) > READ_BYTES
    return _fenced(f"{conn.name}:{real}", {"cluster": conn.name, "path": real,
                                           "truncated": truncated,
                                           "text": text[:READ_BYTES]})


def fetch(executor, perms, *, cluster: Any, path: Any, dest: Any = None) -> str:
    from delfin.remote.errors import LoginRequired, RemoteError
    conn, refusal = _connect(perms, cluster, "cluster_fetch")
    if refusal:
        return refusal
    place, refusal = _place(executor, perms, conn, path, "cluster_fetch",
                            {"cluster": cluster, "path": path, "dest": dest})
    if refusal:
        return refusal
    real, kind = place
    target = str(dest or f"cluster/{conn.name}/{PurePosixPath(real).name}")
    local, error = executor._resolve_in_workspace(target, perms)
    if error:
        return _error("bad_destination", error)
    rel = os.path.relpath(local, perms.workspace)
    if local.exists():
        return _error("exists", f"{rel} exists already; name another dest")
    try:
        sized = conn.run(["du", "-sb", "--", real], timeout=120)
    except LoginRequired:
        return _login_needed(perms, conn.name, "cluster_fetch")
    except RemoteError as exc:
        return _error(exc.code, str(exc))
    try:
        size = int((sized.stdout or "0").split()[0])
    except (ValueError, IndexError):
        return _error("fetch_failed", f"the size of {real} could not be read")
    if size > FETCH_BYTES:
        return _error("too_large", f"{real} holds {size // (1024 * 1024)} MB; a fetch carries at "
                      f"most {FETCH_BYTES // (1024 * 1024)} MB (read parts with cluster_read)")
    # What write_file may write, a fetch may: the same gate, which asks in
    # the default mode. diff_approval stages a write tool's change for the
    # user; a fetch stages nothing, so it is asked there instead.
    from delfin.agent import api_client as A
    args = {"cluster": conn.name, "from": real, "path": target, "size": _size(size)}
    refused = executor._gate_write_path(target, perms, "cluster_fetch", args)
    if refused:
        return _error("denied" if refused.startswith("user denied") else "approval_timed_out"
                      if "TIMED OUT" in refused else "bad_destination", refused)
    if A.effective_mode(perms) == "diff_approval":
        refusal = _ask(executor, perms, "cluster_fetch", args,
                       executor._build_change_preview("cluster_fetch", args, local))
        if refusal:
            return refusal
    local.parent.mkdir(parents=True, exist_ok=True)
    try:
        # A folder's own copies of private names and keys stay behind too.
        got = (conn.download(f"{real}/", f"{local}/", extra=[
                   f"--exclude={name}" for name in _PRIVATE + _KEYS]) if kind == "dir"
               else conn.download(real, str(local)))
    except LoginRequired:
        return _login_needed(perms, conn.name, "cluster_fetch")
    except RemoteError as exc:
        return _error(exc.code, str(exc))
    _audit(perms, "cluster_fetch", "ok" if got.ok else "error", path=str(local),
           extra={"cluster": conn.name, "from": real, "bytes": size})
    if not got.ok:
        return _error("fetch_failed", (got.stderr or "").strip()[:300])
    return json.dumps({"cluster": conn.name, "from": real, "to": rel, "bytes": size})


def upload(executor, perms, *, cluster: Any, src: Any, dest: Any = None) -> str:
    from delfin.remote.errors import LoginRequired, RemoteError
    conn, refusal = _connect(perms, cluster, "cluster_upload")
    if refusal:
        return refusal
    local, error = executor._resolve_in_workspace(str(src or ""), perms, for_read=True)
    if error:
        return _error("bad_source", error)
    if not local.exists():
        return _error("bad_source", f"no {src} in the workspace")
    files, why = _upload_files(executor, perms, local)
    if why:
        return _error("refused", why)
    remote_rel, why = _under_work_dir(conn, dest or local.name)
    if why:
        return _error("bad_destination", why)
    # Where it really leads: a link under the work dir would carry the files
    # elsewhere (rsync follows a linked destination folder).
    try:
        got = conn.run(["bash", "-c", _RESOLVE, "delfin-resolve", conn.cluster.work_dir,
                        remote_rel, "-m"], timeout=60)
    except LoginRequired:
        return _login_needed(perms, conn.name, "cluster_upload")
    except RemoteError as exc:
        return _error(exc.code, str(exc))
    lines = (got.stdout or "").splitlines()
    if not got.ok or len(lines) < 4:
        return _error("upload_failed", f"{remote_rel} could not be resolved on {conn.name}")
    _, work, real, _ = lines[:4]
    if not _inside(real, work) or real == work:
        return _error("bad_destination", f"{remote_rel} leads to {real} on {conn.name}, not "
                      f"into the work dir {work}; uploads stay in it")
    if not _DEST.fullmatch(real):
        return _error("bad_destination", f"{real} on {conn.name} holds characters a transfer "
                      "would split; set a plainer work_dir")
    size = sum(f.stat().st_size for f in files)
    shown = [str(f.relative_to(local)) if local.is_dir() else f.name for f in files]
    args = {"cluster": conn.name, "src": str(src), "dest": remote_rel}
    preview = (f"upload {len(files)} file{'s' if len(files) != 1 else ''} "
               f"({_size(size)}) to {conn.name}:{remote_rel}\nfrom {local}\n"
               + "\n".join(f"  {s}" for s in shown[:20])
               + (f"\n  … {len(shown) - 20} more" if len(shown) > 20 else ""))
    refusal = _ask(executor, perms, "cluster_upload", args, preview)
    from delfin.agent.api_client import _record_security_event
    _record_security_event("egress", "cluster_upload",
                           f"{len(files)} files to {conn.name}", blocked=bool(refusal))
    said = f"upload {len(files)} files ({_size(size)}) to {conn.name}:{remote_rel}"
    if refusal:
        _audit(perms, "cluster_upload", "denied", command=said, path=str(local),
               extra={"cluster": conn.name, "dest": remote_rel, "files": len(files)})
        return refusal
    try:
        conn.run(["mkdir", "-p", "--", str(PurePosixPath(real).parent)])
        got = (conn.upload(f"{local}/", f"{real}/") if local.is_dir()
               else conn.upload(str(local), real))
    except LoginRequired:
        return _login_needed(perms, conn.name, "cluster_upload")
    except RemoteError as exc:
        return _error(exc.code, str(exc))
    _audit(perms, "cluster_upload", "ok" if got.ok else "error", command=said, path=str(local),
           extra={"cluster": conn.name, "dest": remote_rel, "files": len(files),
                  "bytes": size, "approved": True})
    if not got.ok:
        return _error("upload_failed", (got.stderr or "").strip()[:300])
    return json.dumps({"cluster": conn.name, "to": remote_rel, "files": len(files),
                       "bytes": size})


def submit(executor, perms, *, cluster: Any, jobs: Any) -> str:
    from delfin.agent import managed_policy
    from delfin.remote import jobs as remote_jobs
    from delfin.remote.errors import LoginRequired, RemoteError
    conn, refusal = _connect(perms, cluster, "cluster_submit")
    if refusal:
        return refusal
    if not isinstance(jobs, list):
        return _error("bad_argument", "jobs is a list of {folder, kind, input, cores, mem_gb, "
                      "time, partition, name}")
    try:
        if not conn.is_live():
            return _login_needed(perms, conn.name, "cluster_submit")
        batch = remote_jobs.plan(conn, perms.workspace, jobs,
                                 refuse=_file_refusal(executor, perms),
                                 max_core_hours=managed_policy.cluster_core_hours_ceiling(),
                                 auto_install=managed_policy.cluster_auto_install())
    except LoginRequired:
        return _login_needed(perms, conn.name, "cluster_submit")
    except RemoteError as exc:
        return _error(exc.code, str(exc))
    if not batch.ok:
        # As it was asked: the jobs that passed the checks are not the batch.
        asked = [str(j.get("folder") or "") if isinstance(j, dict) else "?" for j in jobs]
        _audit(perms, "cluster_submit", "denied",
               command=f"submit {len(jobs)} job(s) to {conn.name}: " + ", ".join(asked),
               reason="; ".join(batch.problems)[:300],
               extra={"cluster": conn.name, "folders": asked})
        return _error("batch_refused", "the batch was not submitted, and nobody was asked: "
                      + "; ".join(batch.problems))
    folders = [str(j.folder) for j in batch.jobs]
    said = f"submit {len(batch.jobs)} job(s) to {conn.name}: " + ", ".join(
        os.path.basename(f) for f in folders)
    args = {"cluster": conn.name, "jobs": jobs}
    refusal = _ask(executor, perms, "cluster_submit", args, batch.preview())
    if refusal:
        _audit(perms, "cluster_submit", "denied", command=said,
               extra={"cluster": conn.name, "folders": folders})
        return refusal
    try:
        done = remote_jobs.submit(conn, batch, session_id=getattr(perms, "task_session_id", ""),
                                  workspace=perms.workspace)
    except LoginRequired:
        return _login_needed(perms, conn.name, "cluster_submit")
    except remote_jobs.PartialSubmit as exc:
        # The jobs submitted before the failure run: watched like any, and
        # named, so that they are not submitted twice.
        ids = [d["job_id"] for d in exc.done]
        _audit(perms, "cluster_submit", "error", command=f"{said} -> jobs {', '.join(ids)}",
               reason=str(exc)[:300], extra={"cluster": conn.name, "folders": folders,
                                             "approved": True, "job_ids": ids})
        _watch_submitted(perms, conn.name, exc.done)
        return json.dumps({"error": str(exc), "code": exc.code, "cluster": conn.name,
                           "submitted": _submitted(perms, exc.done),
                           "next": ("the submitted jobs run and are watched: do NOT submit "
                                    "them again. Tell the user which were not submitted, and "
                                    "why.")}, ensure_ascii=False)
    except RemoteError as exc:
        _audit(perms, "cluster_submit", "error", command=said, reason=str(exc)[:300],
               extra={"cluster": conn.name, "folders": folders, "approved": True})
        return _error(exc.code, str(exc))
    said += " -> jobs " + ", ".join(d["job_id"] for d in done)
    _audit(perms, "cluster_submit", "ok", command=said, extra={
        "cluster": conn.name, "folders": folders, "approved": True,
        "job_ids": [d["job_id"] for d in done], "commit": batch.commit,
        "core_hours": batch.core_hours})
    watched = _watch_submitted(perms, conn.name, done)
    return json.dumps({"cluster": conn.name, "submitted": _submitted(perms, done),
        "notes": batch.notes,
        "next": (("DELFIN's job watcher follows these jobs: when one ends, its key results "
                  "come back to results_to and you are told in a later turn. END your turn; "
                  "do not poll with cluster_status.") if watched else
                 "cluster_status shows the jobs; cluster_read and cluster_fetch reach their "
                 "folders (remote_dir)")})


def _submitted(perms, done: list[dict]) -> list[dict]:
    """The submitted jobs as the model reads them: id, name, both folders and
    where the watcher brings the results."""
    from delfin.remote import watch
    return [{"job_id": d["job_id"], "name": d["name"], "remote_dir": d["remote_dir"],
             "folder": os.path.relpath(d["folder"], perms.workspace),
             "results_to": os.path.relpath(watch.results_dir(d), perms.workspace)}
            for d in done]


def _watch_submitted(perms, cluster: str, done: list[dict]) -> bool:
    """Each submitted job as a watch of this session's (job_monitor's
    ``cluster`` kind): the watcher's outcome reaches the session's next turn,
    or wakes it. Whether all were registered; never raises."""
    from delfin.agent import job_monitor
    ok = True
    for job in done:
        what = "ORCA run" if job.get("kind") == "orca" else "DELFIN workflow"
        try:
            job_monitor.register_agent_job(
                perms.workspace, f"{cluster}:{job['job_id']}",
                description=f"{job['name']}: {what} on {cluster}",
                extra={"kind": "cluster", "cluster": cluster, "remote_id": str(job["job_id"]),
                       "session_id": str(getattr(perms, "task_session_id", "") or "")})
        except Exception:
            ok = False
    return ok


def cancel(executor, perms, *, cluster: Any, job_ids: Any) -> str:
    from delfin.remote import jobs as remote_jobs
    from delfin.remote.errors import LoginRequired, RemoteError
    conn, refusal = _connect(perms, cluster, "cluster_cancel")
    if refusal:
        return refusal
    try:
        ids = remote_jobs.cancellable(conn, job_ids if isinstance(job_ids, list) else [])
    except LoginRequired:
        return _login_needed(perms, conn.name, "cluster_cancel")
    except RemoteError as exc:
        return _error(exc.code, str(exc))
    args = {"cluster": conn.name, "job_ids": ids}
    refusal = _ask(executor, perms, "cluster_cancel", args,
                   f"cancel job{'s' if len(ids) != 1 else ''} {', '.join(ids)} on {conn.name} "
                   f"({conn.cluster.target})")
    said = f"scancel {' '.join(ids)} on {conn.name}"
    if refusal:
        _audit(perms, "cluster_cancel", "denied", command=said,
               extra={"cluster": conn.name, "job_ids": ids})
        return refusal
    try:
        cancelled = remote_jobs.cancel(conn, ids)
    except LoginRequired:
        return _login_needed(perms, conn.name, "cluster_cancel")
    except RemoteError as exc:
        _audit(perms, "cluster_cancel", "error", command=said, reason=str(exc)[:300],
               extra={"cluster": conn.name, "job_ids": ids, "approved": True})
        return _error(exc.code, str(exc))
    _audit(perms, "cluster_cancel", "ok", command=said,
           extra={"cluster": conn.name, "job_ids": cancelled, "approved": True})
    return json.dumps({"cluster": conn.name, "cancelled": cancelled})


def run(executor, perms, *, cluster: Any, command: Any, timeout: Any = None) -> str:
    from delfin.agent import api_client as A
    from delfin.remote.errors import LoginRequired, RemoteError
    cmd = str(command or "").strip()
    if not cmd or len(cmd) > 4000:
        return _error("bad_argument", "command is one shell command line, at most 4000 characters")
    conn, refusal = _connect(perms, cluster, "cluster_run")
    if refusal:
        return refusal
    refused = _run_refusal(executor, perms, conn, cmd)
    if refused:
        _audit(perms, "cluster_run", "denied", command=cmd, reason=refused[1][:300],
               extra={"cluster": conn.name})
        return _error(*refused)
    egress = executor._scan_bash_egress(cmd, perms)
    if egress:
        A._record_security_event("egress", "cluster_run", f"{conn.name}: {egress}",
                                 blocked=False)
    try:
        limit = int(timeout) if timeout not in (None, "") else RUN_TIMEOUT
    except (TypeError, ValueError):
        return _error("bad_argument", "timeout is a number of seconds")
    limit = max(5, min(limit, RUN_TIMEOUT_MAX))
    args = {"cluster": conn.name, "command": cmd}
    preview = (f"$ {cmd}\n(on {conn.name}: {conn.cluster.target}, in the home directory; "
               f"stops after {limit} s)"
               + (f"\n[EGRESS] {egress}: data may leave the cluster" if egress else ""))
    refusal = _ask(executor, perms, "cluster_run", args, preview)
    if refusal:
        _audit(perms, "cluster_run", "denied", command=cmd, extra={"cluster": conn.name})
        return refusal
    try:
        got = conn.run(["bash", "-c", cmd], timeout=limit)
    except LoginRequired:
        return _login_needed(perms, conn.name, "cluster_run")
    except RemoteError as exc:
        _audit(perms, "cluster_run", "error", command=cmd, reason=str(exc)[:300],
               extra={"cluster": conn.name, "approved": True})
        return _error(exc.code, str(exc))
    _audit(perms, "cluster_run", "ok" if got.ok else "error", command=cmd,
           extra={"cluster": conn.name, "exit_code": got.exit_code, "approved": True})
    return _fenced(f"{conn.name}: {cmd[:60]}", {
        "cluster": conn.name, "exit_code": got.exit_code,
        "stdout": _tail(got.stdout), "stderr": _tail(got.stderr),
        **({"note": "stopped waiting after the timeout; the command may still run"}
           if got.exit_code == 124 else {})})


# --- shared ------------------------------------------------------------------------


def propose(executor, perms, *, name: Any = None, preset: Any = None, user: Any = None,
            host: Any = None, partition: Any = None) -> str:
    from delfin.agent import managed_policy
    from delfin.remote import setup
    if isinstance(name, str) and name.strip() and not managed_policy.cluster_allowed(name.strip()):
        return _error("not_allowed", f"the managed policy does not allow the cluster {name.strip()}")
    got = setup.propose(name=name, preset=preset, user=user, host=host, partition=partition,
                        session_id=getattr(perms, "task_session_id", "") or "",
                        workspace=getattr(perms, "workspace", None))
    if not got.get("ok"):
        return _error("not_proposed", str(got.get("error") or "the proposal was not recorded"))
    _audit(perms, "cluster_propose", "proposed", path=str(name or ""))
    return json.dumps({"ok": True, "proposed": str(name).strip(), "presets": got.get("presets", []),
                       "next": "the user checks the form in the DELFIN window (NEEDS YOU) and "
                               "sends it; nothing is added until then. Do not ask for a "
                               "password; end your turn or go on with other work."},
                      ensure_ascii=False)


def _error(code: str, message: str) -> str:
    return json.dumps({"error": message, "code": code}, ensure_ascii=False)


def _fenced(source: str, payload: dict) -> str:
    """The result as untrusted text: everything in it may come from the cluster."""
    from delfin.agent import untrusted
    return untrusted.wrap(f"cluster {source}"[:80], json.dumps(payload, ensure_ascii=False))


def _tail(text: Optional[str]) -> str:
    data = text or ""
    return data if len(data) <= OUTPUT_BYTES else "…" + data[-OUTPUT_BYTES:]


def _size(n: int) -> str:
    return f"{n / (1024 * 1024):.1f} MB" if n >= 1024 * 1024 else f"{max(1, n // 1024)} KB"


def _connect(perms, cluster: Any, tool: str):
    """The connection to ``cluster``, or the refusal the model reads."""
    from delfin.agent import managed_policy
    from delfin.remote import connection
    from delfin.remote.errors import RemoteError
    name = str(cluster or "").strip()
    if not name:
        return None, _error("bad_argument", "name a cluster (cluster_status lists them)")
    if not managed_policy.cluster_allowed(name):
        _audit(perms, tool, "denied", reason="managed policy", extra={"cluster": name})
        return None, _error("not_allowed", f"a managed policy does not allow the cluster {name}")
    try:
        return connection.connect(name), None
    except RemoteError as exc:
        return None, _error(exc.code, str(exc))


def _login_needed(perms, name: str, tool: str) -> str:
    """Nobody is logged in: a login item in NEEDS YOU, and what to tell the user."""
    from delfin.remote import attention as remote_attention
    remote_attention.ask_for_login(name, reason=f"{tool} needs a login to {name}",
                                   session_id=getattr(perms, "task_session_id", "") or "",
                                   workspace=getattr(perms, "workspace", None))
    return _error("login_required", (
        f"nobody is logged in to {name}. Tell the user: they log in from the DELFIN window "
        f"(/clusters, or the login item now in NEEDS YOU) or with `ctun login -t {name}` in a "
        "terminal. Never ask them for the password here, and do not try other ways in."))


def _ask(executor, perms, tool: str, args: dict, preview: str) -> Optional[str]:
    """None when a person said yes; else the refusal the model reads. Asked
    in every mode: these act on the user's account at a computing centre."""
    from delfin.agent import api_client as A
    if perms is None or perms.confirm_callback is None:
        return _error("needs_approval", (
            f"{tool} is approved by a person every time, and nobody can be asked in this "
            f"session (mode={getattr(perms, 'mode', '')}); it is refused. Tell the user what it "
            "would do."))
    signature = A._tool_action_signature(tool, args)
    if perms.action_denied_earlier("tool", signature):
        return _error("denied", f"the user refused this {tool} earlier in this session; do not "
                      "ask again, ask them what to do instead")
    try:
        ok = bool(A._ask_confirm(perms, tool, args, preview))
    except Exception as exc:
        return _error("approval_failed", f"the approval dialog failed: {exc}")
    if ok:
        return None
    if executor._confirm_timed_out(perms):
        return _error("approval_timed_out", (
            f"the approval of {tool} timed out -- the user is away; this is NOT a refusal. Do "
            "not try it again now; say what is waiting for them."))
    perms.record_denied_action("tool", signature)
    return _error("denied", f"the user refused {tool}. Do NOT retry it or reach the same "
                  "effect another way." + executor._said_why(perms))


def _audit(perms, tool: str, decision: str, *, command: str = "", path: str = "",
           reason: str = "", extra: Optional[dict] = None) -> None:
    try:
        from delfin.agent import audit_log
        audit_log.append(audit_log.make_record(
            tool=tool, decision=decision, mode=str(getattr(perms, "mode", "") or ""),
            command=command, path=path, reason=reason,
            session_id=str(getattr(perms, "task_session_id", "") or ""), extra=extra or {}))
    except Exception:
        pass  # the audit log never blocks what it records (audit_log.append)


def _place(executor, perms, conn, path: Any, tool: str, args: dict):
    """((real path, "file"|"dir"), None) for a path on the cluster that may be
    read -- asked unless it lies in the work dir -- or (None, refusal)."""
    from delfin.remote.errors import LoginRequired, RemoteError
    text = str(path or "").strip()
    if not text or "\n" in text or "\x00" in text:
        return None, _error("bad_argument", "path is one path on the cluster")
    try:
        got = conn.run(["bash", "-c", _RESOLVE, "delfin-resolve", conn.cluster.work_dir, text],
                       timeout=60)
    except LoginRequired:
        return None, _login_needed(perms, conn.name, tool)
    except RemoteError as exc:
        return None, _error(exc.code, str(exc))
    lines = (got.stdout or "").splitlines()
    if got.exit_code == 3 or len(lines) < 4:
        return None, _error("not_found", f"no {text} on {conn.name}")
    home, work, real, kind = lines[:4]
    if any(_inside(real, f"{home}/{name}") for name in _PRIVATE):
        return None, _error("refused", f"{real} is private to the user's account on {conn.name}; "
                            "it is never read or fetched")
    in_home = real[len(home):].lstrip("/") if _inside(real, home) else ""
    if (in_home and perms.matches_path_deny(in_home)) or perms.matches_path_deny(real):
        return None, _error("refused", f"{real} matches the secret deny list; it is never read "
                            "or fetched")
    if (_inside(real, _homes(home)) and not _inside(real, home)
            and not _inside(real, work)):
        return None, _error("refused", f"{real} is another user's; it is never read or fetched")
    if tool == "cluster_fetch" and kind == "dir" and _inside(home, real):
        # The home, or a folder holding it: a fetch would carry ~/.ssh and
        # the rest of _PRIVATE along.
        return None, _error("refused", f"{real} holds the account's private files (~/.ssh "
                            "and the like); fetch a folder below them, such as a job folder")
    if not _inside(real, work):
        verb = "read" if tool == "cluster_read" else "fetch"
        refusal = _ask(executor, perms, tool, args,
                       f"{verb} {real} on {conn.name} ({conn.cluster.target})\n"
                       f"(outside DELFIN's work dir {work})")
        _audit(perms, tool, "denied" if refusal else "ok", path=real,
               extra={"cluster": conn.name, "outside_work_dir": True})
        if refusal:
            return None, refusal
    return (real, kind), None


def _inside(path: str, root: str) -> bool:
    return path == root or path.startswith(root.rstrip("/") + "/")


def _homes(home: str) -> str:
    """Where the cluster keeps its users' homes: the nearest folder above the
    home named home (/home, /pfs/data5/home), else the home's parent."""
    parts = PurePosixPath(home).parts
    for i in range(len(parts) - 2, 0, -1):
        if parts[i] == "home":
            return str(PurePosixPath(*parts[:i + 1]))
    return str(PurePosixPath(home).parent)


def _under_work_dir(conn, dest: Any) -> tuple[str, str]:
    """``dest`` as a path under the cluster's work dir, or a reason it is not."""
    text = str(dest or "").strip().replace("\\", "/")
    parts = PurePosixPath(text).parts
    if (not _DEST.fullmatch(text) or text.startswith("/") or ".." in parts
            or not parts):
        return "", (f"dest is a folder under {conn.cluster.work_dir}, written relative to it "
                    f"with letters, digits, . _ - and / (not {text!r})")
    return f"{conn.cluster.work_dir.rstrip('/')}/{PurePosixPath(text)}", ""


def _file_refusal(executor, perms):
    """Why a workspace file may not leave this machine, or None."""
    from delfin.agent import api_client as A

    def refuse(path: Path) -> Optional[str]:
        if A._is_secret_path(perms, path):
            return "matches the secret deny list"
        if A._in_runtime_masked_dir(path):
            return "is DELFIN's own state"
        return None
    return refuse


def _upload_files(executor, perms, local: Path) -> tuple[list[Path], str]:
    refuse = _file_refusal(executor, perms)
    if local.is_file():
        candidates = [local]
    else:
        candidates = []
        for root, dirs, names in os.walk(local):
            for name in dirs + names:
                if (Path(root) / name).is_symlink():
                    return [], (f"{(Path(root) / name).relative_to(local)} is a link; uploads "
                                "copy no links")
            candidates += [Path(root) / n for n in names]
    if local.is_symlink():
        return [], f"{local.name} is a link; uploads copy no links"
    total = 0
    for path in candidates:
        why = refuse(path)
        if why:
            return [], f"{path.name} {why}; it does not leave this machine"
        total += path.stat().st_size
        if len(candidates) > UPLOAD_FILES or total > UPLOAD_BYTES:
            return [], (f"more than {UPLOAD_FILES} files or "
                        f"{UPLOAD_BYTES // (1024 * 1024)} MB in one upload")
    if not candidates:
        return [], f"{local.name} holds no file"
    return candidates, ""


def _run_refusal(executor, perms, conn, cmd: str) -> Optional[tuple[str, str]]:
    """(code, why) when cluster_run refuses ``cmd`` before anyone is asked."""
    from delfin.agent import api_client as A
    words = _command_words(cmd)
    if _runs_any(words, list(_SUBMIT) + _guard_commands(conn)):
        return "use_cluster_submit", ("cluster_run does not submit jobs (sbatch, srun, salloc, "
                                      "also inside bash -c): cluster_submit checks the budget "
                                      "and asks for the batch")
    if _runs_any(words, list(_CANCEL)):
        return "use_cluster_cancel", ("cluster_run does not cancel jobs: cluster_cancel "
                                      "cancels the user's own")
    hit = perms.matches_bash_deny(A._prose_blanked(cmd))
    if hit:
        return "denied", (f"the command matches a deny pattern ({hit}); it is refused on a "
                          "cluster as it is here")
    secret = executor._bash_denied_path(cmd, perms)
    if secret:
        return "denied", f"the command names a secret path ({secret}); refused"
    if A._prints_a_credential(cmd):
        return "denied", "the command would print a credential; refused"
    for site in ("bash", "cluster_run"):
        rule = A._permission_rule_hit(perms, site, {"command": cmd})
        if rule and rule[0] == "deny":
            return "denied", f"a deny rule of {site} covers this command; refused"
    return None


def _command_words(cmd: str) -> list[str]:
    """``cmd`` as argv for ctun's command reader: the line itself is a bash -c script."""
    return ["bash", "-c", cmd]


def _runs_any(words: list[str], names: list[str]) -> bool:
    """Whether the command line may run one of ``names`` anywhere a command can
    start (ctun's every-command reading, which also counts a command nested
    too deeply to follow)."""
    from cluster_tunnel import budget as ctun_budget
    return ctun_budget.is_guarded(words, names, scope="every-command")


def _guard_commands(conn) -> list[str]:
    """The cluster's own guard commands in ctun's config (aslurm, ...)."""
    try:
        from cluster_tunnel import api as ctun_api
        entry = ctun_api.load_config(conn.cluster.ctun_config).clusters[conn.name]
        return list(entry.budget.guard_commands) if entry.budget else []
    except Exception:
        return []
