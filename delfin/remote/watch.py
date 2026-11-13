"""DELFIN's cluster jobs followed to their end: the results brought back, the sessions told.

Input: the jobs DELFIN submitted (``delfin.remote.jobs.records``), each
cluster's Slurm through its tunnel -- one ctun probe per cluster and tick,
and ``delfin.remote.slurm`` for a job the probe no longer lists -- and the
user's login, which is never asked for here: a cluster without one gets a
login item in NEEDS YOU. Output: ``~/.delfin/remote_watch.json``, what the
watcher found of each job (``outcome`` and ``followed`` read it); an ended
job's key results in the workspace; a notice in NEEDS YOU; and, per
``tick``, the jobs that ended, for the agent service to wake their sessions.

Semantics:

* ``tick`` follows each job that has not ended. Per cluster with a live
  tunnel: one probe of the user's jobs (active, and ended since the oldest
  open job was submitted); a job the probe does not list is asked for by id.
  A cluster without a live tunnel is left alone, and NEEDS YOU is asked for
  a login at most every LOGIN_ASK_EVERY_S while jobs wait there; a cluster
  the managed policy keeps out is not asked. Either way its jobs are given
  up after GIVE_UP_DAYS. Ticks run one at a time across processes (a lock).
* A job has ended at a terminal Slurm state. Its key results then come back
  into ``<job folder>/<cluster>-<id>/`` in the workspace: the files ``KEEP``
  names, by importance, the folder's own before its subfolders', at most
  FILE_BYTES each and FETCH_BYTES in all. Wave functions, densities and
  scratch (``SKIP``) stay on the cluster, as does the job folder there; no
  link comes along (rsync -r copies none). The results folder is marked
  (``jobs.RESULTS_MARKER``): a resubmit of the job folder leaves it out.
* The outcome is judged on the fetched copy as on a local run
  (``api.calculation_status``: exit marker, run log, DELFIN_Data, ORCA's
  termination line) together with Slurm's state: ok when Slurm says
  COMPLETED and DELFIN's files say it succeeded, failed when either says it
  failed, unverified when Slurm says COMPLETED and the files cannot say. A
  job not ok gets the error signatures and job_fix's assessment of the copy.
* A fetch that fails (the tunnel went, a disk is full, the folder is gone)
  is tried again later, spread out over about a day (``_retry_after``),
  FETCH_TRIES times in all; then the job is reported with the reason and
  its results left on the cluster.
* Each ended job is in exactly one tick's ``ended``; ``outcome`` serves it
  to its session from then on (job_monitor's cluster watches). A job about
  which no answer came for GIVE_UP_DAYS is given up and reported so.
* A batch that waits for its install (``jobs.queue_after_install``) is
  followed with the cluster's jobs, ``batch-<id>`` its id: while the
  install runs it waits; once the install is complete it is submitted --
  exactly the approved files, the budget asked again (``jobs.submit``),
  nobody asked -- and its jobs are followed and watched for its session
  from then on. When the install failed or was stopped, or the files or the
  budget no longer allow the batch, nothing is submitted: the batch's
  outcome says why, a notice goes to NEEDS YOU, and it is in the tick's
  ``ended``, which wakes the session. A batch that went in is reported to
  its session at its next look and wakes nobody. A tick after the PC was
  off does the same with what it finds.
"""

from __future__ import annotations

import contextlib
import fnmatch
import json
import os
import re
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath
from typing import Any, Iterator, Optional

from . import jobs as _jobs
from .errors import LoginRequired, RemoteError

#: How long a cluster's probe answers for: a tick sooner asks no cluster
#: (the service's ticks pass their own interval's, at most this).
PROBE_EVERY_S = 55
FILE_BYTES = 100 * 1024 * 1024
FETCH_BYTES = 500 * 1024 * 1024
#: Files a job folder may list for a fetch (a runaway job writes thousands);
#: its own files are listed first, its subfolders' after them.
LIST_FILES = 5000
#: Fetches of an ended job's results, spread out (``_retry_after``): about a
#: day from the first to the last, so an outage of the cluster or the disk
#: does not lose them.
FETCH_TRIES = 8
GIVE_UP_DAYS = 30
#: How often NEEDS YOU is asked for a login a cluster's open jobs wait for.
LOGIN_ASK_EVERY_S = 12 * 3600

#: Slurm's states of a job that has ended (squeue and sacct names).
ENDED = frozenset({"COMPLETED", "FAILED", "TIMEOUT", "CANCELLED", "OUT_OF_MEMORY",
                   "NODE_FAIL", "PREEMPTED", "BOOT_FAIL", "DEADLINE", "REVOKED",
                   "SPECIAL_EXIT"})

#: What comes back of an ended job, the most important first.
KEEP = ("DELFIN_Data.json", "DELFIN*.txt", ".exit_code_*", "delfin_*.out", "delfin_*.err",
        "delfin_run.log", ".delfin_last_run.json", "slurm_reason.txt", "CONTROL.txt",
        "input.txt", "*.inp", "*.xyz", "*.out", "*.log", "*.err", "*.json", "*.txt",
        "*.hess", "*.docx", "*.png", "*.svg", "*.pdf", "*.csv")
#: Never, whatever KEEP says: wave functions, densities, scratch.
SKIP = ("*.gbw", "*.densities", "*.densitiesinfo", "*.tmp", "*.tmp.*", "*.cube", "*.ges",
        "*.uno", "*.qro", "*.unso", "*.bas*")

#: Runs on the cluster in the job's folder: its regular files, size and
#: path, NUL-separated (a name may hold anything but NUL).
_LIST = r'''cd -- "$1" 2>/dev/null || exit 3
{ find . -maxdepth 1 -type f -printf '%s\t%P\0'
  find . -mindepth 2 -maxdepth 4 -type f -printf '%s\t%P\0'; } 2>/dev/null | head -z -n "$2"
'''

_NAME = re.compile(r"[^A-Za-z0-9._-]+")
#: A path the fetch takes: plain names, which DELFIN's and ORCA's outputs have.
_PLAIN = re.compile(r"[\w.+%,=@-]+(?:/[\w.+%,=@-]+)*")


@dataclass
class Tick:
    """What one tick found."""

    ended: list[dict] = field(default_factory=list)   # ``outcome`` of each job that ended now
    following: int = 0           # jobs queued or running on a cluster with a live tunnel
    waiting_for_login: list[str] = field(default_factory=list)


class _FetchFailed(Exception):
    pass


# --- what sessions read ------------------------------------------------------------


def outcome(cluster: Any, job_id: Any) -> Optional[dict]:
    """What the watcher made of an ended job, once its results came back (or
    could not); None while it is queued, running or being fetched."""
    entry = (_load().get("jobs") or {}).get(f"{cluster}:{job_id}") or {}
    return dict(entry["outcome"]) if entry.get("outcome") else None


def seen_state(cluster: Any, job_id: Any) -> str:
    """The job's Slurm state as the watcher last saw it; "" before then."""
    entry = (_load().get("jobs") or {}).get(f"{cluster}:{job_id}") or {}
    return str(entry.get("state") or "")


def followed(cluster: Any, job_id: Any) -> bool:
    """Whether the job is still to be reported: until the watcher has its
    outcome (or gave up on it, after GIVE_UP_DAYS), also before its first
    tick saw the job."""
    return outcome(cluster, job_id) is None


# --- the tick ----------------------------------------------------------------------


def tick(*, now: Optional[float] = None, probe_every_s: Optional[float] = None) -> Tick:
    """Follow every open job one step; the jobs that ended in this step."""
    now = time.time() if now is None else now
    probe_every_s = PROBE_EVERY_S if probe_every_s is None else probe_every_s
    out = Tick()
    with _locked() as held:
        if not held:
            return out             # another process ticks right now
        state = _load()
        known = state.setdefault("jobs", {})
        open_by_cluster: dict[str, list[dict]] = {}
        for record in _jobs.records():
            key = _key(record)
            if key and not (known.get(key) or {}).get("outcome"):
                open_by_cluster.setdefault(str(record["cluster"]), []).append(record)
        batches_by_cluster: dict[str, list[dict]] = {}
        for entry in _jobs.pending():
            batches_by_cluster.setdefault(str(entry["cluster"]), []).append(entry)
        probes = state.setdefault("probed", {})
        for name in sorted(set(open_by_cluster) | set(batches_by_cluster)):
            records = open_by_cluster.get(name, [])
            batches = batches_by_cluster.get(name, [])
            if now - float(probes.get(name) or 0) < probe_every_s:
                out.following += len(batches) + sum(
                    1 for r in records if (known.get(_key(r)) or {}).get("state") not in ENDED)
                continue
            if _follow(name, records, state, now, out, batches=batches):
                probes[name] = now
        _save(state)
    return out


def _follow(name: str, records: list[dict], state: dict, now: float, out: Tick, *,
            batches: Optional[list[dict]] = None) -> bool:
    """One step for a cluster's open jobs and the batches that wait for its
    install; whether its Slurm was asked. A cluster that cannot be asked
    keeps its jobs and batches until they are old enough to give up."""
    from delfin.agent import managed_policy
    from . import attention as remote_attention
    from . import connection, slurm
    known = state.setdefault("jobs", {})
    batches = list(batches or [])

    def give_up(why: str) -> tuple[list[dict], list[dict]]:
        return ([r for r in records if not _give_up_if_old(r, known, now, out, why)],
                [b for b in batches if not _give_up_batch_if_old(b, known, now, out, why)])

    if not managed_policy.cluster_allowed(name):
        give_up(f"a managed policy keeps {name} out")
        return False
    try:
        conn = connection.connect(name)
        live = conn.is_live()
    except RemoteError as exc:
        # Not usable (its settings changed): asked again next tick.
        give_up(f"{name} is not usable ({exc})")
        return False
    if not live:
        waiting, waiting_batches = give_up(f"nobody logged in to {name}")
        asked = state.setdefault("login_asked", {})
        if (waiting or waiting_batches) and now - float(asked.get(name) or 0) >= LOGIN_ASK_EVERY_S:
            # Once in a while, not every tick: a login item the user put
            # aside does not come back a minute later.
            asked[name] = now
            said = []
            if waiting:
                said.append(f"{len(waiting)} DELFIN job(s) on {name} are followed; their results "
                            "come back and their sessions are told once you log in")
            if waiting_batches:
                said.append(f"{len(waiting_batches)} approved batch(es) wait for DELFIN's install "
                            f"on {name}; they are submitted once you log in and it is done")
            remote_attention.ask_for_login(name, reason="; ".join(said))
        if waiting or waiting_batches:
            out.waiting_for_login.append(name)
        return False
    for entry in batches:
        if not _give_up_batch_if_old(entry, known, now, out, "its install never completed"):
            try:
                _resolve_batch(conn, entry, known, now, out)
            except LoginRequired:
                out.waiting_for_login.append(name)
                return False
    if not records:
        return True
    oldest = min(float(r.get("submitted_at") or now) for r in records)
    since = int(min(GIVE_UP_DAYS * 1440, (now - oldest) / 60 + 60))
    try:
        listing = {str(j.jobid): j for j in conn.jobs(since_minutes=since)}
    except LoginRequired:
        out.waiting_for_login.append(name)
        return False
    except RemoteError:
        return False
    for record in records:
        key = _key(record)
        seen = known.setdefault(key, {"cluster": name, "job_id": str(record["job_id"])})
        job = listing.get(str(record["job_id"]))
        said = _slurm_state(job.state) if job is not None else ""
        if not said and not seen.get("ended_at"):
            try:
                said = _slurm_state(slurm.job_fields(conn, str(record["job_id"]))
                                    .get("JobState", ""))
            except LoginRequired:
                out.waiting_for_login.append(name)
                return False
            except RemoteError:
                said = ""
        if said:
            seen["state"], seen["seen_at"] = said, now
        if said in ENDED or seen.get("ended_at"):
            seen.setdefault("ended_at", now)
            if float(seen.get("next_fetch_at") or 0) > now:
                continue           # a failed fetch waits for its next try
            ended = _resolve(conn, record, seen, now)
            if ended is not None:
                out.ended.append(ended)
        elif said:
            out.following += 1
        else:
            _give_up_if_old(record, known, now, out, f"{name} has said nothing of the job")
    return True


def _give_up_if_old(record: dict, known: dict, now: float, out: Tick, why: str) -> bool:
    """Give a job up -- reported as such -- once it is GIVE_UP_DAYS old;
    whether it was."""
    if now - float(record.get("submitted_at") or now) <= GIVE_UP_DAYS * 86400:
        return False
    seen = known.setdefault(_key(record), {"cluster": str(record["cluster"]),
                                           "job_id": str(record["job_id"])})
    seen["outcome"] = _outcome(record, seen, state="UNKNOWN", ok=False, judged={
        "verdict": "failed", "status": "unknown",
        "outcome": f"given up after {GIVE_UP_DAYS} days: {why}"})
    _notify(record, seen["outcome"])
    out.ended.append(seen["outcome"])
    return True


def _resolve(conn, record: dict, seen: dict, now: float) -> Optional[dict]:
    """Fetch an ended job's results and judge them; its outcome, or None
    while a failed fetch is to be tried again."""
    seen["fetch_tries"] = int(seen.get("fetch_tries") or 0) + 1
    try:
        fetched = _fetch(conn, record)
        seen.pop("fetch_error", None)
    except LoginRequired:
        seen["fetch_tries"] -= 1    # no attempt: tried again once logged in
        return None
    except (_FetchFailed, RemoteError, OSError) as exc:
        seen["fetch_error"] = str(exc)[:300]
        if seen["fetch_tries"] < FETCH_TRIES:
            seen["next_fetch_at"] = now + _retry_after(seen["fetch_tries"])
            return None
        fetched = None
    judged = _judge(record, seen, fetched)
    seen["outcome"] = _outcome(record, seen, state=str(seen.get("state") or "?"),
                               ok=judged.pop("ok"), judged=judged, fetched=fetched)
    _notify(record, seen["outcome"])
    return seen["outcome"]


def _retry_after(tries: int) -> float:
    """Seconds to the next fetch after ``tries`` failed ones: 1, 3, 9 ...
    minutes, at most twelve hours."""
    return float(min(12 * 3600, 60 * 3 ** max(0, tries - 1)))


# --- a batch that waits for its install ------------------------------------------


def _batch_key(entry: dict) -> str:
    return f"{entry['cluster']}:batch-{entry['id']}"


def _resolve_batch(conn, entry: dict, known: dict, now: float, out: Tick) -> None:
    """One step for a batch that waits for its install: wait, submit it, or
    say why it was not. Raises LoginRequired when the tunnel went."""
    from . import install as _install
    try:
        current = _install.state(conn, entry["commit"])
    except LoginRequired:
        raise
    except RemoteError:
        out.following += 1
        return                              # asked again next tick
    if current == "running":
        out.following += 1
        return
    done: list[dict] = []
    if current == "complete":
        batch = _jobs.batch_from_dict(entry["batch"])
        batch.install = "complete"
        try:
            done = _jobs.submit(conn, batch, session_id=str(entry.get("session_id") or ""),
                                workspace=entry.get("workspace") or None)
            why = ""
        except LoginRequired:
            raise
        except _jobs.PartialSubmit as exc:
            done, why = list(exc.done), str(exc)
        except RemoteError as exc:
            why = str(exc)
    elif current == "stopped":
        why = f"the install of DELFIN {entry['commit'][:12]} does not go on: {_install.STOPPED_HINT}"
    else:
        try:
            tail = _install.log_tail(conn, entry["commit"]).strip().splitlines()[-3:]
        except RemoteError:
            tail = []
        why = (f"the install of DELFIN {entry['commit'][:12]} on {conn.name} {current}"
               + (": " + " | ".join(t.strip() for t in tail) if tail else "")
               + "; see /clusters (space)")
    ids = [str(d["job_id"]) for d in done]
    for record in done:
        _watch_for_session(record)
    _jobs.settle_pending(entry["id"], state="submitted" if done and not why else "not_submitted",
                         job_ids=ids, reason=why, settled_at=now)
    names = ", ".join(entry.get("jobs") or []) or "the batch"
    if done and not why:
        text = (f"submitted once DELFIN {entry['commit'][:12]} was installed: job(s) "
                f"{', '.join(ids)} ({names}); each is followed, and told when it ends")
        verdict, state = "ok", "SUBMITTED"
    elif done:
        text = f"partly submitted: job(s) {', '.join(ids)} went in and are followed -- {why}"
        verdict, state = "failed", "PARTLY_SUBMITTED"
    else:
        text = f"not submitted ({names}): {why}"
        verdict, state = "failed", "NOT_SUBMITTED"
    record = {"cluster": entry["cluster"], "job_id": f"batch-{entry['id']}",
              "name": names, "kind": "batch", "remote_dir": "",
              "workspace": entry.get("workspace") or "",
              "session_id": entry.get("session_id") or ""}
    seen = known.setdefault(_batch_key(entry), {"cluster": entry["cluster"],
                                                "job_id": record["job_id"]})
    seen["state"], seen["ended_at"] = state, now
    seen["outcome"] = _outcome(record, seen, state=state, ok=verdict == "ok", judged={
        "verdict": verdict, "status": "submitted" if verdict == "ok" else "not_submitted",
        "outcome": text})
    seen["outcome"]["job_ids"] = ids
    if verdict != "ok":
        # Nothing (or not all) went in: the person hears of it, and the
        # session is woken to say so. A batch that went in wakes nobody; its
        # session reads it at its next look, and its jobs when they end.
        _notify(record, seen["outcome"])
        out.ended.append(seen["outcome"])


def _give_up_batch_if_old(entry: dict, known: dict, now: float, out: Tick, why: str) -> bool:
    """A batch that waited GIVE_UP_DAYS: given up, reported as not submitted."""
    if now - float(entry.get("queued_at") or now) <= GIVE_UP_DAYS * 86400:
        return False
    text = f"given up after {GIVE_UP_DAYS} days: {why}"
    _jobs.settle_pending(entry["id"], state="not_submitted", job_ids=[], reason=text,
                         settled_at=now)
    record = {"cluster": entry["cluster"], "job_id": f"batch-{entry['id']}",
              "name": ", ".join(entry.get("jobs") or []) or "the batch", "kind": "batch",
              "remote_dir": "", "workspace": entry.get("workspace") or "",
              "session_id": entry.get("session_id") or ""}
    seen = known.setdefault(_batch_key(entry), {"cluster": entry["cluster"],
                                                "job_id": record["job_id"]})
    seen["outcome"] = _outcome(record, seen, state="NOT_SUBMITTED", ok=False, judged={
        "verdict": "failed", "status": "not_submitted", "outcome": f"not submitted: {text}"})
    _notify(record, seen["outcome"])
    out.ended.append(seen["outcome"])
    return True


def _watch_for_session(record: dict) -> None:
    """A job the watcher submitted, watched for the session that approved its
    batch (as cluster_submit watches its own). Never raises."""
    try:
        from delfin.agent import job_monitor
        workspace = str(record.get("workspace") or "")
        if not workspace:
            return
        what = "ORCA run" if record.get("kind") == "orca" else "DELFIN workflow"
        job_monitor.register_agent_job(
            workspace, f"{record['cluster']}:{record['job_id']}",
            description=f"{record['name']}: {what} on {record['cluster']}",
            extra={"kind": "cluster", "cluster": str(record["cluster"]),
                   "remote_id": str(record["job_id"]),
                   "session_id": str(record.get("session_id") or "")})
    except Exception:
        pass


# --- the fetch ---------------------------------------------------------------------


def results_dir(record: dict) -> Path:
    """Where an ended job's results land: beside its inputs, in a folder of
    its own (``<job folder>/<cluster>-<id>``); under ``cluster/<cluster>/``
    of the workspace when its folder is gone."""
    workspace = Path(str(record.get("workspace") or "")).expanduser()
    folder = Path(str(record.get("folder") or "")).expanduser()
    name = f"{_NAME.sub('-', str(record['cluster']))}-{_NAME.sub('-', str(record['job_id']))}"
    if folder.is_dir() and _within(folder, workspace):
        return folder / name
    job = _NAME.sub("-", str(record.get("name") or "job")).strip("-.") or "job"
    return workspace / "cluster" / _NAME.sub("-", str(record["cluster"])) / f"{job}-{name}"


def _fetch(conn, record: dict) -> dict:
    workspace = Path(str(record.get("workspace") or "")).expanduser()
    if not workspace.is_absolute() or not workspace.is_dir():
        raise _FetchFailed(f"the workspace {workspace} is gone")
    if workspace.resolve() in (Path.home().resolve(), Path("/")):
        # No session works there (KitToolPermissions refuses it), so no
        # record DELFIN wrote names it.
        raise _FetchFailed(f"{workspace} is no workspace results are written into")
    dest = results_dir(record)
    if dest.is_symlink() or (dest.exists() and not dest.is_dir()):
        raise _FetchFailed(f"{dest} is not a folder DELFIN made")
    remote = _on_cluster(str(record.get("remote_dir") or ""))
    listed = conn.run(["bash", "-c", _LIST, "delfin-list", remote, str(LIST_FILES)],
                      timeout=120)
    if listed.exit_code == 3:
        raise _FetchFailed(f"the job folder {record.get('remote_dir')} is gone on {conn.name}")
    if not listed.ok:
        raise _FetchFailed(f"the job folder could not be listed: {(listed.stderr or '')[:200]}")
    chosen, left = _choose(listed.stdout or "")
    dest.mkdir(parents=True, exist_ok=True)
    # What marks it as fetched results: cluster_submit never uploads it.
    (dest / _jobs.RESULTS_MARKER).write_text(json.dumps({
        "cluster": str(record["cluster"]), "job_id": str(record["job_id"]),
        "remote_dir": str(record.get("remote_dir") or ""), "fetched_at": time.time()}),
        encoding="utf-8")
    if chosen:
        fd, names = tempfile.mkstemp(prefix="delfin-fetch-", suffix=".list")
        try:
            with os.fdopen(fd, "wb") as fh:
                fh.write(b"".join(p.encode("utf-8", "surrogateescape") + b"\0"
                                  for p, _ in chosen))
            got = conn.download(f"{remote.rstrip('/')}/", f"{dest}/",
                                extra=["--times", "--from0", f"--files-from={names}"])
        finally:
            os.unlink(names)
        # 23/24: some files could not be read or vanished; the rest came.
        if not got.ok and got.exit_code not in (23, 24):
            raise _FetchFailed(f"the copy failed: {(got.stderr or '').strip()[:200]}")
    came = [(p, s) for p, s in chosen if (dest / p).is_file()]
    missing = [p for p, s in chosen if (p, s) not in came]
    return {"dest": dest, "files": len(came), "bytes": sum(s for _, s in came),
            "left": left + missing}


def _choose(listing: str) -> tuple[list[tuple[str, int]], list[str]]:
    """The files of a job folder's listing to bring back, by importance, and
    the ones KEEP names that stay for their size."""
    entries = []
    for item in listing.split("\0"):
        size, _, path = item.partition("\t")
        parts = PurePosixPath(path).parts
        if not size.isdigit() or not _PLAIN.fullmatch(path) or ".." in parts:
            continue
        base = parts[-1]
        if any(fnmatch.fnmatchcase(base, g) for g in SKIP):
            continue
        rank = next((i for i, g in enumerate(KEEP) if fnmatch.fnmatchcase(base, g)), None)
        if rank is not None:
            entries.append((rank, len(parts), int(size), path))
    chosen, left, total = [], [], 0
    for rank, depth, size, path in sorted(entries):
        if size > FILE_BYTES or total + size > FETCH_BYTES:
            left.append(path)
            continue
        chosen.append((path, size))
        total += size
    return chosen, left


# --- the judgement -----------------------------------------------------------------


def _judge(record: dict, seen: dict, fetched: Optional[dict]) -> dict:
    """ok, status and outcome from Slurm's state and the fetched copy; the
    signatures and job_fix's assessment of a failed one."""
    state = str(seen.get("state") or "")
    status, phrase, exit_code = "unknown", f"Slurm: {state or 'no state'}", None
    if fetched is not None:
        try:
            from delfin.api import calculation_status
            st = calculation_status(str(fetched["dest"]))
            status, phrase = st.state, st.outcome
            for ev in st.evidence:
                code = re.match(r"exit code (-?\d+)", str(ev.get("says") or ""))
                if code and str(ev.get("source", "")).startswith(".exit_code_"):
                    exit_code = int(code.group(1))
        except Exception:
            pass
    # Slurm's COMPLETED says the job script ended with 0; DELFIN's own files
    # say whether the run succeeded. Without them it is not called a success.
    if state != "COMPLETED" or status == "failed":
        verdict = "failed"
    elif status in ("succeeded", "finished"):
        verdict = "ok"
    else:
        verdict = "unverified"
        phrase = (f"unverified: Slurm says COMPLETED, DELFIN's own files do not say it "
                  f"succeeded ({phrase})" if fetched is not None else
                  "unverified: Slurm says COMPLETED, and its results could not be checked")
    ok = verdict == "ok"
    judged = {"ok": ok, "verdict": verdict, "status": status, "outcome": phrase,
              "exit_code": exit_code, "signatures": [], "assessment": ""}
    if not ok and fetched is not None:
        try:
            from delfin.agent import job_fix, job_monitor
            sigs = job_monitor.scan_error_signatures(str(fetched["dest"]))
            fix = job_fix.assess(f"{record['cluster']}:{record['job_id']}",
                                 str(fetched["dest"]), sigs)
            judged["signatures"] = list(sigs)
            judged["assessment"] = f"{fix.summary} {fix.recommendation}".strip()[:400]
        except Exception:
            pass
    return judged


def _outcome(record: dict, seen: dict, *, state: str, ok: bool, judged: dict,
             fetched: Optional[dict] = None) -> dict:
    workspace = Path(str(record.get("workspace") or "")).expanduser()
    results = ""
    if fetched is not None:
        try:
            results = os.path.relpath(fetched["dest"], workspace)
        except ValueError:
            results = str(fetched["dest"])
    return {
        "key": f"{record['cluster']}:{record['job_id']}", "cluster": str(record["cluster"]),
        "job_id": str(record["job_id"]), "name": str(record.get("name") or ""),
        "kind": str(record.get("kind") or ""), "state": state, "ok": bool(ok),
        "verdict": judged.get("verdict") or ("ok" if ok else "failed"),
        "status": judged.get("status", "unknown"), "outcome": judged.get("outcome", ""),
        "exit_code": judged.get("exit_code"), "signatures": judged.get("signatures") or [],
        "assessment": judged.get("assessment", ""), "results": results,
        "files": int((fetched or {}).get("files") or 0),
        "left_on_cluster": len((fetched or {}).get("left") or []),
        "fetch_error": str(seen.get("fetch_error") or "") if fetched is None else "",
        "remote_dir": str(record.get("remote_dir") or ""),
        "session_id": str(record.get("session_id") or ""), "workspace": str(workspace),
        "ended_at": float(seen.get("ended_at") or time.time())}


def summary(ended: dict) -> str:
    """One line of what became of a job, as a session or a notice says it."""
    parts = [f"{ended.get('outcome') or ended.get('state')}"]
    if ended.get("results"):
        parts.append(f"results in {ended['results']} ({ended.get('files', 0)} files"
                     + (f"; {ended['left_on_cluster']} too large, left on the cluster"
                        if ended.get("left_on_cluster") else "") + ")")
    elif ended.get("fetch_error"):
        parts.append(f"its results could not be brought back ({ended['fetch_error']}); "
                     f"they stay on {ended.get('cluster')} in {ended.get('remote_dir')}")
    if ended.get("assessment"):
        parts.append(f"job_fix: {ended['assessment']}")
    return " -- ".join(parts)


def _notify(record: dict, ended: dict) -> None:
    """The job's notice in NEEDS YOU (run_finished or run_failed). Never raises."""
    try:
        from delfin.agent.attention import emit_attention
        verdict = ended.get("verdict") or ("ok" if ended["ok"] else "failed")
        said = {"ok": f"finished ({ended['state']})",
                "unverified": f"finished ({ended['state']}, unverified)"}.get(
                    verdict, f"failed ({ended['state']})")
        title = f"Job {ended['key']} {said}"
        if ended.get("kind") == "batch":
            title = (f"A batch for {ended['cluster']} was partly submitted"
                     if ended["state"] == "PARTLY_SUBMITTED" else
                     f"A batch for {ended['cluster']} was not submitted")
        emit_attention("run_failed" if verdict == "failed" else "run_finished",
                       session_id=ended["session_id"],
                       title=title,
                       detail=f"{ended['name'] or 'job'}: {summary(ended)}",
                       workspace=ended["workspace"] or None)
    except Exception:
        pass


# --- small helpers -----------------------------------------------------------------


def _key(record: dict) -> str:
    cluster, job_id = str(record.get("cluster") or ""), str(record.get("job_id") or "")
    return f"{cluster}:{job_id}" if cluster and job_id else ""


def _slurm_state(text: Any) -> str:
    """The state's word (``CANCELLED by 1234`` is CANCELLED), or ""."""
    word = str(text or "").strip().split(" ")[0].rstrip("+").upper()
    return word if re.fullmatch(r"[A-Z_]{2,32}", word) else ""


def _on_cluster(path: str) -> str:
    """A cluster path as rsync and a shell started in the home take it."""
    if path in ("~", "~/"):
        return "."
    return path[2:] if path.startswith("~/") else path


def _within(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
        return True
    except (OSError, ValueError):
        return False


def _state_path() -> Path:
    # Resolved per call, so the suite's state redirect reaches it (state_paths).
    return Path.home() / ".delfin" / "remote_watch.json"


def _lock_path() -> Path:
    return Path.home() / ".delfin" / "remote_watch.lock"


def _load() -> dict:
    try:
        data = json.loads(_state_path().read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {}
    except (OSError, ValueError):
        return {}


def _save(state: dict) -> None:
    path = _state_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".remote_watch-", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            json.dump(state, fh, indent=1, sort_keys=True)
        os.replace(tmp, path)
    except BaseException:
        with contextlib.suppress(OSError):
            os.unlink(tmp)
        raise


@contextlib.contextmanager
def _locked() -> Iterator[bool]:
    """Held by one tick at a time across processes; False when another holds it."""
    import fcntl
    path = _lock_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    fd = os.open(path, os.O_RDWR | os.O_CREAT, 0o600)
    try:
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            yield False
            return
        try:
            yield True
        finally:
            fcntl.flock(fd, fcntl.LOCK_UN)
    finally:
        os.close(fd)
