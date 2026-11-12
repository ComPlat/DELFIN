"""The clusters and the jobs a window shows, and the login it hands to ctun.

Input: DELFIN's clusters (``delfin.remote.clusters``) and ctun's state; the
agent's background work per session workspace (``background_view``); the
SLURM jobs DELFIN's job monitor watches (``job_monitor.load_watched``).
Output: plain dicts for the service's ``clusters.list`` and ``jobs.list``
results; ``login`` and ``logout`` for ``cluster.login`` / ``cluster.logout``.
Runs in the supervisor, which calls it through ``asyncio.to_thread``.

Semantics:

* A cluster row: name, target, state -- ``refused`` (with its problems),
  ``logged_out`` or ``logged_in`` -- what a login asks for, and while logged
  in the session budget and the number of the user's active jobs. The login
  state is the local control socket; the budget and the job list are asked
  over the tunnel at most once a minute per cluster, since a window polls
  every few seconds.
* While logged in, a row also says whether ORCA is where DELFIN's installer
  looks for it (``delfin.remote.setup.orca_status``, asked at most every ten
  minutes): ORCA is the user's own, so a row without it says what to do.
* ``login``: the password and OTP go to ``delfin.remote.connection``, which
  hands them to ctun, and nowhere else. They are in no return value, log
  line or error text: a failure is reported by its code and a message that
  names neither (ctun's LoginError never does). The supervisor takes this
  request only from a window, never from an agent's process (``ACTING``),
  and never forwards it to a worker.
* ``jobs``: three kinds of group -- the agent's background shells, watches,
  sub-agents and wake-ups, one group per session; the jobs DELFIN's job
  monitor watches; the user's active jobs on each logged-in cluster, all of
  them (also those submitted by hand), asked at most twice a minute per
  cluster, after the batches that wait there for an install (also on a
  cluster nobody is logged in to: they wait for the login). A source that
  cannot be read is named in ``errors``, never shown as empty.
* Everything shown comes from places other processes write (a cluster's job
  names, a watch file in a workspace), so every string is cut to one line
  of printable characters: no escape sequence reaches the terminal.
"""

from __future__ import annotations

import re
import threading
import time
from typing import Any, Callable, Iterable, Optional

#: Seconds between two probes over a cluster's tunnel for its row (budget,
#: job count) and for /jobs (the job list).
ROW_PROBE_S = 60.0
JOBS_PROBE_S = 30.0
#: Seconds between two looks for ORCA on a logged-in cluster.
ORCA_PROBE_S = 600.0

_CONTROL = re.compile("[\x00-\x1f\x7f-\x9f  ‪-‮⁦-⁩]")
# Terminal escape sequences, taken out whole (a coloured log line reads as
# its words): CSI (ESC [ ... final byte), OSC (ESC ] ... BEL or ESC \),
# and any other ESC pair.
_ESCAPE = re.compile("\x1b\\[[0-?]*[ -/]*[@-~]|\x1b\\][^\x07\x1b]*(?:\x07|\x1b\\\\)?|\x1b[@-_]?")
_lock = threading.Lock()
_cache: dict[tuple[str, str], tuple[float, Any]] = {}


def clean(value: Any, limit: int = 160) -> str:
    """``value`` as one line of printable text, at most ``limit`` characters."""
    text = _CONTROL.sub(" ", _ESCAPE.sub("", str(value if value is not None else "")))
    text = " ".join(text.split())
    return text if len(text) <= limit else text[: limit - 1] + "…"


def clean_line(value: Any, limit: int = 300) -> str:
    """``clean`` for a line shown as it is (an output tail, an indented
    list): the same characters taken out, the spacing kept."""
    text = _ESCAPE.sub("", str(value if value is not None else "")).expandtabs(4)
    text = _CONTROL.sub(" ", text).rstrip()
    return text if len(text) <= limit else text[: limit - 1] + "…"


def _cached(key: tuple[str, str], max_age: float, probe: Callable[[], Any],
            now: Optional[float] = None) -> Any:
    now = time.monotonic() if now is None else now
    with _lock:
        hit = _cache.get(key)
        if hit is not None and now - hit[0] < max_age:
            return hit[1]
    value = probe()
    with _lock:
        _cache[key] = (now, value)
    return value


def forget(name: Optional[str] = None) -> None:
    """Drop the cached probes (of one cluster), after a login or logout."""
    with _lock:
        for key in [k for k in _cache if name is None or k[1] == name]:
            del _cache[key]


# --- clusters ------------------------------------------------------------------------


def clusters(*, settings: Optional[dict] = None) -> dict:
    """``{"clusters": [row, ...]}`` for every cluster in remote.clusters."""
    from delfin.remote import clusters as registry
    from delfin.remote import connection
    from delfin.remote.errors import RemoteError

    rows = []
    for name, cluster in registry.load(settings=settings).items():
        row = {"name": clean(name, 64), "target": clean(cluster.target, 120),
               "state": "refused", "problems": [clean(p, 300) for p in cluster.problems],
               "requiresPassword": bool(cluster.requires_password),
               "requiresOtp": bool(cluster.requires_otp),
               "workDir": clean(cluster.work_dir, 200),
               "defaultPartition": clean(cluster.default_partition or "", 64) or None,
               "budget": None, "jobs": None, "orca": None, "error": None}
        rows.append(row)
        if not cluster.usable:
            continue
        try:
            conn = connection.connect(name, settings=settings)
            live = conn.is_live()
        except RemoteError as exc:
            row["error"] = clean(f"{exc.code}: {exc}", 300)
            continue
        row["state"] = "logged_in" if live else "logged_out"
        if not live:
            continue
        try:
            probe = _cached(("row", name), ROW_PROBE_S, lambda: _row_probe(conn))
        except RemoteError as exc:
            row["error"] = clean(f"{exc.code}: {exc}", 300)
            continue
        row["budget"], row["jobs"] = probe
        row["orca"] = _orca(name, conn)
    return {"clusters": rows}


def _orca(name: str, conn) -> Optional[dict]:
    """Where ORCA is on a logged-in cluster (cached); None when not asked."""
    from delfin.remote import setup
    try:
        return _cached(("orca", name), ORCA_PROBE_S, lambda: setup.orca_status(conn))
    except Exception:  # noqa: BLE001 - not known now: the row says nothing of it
        return None


def _row_probe(conn) -> tuple[Optional[dict], int]:
    decision = conn.budget()
    budget = None
    if decision.limit is not None:
        budget = {"used": decision.used, "limit": decision.limit,
                  "unit": clean(decision.unit, 16)}
    active = [j for j in conn.jobs() if getattr(j, "source", "") == "active"]
    return budget, len(active)


def login(name: str, password: str, otp: str, *, settings: Optional[dict] = None,
          timeout: int = 120) -> dict:
    """``{"ok", "opened", "code", "error"}``; never the password or the OTP."""
    from delfin.remote import connection
    from delfin.remote.errors import RemoteError

    try:
        conn = connection.connect(name, settings=settings)
        opened = conn.login(password, otp, timeout=timeout)
    except RemoteError as exc:
        return {"ok": False, "opened": False, "code": exc.code,
                "error": _scrub(str(exc), password, otp)}
    except Exception as exc:  # noqa: BLE001 - the window gets a word, not a traceback
        return {"ok": False, "opened": False, "code": "internal_error",
                "error": _scrub(type(exc).__name__, password, otp)}
    finally:
        forget(name)
    from delfin.remote import attention as remote_attention
    remote_attention.logged_in(name)
    return {"ok": True, "opened": opened, "code": None, "error": None}


def logout(name: str, *, settings: Optional[dict] = None) -> dict:
    from delfin.remote import connection
    from delfin.remote.errors import RemoteError

    try:
        closed = connection.connect(name, settings=settings).logout()
    except RemoteError as exc:
        return {"ok": False, "closed": False, "code": exc.code, "error": clean(str(exc), 300)}
    finally:
        forget(name)
    return {"ok": True, "closed": closed, "code": None, "error": None}


def _scrub(text: str, *secrets: str) -> str:
    """``text`` cleaned, with any secret in it masked: ctun names none, and
    this keeps it so whatever a later release says."""
    out = clean(text, 300)
    for secret in secrets:
        if secret:
            out = out.replace(secret, "***")
    return out


# --- jobs ----------------------------------------------------------------------------


def jobs(sessions: Iterable[dict], *, settings: Optional[dict] = None) -> dict:
    """``{"groups": [...], "errors": [...]}`` for the /jobs screen.

    ``sessions``: the service's sessions, each with ``sessionId``,
    ``workspace`` and ``title``. A group: ``source`` (agent | monitor |
    cluster), ``title``, ``rows`` (each with ``id``, ``kind``, ``label``,
    ``state``, ``detail``) and, for a cluster, ``cluster``."""
    groups: list[dict] = []
    errors: list[str] = []
    try:
        from delfin.agent import background_view
        for s in sessions:
            workspace = str(s.get("workspace") or "")
            if not workspace:
                continue
            view = background_view.collect(workspace, session_id=str(s.get("sessionId") or ""))
            listed = background_view.rows(view)
            errors += [clean(f"{s.get('title') or workspace}: {r['label']} ({r['detail']})", 200)
                       for r in listed if r.get("group") == "errors"]
            rows = [{"id": clean(r["id"], 64), "kind": clean(r["kind"], 40),
                     "label": clean(r["label"], 120), "state": "",
                     "detail": clean(r["detail"], 160), "group": clean(r.get("group"), 16)}
                    for r in listed if r.get("group") != "errors"]
            if rows:
                groups.append({"source": "agent", "session": clean(s.get("sessionId"), 64),
                               "title": clean(s.get("title") or workspace, 80),
                               "rows": rows})
    except Exception as exc:  # noqa: BLE001
        errors.append(clean(f"the agents' background work could not be read: {exc}", 200))
    try:
        from delfin.agent import job_monitor
        watched = job_monitor.load_watched().get("jobs") or {}
        rows = [{"id": clean(jid, 64), "kind": "SLURM",
                 "label": clean((rec or {}).get("folder") or jid, 120),
                 "state": clean((rec or {}).get("last_state") or "", 24),
                 "detail": "watched by DELFIN's job monitor"}
                for jid, rec in watched.items()]
        if rows:
            groups.append({"source": "monitor", "title": "DELFIN's job monitor", "rows": rows})
    except Exception as exc:  # noqa: BLE001
        errors.append(clean(f"the job monitor's list could not be read: {exc}", 200))
    groups += _cluster_groups(settings, errors)
    return {"groups": groups, "errors": errors}


def _cluster_groups(settings: Optional[dict], errors: list[str]) -> list[dict]:
    from delfin.remote import clusters as registry
    from delfin.remote import connection
    from delfin.remote.errors import RemoteError

    out = []
    try:
        listed = registry.load(settings=settings)
    except Exception as exc:  # noqa: BLE001
        errors.append(clean(f"the clusters could not be read: {exc}", 200))
        return out
    waiting = _waiting_batches(errors)
    for name, cluster in listed.items():
        batches = waiting.get(name, [])
        if not cluster.usable and not batches:
            continue
        conn, listing = None, []
        try:
            if cluster.usable:
                conn = connection.connect(name, settings=settings)
                if not conn.is_live():
                    conn = None
                else:
                    listing = _cached(("jobs", name), JOBS_PROBE_S, conn.jobs)
        except RemoteError as exc:
            errors.append(clean(f"{name}: {exc.code}: {exc}", 200))
            conn, listing = None, []
            if not batches:
                continue
        if conn is None and not batches:
            continue
        rows = [_batch_row(name, entry, conn) for entry in batches]
        rows += [{"id": clean(j.jobid, 32), "kind": clean(j.partition, 24),
                  "label": clean(j.name, 120), "state": clean(j.state, 24),
                  "detail": clean(f"{j.elapsed} of {j.timelimit} · {j.nodes} node(s)", 80)}
                 for j in listing if getattr(j, "source", "") == "active"]
        finished = sum(1 for j in listing if getattr(j, "source", "") != "active")
        out.append({"source": "cluster", "cluster": clean(name, 64),
                    "title": clean(f"{name} ({cluster.target})", 120), "rows": rows,
                    "finished": finished})
    return out


# --- batches that wait for their install ---------------------------------------------
#
# What cluster_submit kept because the cluster did not have the PC's commit
# yet (delfin.remote.jobs.queue_after_install): a row in its cluster's group
# until the job watcher has submitted it -- then its Slurm jobs take its
# place -- or settled it as not submitted, which NEEDS YOU reports. A row's
# id is "batch-<id>", which job.detail takes as a batch, never as a Slurm id.

BATCH_PREFIX = "batch-"
_BATCH_ID = re.compile(r"batch-[0-9a-f]{4,32}")

#: An install's state, as the row says it: (state word, what the detail says).
_INSTALL_LOOK = {
    "running": ("INSTALLING", "installing {sha} ({minutes} min)"),
    "complete": ("SUBMITTING", "{sha} installed; it goes in at the watcher's next look"),
    "stopped": ("STOPPED", "the install of {sha} was stopped (hpc_guard); not submitted"),
}


def _waiting_batches(errors: list[str]) -> dict[str, list[dict]]:
    """The kept batches still waiting, per cluster, oldest first."""
    try:
        from delfin.remote import jobs as remote_jobs
        out: dict[str, list[dict]] = {}
        for entry in remote_jobs.pending():
            out.setdefault(str(entry.get("cluster") or ""), []).append(entry)
        return out
    except Exception as exc:  # noqa: BLE001
        errors.append(clean(f"the batches waiting for an install could not be read: {exc}", 200))
        return {}


def _install_state(name: str, conn, sha: str) -> str:
    """complete | running | stopped | failed | ..., asked at most once per probe
    interval; "" when it could not be asked."""
    from delfin.remote import install
    from delfin.remote.errors import RemoteError
    try:
        return _cached((f"install:{sha}", name), JOBS_PROBE_S, lambda: install.state(conn, sha))
    except RemoteError:
        return ""


def _batch_row(name: str, entry: dict, conn) -> dict:
    sha = str(entry.get("commit") or "")[:8]
    count = len(entry.get("jobs") or [])
    label = f"waiting for install · {count} job{'' if count == 1 else 's'}"
    minutes = max(0, int((time.time() - float(entry.get("queued_at") or time.time())) // 60))
    if conn is None:
        state, detail = "LOGIN", f"waits for a login to {name}, then for the install of {sha}"
    else:
        current = _install_state(name, conn, str(entry.get("commit") or ""))
        if not current:
            state, detail = "", f"the install of {sha} could not be asked"
        elif current in _INSTALL_LOOK:
            state, text = _INSTALL_LOOK[current]
            detail = text.format(sha=sha, minutes=minutes)
        else:
            state, detail = "FAILED", f"the install of {sha} {current}; not submitted"
    return {"id": clean(BATCH_PREFIX + str(entry.get("id") or ""), 40), "kind": "BATCH",
            "label": clean(label, 120), "state": state, "detail": clean(detail, 160)}


def _batch_detail(name: str, row_id: str, settings: Optional[dict]) -> dict:
    """A waiting batch's details: the approved jobs, the commit, the install's
    state and the last lines of its log (when the cluster is logged in)."""
    from delfin.remote import connection
    from delfin.remote import install
    from delfin.remote import jobs as remote_jobs
    from delfin.remote.errors import RemoteError

    if not _BATCH_ID.fullmatch(row_id):
        return {"sections": [], "error": "not a batch"}
    entry = next((e for e in remote_jobs.pending(cluster=name)
                  if BATCH_PREFIX + str(e.get("id") or "") == row_id), None)
    if entry is None:
        return {"sections": [], "error": "that batch no longer waits: see its jobs, or NEEDS YOU"}
    sha = str(entry.get("commit") or "")
    jobs = (entry.get("batch") or {}).get("jobs") or []
    state, tail = "", []
    try:
        conn = connection.connect(name, settings=settings)
        if conn.is_live():
            state = install.state(conn, sha)
            tail = (install.log_tail(conn, sha) or "").splitlines()[-TAIL_LINES:]
        else:
            state = "waits for a login"
    except RemoteError as exc:
        state = f"not asked: {exc.code}"
    sections = [
        _section("Batch", [("Cluster", name), ("Commit", sha[:12]), ("Install", state),
                           ("Waiting since", _when(entry.get("queued_at"))),
                           ("Jobs", str(len(jobs)))]),
        _section("Approved jobs", lines=[
            f"{j.get('folder', '')} · {j.get('cores', '?')} cores · {j.get('mem_mb', '?')} MB"
            f" · {j.get('time', '?')}" for j in jobs] or ["(none)"]),
    ]
    if tail:
        sections.append(_section(f"Install log, last {TAIL_LINES} lines", lines=tail))
    return {"sections": sections, "error": None}


# --- details: what Space on an entry shows -------------------------------------------
#
# One answer shape for every kind of entry: sections, each a title, label
# and value pairs and lines (an output tail). A window draws them alike.
# Asked when an entry is opened and again while it stays open, so each is
# reused for a while: a job's for a few seconds, a cluster's for half a
# minute, its host key (a fresh connection to the server) for five.

DETAIL_S = 5.0
CLUSTER_DETAIL_S = 30.0
HOST_KEY_S = 300.0
TAIL_LINES = 20
_TAIL_BYTES = 16384
_JOB_ID = re.compile(r"[0-9]+(?:_[0-9]+)?")


def _section(title: str, pairs: Iterable[tuple[str, Any]] = (),
             lines: Iterable[Any] = ()) -> dict:
    return {"title": clean(title, 60),
            "pairs": [[clean(k, 40), clean(v, 300)] for k, v in pairs if v not in (None, "")],
            "lines": [clean_line(line) for line in lines]}


def _when(epoch: Any) -> str:
    try:
        return time.strftime("%Y-%m-%d %H:%M", time.localtime(float(epoch)))
    except (TypeError, ValueError, OverflowError):
        return ""


def cluster_detail(name: str, *, settings: Optional[dict] = None) -> dict:
    """``{"sections": [...], "error"}``: the session and its budget, DELFIN's
    installs there with what they lack, and the host key."""
    from cluster_tunnel import api as ctun_api
    from delfin.remote import clusters as registry
    from delfin.remote import connection
    from delfin.remote.errors import RemoteError

    cluster = registry.load(settings=settings).get(name)
    if cluster is None:
        return {"sections": [], "error": clean(f"{name} is not one of DELFIN's clusters", 120)}
    if not cluster.usable:
        sections = [_section("Refused", lines=cluster.problems)]
        if cluster.target:
            # ctun knows it: the key its server presents, to compare with the
            # fingerprint its operators publish before `ctun hostkey --pin`.
            where = {"config_path": cluster.ctun_config}
            sections.append(_cached(("hostkey", name), HOST_KEY_S,
                                    lambda: _host_key_section(name, ctun_api, where)))
        return {"sections": sections, "error": None}
    try:
        conn = connection.connect(name, settings=settings)
        live = conn.is_live()
    except RemoteError as exc:
        return {"sections": [], "error": clean(f"{exc.code}: {exc}", 300)}
    where = {"config_path": cluster.ctun_config}
    sections = [_cached(("session", name), CLUSTER_DETAIL_S,
                        lambda: _session_section(conn, cluster, live, ctun_api, where))]
    if live:
        sections.append(_cached(("installs", name), CLUSTER_DETAIL_S,
                                lambda: _installs_section(conn)))
    sections.append(_cached(("hostkey", name), HOST_KEY_S,
                            lambda: _host_key_section(name, ctun_api, where)))
    return {"sections": sections, "error": None}


def _session_section(conn, cluster, live: bool, ctun_api, where: dict) -> dict:
    pairs: list[tuple[str, Any]] = [("Target", cluster.target)]
    try:
        entry = ctun_api.load_config(cluster.ctun_config).clusters[cluster.name]
        budget = entry.budget
    except Exception:  # noqa: BLE001
        budget = None
    if not live:
        pairs.append(("Login", "nobody is logged in"))
    else:
        try:
            (status,) = ctun_api.status(cluster.name, **where)
            since = _when((status.session or {}).get("start_epoch"))
        except Exception:  # noqa: BLE001
            since = ""
        pairs.append(("Login", f"logged in since {since}" if since else "logged in"))
        try:
            decision = conn.budget()
            if decision.limit is None:
                pairs.append(("Budget", "no session limit set"))
            else:
                left = max(0.0, float(decision.limit) - float(decision.used or 0))
                pairs.append(("Budget", f"{decision.used} of {decision.limit} {decision.unit} used,"
                                        f" {left:g} left"))
        except Exception as exc:  # noqa: BLE001
            pairs.append(("Budget", f"not read: {exc}"))
        try:
            listing = _cached(("jobs", cluster.name), JOBS_PROBE_S, conn.jobs)
            active = sum(1 for j in listing if getattr(j, "source", "") == "active")
            pairs.append(("Jobs", f"{active} active, {len(listing) - active} finished recently"))
        except Exception as exc:  # noqa: BLE001
            pairs.append(("Jobs", f"not read: {exc}"))
        orca = _orca(cluster.name, conn)
        if orca is not None:
            pairs.append(("ORCA", orca["path"] + (f" ({orca['hint']})" if orca["hint"] else "")
                          if orca["found"] else f"not found: {orca['hint']}"))
    if budget is not None:
        pairs.append(("Guard", f"{budget.guard_scope}: {', '.join(budget.guard_commands)}"))
    pairs += [("Work dir", cluster.work_dir), ("Install root", cluster.install_root),
              ("Partition", cluster.default_partition or "the cluster's default")]
    return _section("Session", pairs)


def _installs_section(conn) -> dict:
    from delfin.remote import install
    try:
        here = install.installs(conn)
    except Exception as exc:  # noqa: BLE001
        return _section("DELFIN on it", lines=[f"the installs could not be listed: {exc}"])
    try:
        mine = install.pc_commit().sha
    except Exception:  # noqa: BLE001
        mine = ""
    lines = []
    for sha, state in here.items():
        mark = "  ← this PC's commit" if mine and sha == mine else ""
        lines.append(f"{sha[:12]}  {state}{mark}")
        if state == "complete":
            try:
                lines += [f"    lacks: {w}" for w in install.warnings(conn, sha)]
            except Exception:  # noqa: BLE001
                pass
    if not here:
        lines.append("no install of DELFIN's yet")
    if mine and mine not in here:
        lines.append(f"this PC runs {mine[:12]}: installed with the first submit, once it is published")
    return _section("DELFIN on it", lines=lines)


def _host_key_section(name: str, ctun_api, where: dict) -> dict:
    try:
        status = ctun_api.host_key(name, **where)
    except Exception as exc:  # noqa: BLE001
        return _section("Host key", lines=[f"the server could not be asked: {exc}"])
    pinned = ", ".join(k.fingerprint for k in status.pinned) or "none"
    verdict = ("the server presents the pinned key" if status.matches
               else "the server presents a DIFFERENT key than the pinned one"
               if status.pinned else "not pinned")
    return _section("Host key", [("Presented", f"{status.offered.kind} {status.offered.fingerprint}"),
                                 ("Pinned", pinned), ("Verdict", verdict)])


def job_detail(params: dict, sessions: Iterable[dict], *,
               settings: Optional[dict] = None) -> dict:
    """``{"sections": [...], "error"}`` for one row of jobs.list.

    ``source`` cluster: ``cluster`` and ``id``, which must be one of the
    user's active jobs there -- Slurm's facts, where it runs, and its
    output's tail -- or a batch waiting there for its install
    (``batch-<id>``): its approved jobs, the commit and the install's log. ``source`` agent: ``session`` (an id this service runs:
    its workspace is looked up here, never taken from the window), ``group``
    and ``id`` -- the background item's output tail. ``source`` monitor:
    the watch entry."""
    source = str(params.get("source") or "")
    job_id = str(params.get("id") or "")
    if source == "cluster":
        name = str(params.get("cluster") or "")
        if job_id.startswith(BATCH_PREFIX):
            return _cached(("batch", f"{name}:{job_id}"), DETAIL_S,
                           lambda: _batch_detail(name, job_id, settings))
        return _cached(("job", f"{name}:{job_id}"), DETAIL_S,
                       lambda: _cluster_job_detail(name, job_id, settings))
    if source == "agent":
        return _agent_detail(str(params.get("session") or ""), str(params.get("group") or ""),
                             job_id, sessions)
    if source == "monitor":
        return _monitor_detail(job_id)
    return {"sections": [], "error": clean(f"no such source: {source}", 80)}


def _cluster_job_detail(name: str, job_id: str, settings: Optional[dict]) -> dict:
    from delfin.remote import connection
    from delfin.remote.errors import RemoteError

    if not _JOB_ID.fullmatch(job_id):
        return {"sections": [], "error": "not a Slurm job id"}
    try:
        conn = connection.connect(name, settings=settings)
        listing = _cached(("jobs", name), JOBS_PROBE_S, conn.jobs)
        if not any(str(getattr(j, "jobid", "")) == job_id for j in listing):
            return {"sections": [], "error": clean(f"{job_id} is none of your jobs on {name}", 120)}
        from delfin.remote import slurm
        fields = slurm.job_fields(conn, job_id)
        tail = []
        out = fields.get("StdOut", "")
        if fields.get("JobState", "").startswith("PENDING"):
            tail = ["no output yet: the job has not started"]
        elif out.startswith("/"):
            got = conn.run(["tail", "-c", str(_TAIL_BYTES), "--", out], timeout=30)
            if got.ok:
                tail = (got.stdout or "").splitlines()[-TAIL_LINES:]
            elif "No such file" in (got.stderr or ""):
                tail = ["no output file yet"]
            else:
                tail = [f"not readable: {got.stderr or ''}"]
    except RemoteError as exc:
        return {"sections": [], "error": clean(f"{exc.code}: {exc}", 300)}
    return {"sections": _job_sections(job_id, fields, tail), "error": None}


def _job_sections(job_id: str, f: dict[str, str], tail: list[str]) -> list[dict]:
    state = f.get("JobState", "")
    reason = f.get("Reason", "")
    if reason and reason not in ("None", "(null)"):
        state = f"{state} ({reason})"
    tres = f.get("AllocTRES") or f.get("ReqTRES") or f.get("TRES") or ""
    gpus = re.search(r"gres/gpu(?::[^=,]+)?=(\d+)", tres)
    nodes = f.get("NumNodes", "")
    low, _, high = nodes.partition("-")
    nodes = low if high == low else nodes          # Slurm's "1-1" for a pending job
    resources = ", ".join(x for x in (
        f"{nodes} node(s)" if nodes else "",
        f"{f['NumCPUs']} CPUs" if f.get("NumCPUs") else "",
        f"{gpus.group(1)} GPU(s)" if gpus else "") if x)
    start = f.get("StartTime", "")
    time_used = " of ".join(x for x in (f.get("RunTime", ""), f.get("TimeLimit", "")) if x)
    commit = re.search(r"delfin-commit=([0-9a-f]{7,40})", f.get("Comment", ""))
    err = f.get("StdErr", "")
    sections = [
        _section("Slurm", [("Job", f"{job_id} {f.get('JobName', '')}".strip()), ("State", state),
                           ("Partition", f.get("Partition", "")), ("Resources", resources),
                           ("Submitted", f.get("SubmitTime", "")),
                           ("Started", "" if start in ("", "Unknown", "None") else start),
                           ("Time", time_used)]),
        _section("Where it runs", [("Work dir", f.get("WorkDir", "")),
                                   ("Output", f.get("StdOut", "")),
                                   ("Errors", err if err and err != f.get("StdOut") else ""),
                                   ("DELFIN", f"install {commit.group(1)[:12]}" if commit else "")]),
    ]
    if tail or f.get("StdOut"):
        sections.append(_section(f"Output, last {TAIL_LINES} lines", lines=tail or ["(empty)"]))
    return sections


def _agent_detail(session_id: str, group: str, item_id: str, sessions: Iterable[dict]) -> dict:
    workspace = next((str(s.get("workspace") or "") for s in sessions
                      if str(s.get("sessionId") or "") == session_id), "")
    if not workspace:
        return {"sections": [], "error": "that session does not run in this service"}
    if group not in ("shells", "watches", "agents", "wakeups"):
        return {"sections": [], "error": "no such group"}
    try:
        from delfin.agent import background_view
        text = background_view.peek(workspace, group, item_id, lines=TAIL_LINES)
    except Exception as exc:  # noqa: BLE001
        return {"sections": [], "error": clean(f"not read: {exc}", 200)}
    lines = str(text or "").splitlines()[-TAIL_LINES:]
    return {"sections": [_section(f"Output, last {TAIL_LINES} lines", lines=lines or ["(empty)"])],
            "error": None}


def _monitor_detail(job_id: str) -> dict:
    try:
        from delfin.agent import job_monitor
        rec = (job_monitor.load_watched().get("jobs") or {}).get(job_id)
    except Exception as exc:  # noqa: BLE001
        return {"sections": [], "error": clean(f"not read: {exc}", 200)}
    if not isinstance(rec, dict):
        return {"sections": [], "error": "no longer watched"}
    return {"sections": [_section("Watched", [("Job", job_id), ("Folder", rec.get("folder")),
                                              ("Last state", rec.get("last_state")),
                                              ("Since", _when(rec.get("added_at")))])],
            "error": None}
