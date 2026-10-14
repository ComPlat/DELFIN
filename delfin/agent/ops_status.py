"""``delfin-agent ops status``: a read-only overview of DELFIN's own machinery.

Input: one checkout of a repository (default: the current directory, else
the checkout this package runs from), the process table (``/proc``, rows of
the current uid only), the user's state directory (``~/.delfin``).

Output: a dict (``SCHEMA``) with five parts, and a text rendering of it.

* ``checkouts`` -- the main checkout and every worktree from ``git worktree
  list --porcelain``: branch, or detached with the branch it was detached
  from (HEAD reflog) and the checkout that holds that branch; clean/dirty;
  ahead/behind the upstream as recorded locally (no fetch); the kind
  (main, session worktree under ``.delfin/worktrees``, agent worktree --
  ``delfin-wt-*`` or a ``worktrees`` directory under a hidden directory --
  or other); whether a dashboard server runs from it; stale (directory
  missing or prunable).
* ``dashboards`` -- jupyter-server/voila processes of this uid: port, cwd,
  the checkout that contains the cwd, ``PYTHONPATH``, the reload grace
  (culler flag or ``DELFIN_SESSION_GRACE_SECONDS``), kernel processes
  below the server with their kept-session and running-turn records.
* ``sessions`` -- writer locks in the sessions directory (live, stale,
  elsewhere, with the holder's kernel id when the holder is a kernel),
  kept-session records, running-turn records.
* ``leftovers`` -- index.lock files in these checkouts with the git
  processes of this uid working there, stale session locks, stale and
  orphaned agent worktrees (idle, clean, merged), idle agent worktrees
  that still hold work.
* ``findings`` -- one ``problem -> safe next step`` per defect.

Semantics: nothing is written, signalled or fetched. git runs with
``--no-optional-locks`` (``git status`` otherwise refreshes the index and
writes it), and only reading subcommands. Liveness of a pid is a lookup in
the process table that was read, never ``kill(pid, 0)``. Of a process's
environment only three named variables are read; the rest (tokens, keys) is
never kept. The next steps are printed, never run: acting stays with the
user or the agent's gated tools.

Limits: a shell inside a pid namespace (the agent's sandbox) sees only its
own processes; ``process_view`` says ``namespaced`` then and a finding
names it. The number of browser windows per kernel is known only to the
server and is not reported.

Rejected: asking the running server's REST API for kernels and
connections -- it needs the server token, and a status command that reads
tokens widens what an auto-allowed command can reach.
"""

from __future__ import annotations

import json
import os
import re
import socket
import subprocess
import time
from pathlib import Path
from typing import Any, Iterable, Optional

SCHEMA = "delfin.ops_status/1"

#: Mirrors session_store._LOCK_MAX_AGE_S: a lock older than this is stale.
LOCK_MAX_AGE_S = 3600.0
#: Mirrors resume_server.GRACE_SECONDS. Not imported: importing the
#: dashboard package loads the widget stack, which a status line does not need.
GRACE_DEFAULT_S = 90.0
#: Mirrors turn_record.FRESH_SECONDS.
TURN_FRESH_S = 120.0
#: Environment variables read from a dashboard server; nothing else is.
_ENV_KEYS = ("PYTHONPATH", "DELFIN_SESSION_GRACE_SECONDS",
             "DELFIN_SESSION_DROPPED_SECONDS")

_KERNEL_FILE_RE = re.compile(r"kernel-([A-Za-z0-9-]+)\.json$")
_CHECKOUT_MOVE_RE = re.compile(r"^checkout: moving from (\S+) to (\S+)$")


# ---------------------------------------------------------------------------
# process table
# ---------------------------------------------------------------------------

def _read(path: Path) -> bytes:
    try:
        return path.read_bytes()
    except OSError:
        return b""


def read_process_table(proc_root: Path | str = "/proc",
                       uid: Optional[int] = None) -> list[dict]:
    """Rows ``{pid, ppid, name, argv, cwd}`` of processes owned by *uid*.

    The owner is the real uid in ``status``. Unreadable rows are skipped.
    """
    root = Path(proc_root)
    me = os.getuid() if uid is None else int(uid)
    rows = []
    try:
        entries = list(root.iterdir())
    except OSError:
        return rows
    for entry in entries:
        if not entry.name.isdigit():
            continue
        status = _read(entry / "status").decode("utf-8", "replace")
        fields = {}
        for line in status.splitlines():
            key, _, value = line.partition(":")
            fields[key.strip()] = value.strip()
        try:
            owner = int(fields.get("Uid", "").split()[0])
        except (IndexError, ValueError):
            continue
        if owner != me:
            continue
        argv = [a.decode("utf-8", "replace")
                for a in _read(entry / "cmdline").split(b"\0") if a]
        try:
            cwd = os.readlink(entry / "cwd")
        except OSError:
            cwd = ""
        try:
            ppid = int(fields.get("PPid", "0") or 0)
        except ValueError:
            ppid = 0
        rows.append({"pid": int(entry.name), "ppid": ppid,
                     "name": fields.get("Name", ""), "argv": argv,
                     "cwd": cwd, "_dir": entry})
    rows.sort(key=lambda r: r["pid"])
    return rows


def _named_environ(row: dict) -> dict:
    """Only the variables in ``_ENV_KEYS`` from a process's environment."""
    out = {}
    for item in _read(Path(row["_dir"]) / "environ").split(b"\0"):
        key, sep, value = item.partition(b"=")
        name = key.decode("utf-8", "replace")
        if sep and name in _ENV_KEYS:
            out[name] = value.decode("utf-8", "replace")
    return out


def _process_view(proc_root: Path) -> str:
    """``full``, ``namespaced`` (nested pid namespace) or ``unknown``."""
    text = _read(proc_root / "self" / "status").decode("utf-8", "replace")
    for line in text.splitlines():
        if line.startswith("NSpid:"):
            return "namespaced" if len(line.split()[1:]) > 1 else "full"
    return "unknown"


def _server_kind(argv: list[str]) -> str:
    """``voila`` / ``jupyter`` for a dashboard server's argv, else ""."""
    head = argv[:4]
    for i, tok in enumerate(head):
        base = os.path.basename(tok)
        nxt = head[i + 1] if i + 1 < len(head) else ""
        if base == "voila":
            return "voila"
        if base in ("jupyter-server", "jupyter-lab", "jupyter-notebook"):
            return "jupyter"
        if base == "jupyter" and nxt in ("server", "lab", "notebook"):
            return "jupyter"
        if tok == "-m":
            if nxt == "voila":
                return "voila"
            if nxt in ("jupyter_server", "jupyterlab", "notebook"):
                return "jupyter"
            if nxt == "jupyter" and i + 2 < len(argv) and argv[i + 2] in (
                    "server", "lab", "notebook"):
                return "jupyter"
    return ""


def _kernel_id(argv: list[str]) -> str:
    if not any("ipykernel" in a for a in argv[:4]):
        return ""
    for tok in argv:
        m = _KERNEL_FILE_RE.search(tok)
        if m:
            return m.group(1)
    return ""


def _flag_value(argv: list[str], names: Iterable[str]) -> str:
    names = tuple(names)
    for i, tok in enumerate(argv):
        for name in names:
            if tok.startswith(name + "="):
                return tok.split("=", 1)[1]
            if tok == name and i + 1 < len(argv):
                return argv[i + 1]
    return ""


def _is_git(row: dict) -> bool:
    argv0 = os.path.basename(row["argv"][0]) if row["argv"] else ""
    return argv0 == "git" or argv0.startswith("git-") or row["name"] == "git"


# ---------------------------------------------------------------------------
# git, read-only
# ---------------------------------------------------------------------------

def _git(path: Path | str, *args: str, timeout: float = 10.0) -> Optional[str]:
    """stdout of a READING git command, or None. Never takes a lock."""
    env = dict(os.environ)
    env["GIT_OPTIONAL_LOCKS"] = "0"
    env["GIT_TERMINAL_PROMPT"] = "0"
    try:
        out = subprocess.run(
            ["git", "--no-optional-locks", "-C", str(path), *args],
            capture_output=True, text=True, timeout=timeout, env=env,
            stdin=subprocess.DEVNULL)
    except (OSError, subprocess.SubprocessError):
        return None
    return out.stdout if out.returncode == 0 else None


def _parse_porcelain(raw: str) -> list[dict]:
    entries, cur = [], {}
    for line in raw.splitlines() + [""]:
        if not line.strip():
            if cur.get("worktree"):
                entries.append(cur)
            cur = {}
            continue
        key, _, value = line.partition(" ")
        cur[key] = value if value else True
    return entries


def _status(path: Path) -> dict:
    raw = _git(path, "status", "--porcelain=v2", "--branch")
    if raw is None:
        return {}
    info = {"dirty": False, "changed": 0, "untracked": 0,
            "upstream": "", "ahead": None, "behind": None}
    for line in raw.splitlines():
        if line.startswith("# branch.upstream "):
            info["upstream"] = line.split(" ", 2)[2]
        elif line.startswith("# branch.ab "):
            parts = line.split()
            info["ahead"], info["behind"] = int(parts[2]), -int(parts[3])
        elif line.startswith("? "):
            info["untracked"] += 1
        elif line and not line.startswith("#") and not line.startswith("!"):
            info["changed"] += 1
    info["dirty"] = info["changed"] > 0
    return info


def _detached_from(path: Path, head: str) -> str:
    """The branch a detached checkout was detached at, from its reflog."""
    raw = _git(path, "reflog", "show", "--format=%gs", "-n", "50", "HEAD") or ""
    for subject in raw.splitlines():
        m = _CHECKOUT_MOVE_RE.match(subject.strip())
        if not m:
            continue
        target = m.group(2)
        if _git(path, "show-ref", "--verify", "--quiet",
                f"refs/heads/{target}") is not None:
            return target
        pointing = (_git(path, "for-each-ref", "--points-at", target,
                         "--format=%(refname:short)", "refs/heads") or "")
        names = pointing.split()
        if names:
            return names[0]
        break
    pointing = (_git(path, "for-each-ref", "--points-at", head,
                     "--format=%(refname:short)", "refs/heads") or "").split()
    return pointing[0] if len(pointing) == 1 else ""


def _kind(path: str, is_main: bool) -> str:
    if is_main:
        return "main"
    parts = Path(path).parts
    if "/.delfin/worktrees/" in path + "/":
        return "session"
    if Path(path).name.startswith("delfin-wt-"):
        return "agent"
    for i, part in enumerate(parts[:-1]):
        if part == "worktrees" and i > 0 and parts[i - 1].startswith(".") \
                and parts[i - 1] != ".delfin" and parts[i - 1] != ".git":
            return "agent"
    return "other"


def _inside(child: str, parent: str) -> bool:
    if not child or not parent:
        return False
    child = child.rstrip("/") or "/"
    parent = parent.rstrip("/") or "/"
    return child == parent or child.startswith(parent + "/")


def _owner_checkout(path: str, checkouts: list[str]) -> str:
    """The checkout that most narrowly contains *path*, or ""."""
    best = ""
    for c in checkouts:
        if _inside(path, c) and len(c) > len(best):
            best = c
    return best


def _repo_root(start: Path) -> Optional[Path]:
    out = _git(start, "rev-parse", "--path-format=absolute", "--git-common-dir")
    return Path(out.strip()) if out else None


# ---------------------------------------------------------------------------
# collect
# ---------------------------------------------------------------------------

def _default_repo() -> Path:
    cwd = Path.cwd()
    if _repo_root(cwd) is not None:
        return cwd
    return Path(__file__).resolve().parents[2]


def _read_json(path: Path) -> Optional[dict]:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return data if isinstance(data, dict) else None


def _sessions_dir(home: Optional[Path]) -> Path:
    if home is not None:
        return Path(home) / ".delfin" / "agent_sessions"
    try:
        from . import session_store
        return Path(session_store._SESSIONS_DIR)
    except Exception:
        return Path.home() / ".delfin" / "agent_sessions"


def _record_dir(env: str, home: Optional[Path], name: str) -> Path:
    if home is None and os.environ.get(env):
        return Path(os.environ[env])
    return Path(home if home is not None else Path.home()) / ".delfin" / name


def collect(*, repo: Path | str | None = None,
            proc_root: Path | str = "/proc", uid: Optional[int] = None,
            home: Path | str | None = None, host: Optional[str] = None,
            now: Optional[float] = None) -> dict:
    """The status dict described in the module docstring. Never raises on a
    missing piece: what cannot be read is left empty."""
    now = time.time() if now is None else float(now)
    host = socket.gethostname() if host is None else host
    proc_root = Path(proc_root)
    home_p = Path(home) if home is not None else None
    procs = read_process_table(proc_root, uid)
    by_pid = {p["pid"]: p for p in procs}
    findings: list[dict] = []

    def find(kind: str, problem: str, step: str) -> None:
        findings.append({"kind": kind, "problem": problem, "next_step": step})

    view = _process_view(proc_root)
    if view == "namespaced":
        find("process_view",
             "this shell sees only its own pid namespace; dashboard servers "
             "and git processes outside it are not listed",
             "run delfin-agent ops status in a terminal outside the sandbox")

    # -- checkouts ---------------------------------------------------------
    start = Path(repo) if repo is not None else _default_repo()
    common = _repo_root(start)
    rows: list[dict] = []
    main_path = ""
    if common is not None:
        raw = _git(start, "worktree", "list", "--porcelain") or ""
        for i, e in enumerate(_parse_porcelain(raw)):
            path = str(e["worktree"])
            branch = str(e.get("branch", "") or "").replace("refs/heads/", "")
            rows.append({
                "path": path, "kind": _kind(path, i == 0),
                "head": str(e.get("HEAD", "") or ""), "branch": branch,
                "detached": bool(e.get("detached")), "bare": bool(e.get("bare")),
                "exists": Path(path).is_dir(),
                "prunable": e.get("prunable") if e.get("prunable") else "",
                "locked": e.get("locked") if e.get("locked") else "",
            })
        main_path = rows[0]["path"] if rows else ""
    paths = [r["path"] for r in rows]
    branch_holder = {r["branch"]: r["path"] for r in rows if r["branch"]}

    # -- dashboards (needed by the checkout rows) --------------------------
    note = _read_json(Path(home_p if home_p is not None else Path.home())
                      / ".delfin" / "dashboard_here.json") or {}
    kept_dir = _record_dir("DELFIN_SESSION_RECORD_DIR", home_p, "kept_sessions")
    turn_dir = _record_dir("DELFIN_TURN_RECORD_DIR", home_p, "running_turns")
    kept = []
    for f in sorted(kept_dir.glob("*.json")) if kept_dir.is_dir() else []:
        rec = _read_json(f) or {}
        kept.append({"session_name": str(rec.get("session_name") or f.stem),
                     "kernel_id": str(rec.get("kernel_id") or ""),
                     "pid": int(rec.get("pid") or 0),
                     "host": str(rec.get("host") or ""),
                     "updated_at": float(rec.get("updated_at") or 0),
                     "path": str(f)})
    turns = {}
    for f in sorted(turn_dir.glob("*.json")) if turn_dir.is_dir() else []:
        rec = _read_json(f) or {}
        kid = str(rec.get("kernel_id") or f.stem)
        try:
            age = now - f.stat().st_mtime
        except OSError:
            age = float("inf")
        turns[kid] = {"kernel_id": kid, "fresh": age < TURN_FRESH_S,
                      "path": str(f)}
    kept_by_kernel = {k["kernel_id"]: k for k in kept if k["kernel_id"]}

    kernels_by_pid = {p["pid"]: _kernel_id(p["argv"]) for p in procs}
    dashboards = []
    for p in procs:
        kind = _server_kind(p["argv"])
        if not kind:
            continue
        env = _named_environ(p)
        port = _flag_value(p["argv"], ("--port", "--ServerApp.port",
                                       "--Voila.port"))
        if not port and int(note.get("pid") or 0) == p["pid"] \
                and str(note.get("host") or host) == host:
            port = str(note.get("port") or "")
        grace = _flag_value(p["argv"],
                            ("--MappingKernelManager.cull_idle_timeout",)) \
            or env.get("DELFIN_SESSION_GRACE_SECONDS", "")
        kernels = []
        for k in procs:
            kid = kernels_by_pid.get(k["pid"])
            if not kid:
                continue
            up, hops = k["ppid"], 0
            while up and up != p["pid"] and hops < 4:
                up = by_pid.get(up, {}).get("ppid", 0)
                hops += 1
            if up != p["pid"]:
                continue
            kernels.append({
                "kernel_id": kid, "pid": k["pid"],
                "kept_as": kept_by_kernel.get(kid, {}).get("session_name", ""),
                "turn_running": bool(turns.get(kid, {}).get("fresh")),
            })
        cwd = p["cwd"]
        dashboards.append({
            "pid": p["pid"], "kind": kind,
            "port": int(port) if str(port).isdigit() else None,
            "cwd": cwd, "cwd_deleted": cwd.endswith(" (deleted)"),
            "checkout": _owner_checkout(cwd, paths),
            "root_dir": _flag_value(p["argv"], ("--ServerApp.root_dir",)),
            "pythonpath": [x for x in env.get("PYTHONPATH", "").split(":") if x],
            "grace_s": float(grace) if _is_number(grace) else None,
            "kernels": kernels,
        })
    dash_by_checkout: dict[str, list[int]] = {}
    for d in dashboards:
        if d["checkout"]:
            dash_by_checkout.setdefault(d["checkout"], []).append(d["pid"])
        if d["cwd_deleted"]:
            find("dashboard",
                 f"dashboard server pid {d['pid']} runs from a deleted "
                 f"directory {d['cwd'][:-len(' (deleted)')]}",
                 "stop it from its terminal (Ctrl+C) and start delfin-voila "
                 "again from an existing checkout")

    # -- per-checkout git detail ------------------------------------------
    git_rows = [p for p in procs if _is_git(p)]

    def _git_pids_in(checkout: str, gitdir: str) -> list[int]:
        out = []
        for g in git_rows:
            places = [g["cwd"]]
            argv = g["argv"]
            for i, tok in enumerate(argv[:-1]):
                if tok == "-C":
                    c = argv[i + 1]
                    places.append(c if os.path.isabs(c)
                                  else os.path.join(g["cwd"], c))
            if any(_owner_checkout(pl, paths) == checkout or _inside(pl, gitdir)
                   for pl in places if pl):
                out.append(g["pid"])
        return out

    index_locks = []
    stale_worktrees = []
    orphaned = []
    idle_with_work = []
    process_cwds = [p["cwd"] for p in procs]
    main_branch = rows[0]["branch"] if rows else ""
    for r in rows:
        r["dashboards"] = dash_by_checkout.get(r["path"], [])
        r["stale"] = (not r["exists"]) or bool(r["prunable"])
        r.update({"dirty": None, "changed": 0, "untracked": 0, "upstream": "",
                  "ahead": None, "behind": None, "detached_from": "",
                  "held_by": "", "same_commit_as_branch": None,
                  "branch_head": ""})
        if r["stale"]:
            stale_worktrees.append(r["path"])
            find("worktree",
                 f"worktree {r['path']} is "
                 + ("gone from disk" if not r["exists"] else "prunable")
                 + " but still registered",
                 f"git -C {main_path} worktree prune (removes only the "
                 "registration of worktrees whose directory is gone)")
            continue
        if r["bare"]:
            continue
        r.update(_status(Path(r["path"])))
        if r["detached"]:
            frm = _detached_from(Path(r["path"]), r["head"])
            r["detached_from"] = frm
            holder = branch_holder.get(frm, "")
            r["held_by"] = holder if holder != r["path"] else ""
            if frm:
                tip = (_git(r["path"], "rev-parse", "--verify", "--quiet",
                            f"refs/heads/{frm}") or "").strip()
                r["branch_head"] = tip
                r["same_commit_as_branch"] = bool(tip) and tip == r["head"]
            if frm and r["held_by"]:
                who = (f"this dashboard (pid {r['dashboards'][0]})"
                       if r["dashboards"] else r["path"])
                if r["same_commit_as_branch"]:
                    at = f"the same commit {r['head'][:7]}"
                else:
                    at = (f"{r['head'][:7]} while {frm} is at "
                          f"{r['branch_head'][:7]}")
                holder_row = next((x for x in rows if x["path"] == r["held_by"]),
                                  {})
                pre = (f"git -C {main_path} worktree prune, then "
                       if holder_row.get("stale") else
                       f"git -C {r['held_by']} checkout --detach, then ")
                find("held_branch",
                     f"{frm} is held by {r['held_by']}; {who} runs detached "
                     f"at {at}",
                     f"to take the name here: {pre}git -C {r['path']} "
                     f"switch {frm}")
        gitdir = (_git(r["path"], "rev-parse", "--absolute-git-dir") or "").strip()
        r["gitdir"] = gitdir
        if gitdir:
            lock = Path(gitdir) / "index.lock"
            if lock.exists():
                try:
                    age = max(0.0, now - lock.stat().st_mtime)
                except OSError:
                    age = None
                pids = _git_pids_in(r["path"], gitdir)
                index_locks.append({"path": str(lock), "checkout": r["path"],
                                    "age_s": age, "git_pids": pids,
                                    "stale": not pids})
                if pids:
                    find("index_lock",
                         f"{lock}: a git process of yours (pid "
                         f"{', '.join(map(str, pids))}) is running in "
                         f"{r['path']}",
                         "wait for it to finish; do not remove the lock")
                else:
                    find("index_lock",
                         f"stale {lock} ({_age(age)} old, no git process of "
                         f"yours in {r['path']}); git commands there fail",
                         f"rm {lock}")
        if r["kind"] == "agent" and not any(
                _owner_checkout(c, paths) == r["path"] for c in process_cwds):
            # Orphaned only when removing it loses nothing: no process has
            # its working directory there, nothing uncommitted or
            # untracked, HEAD contained in the main branch. An idle
            # worktree that still holds work is listed, not flagged: an
            # agent whose shell runs elsewhere looks idle while it works.
            merged = bool(main_branch) and _git(
                r["path"], "merge-base", "--is-ancestor", "HEAD",
                f"refs/heads/{main_branch}") is not None
            if not r["dirty"] and not r["untracked"] and merged:
                orphaned.append(r["path"])
                find("orphaned_worktree",
                     f"agent worktree {r['path']} has no process of yours "
                     f"working in it and nothing that is not in "
                     f"{main_branch}",
                     f"git -C {main_path} worktree remove {r['path']}")
            else:
                idle_with_work.append(r["path"])
        if (r["branch"] and r["behind"] and not r["dirty"]
                and r["upstream"] and not r["ahead"]):
            find("behind",
                 f"{r['path']} ({r['branch']}) is {r['behind']} commit(s) "
                 f"behind {r['upstream']} as of the last fetch",
                 f"git -C {r['path']} pull --ff-only")

    # -- sessions ----------------------------------------------------------
    sdir = _sessions_dir(home_p)
    locks = []
    for f in sorted(sdir.glob("*.lock")) if sdir.is_dir() else []:
        rec = _read_json(f)
        pid = int((rec or {}).get("pid") or 0)
        ts = float((rec or {}).get("ts") or 0)
        lhost = str((rec or {}).get("host") or "")
        age = now - ts if ts else None
        elsewhere = bool(lhost) and lhost != host
        if not rec or not pid:
            state, why = "stale", "holds no readable pid"
        elif age is not None and age >= LOCK_MAX_AGE_S:
            state, why = "stale", f"older than {int(LOCK_MAX_AGE_S)} s"
        elif elsewhere:
            state, why = "elsewhere", f"held on {lhost}; only its age can be judged here"
        elif pid not in by_pid:
            state, why = "stale", f"holder pid {pid} is not running here"
        else:
            state, why = "live", ""
        locks.append({"session_id": f.stem, "path": str(f), "pid": pid,
                      "host": lhost, "age_s": age, "state": state, "why": why,
                      "kernel_id": kernels_by_pid.get(pid, "") if state == "live"
                      else ""})
    stale_locks = [lk["path"] for lk in locks if lk["state"] == "stale"]
    for lk in locks:
        if lk["state"] == "stale":
            find("session_lock",
                 f"stale writer lock for session {lk['session_id']} ({lk['why']})",
                 f"nothing required: the next save breaks it; or rm {lk['path']}")
    live_kernels = {k for k in kernels_by_pid.values() if k}
    for k in kept:
        here = not k["host"] or k["host"] == host
        k["alive"] = (k["kernel_id"] in live_kernels or k["pid"] in by_pid) \
            if here else None
        if here and k["alive"] is False:
            find("kept_session",
                 f"kept session {k['session_name']} names kernel "
                 f"{k['kernel_id'] or '?'}, which is not running here",
                 "nothing to resume there; the dashboard removes the record "
                 f"when it next starts, or rm {k['path']}")
    return {
        "schema": SCHEMA,
        "generated_at": now,
        "host": host,
        "repository": main_path,
        "process_view": view,
        "checkouts": rows,
        "dashboards": dashboards,
        "sessions": {
            "dir": str(sdir),
            "locks": locks,
            "open": [lk["session_id"] for lk in locks
                     if lk["state"] in ("live", "elsewhere")],
            "kept": kept,
            "running_turns": list(turns.values()),
            "reload_grace_default_s": GRACE_DEFAULT_S,
        },
        "leftovers": {
            "index_locks": index_locks,
            "stale_session_locks": stale_locks,
            "stale_worktrees": stale_worktrees,
            "orphaned_agent_worktrees": orphaned,
            "idle_agent_worktrees_with_work": idle_with_work,
        },
        "findings": findings,
    }


def _is_number(text: Any) -> bool:
    try:
        float(text)
        return True
    except (TypeError, ValueError):
        return False


def _age(seconds: Optional[float]) -> str:
    if seconds is None:
        return "?"
    if seconds < 120:
        return f"{int(seconds)} s"
    if seconds < 7200:
        return f"{int(seconds // 60)} min"
    if seconds < 172800:
        return f"{int(seconds // 3600)} h"
    return f"{int(seconds // 86400)} d"


# ---------------------------------------------------------------------------
# text
# ---------------------------------------------------------------------------

def _checkout_line(r: dict) -> str:
    if r["stale"]:
        state = "STALE (" + ("directory missing" if not r["exists"]
                             else f"prunable: {r['prunable']}") + ")"
        return f"  {r['kind']:<8} {r['path']}  {state}"
    if r["bare"]:
        return f"  {r['kind']:<8} {r['path']}  (bare)"
    if r["detached"]:
        where = f"detached @ {r['head'][:7]}"
        if r["detached_from"]:
            held = f", held by {r['held_by']}" if r["held_by"] else ""
            where += f" (from {r['detached_from']}{held})"
    else:
        where = r["branch"]
    clean = "dirty" if r["dirty"] else "clean"
    if r["untracked"]:
        clean += f", {r['untracked']} untracked"
    up = ""
    if r["upstream"]:
        up = f"  ahead {r['ahead']}, behind {r['behind']} vs {r['upstream']}"
    dash = (f"  [dashboard pid {', '.join(map(str, r['dashboards']))}]"
            if r["dashboards"] else "")
    return f"  {r['kind']:<8} {r['path']}  {where}  {clean}{up}{dash}"


def format_text(st: dict) -> str:
    out = [f"DELFIN ops status -- host {st['host']}, "
           f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(st['generated_at']))}"
           " (read-only; nothing was changed)",
           f"repository {st['repository'] or '(none found)'}", "",
           "Checkouts"]
    out += [_checkout_line(r) for r in st["checkouts"]] or ["  (none)"]
    out += ["", f"Dashboard servers (yours; process view: {st['process_view']})"]
    if not st["dashboards"]:
        out.append("  (none running)")
    for d in st["dashboards"]:
        grace = (f"reload grace {d['grace_s']:g} s" if d["grace_s"] is not None
                 else "no reload grace (kernels are not culled)")
        out.append(f"  pid {d['pid']}  {d['kind']}  port {d['port'] or '?'}  "
                   f"cwd {d['cwd']}  checkout {d['checkout'] or '(not one of these)'}"
                   f"  {grace}")
        if d["pythonpath"]:
            out.append(f"      PYTHONPATH {':'.join(d['pythonpath'])}")
        for k in d["kernels"]:
            tags = []
            if k["kept_as"]:
                tags.append(f"kept as '{k['kept_as']}' (exempt from the grace)")
            if k["turn_running"]:
                tags.append("turn running (exempt while it runs)")
            out.append(f"      kernel {k['kernel_id']} pid {k['pid']}"
                       + (f"  {'; '.join(tags)}" if tags else ""))
    s = st["sessions"]
    out += ["", f"Sessions ({s['dir']})"]
    if not s["locks"] and not s["kept"]:
        out.append("  (no writer locks, no kept sessions)")
    for lk in s["locks"]:
        extra = f"  kernel {lk['kernel_id']}" if lk["kernel_id"] else ""
        why = f"  -- {lk['why']}" if lk["why"] else ""
        out.append(f"  lock {lk['session_id']}  {lk['state']}  pid {lk['pid']}"
                   f"@{lk['host'] or '?'}  age {_age(lk['age_s'])}{extra}{why}")
    for k in s["kept"]:
        alive = {True: "alive", False: "gone", None: "elsewhere"}[k["alive"]]
        out.append(f"  kept {k['session_name']}  kernel {k['kernel_id']}  {alive}")
    for t in s["running_turns"]:
        out.append(f"  turn running in kernel {t['kernel_id']}"
                   + ("" if t["fresh"] else " (record stale)"))
    left = st["leftovers"]
    if left.get("idle_agent_worktrees_with_work"):
        out += ["", "Idle agent worktrees that still hold work (not flagged)"]
        out += [f"  {p}" for p in left["idle_agent_worktrees_with_work"]]
    out += ["", "Findings"]
    if not st["findings"]:
        out.append("  nothing to report")
    for f in st["findings"]:
        out.append(f"  ! {f['problem']} → {f['next_step']}")
    return "\n".join(out)


def main(*, json_out: bool = False, repo: str = "") -> int:
    st = collect(repo=repo or None)
    if json_out:
        print(json.dumps(st, indent=1, default=str))
    else:
        print(format_text(st))
    return 0
