"""A wave: many DELFIN agent sessions on one task list, run and watched.

Where ``swarm.yaml`` describes packages and roles that people start and
review, a wave is the orchestrator's unit: a list of tasks, each one
agent session in a git worktree of its own, started as soon as what it
depends on has ended, watched until it is final. Everything here is a
library function with JSON-able input and output, so the CLI
(``delfin-agent swarm wave``), the MCP server (``swarm/mcp_server.py``)
and the dashboard reach the same code.

A run lives in ``<root>/<run>/``:

  ``tasks.json``         the task list as started (see :func:`validate_tasks`)
  ``state.json``         one entry per started task: pid, worktree, branch,
                         log, start time, the task, how often resumed
  ``<tid>/``             the task's worktree, on branch ``wave/<run>/<tid>``
  ``<tid>.log``          the session's stream-json log
  ``events.jsonl``       what happened, one line each (:func:`observe`)
  ``<tid>.final`` etc.   the stage markers (``swarm/stages.py``)

``root`` is ``agent.swarm.waves_root`` (default ``~/delfin-waves``),
never under ``~/.delfin``, which every member's cage hides.

Each member runs ``python -m delfin.agent.cli chat -p BRIEF`` from the
installed DELFIN -- never from the tree it edits -- in
``bypassPermissions`` with ``--unattended --isolate``, so it runs only
caged, with wall-clock and turn ceilings. Its worktree's push URL leads
nowhere: a member commits, the orchestrator reviews and pushes. On Linux
with a user service manager each member is a scope in one slice with
CPU, memory and task ceilings, at the lowest CPU and I/O priority.

Measured lessons carried over from the waves run by hand
(2026-11-05..24): a lock held for a starter's whole life deadlocked the
stage that resumes members, so every state change takes the lock, reads
the state fresh and writes it back; an ended member whose launcher has
not reaped it is a zombie and counts as ended; a task that waits for
other members starts once all of them ended, or once ``after_min`` of
them ended and ``after_wait_s`` passed; the work of other members is
copied in (``collect``), because their scratch/ is in their worktree.
"""
from __future__ import annotations

import contextlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Iterable

_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,63}$")
_REF = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,63}/[A-Za-z0-9][A-Za-z0-9._-]{0,63}(\.final)?$")
_PUSH_NOWHERE = "file:///nonexistent/wave-members-commit-locally-only"

_DEFAULTS = {
    "waves_root": "~/delfin-waves",
    "max_parallel": 4,
    # Members alive across ALL runs; 0 is no limit. max_parallel is one
    # run's: three runs started side by side ran 21 members (2026-11-24).
    "max_live": 0,
    "provider_caps": {},
    "member_cpu": "200%",
    "member_memory": "12G",
    # 256 was too few: the egress proxy's threads in a busy member failed
    # with "can't start new thread" (two crashes, 2026-11-23/24).
    "member_tasks": 1024,
    "slice": "delfin-swarm.slice",
    "min_available_gb": 8,
    # At the envelope's stop line, not above it: the guard stops the whole
    # slice when the load passes 0.885 x cores, and a start gate left at 1.0
    # made a band (339.8..384 on a 384-core host) where the tick kept
    # launching members onto a host it was about to stop. This number is
    # duplicated rather than imported -- the guard imports this module, so
    # the dependency cannot run the other way -- and the two are pinned
    # together by test_the_start_gate_stops_before_the_envelope_does.
    "max_load_per_core": 0.885,
    "seconds": 5400,
    "turns": 60,
    # A checkout of its own the members run DELFIN from (PYTHONPATH), so a
    # member never runs the tree someone is editing -- with an editable
    # install the installed DELFIN IS a working tree. "" runs the installed
    # one. ``refresh_harness`` puts it on origin/main at each start.
    "harness": "",
}


class WaveError(ValueError):
    """A wave request cannot be served, with the reason in prose."""


# -- settings and places ----------------------------------------------------------

def settings() -> dict:
    """``agent.swarm`` from the user's settings over the defaults."""
    try:
        from delfin.user_settings import load_settings
        own = ((load_settings() or {}).get("agent") or {}).get("swarm") or {}
    except Exception:
        own = {}
    merged = dict(_DEFAULTS)
    merged.update({k: v for k, v in own.items() if k in _DEFAULTS})
    return merged


def root() -> Path:
    path = Path(str(settings()["waves_root"])).expanduser()
    delfin_home = Path.home() / ".delfin"
    if path == delfin_home or delfin_home in path.parents:
        raise WaveError(f"agent.swarm.waves_root ({path}) lies under ~/.delfin, which "
                        "every member's cage hides: the members could not see "
                        "their own worktrees. Choose another directory.")
    return path


def run_dir(run: str) -> Path:
    if not _ID.match(str(run or "")):
        raise WaveError(f"a run name is letters, digits, '.', '_' and '-' (got {run!r})")
    return root() / run


def _state_path(run: str) -> Path:
    return run_dir(run) / "state.json"


def load_state(run: str) -> dict:
    try:
        return json.loads(_state_path(run).read_text())
    except (OSError, ValueError):
        return {"agents": {}}


def _save_state(run: str, st: dict) -> None:
    p = _state_path(run)
    tmp = p.with_suffix(".tmp")
    tmp.write_text(json.dumps(st, indent=1))
    os.replace(tmp, p)


@contextlib.contextmanager
def _state_lock(run: str):
    """The run's state for one change: taken, read fresh, given back. A lock
    held for a starter's whole life deadlocked the stage that resumes its
    members (measured 2026-11-24)."""
    with _flock(run_dir(run) / ".state.lock"):
        yield


@contextlib.contextmanager
def _flock(path: Path, *, wait: bool = True):
    """An exclusive lock on ``path`` for the block; yields False instead of
    waiting when ``wait`` is false and another process holds it."""
    import fcntl
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        try:
            fcntl.flock(fh, fcntl.LOCK_EX | (0 if wait else fcntl.LOCK_NB))
        except BlockingIOError:
            yield False
            return
        try:
            yield True
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)


def runs() -> list[str]:
    """Every run under the root, newest first."""
    try:
        found = [p.parent for p in root().glob("*/state.json")]
    except WaveError:
        return []
    return [p.name for p in sorted(found, key=lambda p: p.stat().st_mtime, reverse=True)]


def runs_oldest_first() -> list[str]:
    """Every run, the one whose first member started earliest first; a run
    with no member yet comes last. The tick starts members in this order, so
    a run that is further along gets a free place before a newer one: by the
    newest state file first, three new runs took every place while an older
    run's integrators waited behind them (2026-11-25)."""
    def age(run: str) -> tuple:
        began = [float(a["started"]) for a in load_state(run).get("agents", {}).values()
                 if isinstance(a, dict) and a.get("started")]
        return (min(began) if began else float("inf"), run)
    return sorted(runs(), key=age)


# -- tasks -----------------------------------------------------------------------

def validate_tasks(tasks: Any) -> list[dict]:
    """The task list, checked. Each task is an object with

      ``id``            letters, digits, '.', '_', '-'; unique
      ``brief``         what the member is asked to do (required)
      ``base``          the git ref its branch starts from (default origin/main)
      ``provider``      ``model``     the endpoint and model (default: the
                        user's agent settings)
      ``seconds``       ``turns``     the session's ceilings
      ``after``         ``run/tid`` or ``run/tid.final`` it waits for
      ``after_min``     ``after_wait_s``   start once that many of them
                        ended and the first of those ended that long ago
      ``collect``       ``run/tid`` whose scratch/ is copied into
                        ``scratch/arms/<tid>/`` at the start
      ``attach``        files copied into ``scratch/``
      ``hide``          paths the member's cage hides besides ~/.delfin
      ``no_realrun``    the stage asks for no real run (a reviewer)

    Raises WaveError with every problem named."""
    if not isinstance(tasks, list) or not tasks:
        raise WaveError("tasks must be a non-empty list of task objects")
    problems, seen = [], set()
    for i, t in enumerate(tasks):
        where = f"task {i}"
        if not isinstance(t, dict):
            problems.append(f"{where}: not an object")
            continue
        tid = str(t.get("id") or "")
        where = f"task {tid or i}"
        if not _ID.match(tid):
            problems.append(f"{where}: id must be letters, digits, '.', '_', '-'")
        elif tid in seen:
            problems.append(f"{where}: id used twice")
        seen.add(tid)
        if not str(t.get("brief") or "").strip():
            problems.append(f"{where}: brief is empty")
        for key in ("after", "collect"):
            refs = t.get(key) or []
            if not isinstance(refs, list) or not all(isinstance(r, str) and _REF.match(r) for r in refs):
                problems.append(f"{where}: {key} must be a list of 'run/tid' (after may end in .final)")
            elif key == "collect" and any(r.endswith(".final") for r in refs):
                problems.append(f"{where}: collect names a member, not its .final")
        for key in ("seconds", "turns", "after_min", "after_wait_s"):
            if key in t and (not isinstance(t[key], int) or t[key] < 0):
                problems.append(f"{where}: {key} must be a whole number >= 0")
        for key in ("attach", "hide"):
            if key in t and not (isinstance(t[key], list) and all(isinstance(x, str) for x in t[key])):
                problems.append(f"{where}: {key} must be a list of paths")
    if problems:
        raise WaveError("; ".join(problems))
    return tasks


# -- processes ---------------------------------------------------------------------

def alive(pid: int) -> bool:
    """Running and not a zombie: a member whose launcher has not reaped it
    yet has ended all the same (a stage stalled on one for an hour,
    2026-11-24)."""
    try:
        os.kill(int(pid), 0)
    except (OSError, ValueError, TypeError):
        return False
    try:
        return Path(f"/proc/{pid}/stat").read_text().rsplit(")", 1)[1].split()[0] != "Z"
    except OSError:
        return True                              # no /proc (macOS): kill(0) answered


def ended(ref: str) -> bool:
    """Whether ``run/tid`` started and no longer runs; ``run/tid.final``
    whether its stage marked it final."""
    run, _, tid = ref.partition("/")
    if tid.endswith(".final"):
        return (run_dir(run) / tid).exists()
    a = load_state(run)["agents"].get(tid)
    return a is not None and not alive(a["pid"])


def _ended_at(ref: str) -> float:
    run, _, tid = ref.partition("/")
    if tid.endswith(".final"):
        p = run_dir(run) / tid
    else:
        a = load_state(run)["agents"].get(tid) or {}
        p = Path(a.get("log") or run_dir(run) / f"{tid}.log")
    try:
        return p.stat().st_mtime
    except OSError:
        return time.time()


def ready(task: dict, now: float | None = None) -> bool:
    """All of ``after`` ended -- or at least ``after_min`` of them, and the
    first of those ended ``after_wait_s`` ago."""
    refs = task.get("after") or []
    if not refs:
        return True
    done = [r for r in refs if ended(r)]
    if len(done) == len(refs):
        return True
    need = int(task.get("after_min") or 0)
    if need and len(done) >= need:
        first = min(_ended_at(r) for r in done)
        return ((now or time.time()) - first) >= int(task.get("after_wait_s") or 0)
    return False


def _room(cfg: dict) -> tuple[bool, str]:
    """Whether this host has room for one more member."""
    try:
        load1 = os.getloadavg()[0]
        cores = os.cpu_count() or 1
        if load1 > cores * float(cfg["max_load_per_core"]):
            return False, f"load {load1:.0f} on {cores} cores"
    except OSError:
        pass
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemAvailable:"):
                gb = int(line.split()[1]) / 1048576
                if gb < float(cfg["min_available_gb"]):
                    return False, f"{gb:.0f} GB available"
    except OSError:
        pass
    return True, ""


def live_members(run_names: Iterable[str] | None = None) -> list[dict]:
    """Live members of the given runs (default: every run)."""
    out = []
    for run in (run_names if run_names is not None else runs()):
        for tid, a in load_state(run)["agents"].items():
            if alive(a["pid"]):
                out.append({"run": run, "id": tid, "pid": a["pid"],
                            "provider": str((a.get("task") or {}).get("provider") or "")})
    return out


def _git(repo: Path, *args: str, cwd: Path | None = None) -> str:
    """git through ``safe_git``: a member's worktree is the member's to write,
    and a program its repository configuration named would run here, outside
    the member's cage, with the person's rights."""
    from delfin.agent import safe_git
    return safe_git.check_output(list(args), str(cwd or repo)).strip()


def _hidden(task: dict) -> list[str]:
    hide = [str(Path(p).expanduser()) for p in task.get("hide") or [] if p]
    hide.append(str(Path.home() / ".delfin"))
    return list(dict.fromkeys(hide))


def _copy_in(src_dir: Path, dest: Path, depth: int = 2, limit: int = 2_000_000) -> None:
    for f in (src_dir.glob("**/*") if src_dir.is_dir() else []):
        rel = f.relative_to(src_dir)
        if (f.is_file() and not f.is_symlink() and len(rel.parts) <= depth
                and rel.parts[0] != "arms" and f.stat().st_size < limit):
            (dest / rel).parent.mkdir(parents=True, exist_ok=True)
            (dest / rel).write_bytes(f.read_bytes())


def member_argv(run: str, task: dict, *, session: str = "", cfg: dict | None = None,
                python: str | None = None) -> list[str]:
    """The command one member runs (without its scope)."""
    cfg = cfg or settings()
    wt = run_dir(run) / task["id"]
    cmd = [python or sys.executable, "-m", "delfin.agent.cli", "chat", "-p", task["brief"],
           "--cwd", str(wt), "--permission-mode", "bypassPermissions", "--unattended",
           "--isolate", "--max-run-seconds", str(int(task.get("seconds", cfg["seconds"]))),
           "--max-run-turns", str(int(task.get("turns", cfg["turns"]))),
           "--output-format", "stream-json"]
    if task.get("provider"):
        cmd += ["--provider", str(task["provider"])]
    if task.get("model"):
        cmd += ["--model", str(task["model"])]
    if session:
        cmd += ["--session", session]
    return cmd


def _scoped(cmd: list[str], cfg: dict) -> list[str]:
    """The member in a scope of the swarm slice where a user service manager
    runs one, at the lowest CPU and I/O priority everywhere.

    ``OOMPolicy=continue``: a command the member runs that outgrows the
    scope's memory is killed, and the member lives on to read exit 137. The
    default (stop) ended the whole scope on that kill, so a full test suite
    run by an integrator killed the integrator (w42i, 2026-11-24)."""
    nice = ["nice", "-n", "15"] + (["ionice", "-c", "3"] if shutil.which("ionice") else [])
    if shutil.which("systemd-run") and os.environ.get("XDG_RUNTIME_DIR"):
        return ["systemd-run", "--user", "--scope", "--quiet", f"--slice={cfg['slice']}",
                "-p", f"CPUQuota={cfg['member_cpu']}", "-p", f"MemoryMax={cfg['member_memory']}",
                "-p", f"TasksMax={int(cfg['member_tasks'])}", "-p", "OOMPolicy=continue",
                *nice, *cmd]
    return [*nice, *cmd]


def _launch(run: str, task: dict, repo: Path, session: str = "") -> dict:
    cfg = settings()
    tid = task["id"]
    d = run_dir(run)
    wt = d / tid
    branch = f"wave/{run}/{tid}"
    if not wt.exists():
        _git(repo, "worktree", "add", "-q", "-b", branch, str(wt), str(task.get("base") or "origin/main"))
    # A member commits; its push goes nowhere (the cage holds no credential
    # either): the orchestrator reviews the branch and pushes it.
    _git(repo, "config", "extensions.worktreeConfig", "true")
    _git(repo, "config", "--worktree", "remote.origin.pushurl", _PUSH_NOWHERE, cwd=wt)
    for src in task.get("attach") or []:
        p = Path(src).expanduser()
        if p.is_file():
            (wt / "scratch").mkdir(exist_ok=True)
            (wt / "scratch" / p.name).write_bytes(p.read_bytes())
    for ref in task.get("collect") or []:
        orun, _, otid = ref.partition("/")
        _copy_in(run_dir(orun) / otid / "scratch", wt / "scratch" / "arms" / otid)
    log = d / f"{tid}.log"
    env = dict(os.environ, DELFIN_SWARM_RUN=run, DELFIN_NO_BROWSER="1",
               DELFIN_AGENT_LABEL=f"{run}/{tid}",
               DELFIN_CAGE_HIDE=os.pathsep.join(_hidden(task)))
    if cfg.get("harness"):
        env["PYTHONPATH"] = str(Path(str(cfg["harness"])).expanduser())
    for name in list(env):
        if name.startswith("VSCODE_") or name in ("BROWSER", "TERM_PROGRAM", "LD_PRELOAD",
                                                  "OPENAI_API_KEY"):
            env.pop(name, None)
    with open(log, "ab") as fh:
        # cwd is the run root, not the worktree: ``python -m`` puts the cwd
        # first on sys.path, and a member must run the installed DELFIN,
        # never the tree it edits.
        proc = subprocess.Popen(_scoped(member_argv(run, task, session=session, cfg=cfg), cfg),
                                cwd=str(root()), env=env, stdin=subprocess.DEVNULL,
                                stdout=fh, stderr=subprocess.STDOUT, start_new_session=True)
    return {"pid": proc.pid, "branch": branch, "worktree": str(wt), "log": str(log),
            "started": time.time(), "task": task}


def refresh_harness(cfg: dict | None = None) -> dict:
    """Put the harness checkout on origin/main. Not under a running member,
    not over local changes, and not onto a main that does not import: then
    the old commit is restored and the start refused (WaveError).

    On 2026-11-23 the harness stood 3,699 commits behind main and a whole
    wave ran without a month of fixes; separate was never meant to mean
    stale."""
    cfg = cfg or settings()
    path = str(cfg.get("harness") or "")
    if not path:
        return {"harness": "", "moved": False, "why": "no harness configured"}
    h = Path(path).expanduser()
    from delfin.agent import safe_git
    users = [m for m in live_members()]
    if users:
        return {"harness": str(h), "moved": False, "why": f"{len(users)} member(s) run from it"}
    if safe_git.out(h, "status", "--porcelain").strip():
        return {"harness": str(h), "moved": False, "why": "local changes in the harness"}
    old = safe_git.out(h, "rev-parse", "HEAD").strip()
    safe_git.run(["fetch", "-q", "origin"], h, timeout=300)
    new = safe_git.out(h, "rev-parse", "origin/main").strip()
    if not new or new == old:
        return {"harness": str(h), "moved": False, "why": f"on origin/main ({old[:9]})"}
    safe_git.run(["checkout", "-q", "--detach", new], h)
    # Imported FROM the harness: without the check a harness missing the
    # package imports the installed DELFIN and the smoke passes on nothing.
    smoke = subprocess.run(
        [sys.executable, "-c", "import os, delfin.agent.cli, delfin.agent.api_client, delfin; "
                               "print(os.path.realpath(delfin.__file__))"],
        cwd=str(h), env=dict(os.environ, PYTHONPATH=str(h)), capture_output=True, text=True)
    if smoke.returncode != 0 or not smoke.stdout.strip().startswith(os.path.realpath(h)):
        safe_git.run(["checkout", "-q", "--detach", old], h)
        raise WaveError(f"the harness on origin/main {new[:9]} does not import; back on {old[:9]} "
                        f"and the start is refused: {(smoke.stderr or smoke.stdout)[-400:]}")
    return {"harness": str(h), "moved": True, "from": old[:9], "to": new[:9]}


def _refuse_on_a_login_node() -> None:
    try:
        from delfin.hpc_guard import on_login_node
        from .job import login_node_opt_in
        if on_login_node() and not login_node_opt_in():
            raise WaveError("a wave does not start on a cluster login node: run it inside a "
                            "SLURM job, or allow it with hpc.swarm_on_login_node")
    except WaveError:
        raise
    except Exception:
        pass


# -- the verbs ---------------------------------------------------------------------

def start(run: str, tasks: list[dict], repo: str | Path, *, parallel: int | None = None,
          wait: bool = True, poll_s: float = 30.0) -> dict:
    """Start the run's tasks as they become ready, at most ``parallel`` at a
    time and within each provider's cap; with ``wait`` until every task has
    started. A task already in the state is not started twice, so a start
    after an interruption continues. Returns ``{"started": [...], "waiting": [...]}``."""
    _refuse_on_a_login_node()
    tasks = validate_tasks(tasks)
    repo = Path(repo).resolve()
    _git(repo, "rev-parse", "--git-dir")
    cfg = settings()
    caps = {str(k): int(v) for k, v in (cfg["provider_caps"] or {}).items()}
    d = run_dir(run)
    d.mkdir(parents=True, exist_ok=True)
    if not (d / "tasks.json").exists() and cfg.get("harness"):
        record_harness = refresh_harness(cfg)        # a new run: the harness on main
    else:
        record_harness = None
    (d / "tasks.json").write_text(json.dumps(tasks, indent=1))
    with _state_lock(run):
        st = load_state(run)
        st.setdefault("repo", str(repo))
        # The run keeps its limit: the tick starts the rest of the queue
        # with no limit of its own, and used the default (2026-11-24).
        if parallel:
            st["parallel"] = int(parallel)
        _save_state(run, st)
    parallel = int(parallel or st.get("parallel") or cfg["max_parallel"])
    max_live = int(cfg.get("max_live") or 0)
    started: list[str] = []
    while True:
        st = load_state(run)
        queue = [t for t in tasks if t["id"] not in st["agents"]]
        if not queue:
            break
        running = [a for a in st["agents"].values() if alive(a["pid"])]
        launched = False
        everyone = live_members()
        if (len(running) < parallel and (not max_live or len(everyone) < max_live)
                and _room(cfg)[0]):
            live = {}
            for m in everyone:
                live[m["provider"]] = live.get(m["provider"], 0) + 1
            pick = next((t for t in queue if ready(t)
                         and live.get(str(t.get("provider") or ""), 0)
                         < caps.get(str(t.get("provider") or ""), 10 ** 6)), None)
            if pick is not None:
                with _state_lock(run):
                    st = load_state(run)
                    if pick["id"] not in st["agents"]:
                        st["agents"][pick["id"]] = _launch(run, pick, repo)
                        _save_state(run, st)
                        started.append(pick["id"])
                        record(run, "started", pick["id"])
                        launched = True
        if not wait:
            if not launched:
                break
            continue
        if not launched:
            time.sleep(poll_s)
    st = load_state(run)
    out = {"run": run, "started": started,
           "waiting": [t["id"] for t in tasks if t["id"] not in st["agents"]]}
    if record_harness:
        out["harness"] = record_harness
    return out


def add_tasks(run: str, tasks: list[dict]) -> dict:
    """Append tasks to a run that exists -- a later stage the orchestrator
    decides on once it has seen the first (reviewers after a lander, a fixer
    after the reviewers). An id the run already has is refused; the tick
    starts the new tasks when they are ready."""
    d = run_dir(run)
    try:
        current = json.loads((d / "tasks.json").read_text())
    except (OSError, ValueError):
        raise WaveError(f"run {run} has no task list; start it first")
    known = {t["id"] for t in current}
    clash = sorted(known & {str(t.get("id")) for t in tasks if isinstance(t, dict)})
    if clash:
        raise WaveError(f"run {run} already has tasks {', '.join(clash)}")
    merged = validate_tasks(current + list(tasks))
    with _state_lock(run):
        tmp = d / "tasks.json.tmp"
        tmp.write_text(json.dumps(merged, indent=1))
        os.replace(tmp, d / "tasks.json")
    for t in tasks:
        record(run, "task_added", t["id"])
    return {"run": run, "added": [t["id"] for t in tasks], "tasks": len(merged)}


def _last_result(log: Path) -> dict:
    try:
        lines = log.read_text(errors="replace").splitlines()
    except OSError:
        return {}
    for line in reversed(lines):
        if line.startswith('{"type": "result"'):
            try:
                return json.loads(line)
            except ValueError:
                return {}
    return {}


def resume(run: str, tid: str, note: str = "", *, fresh: bool = False,
           wait_for_room: bool = True, poll_s: float = 60.0) -> dict:
    """Continue member ``tid``'s session in its worktree, with ``note`` (the
    orchestrator's review, or what the stage found) ahead of its task. A
    member still running is not resumed. ``fresh`` starts a new session
    that is told to read the branch's commits first."""
    with _state_lock(run):
        a = load_state(run)["agents"].get(tid)
    if a is None:
        raise WaveError(f"{run}/{tid} was never started")
    if alive(a["pid"]):
        return {"run": run, "id": tid, "resumed": False, "why": "still running"}
    cfg = settings()
    while wait_for_room and not _room(cfg)[0]:
        time.sleep(poll_s)
    res = _last_result(Path(a["log"]))
    sid = "" if fresh else str(res.get("session_id") or "")
    task = dict(a["task"])
    base_s = int(task.setdefault("base_seconds", task.get("seconds", cfg["seconds"])))
    # The session's run clock continues across a resume: grant the base
    # allotment again on top of what was used.
    task["seconds"] = int(task.get("seconds", base_s)) + base_s
    brief = task.setdefault("orig_brief", task["brief"])
    lead = ("You take over a task in this worktree. Read `git log` and `git status` first: "
            "the commits there are earlier work on it. " if fresh else
            "Continue your task from where you stopped. ")
    task["brief"] = (lead + "Finish it completely in this turn; do not announce next steps, "
                     "do them.\n\n" + (f"Review of your work so far:\n{note.strip()}\n\n" if note.strip() else "")
                     + "Your task was:\n" + brief)
    with _state_lock(run):
        st = load_state(run)
        repo = Path(st.get("repo") or ".")
        cur = st["agents"][tid]
        if alive(cur["pid"]):
            return {"run": run, "id": tid, "resumed": False, "why": "started meanwhile"}
        info = _launch(run, task, repo, session=sid)
        info["resumes"] = int(cur.get("resumes", 0)) + 1
        st["agents"][tid] = info
        _save_state(run, st)
    record(run, "resumed", tid, note=note[:400])
    return {"run": run, "id": tid, "resumed": True, "pid": info["pid"], "session": sid[:8]}


def stop(run: str, tid: str = "") -> dict:
    """End the members of a run (or one member): their scopes' process
    groups get SIGTERM, then SIGKILL. Worktrees, commits and logs stay."""
    import signal
    st = load_state(run)
    targets = [t for t in st["agents"] if not tid or t == tid]
    stopped = []
    for t in targets:
        pid = st["agents"][t]["pid"]
        if not alive(pid):
            continue
        for sig in (signal.SIGTERM, signal.SIGKILL):
            try:
                os.killpg(pid, sig)
            except OSError:
                break
            for _ in range(50):
                if not alive(pid):
                    break
                time.sleep(0.1)
            if not alive(pid):
                break
        stopped.append(t)
        record(run, "stopped", t)
    return {"run": run, "stopped": stopped}


def status(run: str) -> dict:
    """One row per task: state, commits on its branch, uncommitted files,
    whether it wrote scratch/PR.md, how its last session ended, its stage
    markers, minutes since its log last grew."""
    st = load_state(run)
    try:
        tasks = json.loads((run_dir(run) / "tasks.json").read_text())
    except (OSError, ValueError):
        tasks = [a["task"] for a in st["agents"].values()]
    repo = Path(st.get("repo") or ".")
    rows = []
    for t in tasks:
        tid = t["id"]
        a = st["agents"].get(tid)
        d = run_dir(run)
        marks = [m for m in ("asked", "realrun", "final") if (d / f"{tid}.{m}").exists()]
        if a is None:
            rows.append({"id": tid, "state": "waiting" if not ready(t) else "queued",
                         "marks": marks})
            continue
        wt = Path(a["worktree"])
        try:
            commits = int(_git(repo, "rev-list", "--count",
                               f"{t.get('base') or 'origin/main'}..{a['branch']}") or 0)
        except (subprocess.CalledProcessError, OSError, ValueError):
            commits = None
        try:
            dirty = len([x for x in _git(repo, "status", "--porcelain", cwd=wt).splitlines() if x])
        except (subprocess.CalledProcessError, OSError):
            dirty = None
        log = Path(a["log"])
        try:
            idle_min = round((time.time() - log.stat().st_mtime) / 60, 1)
        except OSError:
            idle_min = None
        res = _last_result(log)
        rows.append({"id": tid, "state": "running" if alive(a["pid"]) else "ended",
                     "pid": a["pid"], "branch": a["branch"], "commits": commits, "dirty": dirty,
                     "pr_md": (wt / "scratch" / "PR.md").is_file(),
                     "stop_reason": res.get("stop_reason"), "error": (res.get("error") or "")[:200],
                     "resumes": a.get("resumes", 0), "marks": marks, "idle_min": idle_min,
                     "provider": t.get("provider"), "model": t.get("model")})
    return {"run": run, "root": str(run_dir(run)), "members": rows}


# -- events ------------------------------------------------------------------------

def _events_path(run: str) -> Path:
    return run_dir(run) / "events.jsonl"


def record(run: str, kind: str, tid: str = "", **data: Any) -> dict:
    """Append one event to the run's log and hand it to every subscriber of
    the swarm's event feed (``status.emit_event``)."""
    ev = {"t": round(time.time(), 1), "run": run, "kind": kind, "id": tid, **data}
    with open(_events_path(run), "a") as fh:
        fh.write(json.dumps(ev) + "\n")
    try:
        from . import status as _status
        _status.emit_event({"kind": f"wave_{kind}", "run": run, "member": tid, **data})
    except Exception:
        pass
    return ev


def observe(run: str, *, silent_min: float = 60.0) -> list[dict]:
    """Record what changed since the last look, once each: a member that
    ended (with how), one whose log has not grown for ``silent_min``, a
    traceback in a member's log, a stage marker, every task final.
    Returns the new events.

    One look at a time per run: the tick (a service), ``watch`` and a tick
    run by hand all call it, and two looks at once recorded each event
    twice and collided on the cursor file (seen 2026-11-24)."""
    with _flock(run_dir(run) / ".observe.lock"):
        return _observe(run, silent_min)


def _observe(run: str, silent_min: float) -> list[dict]:
    d = run_dir(run)
    seen_p = d / ".observed.json"
    try:
        seen = json.loads(seen_p.read_text())
    except (OSError, ValueError):
        seen = {}
    new = []
    st = load_state(run)
    for tid, a in st["agents"].items():
        log = Path(a["log"])
        if not alive(a["pid"]):
            key = f"ended:{tid}:{a['pid']}"
            if key not in seen:
                res = _last_result(log)
                wt = Path(a["worktree"])
                new.append(record(run, "ended", tid, stop_reason=res.get("stop_reason"),
                                  error=(res.get("error") or "")[:300],
                                  pr_md=(wt / "scratch" / "PR.md").is_file(),
                                  realrun_asked=(wt / "scratch" / "realrun.json").is_file()))
                seen[key] = time.time()
        else:
            try:
                idle = (time.time() - log.stat().st_mtime) / 60
            except OSError:
                idle = 0
            key = f"silent:{tid}:{a['pid']}:{int(idle // silent_min)}"
            if idle >= silent_min and key not in seen:
                new.append(record(run, "silent", tid, idle_min=round(idle)))
                seen[key] = time.time()
        try:
            n_tb = log.read_text(errors="replace").count("Traceback (most recent call last)")
        except OSError:
            n_tb = 0
        key = f"crash:{tid}:{n_tb}"
        if n_tb and key not in seen:
            new.append(record(run, "traceback", tid, count=n_tb))
            seen[key] = time.time()
    for m in d.glob("*.final"):
        key = f"final:{m.name}"
        if key not in seen:
            new.append(record(run, "final", m.name[:-len(".final")],
                              why=m.read_text(errors="replace").strip()[:200]))
            seen[key] = time.time()
    try:
        tasks = json.loads((d / "tasks.json").read_text())
        if tasks and all((d / f"{t['id']}.final").exists() for t in tasks) and "all_final" not in seen:
            new.append(record(run, "all_final"))
            seen["all_final"] = time.time()
    except (OSError, ValueError, KeyError):
        pass
    tmp = seen_p.with_name(f".observed.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(seen))
    os.replace(tmp, seen_p)
    return new


def events(run: str, since: int = 0) -> dict:
    """The run's events from line ``since`` on, and the cursor to pass next."""
    try:
        lines = _events_path(run).read_text().splitlines()
    except OSError:
        lines = []
    return {"events": [json.loads(x) for x in lines[since:] if x.strip()], "cursor": len(lines)}


def watch(run_names: list[str], *, timeout_s: float = 1800.0, poll_s: float = 20.0,
          cursors: dict | None = None) -> dict:
    """Wait until something happens in any of the runs, then say what:
    ``{"events": [...], "cursors": {...}}``; empty events at the timeout.
    Pass the cursors back to continue where the last answer stopped."""
    cursors = dict(cursors or {})
    deadline = time.time() + max(0.0, timeout_s)
    for run in run_names:
        cursors.setdefault(run, events(run)["cursor"])
    while True:
        out = []
        for run in run_names:
            observe(run)
            got = events(run, int(cursors.get(run, 0)))
            out += got["events"]
            cursors[run] = got["cursor"]
        if out or time.time() >= deadline:
            return {"events": out, "cursors": cursors}
        time.sleep(min(poll_s, max(0.5, deadline - time.time())))
