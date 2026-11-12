"""A swarm runs in a SLURM job: check, submit, watch, stop, and the job body.

A cluster's login nodes are shared by everyone; a swarm is up to eight
agent sessions at once and belongs on a compute node (issue #640). The
person's rule: interactive sessions may stay on the login node, a swarm
always runs in its own SLURM job, at most :data:`MAX_MEMBERS_PER_JOB`
member sessions per job. Nobody has to write an sbatch script for it:

* :func:`preflight` checks a spec before anything is submitted -- the
  spec rules (delfin.agent.swarm.spec), the member count, the git
  repository, branches left over from an earlier run, uncommitted work --
  and says in plain words what will happen (:func:`confirmation_text`);
* :func:`submit` hands the swarm to sbatch through
  ``delfin.scheduler_client`` (the same way ``delfin.dashboard_job`` does),
  with a frozen copy of the checked spec, so an edit made while the job
  waits for a node does not change what starts;
* :func:`jobs` and :func:`cancel` list and stop swarm jobs;
* :func:`run_in_job` is the job body on the compute node: it creates the
  member worktrees (``launch.start_members``), starts one member session
  per tmux window, spaced, and holds the allocation until every member
  session has ended.

State lives under ``~/.delfin/swarm_jobs`` (0700): one record per job, the
frozen spec, the job log and, per member, its first prompt and terminal
log. The job inherits the submitting environment as sbatch's default does
(never ``--export`` with named variables, a known cause of failed job
starts on this cluster).
"""

from __future__ import annotations

import json
import os
import secrets
import shlex
import shutil
import sys
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

#: The most member sessions one swarm job holds (the person's decision,
#: 11.11.: more than eight sessions on one node starve each other).
MAX_MEMBERS_PER_JOB = 8

#: The SLURM job name every swarm job carries.
JOB_NAME = "delfin-swarm"

#: How long a swarm job asks for unless ``hpc.swarm_time`` says otherwise.
DEFAULT_TIME = "24:00:00"

#: Cores and memory per member session when the settings name none, and
#: the floor and ceiling of the whole job.
CPUS_PER_MEMBER = 2
MEM_GB_PER_MEMBER = 8
MIN_CPUS, MAX_CPUS = 4, 16
MIN_MEM_GB, MAX_MEM_GB = 16, 64

#: A job the scheduler does not list yet within this long after its submit
#: is still starting, not gone (squeue lags behind sbatch).
SUBMIT_GRACE_S = 120.0

#: Scheduler argv heads as module constants (the direct-call guard scans
#: call sites); every run goes through ``delfin.scheduler_client``.
_SBATCH = "sbatch"
_SQUEUE_HEAD = ("squeue", "-h", "-u")
_SCANCEL = "scancel"

#: Where a member session waits before its start while the node is short
#: of memory (``limits.memory_gate``).
MEMORY_GATE_GB = 4.0
MEMORY_GATE_WAIT_S = 600.0


class SwarmJobError(RuntimeError):
    """A swarm job cannot be submitted or run, with the reason in prose."""


# --------------------------------------------------------------------------
# Where a start goes
# --------------------------------------------------------------------------

def host_kind(env: dict | None = None) -> str:
    """``login_node``, ``slurm_job`` or ``workstation`` (one definition,
    ``delfin.agent.calc_policy.host_kind``)."""
    from delfin.agent.calc_policy import host_kind as _host_kind
    return _host_kind(env)


#: Opt-in for members on a login node: env DELFIN_SWARM_ON_LOGIN_NODE=1, or
#: the user setting hpc.swarm_on_login_node: true. Off by default.
LOGIN_NODE_OPT_IN_ENV = "DELFIN_SWARM_ON_LOGIN_NODE"


def _user_settings() -> dict:
    try:
        from delfin.hpc_guard import _read_user_settings
        return _read_user_settings()
    except Exception:
        return {}


def login_node_opt_in(env: dict | None = None,
                      settings: dict | None = None) -> bool:
    """Whether the person allowed swarm members on a cluster login node.

    Off by default: a login node submits the swarm as a job. Allowed, the
    members run there as interactive agent sessions do -- and, as those,
    every test, build and calculation they start is sent to SLURM (the
    login-node bash gate holds for every session on a login node)."""
    env = os.environ if env is None else env
    raw = str(env.get(LOGIN_NODE_OPT_IN_ENV, "")).strip().lower()
    if raw:
        return raw in ("1", "true", "yes", "on")
    settings = _user_settings() if settings is None else settings
    hpc = settings.get("hpc") if isinstance(settings, dict) else None
    return (hpc or {}).get("swarm_on_login_node") is True


def login_node_room(members: int, *, settings: dict | None = None,
                    on_host: int | None = None) -> str:
    """Why *members* more sessions may not open on this login node, or "".

    The same cap as interactive sessions (hpc.login_node_sessions, at most
    MAX_MEMBERS_PER_JOB): the login node's per-user memory limit is shared
    by every session the person has there."""
    settings = _user_settings() if settings is None else settings
    hpc = settings.get("hpc") if isinstance(settings, dict) else None
    try:
        limit = int((hpc or {}).get("login_node_sessions") or MAX_MEMBERS_PER_JOB)
    except (TypeError, ValueError):
        limit = MAX_MEMBERS_PER_JOB
    limit = max(1, min(MAX_MEMBERS_PER_JOB, limit))
    if on_host is None:
        try:
            import socket
            from delfin.agent import session_presence
            me = socket.gethostname().split(".")[0]
            on_host = sum(1 for r in session_presence.open_sessions()
                          if str(r.get("host") or "").split(".")[0] == me)
        except Exception:
            on_host = 0
    if on_host + members <= limit:
        return ""
    return (f"{members} member session{'s' if members != 1 else ''} would make "
            f"{on_host + members} agent sessions on this login node; at most "
            f"{limit} run here (hpc.login_node_sessions). Start it as a SLURM "
            "job instead (drop --here), or end sessions first.")


def route(env: dict | None = None, *, prefer_job: bool = False,
          here: bool = False) -> str:
    """``"job"`` when the swarm is submitted as a SLURM job, ``"here"`` when
    it starts in this process as before.

    A login node submits, unless the person allowed members there
    (*here*, or :func:`login_node_opt_in`); their tests and builds then
    still go to SLURM. Inside a job (a dashboard job, an interactive srun)
    the start stays where it is unless *prefer_job* asks for a job of its
    own and sbatch is here. A workstation or a private server (no
    scheduler, or ``hpc.host_kind`` naming it) always starts here."""
    kind = host_kind(env)
    if kind == "login_node":
        return "here" if (here or login_node_opt_in(env)) else "job"
    if prefer_job and kind == "slurm_job" and shutil.which(_SBATCH):
        return "job"
    return "here"


# --------------------------------------------------------------------------
# State
# --------------------------------------------------------------------------

def state_dir() -> Path:
    """``~/.delfin/swarm_jobs``, created 0700. Resolved per call."""
    d = Path.home() / ".delfin" / "swarm_jobs"
    d.mkdir(parents=True, mode=0o700, exist_ok=True)
    try:
        os.chmod(d, 0o700)
    except OSError:
        pass
    return d


def _write_private(path: Path, text: str, mode: int = 0o600) -> None:
    tmp = path.with_name(path.name + ".tmp")
    fd = os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, mode)
    with os.fdopen(fd, "w", encoding="utf-8") as fh:
        fh.write(text)
    os.chmod(tmp, mode)
    tmp.replace(path)


def record_path(job: str) -> Path:
    return state_dir() / f"{job}.json"


def progress_path(job: str) -> Path:
    return state_dir() / f"{job}.progress.json"


def log_path(job: str) -> Path:
    return state_dir() / f"swarm-{job}.log"


def member_dir(job: str) -> Path:
    d = state_dir() / str(job)
    d.mkdir(mode=0o700, exist_ok=True)
    return d


def _read_json(path: Path) -> dict:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


def records() -> list:
    """Every swarm job record, oldest first (progress merged in)."""
    out = []
    for p in sorted(state_dir().glob("*.json")):
        if p.name.endswith(".progress.json"):
            continue
        r = _read_json(p)
        if r.get("job"):
            r["progress"] = _read_json(progress_path(str(r["job"])))
            out.append(r)
    out.sort(key=lambda r: float(r.get("submitted_at") or 0))
    return out


# --------------------------------------------------------------------------
# Resources
# --------------------------------------------------------------------------

def _hpc_settings(settings: dict | None) -> dict:
    if settings is None:
        from delfin import hpc_guard
        settings = hpc_guard._read_user_settings()
    hpc = settings.get("hpc") if isinstance(settings, dict) else None
    return hpc if isinstance(hpc, dict) else {}


def job_resources(n_members: int, settings: dict | None = None) -> dict:
    """Partition, cores, memory and time of a swarm job.

    From the ``hpc`` entry of the user settings (``swarm_partition``,
    ``swarm_cpus``, ``swarm_mem``, ``swarm_time``; the partition falls back
    to ``agent_partition``), else sized to the swarm: two cores and 8 GB
    per member session, between 4 and 16 cores and 16 and 64 GB."""
    hpc = _hpc_settings(settings)
    n = max(1, int(n_members))
    cpus = min(MAX_CPUS, max(MIN_CPUS, CPUS_PER_MEMBER * n))
    mem = min(MAX_MEM_GB, max(MIN_MEM_GB, MEM_GB_PER_MEMBER * n))
    job = {
        "partition": str(hpc.get("swarm_partition")
                         or hpc.get("agent_partition") or "cpu"),
        "cpus": cpus,
        "mem": f"{mem}G",
        "time": DEFAULT_TIME,
    }
    for key in ("cpus", "mem", "time"):
        value = hpc.get(f"swarm_{key}")
        if value not in (None, ""):
            job[key] = int(value) if key == "cpus" else str(value)
    return job


# --------------------------------------------------------------------------
# Checks before a start
# --------------------------------------------------------------------------

def _git(repo: Path, *args: str):
    from delfin.agent.safe_git import run as _safe_git_run
    return _safe_git_run(list(args), cwd=repo)


def repo_of(spec_path: Path | str) -> Path:
    """The git repository the spec belongs to (its folder's top level).
    Raises :class:`SwarmJobError` when the folder is not in one."""
    folder = Path(spec_path).expanduser().resolve().parent
    if not folder.is_dir():
        raise SwarmJobError(f"The folder {folder} does not exist.")
    out = _git(folder, "rev-parse", "--show-toplevel")
    top = (out.stdout or "").strip()
    if out.returncode != 0 or not top:
        raise SwarmJobError(
            f"{folder} is not inside a git repository. A swarm gives every "
            "member its own git worktree and branch, so the spec must live "
            "in the repository the swarm works on.")
    return Path(top)


#: Environment a member's commits take their author and committer from.
IDENTITY_ENV = ("GIT_AUTHOR_NAME", "GIT_AUTHOR_EMAIL",
                "GIT_COMMITTER_NAME", "GIT_COMMITTER_EMAIL")

#: Environment variable that names a swarm run (delfin.agent.swarm_memory).
RUN_ENV = "DELFIN_SWARM_RUN"


def git_identity(repo: Path | str) -> dict:
    """``{"name", "email"}`` the members' commits are made as, or ``{}``.

    What ``git config`` gives in *repo*: the repository's own setting
    first, then the person's global one. Without both, git guesses an
    identity from the login and the host name -- on a compute node that
    puts the node's name into every commit, and the commit is no longer
    the person's. A swarm therefore never starts on a guess."""
    out = {}
    for key in ("name", "email"):
        r = _git(Path(repo), "config", "--get", f"user.{key}")
        value = (r.stdout or "").strip() if r.returncode == 0 else ""
        if not value:
            return {}
        out[key] = value
    return out


def identity_env(identity: dict) -> dict:
    """The GIT_AUTHOR_*/GIT_COMMITTER_* variables for *identity*."""
    if not identity.get("name") or not identity.get("email"):
        return {}
    return {"GIT_AUTHOR_NAME": identity["name"],
            "GIT_COMMITTER_NAME": identity["name"],
            "GIT_AUTHOR_EMAIL": identity["email"],
            "GIT_COMMITTER_EMAIL": identity["email"]}


def _no_identity(repo: Path | str) -> str:
    return (f"No git identity for the members' commits in {repo}: git would "
            "guess one from the login and the compute node's host name. Set "
            "it first: `git config user.name \"Your Name\"` and `git config "
            "user.email you@example.org` (in the repository, or with "
            "--global).")


#: What each member permission mode means, for the confirmation.
MODE_WORDS = {
    "acceptEdits": "writes its package's files without asking",
    "default": "asks before every write",
    "plan": "plans only, writes nothing",
}


def _parse(spec_text: str) -> dict:
    import yaml
    from .spec import SPEC_VERSION, SpecError
    if not isinstance(spec_text, str) or not spec_text.strip():
        raise SpecError("The spec is empty. Load a swarm.yaml or start from "
                        "the example.")
    try:
        data = yaml.safe_load(spec_text)
    except yaml.YAMLError as exc:
        raise SpecError(f"swarm.yaml is not valid YAML: {exc}") from exc
    if not isinstance(data, dict):
        raise SpecError("swarm.yaml must be a mapping at the top level "
                        "(goal:, packages:, ...)")
    data.setdefault("spec_version", SPEC_VERSION)
    return data


def preflight(spec_text: str, spec_path: Path | str, *,
              env: dict | None = None, prefer_job: bool = False,
              route_: str | None = None,
              running: list | None = None,
              settings: dict | None = None) -> dict:
    """Everything that can be checked before a start, in plain words.

    Returns ``{"ok", "errors", "warnings", "route", "repo", "head",
    "members", "resources", "goal", "spacing_s", "steps"}``. ``ok`` is
    False when anything in ``errors`` would make the start fail; warnings
    are worth reading but do not block. Starts nothing, writes nothing.
    *running* is the list :func:`jobs` gave (None: not asked)."""
    from . import launch as _launch
    from .spec import SpecError

    out: dict[str, Any] = {
        "ok": False, "errors": [], "warnings": [], "members": [],
        "route": route_ or route(env, prefer_job=prefer_job),
        "repo": "", "head": "", "resources": {}, "goal": "",
        "spacing_s": 15, "steps": [], "spec_path": str(spec_path or ""),
        "modes": {}, "git_identity": {},
    }
    errors, warnings = out["errors"], out["warnings"]
    try:
        spec = _parse(spec_text)
        plan = _launch.plan_launch(spec)
    except (SpecError, _launch.LaunchError) as exc:
        errors.append(str(exc))
        return out
    members = plan["members"]
    out["members"] = [m["name"] for m in members]
    out["modes"] = {m["name"]: m["permission_mode"] for m in members}
    out["goal"] = plan["goal"]
    out["spacing_s"] = plan["start_spacing_s"]
    if out["route"] == "here" and host_kind(env) == "login_node":
        # Allowed on this login node (hpc.swarm_on_login_node): within the
        # session cap, and every test and build still goes to SLURM.
        why = login_node_room(len(members), settings=settings)
        if why:
            errors.append(why)
        else:
            warnings.append("The members run on this login node; every test, "
                            "build and calculation they start goes to SLURM.")
    try:
        spacing = float(plan["start_spacing_s"])
        if spacing < 0:
            raise ValueError
    except (TypeError, ValueError):
        errors.append("limits.start_spacing_s must be a number of seconds "
                      f"(0 or more), not {plan['start_spacing_s']!r}.")
    n = len(members)
    if out["route"] == "job" and n > MAX_MEMBERS_PER_JOB:
        errors.append(
            f"This swarm has {n} member sessions; one swarm job holds at "
            f"most {MAX_MEMBERS_PER_JOB}. Cut packages or roles (a package "
            "without 'roles' is one builder), or split the packages over "
            "two spec files and start each one as its own swarm.")
    if not spec_path:
        errors.append("Name the spec file (swarm.yaml in the repository "
                      "the swarm works on).")
        return out
    try:
        repo = repo_of(spec_path)
    except SwarmJobError as exc:
        errors.append(str(exc))
        return out
    out["repo"] = str(repo)
    head = _git(repo, "rev-parse", "--short", "HEAD")
    out["head"] = (head.stdout or "").strip() if head.returncode == 0 else ""
    if not out["head"]:
        errors.append(f"{repo} has no commit yet: the members branch off "
                      "HEAD, so commit something first.")
    out["git_identity"] = git_identity(repo)
    if not out["git_identity"]:
        errors.append(_no_identity(repo))
    for m in members:
        wt = repo / m["worktree"]
        if wt.exists():
            warnings.append(
                f"{m['name']}: its worktree {m['worktree']} exists already; "
                "the member continues there (nothing in it is reset).")
            continue
        ref = _git(repo, "show-ref", "--verify", "--quiet",
                   f"refs/heads/{m['branch']}")
        if ref.returncode == 0:
            errors.append(
                f"{m['name']}: the branch {m['branch']} exists already but "
                f"its worktree {m['worktree']} does not, so the start would "
                "fail. Rename that branch (git branch -m) or rename the "
                "package.")
    dirty = _git(repo, "status", "--porcelain", "--untracked-files=no")
    if dirty.returncode == 0 and (dirty.stdout or "").strip():
        warnings.append(
            f"{repo} has uncommitted changes. The members start from the "
            f"last commit ({out['head'] or 'HEAD'}) and do not see them; "
            "commit first if they should.")
    for r in running or []:
        if str(r.get("repo") or "") == str(repo) and \
                str(r.get("goal") or "") == out["goal"]:
            warnings.append(
                f"A swarm job with the same goal already runs from this "
                f"repository (job {r.get('job')}, {r.get('state', '?')}). "
                "Starting again runs its members twice; watch or stop that "
                "one first.")
    if out["route"] == "job":
        out["resources"] = job_resources(n, settings)
    out["ok"] = not errors
    out["steps"] = _steps(out)
    return out


def _modes_line(pre: dict) -> str:
    modes = pre.get("modes") or {}
    parts = [f"{name}: {modes.get(name, 'default')} "
             f"({MODE_WORDS.get(modes.get(name, 'default'), '')})"
             for name in pre["members"]]
    return ("Permission modes, as the spec names them per package (never "
            "bypass): " + "; ".join(parts) + ". A write outside the "
            "package's own files is refused in every mode.")


def approvals_line() -> str:
    """Where an approval a member still asks for is answered."""
    return ("An approval a member still asks for waits for you, since nobody "
            "sits at its terminal: answer it in the dashboard's Swarm Watch "
            "tab, or list it with `delfin-agent approvals` and answer with "
            "`delfin-agent approvals approve <id> --code <code>` (the id and "
            "code also stand in the member's log).")


def _steps(pre: dict) -> list:
    """What will happen, step by step, for the confirmation."""
    n = len(pre["members"])
    names = ", ".join(pre["members"])
    repo = pre["repo"]
    sp = pre["spacing_s"]
    ident = pre.get("git_identity") or {}
    commits = (f"Their commits are made as {ident['name']} <{ident['email']}> "
               "(your git config)." if ident else "")
    if pre["route"] == "job":
        res = pre["resources"]
        return [
            f"Submit one SLURM job '{JOB_NAME}' (partition {res['partition']}, "
            f"{res['cpus']} cores, {res['mem']} memory, at most "
            f"{res['time']}). It waits in the queue until a compute node is "
            "free; nothing runs on this login node.",
            f"On the node: create {n} git worktree(s) under {repo}/.swarm, "
            f"each on its own branch from {pre['head'] or 'HEAD'}.",
            (f"Start {n} member session(s), {sp} s apart: {names}. "
             + commits).strip(),
            _modes_line(pre),
            approvals_line(),
            "The job ends when every member session has ended, at the time "
            "limit, or when you stop it. Nothing is merged or pushed: you "
            "review each branch yourself.",
        ]
    return [
        f"Create {n} git worktree(s) under {repo}/.swarm, each on its own "
        f"branch from {pre['head'] or 'HEAD'}, with the member's write "
        "scope in its settings.",
        f"Members: {names}. Their sessions are not started here: open each "
        "worktree in an agent session (it reads its role prompt from "
        "delfin/agent/swarm/prompts.py).",
        _modes_line(pre),
    ]


def watch_lines(job: str = "") -> list:
    """How to watch and stop a swarm, the same words everywhere."""
    stop = f"delfin-agent swarm stop {job}" if job else \
        "delfin-agent swarm stop JOB"
    return [
        "Watch it: the dashboard's Swarm Watch tab, '/swarm watch' in the "
        "TUI, or `delfin-agent swarm status` (members) and "
        "`delfin-agent swarm jobs` (job, node, time left).",
        f"Stop it: `{stop}` or the Stop button in the dashboard's Swarm "
        "tab. The member sessions end with the job; their worktrees and "
        "branches stay.",
    ]


def confirmation_text(pre: dict) -> str:
    """The text a person reads before saying yes: problems first, then
    exactly what will happen."""
    lines = []
    if pre["errors"]:
        lines.append("Cannot start yet:")
        lines += [f"  - {e}" for e in pre["errors"]]
    if pre["warnings"]:
        lines.append("Please note:")
        lines += [f"  - {w}" for w in pre["warnings"]]
    if pre["ok"]:
        where = ("as a SLURM job" if pre["route"] == "job"
                 else "here, without a SLURM job")
        lines.append(f"Swarm '{pre['goal']}' starts {where}. What happens:")
        lines += [f"  {i}. {s}" for i, s in enumerate(pre["steps"], 1)]
        if pre["route"] == "job":
            lines += [f"  {w}" for w in watch_lines()]
    return "\n".join(lines)


# --------------------------------------------------------------------------
# Scheduler
# --------------------------------------------------------------------------

def submit(spec_text: str, pre: dict, *, run: Callable | None = None,
           cwd_python: str | None = None) -> dict:
    """Submit the checked swarm as a SLURM job. Returns its record.

    *pre* must be a :func:`preflight` result that is ``ok`` with route
    ``job``; the spec text is frozen into the state directory so the job
    starts exactly what was checked. *run(cmd, env)* replaces the real
    submit (tests)."""
    if not pre.get("ok"):
        raise SwarmJobError("The spec has not passed its checks: "
                            + "; ".join(pre.get("errors") or ["not checked"]))
    if pre.get("route") != "job":
        raise SwarmJobError("This start does not go to a SLURM job.")
    if len(pre["members"]) > MAX_MEMBERS_PER_JOB:
        raise SwarmJobError(f"At most {MAX_MEMBERS_PER_JOB} member sessions "
                            "per swarm job.")
    res = pre["resources"]
    repo = pre["repo"]
    state = state_dir()
    identity = pre.get("git_identity") or {}
    if not identity_env(identity):
        raise SwarmJobError(_no_identity(repo))
    token = secrets.token_hex(6)
    frozen = state / f"spec-{token}.yaml"
    _write_private(frozen, spec_text)
    # The identity checked now is the one the members commit as; it goes to
    # the job in a private file, not on the command line the scheduler
    # shows to everybody.
    ident_file = state / f"spec-{token}.identity.json"
    _write_private(ident_file, json.dumps(
        {"name": identity["name"], "email": identity["email"]}))
    python = cwd_python or sys.executable
    inner = [python, "-m", "delfin.agent.swarm.job", "run",
             "--spec", str(frozen), "--repo", repo,
             "--identity", str(ident_file)]
    cmd = [_SBATCH, "--parsable", f"--job-name={JOB_NAME}",
           f"--partition={res['partition']}", f"--time={res['time']}",
           "--nodes=1", "--ntasks=1", f"--cpus-per-task={int(res['cpus'])}",
           f"--mem={res['mem']}", f"--output={state}/swarm-%j.log",
           f"--chdir={repo}",
           "--wrap=exec " + " ".join(shlex.quote(a) for a in inner)]
    env = dict(os.environ)
    if run is None:
        from delfin import scheduler_client
        out = scheduler_client.query(
            scheduler_client.kind_of_argv(cmd), cmd,
            runner=scheduler_client.child_env_runner(env))
        out = "" if out is scheduler_client.THROTTLED else str(out or "")
    else:
        out = str(run(cmd, env) or "")
    jid = out.strip().split(";")[0]
    if not jid.isdigit():
        raise SwarmJobError(
            f"sbatch gave no job id: {out.strip() or 'no output'}. The spec "
            f"is kept at {frozen}.")
    rec = {"job": jid, "goal": pre["goal"], "repo": repo,
           "spec_path": pre.get("spec_path", ""), "spec_copy": str(frozen),
           "identity_file": str(ident_file),
           "members": list(pre["members"]), "resources": res,
           "head": pre.get("head", ""), "submitted_at": time.time(),
           "log": str(log_path(jid))}
    _write_private(record_path(jid), json.dumps(rec, indent=1))
    return rec


def _squeue(runner: Callable | None = None) -> tuple:
    """({job id: (state, node, time left, reason)}, fresh) for the user's
    swarm jobs. ``fresh`` is False when the throttle answered from its
    cache: an answer from before a submit never declares that job ended."""
    from delfin import scheduler_client

    user = os.environ.get("USER") or os.environ.get("LOGNAME") or ""
    argv = [*_SQUEUE_HEAD, user, "-n", JOB_NAME, "-o", "%i|%T|%N|%L|%R"]
    kind = scheduler_client.kind_of_argv(argv)
    out = scheduler_client.query(kind, argv, runner=runner)
    fresh = out is not scheduler_client.THROTTLED and out is not None
    if out is scheduler_client.THROTTLED:
        out = scheduler_client.default_client.last_result(kind, argv)
        if out is scheduler_client.THROTTLED:
            out = None
    found = {}
    for line in str(out or "").splitlines():
        parts = line.strip().split("|")
        if len(parts) >= 4 and parts[0]:
            found[parts[0]] = (parts[1], parts[2], parts[3],
                               parts[4] if len(parts) > 4 else "")
    return found, fresh


def jobs(runner: Callable | None = None) -> list:
    """Records of swarm jobs the scheduler still knows, oldest first, each
    with ``state``, ``node``, ``left`` and ``reason`` filled in. A job is
    dropped only on a fresh answer that does not list it, and never within
    :data:`SUBMIT_GRACE_S` of its submit."""
    recs = [r for r in records() if not r.get("ended_at")]
    if not recs:
        return []  # nothing submitted from here: the scheduler is not asked
    live, fresh = _squeue(runner)
    out = []
    now = time.time()
    for r in recs:
        job = str(r["job"])
        if job in live:
            r["state"], r["node"], r["left"], r["reason"] = live[job]
            out.append(r)
        elif not fresh or now - float(r.get("submitted_at") or 0) < SUBMIT_GRACE_S:
            r["state"], r["node"], r["left"], r["reason"] = "UNKNOWN", "", "-", ""
            out.append(r)
        else:
            _mark_ended(r, now)
    return out


def _mark_ended(rec: dict, now: float) -> None:
    """Keep the record as history, but ask the scheduler about it no more."""
    keep = {k: v for k, v in rec.items() if k != "progress"}
    keep["ended_at"] = now
    try:
        _write_private(record_path(str(rec["job"])), json.dumps(keep, indent=1))
    except OSError:
        pass


def cancel(job_ids: list, *, run: Callable | None = None) -> str:
    """Stop swarm jobs. Returns what happened, in a sentence. Only jobs
    with a record here are stopped (never another job of the user's).
    *run(cmd, env)* replaces the real scancel (tests): it returns the
    output, "" on success, None when it failed."""
    ids = [str(j) for j in job_ids if str(j).isdigit()]
    if not ids:
        return "No swarm job to stop."
    known = {str(r["job"]) for r in records()}
    foreign = [j for j in ids if j not in known]
    if foreign:
        return (f"Job {', '.join(foreign)} is not a swarm job started from "
                "here; nothing was stopped.")
    cmd = [_SCANCEL, *ids]
    if run is None:
        from delfin import scheduler_client
        out = scheduler_client.query(
            scheduler_client.kind_of_argv(cmd), cmd,
            runner=scheduler_client.child_env_runner(dict(os.environ)))
        if out is scheduler_client.THROTTLED:
            return ("The scheduler was asked to stop a job moments ago; try "
                    "again in half a minute.")
        ok = out is not None and not str(out).strip()
    else:
        out = run(cmd, dict(os.environ))
        ok = out is not None and not str(out).strip()
    if not ok:
        return f"scancel {' '.join(ids)} failed; is the job still running?"
    return (f"Stopped swarm job {', '.join(ids)}. Its member sessions end; "
            "their worktrees and branches stay.")


def status_line(rec: dict) -> str:
    """One line on a swarm job for a person."""
    state = str(rec.get("state") or "UNKNOWN")
    prog = rec.get("progress") or {}
    n = len(rec.get("members") or [])
    head = f"Swarm job {rec.get('job')} ({rec.get('goal', '')})"
    if state == "PENDING":
        why = rec.get("reason") or ""
        return (f"{head}: waiting for a free compute node"
                + (f" ({why})" if why and why != "None" else "") + ".")
    if state == "RUNNING":
        started = len(prog.get("started") or [])
        what = (f"{started} of {n} member sessions started"
                if started < n else f"all {n} member sessions started")
        if prog.get("phase") == "worktrees":
            what = "creating the member worktrees"
        return (f"{head}: running on {rec.get('node') or '?'}, {what}, "
                f"{rec.get('left') or '?'} left.")
    if state == "UNKNOWN":
        return f"{head}: submitted, the scheduler has not listed it yet."
    return f"{head}: {state.lower()}."


def format_jobs(rows: list) -> str:
    if not rows:
        return ("No swarm job runs. Start one with: delfin-agent swarm start "
                "[swarm.yaml]")
    out = [status_line(r) for r in rows]
    for r in rows:
        attach = (r.get("progress") or {}).get("attach")
        if attach:
            out.append(f"  job {r['job']} terminals: {attach}")
    out.append("Stop one: delfin-agent swarm stop JOB")
    return "\n".join(out)


# --------------------------------------------------------------------------
# The job body (on the compute node)
# --------------------------------------------------------------------------

def _tmux(*args: str, env: dict | None = None):
    import subprocess
    return subprocess.run(["tmux", *args], capture_output=True, text=True,
                          env=env, check=False)


def _mem_available_gb() -> float:
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) / 1024 / 1024
    except (OSError, ValueError, IndexError):
        pass
    return float("inf")


def member_launch_script(member: dict, repo: Path, prompt_file: Path,
                         python: str | None = None) -> str:
    """The bash script one member session runs in its tmux window."""
    python = python or sys.executable
    wt = repo / member["worktree"]
    argv = [python, "-m", "delfin.agent.cli", "chat", "-n", member["name"],
            "--permission-mode", member["permission_mode"]]
    return ("#!/bin/bash\n"
            f"cd {shlex.quote(str(wt))} || exit 1\n"
            f"exec {' '.join(shlex.quote(a) for a in argv)} "
            f"\"$(cat {shlex.quote(str(prompt_file))})\"\n")


def run_id_for(job: str) -> str:
    """The swarm run id of the members of swarm job *job*."""
    return f"swarm-{job}"


def run_in_job(spec_file: Path | str, repo: Path | str, *,
               identity_file: Path | str | None = None,
               env: dict | None = None,
               tmux: Callable = _tmux,
               sleep: Callable[[float], None] = time.sleep,
               start_members: Callable | None = None,
               poll_s: float = 60.0,
               mem_gb: Callable[[], float] = _mem_available_gb,
               say: Callable[[str], None] = lambda t: print(t, flush=True),
               ) -> int:
    """The body of a swarm job: worktrees, member sessions, then hold.

    Refuses outside a SLURM job: the members must run on a compute node.
    The member sessions are tied to this process (it is their lifeline),
    so they end with the job. Each member gets the git identity recorded
    at submit (*identity_file*; else the repository's git config, and no
    start without one) and the run id ``swarm-<job>``, which
    ``swarm status`` uses to tell this run's members from earlier ones."""
    env = os.environ if env is None else env
    job = str(env.get("SLURM_JOB_ID") or "").strip()
    if not job:
        say("This is the body of a swarm job and runs only inside one. "
            "Start a swarm with: delfin-agent swarm start [swarm.yaml]")
        return 2
    from . import launch as _launch
    from .spec import load_spec
    from delfin.agent import lifeline as _lifeline

    repo = Path(repo)
    os.chdir(repo)
    spec = load_spec(spec_file)
    plan = _launch.plan_launch(spec)
    members = plan["members"]
    if len(members) > MAX_MEMBERS_PER_JOB:
        say(f"Refused: {len(members)} member sessions, at most "
            f"{MAX_MEMBERS_PER_JOB} per swarm job.")
        return 2
    identity = _read_json(Path(identity_file)) if identity_file else {}
    git_env = identity_env(identity) or identity_env(git_identity(repo))
    if not git_env:
        say("Refused: " + _no_identity(repo))
        return 2
    import socket
    node = socket.gethostname().split(".")[0]
    run_id = run_id_for(job)
    progress = {"job": job, "node": node, "phase": "worktrees",
                "run": run_id, "started": [], "started_at": {},
                "at": time.time()}

    def _save() -> None:
        progress["at"] = time.time()
        _write_private(progress_path(job), json.dumps(progress, indent=1))

    _save()
    say(f"swarm job {job} on {node}: {plan['goal']} -- "
        f"{len(members)} member(s)")
    starter = start_members or _launch.start_members
    starter(spec, execute=True, repo=repo)

    # tmux needs the job's own TMPDIR for its socket (10.11.); the member
    # sessions live under this process, never under the login node's.
    tenv = dict(env)
    tenv["TMUX_TMPDIR"] = str(env.get("TMPDIR") or tenv.get("TMUX_TMPDIR")
                              or "/tmp")
    tenv.pop("TMUX", None)
    tenv.pop("TMUX_PANE", None)
    for key in (_lifeline.ENV_PID, _lifeline.ENV_TICKS):
        tenv.pop(key, None)
    os.environ.pop(_lifeline.ENV_PID, None)
    tenv.update(_lifeline.claim_root())
    tenv.update(git_env)
    tenv[RUN_ENV] = run_id
    session = f"delfin-swarm-{job}"
    if tmux("new-session", "-d", "-s", session, "-n", "hub",
            env=tenv).returncode != 0:
        say("tmux did not start; the member sessions cannot run.")
        progress["phase"] = "failed"
        _save()
        return 1
    tmux("set-option", "-t", session, "remain-on-exit", "on", env=tenv)
    sock = f"{tenv['TMUX_TMPDIR']}/tmux-{os.getuid()}/default"
    progress["attach"] = (f"srun --jobid={job} --overlap --pty tmux -S "
                          f"{sock} attach -t {session}")
    progress["phase"] = "members"
    _save()
    mdir = member_dir(job)
    use_script = shutil.which("script") is not None
    spacing = float(plan["start_spacing_s"] or 0)
    for i, m in enumerate(members):
        if i:
            sleep(spacing)
        if plan.get("memory_gate", True):
            waited = 0.0
            while mem_gb() < MEMORY_GATE_GB and waited < MEMORY_GATE_WAIT_S:
                say(f"{m['name']}: waiting for free memory")
                sleep(30.0)
                waited += 30.0
        prompt = mdir / f"{m['name']}.prompt"
        _write_private(prompt, m["prompt_text"])
        launch = mdir / f"{m['name']}.sh"
        _write_private(launch, member_launch_script(m, repo, prompt), 0o700)
        log = mdir / f"{m['name']}.log"
        cmd = (f"script -q -f -c {shlex.quote(str(launch))} "
               f"{shlex.quote(str(log))}" if use_script
               else shlex.quote(str(launch)))
        tmux("new-window", "-t", session, "-n", m["name"], cmd, env=tenv)
        progress["started"].append(m["name"])
        progress["started_at"][m["name"]] = time.time()
        _save()
        say(f"{m['name']}: started in {m['worktree']} (log {log})")
    progress["phase"] = "running"
    _save()
    names = {m["name"] for m in members}
    while True:
        sleep(poll_s)
        out = tmux("list-windows", "-t", session, "-F",
                   "#{window_name} #{pane_dead}", env=tenv)
        if out.returncode != 0:
            break
        alive = [line.split()[0] for line in out.stdout.splitlines()
                 if line.split() and line.split()[0] in names
                 and line.split()[-1] == "0"]
        if not alive:
            break
    tmux("kill-session", "-t", session, env=tenv)
    progress["phase"] = "ended"
    _save()
    say(f"swarm job {job}: every member session has ended")
    return 0


def main(argv: list | None = None) -> int:
    import argparse
    parser = argparse.ArgumentParser(
        prog="python -m delfin.agent.swarm.job",
        description="The body of a swarm job (submitted by "
                    "`delfin-agent swarm start`); not for direct use.")
    sub = parser.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="Run the swarm inside this SLURM job")
    r.add_argument("--spec", required=True)
    r.add_argument("--repo", required=True)
    r.add_argument("--identity", default="",
                   help="Private file with the git identity checked at submit")
    args = parser.parse_args(argv)
    return run_in_job(args.spec, args.repo,
                      identity_file=args.identity or None)


if __name__ == "__main__":
    raise SystemExit(main())
