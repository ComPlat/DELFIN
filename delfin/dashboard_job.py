"""``delfin-voila --job``: the dashboard runs in a SLURM job.

A cluster's login nodes are shared by everyone; the dashboard's server, its
kernels and the agent sessions inside them are real work and belong on a
compute node (issue #640). On a login node ``delfin-voila`` runs the
dashboard itself with ONE agent session; everything beyond that runs here:

    delfin-voila --job          start a dashboard job, or reconnect to the one running
    delfin-voila --job --keep   the same, every session it opens is kept
    delfin-voila --list       the dashboard jobs that run, with node and time left
    delfin-voila --attach [JOB]   reconnect (the newest one without JOB)
    delfin-voila --stop [JOB]     end a dashboard job
    delfin-voila --new        start another one although one runs

The job plays the part tmux played on a login node: it is not tied to a
terminal, so closing the terminal, the browser or the tunnel ends nothing.
The login node only keeps an ``ssh -L`` tunnel to the job's node (no load);
an editor that forwards ports (VS Code) forwards it on to the laptop as before.

State lives under ``~/.delfin/dashboard_jobs`` (0700): one record per job and
its access token (0600). The token reaches the job through the environment
(sbatch passes the whole environment, never ``--export`` with named
variables), so it appears in no command line and no shared log. The
scheduler is asked through ``scheduler_client`` (throttled, one query kind).
"""

from __future__ import annotations

import json
import os
import re
import secrets
import shlex
import shutil
import socket
import subprocess
import sys
import time
from collections.abc import Callable, Sequence
from pathlib import Path

#: The SLURM job name every dashboard job carries.
JOB_NAME = "delfin-dashboard"

#: How long a dashboard job asks for unless the user settings say otherwise
#: (``hpc.dashboard_time``). Longer than an agent job: it is the place the
#: kept sessions live.
DEFAULT_TIME = "24:00:00"

#: The line the launcher prints once Voila is about to serve.
_START_RE = re.compile(r"Starting DELFIN Dashboard on (https?)://[^:\s]+:(\d+)")

#: Flags that only make sense for the launcher on the login node, or that
#: the job sets itself.
_LOCAL_ONLY = {"--list", "--attach", "--stop", "--new", "--keep", "--job",
               "--open-browser", "--no-browser", "--strict-port"}
_LOCAL_WITH_VALUE = {"--port", "--token", "--ip"}


def state_dir() -> Path:
    """``~/.delfin/dashboard_jobs``, created 0700. Resolved per call."""
    d = Path.home() / ".delfin" / "dashboard_jobs"
    d.mkdir(parents=True, mode=0o700, exist_ok=True)
    try:
        os.chmod(d, 0o700)
    except OSError:
        pass
    return d


def _write_private(path: Path, text: str) -> None:
    tmp = path.with_name(path.name + ".tmp")
    fd = os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
    with os.fdopen(fd, "w", encoding="utf-8") as fh:
        fh.write(text)
    os.chmod(tmp, 0o600)
    tmp.replace(path)


def job_spec(settings: dict | None = None) -> dict:
    """Partition, cpus and memory as the agent jobs (``hpc.agent_*``); time
    from ``hpc.dashboard_time`` or :data:`DEFAULT_TIME`."""
    from delfin import hpc_guard

    if settings is None:
        settings = hpc_guard._read_user_settings()
    job = hpc_guard.job_settings(settings)
    hpc = settings.get("hpc") if isinstance(settings, dict) else None
    hpc = hpc if isinstance(hpc, dict) else {}
    job["time"] = str(hpc.get("dashboard_time") or DEFAULT_TIME)
    return job


def passthrough_args(argv: Sequence[str]) -> list:
    """The launcher arguments the dashboard inside the job should get."""
    out: list = []
    skip = False
    args = list(argv)
    for i, a in enumerate(args):
        if skip:
            skip = False
            continue
        name = a.split("=", 1)[0]
        if name in _LOCAL_WITH_VALUE:
            skip = "=" not in a
            continue
        if name in ("--attach", "--stop"):
            # optional value: drop a following job id too
            if "=" not in a and i + 1 < len(args) and args[i + 1].isdigit():
                skip = True
            continue
        if name in _LOCAL_ONLY:
            continue
        out.append(a)
    return out


# --------------------------------------------------------------------------
# Records
# --------------------------------------------------------------------------

def _record_path(job: str) -> Path:
    return state_dir() / f"{job}.json"


def token_path(job: str) -> Path:
    return state_dir() / f"{job}.token"


def log_path(job: str) -> Path:
    return state_dir() / f"{job}.log"


def records() -> list:
    out = []
    for p in sorted(state_dir().glob("*.json")):
        try:
            r = json.loads(p.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        if isinstance(r, dict) and r.get("job"):
            out.append(r)
    out.sort(key=lambda r: float(r.get("submitted_at") or 0))
    return out


def _forget(job: str) -> None:
    """Drop the token of a job that has ended (the record stays as history)."""
    try:
        token_path(job).unlink()
    except OSError:
        pass


# --------------------------------------------------------------------------
# Scheduler
# --------------------------------------------------------------------------

#: Scheduler argv heads as module constants (the direct-call guard scans
#: call sites); every run goes through ``delfin.scheduler_client``.
_SQUEUE_HEAD = ("squeue", "-h", "-u")
_SBATCH = "sbatch"

#: A job the scheduler does not list yet within this long after its submit
#: is still starting, not gone (squeue lags behind sbatch).
SUBMIT_GRACE_S = 120.0


def _squeue(runner: Callable | None = None) -> tuple:
    """({job id: (state, node, time left)}, fresh) for the user's dashboard
    jobs. ``fresh`` is False when the throttle answered from its cache: an
    answer from before a submit must never declare that job ended."""
    from delfin import scheduler_client

    user = os.environ.get("USER") or os.environ.get("LOGNAME") or ""
    argv = [*_SQUEUE_HEAD, user, "-n", JOB_NAME, "-o", "%i|%T|%N|%L"]
    out = scheduler_client.query("squeue", argv, runner=runner)
    fresh = out is not scheduler_client.THROTTLED and out is not None
    if out is scheduler_client.THROTTLED:
        out = scheduler_client.default_client.last_result("squeue", argv)
    jobs = {}
    for line in str(out or "").splitlines():
        parts = line.strip().split("|")
        if len(parts) == 4 and parts[0]:
            jobs[parts[0]] = (parts[1], parts[2], parts[3])
    return jobs, fresh


def running(runner: Callable | None = None) -> list:
    """Records of dashboard jobs the scheduler still knows, newest last,
    each with ``state``, ``node`` and ``left`` filled in. A job is only
    dropped (and its token with it) on a fresh answer that does not list
    it, and never within :data:`SUBMIT_GRACE_S` of its submit."""
    live, fresh = _squeue(runner)
    out = []
    now = time.time()
    for r in records():
        job = str(r["job"])
        if job in live:
            r["state"], r["node"], r["left"] = live[job]
            out.append(r)
        elif (not fresh
              or now - float(r.get("submitted_at") or 0) < SUBMIT_GRACE_S):
            r["state"], r["node"], r["left"] = "UNKNOWN", "", "-"
            out.append(r)
        else:
            _forget(job)
    return out


def submit(argv: Sequence[str], *, keep: bool, port: int,
           run: Callable | None = None, cwd: str | None = None) -> str:
    """Submit a dashboard job. Returns the job id."""
    job = job_spec()
    token = secrets.token_urlsafe(32)
    cwd = cwd or os.getcwd()
    inner = [sys.executable, "-m", "delfin.cli_voila", "--no-browser",
             f"--port={int(port)}", *passthrough_args(argv)]
    if keep:
        inner.append("--keep")
    pending = state_dir() / f"pending-{os.getpid()}"
    env = dict(os.environ)
    env["DELFIN_VOILA_TOKEN"] = token
    env["DELFIN_LAUNCH_CWD"] = cwd
    cmd = [_SBATCH, "--parsable", f"--job-name={JOB_NAME}",
           f"--partition={job['partition']}", f"--time={job['time']}",
           "--nodes=1", "--ntasks=1", f"--cpus-per-task={int(job['cpus'])}",
           f"--mem={job['mem']}", f"--output={pending}.log",
           f"--chdir={cwd}",
           "--wrap=exec " + " ".join(shlex.quote(a) for a in inner)]
    if run is None:
        from delfin import scheduler_client
        out = scheduler_client.query(
            scheduler_client.kind_of_argv(cmd), cmd,
            runner=scheduler_client.child_env_runner(env))
        out = "" if out is scheduler_client.THROTTLED else str(out or "")
    else:
        out = run(cmd, env)
    jid = out.strip().split(";")[0]
    if not jid.isdigit():
        raise RuntimeError(f"sbatch gave no job id: {out.strip() or 'no output'}")
    _write_private(token_path(jid), token)
    _write_private(_record_path(jid), json.dumps({
        "job": jid, "keep": bool(keep), "submitted_at": time.time(),
        "cwd": cwd, "log": f"{pending}.log", "spec": job,
        "args": passthrough_args(argv)}, indent=1))
    return jid


def cancel(job: str, run: Callable = subprocess.run) -> bool:
    proc = run(["scancel", str(job)], capture_output=True, text=True)
    if proc.returncode == 0:
        _forget(str(job))
        return True
    return False


# --------------------------------------------------------------------------
# Reaching the job
# --------------------------------------------------------------------------

def served_port(rec: dict) -> tuple:
    """(scheme, port) the dashboard in the job serves on, or ("", 0) while
    it has not started yet."""
    try:
        text = Path(rec.get("log") or log_path(rec["job"])).read_text(
            encoding="utf-8", errors="replace")
    except OSError:
        return "", 0
    m = _START_RE.search(text)
    return (m.group(1), int(m.group(2))) if m else ("", 0)


def dashboard_path(node: str) -> str:
    """The render path from the note the dashboard leaves, if it is the
    note of the job's node; "/" otherwise (Voila redirects)."""
    try:
        from delfin.agent import where
        note = where.dashboard()
    except Exception:
        return "/"
    host = str(note.get("host") or "").split(".")[0]
    if node and host == node.split(".")[0] and note.get("url"):
        return str(note["url"])
    return "/"


def free_local_port(start: int) -> int:
    for port in range(start, start + 100):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            try:
                s.bind(("127.0.0.1", port))
            except OSError:
                continue
            return port
    raise RuntimeError(f"no free local port in {start}-{start + 99}")


def node_reachable(node: str, run: Callable = subprocess.run) -> bool:
    try:
        proc = run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10",
                    "-o", "StrictHostKeyChecking=accept-new", node, "true"],
                   capture_output=True, text=True, timeout=30)
    except Exception:
        return False
    return proc.returncode == 0


def enable_node_ssh(*, home: Path | None = None,
                    run: Callable = subprocess.run) -> str:
    """Let the user's own key open their own job nodes: generate a key when
    there is none, add its public half to ``authorized_keys`` (appended,
    nothing removed). Returns what was done."""
    ssh = (home or Path.home()) / ".ssh"
    ssh.mkdir(mode=0o700, exist_ok=True)
    key = ssh / "id_ed25519"
    done = []
    if not key.exists():
        run(["ssh-keygen", "-q", "-t", "ed25519", "-N", "", "-f", str(key),
             "-C", "delfin-node-access"], check=True)
        done.append(f"created {key}")
    pub = Path(str(key) + ".pub").read_text(encoding="utf-8").strip()
    auth = ssh / "authorized_keys"
    have = auth.read_text(encoding="utf-8") if auth.exists() else ""
    if pub not in have:
        with open(auth, "a", encoding="utf-8") as fh:
            if have and not have.endswith("\n"):
                fh.write("\n")
            fh.write(pub + "\n")
        done.append(f"added your public key to {auth}")
    os.chmod(auth, 0o600)
    return "; ".join(done) or "already set up"


def _ask(question: str) -> bool:
    if not sys.stdin.isatty():
        return False
    try:
        answer = input(f"{question} [Y/n] ").strip().lower()
    except EOFError:
        return False
    return answer in ("", "y", "yes", "j", "ja")


def _say(msg: str = "") -> None:
    print(msg, flush=True)


def wait_until_served(rec: dict, *, poll: float = 5.0, limit: float = 7200.0,
                      runner: Callable | None = None) -> dict | None:
    """Wait for the job to run and its dashboard to print its port."""
    t0 = time.time()
    shown = ""
    while time.time() - t0 < limit:
        live = {r["job"]: r for r in running(runner)}
        cur = live.get(rec["job"])
        if cur is None:
            _say(f"\nDashboard job {rec['job']} has ended -- see "
                 f"{rec.get('log')}")
            return None
        if cur["state"] != shown and cur["state"] != "UNKNOWN":
            _say(f"  job {cur['job']}: {cur['state'].lower()}"
                 + (f" on {cur['node']}" if cur["state"] == "RUNNING" else
                    " (waiting for a free node)"))
            shown = cur["state"]
        if cur["state"] == "RUNNING" and served_port(cur)[1]:
            return cur
        time.sleep(poll)
    return None


def connect(rec: dict, *, port: int, open_browser: bool | None = None,
            runner: Callable | None = None) -> int:
    """Wait for *rec*'s dashboard, tunnel to it, print the link. Blocks
    until Ctrl+C; the dashboard keeps running after that."""
    cur = wait_until_served(rec, runner=runner)
    if cur is None:
        return 1
    node = cur["node"]
    scheme, remote = served_port(cur)
    if not node_reachable(node):
        _say("\nThe tunnel needs ssh from this login node to the job's node "
             f"({node}), and it was refused.")
        if _ask("Allow your own ssh key to open your own job nodes "
                "(adds it to ~/.ssh/authorized_keys)?"):
            _say("  " + enable_node_ssh())
        if not node_reachable(node):
            _say("Still no ssh to the node. The dashboard runs in job "
                 f"{cur['job']}; ask your cluster team whether ssh to nodes "
                 "of your own running jobs is allowed.")
            return 1
    local = free_local_port(port)
    try:
        token = token_path(cur["job"]).read_text(encoding="utf-8").strip()
    except OSError:
        token = ""
    path = dashboard_path(node)
    url = f"{scheme}://localhost:{local}{path}"
    tunnel = subprocess.Popen(
        ["ssh", "-N", "-o", "ExitOnForwardFailure=yes",
         "-o", "ServerAliveInterval=30", "-o", "BatchMode=yes",
         "-L", f"{local}:localhost:{remote}", node])
    time.sleep(1.5)
    if tunnel.poll() is not None:
        _say(f"The tunnel to {node} did not come up (ssh exit {tunnel.returncode}).")
        return 1
    _say("")
    _say(f"DELFIN dashboard (job {cur['job']} on {node}, {cur['left']} left"
         + (", sessions kept" if cur.get("keep") else "") + ")")
    if sys.stdout.isatty() and token:
        _say(f"  Open: {url}?token={token}")
    else:
        _say(f"  Open: {url}  (token: {token_path(cur['job'])})")
    _say("  Ctrl+C closes only this tunnel -- the dashboard keeps running.")
    _say("  Back in later:  delfin-voila --attach     All of them:  delfin-voila --list")
    _say("  End it:         delfin-voila --stop")
    if open_browser and token:
        opener = os.environ.get("BROWSER")
        if opener:
            try:
                subprocess.Popen([opener, f"{url}?token={token}"],
                                 stdout=subprocess.DEVNULL,
                                 stderr=subprocess.DEVNULL)
            except OSError:
                pass
    try:
        return int(tunnel.wait() or 0)
    except KeyboardInterrupt:
        tunnel.terminate()
        _say(f"\nTunnel closed. The dashboard keeps running in job {cur['job']}"
             " -- back in with: delfin-voila --attach")
        return 0


# --------------------------------------------------------------------------
# Command
# --------------------------------------------------------------------------

def format_list(rows: list) -> str:
    if not rows:
        return "No dashboard job runs. Start one with: delfin-voila --job  (or --job --keep)"
    out = ["JOB        STATE      NODE         LEFT         KEPT  STARTED"]
    for r in rows:
        started = time.strftime("%d.%m. %H:%M",
                                time.localtime(float(r.get("submitted_at") or 0)))
        out.append(f"{r['job']:<10} {r['state']:<10} {r['node'] or '-':<12} "
                   f"{r['left']:<12} {'yes' if r.get('keep') else 'no':<5} "
                   f"{started}")
    out.append("Reconnect: delfin-voila --attach [JOB]    End: delfin-voila --stop [JOB]")
    return "\n".join(out)


def _pick(rows: list, job: str) -> dict | None:
    if job and job not in ("latest", "all"):
        return next((r for r in rows if r["job"] == job), None)
    return rows[-1] if rows else None


def main(args, argv: Sequence[str]) -> int:
    """Run ``delfin-voila`` from a login node (or ``--list``/``--attach``/
    ``--stop`` from anywhere)."""
    if not shutil.which("sbatch"):
        _say("No SLURM here (sbatch not found).")
        return 1
    rows = running()
    if args.list:
        _say(format_list(rows))
        return 0
    if args.stop is not None:
        if args.stop == "all":
            targets = rows
        else:
            one = _pick(rows, args.stop)
            if args.stop == "latest" and len(rows) > 1:
                _say(format_list(rows))
                _say("Several dashboards run -- name one: delfin-voila --stop JOB"
                     "  (or --stop all)")
                return 1
            targets = [one] if one else []
        if not targets:
            _say("No such dashboard job.\n" + format_list(rows))
            return 1
        for r in targets:
            _say(f"job {r['job']}: " + ("ended" if cancel(r["job"]) else "scancel failed"))
        return 0
    open_browser = bool(getattr(args, "open_browser", None))
    if args.attach is not None:
        rec = _pick(rows, args.attach)
        if rec is None:
            _say("No such dashboard job.\n" + format_list(rows))
            return 1
        return connect(rec, port=args.port, open_browser=open_browser)
    if rows and not args.new:
        rec = rows[-1]
        _say(f"A dashboard already runs (job {rec['job']}) -- reconnecting. "
             "Another one: delfin-voila --job --new")
        if args.keep and not rec.get("keep"):
            _say("  (it was started without --keep; arm 'Keep session' in its "
                 "header, or start a kept one with --job --new --keep)")
        return connect(rec, port=args.port, open_browser=open_browser)
    try:
        jid = submit(argv, keep=bool(args.keep), port=args.port)
    except Exception as exc:
        _say(f"Could not submit the dashboard job: {exc}")
        return 1
    spec = job_spec()
    _say(f"The dashboard runs as SLURM job {jid} "
         f"({spec['partition']}, {spec['cpus']} cores, {spec['mem']}, "
         f"{spec['time']}){', every session kept' if args.keep else ''}.")
    rec = next((r for r in records() if r["job"] == jid), {"job": jid})
    return connect(rec, port=args.port, open_browser=open_browser)
