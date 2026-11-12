"""Keep DELFIN's heavy parts off shared HPC login nodes.

A cluster login node is shared by every user of the machine. On 2026-11-10
the cluster operators suspended an account because DELFIN agent sessions,
the agent service, dashboards and test runs filled a login node with tasks
(GitHub issue #640). The rule since then is enforced here, not left to
whoever starts something:

* agent sessions, the agent service, the dashboard server, the MCP servers
  and test runs do not start on a login node; inside a SLURM job they run
  as before, and on a workstation (no scheduler) nothing changes;
* an interactive session asked for on a login node is started on a compute
  node instead, through ``srun --pty`` (:func:`run_on_compute_node`);
* :func:`watch` is the net under all of it: it samples this user's own
  processes on the login node and stops (SIGSTOP, by PID) a DELFIN process
  that keeps the user above the core budget. It is not started by hand:
  every ``delfin-agent`` call on a login node makes sure one watcher runs
  on that host (:func:`ensure_watcher`, one per user and host).

Where a host is a login node is decided in one place,
``delfin.agent.calc_policy.host_kind`` (no job id, scheduler commands on
PATH). There is deliberately no environment switch that turns the refusal
off: a session must not be able to talk its way back onto the login node.
"""
from __future__ import annotations

import json
import os
import shutil
import signal
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Optional, Sequence

#: Cores this user's processes may keep busy on a login node, together.
LOGIN_CORE_BUDGET = 2.0

#: Where the compute-node session goes when the settings name nothing. The
#: dev partitions cap the time at 30 minutes, too short for a session.
DEFAULT_JOB = {
    "partition": "cpu",
    "time": "08:00:00",
    "cpus": 4,
    "mem": "16G",
}

#: The srun binary, named once (it is not one of the throttled query
#: commands of delfin.scheduler_client; it starts a job step).
_SRUN = "srun"


class LoginNodeRefused(RuntimeError):
    """A heavy part of DELFIN was asked to start on a shared login node."""


def on_login_node(env: Optional[dict] = None) -> bool:
    """True on a cluster login node: no SLURM job, scheduler on PATH."""
    from delfin.agent.calc_policy import host_kind
    return host_kind(env) == "login_node"


def refusal_text(what: str, how: str = "") -> str:
    """The one message every refusal prints: what, why, what to do."""
    how = how or ("start it inside a SLURM job (sbatch, or srun --pty for "
                  "an interactive shell on a compute node)")
    return (f"DELFIN does not run {what} on a login node. Login nodes are "
            "shared by every user of the cluster; heavy work there breaks "
            f"the cluster's rules and gets the account suspended. Instead: "
            f"{how}.")


#: What runs on a login node after all (the person's decision, 11.11.:
#: users have to get in again). An interactive session, its agent service
#: and the small stdio servers it starts are light -- they mostly wait for
#: the model -- and what they would RUN that is heavy (tests, builds,
#: installs) the bash gate refuses there and sends to SLURM
#: (calc_policy.login_node_heavy). Swarms, waves, the schedule daemon,
#: benchmarks and the pre-PR check stay off the login node.
#: DELFIN_LOGIN_NODE_STRICT=1 refuses these too.
LOGIN_NODE_INTERACTIVE = frozenset({
    "the agent service",
    "the DELFIN tools server",
    "the DELFIN ops server",
    "the DELFIN docs server",
})


def refuse_on_login_node(what: str, how: str = "", *,
                         env: Optional[dict] = None) -> None:
    """Raise :class:`LoginNodeRefused` when this is a login node -- except for
    the interactive pieces in :data:`LOGIN_NODE_INTERACTIVE`."""
    e = os.environ if env is None else env
    strict = str(e.get("DELFIN_LOGIN_NODE_STRICT", "")).strip().lower() in (
        "1", "true", "yes", "on")
    if what in LOGIN_NODE_INTERACTIVE and not strict:
        return
    if on_login_node(env):
        raise LoginNodeRefused(refusal_text(what, how))


def exit_on_login_node(what: str, how: str = "", *,
                       env: Optional[dict] = None) -> None:
    """For a command-line entry point: print the refusal and exit 2."""
    try:
        refuse_on_login_node(what, how, env=env)
    except LoginNodeRefused as exc:
        print(str(exc), file=sys.stderr)
        raise SystemExit(2) from None


# --------------------------------------------------------------------------
# An interactive session on a compute node
# --------------------------------------------------------------------------

def job_settings(settings: Optional[dict] = None) -> dict:
    """The compute-node job for a session: DEFAULT_JOB under the ``hpc``
    entry of the user settings (``agent_partition``, ``agent_time``,
    ``agent_cpus``, ``agent_mem``). Read only; a broken file gives the
    defaults."""
    if settings is None:
        settings = _read_user_settings()
    hpc = settings.get("hpc") if isinstance(settings, dict) else None
    hpc = hpc if isinstance(hpc, dict) else {}
    job = dict(DEFAULT_JOB)
    for key in DEFAULT_JOB:
        value = hpc.get(f"agent_{key}")
        if value not in (None, ""):
            job[key] = value
    return job


def _read_user_settings() -> dict:
    try:
        from delfin.user_settings import get_settings_path
        data = json.loads(Path(get_settings_path()).read_text(encoding="utf-8"))
    except Exception:
        return {}
    return data if isinstance(data, dict) else {}


def compute_node_argv(argv: Sequence[str], *, job: Optional[dict] = None,
                      name: str = "delfin-agent") -> list:
    """The srun command that runs *argv* interactively on a compute node.

    The environment goes along as srun's default does (all of it); there is
    no ``--export`` with named variables (a known cause of failed job
    starts on this cluster)."""
    job = job_settings() if job is None else job
    return [_SRUN, "--pty",
            f"--partition={job['partition']}",
            f"--time={job['time']}",
            "--nodes=1", "--ntasks=1",
            f"--cpus-per-task={int(job['cpus'])}",
            f"--mem={job['mem']}",
            f"--job-name={name}",
            *argv]


def run_on_compute_node(argv: Sequence[str], what: str, *,
                        env: Optional[dict] = None,
                        execvp: Callable = os.execvp,
                        isatty: Callable[[], bool] = lambda: sys.stdin.isatty(),
                        ) -> None:
    """On a login node, replace this process with *argv* under ``srun --pty``.

    Returns (doing nothing) off a login node. A non-interactive caller (no
    terminal: a script, a cron job, an agent) is refused instead -- it
    belongs in an sbatch job, not in a blocking interactive step."""
    if not on_login_node(env):
        return
    if not isatty() or not shutil.which(_SRUN):
        exit_on_login_node(what, "submit it with sbatch", env=env)
        return
    cmd = compute_node_argv(list(argv))
    print(f"{what}: login node -- starting it on a compute node "
          f"({' '.join(cmd[:8])} ...); it waits for a free node.",
          file=sys.stderr)
    execvp(cmd[0], cmd)


#: delfin-agent subcommands (by handler) that run agent turns or tests and
#: move onto a compute node when started from a terminal on a login node.
#: An interactive session (chat, run, the TUI) stays on the login node
#: since 11.11. (LOGIN_NODE_INTERACTIVE); what it runs heavy goes to SLURM
#: through the bash gate.
MOVED_COMMANDS = {
    "_run_prcheck_subcommand": "the pre-PR check (it runs tests)",
    "cmd_bench": "a benchmark run",
}


def keep_cli_off_login_node(func, args, argv: Sequence[str], *,
                            env: Optional[dict] = None,
                            execvp: Callable = os.execvp,
                            isatty: Optional[Callable[[], bool]] = None) -> None:
    """The delfin-agent entry point's guard, called once its arguments are
    parsed. A session, a turn, the TUI, a pre-PR check or a benchmark moves
    onto a compute node (srun --pty) from a terminal and is refused without
    one. Starting the schedule daemon is refused: it outlives the command,
    so it belongs in an sbatch job. Starting a swarm is handed over to the
    swarm command, which submits it as a SLURM job itself
    (delfin.agent.swarm.job; its members never start on a login node:
    launch.start_members refuses there). Everything else (status, lists,
    approvals, watch, where, ...) stays: it is light."""
    try:
        ensure_watcher(env=env)
    except Exception:
        pass  # the watcher is a net; it never breaks the command
    name = getattr(func, "__name__", "")
    # `delfin-agent swarm start` on a login node: no refusal here -- the
    # swarm command checks the spec and submits the swarm as a SLURM job
    # (delfin.agent.swarm.job.route says "job" on every login node).
    if name == "cmd_scheduler" and getattr(args, "scheduler_action", "") == "start":
        exit_on_login_node("the schedule daemon",
                           "run it inside an sbatch job", env=env)
    what = MOVED_COMMANDS.get(name)
    if what:
        kw = {} if isatty is None else {"isatty": isatty}
        run_on_compute_node([sys.executable, "-m", "delfin.agent.cli", *argv],
                            what, env=env, execvp=execvp, **kw)


# --------------------------------------------------------------------------
# The load watcher
# --------------------------------------------------------------------------

@dataclass
class Proc:
    pid: int
    cpu_s: float      # utime + stime so far
    cmdline: str


def own_processes(uid: Optional[int] = None,
                  proc: Path = Path("/proc")) -> list:
    """This user's processes with their CPU seconds so far (from /proc)."""
    uid = os.getuid() if uid is None else uid
    tick = os.sysconf("SC_CLK_TCK") if hasattr(os, "sysconf") else 100
    out = []
    for entry in proc.iterdir():
        if not entry.name.isdigit():
            continue
        try:
            if entry.stat().st_uid != uid:
                continue
            fields = (entry / "stat").read_text().rsplit(")", 1)[1].split()
            cpu = (int(fields[11]) + int(fields[12])) / tick
            cmd = (entry / "cmdline").read_bytes().replace(b"\0", b" ")
        except (OSError, IndexError, ValueError):
            continue
        out.append(Proc(int(entry.name), cpu, cmd.decode(errors="replace").strip()))
    return out


def busy_cores(before: Iterable[Proc], after: Iterable[Proc],
               seconds: float) -> dict:
    """Cores each pid kept busy between two samples, by pid."""
    first = {p.pid: p.cpu_s for p in before}
    return {p.pid: (p.cpu_s - first[p.pid]) / seconds
            for p in after if p.pid in first and seconds > 0}


def is_delfin(cmdline: str) -> bool:
    """A process DELFIN started (its interpreter runs delfin code or a
    delfin entry point), never the person's shell or editor."""
    return "delfin" in cmdline and "hpc_guard" not in cmdline


def watch(*, budget: float = LOGIN_CORE_BUDGET, interval: float = 30.0,
          strikes: int = 2, rounds: Optional[int] = None,
          sample: Callable[[], list] = own_processes,
          stop: Callable[[int], None] = lambda pid: os.kill(pid, signal.SIGSTOP),
          note: Callable[[str], None] = lambda text: print(text, flush=True),
          sleep: Callable[[float], None] = time.sleep,
          env: Optional[dict] = None) -> int:
    """Sample this user's load on the login node; stop what keeps it high.

    Every *interval* seconds the cores this user's processes keep busy are
    added up. Over *budget* for *strikes* samples in a row, the busiest
    DELFIN process is stopped with SIGSTOP (by its PID; ``kill -CONT <pid>``
    resumes it inside a job) and the stop is reported through *note*.
    Returns the number of processes stopped. Off a login node it does
    nothing."""
    if not on_login_node(env):
        note("hpc_guard watch: not a login node, nothing to watch")
        return 0
    stopped = 0
    over = 0
    prev = sample()
    done = 0
    while rounds is None or done < rounds:
        sleep(interval)
        cur = sample()
        load = busy_cores(prev, cur, interval)
        prev = cur
        done += 1
        total = sum(load.values())
        if total <= budget:
            over = 0
            continue
        over += 1
        if over < strikes:
            continue
        names = {p.pid: p.cmdline for p in cur}
        candidates = sorted((pid for pid in load if is_delfin(names.get(pid, ""))),
                            key=lambda pid: load[pid], reverse=True)
        if not candidates:
            note(f"hpc_guard: {total:.1f} cores busy on the login node, none "
                 "of it DELFIN -- look with htop")
            over = 0
            continue
        pid = candidates[0]
        try:
            stop(pid)
        except OSError as exc:
            note(f"hpc_guard: could not stop {pid}: {exc}")
            continue
        stopped += 1
        over = 0
        note(f"hpc_guard: {total:.1f} cores busy on the login node (budget "
             f"{budget:g}); stopped pid {pid} ({load[pid]:.1f} cores): "
             f"{names.get(pid, '')[:160]}")
    return stopped


# --------------------------------------------------------------------------
# One watcher per user and login node, started by DELFIN itself
# --------------------------------------------------------------------------

def _watch_files() -> tuple:
    """The lock and log of this host's watcher (HOME is shared between the
    login nodes, so the host name is part of the name)."""
    import socket
    base = Path.home() / ".delfin"
    host = socket.gethostname().split(".")[0] or "host"
    return base / f"hpc_watch.{host}.lock", base / f"hpc_watch.{host}.log"


def _try_lock(path: Path):
    """An open file holding an exclusive lock on *path*, or None when
    another process holds it."""
    import fcntl
    path.parent.mkdir(parents=True, exist_ok=True)
    fh = open(path, "a+")
    try:
        fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        fh.close()
        return None
    return fh


def watcher_running() -> bool:
    """True when this host's watcher holds its lock."""
    fh = _try_lock(_watch_files()[0])
    if fh is None:
        return True
    fh.close()
    return False


def _spawn_watcher() -> None:
    """Start the watcher detached (own session, no terminal), its output in
    the host's log file."""
    import subprocess
    log = _watch_files()[1]
    log.parent.mkdir(parents=True, exist_ok=True)
    with open(log, "a") as out:
        subprocess.Popen([sys.executable, "-m", "delfin.hpc_guard", "watch"],
                         stdin=subprocess.DEVNULL, stdout=out,
                         stderr=subprocess.STDOUT, close_fds=True,
                         start_new_session=True)


def ensure_watcher(*, env: Optional[dict] = None,
                   running: Optional[Callable[[], bool]] = None,
                   spawn: Optional[Callable[[], None]] = None) -> bool:
    """On a login node, start this host's watcher unless it runs already.
    Returns True when one was started. Off a login node nothing happens."""
    if not on_login_node(env):
        return False
    if spawn is None and "PYTEST_CURRENT_TEST" in os.environ:
        return False  # a test that simulates a login node starts no daemon
    running = running or watcher_running
    spawn = spawn or _spawn_watcher
    if running():
        return False
    spawn()
    return True


def main(argv: Optional[Sequence[str]] = None) -> int:
    import argparse
    parser = argparse.ArgumentParser(
        prog="python -m delfin.hpc_guard",
        description="Keep DELFIN off shared login nodes: check this host, or "
                    "watch this user's load and stop DELFIN processes that "
                    "keep it over the core budget.")
    sub = parser.add_subparsers(dest="cmd", required=True)
    sub.add_parser("check", help="Say whether this host is a login node")
    w = sub.add_parser("watch", help="Watch the load (until Ctrl+C)")
    w.add_argument("--budget", type=float, default=LOGIN_CORE_BUDGET)
    w.add_argument("--interval", type=float, default=30.0)
    args = parser.parse_args(argv)
    if args.cmd == "check":
        print("login node" if on_login_node() else "not a login node")
        return 0
    if not on_login_node():
        print("hpc_guard watch: not a login node, nothing to watch")
        return 0
    lock = _try_lock(_watch_files()[0])
    if lock is None:
        print("hpc_guard watch: a watcher already runs on this host")
        return 0
    try:
        watch(budget=args.budget, interval=args.interval)
    except KeyboardInterrupt:
        pass
    finally:
        lock.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
