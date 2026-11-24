"""A program DELFIN starts listens on the loopback only.

Input: the argv of a program DELFIN is about to start on this machine (an
ORCA run, a check of an input, a saddle search). Output: the argv to start
instead -- the same program under the socket guard's inbound rule
(``delfin/agent/landlock_exec.py --fs 0 --net open``), or the argv
unchanged where that rule does not apply or cannot hold.

Why. A parallel ORCA runs on MPI, and MPI listens for its ranks on every
interface of the host: measured 2026-11-24, OpenMPI's TCP transport on
0.0.0.0:1024-1025 and [::]:1024-1025, mpirun on four more wildcard ports,
for a two-process water single point. On a host without a firewall --
a shared server, a workstation, a cluster's login node -- those ports are
reachable from the network for as long as the run lasts. The interface
settings of OpenMPI and PMIx do not keep mpirun off the wildcard, and they
are one MPI's settings; the rule here holds at the system call, for any
program and any MPI. Under it every bind() and listen() of the program is
carried out by the guard's supervisor, on the loopback: the ranks of one
machine still reach each other (shared memory, the loopback), nothing
reaches them from elsewhere. Several parallel runs at once take the next
free ports, as they do on the wildcard.

When the argv stays as it is:

  * inside a SLURM allocation of several nodes the rule is ``--inbound
    any``: the job's ranks on the other nodes must connect, and the
    cluster network is the job's;
  * a program already under the guard (an agent's command running a DELFIN
    workflow; ``DELFIN_INBOUND_GUARD=held``, set by the guard) is held by
    it -- the kernel takes one notifying filter per process tree;
  * ``DELFIN_INBOUND_GUARD=off`` switches it off, for a site whose MPI
    needs another interface on one node;
  * on a host where the guard cannot run (no seccomp user notification, a
    ptrace policy that keeps a parent from a child's descriptors, macOS):
    a one-time probe per process finds out, and the program runs as
    before, rather than not at all.

Stopping is unchanged for the callers: the program keeps its process group
(the guard's supervisor forks it in the same group), a signal to the
supervisor is passed on to it, and it ends if the supervisor ends.
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

#: The guard's own entry point, by path: under ``python -I`` nothing in the
#: working directory or the environment changes what runs first.
_HELPER = Path(__file__).resolve().parent / "agent" / "landlock_exec.py"
_MARKER = "DELFIN_INBOUND_GUARD"
_usable: bool | None = None


def multi_node_allocation(env=None) -> bool:
    """Whether ``env`` (default: this process's) is inside a SLURM
    allocation of more than one node (``SLURM_JOB_NUM_NODES``, or
    ``SLURM_NNODES`` on older releases)."""
    env = os.environ if env is None else env
    if not env.get("SLURM_JOB_ID"):
        return False
    try:
        return int(env.get("SLURM_JOB_NUM_NODES") or env.get("SLURM_NNODES") or 1) > 1
    except ValueError:
        return False


def _wrap(argv: list[str], inbound: str) -> list[str]:
    return [sys.executable, "-I", str(_HELPER), "--fs", "0", "--net", "open",
            "--allow-socket", "/", "--inbound", inbound, "--", *argv]


#: Exits 0 only when the guard itself holds the probe: its supervisor sets
#: the marker in the one path where every bind is its own. A helper that
#: fell back to a lesser filter (no user notification on this kernel) runs
#: the probe too, but without the marker -- and that fallback refuses Unix
#: sockets outright, which an MPI run needs.
_PROBE = ("import os, sys; "
          f"sys.exit(0 if os.environ.get({_MARKER!r}) == 'held' else 3)")


def usable() -> bool:
    """Whether the guard holds a program on this host: probed by running a
    probe under it (about 50 ms). A probe that answered is kept for the
    process; one that could not run at all (an exception) is not, so a
    passing failure costs this one start its guard, not every later one."""
    global _usable
    if _usable is not None:
        return _usable
    if not (_HELPER.is_file() and sys.platform.startswith("linux")):
        _usable = False
        return False
    try:
        rc = subprocess.run(_wrap([sys.executable, "-I", "-c", _PROBE], "loopback"),
                            capture_output=True, timeout=30,
                            stdin=subprocess.DEVNULL).returncode
    except Exception:
        return False
    _usable = rc == 0
    return _usable


def guarded(argv, env=None) -> list[str]:
    """``argv`` under the inbound rule, or unchanged (see the module notes).

    ``env`` is the environment the program will get (default: this
    process's); the marker, the switch and the SLURM allocation are read
    from it."""
    argv = [str(a) for a in argv]
    env = os.environ if env is None else env
    if env.get(_MARKER, "").strip().lower() in ("held", "off", "0", "false", "no"):
        return argv
    if not usable():
        return argv
    return _wrap(argv, "any" if multi_node_allocation(env) else "loopback")
