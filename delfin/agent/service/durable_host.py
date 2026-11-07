"""Durable host: the agent supervisor as a SLURM job.

Package F1 (durable agent sessions). The supervisor/service normally run
on the login node. Targets 2-4 make it survive, and keep running, when the
host is a compute node reached only through the shared filesystem.

Nothing here opens a TCP port. All cross-node state travels through the
shared filesystem that already exists: the session store
(``session_store._SESSIONS_DIR``), pause flags (``session_pause``), approval
rooms (``file_confirm``), handoff/step journals. This module only:

* keeps a *supervisor heartbeat* so a second supervisor for the same user
  refuses to start while the first is alive (Target 2);
* renders a SLURM job script that starts the **same** supervisor with the
  **same** settings as the login path — never a wider permission and never
  as root (Target 4);
* documents the filesystem-only channel the login-node CLI uses to reach a
  session whose host is a compute node (Target 3).

Targets are red-tested per battery in ``tests/test_f1_adversarial_*.py``.
"""

from __future__ import annotations

import json
import re
import sys
import time
from pathlib import Path
from typing import Optional

from .. import session_store


# --------------------------------------------------------------------------
# Target 2 — supervisor heartbeat on the shared filesystem
# --------------------------------------------------------------------------

HEARTBEAT_FILE = "supervisor_heartbeat.json"
"""Basename of the supervisor heartbeat, stored under the sessions dir."""

FRESH_WINDOW_S = 90
"""How long a heartbeat is considered fresh (no wall-clock sleep in tests:
the fresh/stale boundary is decided only on the record's ``ts``)."""


def _heartbeat_path(root: Optional[str]) -> Path:
    """Resolve the heartbeat record path on the shared filesystem.

    The authoritative location is ``session_store._SESSIONS_DIR`` -- the
    same shared dir every supervisor and the login-node CLI already use,
    and the one the test suite redirects per run (``state_paths``). ``root``
    is accepted for call-site compatibility but is NOT trusted: a second
    supervisor on another node must find the file the first wrote, so the
    path is resolved from the session store every call, never from a local
    rundir.
    """
    return Path(session_store._SESSIONS_DIR) / HEARTBEAT_FILE


def read_supervisor_heartbeat(root: Optional[str] = None) -> Optional[dict]:
    """Return the heartbeat record (its clock ``ts``, pid, host, sessions)
    or None when absent/corrupt. Never raises."""
    try:
        record = json.loads(_heartbeat_path(root).read_text(encoding="utf-8"))
        float(record["ts"])
        return record
    except (OSError, ValueError, TypeError, KeyError):
        return None


def heartbeat_fresh(root: Optional[str] = None) -> bool:
    """True while a supervisor heartbeat is fresh (record ts within
    ``FRESH_WINDOW_S`` of now). Corrupt or absent => not fresh, so a second
    supervisor may start. The boundary moves only with the record's clock.
    """
    record = read_supervisor_heartbeat(root)
    if record is None:
        return False
    ts = float(record["ts"])
    return (time.time() - ts) <= FRESH_WINDOW_S


def write_supervisor_heartbeat(
    sessions: list,
    *,
    pid: int,
    host: Optional[str] = None,
    root: Optional[str] = None,
) -> Path:
    """Write a fresh heartbeat record for the supervised sessions. The real
    writer every durable-supervisor path (and every test) uses."""
    if host is None:
        try:
            import socket as _socket  # only to name the host, never to bind

            host = _socket.gethostname()
        except Exception:  # pragma: no cover - host naming failure is benign
            host = "unknown"
    record = {"pid": pid, "ts": time.time(), "host": host,
              "sessions": sorted(sessions)}
    p = _heartbeat_path(root)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(record), encoding="utf-8")
    return p


# --------------------------------------------------------------------------
# Target 4 — the rendered durable SLURM job (same supervisor, same sandbox)
# --------------------------------------------------------------------------

# Allow-lists mirror the discipline of slurm_tests.build_job (itself the
# model for a rendered, tested job script): refuse anything that could bleed
# into the rendered walltime or job name.
_TIME_RE = re.compile(r"^\d{2}:\d{2}:\d{2}$")
_VENV_RE = re.compile(r"^[A-Za-z0-9_/\-\. ]+$")
_NAME_RE = re.compile(r"^[A-Za-z0-9_\-]+$")


def _require(value: str, pattern: re.Pattern, what: str) -> str:
    if not value or not pattern.match(value):
        raise ValueError(f"refusing to render durable job: bad {what} {value!r}")
    return value


def render_durable_job(
    walltime: str = "12:00:00",
    *,
    partition: Optional[str] = None,
    qos: Optional[str] = None,
    signal_minutes_before: int = 5,
) -> str:
    """Render a SLURM job that runs the SAME supervisor on the node.

    Same supervisor command as ``launcher.start_service``
    (``python -m delfin.agent.service.supervisor``), the same venv
    interpreter that renders it, ``--export=ALL`` EXACTLY (no explicit var
    list), a hard ``set -euo pipefail`` guard, and a ``--signal`` well before
    the walltime end so the supervisor checkpoints and can requeue. No sudo,
    no root escalation, no privilege-related env is planted.
    """
    _require(walltime, _TIME_RE, "walltime")
    if signal_minutes_before < 0 or signal_minutes_before > 60:
        raise ValueError("signal_minutes_before must be 0..60")
    # The signal must fire before the walltime end; parse hh:mm:ss.
    hh, mm, _ss = (int(x) for x in walltime.split(":"))
    total_minutes = hh * 60 + mm
    if total_minutes - signal_minutes_before <= 0:
        raise ValueError("walltime too short for the pre-end signal")

    signal_before = f"{signal_minutes_before}:00" if signal_minutes_before else "0"

    lines = [
        "#!/usr/bin/env bash",
        "#SBATCH --export=ALL",
        f"#SBATCH --time={walltime}",
        f"#SBATCH --signal=INT@{signal_before}",
        "#SBATCH --job-name=delfin-durable",
        "#SBATCH --output=delfin-durable.%j.out",
        "#SBATCH --error=delfin-durable.%j.err",
        "set -euo pipefail",
    ]
    if partition is not None:
        _require(partition, _NAME_RE, "partition")
        lines.append(f"#SBATCH --partition={partition}")
    if qos is not None:
        _require(qos, _NAME_RE, "qos")
        lines.append(f"#SBATCH --qos={qos}")

    lines.append("")
    lines.append("# Keep the supervisor alive across requeues that follow a")
    lines.append("# pre-end --signal: the job only ends when the supervisor")
    lines.append("# has checkpointed and asked to stop.")
    lines.append("")
    # Same interpreter as the renderer (same venv, same settings).
    lines.append(f'exec "{sys.executable}" -m delfin.agent.service.supervisor')
    lines.append("")
    return "\n".join(lines)


# --------------------------------------------------------------------------
# Cross-node notes (Target 3) — maintained here as the shared surface
# --------------------------------------------------------------------------

__doc__ += """

Cross-node reachability (Target 3)
----------------------------------
A supervisor running on a compute node is NOT reachable through a unix
socket from the login node (a unix socket is per-host). DELFIN already
reaches a running session through the shared filesystem, and the channel
here is filesystem-only — it never offers the socket ``bind``, a
socketserver ``TCPServer`` listener, or an ``http`` server (no TCP port is
opened). The durable run uses the same channels in exactly the same way:

* **session state** — the step journal, resume note and checkpointed turn
  live under ``session_store._SESSIONS_DIR`` on the shared fs; ``delfin-agent
  service status`` reads them from the login node as from any node.
* **approvals** — ``FileConfirmBroker`` writes a request record that names
  the asking host (``file_confirm`` records ``host``; ``delfin-agent
  approvals ls/approve`` read that directory on the shared fs).
* **messages / handoff** — the handoff file and inbox flags are plain files
  under ``_SESSIONS_DIR``, written by one side and read by the other.
* **pause** — ``session_pause`` writes a flag file; ``delfin-agent pause``
  from the login node toggles the same file for the compute-node session.

The login-node CLI stays the same binary; only the host of the supervisor
moves. No port, no token, no wider permission.
"""
