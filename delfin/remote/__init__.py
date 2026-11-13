"""DELFIN on HPC clusters, driven from the user's PC through cluster-tunnel (ctun).

ctun keeps one logged-in ssh master per cluster: a person logs in once with
password and OTP, and every later command reuses that connection. The
master is the user's whole identity on the cluster, so only DELFIN's own
process drives it (``connection``), for the clusters the user named in
``~/.delfin_settings.json`` and only once ctun has their host key pinned
(``clusters``). A command an agent runs reaches none of it: ctun's
directories are masked in every cage and denied to every file tool
(``state_dirs``), and ``connection.connect`` refuses inside an agent session.
"""

from __future__ import annotations

import os
from pathlib import Path

__all__ = ["state_dirs"]


def state_dirs() -> list[str]:
    """ctun's directories, created now and resolved: its config directory
    (config.yaml, the budget scripts, the pinned host keys), its cache
    directory (control sockets, sessions with their budget limit, the
    command log) and a socket directory configured elsewhere.

    Created so that a cage masks them from the start: an agent that could
    plant ~/.config/cluster-tunnel/config.yaml before the user writes one
    would choose the ssh options DELFIN's own process runs with. Never
    raises; an unusable location is left out.
    """
    out: list[str] = []
    try:
        from cluster_tunnel import paths as ctun_paths
    except Exception:
        return out
    dirs: list[Path] = []
    try:
        dirs += [ctun_paths.config_dir(), ctun_paths.cache_dir()]
    except Exception:
        return out
    try:
        from cluster_tunnel import api as ctun_api
        socket_dir = ctun_api.load_config().defaults.socket_dir
        if socket_dir:
            dirs.append(Path(socket_dir).expanduser())
    except Exception:
        pass  # no config yet, or a broken one: the default locations stand
    home = Path.home()
    # Masking one of these would hide far more than ctun's state.
    never = {Path("/"), Path("/tmp"), Path("/run"), home}
    try:
        never.add(home.resolve())
    except OSError:
        pass
    for d in dirs:
        try:
            d.mkdir(parents=True, exist_ok=True, mode=0o700)
            real = d.resolve()
        except OSError:
            continue
        if real in never or not real.is_dir():
            continue
        if any(real == Path(o) or Path(o) in real.parents for o in out):
            continue  # inside a directory already listed
        out.append(str(real))
    return out


def _inside_agent_session() -> bool:
    try:
        from delfin.agent.approval_guard import inside_agent_session
        return inside_agent_session()
    except Exception:
        return bool(os.environ.get("DELFIN_AGENT_SESSION"))
