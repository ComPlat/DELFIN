"""Command isolation on macOS: Seatbelt (``sandbox-exec``).

macOS has neither bubblewrap nor Landlock nor seccomp. It has the Seatbelt
sandbox, which the system and coding-agent sandboxes use: a profile in the
Sandbox Profile Language, applied by ``sandbox-exec`` to the command and
everything it starts. The same promise as on Linux:

* writes only beneath the workspace roots and a private temp directory;
* the credential locations unreadable;
* Unix sockets only beneath an allow-list (``deny``-listed session doors in
  the attended profile);
* with a restricted network, TCP only to the egress proxy's loopback port.

Seatbelt cannot filter a remote IP address other than ``localhost``, so the
cloud metadata address is not refused by address here; in the restricted
modes it is unreachable like every other destination.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from typing import Iterable, Optional

_SANDBOX_EXEC = "/usr/bin/sandbox-exec"
_FUNCTIONAL: Optional[bool] = None


def available() -> bool:
    """Whether ``sandbox-exec`` runs a profile here. Probed once."""
    global _FUNCTIONAL
    if _FUNCTIONAL is not None:
        return _FUNCTIONAL
    ok = False
    if sys.platform == "darwin" and (os.path.exists(_SANDBOX_EXEC) or shutil.which("sandbox-exec")):
        try:
            r = subprocess.run([_SANDBOX_EXEC, "-p", "(version 1)(allow default)", "/usr/bin/true"],
                               capture_output=True, timeout=10)
            ok = r.returncode == 0
        except (OSError, subprocess.TimeoutExpired):
            ok = False
    _FUNCTIONAL = ok
    return ok


def _q(path: str) -> str:
    return '"' + str(path).replace("\\", "\\\\").replace('"', '\\"') + '"'


def _real(paths: Iterable[str]) -> list[str]:
    out = []
    for p in paths:
        if not p:
            continue
        real = os.path.realpath(p)
        if real not in out:
            out.append(real)
    return out


def profile(*, write_roots: Iterable[str] = (), hide: Iterable[str] = (),
            allow_sockets: Iterable[str] = (), deny_sockets: Iterable[str] = (),
            restrict_files: bool = True, net_mode: str = "open",
            proxy_port: int = 0, write_deny: Iterable[str] = ()) -> str:
    """The Seatbelt profile. ``restrict_files`` False is the attended
    profile: the filesystem as it is, the session doors still shut."""
    lines = ["(version 1)", "(allow default)"]
    if restrict_files:
        roots = _real(write_roots)
        lines.append("(deny file-write*)")
        allowed = " ".join(f"(subpath {_q(r)})" for r in roots)
        lines.append("(allow file-write* " + allowed +
                     ' (literal "/dev/null") (literal "/dev/zero")'
                     ' (literal "/dev/dtracehelper") (regex #"^/dev/tty")'
                     ' (regex #"^/dev/fd/"))')
        denied = _real(write_deny)
        if denied:
            # After the allow, so it wins: git hooks stay readable, a new
            # one cannot be written.
            lines.append("(deny file-write* " +
                         " ".join(f"(subpath {_q(d)})" for d in denied) + ")")
    hidden = _real(hide)
    if hidden:
        lines.append("(deny file-read* file-write* " +
                     " ".join(f"(subpath {_q(h)})" for h in hidden) + ")")
    sockets = _real(allow_sockets)
    if sockets:
        lines.append("(deny network-outbound (remote unix-socket))")
        lines.append("(allow network-outbound " +
                     " ".join(f"(remote unix-socket (subpath {_q(s)}))" for s in sockets) +
                     ' (remote unix-socket (path-literal "/private/var/run/mDNSResponder")))')
    patterns = [d[len("regex:"):] for d in deny_sockets if str(d).startswith("regex:")]
    doors = _real(d for d in deny_sockets if not str(d).startswith("regex:"))
    if doors or patterns:
        lines.append("(deny network-outbound " +
                     " ".join(f"(remote unix-socket (subpath {_q(d)}))" for d in doors) +
                     "".join(f' (remote unix-socket (regex #"{pat}"))' for pat in patterns) + ")")
    if net_mode in ("proxy", "none"):
        lines.append("(deny network-outbound (remote ip))")
        if net_mode == "proxy" and proxy_port:
            lines.append(f'(allow network-outbound (remote ip "localhost:{int(proxy_port)}"))')
    # Another process's arguments and environment are not the command's to read.
    lines.append("(deny process-info* (target others))")
    return "\n".join(lines) + "\n"


def argv(plain: list[str], prof: str) -> list[str]:
    return [_SANDBOX_EXEC, "-p", prof, *plain]


def session_doors() -> list[str]:
    """macOS sockets through which a command reaches the user's sessions."""
    uid = os.getuid()
    home = os.path.expanduser("~")
    doors = [f"/private/tmp/tmux-{uid}", os.path.join(home, ".ssh"),
             r"regex:^/private/tmp/com\.apple\.launchd\.", "/var/run/docker.sock",
             os.path.join(home, ".docker", "run"), "/private/tmp/.X11-unix"]
    agent = os.environ.get("SSH_AUTH_SOCK", "")
    if agent:
        doors.append(os.path.dirname(agent))
    return doors
