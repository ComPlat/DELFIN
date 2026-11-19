"""The messenger without root: ``delfin-agent messenger setup --user``.

For hosts where the person has no root (cluster login nodes). The same
three programs as scripts/pi/pi_messenger.sh -- the pinned simplex-chat
daemon, the bridge (python -m delfin.agent.service.bridge) and the gateway
(python -m delfin.messenger run) -- run as the user:

* as systemd --user units when ``systemctl --user`` answers, else as
  supervised background processes with a pid file each (``pid start-ticks``;
  ``stop --user`` ends a process only when both still match);
* at low priority (nice 10, the idle I/O class), without whisper, and with
  an idle exit (``--idle-exit HOURS``, default 12; 0: never);
* state under ``$XDG_DATA_HOME/delfin-messenger`` (0700), units under
  ``$XDG_CONFIG_HOME/systemd/user``.

What changes against the root install, and is said before anything starts:

* the daemon's API listens on 127.0.0.1 at a random port. The root install
  lets only the gateway's user reach it (an nftables rule); here no rule
  can, so ANY local user of the host can connect to it and read the bot's
  chats and write as the bot. Setup runs only with
  ``--accept-loopback-risk``. (Approvals are not reachable that way: they
  go through the bridge's Unix socket, below.)
* the gateway, the bridge and the agent service are one user. The bridge
  admits its own uid only for the gateway's recorded process (pid and
  start time from the gateway's pid file, its environment readable and
  without the agent marker): bridge.rootless_peer_ok. The boundary against
  the agent's own commands is the agent's command cage (bubblewrap, the
  socket guard), as in people mode; without a cage, a command you approve
  runs as you and can reach what you can.

Rejected: an unpinned simplex-chat build (refused, the digests are the
installer's, pinned by a test); a private network namespace for the daemon
(it needs the network, and slirp4netns or pasta is not on every host).
"""
from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import signal
import socket
import subprocess
import sys
from pathlib import Path

#: The pinned simplex-chat release and its digests (scripts/pi/pi_messenger.sh;
#: test_the_messenger_runs_without_root keeps the two copies equal).
SIMPLEX_VERSION = "v7.0.3"
PINNED = frozenset({
    "895fb14cfaa662d1c0947f2f871141fc89680fe3e70bff64c366cbe9e59aa4f0",
    "2d2e62351f11bc51ae659618584722b38ea5a6796b806c9a388ee6e3124c90ac",
    "5afb1d25efe5ccf564a1ab124bc7f410a7a73171c974a4d8b7f4d8f2e3d62e77",
    "91f332edc813cab78a848da7ba75ed91bcb8d09abf1fc92eca74330dc29ec970",
})
UNITS = ("delfin-messenger-chat.service", "delfin-messenger-bridge.service",
         "delfin-messenger-gateway.service")
_ORDER = (("daemon", UNITS[0]), ("bridge", UNITS[1]), ("gateway", UNITS[2]))

RISK = (
    "Rootless mode: the SimpleX daemon's API listens on 127.0.0.1 without a firewall rule, "
    "so other users of this host can connect to it, read the bot's chats and write as the bot. "
    "Approvals stay behind the bridge's Unix socket (0600, your gateway process only). "
    "Run it only where you accept that; then add --accept-loopback-risk.")
LOGIN_NODE = (
    "On a cluster login node: some sites forbid long-running background programs there. "
    "Check your site's rules first. The messenger runs at low priority (nice 10, idle I/O), "
    "polls every 5 s, and exits after --idle-exit hours without a message.")


def state_dir() -> Path:
    base = os.environ.get("XDG_DATA_HOME", "").strip() or str(Path.home() / ".local" / "share")
    return Path(base) / "delfin-messenger"


def unit_dir() -> Path:
    base = os.environ.get("XDG_CONFIG_HOME", "").strip() or str(Path.home() / ".config")
    return Path(base) / "systemd" / "user"


def free_port() -> int:
    """A free TCP port on 127.0.0.1, chosen by the kernel."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("127.0.0.1", 0))
        return int(s.getsockname()[1])


def start_ticks(pid) -> str:
    """A process's start time (clock ticks since boot), '' where unknown."""
    from .pairing import _start_time
    return _start_time(pid) or ""


def _environment(port, bot_name):
    sd = state_dir()
    return {
        "DELFIN_MESSENGER_PORT": str(port),
        "DELFIN_MESSENGER_ROOTLESS": "1",
        "DELFIN_MESSENGER_DIR": str(sd / "gateway"),
        "DELFIN_MESSENGER_SOCKET": str(sd / "run" / "bridge.sock"),
        "DELFIN_MESSENGER_GATEWAY_PIDFILE": str(sd / "run" / "gateway.pid"),
        "DELFIN_MESSENGER_LABEL": str(bot_name).removeprefix("DELFIN-") or "rootless",
    }


def plan(*, simplex_bin, port, bot_name, python=None, idle_exit_s=0) -> dict:
    """What runs: ``env`` (the messenger's settings), ``processes`` (argv
    of daemon, bridge, gateway) and ``units`` (the systemd --user files)."""
    python = python or sys.executable
    sd = state_dir()
    env = _environment(port, bot_name)
    processes = {
        "daemon": [str(simplex_bin), "-d", str(sd / "chat" / "simplex.db"), "--chat-server-port", str(port),
                   "--files-folder", str(sd / "files"), "--auto-accept-files", "0",
                   "--create-bot-display-name", str(bot_name), "--create-bot-allow-files"],
        "bridge": [python, "-m", "delfin.agent.service.bridge"],
        "gateway": [python, "-m", "delfin.messenger", "run"]
                   + (["--idle-exit", str(int(idle_exit_s))] if idle_exit_s else []),
    }
    units = {}
    after = {"daemon": "network-online.target", "bridge": UNITS[0], "gateway": f"{UNITS[0]} {UNITS[1]}"}
    for name, unit in _ORDER:
        lines = ["[Unit]", f"Description=DELFIN messenger ({name}, rootless)", f"After={after[name]}",
                 "", "[Service]"]
        lines += [f"Environment={k}={v}" for k, v in sorted(env.items())]
        lines += ["ExecStart=" + " ".join(_quote(a) for a in processes[name]),
                  # An idle exit is not a failure: the gateway stays down until started.
                  "Restart=on-failure", "RestartSec=10", "Nice=10", "IOSchedulingClass=idle",
                  "UMask=0077", "", "[Install]", "WantedBy=default.target", ""]
        units[unit] = "\n".join(lines)
    return {"env": env, "processes": processes, "units": units}


def _quote(arg: str) -> str:
    return arg if arg and all(c.isalnum() or c in "/._-=:" for c in arg) else '"' + arg.replace('"', '\\"') + '"'


def systemd_user(run=subprocess.run) -> bool:
    """Whether ``systemctl --user`` answers (a user manager is running)."""
    try:
        done = run(["systemctl", "--user", "is-system-running"], capture_output=True, text=True,
                   timeout=10, check=False)
    except (OSError, subprocess.SubprocessError):
        return False
    return done.returncode == 0 or (done.stdout or "").strip() in {"running", "degraded"}


def _spawn(name, argv, env, log):
    """One supervised background process (its own session, output to
    ``log``); the pid."""
    with open(log, "ab") as out:
        proc = subprocess.Popen(argv, env=env, stdin=subprocess.DEVNULL, stdout=out,
                                stderr=subprocess.STDOUT, start_new_session=True)
    return proc.pid


def _low(argv):
    """``argv`` at low priority: nice 10 and, where ionice exists, idle I/O."""
    prefix = ["nice", "-n", "10"]
    if shutil.which("ionice"):
        prefix += ["ionice", "-c", "3"]
    return prefix + list(argv)


def _digest(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def setup(argv, *, out=print, run=subprocess.run, spawn=_spawn, ask=input) -> int:
    """``delfin-agent messenger setup --user``: check, install, start, pair."""
    parser = argparse.ArgumentParser(prog="delfin-agent messenger setup --user")
    parser.add_argument("--user", action="store_true")
    parser.add_argument("--simplex-chat", default="", help="the pinned simplex-chat binary")
    parser.add_argument("--port", type=int, default=0)
    parser.add_argument("--name", default="")
    parser.add_argument("--idle-exit", type=float, default=12.0, help="hours without a message; 0: never")
    parser.add_argument("--accept-loopback-risk", action="store_true")
    parser.add_argument("--no-pair", action="store_true")
    args = parser.parse_args([a for a in argv if a != "--user"])
    from .pairing import _from_an_agent
    if _from_an_agent():
        out("Refused: an agent session cannot set up the messenger; run it in your own terminal.")
        return 2
    out(LOGIN_NODE)
    if not args.accept_loopback_risk:
        out(RISK)
        return 2
    binary = args.simplex_chat or shutil.which("simplex-chat") or ""
    if not binary or not Path(binary).is_file():
        out(f"simplex-chat {SIMPLEX_VERSION} not found. Download the build for this machine from "
            f"github.com/simplex-chat/simplex-chat/releases/tag/{SIMPLEX_VERSION} and pass "
            "--simplex-chat PATH; its SHA-256 must be one of the pinned digests.")
        return 2
    digest = _digest(binary)
    if digest not in PINNED:
        out(f"Refused: {binary} has SHA-256 {digest}, not one of the pinned {SIMPLEX_VERSION} digests.")
        return 2
    port = args.port or free_port()
    name = args.name or "DELFIN-" + (socket.gethostname().split(".")[0] or "rootless")
    p = plan(simplex_bin=binary, port=port, bot_name=name, idle_exit_s=int(args.idle_exit * 3600))
    sd = state_dir()
    for sub in ("", "chat", "files", "gateway", "run", "log"):
        (sd / sub).mkdir(parents=True, exist_ok=True)
        os.chmod(sd / sub, 0o700)
    if systemd_user(run):
        ud = unit_dir()
        ud.mkdir(parents=True, exist_ok=True)
        for unit, text in p["units"].items():
            (ud / unit).write_text(text)
        run(["systemctl", "--user", "daemon-reload"], check=False, timeout=30)
        run(["systemctl", "--user", "enable", "--now", *UNITS], check=False, timeout=60)
        out(f"Started as systemd --user units ({', '.join(UNITS)}); the API on 127.0.0.1:{port}.")
    else:
        env = {**os.environ, **p["env"]}
        for name_, _unit in _ORDER:
            pid = spawn(name_, _low(p["processes"][name_]), env, sd / "log" / f"{name_}.log")
            if name_ != "gateway":
                # The gateway writes its own (delfin/messenger/__main__.py).
                (sd / "run" / f"{name_}.pid").write_text(f"{pid} {start_ticks(pid)}\n")
        out(f"Started in the background (pid files in {sd / 'run'}); the API on 127.0.0.1:{port}. "
            "Stop: delfin-agent messenger stop --user.")
    if not args.no_pair:
        from .pairing import cli_main
        return cli_main(["pair"], ask=ask, out=out)
    return 0


def stop(*, out=print, kill=os.kill, run=subprocess.run) -> int:
    """``delfin-agent messenger stop --user``: the units, or each recorded
    process whose pid AND start time still match (a reused pid is left)."""
    if systemd_user(run) and any((unit_dir() / u).exists() for u in UNITS):
        run(["systemctl", "--user", "stop", *UNITS], check=False, timeout=60)
        out("Stopped the messenger's user units.")
        return 0
    for name, _unit in reversed(_ORDER):
        path = state_dir() / "run" / f"{name}.pid"
        try:
            pid_text, ticks = (path.read_text().split() + ["", ""])[:2]
            pid = int(pid_text)
        except (OSError, ValueError):
            continue
        if ticks and start_ticks(pid) == ticks:
            try:
                kill(pid, signal.SIGTERM)
                out(f"Stopped {name} ({pid}).")
            except OSError:
                pass
        path.unlink(missing_ok=True)
    return 0
