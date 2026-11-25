"""The wave's tick: one idempotent pass that keeps every run moving and the
host safe.

Each tick, for every run under the waves root:

  * starts the tasks that became ready, within the limits (``wave.start``
    without waiting), so no long-lived starter can die and leave a queue;
  * advances the stage of a run started with stages (``stages.advance``);
  * records what changed (``wave.observe``): a member ended, went silent,
    left a traceback, became final;

and for the host:

  * leaves the envelope two ticks in a row (load above
    ``guard_max_load_per_core`` x cores, or less than
    ``guard_min_available_gb``): every member is stopped;
  * a process a member started listens beyond the loopback: every member is
    stopped. The socket guard binds every member's socket to the loopback,
    so this should never happen; it is the check that the guard held
    (2026-11-24: an agent's MPI listened on 0.0.0.0 for three minutes, and
    this check, run outside DELFIN then, ended it).

A stop writes ``<root>/EMERGENCY_STOPPED``; while it exists a tick only
records, and starts or resumes nothing, until a person clears it
(``clear_stop``). Alerts go to the attention inbox: a stop as a security
alert, a stalled member as a session state.

The tick runs in one of three ways, whichever the host has:

  * a systemd user service (``install_unit``): ``wave tick --every 60``
    from the harness, restarted when it fails, started again at boot;
  * the scheduler daemon (``install_entry``): an interval entry the daemon
    runs as a deterministic pass, never a model turn;
  * ``delfin-agent swarm wave tick --every 60`` under any supervisor
    (a cluster's tmux, a batch job), or ``wave tick`` once from a timer.
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

from . import stages, wave

PROMPT_MARKER = "[swarm-guard]"
_GUARD_DEFAULTS = {"guard_max_load_per_core": 1.5, "guard_min_available_gb": 4}


def _flag() -> Path:
    return wave.root() / "EMERGENCY_STOPPED"


def _state_file() -> Path:
    return wave.root() / ".guard.json"


def _load() -> dict:
    try:
        return json.loads(_state_file().read_text())
    except (OSError, ValueError):
        return {}


def _save(st: dict) -> None:
    p = _state_file()
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_name(f"{p.name}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(st))
    os.replace(tmp, p)


def _attention(kind: str, title: str, detail: str) -> None:
    try:
        from delfin.agent.attention import emit_attention
        emit_attention(kind, title=title, detail=detail[:400], workspace=str(wave.root()))
    except Exception:
        pass


def _limits() -> dict:
    try:
        from delfin.user_settings import load_settings
        own = ((load_settings() or {}).get("agent") or {}).get("swarm") or {}
    except Exception:
        own = {}
    return {k: own.get(k, v) for k, v in _GUARD_DEFAULTS.items()}


def outside_envelope() -> str:
    """Why the host is outside the envelope now, or ""."""
    lim = _limits()
    try:
        load1, cores = os.getloadavg()[0], os.cpu_count() or 1
        if load1 > cores * float(lim["guard_max_load_per_core"]):
            return f"load {load1:.0f} on {cores} cores"
    except OSError:
        pass
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemAvailable:"):
                gb = int(line.split()[1]) / 1048576
                if gb < float(lim["guard_min_available_gb"]):
                    return f"{gb:.1f} GB available"
    except OSError:
        pass
    return ""


def _children(pid: int) -> set:
    kids: dict = {}
    for name in os.listdir("/proc") if os.path.isdir("/proc") else []:
        if name.isdigit():
            try:
                ppid = int(Path(f"/proc/{name}/stat").read_text().rsplit(")", 1)[1].split()[1])
            except (OSError, IndexError, ValueError):
                continue
            kids.setdefault(ppid, []).append(int(name))
    out, todo = set(), [pid]
    while todo:
        p = todo.pop()
        if p not in out:
            out.add(p)
            todo += kids.get(p, [])
    return out


def listeners_beyond_loopback(pids: set) -> list[str]:
    """``ss`` lines of sockets the given processes listen on beyond the
    loopback ([] when ``ss`` is not there)."""
    out = []
    for flag in ("-tlnpH", "-ulnpH"):
        try:
            text = subprocess.run(["ss", flag], capture_output=True, text=True, timeout=20).stdout
        except (OSError, subprocess.SubprocessError):
            return []
        for line in text.splitlines():
            cols = line.split()
            if len(cols) < 4 or not {int(x) for x in re.findall(r"pid=(\d+)", line)} & pids:
                continue
            addr = cols[3]
            if not addr.startswith(("127.", "[::1]", "::1", "[::ffff:127.")):
                out.append(f"{flag[1]} {addr}")
    return out


def emergency_stop(why: str) -> dict:
    """Stop every member of every run, write the flag, alert."""
    stopped = {}
    for run in wave.runs():
        got = wave.stop(run)["stopped"]
        if got:
            stopped[run] = got
        wave.record(run, "emergency_stop", why=why)
    try:
        slice_name = wave.settings()["slice"]
        _systemctl("stop", slice_name)
    except (OSError, subprocess.SubprocessError):
        pass
    flag = _flag()
    flag.parent.mkdir(parents=True, exist_ok=True)
    flag.write_text(f"{time.strftime('%Y-%m-%d %H:%M:%S %Z')} {why}\n")
    _attention("security_alert", "Swarm stopped", why)
    return {"stopped": stopped, "why": why}


def clear_stop() -> bool:
    """Remove the flag a stop wrote (a person's decision). True when one was there."""
    try:
        _flag().unlink()
        return True
    except FileNotFoundError:
        return False


def tick(*, start_tasks: bool = True) -> dict:
    """One pass; returns what it did and saw. One pass at a time: while
    another runs (the service's, one run by hand), this one returns
    ``{"busy": True}`` at once and does nothing."""
    with wave._flock(wave.root() / ".tick.lock", wait=False) as mine:
        if not mine:
            return {"runs": {}, "stopped": None, "busy": True}
        return _tick(start_tasks)


def _tick(start_tasks: bool) -> dict:
    did: dict = {"runs": {}, "stopped": None}
    st = _load()
    flag = _flag()
    if flag.exists():
        for run in wave.runs():
            did["runs"][run] = {"events": wave.observe(run)}
        did["held"] = flag.read_text().strip()
        return did
    why = outside_envelope()
    st["outside"] = (st.get("outside", 0) + 1) if why else 0
    if st["outside"] >= 2:
        did["stopped"] = emergency_stop(f"the host left the envelope on two ticks: {why}")
        _save(st)
        return did
    members = wave.live_members()
    mine: set = set()
    for m in members:
        mine |= _children(int(m["pid"]))
    loud = listeners_beyond_loopback(mine) if mine else []
    if loud:
        did["stopped"] = emergency_stop("a member's process listens beyond the loopback: "
                                        + "; ".join(loud[:3]))
        _save(st)
        return did
    for run in wave.runs_oldest_first():
        rd: dict = {}
        state = wave.load_state(run)
        # Only a run started through this interface (wave_start sets it) is
        # driven: the root may hold runs another tool started, and their
        # queued tasks are not this tick's to launch.
        if start_tasks and state.get("driven") is True:
            try:
                tasks = json.loads((wave.run_dir(run) / "tasks.json").read_text())
                waiting = [t for t in tasks if t["id"] not in state["agents"]]
                if waiting:
                    rd["started"] = wave.start(run, tasks, state.get("repo") or ".", wait=False)["started"]
            except (OSError, ValueError, wave.WaveError, subprocess.CalledProcessError) as exc:
                rd["start_error"] = str(exc)[:300]
        if state.get("stages"):
            try:
                rd["stage"] = stages.advance(run)
            except wave.WaveError as exc:
                rd["stage_error"] = str(exc)[:300]
        ev = wave.observe(run)
        rd["events"] = ev
        for e in ev:
            if e["kind"] in ("silent", "traceback"):
                _attention("session_state", f"Swarm {run}/{e['id']}: {e['kind']}",
                           json.dumps({k: v for k, v in e.items() if k not in ("t", "run")}))
        did["runs"][run] = rd
    _save(st)
    return did


# -- the scheduler entry -------------------------------------------------------------

def install_entry(every_minutes: int = 1) -> dict:
    """Add (or replace) the guard's interval entry in the schedule the
    daemon runs. Idempotent; other entries stay as they are."""
    from delfin.agent.scheduler import Scheduler
    every = int(every_minutes) * 60
    if every < 60:
        raise ValueError("every_minutes must be >= 1")
    sched = Scheduler(path=Path.home() / ".delfin" / "cron.json")
    remove_entry(sched)
    root = wave.root()
    root.mkdir(parents=True, exist_ok=True)
    entry = sched.schedule_interval(every_seconds=every,
                                    prompt=f"{PROMPT_MARKER} tick of every wave, deterministic",
                                    reason="swarm guard, installed by the user",
                                    workspace=str(root))
    return {"id": entry.id, "every_seconds": every}


def remove_entry(sched=None) -> bool:
    from delfin.agent.scheduler import Scheduler
    sched = sched or Scheduler(path=Path.home() / ".delfin" / "cron.json")
    removed = False
    for ent in list(sched.list_entries()):
        if PROMPT_MARKER in ent.prompt:
            sched.delete(ent.id)
            removed = True
    return removed


# -- the tick as a loop and as a service ----------------------------------------------

UNIT = "delfin-wave-tick.service"


def loop(every_s: float = 60, *, ticks: int = 0, sleep=time.sleep, log=print) -> int:
    """Tick every ``every_s`` seconds, one line per tick; ``ticks`` > 0 ends
    after that many. A tick that raises is logged and the next one runs: the
    loop is the part that must not end. Returns the number of ticks run."""
    n = 0
    while True:
        try:
            did = tick()
            runs = did.get("runs") or {}
            moved = {r: sorted({e.get("kind", "?") for e in d.get("events") or []}
                               | ({"started"} if d.get("started") else set()))
                     for r, d in runs.items()}
            log(json.dumps({"t": time.strftime("%Y-%m-%d %H:%M:%S"), "runs": len(runs),
                            "moved": {r: k for r, k in moved.items() if k},
                            "held": bool(did.get("held")), "busy": bool(did.get("busy")),
                            "stopped": bool(did.get("stopped"))}))
        except Exception as exc:
            log(f"tick failed: {type(exc).__name__}: {str(exc)[:300]}")
        n += 1
        if ticks and n >= ticks:
            return n
        sleep(every_s)


def _unit_dir() -> Path:
    base = os.environ.get("XDG_CONFIG_HOME") or str(Path.home() / ".config")
    return Path(base) / "systemd" / "user"


def _q(text: str) -> str:
    """One systemd word: quoted, with % written as %%."""
    return '"' + str(text).replace("%", "%%").replace('"', '\\"') + '"'


def unit_text(every_s: int = 60, python: str = "") -> str:
    """The user service that runs ``wave tick --every`` from the harness (or
    the installed DELFIN when none is configured)."""
    cfg = wave.settings()
    harness = str(cfg.get("harness") or "")
    env = ["PYTHONUNBUFFERED=1"]
    if harness:
        env.append("PYTHONPATH=" + str(Path(harness).expanduser()))
    argv = [python or sys.executable, "-P", "-W", "ignore", "-m", "delfin.agent.cli",
            "swarm", "wave", "tick", "--every", str(int(every_s))]
    return "\n".join([
        "[Unit]",
        "Description=DELFIN waves: the tick (start, advance, record, guard)",
        "",
        "[Service]",
        "Type=simple",
        "WorkingDirectory=" + str(wave.root()).replace("%", "%%"),
        *("Environment=" + _q(e) for e in env),
        "ExecStart=" + " ".join(_q(a) for a in argv),
        "Restart=on-failure",
        "RestartSec=30",
        "Nice=10",
        "IOSchedulingClass=idle",
        "NoNewPrivileges=yes",
        "",
        "[Install]",
        "WantedBy=default.target",
        "",
    ])


def _systemctl(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["systemctl", "--user", *args], capture_output=True, text=True, timeout=60)


def install_unit(every_s: int = 60) -> dict:
    """Write the user service, enable and start it (idempotent: a running
    one is restarted on the new text). Raises WaveError where the host has
    no systemd user manager; ``wave tick --every`` runs anywhere."""
    if int(every_s) < 10:
        raise ValueError("every_s must be >= 10")
    try:
        probe = _systemctl("show-environment")
    except (OSError, subprocess.SubprocessError) as exc:
        probe = None
        why = f"{type(exc).__name__}"
    else:
        why = (probe.stderr or "").strip()[:200]
    if probe is None or probe.returncode != 0:
        raise wave.WaveError(f"no systemd user manager here ({why}); run "
                             "`delfin-agent swarm wave tick --every 60` under any supervisor, "
                             "or `wave guard install` for the scheduler daemon")
    path = _unit_dir() / UNIT
    path.parent.mkdir(parents=True, exist_ok=True)
    text = unit_text(every_s)
    changed = not path.exists() or path.read_text() != text
    path.write_text(text)
    wave.root().mkdir(parents=True, exist_ok=True)
    _systemctl("daemon-reload")
    on = _systemctl("enable", "--now", UNIT)
    if changed:
        _systemctl("restart", UNIT)
    active = _systemctl("is-active", UNIT).stdout.strip()
    return {"unit": str(path), "every_s": int(every_s), "active": active,
            "error": (on.stderr or "").strip()[:200] if on.returncode else ""}


def remove_unit() -> bool:
    """Stop, disable and delete the user service. True when one was there."""
    path = _unit_dir() / UNIT
    if not path.exists():
        return False
    try:
        _systemctl("disable", "--now", UNIT)
    except (OSError, subprocess.SubprocessError):
        pass
    path.unlink()
    try:
        _systemctl("daemon-reload")
    except (OSError, subprocess.SubprocessError):
        pass
    return True
