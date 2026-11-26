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
    ``min_available_gb()`` -- the larger of a floor and a share of the
    host's total): every member is stopped;
  * a process a member started listens beyond the loopback: every member is
    stopped. The socket guard binds every member's socket to the loopback,
    so this should never happen; it is the check that the guard held
    (2026-11-24: an agent's MPI listened on 0.0.0.0 for three minutes, and
    this check, run outside DELFIN then, ended it).
  * a member sits stalled -- ended 15+ minutes ago with nothing done by its
    stage, alive with a log that has not grown for 60 minutes, or a real run
    past its ceiling -- is raised once through the attention path
    (``find_stalls``);
  * what the members ran into is counted as it appears (``friction.sweep``,
    each log read from where the last tick stopped) into
    ``<root>/.friction.jsonl``, and ``<root>/friction.json`` is rewritten
    at most once an hour.

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

from . import friction, stages, wave

PROMPT_MARKER = "[swarm-guard]"

#: The envelope the tick holds the host inside, and why these numbers.
#:
#: The operator script the tick replaced (a cron line every five minutes,
#: until 2026-11-24) stopped everything at ``load1 > 340`` or
#: ``MemAvailable < 150 GB``, read twice in a row. The first defaults below
#: were looser than that on every machine: 1.5 x cores is 576 on a 384-core
#: host (1.7x the script's 340) and 4 GB is 2% of its 150 GB -- a shared
#: host was left unguarded, which is the one job this guard has.
#:
#: 0.885 load per core and ``max(4 GB, 10 % of total)`` reproduce the old
#: protection on that host, measured on it (384 cores, MemTotal 1511.3 GiB,
#: ``os.cpu_count`` 384): 0.885 x 384 = 339.8 against the script's 340, and
#: 10 % of 1511.3 = 151.1 GiB against its 150. It scales with the machine so
#: a small one is not stopped at idle: on 8 cores the load stop is 7.1, and
#: the 4 GB floor -- the number the tick shipped with -- decides below 40 GB
#: of total, so a 16 GB laptop keeps the stop it already had and only hosts
#: above 40 GB get a stricter one. An 8 GB floor would have been 50 % of a
#: 16 GB laptop, far heavier in proportion than 150 GB of 1.5 TB, which is
#: the stopped-at-idle failure this is meant to avoid. The script's numbers
#: were absolute and belonged to one machine; per core and per total, one set
#: of defaults protects a laptop and the cluster host alike.
#:
#: Both stay a setting (``agent.swarm.guard_max_load_per_core`` /
#: ``guard_min_memory_frac`` / ``guard_min_memory_floor_gb``), so a host
#: that genuinely wants the looser 1.5 keeps it by name.
#:
#: The load stop is the envelope's only LOAD number, and ``wave._room``
#: holds its own ``max_load_per_core`` for a different question: how many
#: members may *start*. The two are deliberately equal (both 0.885), so a
#: host the envelope is about to stop can never be one the tick keeps
#: launching onto -- there is no band between them. They are duplicated
#: rather than imported because the guard imports this module and the
#: dependency cannot run the other way; they are pinned together by
#: ``test_the_start_gate_stops_before_the_envelope_does``. A start gate
#: left above the stop line would mean launching work onto a host the
#: envelope condemns, which is what that test exists to prevent.
_GUARD_DEFAULTS = {"guard_max_load_per_core": 0.885,
                   "guard_min_memory_frac": 0.10,
                   "guard_min_memory_floor_gb": 4.0,
                   # A member that ended 15+ minutes ago whose stage did
                   # nothing with it: the stage resumes on endpoint errors
                   # only after a 900 s pause, so the floor sits just past
                   # that pause -- below it every resumed member would
                   # raise an alarm on its way to being resumed.
                   "stall_ended_min": 15.0,
                   # A live member whose log has not grown in an hour.
                   # ``wave.observe`` records a "silent" event from 60
                   # minutes; a member that is merely thinking for a few
                   # minutes is normal, so this is an attention item, not a
                   # stop.
                   "stall_silent_min": 60.0,
                   # A real run (base + branch, each probe) that is still
                   # going after 2 h 30: each side carries a 1500 s ceiling
                   # of its own, so this is several probes deep, or the
                   # thread that started it died with the marker left
                   # behind (2026-11-24: an integrator's real run held a
                   # wave for hours and nothing said so).
                   "stall_realrun_min": 150.0,
                   # The friction ledger appends every tick (it reads each
                   # log from the byte it stopped at); the grouped report is
                   # rewritten at most this often.
                   "friction_report_min": 60.0}


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


def _now() -> float:
    """This pass's idea of now. One seam so a test can move the clock the
    stall and report-age checks read without touching the members' timing."""
    return time.time()


def _attention(kind: str, title: str, detail: str) -> None:
    try:
        from delfin.agent.attention import emit_attention
        emit_attention(kind, title=title, detail=detail[:400], workspace=str(wave.root()))
    except Exception:
        pass


def _limits() -> dict:
    """``agent.swarm`` over the defaults, with one translation: a host that
    already set the old absolute ``guard_min_available_gb`` keeps it, as the
    floor -- an override a person wrote must not stop working because the
    ceiling beside it became a share."""
    try:
        from delfin.user_settings import load_settings
        own = ((load_settings() or {}).get("agent") or {}).get("swarm") or {}
    except Exception:
        own = {}
    lim = {k: own.get(k, v) for k, v in _GUARD_DEFAULTS.items()}
    if "guard_min_available_gb" in own and "guard_min_memory_floor_gb" not in own:
        lim["guard_min_memory_floor_gb"] = own["guard_min_available_gb"]
    return lim


def _meminfo_kb() -> dict:
    """``MemTotal`` and ``MemAvailable`` of this host in kB ({} when /proc
    cannot say)."""
    try:
        out = {}
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith(("MemTotal:", "MemAvailable:")):
                parts = line.split()
                out[parts[0][:-1]] = int(parts[1])
        return out
    except (OSError, IndexError, ValueError):
        return {}


def host_readings() -> tuple:
    """``(load1, cores, available_gb, total_gb)`` as this host reports them.

    One seam on purpose: what decides a stop is the arithmetic in
    :func:`outside_envelope`, and a test that fakes this can pin that
    arithmetic on a 384-core host and on a laptop without either reaching
    /proc."""
    try:
        load1 = os.getloadavg()[0]
    except OSError:
        load1 = 0.0
    info = {k: v / 1048576 for k, v in _meminfo_kb().items()}
    # None means "the host did not say", which is not the same as 0.0: an
    # available memory of zero is total exhaustion and must stop, so a test
    # of truthiness here would fail open on exactly the host in trouble.
    return (load1, os.cpu_count() or 1, info.get("MemAvailable"), info.get("MemTotal"))


def min_available_gb(lim: dict | None = None, total_gb: float | None = None) -> float:
    """The memory floor of the envelope: the larger of a floor and a share
    of this host's total, so the same default stops a 1.5 TB host at 151 GB
    and a 16 GB one at the 4 GB floor. With no MemTotal to read, the floor
    alone."""
    lim = lim or _limits()
    if total_gb is None:
        total_gb = host_readings()[3]
    if total_gb is None:
        return float(lim["guard_min_memory_floor_gb"])
    return max(float(lim["guard_min_memory_floor_gb"]),
               float(lim["guard_min_memory_frac"]) * float(total_gb))


def outside_envelope() -> str:
    """Why the host is outside the envelope now, or ""."""
    lim = _limits()
    load1, cores, avail_gb, total_gb = host_readings()
    if load1 > cores * float(lim["guard_max_load_per_core"]):
        return f"load {load1:.0f} on {cores} cores"
    if avail_gb is not None and avail_gb < min_available_gb(lim, total_gb):
        return f"{avail_gb:.1f} GB available"
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
    lim = _limits()
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
            # "silent" is NOT alerted here any more: find_stalls below owns
            # every stall and deduplicates it, and alerting the same quiet
            # member from both paths sent it twice per tick. The event is
            # still recorded in the run's log by observe.
            if e["kind"] == "traceback":
                _attention("session_state", f"Swarm {run}/{e['id']}: {e['kind']}",
                           json.dumps({k: v for k, v in e.items() if k not in ("t", "run")}))
        did["runs"][run] = rd
    found = find_stalls(lim)
    # Called every tick, not only when something is stalled: the call is what
    # forgets the stalls that have cleared, and the tick where the last one
    # clears is the only tick that can forget it. Under `if found:` that tick
    # was skipped, so a member that stalled, recovered and stalled again was
    # never reported the second time -- its key stayed in the state and
    # deduplicated it forever.
    did["stalls"] = raise_stalls(found, st.setdefault("stalls", {}))
    did["friction"] = _ledger(lim)
    _save(st)
    return did


def _ledger(lim: dict) -> dict:
    """Count what the members ran into, once per tick, cheaply.

    ``friction.sweep`` reads each member log from the byte the last sweep
    stopped at, so a tick over a full waves root costs a stat per log plus
    the bytes written since the last one; the grouped report is rewritten
    at most every ``friction_report_min`` minutes (it re-reads the whole
    ledger, and per-model rates do not move every minute). The age that
    decides it is the report file's own mtime, not a counter in the tick's
    state: a report that stopped being rewritten because the root went
    quiet still has to catch up, and a restarted tick inherits the truth
    from disk. A ledger that cannot be written (a full root, a locked file)
    is said in the tick's line and tried again next tick -- the count is the
    point, not the alert."""
    try:
        out = friction.sweep()
    except Exception as exc:
        return {"error": f"{type(exc).__name__}: {str(exc)[:160]}"}
    every = float(lim["friction_report_min"]) * 60
    report = wave.root() / "friction.json"
    try:
        age = _now() - report.stat().st_mtime
    except OSError:
        age = every + 1.0                             # never written
    if age >= every:
        try:
            out["report"] = friction.write_report()["records"]
        except Exception as exc:
            out["report_error"] = str(exc)[:160]
    return out


# -- stalls ---------------------------------------------------------------------------

#: A member can sit in three ways, and until the tick took over none of them
#: was said aloud: a member that ended and whose stage did nothing with it
#: (the stage resumes on endpoint errors and asks for real runs; a member
#: that falls through every branch simply stops being anyone's business), a
#: member alive but writing nothing for an hour, and a real run whose thread
#: died and left its ``<tid>.realrun_running`` marker behind -- the stage
#: then refuses the member forever, because the marker says a run is in
#: flight. ``observe`` records the first two as events in the run's own log,
#: which nobody read, so a wave could hold for hours with the tick passing
#: over it happily. Each one goes to the attention inbox ONCE per member and
#: kind (``stalls`` in the tick's state), not once per tick; a member that
#: stops being stalled is forgotten, so it can be reported again later.
def find_stalls(lim: dict | None = None, now: float | None = None) -> list[dict]:
    """Every member sitting stalled, each with ``run``, ``id``, ``kind``,
    ``min`` and ``why``. Reads only; raises nothing."""
    lim = lim or _limits()
    now = _now() if now is None else now
    out = []
    for run in wave.runs():
        try:
            state = wave.load_state(run)
        except wave.WaveError:
            continue
        d = wave.run_dir(run)
        try:
            task_ids = [t["id"] for t in json.loads((d / "tasks.json").read_text())]
        except (OSError, ValueError, KeyError, TypeError):
            task_ids = []
        agents = state["agents"]
        for tid in (task_ids or list(agents)):
            a = agents.get(tid)
            if a is None or (d / f"{tid}.final").exists():
                continue                              # never started, or done
            at = {"run": run, "id": tid}
            # The real run first, and before the member's own state: the
            # stage starts it for a member that has ENDED, so "ended" must
            # not swallow it -- and while its marker stands, the stage
            # leaves the member alone, which is exactly when a stale marker
            # (a thread that died with it) holds the wave forever.
            marker = d / f"{tid}.realrun_running"
            if marker.exists():
                try:
                    who = int(marker.read_text().strip())
                    age = (now - marker.stat().st_mtime) / 60
                except (OSError, ValueError):
                    out.append({**at, "kind": "realrun", "min": 0,
                                "why": "its real run's marker cannot be read"})
                    continue
                if not wave.alive(who):
                    out.append({**at, "kind": "realrun", "min": round(age),
                                "why": "its real run's marker is held by a process "
                                       "that is gone"})
                elif age >= float(lim["stall_realrun_min"]):
                    out.append({**at, "kind": "realrun", "min": round(age),
                                "why": "a real run older than its ceiling"})
                continue
            if not wave.alive(a["pid"]):
                # "Ended with nothing done by its stage" presupposes that a
                # stage of OURS owns the member. Every other action this tick
                # takes is guarded the same way (see the start block below):
                # the root may hold runs another tool started, and their
                # members are not this tick's to judge. Measured on arm a2
                # without this gate: a foreign run and a stages-off run both
                # put an `ended` item in the inbox, and so did a member whose
                # stage had already acted on it (`<tid>.asked`).
                if state.get("driven") is not True:
                    continue
                marks = [m for m in ("asked", "realrun", "final")
                         if (d / f"{tid}.{m}").exists()]
                if marks:
                    continue                          # the stage did act on it
                quiet = (now - wave._ended_at(f"{run}/{tid}")) / 60
                if quiet >= float(lim["stall_ended_min"]):
                    out.append({**at, "kind": "ended", "min": round(quiet),
                                "why": "ended and the stage did nothing with it"})
                continue
            try:
                last = Path(a["log"]).stat().st_mtime
            except OSError:
                last = 0
            # A resume keeps the member's log (``wave._launch`` opens it "ab",
            # which does not touch its mtime) and gives it a fresh pid and a
            # fresh ``started``. Judged on the log alone, a member resumed a
            # second ago after two hours of work reads as "silent for 120
            # minutes" -- measured: a member with ``started = now`` and a
            # 90-minute-old log was reported silent 90. The session is the
            # newer of the two, so the newer of the two is what "has not grown"
            # is measured from.
            started = a.get("started")
            if isinstance(started, (int, float)) and started > last:
                last = float(started)
            idle = (now - last) / 60
            if idle >= float(lim["stall_silent_min"]):
                out.append({**at, "kind": "silent", "min": round(idle),
                            "why": "alive, nothing written to its log"})
    return out


def raise_stalls(found: list[dict], seen: dict) -> int:
    """Send each newly found stall through the attention path once, and
    forget the ones no longer seen so they can be raised again. Returns how
    many were sent this pass."""
    raised = 0
    for s in found:
        key = f"{s['run']}/{s['id']}:{s['kind']}"
        if key in seen:
            continue
        seen[key] = s["min"]
        raised += 1
        _attention("session_state", f"Swarm {s['run']}/{s['id']} stalled ({s['kind']})",
                   json.dumps({k: v for k, v in s.items() if k != "run"}))
    live = {f"{s['run']}/{s['id']}:{s['kind']}" for s in found}
    for key in [k for k in seen if k not in live]:
        seen.pop(key, None)
    return raised


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
                            "stopped": bool(did.get("stopped")),
                            "stalls": did.get("stalls") or 0,
                            "friction": (did.get("friction") or {}).get("appended", 0)}))
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
