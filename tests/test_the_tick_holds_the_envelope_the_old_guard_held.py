"""The tick holds the server envelope the old guard held, and keeps the ledger.

The tick (``delfin/agent/swarm/guard.py``, run by the systemd unit it installs)
replaced an operator script that ran every five minutes outside DELFIN. Three of
its duties did not move over:

  envelope  the script stopped everything at ``load1 > 340`` or
            ``MemAvailable < 150 GB``, read twice in a row. The tick's defaults
            stopped at 1.5 x cores (576 on a 384-core host) and 4 GB -- looser
            than the script on every machine, on the one host that is shared.
  ledger    ``friction.py`` can count what members ran into, but only when a
            person asks; nothing ran it unattended.
  stalls    a member that ended and whose stage did nothing with it, a member
            alive with a log that stopped growing, a real run that outlived its
            thread -- all three were recorded in the run's own ``events.jsonl``,
            which nothing read.

Every case here drives ``guard.tick`` (the entry point the service, the
scheduler entry and ``wave tick`` all reach) with a fake host -- load, memory,
members, logs -- and no case reads or writes the live waves root or the real
slice: the ``host`` fixture points the root at a tmp dir and names the slice
after the test process, as the wave suite does.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

from delfin.agent.swarm import friction, guard, mcp_server as M, wave

_MEMBER = r'''
import json, os, subprocess, sys, time
args = sys.argv[1:]
brief = args[args.index("-p") + 1]
wt = args[args.index("--cwd") + 1]
os.chdir(wt)
for word in brief.split():
    if word.startswith("sleep:"):
        time.sleep(float(word[6:]))
print(json.dumps({"type": "result", "stop_reason": "end_of_turn",
                  "tool_calls": [1, 2, 3]}))
'''


@pytest.fixture
def host(tmp_path, monkeypatch):
    """A waves root of its own, a git repository, stand-in members, and a host
    the test decides: ``set_host`` gives the load and memory the envelope
    sees, ``at`` moves the clock the stall checks read."""
    root = tmp_path / "waves"
    repo = tmp_path / "repo"
    repo.mkdir()
    for args in (["init", "-q", "-b", "main"], ["config", "user.email", "t@example.org"],
                 ["config", "user.name", "t"]):
        subprocess.run(["git", *args], cwd=repo, check=True)
    (repo / "README").write_text("r\n")
    subprocess.run(["git", "add", "."], cwd=repo, check=True)
    subprocess.run(["git", "commit", "-qm", "init"], cwd=repo, check=True)
    member = tmp_path / "member.py"
    member.write_text(_MEMBER)
    # The shipped slice and root belong to the live swarm: an emergency stop
    # in a test must never reach them (measured 2026-11-24).
    cfg = dict(wave._DEFAULTS, waves_root=str(root), min_available_gb=0,
               max_load_per_core=1e9, slice=f"delfin-test-{os.getpid()}.slice")
    monkeypatch.setattr(wave, "settings", lambda: dict(cfg))
    # The envelope limits read agent.swarm through the real _limits(), so the
    # settings are what this dict holds (the shipped defaults unless a case
    # names another value) -- patching the source rather than the reader
    # keeps the retired-key path under test too.
    own: dict = {}
    import delfin.user_settings as _us
    monkeypatch.setattr(_us, "load_settings", lambda: {"agent": {"swarm": own}})
    monkeypatch.setattr(guard, "_systemctl", lambda *a: subprocess.CompletedProcess(a, 0, "", ""))
    monkeypatch.setattr(wave, "_scoped", lambda cmd, cfg: cmd)
    monkeypatch.setattr(wave, "_refuse_on_a_login_node", lambda: None)
    monkeypatch.setattr(wave, "member_argv", lambda run, task, session="", cfg=None, python=None:
                        [sys.executable, str(member), "-p", task["brief"],
                         "--cwd", str(wave.run_dir(run) / task["id"])])
    seen: list = []
    monkeypatch.setattr(guard, "_attention",
                        lambda kind, title, detail: seen.append((kind, title, detail)))
    # The host this test decides. Fake it at BOTH seams: the one the guard
    # reads through (absent on the base, hence raising=False) and the raw
    # system call the base reads directly, so a case says the same thing on
    # either tree instead of depending on the machine it runs on (rule 6).
    reading = {"load1": 1.0, "cores": 384, "avail_gb": 900.0, "total_gb": 1511.3}
    monkeypatch.setattr(guard.os, "getloadavg",
                        lambda: (reading["load1"], reading["load1"], reading["load1"]))
    monkeypatch.setattr(guard.os, "cpu_count", lambda: reading["cores"])

    def meminfo_text() -> str:
        """This host's /proc/meminfo. avail_gb None leaves MemAvailable out,
        which is the host that cannot say -- not the same as 0."""
        text = f"MemTotal:       {int(reading['total_gb'] * 1048576)} kB\n"
        if reading["avail_gb"] is not None:
            text += f"MemAvailable:   {int(reading['avail_gb'] * 1048576)} kB\n"
        return text

    # /proc/meminfo itself is faked, not only the reader around it: the base
    # reads the file directly, so faking only the reader would leave the base
    # judging this host's real memory (rule 6).
    real_path = guard.Path

    class _HostFile:
        """A path like ``Path``, except that /proc/meminfo answers from this
        dict. A plain stand-in rather than a Path subclass: subclassing
        pathlib differs between 3.11 and 3.12 and a test may not care."""

        def __init__(self, *parts):
            self._p = real_path(*parts)

        def read_text(self, *args, **kwargs):
            return (meminfo_text() if str(self._p) == "/proc/meminfo"
                    else self._p.read_text(*args, **kwargs))

        def __getattr__(self, name):
            return getattr(self._p, name)

        def __str__(self):
            return str(self._p)

    monkeypatch.setattr(guard, "Path", lambda *p: _HostFile(*p))
    # Which pids this host says are running. An ended member gets a pid that
    # is not in here: pid 0 would not do, because os.kill(0, 0) signals the
    # whole process group and wave.alive(0) answers True.
    alive: set = set()
    monkeypatch.setattr(wave, "alive", lambda pid: int(pid) in alive)
    # The two host seams _tick reaches through _children/pids, faked beside the
    # others. Without them a case that uses a LIVE pid (os.getpid()) makes the
    # tick read this machine: listeners_beyond_loopback() shells out to the real
    # `ss`, and on a host where this process listens beyond the loopback the
    # tick emergency-stops instead of doing what the case asserts (measured:
    # a fake `ss` line for os.getpid() turns the silent-member case into
    # "a member's process listens beyond the loopback"). A pid no process holds
    # answers the same way on every host.
    monkeypatch.setattr(guard, "_children", lambda pid: {int(pid)})
    monkeypatch.setattr(guard, "listeners_beyond_loopback", lambda pids: [])
    # When a member ended == when its log last changed, as in the real code.
    monkeypatch.setattr(wave, "_ended_at", lambda ref: Path(
        wave.load_state(ref.partition("/")[0])["agents"][ref.partition("/")[2]]["log"]
    ).stat().st_mtime)
    real_time = time.time
    clock = {"now": real_time(), "moving": False}

    def advance(d: float) -> None:
        """Move the guard's clock -- only for the cases that wait on an age.

        The ledger-cost case must NOT move it: there the report is a few
        minutes old and the tick is timed for real, so a moved clock would
        make every pass rewrite the report and measure the wrong thing."""
        clock["moving"] = True
        clock["now"] += d

    monkeypatch.setattr(guard, "_now", lambda: clock["now"], raising=False)
    return {"root": root, "repo": repo, "seen": seen, "limits": own,
            "host": reading, "clock": clock, "advance": advance, "alive": alive}


def _task(tid, brief="work", **kw):
    return {"id": tid, "brief": brief, "base": "main", **kw}


def _agent_entry(run, tid, log: Path, *, started=None):
    """One member in state.json, ended (its pid is one no process holds).

    ``driven`` is set the way ``h_start`` sets it: a run the tick drives is the
    only kind whose ended member the tick judges. A test that leaves it out is
    describing a run another tool started."""
    st = wave.load_state(run)
    st["driven"] = True
    st.setdefault("agents", {})[tid] = {
        "pid": 4242, "branch": f"wave/{run}/{tid}", "worktree": str(wave.run_dir(run) / tid),
        "log": str(log), "started": started if started is not None else time.time(),
        "task": _task(tid)}
    wave._save_state(run, st)


# -- 1. the envelope -------------------------------------------------------------------

def test_the_old_scripts_stop_line_is_the_default_on_a_384_core_host(host):
    """If the envelope is the old one, a 384-core host at load 341 -- which the
    operator script stopped at -- must be outside it. On the shipped 1.5 x
    cores the stop was 576, so load 341 read as healthy."""
    host["host"].update(load1=341.0, cores=384)
    assert guard.outside_envelope() == "load 341 on 384 cores"


def test_a_host_just_under_the_old_stop_line_is_left_alone(host):
    """Control for the line above: 339 on 384 cores is inside the envelope, so
    the test below is the arithmetic and not a guard that stops everything."""
    host["host"].update(load1=339.0, cores=384)
    assert guard.outside_envelope() == ""


def test_the_memory_floor_scales_with_the_host(host):
    """The share decides on a big host (10 % of 1511.3 = 151.1 GB, against the
    script's 150); below 40 GB of total the 4 GB floor the tick already
    shipped decides, so a 16 GB laptop keeps the stop it had.

    Asked through ``outside_envelope`` -- the function the tick calls -- and
    not through the helper that computes the floor: a test that called the
    helper would be red on the base for ``AttributeError`` (the helper did not
    exist) instead of for the bug, which is that 140 GB free on a 1.5 TB host
    was judged healthy."""
    # 140 GB free of 1511.3: under the 151.1 share, far over the flat 4 GB.
    host["host"].update(cores=384, total_gb=1511.3, avail_gb=140.0)
    assert guard.outside_envelope() == "140.0 GB available"
    # The same 140 GB on a 16 GB host is plenty: the floor there is 4 GB.
    host["host"].update(cores=8, total_gb=16.0, avail_gb=14.0)
    assert guard.outside_envelope() == ""


def test_low_available_memory_stops_a_large_host(host):
    """120 GB free on a 1.5 TB host is 8 % -- under the 10 % share, and 30 GB
    under the old script's floor. A flat 4 GB default called this healthy."""
    host["host"].update(cores=384, total_gb=1511.3, avail_gb=120.0)
    assert guard.outside_envelope() == "120.0 GB available"


def test_a_small_host_with_6_gb_free_is_inside_the_envelope(host):
    """Control: 6 GB free of 16 GB is normal on a laptop and must not stop
    anything. An absolute 150 GB floor would have stopped it."""
    host["host"].update(cores=8, total_gb=16.0, avail_gb=6.0)
    assert guard.outside_envelope() == ""


def test_no_memory_left_at_all_is_outside_the_envelope(host):
    """A host with 0 GB available is the one in trouble. `if avail_gb and
    ...` -- truthiness on a reading -- would read it as healthy and fail
    open, which is the wrong direction for a brake."""
    host["host"].update(cores=384, total_gb=1511.3, avail_gb=0.0)
    assert guard.outside_envelope() == "0.0 GB available"


def test_a_host_that_cannot_read_its_memory_is_not_stopped_for_it(host):
    """Control for the case above: no MemAvailable to read at all is unknown,
    not zero, and an unknown number must not pull the brake either."""
    host["host"].update(cores=384, total_gb=1511.3, avail_gb=None)
    assert guard.outside_envelope() == ""


def test_the_envelope_limits_stay_settings(host):
    """Both numbers are overridable in agent.swarm, so a host that wants the
    looser 1.5 per core says so by name."""
    host["host"].update(load1=400.0, cores=384)
    host["limits"]["guard_max_load_per_core"] = 1.5
    assert guard.outside_envelope() == ""
    host["limits"]["guard_max_load_per_core"] = 0.885
    assert guard.outside_envelope().startswith("load 400")


def test_two_ticks_outside_the_new_envelope_stop_everything(host):
    """The duty end to end: no member stops on the first tick outside, the
    second one stops them and writes the hold flag."""
    M.h_start("r1", [_task("a", "sleep:30")], str(host["repo"]))
    host["alive"].add(wave.load_state("r1")["agents"]["a"]["pid"])
    host["host"]["load1"] = 400.0
    assert guard.tick()["stopped"] is None
    out = guard.tick()["stopped"]
    assert out["stopped"] == {"r1": ["a"]}
    assert (host["root"] / "EMERGENCY_STOPPED").exists()
    assert guard.clear_stop() is True


def test_an_in_envelope_host_is_not_stopped(host):
    """Control for the two tests above: the same member, the same ticks, a host
    at 0.41 load per core and 898 GB free runs on."""
    M.h_start("r2", [_task("a", "sleep:30")], str(host["repo"]))
    host["alive"].add(wave.load_state("r2")["agents"]["a"]["pid"])
    host["host"].update(load1=157.0, cores=384, avail_gb=898.0, total_gb=1511.3)
    assert guard.tick()["stopped"] is None
    assert guard.tick()["stopped"] is None
    assert wave.live_members(["r2"]), "the member must not have been stopped"
    wave.stop("r2")


# -- 2. the friction ledger -------------------------------------------------------------

def test_the_tick_counts_friction_without_anyone_asking(host):
    """The ledger is the tick's duty now: a member log that appears is counted
    into <root>/.friction.jsonl by the next tick, with no person and no tool
    call in between."""
    run = "f1"
    wave.run_dir(run).mkdir(parents=True)
    wave._save_state(run, {"agents": {}})
    log = wave.run_dir(run) / "a.log"
    log.write_text(json.dumps({"type": "notice",
                               "text": "Plan mode is on; edits refused"}) + "\n")
    _agent_entry(run, "a", log)
    did = guard.tick()
    assert did["friction"]["appended"] == 1
    recs = friction.read_ledger()
    assert [(r["kind"], r["code"], r["run"], r["id"]) for r in recs] == \
        [("refusal", "plan_mode", run, "a")]
    assert recs[0]["src"].startswith("member:") and recs[0]["src"].endswith(f"{run}/a")


def test_a_second_tick_appends_only_what_is_new(host):
    """If the ledger re-read whole logs, the same denial would be counted
    twice every tick and the counts would mean nothing."""
    run = "f2"
    wave.run_dir(run).mkdir(parents=True)
    wave._save_state(run, {"agents": {}})
    log = wave.run_dir(run) / "a.log"
    line = json.dumps({"type": "permission_denied", "tool": "bash",
                       "reason_code": "allow_list", "text": "not on the auto-allow list"})
    with open(log, "w") as fh:
        fh.write(line + "\n")
    _agent_entry(run, "a", log)
    assert guard.tick()["friction"]["appended"] == 1
    assert guard.tick()["friction"]["appended"] == 0
    with open(log, "a") as fh:
        fh.write(line + "\n")
    assert guard.tick()["friction"]["appended"] == 1
    assert len(friction.read_ledger()) == 2


def test_a_part_line_is_not_counted_until_it_is_whole(host):
    """A stream-json log is read while the member writes it. A half line that
    cannot be classified must be re-read next tick, not skipped forever: the
    cursor stays behind it."""
    run = "f3"
    wave.run_dir(run).mkdir(parents=True)
    wave._save_state(run, {"agents": {}})
    log = wave.run_dir(run) / "a.log"
    whole = json.dumps({"type": "notice", "text": "Plan mode is on"}) + "\n"
    part = '{"type": "notice", "text": "Context com'
    with open(log, "w") as fh:
        fh.write(whole + part)
    _agent_entry(run, "a", log)
    assert guard.tick()["friction"]["appended"] == 1
    with open(log, "a") as fh:
        fh.write('pacted"}' + "\n")
    assert guard.tick()["friction"]["appended"] == 1
    assert [r["code"] for r in friction.read_ledger()] == ["plan_mode", "compacted"]


def test_the_tick_rewrites_the_report_at_most_hourly(host):
    """The grouped report re-reads the whole ledger; the tick writes it once an
    hour and the ledger every tick, so a busy root does not pay for a report
    every minute."""
    run = "f4"
    wave.run_dir(run).mkdir(parents=True)
    wave._save_state(run, {"agents": {}})
    log = wave.run_dir(run) / "a.log"
    _agent_entry(run, "a", log)
    line = json.dumps({"type": "notice", "text": "Plan mode is on"}) + "\n"
    report = host["root"] / "friction.json"

    with open(log, "w") as fh:
        fh.write(line)
    guard.tick()
    first = json.loads(report.read_text())
    assert first["records"] == 1

    with open(log, "a") as fh:
        fh.write(line)
    host["advance"](120)
    guard.tick()
    assert len(friction.read_ledger()) == 2, "the ledger keeps up every tick"
    assert json.loads(report.read_text())["records"] == 1, "the report waits its hour"

    host["advance"](3700)
    guard.tick()
    assert json.loads(report.read_text())["records"] == 2


def test_the_ledger_cost_stays_small(host):
    """What the tick pays per pass, measured: reading each log from where the
    last pass stopped, not from the top. On a root with a finished 400 KB log
    the second pass must not read those bytes again."""
    run = "f5"
    wave.run_dir(run).mkdir(parents=True)
    wave._save_state(run, {"agents": {}})
    log = wave.run_dir(run) / "a.log"
    line = json.dumps({"type": "notice", "text": "waiting for a slot on kit"}) + "\n"
    with open(log, "w") as fh:
        fh.write(line * 4000)                        # ~330 KB of settled history
    _agent_entry(run, "a", log)
    t0 = time.perf_counter()
    first = guard.tick()["friction"]
    cold = time.perf_counter() - t0
    t0 = time.perf_counter()
    second = guard.tick()["friction"]
    warm = time.perf_counter() - t0
    assert first["appended"] == 4000 and second["appended"] == 0
    assert second["sources"] == first["sources"]
    assert warm < cold, f"a settled root must cost less than a first read ({warm:.4f} vs {cold:.4f})"


def test_an_old_absolute_memory_override_still_decides(host, monkeypatch):
    """A host that set the retired ``guard_min_available_gb`` keeps it. A
    changed default is a changed meaning for every caller: silently dropping
    that key would hand a host a looser stop than the one it asked for."""
    import delfin.user_settings as us
    monkeypatch.setattr(us, "load_settings",
                        lambda: {"agent": {"swarm": {"guard_min_available_gb": 200}}})
    host["host"].update(cores=384, total_gb=1511.3, avail_gb=160.0)
    assert guard._limits()["guard_min_memory_floor_gb"] == 200
    assert guard.outside_envelope() == "160.0 GB available"
    # ...and the new key wins over the retired one when both are present.
    monkeypatch.setattr(us, "load_settings",
                        lambda: {"agent": {"swarm": {"guard_min_available_gb": 200,
                                                     "guard_min_memory_floor_gb": 4}}})
    assert guard.outside_envelope() == ""


def test_the_wave_tick_command_keeps_the_ledger(host, capsys):
    """Reachability: the same duty through the command a person or a timer
    types, not only through the Python function. If ``wave tick`` did not
    reach the ledger, the counts would exist only for a caller inside
    DELFIN -- the defect that left friction uncounted in the first place."""
    from delfin.agent import cli

    run = "f6"
    wave.run_dir(run).mkdir(parents=True)
    wave._save_state(run, {"agents": {}})
    log = wave.run_dir(run) / "a.log"
    log.write_text(json.dumps({"type": "notice", "text": "Plan mode is on"}) + "\n")
    _agent_entry(run, "a", log)
    assert cli.main(["swarm", "wave", "tick"]) == 0
    assert json.loads(capsys.readouterr().out)["friction"]["appended"] == 1
    assert len(friction.read_ledger()) == 1


def test_the_tick_service_line_says_what_it_found(host, monkeypatch):
    """The service's one line per tick is what a person reads in the journal:
    a stall or a count that never reaches it is invisible on a host where
    nobody opens the waves root."""
    run = "f7"
    wave.run_dir(run).mkdir(parents=True)
    (wave.run_dir(run) / "tasks.json").write_text(json.dumps([_task("a")]))
    log = wave.run_dir(run) / "a.log"
    log.write_text('{"type": "result", "stop_reason": "end_of_turn", "tool_calls": 1}\n')
    _agent_entry(run, "a", log)
    lines: list = []
    monkeypatch.setattr(guard, "_attention", lambda *a, **k: None)
    stamp = time.time()
    os.utime(log, (stamp, stamp))
    host["advance"](20 * 60)                            # ended 20 min ago
    guard.loop(60, ticks=1, sleep=lambda s: None, log=lines.append)
    row = json.loads(lines[-1])
    assert row["stalls"] >= 1, f"the journal line must carry the stall: {row}"
    assert "friction" in row


# -- 3. stalls --------------------------------------------------------------------------

def test_a_member_ended_20_minutes_ago_with_nothing_done_is_raised_once(host):
    """The member ended; its stage has no branch left to run (it is not final,
    so nothing more will happen to it). Nobody said so. 15 minutes is past the
    stage's 900 s endpoint-error pause, so this is not a member about to be
    resumed."""
    run = "s1"
    wave.run_dir(run).mkdir(parents=True)
    (wave.run_dir(run) / "tasks.json").write_text(json.dumps([_task("a")]))
    log = wave.run_dir(run) / "a.log"
    log.write_text(json.dumps({"type": "result", "stop_reason": "end_of_turn",
                               "tool_calls": 3}) + "\n")
    _agent_entry(run, "a", log)
    # The log keeps its real mtime and the guard's clock moves past it: the
    # member ended 20 minutes ago, measured on the clock the guard reads.
    stamp = time.time()
    os.utime(log, (stamp, stamp))
    host["advance"](20 * 60)
    assert guard.tick()["stalls"] == 1
    kind, title, detail = host["seen"][0]
    assert kind == "session_state" and "s1/a" in title and "ended" in title
    assert "stalled (" in title, f"the tick's own wording, not observe()'s: {title}"
    assert json.loads(detail)["min"] >= 19
    assert "s1/a:ended" in guard._load()["stalls"]
    del host["seen"][:]
    host["advance"](600)
    assert guard.tick().get("stalls") in (None, 0)
    assert host["seen"] == [], "the same stall must not be raised twice"


def test_a_member_that_ended_a_minute_ago_is_not_a_stall(host):
    """Control: a member one minute past its end is mid-stage -- the stage may
    still be resuming it. Nothing goes to the inbox."""
    run = "s2"
    wave.run_dir(run).mkdir(parents=True)
    (wave.run_dir(run) / "tasks.json").write_text(json.dumps([_task("a")]))
    log = wave.run_dir(run) / "a.log"
    log.write_text(json.dumps({"type": "result", "stop_reason": "end_of_turn",
                               "tool_calls": 3}) + "\n")
    _agent_entry(run, "a", log)
    guard.tick()
    assert host["seen"] == []


def test_a_member_whose_real_run_marker_outlived_its_process_is_raised(host):
    """The stage keeps a member out of every later pass while
    ``<tid>.realrun_running`` exists. If the process named there died, the
    marker holds the member forever and the wave never finishes."""
    run = "s3"
    wave.run_dir(run).mkdir(parents=True)
    (wave.run_dir(run) / "tasks.json").write_text(json.dumps([_task("a")]))
    log = wave.run_dir(run) / "a.log"
    log.write_text('{"type": "result", "stop_reason": "end_of_turn", "tool_calls": 1}\n')
    _agent_entry(run, "a", log)
    marker = wave.run_dir(run) / "a.realrun_running"
    marker.write_text("999999")                      # a pid that is not running
    assert guard.tick()["stalls"] == 1
    assert "real run" in host["seen"][0][2]


def test_a_live_member_with_a_silent_log_is_raised(host):
    """Alive, and nothing written for over an hour.

    ``started`` is as old as the log on purpose: a member resumed seconds ago
    with an hour-old log is NOT this case (see
    ``test_a_member_the_stage_just_resumed_is_not_silent``), so the session must
    have begun as long ago as its silence to be a stall."""
    run = "s4"
    wave.run_dir(run).mkdir(parents=True)
    (wave.run_dir(run) / "tasks.json").write_text(json.dumps([_task("a")]))
    log = wave.run_dir(run) / "a.log"
    log.write_text("{}\n")
    old = time.time() - 90 * 60
    st = wave.load_state(run)
    host["alive"].add(4242)                                     # alive
    st["agents"] = {"a": {"pid": 4242,
                          "branch": f"wave/{run}/a",
                          "worktree": str(wave.run_dir(run) / "a"),
                          "log": str(log), "started": old,
                          "task": _task("a")}}
    wave._save_state(run, st)
    os.utime(log, (old, old))
    assert guard.tick()["stalls"] == 1
    assert "stalled (silent)" in host["seen"][0][1]


def test_a_task_that_was_never_started_is_not_a_stall(host):
    """Control: a queued task has no member and no log at all. It waits for
    what it runs after; it is not stuck."""
    run = "s5"
    wave.run_dir(run).mkdir(parents=True)
    (wave.run_dir(run) / "tasks.json").write_text(json.dumps([_task("a")]))
    wave._save_state(run, {"agents": {}, "driven": True})
    guard.tick()
    assert host["seen"] == []


def test_a_member_whose_stage_already_acted_is_not_a_stall(host):
    """Control: the member ended 20 minutes ago, but its stage DID act on it --
    it asked for a real run and is waiting for one. "Ended with nothing done by
    its stage" is false, and the inbox must stay empty. Without this check the
    wave raises a stall for every member it is correctly waiting on."""
    run = "s6"
    wave.run_dir(run).mkdir(parents=True)
    (wave.run_dir(run) / "tasks.json").write_text(json.dumps([_task("a")]))
    log = wave.run_dir(run) / "a.log"
    log.write_text(json.dumps({"type": "result", "stop_reason": "end_of_turn",
                               "tool_calls": 3}) + "\n")
    _agent_entry(run, "a", log)
    (wave.run_dir(run) / "a.asked").write_text("asked for a real run\n")
    stamp = time.time()
    os.utime(log, (stamp, stamp))
    host["advance"](20 * 60)
    guard.tick()
    assert host["seen"] == []


def test_a_run_another_tool_started_gets_no_ended_stall(host):
    """Control: the root may hold runs this tick does not drive -- the start
    block guards every action it takes with ``driven``, and "ended and its
    stage did nothing" is the same kind of claim: a member that ended is the
    stage's to resume, and a foreign run keeps its own stage.

    The gate is on ``ended`` only, and deliberately so: ``silent`` and
    ``realrun`` describe a process that is already gone or wedged in THIS
    root, which is true whoever started it. This test therefore pins one
    kind, not a blanket silence -- a foreign run whose member is still alive
    and whose log stopped growing is still raised as ``silent``, as it is on
    origin/main."""
    run = "s7"
    wave.run_dir(run).mkdir(parents=True)
    (wave.run_dir(run) / "tasks.json").write_text(json.dumps([_task("a")]))
    log = wave.run_dir(run) / "a.log"
    log.write_text(json.dumps({"type": "result", "stop_reason": "end_of_turn",
                               "tool_calls": 3}) + "\n")
    _agent_entry(run, "a", log)
    st = wave.load_state(run)
    st.pop("driven", None)                       # another tool started it
    wave._save_state(run, st)
    stamp = time.time()
    os.utime(log, (stamp, stamp))
    host["advance"](20 * 60)
    guard.tick()
    assert host["seen"] == []


def test_a_foreign_runs_live_member_is_still_called_silent(host):
    """The other half of the gate above, so its scope cannot drift: ``ended``
    is withheld from a run this tick does not drive, ``silent`` is not. A
    member that is still alive with a log that stopped growing is a wedged
    process in this root whoever started it, and origin/main raises it too --
    this is the case that would fail if the ``driven`` gate were widened to
    swallow every kind."""
    run = "s8"
    wave.run_dir(run).mkdir(parents=True)
    (wave.run_dir(run) / "tasks.json").write_text(json.dumps([_task("a")]))
    log = wave.run_dir(run) / "a.log"
    log.write_text("")
    stamp = time.time() - 2 * 3600
    os.utime(log, (stamp, stamp))
    _agent_entry(run, "a", log, started=stamp)     # started long ago too
    st = wave.load_state(run)
    st.pop("driven", None)
    st["agents"]["a"]["pid"] = 777
    wave._save_state(run, st)
    host["alive"].add(777)                         # the host says it runs
    host["advance"](120 * 60)
    guard.tick()
    assert len(host["seen"]) == 1, host["seen"]
    kind, title, _ = host["seen"][0]
    # The behaviour, not the wording: one session_state alert about this
    # member's silence. The title differs between the two trees ("silent" on
    # origin/main, "stalled (silent)" here), and a control that pins the
    # wording would be green on one tree only.
    assert kind == "session_state", kind
    assert "s8/a" in title and "silent" in title.lower(), title


def test_a_member_the_stage_just_resumed_is_not_silent(host):
    """The stage resumes a member that fell through its branches, and a resume
    keeps the old log: ``wave._launch`` opens it "ab", which does not touch its
    mtime, and gives the member a fresh pid and ``started``. Judged on the log
    alone the member reads as silent for as long as its previous session ran --
    measured on the arm's code: a member with ``started`` = now and a
    90-minute-old log was reported "silent, 90 min" in the tick right after its
    resume, i.e. a false alarm that arrives exactly when the wave is
    recovering."""
    run = "s8"
    wave.run_dir(run).mkdir(parents=True)
    (wave.run_dir(run) / "tasks.json").write_text(json.dumps([_task("a")]))
    log = wave.run_dir(run) / "a.log"
    log.write_text(json.dumps({"type": "result", "stop_reason": "end_of_turn",
                               "tool_calls": 3}) + "\n")
    _agent_entry(run, "a", log, started=time.time())   # resumed seconds ago
    old = time.time() - 90 * 60
    os.utime(log, (old, old))                          # the old session's log
    host["alive"].add(4242)                            # and it is running now
    guard.tick()
    assert host["seen"] == []


def test_a_member_silent_past_the_old_session_is_still_raised(host):
    """Control for the case above: the start time must not hide a member that
    really is quiet. Alive, its log 90 minutes old AND started 90 minutes ago
    -- nothing has happened since the session began -- is a stall."""
    run = "s9"
    wave.run_dir(run).mkdir(parents=True)
    (wave.run_dir(run) / "tasks.json").write_text(json.dumps([_task("a")]))
    log = wave.run_dir(run) / "a.log"
    log.write_text("{}\n")
    old = time.time() - 90 * 60
    os.utime(log, (old, old))
    _agent_entry(run, "a", log, started=old)
    host["alive"].add(4242)
    assert guard.tick()["stalls"] == 1
    assert "stalled (silent)" in host["seen"][0][1]


def test_a_stall_that_clears_is_forgotten_so_it_can_be_raised_again(host):
    """``raise_stalls`` documents that it "forgets the ones no longer seen so
    they can be raised again".

    The forget loop only ran when ``found`` was non-empty (``if found:``
    around the call in ``tick``), so on the very tick where the last stall
    clears -- ``found == []``, the only tick that can forget the last one --
    it never ran. A member that stalled, recovered and stalled again was
    therefore never reported the second time: the key stayed in the state
    and deduplicated it forever.
    """
    run = "s10"
    wave.run_dir(run).mkdir(parents=True)
    (wave.run_dir(run) / "tasks.json").write_text(json.dumps([_task("a")]))
    log = wave.run_dir(run) / "a.log"
    log.write_text("{}\n")
    old = time.time() - 90 * 60
    os.utime(log, (old, old))
    _agent_entry(run, "a", log, started=old)
    host["alive"].add(4242)

    # 1. it is silent past the 60-minute line: raised, and remembered.
    assert guard.tick()["stalls"] == 1
    assert "s10/a:silent" in guard._load()["stalls"]

    # 2. it writes to its log again and the guard's clock moves on: the stall
    #    is gone, so the tick forgets it.
    stamp = time.time()
    os.utime(log, (stamp, stamp))
    state = wave.load_state(run)
    state["agents"]["a"]["started"] = stamp
    wave._save_state(run, state)
    host["advance"](5 * 60)
    assert guard.tick().get("stalls") in (None, 0)
    assert "s10/a:silent" not in guard._load()["stalls"], (
        "a cleared stall must be forgotten, or it can never be raised again")

    # 3. it goes quiet again: the same stall is raised a second time.
    host["advance"](90 * 60)
    assert guard.tick()["stalls"] == 1, (
        "the second stall of the same member must be reported, not "
        "deduplicated against the first")
    assert "s10/a:silent" in guard._load()["stalls"]


def test_the_start_gate_stops_before_the_envelope_does(host):
    """A host must stop launching members before it stops everything.

    The envelope stops the whole slice (``emergency_stop``: every member of
    every run, including a foreign tool's) once it is outside for two ticks.
    ``wave._room`` decides how many members start, and it has its own load
    number on purpose -- a different question, as the guard's constants say.
    But if the start gate stays open past the stop line, there is a band of
    host load in which the tick keeps launching members onto a host the
    envelope has already condemned, and then stops them all. On base the
    order was the other way (stop 1.5 x cores, start 1.0 x cores); the 0.885
    stop line inverted it, and nothing pinned the relation.

    Not pinned here: any particular pair of numbers. Only the order that
    makes the two knobs safe together.
    """
    stop = guard._GUARD_DEFAULTS["guard_max_load_per_core"]
    start = wave._DEFAULTS["max_load_per_core"]
    assert start <= stop, (
        f"wave starts members up to load {start} x cores while the envelope "
        f"stops everything above {stop} x cores: between the two the tick "
        f"launches onto a host it is about to stop. Keep the start gate at or "
        f"below the stop line.")


def test_the_ledger_cursor_and_report_are_owner_only(host):
    """The friction files quote member logs, so they are not world-readable.

    ``.friction.jsonl`` records snippets of what each member's log said, and
    it sits beside the wave state -- which on the cluster is a shared path.
    The three files were written with the process umask (0644 on a normal
    host), so any other account could read them. DELFIN's other writers of
    this kind (attention.py, audit_log.py, approval_code.py) create 0o600.
    """
    import stat

    from delfin.agent.swarm import friction

    run = "s7"
    wave.run_dir(run).mkdir(parents=True)
    (wave.run_dir(run) / "tasks.json").write_text(json.dumps([_task("a")]))
    log = wave.run_dir(run) / "a.log"
    log.write_text(json.dumps({"type": "result", "stop_reason": "end_of_turn",
                               "tool_calls": 1}) + "\n")
    _agent_entry(run, "a", log)
    host["alive"].add(4242)

    friction.sweep()
    friction.write_report()
    for path in (friction._ledger_path(), friction._cursor_path(),
                 wave.root() / "friction.json"):
        assert path.exists(), f"{path} was not written"
        mode = stat.S_IMODE(path.stat().st_mode)
        assert mode == 0o600, f"{path.name} is {oct(mode)}, not 0o600"


def test_the_ledger_is_bounded_so_a_long_lived_root_cannot_grow_it(host):
    """The ledger only appended; nothing removed anything.

    Each tick appends the lines its members wrote since the last one, and
    ``write_report`` re-reads the whole file -- on a host running waves for
    months both the file and the per-report cost grew without bound. The
    report path keeps the newest lines and says how many it dropped.
    """
    from delfin.agent.swarm import friction

    led = friction._ledger_path()
    led.parent.mkdir(parents=True, exist_ok=True)
    keep = 10
    total = 25
    led.write_text("\n".join(
        json.dumps({"run": "r", "id": f"m{i}", "src": "member:r/m",
                    "kind": "end", "code": "ok", "detail": "", "model": ""})
        for i in range(total)) + "\n")

    # The oldest records go, the newest stay, and the caller is told.
    dropped = friction._prune_ledger(keep=keep)
    assert dropped == total - keep
    lines = [line for line in led.read_text().splitlines() if line.strip()]
    assert len(lines) == keep, f"the ledger is not bounded: {len(lines)} lines"
    tail = json.loads(lines[-1])
    assert tail["id"] == f"m{total - 1}", "the NEWEST records must be kept"

    # Below the cap it is a no-op, so a normal root keeps everything.
    assert friction._prune_ledger(keep=keep) == 0
    assert len([l for l in led.read_text().splitlines() if l.strip()]) == keep


def test_the_guard_prose_does_not_claim_a_start_gate_above_the_stop_line():
    """The guard's own comment must describe the numbers the code ships.

    A real run of the swarm tick read `guard.py`'s header and repeated it as
    fact: that `wave._room` holds "its own max_load_per_core (default 1.0)"
    so that the start gate sits above the stop line. That was true before
    the start gate was brought down to the envelope's stop line; after it,
    both constants are 0.885 and the comment was the only thing still
    claiming otherwise. A comment that contradicts its own module is worse
    than no comment: it is read as documentation.
    """
    import inspect
    import re

    from delfin.agent.swarm import guard, wave

    header = inspect.getsource(guard)
    start = wave._DEFAULTS["max_load_per_core"]
    stop = guard._GUARD_DEFAULTS["guard_max_load_per_core"]

    # The copied-verbatim falsehood: a literal 1.0 attributed to wave.
    assert "default 1.0" not in header, (
        "guard.py still says wave's start gate defaults to 1.0, but it is "
        f"{start}"
    )
    # And whatever the prose says about the two numbers, the code agrees:
    # the start gate is at the stop line, not above it.
    assert start <= stop, f"start gate {start} is above the stop line {stop}"
    assert start == stop, (
        "the header argues the two knobs are deliberately different, but "
        f"they are now identical ({start}): either restore the argument or "
        "state the new relationship"
    )
    # The sentence that carried the falsehood named a stale value.
    m = re.search(r"``max_load_per_core`` \(default ([0-9.]+)\)", header)
    if m:
        assert float(m.group(1)) == start, (
            f"guard.py names wave's default as {m.group(1)}, but it is {start}"
        )
