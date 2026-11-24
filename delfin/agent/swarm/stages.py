"""The stage after a member ends: its real runs, the resume that reads them,
and the mark that it is final.

For each member of a run whose session has ended and that is not final:

  0. a member that ended on an endpoint error without writing PR.md is
     resumed after a pause, up to three times (``<tid>.error_resumes``);
  1. a reviewer (``no_realrun``) is final at once;
  2. a member that already had its real runs read (``<tid>.realrun``) is
     final ("after its real run");
  3. a member that wrote ``scratch/realrun.json`` gets them
     (``swarm/realrun.py``, outside its cage) and is resumed with a note
     pointing at ``scratch/realrun/summary.md``;
  4. a member that ended without asking is resumed once with a note asking
     for them (``<tid>.asked``); ending a second time without asking makes
     it final, and says so.

Each step leaves a marker next to the member's log, so a stage that was
restarted continues where it stopped and never repeats a step. A later
task's ``after: ["run/tid.final"]`` waits for the third marker.
"""
from __future__ import annotations

import os
import threading
import time
from pathlib import Path

from . import realrun, wave

MAX_PARALLEL_RUNS = 6
#: A member whose session ended on an endpoint error (a 503, a quota, a
#: timeout the retries did not outlast) and wrote no PR.md is resumed after
#: this pause, at most ``MAX_ERROR_RESUMES`` times: a provider outage is not
#: the member's result (2026-11-24: Blablador answered 503 for two models).
ERROR_RESUME_AFTER_S = 900
MAX_ERROR_RESUMES = 3

_READ_YOUR_RUNS = """
Your real runs are done: scratch/realrun/summary.md, and one log per probe
and side in scratch/realrun/. They were made outside your cage, on your
base and on your branch HEAD, same prompt and model, a fresh HOME each.

Judge them against the prediction you wrote before they ran. If the
branch is not better than the base where your task says it must be, say
so in the first line of scratch/PR.md -- that is a result, not a failure
to hide. If a run shows a defect, fix it, commit, and say that the fix has
not had a real run of its own. Then finish scratch/PR.md with the real
runs table, what it shows, and every suite you ran with its numbers. Do
not write scratch/realrun.json again.
"""

_ASK_FOR_RUNS = """
You ended without scratch/realrun.json. Nothing lands without a real run.
Write it now: up to three probes, each a prompt a person would type that
shows your change working ({"probes": [{"name": "...", "prompt": "...",
"workspace": "" or "repo", "turns": 20}]}). Write your prediction for each
probe into scratch/PR.md first, then end your turn: the runs are made for
you, and you are resumed to read them.
"""


def _in_flight(run: str, tid: str) -> bool:
    """A real run another process started and that still runs: a tick from
    a timer is a process of its own, and its in-memory ``busy`` is gone by
    the next one."""
    p = wave.run_dir(run) / f"{tid}.realrun_running"
    try:
        return wave.alive(int(p.read_text().strip()))
    except (OSError, ValueError):
        return False


def _mark(run: str, tid: str, kind: str, text: str = "") -> None:
    (wave.run_dir(run) / f"{tid}.{kind}").write_text(text)


def advance(run: str, busy: set | None = None, *, start_real_run=None) -> list[dict]:
    """One pass over the run; returns what it did. ``start_real_run(run,
    tid)`` starts the real runs and the resume after them (default: a
    thread); ``busy`` holds the members whose real runs are in flight."""
    busy = busy if busy is not None else set()
    d = wave.run_dir(run)
    st = wave.load_state(run)
    try:
        import json
        tasks = json.loads((d / "tasks.json").read_text())
    except (OSError, ValueError):
        tasks = [a["task"] for a in st["agents"].values()]
    done = []
    for task in tasks:
        tid = task["id"]
        if (d / f"{tid}.final").exists():
            continue
        a = st["agents"].get(tid)
        if a is None or wave.alive(a["pid"]) or tid in busy or _in_flight(run, tid):
            continue
        wt = Path(a["worktree"])
        res = wave._last_result(Path(a["log"]))
        if res.get("stop_reason") == "error" and not (wt / "scratch" / "PR.md").is_file():
            counter = d / f"{tid}.error_resumes"
            n = int(counter.read_text() or 0) if counter.exists() else 0
            try:
                quiet_s = time.time() - Path(a["log"]).stat().st_mtime
            except OSError:
                quiet_s = 0.0
            if n < MAX_ERROR_RESUMES:
                if quiet_s >= ERROR_RESUME_AFTER_S:
                    counter.write_text(str(n + 1))
                    wave.resume(run, tid, "Your session ended on an endpoint error "
                                f"({(res.get('error') or '')[:200]}). Continue where you stopped.")
                    done.append({"id": tid, "did": "resumed after an endpoint error", "n": n + 1})
                continue
        if task.get("no_realrun"):
            _mark(run, tid, "final", "ended; no real run asked of it\n")
            done.append({"id": tid, "did": "final", "why": "no real run asked of it"})
        elif (d / f"{tid}.realrun").exists():
            _mark(run, tid, "final", "after its real run\n")
            done.append({"id": tid, "did": "final", "why": "after its real run"})
        elif (wt / "scratch" / "realrun.json").is_file():
            busy.add(tid)
            # Written before the run starts, so no other tick can start it twice.
            (d / f"{tid}.realrun_running").write_text(str(os.getpid()))
            (start_real_run or _start_real_run)(run, tid, busy)
            done.append({"id": tid, "did": "real run"})
        elif not (d / f"{tid}.asked").exists():
            _mark(run, tid, "asked")
            wave.resume(run, tid, _ASK_FOR_RUNS)
            done.append({"id": tid, "did": "asked for a real run"})
        else:
            _mark(run, tid, "final", "ended twice without asking for a real run\n")
            done.append({"id": tid, "did": "final", "why": "ended twice without asking"})
    return done


_gate = threading.Semaphore(MAX_PARALLEL_RUNS)


def _real_run_then_resume(run: str, tid: str, busy: set) -> None:
    running = wave.run_dir(run) / f"{tid}.realrun_running"
    try:
        with _gate:
            try:
                res = realrun.run(run, tid)
                note = _READ_YOUR_RUNS
                rc = "0"
            except Exception as exc:          # a refused request is the member's to fix
                note = (f"Your scratch/realrun.json could not be run: {exc}. Fix it and end "
                        "your turn; the runs are made for you.")
                rc = "1"
                res = None
        if res is not None:
            _mark(run, tid, "realrun", rc)
        else:
            # A request that cannot be run is the member's to fix -- twice.
            failed = wave.run_dir(run) / f"{tid}.realrun_refused"
            n = int(failed.read_text() or 0) + 1 if failed.exists() else 1
            failed.write_text(str(n))
            if n >= 2:
                _mark(run, tid, "final", "its real-run request could not be run twice\n")
                return
        wave.resume(run, tid, note)
    finally:
        running.unlink(missing_ok=True)
        busy.discard(tid)


def _start_real_run(run: str, tid: str, busy: set) -> None:
    threading.Thread(target=_real_run_then_resume, args=(run, tid, busy), daemon=False).start()


def serve(run: str, *, every_s: float = 60.0, until_final: bool = True) -> int:
    """Advance the run every ``every_s`` until every task is final."""
    busy: set = set()
    while True:
        try:
            advance(run, busy)
        except wave.WaveError:
            pass
        st = wave.status(run)
        if until_final and st["members"] and all("final" in m["marks"] for m in st["members"]) \
                and not busy:
            return 0
        time.sleep(every_s)
