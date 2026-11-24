"""Two fast background jobs are both announced with the exit code they had.

WHY THIS FILE EXISTS. Two independent real runs, on the base and on the
branch, each started two quick ``echo`` jobs and reported the SAME anomaly
without being asked to look for it: the completion notice said

    [background] job bee9bd0c (Print two in background) finished with exit
    code None. Read its output with bash_output before relying on it.

while ``bash_output`` on that very job reported ``exit_code: 0``. The
notice is what a waking session reasons from; "exit code None" is not a
weaker statement than "0", it is a WRONG one -- it says the status could
not be recovered when in fact the job had exited normally and the process
that started it was holding the real code.

The cause is in ``drain_finished_events``. A record with no ``finished_at``
whose pid is no longer alive takes the orphan branch:

    finished_at = rec["finished_at"] = now

and the event is then built with ``"exit_code": rec.get("exit_code")``.
Between a fast job's child dying and its watchdog thread reaching
``_update_job_record(... exit_code=rc_final, finished_at=...)`` there is a
window in which the pid is gone and the code is not yet written. A drain
that lands in that window announces the code as None and claims the event,
and the real code arrives in the record immediately afterwards -- too late,
the notice has already been composed from the stale read.

So the two statements contradict each other on a path the user reads. This
file measures the notice, not ``bash_output``, because the notice is the
thing that was wrong.
"""

from __future__ import annotations

import json
import tempfile
import time
from pathlib import Path

import pytest

from delfin.agent import bash_jobs
from delfin.agent.api_client import KitToolPermissions, _doc_executor


def _call(tool, args, perms):
    raw = str(_doc_executor.execute(tool, args, perms))
    body = raw
    if body.lstrip().startswith("[UNTRUSTED"):
        body = "\n".join(l for l in body.split("\n")
                         if not l.startswith(("[UNTRUSTED", "[END")))
    try:
        return json.loads(body.strip()), raw
    except Exception:
        return {"_raw": raw}, raw


@pytest.fixture
def box():
    with tempfile.TemporaryDirectory(prefix="f2-nonexit-") as d:
        yield Path(d)


@pytest.fixture
def perms(box):
    return KitToolPermissions(mode="bypassPermissions", workspace=str(box))


def _ids(out):
    found = []
    for key in ("job_id", "id"):
        if isinstance(out.get(key), str):
            found.append(out[key])
    for key in ("jobs", "started"):
        for entry in out.get(key) or []:
            if isinstance(entry, dict):
                for k in ("job_id", "id"):
                    if isinstance(entry.get(k), str):
                        found.append(entry[k])
    return found


def test_both_of_two_fast_jobs_are_announced_with_their_real_exit_code(
        box, perms):
    """The notice must not say None for a job that exited 0.

    The input that makes this fail is two ``echo`` jobs started back to
    back and drained at a turn boundary: the second one is still in the
    child-died-code-not-yet-written window when the drain runs.
    """
    started = []
    # Six, not two: the window between a child dying and its watchdog
    # writing the code is short, and two jobs missed it in every attempt
    # while six found it in 6 of 6 (11 of 36 jobs announced None before
    # the fix). A test that races once proves nothing on the day it wins.
    for i in range(6):
        out, raw = _call("bash_background",
                         {"command": f"echo {i}",
                          "description": f"Print {i} in background"},
                         perms)
        ids = _ids(out)
        assert ids, f"bash_background started no job: {raw[:400]}"
        started.append(ids[0])

    # A turn boundary: the drain that composes the notice. It CLAIMS what it
    # returns, so a second drain does not repeat it -- the loop therefore
    # accumulates instead of overwriting what an earlier drain gave.
    deadline = time.time() + 30.0
    by_id = {}
    while time.time() < deadline and len(by_id) < len(started):
        for e in bash_jobs.drain_finished_events(str(box)) or ():
            by_id[e.get("job_id")] = e
        if len(by_id) >= len(started):
            break
        time.sleep(0.005 if not by_id else 0.05)

    assert len(by_id) >= len(started), (
        f"expected both jobs to be announced, got {len(by_id)}: "
        f"{sorted(by_id)}")
    for jid in started:
        assert jid in by_id, (
            f"job {jid} started through the dispatch was never announced; "
            f"announced: {sorted(by_id)}")
        announced = by_id[jid].get("exit_code")
        # The literal code is host-dependent (a sandbox that refuses the
        # command exits 126, a plain host 0), so this asserts the NOTICE
        # against the JOB's own record rather than against a number. What
        # makes it fail is the defect: announced None while the record
        # holds the real code.
        reported, raw_status = _call("bash_output", {"job_id": jid}, perms)
        assert announced == reported.get("exit_code"), (
            f"job {jid} was announced with exit_code={announced!r} while "
            f"bash_output reports {reported.get('exit_code')!r}: the notice "
            f"a waking session reasons from contradicts the job's own "
            f"record ({raw_status[:300]})")
        assert announced is not None, (
            f"job {jid} was announced with exit_code=None -- read as "
            f"'status unknown' -- although it had ended and its record "
            f"held a code")


def test_the_notice_never_says_none_when_the_process_holds_a_code(
        box, perms):
    """Deterministic form of the race above, which only lands half the time.

    The window (child reaped, watchdog not yet written) is real but short:
    on the unfixed tree this file's first test is red in only 6 of 12 runs,
    which is a weak red. So the state is CONSTRUCTED instead of raced: a
    real job is started and really finishes, and then the persisted record
    is rewound to what it looks like in that window -- ``finished_at`` and
    ``exit_code`` cleared, the pid already dead -- while this process still
    holds the code in memory. That is exactly the state the defect
    mishandled, reached without waiting for it.

    What makes this fail is the announced ``None``: the drain must ask the
    process that started the job rather than announce a status it has.
    """
    out, raw = _call("bash_background",
                     {"command": "echo hi; exit 3",
                      "description": "a job that ends non-zero"}, perms)
    ids = _ids(out)
    assert ids, f"bash_background started no job: {raw[:400]}"
    jid = ids[0]

    # Let it really finish, so the in-memory job knows its code.
    known = None
    for _ in range(600):
        job = bash_jobs._REGISTRY._jobs.get(jid)
        if job is not None:
            known = job.poll()
            if known is not None:
                break
        time.sleep(0.02)
    assert known is not None, (
        f"job {jid} never finished in this process, so its code cannot be "
        f"held; the test cannot construct the window")

    # Rewind the persisted record to the mid-window state.
    reg = bash_jobs._load_registry_file(str(box))
    rec = reg["jobs"][jid]
    rec["exit_code"] = None
    rec["finished_at"] = None
    rec["claimed_at"] = 0.0
    rec["acknowledged"] = False
    bash_jobs._atomic_write_json(bash_jobs._registry_path(str(box)), reg)

    events = [e for e in bash_jobs.drain_finished_events(str(box)) or ()
              if e.get("job_id") == jid]

    # Either the drain waits (no event yet, which is the honest "not known
    # yet") or it reports the code it holds -- never None.
    for e in events:
        assert e.get("exit_code") == known, (
            f"the drain announced exit_code={e.get('exit_code')!r} for a job "
            f"this process knows ended {known}: a notice that says the status "
            f"is unknown while it is in hand")


def test_the_announced_code_and_the_recorded_code_agree(box, perms):
    """The notice and ``bash_output`` may not contradict each other.

    The control for the test above: it asserts the relationship, not a
    literal, so it also holds for a non-zero code.
    """
    out, raw = _call("bash_background",
                     {"command": "echo hi; exit 3",
                      "description": "Fail with 3"}, perms)
    ids = _ids(out)
    assert ids, f"bash_background started no job: {raw[:400]}"
    jid = ids[0]

    deadline = time.time() + 30.0
    event = None
    while time.time() < deadline:
        for e in bash_jobs.drain_finished_events(str(box)) or ():
            if e.get("job_id") == jid:
                event = e
                break
        if event is not None:
            break
        time.sleep(0.01)

    assert event is not None, f"job {jid} was never announced"
    reported, raw_status = _call("bash_output", {"job_id": jid}, perms)
    assert reported.get("exit_code") == event.get("exit_code"), (
        f"the notice said exit_code={event.get('exit_code')!r} and "
        f"bash_output said {reported.get('exit_code')!r}: {raw_status[:300]}")
