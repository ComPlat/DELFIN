"""Two processes wrote the same file and one of the writes was the loss.

The defect, measured before this existed: three stores are a full-file
load → mutate → write cycle guarded by a ``threading.Lock``, and each of
their docstrings claimed cross-process safety on the grounds that the
WRITE is atomic. ``os.replace`` stops a reader observing a torn file. It
has nothing to say about two processes that both read, both mutate their
own copy, and both write back — the second write IS the first one's
update, gone.

  * the attention inbox: four processes emitting forty events each left
    48 of 160 records, and the ones that vanish are exactly the ones
    telling the user the agent is blocked on them;
  * the per-user job index: a lost update leaves a live job with no route
    back to its registry, so ``bash_status`` from another process answers
    "unknown job_id" about a job that is running fine;
  * the workspace job registry: the same race between the watchdog
    writing the exit code and a drain setting the acknowledged flag drops
    whichever loses.

``fcntl.flock`` was already used elsewhere in this codebase. What is
pinned here: N processes times M writes leaves exactly N*M records, in
all three, and the ordinary single-process behaviour is unchanged.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from delfin.agent import attention, bash_jobs


PROCS, WRITES = 4, 40

_CHILD = """
import sys
sys.path.insert(0, {root!r})
import pathlib
home = pathlib.Path(sys.argv[1])
pathlib.Path.home = staticmethod(lambda: home)
tag = sys.argv[2]
n = int(sys.argv[3])
{body}
"""

_EMIT = """
from delfin.agent import attention
attention._deliver = lambda event: None
for i in range(n):
    attention.emit_attention("run_finished", title=f"{tag}-{i}")
"""

_INDEX = """
from delfin.agent import bash_jobs
bash_jobs._INDEX_PATH = home / ".delfin" / "bash_jobs_index.json"
bash_jobs._INDEX_PATH.parent.mkdir(parents=True, exist_ok=True)
for i in range(n):
    bash_jobs._note_job_workspace(f"{tag}{i:03d}", "/tmp/ws")
"""

_REGISTRY = """
import time
from delfin.agent import bash_jobs
for i in range(n):
    bash_jobs._persist_job_start(str(home), {
        "job_id": f"{tag}{i:03d}", "pid": 0, "started_at": time.time(),
        "command": "true", "acknowledged": False, "finished_at": None,
    })
"""


def _run_concurrently(body: str, home: Path) -> None:
    """PROCS real processes, all writing the same file at once."""
    root = str(Path(__file__).resolve().parent.parent)
    code = _CHILD.format(root=root, body=body)
    children = [
        subprocess.Popen([sys.executable, "-c", code, str(home),
                          f"p{k}", str(WRITES)])
        for k in range(PROCS)
    ]
    for child in children:
        assert child.wait(timeout=120) == 0


# ---------------------------------------------------------------------------
# N processes x M writes = N*M records, in each of the three
# ---------------------------------------------------------------------------

def test_no_attention_event_is_lost_between_processes(tmp_path):
    _run_concurrently(_EMIT, tmp_path)
    lines = [ln for ln in
             (tmp_path / ".delfin" / "attention_inbox.jsonl")
             .read_text(encoding="utf-8").splitlines() if ln.strip()]
    assert len(lines) == PROCS * WRITES
    assert len({json.loads(ln)["id"] for ln in lines}) == PROCS * WRITES


def test_no_job_stays_unaddressable_after_a_concurrent_index_write(tmp_path):
    _run_concurrently(_INDEX, tmp_path)
    data = json.loads(
        (tmp_path / ".delfin" / "bash_jobs_index.json").read_text("utf-8"))
    assert len(data["jobs"]) == PROCS * WRITES


def test_no_job_record_is_lost_from_the_workspace_registry(tmp_path):
    _run_concurrently(_REGISTRY, tmp_path)
    data = json.loads(
        (tmp_path / ".delfin" / "bash_jobs.json").read_text("utf-8"))
    assert len(data["jobs"]) == PROCS * WRITES


# ---------------------------------------------------------------------------
# The lock itself: exclusive, bounded, and never the thing that breaks
# ---------------------------------------------------------------------------

def test_the_lock_is_exclusive_across_processes(monkeypatch, tmp_path):
    monkeypatch.setattr(bash_jobs, "_LOCK_TIMEOUT_S", 0.3)
    guarded = tmp_path / "state.json"
    holder = subprocess.Popen([
        sys.executable, "-c",
        "import sys, time; sys.path.insert(0, %r)\n"
        "from delfin.agent.bash_jobs import cross_process_lock\n"
        "from pathlib import Path\n"
        "with cross_process_lock(Path(sys.argv[1])):\n"
        "    print('held', flush=True); time.sleep(30)\n"
        % str(Path(__file__).resolve().parent.parent),
        str(guarded)], stdout=subprocess.PIPE, text=True)
    try:
        assert holder.stdout.readline().strip() == "held"
        import time
        started = time.monotonic()
        with bash_jobs.cross_process_lock(guarded):
            waited = time.monotonic() - started
        # It waited out the deadline rather than taking a held lock, and
        # then went on rather than hanging: a bookkeeping file must never
        # become a way to stop the agent.
        assert waited >= bash_jobs._LOCK_TIMEOUT_S
    finally:
        holder.kill()
        holder.wait(timeout=10)


def test_taking_the_lock_never_disturbs_the_file_it_guards(tmp_path):
    guarded = tmp_path / "state.json"
    guarded.write_text('{"jobs": {}}', encoding="utf-8")
    with bash_jobs.cross_process_lock(guarded):
        pass
    assert guarded.read_text(encoding="utf-8") == '{"jobs": {}}'


def test_an_unusable_lock_path_does_not_break_the_caller(tmp_path):
    """Best-effort: the cycle runs whether or not the lock could be
    taken, exactly as it did before this existed."""
    ran = []
    with bash_jobs.cross_process_lock(tmp_path / "no" / "such" / "dir" / "\0"):
        ran.append(True)
    assert ran == [True]


# ---------------------------------------------------------------------------
# Single-process behaviour is unchanged
# ---------------------------------------------------------------------------

def test_the_inbox_still_answers_the_way_it_did(monkeypatch, tmp_path):
    monkeypatch.setattr(Path, "home", staticmethod(lambda: tmp_path))
    monkeypatch.setattr(attention, "_deliver", lambda event: None)
    event_id = attention.emit_attention(
        "question_pending", session_id="s1", title="which folder?")
    assert [ev["id"] for ev in attention.list_pending()] == [event_id]
    assert attention.resolve(event_id, answer="the 2026 one") is True
    drained = attention.drain_resolved("s1")
    assert [ev["answer"] for ev in drained] == ["the 2026 one"]
    assert attention.drain_resolved("s1") == []     # exactly once
    assert attention.list_pending() == []


def test_a_job_record_still_round_trips(tmp_path):
    import time

    now = time.time()
    bash_jobs._persist_job_start(str(tmp_path), {
        "job_id": "abc123", "pid": 0, "started_at": now,
        "command": "true", "acknowledged": False, "finished_at": None})
    bash_jobs._update_job_record(tmp_path, "abc123", exit_code=0,
                                 finished_at=now + 1)
    rec = bash_jobs._load_registry_file(tmp_path)["jobs"]["abc123"]
    assert rec["exit_code"] == 0 and rec["finished_at"] == now + 1
    events = bash_jobs.drain_finished_events(tmp_path)
    assert [ev["job_id"] for ev in events] == ["abc123"]
    assert bash_jobs.drain_finished_events(tmp_path) == []   # exactly once
