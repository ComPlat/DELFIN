"""Background work: bounded, and watched without being asked to poll.

Three gaps, all from the background-and-scheduling review.

NO CONCURRENCY CAP. Background SUB-agents have had an eight-slot
semaphore since they were built; background COMMANDS had none at all. On
a laptop an unbounded fan-out is untidy. On a shared cluster node it is
other people's CPU, and the agent is the one participant that can start
hundreds without noticing.

SUBMITTED JOBS WERE WATCHED BY NOBODY. The auto-watch helper is wired
only into the dashboard's Submit tab, never into the agent's own
``sbatch``. ``watch_job`` appears once in the whole prompt pack and has
zero recorded uses — asking the model to remember a second call after
every submit is exactly the kind of rule this project has learned does
not bind. The job id is in the output; the framework can read it.

THE WAKE-UP CEILING WAS AN HOUR. A six-hour calculation could not be
waited on in one call. And the entry recorded no workspace, so the
scheduler stored the process's cwd — under Voila the notebook's
directory, not the agent's — and the daemon then disabled the entry
because the path did not match. A wake-up that silently never fires is
worse than one that was refused.
"""

from __future__ import annotations

import json
import pathlib
import tempfile

import pytest

from delfin.agent import api_client as A
from delfin.agent import bash_jobs


# ---------------------------------------------------------------------------
# The cap
# ---------------------------------------------------------------------------

def test_there_is_a_default_cap():
    assert bash_jobs._DEFAULT_MAX_BG_JOBS >= 1
    assert bash_jobs._max_background_jobs() >= 1


def test_the_cap_is_configurable(monkeypatch):
    monkeypatch.setattr(
        "delfin.user_settings.load_settings",
        lambda *a, **kw: {"agent": {"max_background_jobs": 9}})
    assert bash_jobs._max_background_jobs() == 9


def test_a_nonsense_setting_does_not_disable_the_cap(monkeypatch):
    """A cap that can be switched off by a typo is not a cap."""
    for junk in ("", "many", None, 0, -3):
        monkeypatch.setattr(
            "delfin.user_settings.load_settings",
            lambda *a, **kw: {"agent": {"max_background_jobs": junk}})
        assert bash_jobs._max_background_jobs() >= 1


def test_unreadable_settings_fall_back_to_the_default(monkeypatch):
    monkeypatch.setattr(
        "delfin.user_settings.load_settings",
        lambda *a, **kw: (_ for _ in ()).throw(OSError("no settings")))
    assert bash_jobs._max_background_jobs() == bash_jobs._DEFAULT_MAX_BG_JOBS


def test_the_registry_refuses_rather_than_queues(monkeypatch):
    """A queued job looks started to the model, which then waits for output
    that cannot arrive. Refusing with the count tells it what to do."""
    reg = bash_jobs._Registry.__new__(bash_jobs._Registry)
    monkeypatch.setattr(reg, "count_running", lambda: 4, raising=False)
    monkeypatch.setattr(bash_jobs, "_max_background_jobs", lambda: 4)
    with pytest.raises(ValueError) as excinfo:
        reg.start(command="sleep 1", cwd="/tmp")
    message = str(excinfo.value)
    assert "cap is 4" in message
    assert "bash_status" in message, (
        "the refusal must name the action that clears it")


def test_counting_uses_poll_not_a_stored_flag():
    """A process that died without anyone asking would otherwise hold a
    slot forever, and a cap that leaks slots stops the agent working for a
    reason nobody can act on."""
    source = pathlib.Path(bash_jobs.__file__).read_text(encoding="utf-8")
    body = source[source.index("def count_running"):]
    body = body[:body.index("def list_jobs")]
    assert "poll()" in body


# ---------------------------------------------------------------------------
# Submitted SLURM jobs
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("stdout,expected", [
    ("Submitted batch job 123456\n", ["123456"]),
    ("prep done\nSubmitted batch job 7\nSubmitted batch job 8\n", ["7", "8"]),
    ("nothing was submitted here", []),
    ("", []),
])
def test_every_submitted_job_id_is_found(stdout, expected):
    assert A._SBATCH_SUBMITTED_RE.findall(stdout) == expected


def test_submitted_jobs_are_registered(monkeypatch):
    registered: list[tuple] = []
    monkeypatch.setattr(
        "delfin.agent.job_monitor.register_agent_job",
        lambda ws, jid, description="": registered.append((str(ws), jid)))
    ws = pathlib.Path(tempfile.mkdtemp(prefix="ws_"))
    perms = A.KitToolPermissions(workspace=ws)
    ids = A._auto_watch_submitted_jobs("Submitted batch job 4242\n", perms)
    assert ids == ["4242"]
    assert registered == [(str(ws), "4242")]


def test_registration_failure_never_breaks_the_command(monkeypatch):
    """Watching is a convenience; the command already ran."""
    monkeypatch.setattr(
        "delfin.agent.job_monitor.register_agent_job",
        lambda *a, **kw: (_ for _ in ()).throw(OSError("registry down")))
    perms = A.KitToolPermissions(
        workspace=pathlib.Path(tempfile.mkdtemp(prefix="ws_")))
    assert A._auto_watch_submitted_jobs("Submitted batch job 1\n", perms) == []


def test_no_permissions_means_no_registration():
    assert A._auto_watch_submitted_jobs("Submitted batch job 1\n", None) == []


def test_the_result_tells_the_model_not_to_poll():
    """Otherwise it loops on squeue, which is what the watch exists to
    replace."""
    source = pathlib.Path(A.__file__).read_text(encoding="utf-8")
    assert "Do not\n                \"loop on squeue" in source or \
           "loop on squeue" in source


# ---------------------------------------------------------------------------
# The wake-up
# ---------------------------------------------------------------------------

def test_a_wake_up_can_span_a_long_calculation():
    source = pathlib.Path(A.__file__).read_text(encoding="utf-8")
    # The comment explaining the change sits between the key and the value,
    # so a 400-char window landed inside the prose rather than on the number.
    block = source[source.index('"delay_seconds": {'):]
    assert '"maximum": 86400' in block[:900], (
        "the ceiling is back below a day, so a six-hour job cannot be "
        "waited on in one call")


def test_the_wake_up_records_the_agents_workspace():
    ws = pathlib.Path(tempfile.mkdtemp(prefix="ws_"))
    perms = A.KitToolPermissions(workspace=ws)
    executor = A._DocToolExecutor.__new__(A._DocToolExecutor)
    out = json.loads(executor._execute_scheduler(
        "schedule_wakeup",
        {"delay_seconds": 7200, "prompt": "check the job", "reason": "long run"},
        perms))
    assert out.get("status") == "ok"

    from delfin.agent import scheduler as sched
    entries = [e for e in sched.get_scheduler().list_entries()
               if e.id == out.get("id")]
    try:
        assert entries, "the entry was not stored"
        assert str(entries[0].workspace) == str(ws), (
            "the scheduler recorded the process cwd again; the daemon "
            "disables entries whose path does not match, so the wake-up "
            "would silently never fire")
    finally:
        sched.get_scheduler().delete(out.get("id"))
