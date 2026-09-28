"""A finished shell wakes the session that started it, and no other.

Input: the job registry and the waking session's id. Output: the finished
shells belonging to that session.

Measured 2026-09-28, in three bug reports from one wave: a benchmark shell
(`bash bench_trial/run_trial4_login.sh`) finished, and three sessions were
woken. Two of them spent a turn establishing that the job was not theirs;
one wrote to the owner to ask; the owner replied "please ignore". One
event, three turns, two of them pure noise -- and the reports say so in
the sessions' own words ("Watcher-Querverbindung", "stale from an earlier
session context").

The dashboard's wake tick already scoped its two other sources by session
-- watched jobs and background agents -- so the asymmetry was in this one
call. The terminal REPL had the same gap.

A job with no session recorded still reaches everyone: an unowned job is
better announced twice than lost, and jobs predating the field have none.
"""

from __future__ import annotations

import pytest

from delfin.agent import job_wake


class _Job:
    def __init__(self, job_id, code, session_id="", command="a command"):
        self.job_id = job_id
        self._code = code
        self.session_id = session_id
        self.command = command

    def poll(self):
        return self._code


class _Registry:
    def __init__(self, *jobs):
        self._jobs = list(jobs)

    def list_jobs(self, include_finished=True):
        return [j for j in self._jobs
                if include_finished or j.poll() is None]


def _ids(out):
    return sorted(d["job_id"] for d in out)


def test_only_the_owning_session_is_told():
    reg = _Registry(_Job("mine", 0, "session-a"),
                    _Job("theirs", 0, "session-b"))
    assert _ids(job_wake.finished_shells(set(), registry=reg,
                                         session_id="session-a")) == ["mine"]


def test_the_other_session_is_told_about_its_own():
    reg = _Registry(_Job("mine", 0, "session-a"),
                    _Job("theirs", 0, "session-b"))
    assert _ids(job_wake.finished_shells(set(), registry=reg,
                                         session_id="session-b")) == ["theirs"]


def test_a_job_with_no_session_still_reaches_everyone():
    """Jobs from before the field existed, and anything started outside a
    session. Announced twice beats lost."""
    reg = _Registry(_Job("unowned", 0, ""), _Job("theirs", 0, "session-b"))
    assert _ids(job_wake.finished_shells(set(), registry=reg,
                                         session_id="session-a")) == ["unowned"]


def test_without_a_session_nothing_is_filtered():
    """The old signature still means the old thing, so a caller that has
    no session to give is not silently given an empty list."""
    reg = _Registry(_Job("a", 0, "session-a"), _Job("b", 0, "session-b"))
    assert _ids(job_wake.finished_shells(set(), registry=reg)) == ["a", "b"]


def test_a_running_job_is_not_reported():
    reg = _Registry(_Job("running", None, "session-a"))
    assert job_wake.finished_shells(set(), registry=reg,
                                    session_id="session-a") == []


def test_a_job_is_reported_once():
    reg = _Registry(_Job("mine", 0, "session-a"))
    seen: set = set()
    assert len(job_wake.finished_shells(seen, registry=reg,
                                        session_id="session-a")) == 1
    assert job_wake.finished_shells(seen, registry=reg,
                                    session_id="session-a") == []


def test_a_registry_that_raises_wakes_nobody():
    class _Broken:
        def list_jobs(self, include_finished=True):
            raise RuntimeError("registry gone")

    assert job_wake.finished_shells(set(), registry=_Broken(),
                                    session_id="s") == []


# -- both surfaces pass it -------------------------------------------------

@pytest.mark.parametrize("module, needle", [
    ("delfin/dashboard/tab_agent.py", "finished_shells(seen, session_id="),
    ("delfin/agent/repl.py", "job_wake.finished_shells("),
])
def test_the_callers_scope_it(module, needle):
    """A filter no caller uses reports everything to everyone, which is
    where this started."""
    import pathlib

    root = pathlib.Path(job_wake.__file__).resolve().parents[2]
    src = (root / module).read_text(encoding="utf-8")
    assert needle in src, module
    at = src.index(needle)
    assert "session_id=" in src[at:at + 220], (
        f"{module} calls finished_shells without a session")
