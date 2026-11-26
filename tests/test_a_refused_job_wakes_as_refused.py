"""A job the sandbox refused wakes as REFUSED, not as a failed command.

A ``bash_background`` job whose sandbox could not be applied exits 126 with
"the command was not run" on stderr. ``finished_shells`` rendered that as
``exit 126``, so the wake notice read "shell <id> [exit 126] echo hello" --
as if ``echo hello`` had run and failed. An agent told that spends its turn
debugging a command that never ran. The notice owes the watcher's own
finding: the sandbox refused it.
"""

from __future__ import annotations

import tempfile
import types

from delfin.agent import job_wake


def _stub_job(job_id: str, exit_code: int, stderr_text: str, session_id: str = "sess-1"):
    """A duck-typed BashJob over a real kit_bg stderr tempfile.

    ``read_output`` only reads files that look like what the registry
    creates (``is_own_output_file``), so the stderr lives in a real
    ``kit_bg_*.stderr`` tempfile: regular, single-link, owned by this user.
    """
    import os

    fd, path = tempfile.mkstemp(prefix="kit_bg_", suffix=".stderr")
    os.close(fd)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(stderr_text)

    job = types.SimpleNamespace()
    job.job_id = job_id
    job.session_id = session_id
    job.command = "echo hello"
    job.stderr_path = path
    job.stdout_path = ""
    job.poll = lambda: exit_code
    return job


class _RegistryStub:
    def __init__(self, jobs):
        self._jobs = jobs

    def list_jobs(self, include_finished=True):
        return self._jobs


REFUSED_STDERR = (
    "refused: sandbox could not be applied ([Errno 16] seccomp filter: "
    "Device or resource busy); the command was not run."
)


def test_a_refused_job_wakes_as_refused_not_as_exit():
    job = _stub_job("refused1", 126, REFUSED_STDERR)
    seen = set()
    done = job_wake.finished_shells(
        seen, registry=_RegistryStub([job]), session_id="sess-1")
    assert len(done) == 1
    state = done[0]["state"]
    assert "not run" in state.lower() or "refused" in state.lower(), (
        f"the notice must say the command was never run, got: {state!r}")
    assert "126" not in state, (
        f"the bare exit code reads as the command failing: {state!r}")


def test_a_job_that_ran_and_failed_still_says_exit():
    # Control: a command that RAN and failed keeps its exit-code notice.
    job = _stub_job("ranfail1", 1, "ls: cannot access 'x': No such file\n")
    seen = set()
    done = job_wake.finished_shells(
        seen, registry=_RegistryStub([job]), session_id="sess-1")
    assert done[0]["state"] == "exit 1"
    assert done[0]["ok"] is False
