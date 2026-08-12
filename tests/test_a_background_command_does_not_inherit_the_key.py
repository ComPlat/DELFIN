"""bash_background handed the model's own provider key to its child.

Foreground bash was fixed to run on a scrubbed environment. ``bash_jobs
.start`` still took ``os.environ.copy()``, and ``bash_background`` is
model-callable -- so the model CHOOSING to background a command was the
thing that handed the credential over. ``env`` and ``printenv`` are on the
auto-allow list; reaching for them is an ordinary debugging reflex.

The output files are 0600, so no third party reads them. What happens is
that the key reaches ``/tmp`` in clear the moment a command prints its
environment, and that text comes back as a tool result.

Nothing ever deleted those files either: 6740 orphaned ``kit_bg_*`` files
were observed in ``/tmp``. Dropping the registry record is the last moment
anything knows where they are, so that is where they are unlinked.
"""

from __future__ import annotations

import os

import pytest

from delfin.agent import bash_jobs


@pytest.fixture
def started_jobs():
    """Start jobs and take their tempfiles away again afterwards.

    A job's output files outlive it by design -- ``bash_output`` is read
    after it finishes -- and their only cleanup is the registry record's
    seven-day prune. A test's registry lives in ``tmp_path`` and is gone
    long before that, so every job a test starts is an immediate orphan.
    That is a large part of how 6740 of them accumulated.
    """
    jobs = []
    yield jobs
    for job in jobs:
        for p in (job.stdout_path, job.stderr_path):
            try:
                p.unlink()
            except OSError:
                pass


# ---------------------------------------------------------------------------
# The environment the child starts from
# ---------------------------------------------------------------------------

def test_the_child_environment_has_no_provider_key(monkeypatch):
    monkeypatch.setenv("ANTHROPIC_API_KEY", "value-must-not-appear")
    monkeypatch.setenv("KIT_TOOLBOX_API_KEY", "value-must-not-appear")
    monkeypatch.setenv("MY_SERVICE_TOKEN", "value-must-not-appear")
    monkeypatch.setenv("SOME_PASSWORD", "value-must-not-appear")

    env = bash_jobs._base_child_env()

    for name in ("ANTHROPIC_API_KEY", "KIT_TOOLBOX_API_KEY",
                 "MY_SERVICE_TOKEN", "SOME_PASSWORD"):
        assert name not in env, name
    assert "value-must-not-appear" not in "".join(env.values())


def test_the_child_environment_is_the_same_one_foreground_bash_uses():
    """Two scrubbers would drift; this pins that there is one."""
    from delfin.agent import api_client

    assert bash_jobs._base_child_env() == api_client._scrubbed_bash_env()


def test_what_a_command_actually_needs_survives(monkeypatch):
    """Scrubbing that breaks ordinary work would be reverted within a
    week -- a background xtb run needs its PATH and its scratch dir."""
    monkeypatch.setenv("XTBPATH", "/opt/xtb")
    env = bash_jobs._base_child_env()
    for name in ("PATH", "HOME", "XTBPATH"):
        assert name in env, name


def test_a_backgrounded_printenv_cannot_show_the_key(monkeypatch, tmp_path,
                                                     started_jobs):
    """The end-to-end shape of the leak: the model backgrounds a command
    that prints its environment, and the output comes back as a result."""
    monkeypatch.setenv("ANTHROPIC_API_KEY", "value-must-not-appear")
    monkeypatch.setattr(bash_jobs, "_INDEX_PATH", tmp_path / "index.json")
    mgr = bash_jobs._Registry()
    job = mgr.start("printenv", cwd=str(tmp_path), workspace=str(tmp_path))
    started_jobs.append(job)
    job.proc.wait(timeout=30)
    out = job.stdout_path.read_text(encoding="utf-8", errors="replace")
    assert "value-must-not-appear" not in out
    assert "ANTHROPIC_API_KEY" not in out


# ---------------------------------------------------------------------------
# The tempfiles
# ---------------------------------------------------------------------------

def test_pruning_a_record_unlinks_its_output_files(tmp_path):
    stdout = tmp_path / "kit_bg_x.stdout"
    stderr = tmp_path / "kit_bg_x.stderr"
    for p in (stdout, stderr):
        p.write_text("output nobody will ever read\n", encoding="utf-8")

    jobs = {"j1": {"started_at": 0.0,          # far older than the cutoff
                   "stdout_path": str(stdout),
                   "stderr_path": str(stderr)}}

    assert bash_jobs._prune_old_records(jobs, now=1_000_000_000.0)

    assert jobs == {}
    assert not stdout.exists()
    assert not stderr.exists()


def test_a_live_record_keeps_its_output_files(tmp_path):
    """The prune must not reach a job whose output is still being read."""
    import time

    stdout = tmp_path / "kit_bg_live.stdout"
    stdout.write_text("still running\n", encoding="utf-8")
    jobs = {"j2": {"started_at": time.time(), "stdout_path": str(stdout)}}

    assert not bash_jobs._prune_old_records(jobs, now=time.time())

    assert "j2" in jobs
    assert stdout.exists()


def test_pruning_survives_a_file_that_is_already_gone(tmp_path):
    """The expected case, not the exceptional one: a job's files may have
    been cleaned by the OS or by a previous prune."""
    jobs = {"j3": {"started_at": 0.0,
                   "stdout_path": str(tmp_path / "vanished.stdout"),
                   "stderr_path": ""}}

    assert bash_jobs._prune_old_records(jobs, now=1_000_000_000.0)
    assert jobs == {}


def test_a_job_output_file_is_not_world_readable(tmp_path, monkeypatch,
                                                 started_jobs):
    import stat

    monkeypatch.setattr(bash_jobs, "_INDEX_PATH", tmp_path / "index.json")
    mgr = bash_jobs._Registry()
    job = mgr.start("echo hi", cwd=str(tmp_path), workspace=str(tmp_path))
    started_jobs.append(job)
    job.proc.wait(timeout=30)
    mode = stat.S_IMODE(os.stat(job.stdout_path).st_mode)
    assert mode & 0o077 == 0, oct(mode)
