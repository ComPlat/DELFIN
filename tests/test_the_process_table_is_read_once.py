"""Nine readings of the process table, for a question nobody had asked.

Building the dashboard read the system's process list once per tab that lists
jobs. Each reading costs a whole ``ps`` -- which reads all of /proc whatever
filter it is given -- plus a readlink per process: 0.67 s on the 384-core host
this was found on, with 6465 processes. Nine of them is 5.4 s of a 9.2 s cold
start.

The table has exactly one consumer: ``_update_job_status``, and that runs only
for a job whose status is RUNNING. With nothing running, every one of those
readings was taken and thrown away. So the callers ask first.

    9x list_jobs, nothing running   9 readings  ->  0
    server response                 9153 ms     ->  4665 ms
    first widget                    12.27 s     ->  7.62 s

A first attempt cached the reading for two seconds instead. That was worse in
a way that matters: the table is what decides whether a job that no longer
answers to its own pid is still alive through its children, and judging that
against a table from before the job existed marks a live job FAILED -- which
is written to ~/.delfin_jobs.json and does not come back. Asking only when
there is something to ask about has no staleness in it at all.
"""

import json
import os
import subprocess
import sys
import time

import pytest

from delfin.dashboard import tab_calculations_browser as browser
from delfin.dashboard.backend_local import LocalJobBackend


@pytest.fixture
def jobs_file(tmp_path):
    return tmp_path / 'jobs.json'


def _backend(jobs_file, jobs):
    jobs_file.write_text(json.dumps({'_next_job_id': 1, 'jobs': jobs}))
    made = LocalJobBackend(jobs_file=jobs_file)
    made._worker_running = False
    return made


def test_nothing_running_means_nothing_is_asked(jobs_file, tmp_path):
    made = _backend(jobs_file, [{'job_id': 1, 'status': 'COMPLETED',
                                 'job_dir': str(tmp_path)}])
    readings = []
    made._list_processes = lambda: readings.append(1) or []

    for _ in range(9):
        made.list_jobs()

    assert readings == [], 'the table has no consumer when nothing is RUNNING'


def test_a_running_job_is_still_judged_against_a_fresh_reading(jobs_file, tmp_path):
    """Not a cached one. The table decides whether a job that no longer answers
    to its own pid is alive through its children; a stale one marks it FAILED,
    and that is written to disk."""
    made = _backend(jobs_file, [{'job_id': 2, 'status': 'RUNNING',
                                 'job_dir': str(tmp_path)}])
    readings = []
    # A live child of the job, so it stays RUNNING and is asked about again.
    # Any pid but this process's own: the reading skips itself, and a fixed
    # "surely nobody has this" number is not safe -- 999999 is a perfectly
    # ordinary pid where pid_max is four million, which is where this failed.
    other_pid = os.getpid() + 1
    made._list_processes = lambda: (
        readings.append(1) or [(other_pid, other_pid, 'orca', str(tmp_path))])

    made.list_jobs()
    after_one = len(readings)
    made.list_jobs()

    assert json.loads(jobs_file.read_text())['jobs'][0]['status'] == 'RUNNING'
    assert after_one >= 1, 'a running job has to be checked against something'
    assert len(readings) > after_one, 'each round asks again, nothing is kept'


def test_a_live_job_survives_the_round(jobs_file, tmp_path):
    child = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(20)'],
                             cwd=str(tmp_path), start_new_session=True)
    try:
        time.sleep(0.4)
        made = _backend(jobs_file, [{
            'job_id': 3, 'status': 'RUNNING', 'job_dir': str(tmp_path),
            'pid': child.pid, 'pgid': os.getpgid(child.pid)}])

        made.list_jobs()

        assert json.loads(jobs_file.read_text())['jobs'][0]['status'] == 'RUNNING'
    finally:
        child.terminate()
        child.wait()


def test_a_job_whose_processes_are_gone_is_marked(jobs_file, tmp_path):
    child = subprocess.Popen([sys.executable, '-c', 'pass'], cwd=str(tmp_path),
                             start_new_session=True)
    child.wait()
    made = _backend(jobs_file, [{
        'job_id': 4, 'status': 'RUNNING', 'job_dir': str(tmp_path / 'gone'),
        'pid': child.pid, 'pgid': child.pid}])

    made.list_jobs()

    assert json.loads(jobs_file.read_text())['jobs'][0]['status'] == 'FAILED'


def test_an_empty_table_is_an_answer_and_not_a_missing_one(jobs_file, tmp_path):
    """``[]`` means 'nothing of this job is running'. Read as 'no table given'
    -- which is what a falsy test does -- it forks another ``ps`` per job."""
    made = _backend(jobs_file, [])
    readings = []
    made._list_processes = lambda: readings.append(1) or []

    made._job_has_active_processes({'job_dir': str(tmp_path)}, process_table=[])

    assert readings == []


# ---------------------------------------------------------------------------
# and the same question asked three times by three builds of one tab
# ---------------------------------------------------------------------------


def test_the_three_builds_of_the_browser_share_one_reading(monkeypatch):
    """Calculations, Archive and Office are three builds of one builder, and
    each lists its directory as it is created."""
    calls = []

    def fake_run(cmd, *a, **k):
        calls.append(cmd)
        return subprocess.CompletedProcess(cmd, 0, stdout='1 orca\n', stderr='')

    monkeypatch.setattr(browser.subprocess, 'run', fake_run)
    browser._user_process_snapshot[:] = [0.0, None]

    out = [browser._recent_user_processes('someone') for _ in range(3)]

    assert len(calls) == 1, calls
    assert out == ['1 orca\n'] * 3
    assert calls[0][:2] == ['ps', '-u']


def test_a_process_listing_that_fails_is_not_remembered(monkeypatch):
    """A failed reading must not be served for the next two seconds as though
    it were an answer -- the caller treats None as 'could not tell'."""
    monkeypatch.setattr(browser.subprocess, 'run',
                        lambda *a, **k: subprocess.CompletedProcess(a[0], 1, '', ''))
    browser._user_process_snapshot[:] = [0.0, None]

    assert browser._recent_user_processes('someone') is None
    assert browser._user_process_snapshot[1] is None
