"""Nine readings of the process table, for one answer.

Building the dashboard asked the system for the running processes once per tab
that lists jobs. Every one of them wanted the same thing about the same
instant, and each cost a whole ``ps`` -- which reads all of /proc whatever
filter it is given -- plus a readlink per process.

Measured on the 384-core host this was found on, 6465 processes: one reading
0.67 s, the eight that followed it 0.0 s instead of 5.4 s. End to end, on the
running dashboard, the server response went 9153 -> 6344 ms with the backend
snapshot and 6344 -> 4665 ms with the per-user one; time to the first widget
12.27 -> 7.62 s.

A snapshot is already a moment in the past by the time it has been parsed, so
reusing it for a moment longer answers the same question. It is dropped the
instant the backend starts a process itself, so a job that has just been
launched is never judged against a table from before it existed.
"""

import subprocess
import time

import pytest

from delfin.dashboard import tab_calculations_browser as browser
from delfin.dashboard.backend_local import LocalJobBackend


@pytest.fixture
def backend(tmp_path):
    made = LocalJobBackend(jobs_file=tmp_path / 'jobs.json')
    yield made
    made._worker_running = False


def test_the_burst_costs_one_reading(backend, monkeypatch):
    readings = []
    monkeypatch.setattr(backend, '_read_processes',
                        lambda: readings.append(1) or [(1, 1, 'x', None)])

    first = backend._list_processes()
    rest = [backend._list_processes() for _ in range(8)]

    assert len(readings) == 1, 'the same instant was asked for nine times'
    assert all(r is first for r in rest)


def test_a_snapshot_does_not_outlive_its_moment(backend, monkeypatch):
    readings = []
    monkeypatch.setattr(backend, '_read_processes',
                        lambda: readings.append(1) or [])
    monkeypatch.setattr(type(backend), '_PROCESS_SNAPSHOT_TTL', 0.05)

    backend._list_processes()
    time.sleep(0.08)
    backend._list_processes()

    assert len(readings) == 2


def test_starting_a_job_throws_the_snapshot_away(backend, monkeypatch):
    """The next thing that happens after a launch is a status check, and a
    table from before the launch would not have the new process in it."""
    monkeypatch.setattr(backend, '_read_processes', lambda: [])
    backend._list_processes()
    assert backend._process_snapshot is not None

    source = open(backend.__class__.__module__.replace('.', '/') + '.py',
                  encoding='utf-8').read()
    started = source.split('proc = subprocess.Popen(')[1].split("job['status'] = 'RUNNING'")[0]
    assert 'self._process_snapshot = None' in started


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
