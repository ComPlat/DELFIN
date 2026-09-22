"""The Job Status tab shows the user's fairshare and CPU limit, always.

Whether waiting jobs start sooner or later depends on the fairshare, and
whether they start at all on the account's CPU limit; neither was visible,
and the summary disappeared with an empty queue.
"""

from datetime import datetime
from pathlib import Path

import pytest

from delfin.dashboard import backend_slurm
from delfin.dashboard.backend_base import JobBackend, JobInfo
from delfin.dashboard.backend_slurm import SlurmJobBackend
from delfin.dashboard.tab_job_status import build_slurm_job_table, empty_queue_html

NOW = datetime(2026, 9, 14, 14, 0)

SSHARE = 'ka|ka_user|0.000679|0.001409|0.394548|0.481840\n'
SACCTMGR = 'ka|ka||cpu=3840||\n'


class FakeClock:
    def __init__(self):
        self.now = 1000.0

    def __call__(self):
        return self.now


@pytest.fixture
def backend(tmp_path, monkeypatch):
    """A SLURM backend whose sshare and sacctmgr are counted, not executed."""
    calls = []
    outputs = {'sshare': SSHARE, 'sacctmgr': SACCTMGR}

    def fake_run(cmd, **kwargs):
        calls.append(cmd[0])
        output = outputs[cmd[0]]
        if isinstance(output, Exception):
            raise output

        class _Result:
            returncode = 0
            stdout = output
            stderr = ''

        return _Result()

    clock = FakeClock()
    monkeypatch.setattr(backend_slurm.subprocess, 'run', fake_run)
    monkeypatch.setattr(backend_slurm.time, 'monotonic', clock)
    obj = SlurmJobBackend(Path(tmp_path), tool_binaries={'orca': '/opt/orca'})
    obj._calls, obj._outputs, obj._clock = calls, outputs, clock
    return obj


def test_the_backend_reads_fairshare_and_the_cpu_limit(backend):
    standing = backend.account_standing()

    assert standing['account'] == 'ka'
    assert standing['fairshare'] == pytest.approx(0.394548)
    assert standing['level_fs'] == pytest.approx(0.48184)
    assert standing['effective_usage'] / standing['norm_shares'] == pytest.approx(2.075, abs=0.001)
    assert standing['cpu_limit'] == 3840
    assert standing['max_jobs'] is None and standing['max_submit'] is None


def test_the_default_account_is_the_one_shown(backend):
    backend._outputs['sshare'] = ('other|ka_user|0.1|0.05|0.8|1.2\n'
                                  'ka|ka_user|0.000679|0.001409|0.394548|0.481840\n')
    backend._outputs['sacctmgr'] = ('ka|other||cpu=10||\n'
                                    'ka|ka|cpu_il|cpu=100||\n'
                                    'ka|ka||cpu=3840|50|200\n')
    standing = backend.account_standing()

    assert standing['fairshare'] == pytest.approx(0.394548)
    assert standing['cpu_limit'] == 3840, 'the association for all partitions, not one partition'
    assert (standing['max_jobs'], standing['max_submit']) == (50, 200)


def test_every_refresh_may_ask_without_loading_the_cluster(backend):
    for _ in range(100):
        backend.account_standing()
    assert backend._calls.count('sshare') == 1
    assert backend._calls.count('sacctmgr') == 1

    backend._clock.now += backend._FAIRSHARE_CACHE_SECONDS + 1
    backend.account_standing()
    assert backend._calls.count('sshare') == 2, 'fairshare is recomputed every few minutes'
    assert backend._calls.count('sacctmgr') == 1, 'limits change rarely'

    backend._clock.now += backend._LIMITS_CACHE_SECONDS
    backend.account_standing()
    assert backend._calls.count('sacctmgr') == 2


def test_a_failing_query_keeps_the_last_answer_and_is_not_retried_every_tick(backend):
    assert backend.account_standing()['fairshare'] == pytest.approx(0.394548)

    backend._outputs['sshare'] = OSError('sshare: connection timed out')
    backend._clock.now += backend._FAIRSHARE_CACHE_SECONDS + 1
    for _ in range(10):
        assert backend.account_standing()['fairshare'] == pytest.approx(0.394548)
    assert backend._calls.count('sshare') == 2


def test_a_cluster_that_tells_nothing_gives_no_standing(backend):
    backend._outputs['sshare'] = OSError('sshare: Access/permission denied')
    backend._outputs['sacctmgr'] = ''
    assert backend.account_standing() is None


def test_other_backends_have_no_standing():
    from delfin.dashboard.backend_local import LocalJobBackend

    assert LocalJobBackend.account_standing is JobBackend.account_standing


STANDING = {'account': 'ka', 'fairshare': 0.394548, 'level_fs': 0.48184,
            'norm_shares': 0.000679, 'effective_usage': 0.001409,
            'cpu_limit': 3840, 'max_jobs': None, 'max_submit': None}


def test_the_summary_shows_fairshare_and_the_cpu_limit():
    job = JobInfo(job_id='901', name='Opt Fe', status='R',
                  extra={'cpus': '40', 'reason': 'uc2n872', 'time_used': '1:00', 'time_limit': '1:00:00'})
    html = build_slurm_job_table([job], [], now=NOW, standing=STANDING)

    assert '0.39</div><div class="djs-stat-label">Fairshare · high is good' in html
    assert 'color:#b26a00;">0.39' in html, 'below one\'s share is marked'
    assert 'Lately you used 2.1× your share.' in html
    assert '0.5 means you used exactly your share' in html
    assert '40<span class="djs-stat-suffix"> / 3840</span>' in html
    assert '40 of the 3840 CPUs you may use at once (account ka)' in html


def test_fairshare_is_in_view_with_an_empty_queue():
    html = empty_queue_html('all', standing=STANDING, now=NOW)
    assert 'No jobs in the queue' in html
    assert 'Fairshare · high is good' in html
    assert '/ 3840' in html


def test_without_a_standing_the_summary_is_as_before():
    html = empty_queue_html('all', now=NOW)
    assert 'Fairshare' not in html and 'djs-stat-suffix"' not in html
    assert '>0</div><div class="djs-stat-label">CPUs in use' in html


def test_a_high_fairshare_reads_as_good_and_a_low_one_as_bad():
    good = empty_queue_html(standing={'fairshare': 0.8}, now=NOW)
    bad = empty_queue_html(standing={'fairshare': 0.05}, now=NOW)
    assert 'color:#2e7d32;">0.80' in good
    assert 'color:#d32f2f;">0.05' in bad
