"""Periodic widgets must not fork one ``squeue`` per tick.

The dashboard has several independent refresh timers (the agent activity
panel ticks every 5 s, the job-event watcher every 30 s). Each one used to
run its own ``squeue``, which puts hundreds of queries per hour on the
cluster controller -- a shared, single-instance service -- for data that
changes on the scale of a minute. ``list_jobs()`` now serves all callers
from one query per interval.
"""

from pathlib import Path

import pytest

from delfin.dashboard import backend_slurm
from delfin.dashboard.backend_slurm import SlurmJobBackend


SQUEUE_OUTPUT = (
    '  JOBID   PARTITION   NAME   USER  ST  TIME  NODES  NODELIST\n'
    '  123456      single  water   foo   R  1:02      1  uc3n123\n'
)


class FakeClock:
    """Monotonic clock we advance by hand, so tests never sleep."""

    def __init__(self):
        self.now = 1000.0

    def __call__(self):
        return self.now

    def advance(self, seconds):
        self.now += seconds


@pytest.fixture
def backend(tmp_path, monkeypatch):
    """A SLURM backend whose squeue calls are counted, not executed."""
    calls = []

    def fake_run(cmd, **kwargs):
        calls.append(list(cmd))

        class _Result:
            returncode = 0
            stdout = SQUEUE_OUTPUT if cmd[0] == 'squeue' else ''
            stderr = ''

        return _Result()

    monkeypatch.setattr(backend_slurm.subprocess, 'run', fake_run)
    monkeypatch.delenv('DELFIN_SQUEUE_MIN_INTERVAL', raising=False)
    # A non-empty tool_binaries mapping skips the auto-detection scan.
    obj = SlurmJobBackend(Path(tmp_path), tool_binaries={'orca': '/opt/orca'})
    obj._calls = calls
    return obj


def _squeue_calls(backend_obj):
    return [c for c in backend_obj._calls if c and c[0] == 'squeue']


# --------------------------------------------------------------------------
# The rate limit itself

def test_first_call_queries_slurm(backend):
    jobs = backend.list_jobs()
    assert len(jobs) == 1
    assert jobs[0].job_id == '123456'
    assert len(_squeue_calls(backend)) == 1


def test_repeated_calls_within_the_interval_share_one_query(backend):
    for _ in range(20):
        assert len(backend.list_jobs()) == 1
    assert len(_squeue_calls(backend)) == 1


def test_query_repeats_once_the_interval_has_passed(backend, monkeypatch):
    clock = FakeClock()
    monkeypatch.setattr(backend_slurm.time, 'monotonic', clock)

    backend.list_jobs()
    clock.advance(backend._JOBS_CACHE_SECONDS - 0.1)
    backend.list_jobs()
    assert len(_squeue_calls(backend)) == 1, 'refreshed too early'

    clock.advance(0.2)
    backend.list_jobs()
    assert len(_squeue_calls(backend)) == 2


def test_five_second_polling_for_an_hour_stays_bounded(backend, monkeypatch):
    """The regression in numbers: 720 ticks/h must not be 720 queries/h."""
    clock = FakeClock()
    monkeypatch.setattr(backend_slurm.time, 'monotonic', clock)

    for _ in range(720):  # one hour at the activity panel's 5 s tick
        backend.list_jobs()
        clock.advance(5.0)

    expected = 3600 / backend._JOBS_CACHE_SECONDS
    assert len(_squeue_calls(backend)) <= expected + 1
    assert len(_squeue_calls(backend)) < 200


def test_force_bypasses_the_cache(backend):
    backend.list_jobs()
    backend.list_jobs(force=True)
    backend.list_jobs(force=True)
    assert len(_squeue_calls(backend)) == 3


def test_interval_is_configurable_per_site(backend, monkeypatch):
    clock = FakeClock()
    monkeypatch.setattr(backend_slurm.time, 'monotonic', clock)
    monkeypatch.setenv('DELFIN_SQUEUE_MIN_INTERVAL', '120')

    backend.list_jobs()
    clock.advance(60.0)
    backend.list_jobs()
    assert len(_squeue_calls(backend)) == 1

    clock.advance(61.0)
    backend.list_jobs()
    assert len(_squeue_calls(backend)) == 2


def test_unparsable_interval_falls_back_to_the_default(backend, monkeypatch):
    monkeypatch.setenv('DELFIN_SQUEUE_MIN_INTERVAL', 'soon')
    assert backend._jobs_cache_seconds() == backend._JOBS_CACHE_SECONDS


def test_a_backwards_clock_does_not_freeze_the_cache(backend, monkeypatch):
    """Guard against a negative age wedging the cache forever."""
    clock = FakeClock()
    monkeypatch.setattr(backend_slurm.time, 'monotonic', clock)
    backend.list_jobs()
    clock.now -= 10_000.0
    backend.list_jobs()
    assert len(_squeue_calls(backend)) == 2


# --------------------------------------------------------------------------
# Freshness where the user expects it

def test_cancel_makes_the_next_listing_fresh(backend):
    backend.list_jobs()
    ok, _msg = backend.cancel_job('123456')
    assert ok
    backend.list_jobs()
    assert len(_squeue_calls(backend)) == 2


def test_callers_cannot_mutate_the_cached_list(backend):
    first = backend.list_jobs()
    first.clear()
    assert len(backend.list_jobs()) == 1


def test_an_idle_queue_is_reported_as_empty(backend, monkeypatch):
    """Header-only output means "no jobs", which is a valid answer to cache."""
    def empty_run(cmd, **kwargs):
        backend._calls.append(list(cmd))

        class _Result:
            returncode = 0
            stdout = '  JOBID PARTITION NAME USER ST TIME NODES NODELIST\n'
            stderr = ''

        return _Result()

    monkeypatch.setattr(backend_slurm.subprocess, 'run', empty_run)
    assert backend.list_jobs() == []
    assert backend.list_jobs() == []
    assert len(_squeue_calls(backend)) == 1


def test_a_failing_squeue_is_rate_limited_too(backend, monkeypatch):
    """A struggling controller is exactly when retrying every tick is worst."""
    clock = FakeClock()
    monkeypatch.setattr(backend_slurm.time, 'monotonic', clock)

    def failing_run(cmd, **kwargs):
        backend._calls.append(list(cmd))
        raise OSError('squeue: connection timed out')

    monkeypatch.setattr(backend_slurm.subprocess, 'run', failing_run)
    for _ in range(10):
        backend.list_jobs()
        clock.advance(5.0)

    assert len(_squeue_calls(backend)) <= 3


def test_a_failing_squeue_keeps_the_last_good_answer(backend, monkeypatch):
    """The job list must not flash empty because one query hiccupped."""
    clock = FakeClock()
    monkeypatch.setattr(backend_slurm.time, 'monotonic', clock)
    assert len(backend.list_jobs()) == 1

    def failing_run(cmd, **kwargs):
        backend._calls.append(list(cmd))
        raise OSError('squeue: connection timed out')

    monkeypatch.setattr(backend_slurm.subprocess, 'run', failing_run)
    clock.advance(backend._JOBS_CACHE_SECONDS + 1.0)
    assert len(backend.list_jobs()) == 1, 'stale-but-true beats a false empty'


# --------------------------------------------------------------------------
# Interface + regression guards

def test_local_backend_accepts_the_same_signature():
    """force= must be safe to pass whatever backend is active."""
    import inspect

    from delfin.dashboard.backend_local import LocalJobBackend

    sig = inspect.signature(LocalJobBackend.list_jobs)
    assert 'force' in sig.parameters
    assert sig.parameters['force'].default is False


@pytest.mark.parametrize(
    'module_name',
    ['tab_agent_activity', 'tab_agent', 'tab_calculations_browser'],
)
def test_periodic_widgets_never_force_a_refresh(module_name):
    """A timer that forces the refresh would defeat the rate limit."""
    import delfin.dashboard as dash

    source = Path(dash.__file__).with_name(f'{module_name}.py')
    text = source.read_text(encoding='utf-8')
    assert 'list_jobs(force=True)' not in text, (
        f'{module_name} polls on a timer and must use the shared cache'
    )
