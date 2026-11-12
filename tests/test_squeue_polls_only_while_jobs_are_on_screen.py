"""squeue is asked on a timer only while a view of jobs is on screen.

Found in the bug hunt of 11.11. (job 8543241): with the Submit Job or the
Agent tab in front, ``squeue`` still ran about every 25 s for the whole life
of the dashboard (the Agent Activity panel ticks every 5 s, the agent's job
watcher every 30 s, and both ask the backend). VS Code and the GitHub UI
pause polling for views that are not on screen.

Now the dashboard tells the backend which tab is on screen
(``DashboardContext.jobs_view_visible``). With Job Status or Agent Activity
in front, periodic callers keep the normal 25 s throttle. Anywhere else they
get the last answer, and squeue is asked at most every 5 minutes (so a job
that ends is still noticed, later), or never with
``DELFIN_SQUEUE_BACKGROUND_INTERVAL=0``. Refresh (``force=True``) always asks.
"""

from __future__ import annotations

from pathlib import Path

import ipywidgets as widgets
import pytest

from delfin.dashboard import backend_slurm
from delfin.dashboard.backend_slurm import SlurmJobBackend
from delfin.dashboard.context import DashboardContext

SQUEUE_OUTPUT = (
    '  JOBID   PARTITION   NAME   USER  ST  TIME  NODES  NODELIST\n'
    '  123456      single  water   foo   R  1:02      1  node123\n'
)
TITLES = ['DELFIN Agent', 'Submit Job', 'Job Status', 'Agent Activity']


class FakeClock:
    def __init__(self):
        self.now = 1000.0

    def __call__(self):
        return self.now


@pytest.fixture
def board(tmp_path, monkeypatch):
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
    monkeypatch.delenv('DELFIN_SQUEUE_BACKGROUND_INTERVAL', raising=False)
    clock = FakeClock()
    monkeypatch.setattr(backend_slurm.time, 'monotonic', clock)
    backend = SlurmJobBackend(Path(tmp_path), tool_binaries={'orca': '/opt/orca'})
    ctx = DashboardContext(backend=backend)
    ctx.tabs_widget = widgets.Tab(children=[widgets.VBox() for _ in TITLES])
    ctx.tab_indices = {t: i for i, t in enumerate(TITLES)}
    ctx.tab_indices_title_of = lambda i: TITLES[i] if i is not None and 0 <= i < len(TITLES) else None

    def show(title):
        ctx.tabs_widget.selected_index = TITLES.index(title)

    def squeues():
        return len([c for c in calls if c and c[0] == 'squeue'])

    def poll_for(seconds, every=5.0):
        t = 0.0
        while t < seconds:
            backend.list_jobs()
            clock.now += every
            t += every

    return ctx, backend, show, squeues, poll_for


def test_the_context_hands_its_visibility_to_the_backend(board):
    ctx, backend, show, _squeues, _poll = board
    assert backend.jobs_view_visible == ctx.jobs_view_visible
    show('Submit Job')
    assert ctx.jobs_view_visible() is False
    show('Job Status')
    assert ctx.jobs_view_visible() is True
    show('Agent Activity')
    assert ctx.jobs_view_visible() is True


def test_no_tab_bar_means_visible_as_before():
    assert DashboardContext().jobs_view_visible() is True


def test_job_status_on_screen_keeps_the_normal_pace(board):
    _ctx, _backend, show, squeues, poll_for = board
    show('Job Status')
    poll_for(600)  # ten minutes of a 5 s tick
    assert 20 <= squeues() <= 25


def test_another_tab_on_screen_polls_at_the_background_pace(board):
    _ctx, backend, show, squeues, poll_for = board
    show('Submit Job')
    backend.list_jobs()  # the first answer (start-up) is always asked for
    before = squeues()
    poll_for(600)
    assert squeues() - before <= 2  # 300 s background pace, not 25 s
    assert backend.list_jobs()[0].job_id == '123456'  # the cache still answers


def test_background_polling_can_be_switched_off(board, monkeypatch):
    _ctx, backend, show, squeues, poll_for = board
    monkeypatch.setenv('DELFIN_SQUEUE_BACKGROUND_INTERVAL', '0')
    show('DELFIN Agent')
    backend.list_jobs()
    before = squeues()
    poll_for(3600)
    assert squeues() == before


def test_refresh_always_asks(board):
    _ctx, backend, show, squeues, _poll = board
    show('Submit Job')
    backend.list_jobs()
    before = squeues()
    backend.list_jobs(force=True)
    assert squeues() == before + 1


def test_opening_a_job_view_resumes_at_once(board):
    ctx, backend, show, squeues, poll_for = board
    show('Submit Job')
    backend.list_jobs()
    poll_for(60)
    before = squeues()
    show('Agent Activity')
    backend.list_jobs()  # cache older than 25 s: asks now
    assert squeues() == before + 1


def test_a_backend_without_a_dashboard_is_unchanged(tmp_path, monkeypatch):
    monkeypatch.setattr(backend_slurm.subprocess, 'run', lambda cmd, **kw: type(
        'R', (), {'returncode': 0, 'stdout': SQUEUE_OUTPUT, 'stderr': ''})())
    backend = SlurmJobBackend(Path(tmp_path), tool_binaries={'orca': '/opt/orca'})
    assert backend.jobs_view_visible is None
    assert backend._jobs_view_open() is True


def test_the_activity_panel_refreshes_only_while_it_is_shown():
    source = (Path(__file__).resolve().parents[1] / 'delfin' / 'dashboard'
              / 'tab_agent_activity.py').read_text(encoding='utf-8')
    loop = source[source.index('def _live_loop'):source.index('def _reload')]
    assert 'selected_tab_title' in loop
    assert '"Agent Activity"' in loop
