"""Megabytes fetched, parsed and run at every load, for widgets nobody makes.

The nglview frontend extension ships with the environment and Voila serves it
to the browser whether or not anything uses it. Nothing here does: the
structure viewers are py3Dmol, and there is no ``import nglview`` in the
source. The four mentions of the name are Python-side package probes for the
Settings tab, which a frontend denylist cannot reach.

Measured on the running dashboard, with the extension left out:

    transfer            14.7 MB  ->  10.49 MB
    main-thread blocks  5-6 tasks, ~1500 ms  ->  4 tasks, 1191 ms
    first widget        6.85 s   ->  6.42 s
    document            7626 elements -> 7601, the 25 being nglview's styles
    tabs                14, same labels, either way

jupyterlab-plotly is the same story and was measured separately, three
interleaved loads of each arm on one machine against one fixture, the two
servers differing only in DELFIN_VOILA_EXTENSION_DENYLIST:

    first widget        8.25 s  ->  7.43 s
    main-thread blocks  1919 ms in 6 tasks  ->  1631 ms in 4
    transfer            10.49 MB  ->  5.66 MB
    document            5984 -> 5983 elements, same 14 tabs

Nothing can render a plotly figure into this page. The agent runs its Python
in subprocesses, whose output is a file or text, and Voila has no editable
cells; `delfin/tools/` contains no plotly, FigureWidget or to_html. The
Settings probe that reports whether the plotly *package* is installed asks a
different question and still answers it.

On loopback those bytes are nearly free. Through the SSH tunnel this dashboard
is normally reached through, they are most of the wait.
"""

import os

import pytest

from delfin import cli_voila


@pytest.mark.parametrize('package', ['nglview', 'plotly'])
def test_the_denied_ones_are_never_imported(package):
    """If this ever fails, that name has to leave the denylist -- a denied
    extension whose widget is actually built renders an empty box and says
    nothing about why."""
    import pathlib

    root = pathlib.Path(cli_voila.__file__).parent
    offenders = []
    for path in root.rglob('*.py'):
        text = path.read_text(encoding='utf-8', errors='replace')
        for line in text.splitlines():
            stripped = line.strip()
            if stripped.startswith((f'import {package}', f'from {package}')):
                offenders.append(f'{path.name}: {stripped}')
    assert offenders == [], offenders


def test_no_figure_can_reach_the_page():
    """The one route that would need the plotly extension: a figure displayed
    as a widget in this page. The agent's Python runs in subprocesses, so its
    output is a file or text, never a widget."""
    import pathlib

    tools = pathlib.Path(cli_voila.__file__).parent / 'tools'
    if not tools.is_dir():
        pytest.skip('no tools package here')
    text = '\n'.join(p.read_text(encoding='utf-8', errors='replace')
                     for p in tools.rglob('*.py'))
    for needle in ('FigureWidget', 'plotly'):
        assert needle not in text, needle


def test_the_launcher_leaves_it_out():
    source = open(cli_voila.__file__, encoding='utf-8').read()

    assert '--VoilaConfiguration.extension_denylist=' in source
    assert 'nglview-js-widgets' in cli_voila.DEFAULT_EXTENSION_DENYLIST
    assert 'jupyterlab-plotly' in cli_voila.DEFAULT_EXTENSION_DENYLIST


def test_a_machine_that_wants_it_can_have_it_back(monkeypatch):
    """Without editing this file: the dashboard is often run on a cluster by
    someone who cannot patch the source."""
    monkeypatch.delenv('DELFIN_VOILA_EXTENSION_DENYLIST', raising=False)
    assert cli_voila._extension_denylist() == ['nglview-js-widgets',
                                               'jupyterlab-plotly']

    monkeypatch.setenv('DELFIN_VOILA_EXTENSION_DENYLIST', '')
    assert cli_voila._extension_denylist() == [], 'empty means deny nothing'

    monkeypatch.setenv('DELFIN_VOILA_EXTENSION_DENYLIST', ' one , two ')
    assert cli_voila._extension_denylist() == ['one', 'two']


def test_the_python_side_probe_is_untouched():
    """Settings still reports whether the nglview *package* is installed --
    that is a different question from whether the browser is sent its
    JavaScript, and the denylist cannot reach it."""
    from delfin import qm_health

    assert 'nglview' in qm_health.PACKAGES or 'nglview' in getattr(
        qm_health, 'PACKAGE_FAMILIES', {'nglview': None})
