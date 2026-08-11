"""4.55 MB fetched, parsed and run at every load, for a widget nobody makes.

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

On loopback those bytes are nearly free. Through the SSH tunnel this dashboard
is normally reached through, they are most of the wait.
"""

import os

import pytest

from delfin import cli_voila


def test_nglview_is_never_imported():
    """If this ever fails, the denylist has to go -- a denied extension whose
    widget is actually built renders "model not found" and says nothing."""
    import pathlib

    root = pathlib.Path(cli_voila.__file__).parent
    offenders = []
    for path in root.rglob('*.py'):
        text = path.read_text(encoding='utf-8', errors='replace')
        for line in text.splitlines():
            stripped = line.strip()
            if stripped.startswith(('import nglview', 'from nglview')):
                offenders.append(f'{path.name}: {stripped}')
    assert offenders == [], offenders


def test_the_launcher_leaves_it_out():
    source = open(cli_voila.__file__, encoding='utf-8').read()

    assert '--VoilaConfiguration.extension_denylist=' in source
    assert 'nglview-js-widgets' in cli_voila.DEFAULT_EXTENSION_DENYLIST


def test_a_machine_that_wants_it_can_have_it_back(monkeypatch):
    """Without editing this file: the dashboard is often run on a cluster by
    someone who cannot patch the source."""
    monkeypatch.delenv('DELFIN_VOILA_EXTENSION_DENYLIST', raising=False)
    assert cli_voila._extension_denylist() == ['nglview-js-widgets']

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
