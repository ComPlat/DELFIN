"""Lazy tab building: at startup only the first visible tab is built; every
other tab builds on its first open and then keeps its state.

The dashboard must not lie: a tab that is on the bar but not yet built is
still reachable, opens on first click, and a second open shows the same
widget object that was built the first time.
"""
import pytest


def _spec(tab_id, title, order, built=None, build_calls=None):
    def build():
        build_calls.append(tab_id)
        return built[tab_id]
    return {
        'id': tab_id,
        'title': title,
        'widget': None,
        'build': build,
        'default_order': order,
        'default_visible': True,
        'available': True,
        'fixed': False,
        'reason': '',
    }


@pytest.fixture()
def widgets():
    return {key: {'mark': key} for key in ('submit', 'recalc', 'agent')}


@pytest.fixture()
def specs(widgets):
    build_calls = []
    return [
        _spec('submit', 'Submit Job', 10, widgets, build_calls),
        _spec('recalc', 'Recalc', 20, widgets, build_calls),
        _spec('agent', 'Agent', 30, widgets, build_calls),
    ], build_calls


def test_startup_builds_only_the_first_visible_tab(specs):
    from delfin.dashboard import tab_schedule
    spec_list, build_calls = specs
    startup = tab_schedule.startup_specs_to_build(spec_list)
    assert [s['id'] for s in startup] == ['submit']
    assert build_calls == []


def test_open_spec_builds_once_and_keeps_state(specs):
    from delfin.dashboard import tab_schedule
    spec_list, build_calls = specs
    _, _submit = spec_list[0], spec_list[1]
    first = tab_schedule.open_spec(spec_list[1])
    second = tab_schedule.open_spec(spec_list[1])
    assert first is second
    assert first == {'mark': 'recalc'}
    assert build_calls == ['recalc']


def test_open_spec_survives_a_failing_build(specs):
    from delfin.dashboard import tab_schedule
    spec_list, _ = specs
    broken = dict(spec_list[2])
    broken['build'] = lambda: (_ for _ in ()).throw(RuntimeError('no backend'))
    assert tab_schedule.open_spec(broken) is None


def test_open_spec_is_idempotent_after_an_explicit_widget(specs):
    from delfin.dashboard import tab_schedule
    spec_list, build_calls = specs
    widget = {'mark': 'prebuilt'}
    spec_list[0]['widget'] = widget
    assert tab_schedule.open_spec(spec_list[0]) is widget
    assert build_calls == []


def test_the_viewer_package_is_not_imported_by_the_dashboard_modules():
    """The startup path must not pull py3Dmol in: importing the viewer and
    the structure editor leaves it unloaded until a viewer is built, and the
    first use installs the script-safe argument rule."""
    import subprocess
    import sys
    code = (
        "import sys\n"
        "import delfin.dashboard.molecule_viewer as mv\n"
        "import delfin.dashboard.structure_editor\n"
        "assert 'py3Dmol' not in sys.modules, 'imported at startup'\n"
        "m = mv.py3dmol()\n"
        "assert m.tostr._delfin_script_safe\n"
    )
    done = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert done.returncode == 0, done.stderr[-800:]
