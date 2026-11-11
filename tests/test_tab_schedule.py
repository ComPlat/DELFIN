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


def test_the_viewer_package_is_not_imported_by_the_dashboard_modules(tmp_path):
    """The startup path must not pull py3Dmol in: importing the viewer and
    the structure editor leaves it unloaded until a viewer is built, and the
    first use installs the script-safe argument rule."""
    import subprocess
    import sys

    from conftest import child_env
    code = (
        "import sys\n"
        "import delfin.dashboard.molecule_viewer as mv\n"
        "import delfin.dashboard.structure_editor\n"
        "assert 'py3Dmol' not in sys.modules, 'imported at startup'\n"
        "m = mv.py3dmol()\n"
        "assert m.tostr._delfin_script_safe\n"
    )
    done = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                          env=child_env(tmp_path))
    assert done.returncode == 0, done.stderr[-800:]


def test_every_visible_tab_is_on_the_bar_and_builds_on_its_first_open(tmp_path):
    """A deferred tab is still on the bar (a placeholder until opened), and
    opening it builds it once: the cold start must not drop the Agent tab,
    the ORCA Builder, ChemDarwin or Job Status from the bar."""
    import subprocess
    import sys

    from conftest import child_env
    code = (
        "import ipywidgets as w\n"
        "import delfin.dashboard as d\n"
        "d.display = lambda *a, **k: None\n"
        "tabs = []\n"
        "init = w.Tab.__init__\n"
        "def record(self, *a, **k):\n"
        "    init(self, *a, **k)\n"
        "    tabs.append(self)\n"
        "w.Tab.__init__ = record\n"
        "d.create_dashboard(backend='local')\n"
        "bar = next(t for t in tabs if 'Submit Job' in t.titles)\n"
        "for title in ('Recalc', 'ORCA Builder', 'ChemDarwin', 'DELFIN Agent', 'Job Status', 'Calculations'):\n"
        "    assert title in bar.titles, (title, bar.titles)\n"
        "i = bar.titles.index('DELFIN Agent')\n"
        "placeholder = bar.children[i]\n"
        "assert not placeholder.children, 'built at start'\n"
        "bar.selected_index = i\n"
        "i = bar.titles.index('DELFIN Agent')\n"
        "assert bar.children[i] is not placeholder, 'not built on open'\n"
        "assert bar.selected_index == i\n"
        "import gc\n"
        "from delfin.dashboard.context import DashboardContext\n"
        "ctx = next(o for o in gc.get_objects() if isinstance(o, DashboardContext))\n"
        "for refs in ('submit_refs', 'orca_builder_refs', 'recalc_refs',\n"
        "             'job_status_refs', 'calc_browser_refs'):\n"
        "    assert getattr(ctx, refs), ('the agent tab drives it', refs)\n"
    )
    done = subprocess.run([sys.executable, "-c", code], capture_output=True,
                          text=True, env=child_env(tmp_path), timeout=600)
    assert done.returncode == 0, done.stderr[-1500:]


# -- warm-up: tabs become ready in idle moments (11.11.) ----------------------

from delfin.dashboard import tab_schedule  # noqa: E402

def _fake_loop():
    calls = []
    return calls, (lambda delay, fn: calls.append((delay, fn)) or True)


def test_warm_up_builds_one_tab_per_idle_step_in_order():
    built = []
    specs = [{"id": "agent", "build": lambda: None},
             {"id": "calc", "build": lambda: None},
             {"id": "done", "widget": object(), "build": lambda: None}]

    def build_now(spec):
        built.append(spec["id"])
        spec["widget"] = object()
        spec.pop("build", None)
    calls, later = _fake_loop()
    assert tab_schedule.warm_up(None, build_now, specs, call_later=later, env={})
    while calls:
        _, fn = calls.pop(0)
        fn()
    assert built == ["agent", "calc"]          # the built one is skipped


def test_a_tab_a_click_built_first_is_not_built_again():
    built = []
    spec = {"id": "agent", "build": lambda: None}
    calls, later = _fake_loop()
    tab_schedule.warm_up(None, lambda s: built.append(s["id"]), [spec, spec],
                         call_later=later, env={})
    spec["widget"] = object()                  # the reader clicked it meanwhile
    while calls:
        calls.pop(0)[1]()
    assert built == []


def test_warm_up_can_be_switched_off():
    calls, later = _fake_loop()
    assert not tab_schedule.warm_up(None, lambda s: None, [{"build": print}],
                                    call_later=later,
                                    env={"DELFIN_DASHBOARD_WARM": "0"})
    assert calls == []


def test_without_a_running_loop_nothing_is_scheduled():
    assert tab_schedule._call_later(0.1, lambda: None) is False
