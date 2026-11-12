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


def test_startup_builds_every_visible_tab(specs):
    """The person (11.11.): every tab is there at once, none waits for its
    first click."""
    from delfin.dashboard import tab_schedule
    spec_list, build_calls = specs
    startup = tab_schedule.startup_specs_to_build(spec_list, env={})
    assert [s['id'] for s in startup] == ['submit', 'recalc', 'agent']
    assert build_calls == []


def test_lazy_startup_builds_only_the_first_visible_tab(specs):
    from delfin.dashboard import tab_schedule
    spec_list, build_calls = specs
    startup = tab_schedule.startup_specs_to_build(
        spec_list, env={'DELFIN_DASHBOARD_LAZY': '1'})
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
    # The lazy path (DELFIN_DASHBOARD_LAZY=1) still builds on first open,
    # and the tab built then carries its own init scripts (self_running_html):
    # sent as a Javascript output after the page was out, they never ran.
    code += (
        "built = bar.children[i]\n"
        "found = []\n"
        "def walk(w):\n"
        "    if 'delfin-tab-js' in str(getattr(w, 'value', '')): found.append(w)\n"
        "    for c in getattr(w, 'children', ()): walk(c)\n"
        "walk(built)\n"
        "assert found, 'the tab built on open carries no init scripts'\n"
    )
    env = dict(child_env(tmp_path), DELFIN_DASHBOARD_LAZY="1")
    done = subprocess.run([sys.executable, "-c", code], capture_output=True,
                          text=True, env=env, timeout=600)
    assert done.returncode == 0, done.stderr[-1500:]


def test_every_visible_tab_is_built_at_start(tmp_path):
    """The default: the agent tab (and every other visible one) is built
    before the page goes out, with its scripts."""
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
        "for title in ('DELFIN Agent', 'Job Status', 'Calculations'):\n"
        "    child = bar.children[bar.titles.index(title)]\n"
        "    assert getattr(child, 'children', None), ('not built at start', title)\n"
    )
    # Building every tab builds the Literature tab, which keeps its folder
    # in the checkout (literature/, git-ignored). The run must leave the
    # checkout as it found it: a folder this test made goes again.
    import shutil
    from pathlib import Path
    lit = Path(__file__).resolve().parent.parent / "literature"
    had_lit = lit.exists()
    try:
        done = subprocess.run([sys.executable, "-c", code], capture_output=True,
                              text=True, env=child_env(tmp_path), timeout=600)
    finally:
        if not had_lit and lit.is_dir():
            shutil.rmtree(lit, ignore_errors=True)
    assert done.returncode == 0, done.stderr[-1500:]


def test_tabs_built_before_shown_wait_for_their_root():
    """A tab built in the background (warm-up, or the agent tab's needs) has
    no DOM until its first open; its bind script gave up after 4 s and the
    splitter and the explorer's clicks stayed dead (11.11., Calculations)."""
    from pathlib import Path
    root = Path(__file__).resolve().parent.parent / "delfin" / "dashboard"
    calc = (root / "tab_calculations_browser.py").read_text(encoding="utf-8")
    lit = (root / "tab_literature.py").read_text(encoding="utf-8")
    assert "__delfinCalcRootWaits" in calc and "initCalcScopeBind(0)" in calc
    assert "__delfinLitRootWaits" in lit and "initLitSplitter(0)" in lit
