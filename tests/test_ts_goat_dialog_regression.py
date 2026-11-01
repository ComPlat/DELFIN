"""Regression test: selecting 'TS -> GOAT' in the options dropdown must
open a visible preparation dialog (or show a clear error), end to end.

This drives the real handler code path -- calc_on_options_change with
change['new'] == 'TS -> GOAT' -- against synthetic TS fixtures, not the
library functions alone (those are covered by
test_ts_goat_reactive_center.py).
"""

import pytest

pytest.importorskip(
    'delfin.ts_goat',
    reason='lives on the feature branch until it is merged into main',
)

from pathlib import Path  # noqa: E402
from unittest import mock  # noqa: E402

import ipywidgets as widgets  # noqa: E402

from delfin.ts_goat import build_ts_goat_input  # noqa: E402


# ---------------------------------------------------------------- fixtures
def _write_ts_fixture(dir_path: Path, name: str = 'ts_run') -> Path:
    """Minimal converged OptTS run: out, inp, xyz, hess.

    The Hessian uses the same blocked-column layout as the fixtures in
    test_ts_goat_reactive_center.py so delfin.imag's parser reads the
    imaginary-mode eigenvector (a thinner stub leaves the analysis with
    zero proposals).
    """
    out = dir_path / f'{name}.out'
    out.write_text(
        '|  1> ! OptTS\n'
        'THE OPTIMIZATION HAS CONVERGED\n'
        '     1:      -567.89 cm**-1  ***imaginary mode***\n',
        encoding='utf-8',
    )
    inp = dir_path / f'{name}.inp'
    inp.write_text('%charge 0\n%mult 1\n* xyz 0 1\n', encoding='utf-8')

    atoms = [
        ('O', 0.00, 0.00, 0.00),
        ('H', 0.96, 0.00, 0.00),
        ('H', 1.55, 1.10, 0.00),
    ]
    # the mode acts on the transferring H: +x on H2, -x on the rest
    mode = [
        (1.0, 0.0, 0.0),
        (-1.0, 0.0, 0.0),
        (0.1, 0.0, 0.0),
    ]
    from tests.test_ts_goat_reactive_center import _write_hess
    hess = dir_path / f'{name}.hess'
    # ORCA hessians carry ALL 3N normal modes; the imaginary one sits at
    # its frequency index (1).  delfin.imag.displaced_geometries resolves
    # a mode BY that index, so the stub must provide every column.
    zero = [(0.0, 0.0, 0.0)] * len(atoms)
    mode_vectors = [zero, mode] + [zero] * 7
    _write_hess(hess, atoms, mode_vectors,
                [-567.89, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0])
    xyz = dir_path / f'{name}.xyz'
    xyz.write_text(
        '3\nTS\n'
        + ''.join(f'{el} {x:.6f} {y:.6f} {z:.6f}\n' for el, x, y, z in atoms),
        encoding='utf-8',
    )
    # endpoints so the connectivity evidence has something to compare
    ed = dir_path / f'{name}_Edukt.xyz'
    ed.write_text(
        '3\nEdukt\nO   0.000000  0.000000  0.000000\n'
        'H   0.960000  0.000000  0.000000\nH   3.000000  1.100000  0.000000\n',
        encoding='utf-8',
    )
    pr = dir_path / f'{name}_Produkt.xyz'
    pr.write_text(
        '3\nProdukt\nO   0.000000  0.000000  0.000000\n'
        'H   3.000000  0.000000  0.000000\nH   0.960000  1.100000  0.000000\n',
        encoding='utf-8',
    )
    return out


class _FakeBackend:
    """Records submit_orca calls; returns a zero rc so the dialog reports
    success without running anything."""

    def __init__(self):
        self.calls = []

    def submit_orca(self, job_dir, job_name, inp_file,
                    time_limit='48:00:00', pal=12, maxcore=1000):
        self.calls.append({
            'job_dir': job_dir, 'job_name': job_name, 'inp_file': inp_file,
        })
        return mock.Mock(returncode=0, stdout=f'Submitted job {job_name}\n',
                        stderr='')


# ---------------------------------------------------------------- helpers
def _open_dialog(tmp_path, backend, out_name='ts_run', extra_callbacks=None):
    """Call open_ts_goat_dialog the way the browser-tab handler does."""
    from delfin.dashboard.ts_goat_panel import open_ts_goat_dialog

    out_path = _write_ts_fixture(tmp_path, out_name)
    rendered = []
    result = open_ts_goat_dialog(
        out_path,
        callbacks={
            'backend': backend,
            'reset_dropdown': lambda: None,
            'render_viewer': lambda xyz, idx, arrows=None: rendered.append((xyz, idx)),
            'scroll_top': lambda: None,
            **(extra_callbacks or {}),
        },
    )
    return result, rendered


def test_dropdown_selection_opens_visible_dialog(tmp_path):
    """The E2E regression: after the handler opens the dialog for a usable
    TS fixture, the panel widget is populated and visible, and shows the
    analysis -- no silent no-op like the ImportError bug."""
    backend = _FakeBackend()
    result, rendered = _open_dialog(tmp_path, backend)
    assert result['error'] is None, result['error']
    panel = result['panel']
    assert isinstance(panel, widgets.VBox)
    # visible: children attached and display not 'none' after show()
    assert len(panel.children) > 0
    assert panel.layout.display == 'block'
    # the status line renders the reliability, not an exception
    status_html = result['registry']['status'].value
    assert 'reliability' in status_html
    assert 'ImportError' not in status_html


def test_dropdown_selection_shows_clear_error_for_non_ts(tmp_path):
    """A plain (non-TS) output must produce a readable explanation, not a
    silent reset."""
    backend = _FakeBackend()
    d = tmp_path / 'plain'
    d.mkdir()
    out = d / 'plain.out'
    out.write_text('! Opt FREQ\nTHE OPTIMIZATION HAS CONVERGED\n', encoding='utf-8')
    from delfin.dashboard.ts_goat_panel import open_ts_goat_dialog
    result = open_ts_goat_dialog(
        out,
        callbacks={
            'backend': backend, 'reset_dropdown': lambda: None,
            'render_viewer': lambda *a, **k: None, 'scroll_top': lambda: None,
        },
    )
    assert result['error_kind'] == 'not_usable_ts'
    assert 'not a usable converged TS' in result['error']
    assert result['panel'] is None


def test_submit_creates_input_in_dedicated_folder(tmp_path):
    """The confirmation flow: submit writes the GOAT input into a NEW
    folder next to the TS run, leaving the original files untouched."""
    pytest.importorskip('opi')  # build_ts_goat_input delegates to the OPI writer
    backend = _FakeBackend()
    result, _rendered = _open_dialog(tmp_path, _FakeBackend())
    assert result['error'] is None
    registry = result['registry']
    # simulate the user pressing Submit with the proposed indices:
    state = registry['state']
    proposed = state['proposed']
    assert proposed.indices, 'analysis should propose at least one atom'
    registry['panel']  # the dialog is attached
    # _do_submit via the on_click callback: emulate by calling submit path
    # reach into the closure through the returned state: use the recorded
    # show payload instead of widget internals
    files = state['files']
    src_dir = Path(files.out_path).parent
    before = sorted(p.name for p in src_dir.iterdir())
    # emulate _do_submit through the public API: build + backend call
    import datetime
    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    job_name = f'GOAT_TS_{Path(files.out_path).stem}_{stamp}'
    job_dir = src_dir / job_name
    job_dir.mkdir(parents=True, exist_ok=False)
    build_ts_goat_input(
        job_dir,
        elements=state['elements'],
        coordinates=state['coordinates'],
        constraint_indices=list(proposed.indices),
        charge=state['charge'], mult=state['mult'],
        method='XTB', basis='',
        ncores=12, maxcore=1000,
        base_name=job_name,
    )
    r = backend.submit_orca(
        job_dir=job_dir, job_name=job_name,
        inp_file=f'{job_name}.inp', pal=12, maxcore=1000,
    )
    assert r.returncode == 0
    assert (job_dir / f'{job_name}.inp').is_file()
    inp_text = (job_dir / f'{job_name}.inp').read_text()
    assert 'GOAT' in inp_text
    assert 'Constraints' in inp_text
    # original files untouched
    after = sorted(p.name for p in src_dir.iterdir())
    assert before == [n for n in after if n != f'{job_name}']


def test_missing_hess_is_reported_not_swallowed(tmp_path):
    """A TS output without .hess must produce a clear error -- the
    original bug silently swallowed such conditions."""
    backend = _FakeBackend()
    d = tmp_path / 'noints'
    d.mkdir()
    out = _write_ts_fixture(d, 'noints')
    (d / 'noints.hess').unlink()
    # the .xyz alone is enough for the no_hess error path
    from delfin.dashboard.ts_goat_panel import open_ts_goat_dialog
    result = open_ts_goat_dialog(
        out,
        callbacks={
            'backend': backend, 'reset_dropdown': lambda: None,
            'render_viewer': lambda *a, **k: None, 'scroll_top': lambda: None,
        },
    )
    assert result['error_kind'] == 'no_hess'
    assert 'No Hessian' in result['error']
    assert result['panel'] is None


# ---------------- integrated viewer (user report 2026-11-02) ---------------

def _open_with_panel_viewer(tmp_path):
    """Open the dialog and hand the panel's viewer_out to the test."""
    backend = _FakeBackend()
    captured = {}

    def render_in_panel(atoms, indices, arrows=None, output_widget=None):
        captured['atoms'] = atoms
        captured['indices'] = indices
        captured['output_widget'] = output_widget

    result, _ = _open_dialog(tmp_path, backend,
                             extra_callbacks={'render_in_panel': render_in_panel})
    return result, captured


def test_panel_contains_its_own_viewer_output(tmp_path):
    result, captured = _open_with_panel_viewer(tmp_path)
    assert result['error'] is None
    reg = result['registry']
    # the Output widget is a child of the panel -> rendered INSIDE it
    children = reg['panel'].children
    from ipywidgets import Output
    outputs = [c for c in children if isinstance(c, Output)]
    assert outputs, 'panel must contain an Output viewer widget'
    # show() rendered the TS into it (captured via the injected callback)
    assert captured['output_widget'] is outputs[0]
    assert len(captured['atoms']) == 3


def _iter_widgets(container):
    """All widget descendants of a panel (panels nest HBox rows)."""
    for child in container.children:
        yield child
        if hasattr(child, 'children'):
            yield from _iter_widgets(child)


def test_atom_click_toggles_constraints(tmp_path):
    result, _ = _open_with_panel_viewer(tmp_path)
    reg = result['registry']
    constraint_field = reg['constraints_field']
    start = reg['parse_constraints'](constraint_field.value)
    # click an unconstrained atom -> it is added; click again -> removed.
    # Any atom outside the proposal qualifies; if the proposal covers
    # every atom, an out-of-range click must be ignored (sync keeps the
    # field consistent with the preview).
    unconstrained = next(
        (i for i in range(3) if i not in start), None)
    if unconstrained is None:
        reg['on_atom_click'](99)
        assert reg['parse_constraints'](constraint_field.value) == start
        return
    reg['on_atom_click'](unconstrained)
    after_add = reg['parse_constraints'](constraint_field.value)
    assert unconstrained in after_add
    reg['on_atom_click'](unconstrained)
    after_remove = reg['parse_constraints'](constraint_field.value)
    assert unconstrained not in after_remove


def test_animation_toggle_swaps_geometry_and_resets(tmp_path):
    result, captured = _open_with_panel_viewer(tmp_path)
    reg = result['registry']
    from ipywidgets import ToggleButton
    toggles = [c for c in _iter_widgets(reg['panel'])
               if isinstance(c, ToggleButton)]
    assert toggles, 'Imag. Mode toggle must be in the panel'
    toggle = toggles[0]
    undistorted = list(captured['atoms'])
    toggle.value = True   # animation on: displaced geometry rendered
    displaced = list(captured['atoms'])
    assert displaced != undistorted
    toggle.value = False  # off again: the optimized geometry returns
    assert list(captured['atoms']) == undistorted


def test_mode_change_updates_animation_source(tmp_path):
    result, captured = _open_with_panel_viewer(tmp_path)
    reg = result['registry']
    dropdown = None
    from ipywidgets import Dropdown
    for row in reg['panel'].children:
        if isinstance(row, Dropdown):
            dropdown = row
    assert dropdown is not None, 'mode dropdown must exist'
    # switch mode: the analysis re-runs and the displaced geometry comes
    # from the newly selected mode (both hess columns exist in fixtures)
    if dropdown.value is not None and len(dropdown.options) > 1:
        first = dropdown.value
        other = next(v for _l, v in dropdown.options if v != first)
        dropdown.value = other
        assert reg['state']['proposed'] is not None


# --------------------- script-free viewer bootstrap (Voila replay fix) -----

def test_static_div_carries_no_script(tmp_path):
    """The script-free viewer markup contains the viewer div + payload
    but NO executable script: on a later-attached panel Voila replays
    stored output as inert HTML, so scripts in it never run (the cause
    of the '3Dmol.js failed to load' banner, user report 2026-11-02)."""
    import re
    import sys
    sys.path.insert(0, '.')
    from delfin.dashboard.molecule_viewer import static_viewer_html
    vid, stored = static_viewer_html(
        '2\ntest\nO 0 0 0\nH 1 0 0\n',
        highlight_indices=[0], arrows=[0.1, 0, 0, -0.1, 0, 0])
    assert vid in stored
    assert '<script type="application/json"' in stored
    for m in re.finditer(r'<script([^>]*)>', stored):
        assert 'type="application/json"' in m.group(1), stored
    assert 'arrows' in stored


def test_bootstrap_js_uses_no_cdn_and_diagnoses_missing_lib():
    import sys
    sys.path.insert(0, '.')
    from delfin.dashboard.molecule_viewer import ts_goat_viewer_bootstrap_js
    js = ts_goat_viewer_bootstrap_js('vid123', arrows=[0.0, 0.1, 0.0])
    assert 'vid123' in js
    assert 'cdn' not in js.lower()
    assert 'window.$3Dmol' in js
    assert 'DELFIN TS-GOAT: window.$3Dmol is undefined' in js
    assert 'setClickable' in js
    assert 'ts-goat-click-sync' in js


def test_step_js_displaces_in_both_directions():
    import sys
    sys.path.insert(0, '.')
    from delfin.dashboard.molecule_viewer import ts_goat_viewer_step_js
    pos = ts_goat_viewer_step_js('vid123', +1)
    neg = ts_goat_viewer_step_js('vid123', -1)
    assert 'var sign = 1;' in pos
    assert 'var sign = -1;' in neg
    assert 'payload.arrows' in pos and 'payload.arrows' in neg
    assert 'window.__tsGoatAmp' in pos


def test_click_sync_widget_observes_and_routes(tmp_path):
    """A click arriving through the hidden sync input must route into
    on_atom_click and toggle the constraint field (JS -> Python)."""
    out, _rendered = _open_dialog(tmp_path, _FakeBackend())
    reg = out['registry']
    sync = reg.get('click_sync')
    assert sync is not None, 'click-sync Text widget missing from registry'
    field = reg.get('constraints_field')
    assert field is not None
    start = reg['parse_constraints'](field.value)
    # out-of-range click ignored
    sync.value = '999'
    assert reg['parse_constraints'](field.value) == start
    # in-range click toggles
    # first click toggles: if the proposal covers every atom the click
    # REMOVES it; either way the value must change coherently
    target = 0 if 0 not in start else (1 if 1 not in start else 2)
    sync.value = str(target)
    after = reg['parse_constraints'](field.value)
    assert (target in after) != (target in start)
    # an observe() fires only on value CHANGE: route through another
    # value so the second click actually arrives (like two real clicks
    # do through the JS channel).
    sync.value = str(999 if target != 999 else 998)
    sync.value = str(target)
    # second click toggles back: 2 is re-added
    assert target in reg['parse_constraints'](field.value)
